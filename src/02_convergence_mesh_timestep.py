#!/usr/bin/env python3
"""
Mesh and time-step convergence study for the WAAM thermal model.

Identical physics to waam_thermal.py. A reduced scenario (3 layers, fixed 10 s
interpass dwell) is run at several cell sizes and several time-step safety
factors, and the two quantities the report relies on -- peak temperature at the
deposition surface and instantaneous melt-pool length -- are recorded.
"""
import numpy as np, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(ROOT, "results", "data")
os.makedirs(OUT, exist_ok=True)

RHO, CP, K = 4430.0, 610.0, 7.2
LATENT, T_SOL, T_LIQ = 286e3, 1877.0, 1933.0
EMISS, SIGMA, H_CONV, T_AMB = 0.4, 5.670e-8, 15.0, 298.0
CURRENT, VOLTAGE, ETA = 160.0, 11.0, 0.62
V_TRAV, R_SRC, W_WALL = 4.5e-3, 3.0e-3, 6.0e-3
H_LAYER, L_WALL, SUB_T = 1.5e-3, 100.0e-3, 6.0e-3
N_LAYER, DWELL = 3, 10.0
Q_ARC = ETA * VOLTAGE * CURRENT
ALPHA = K / (RHO * CP)


def run(dx, cfl):
    dz = dx
    nx = int(round(L_WALL / dx))
    nzs = int(round(SUB_T / dz))
    rpl = int(round(H_LAYER / dz))
    nz = nzs + N_LAYER * rpl
    dt = cfl * dx**2 / ALPHA

    x = (np.arange(nx) + 0.5) * dx
    z = (np.arange(nz) + 0.5) * dz
    T = np.full((nz, nx), T_AMB)
    act = np.zeros((nz, nx), dtype=bool)
    act[:nzs, :] = True
    mon = nx // 2

    def step(Tf, a, q):
        Tp = np.where(a, Tf, T_AMB)
        Tn = np.empty_like(Tp); Tn[:-1] = Tp[1:];  Tn[-1] = Tp[-1]
        Ts = np.empty_like(Tp); Ts[1:]  = Tp[:-1]; Ts[0]  = Tp[0]
        Te = np.empty_like(Tp); Te[:, :-1] = Tp[:, 1:];  Te[:, -1] = Tp[:, -1]
        Tw = np.empty_like(Tp); Tw[:, 1:]  = Tp[:, :-1]; Tw[:, 0]  = Tp[:, 0]
        an  = np.empty_like(a); an[:-1] = a[1:];  an[-1] = False
        as_ = np.empty_like(a); as_[1:] = a[:-1]; as_[0] = False
        ae  = np.empty_like(a); ae[:, :-1] = a[:, 1:];  ae[:, -1] = False
        aw  = np.empty_like(a); aw[:, 1:]  = a[:, :-1]; aw[:, 0]  = False
        Tn = np.where(an, Tn, Tp); Ts = np.where(as_, Ts, Tp)
        Te = np.where(ae, Te, Tp); Tw = np.where(aw, Tw, Tp)
        lap = (Te - 2 * Tp + Tw) / dx**2 + (Tn - 2 * Tp + Ts) / dz**2
        dT = Tp - T_AMB
        rad = EMISS * SIGMA * (Tp**4 - T_AMB**4)
        qs = -(2.0 / W_WALL) * (H_CONV * dT + rad)
        qs += np.where(a & ~an, -(1.0 / dz) * (H_CONV * dT + rad), 0.0)
        c = np.full_like(Tp, CP)
        m = (Tp > T_SOL) & (Tp < T_LIQ)
        c[m] += LATENT / (T_LIQ - T_SOL)
        return np.clip(np.where(a, Tp + dt * (K * lap + q + qs) / (RHO * c),
                                T_AMB), T_AMB, 4000.0)

    peak_top, pool_top = None, None
    for layer in range(N_LAYER):
        r0 = nzs + layer * rpl
        r1 = r0 + rpl
        act[r0:r1, :] = True
        T[r0:r1, :] = T_AMB
        src = r1 - 1
        n_steps = int(np.ceil((L_WALL / V_TRAV) / dt))
        band = slice(max(r0 - int(round(3.0e-3 / dz)), 0), r1)
        zz_b = z[band]
        lp, pool_now = T_AMB, 0.0
        for n in range(n_steps):
            xs = V_TRAV * (n * dt)
            q = np.zeros_like(T)
            xx, zz = np.meshgrid(x, zz_b)
            d2 = (xx - xs)**2 + (zz - z[src])**2
            q[band] = (3.0 * Q_ARC / (np.pi * R_SRC**2 * W_WALL)
                       * np.exp(-3.0 * d2 / R_SRC**2))
            q[~act] = 0.0
            T = step(T, act, q)
            lp = max(lp, T[src, mon])
            if n == n_steps // 2:
                pool_now = float(np.sum(T[src] >= T_LIQ) * dx * 1000.0)
        t0 = 0.0
        while t0 < DWELL:
            T = step(T, act, np.zeros_like(T)); t0 += dt
        if layer == N_LAYER - 1:
            peak_top, pool_top = lp - 273.15, pool_now
    return dict(dx_mm=dx * 1000, cfl=cfl, dt_ms=dt * 1000,
                nx=nx, nz=nz, cells=nx * nz,
                peak_C=round(peak_top, 1), pool_mm=round(pool_top, 2))


print("=== MESH REFINEMENT (cfl = 0.12) ===")
mesh = [run(d, 0.12) for d in (0.75e-3, 0.5e-3, 0.25e-3)]
for r in mesh:
    print(f"  dx={r['dx_mm']:.2f} mm  cells={r['cells']:6d}  dt={r['dt_ms']:.2f} ms"
          f"  peak={r['peak_C']:.1f} C  pool={r['pool_mm']:.2f} mm")

print("=== TIME-STEP REFINEMENT (dx = 0.5 mm) ===")
tstep = [run(0.5e-3, c) for c in (0.12, 0.06, 0.03)]
for r in tstep:
    print(f"  cfl={r['cfl']:.3f}  dt={r['dt_ms']:.2f} ms"
          f"  peak={r['peak_C']:.1f} C  pool={r['pool_mm']:.2f} mm")

ref = mesh[-1]
for r in mesh:
    r["d_peak_pc"] = round(100 * (r["peak_C"] - ref["peak_C"]) / ref["peak_C"], 2)
    r["d_pool_pc"] = round(100 * (r["pool_mm"] - ref["pool_mm"]) / ref["pool_mm"], 2)
reft = tstep[-1]
for r in tstep:
    r["d_peak_pc"] = round(100 * (r["peak_C"] - reft["peak_C"]) / reft["peak_C"], 2)
    r["d_pool_pc"] = round(100 * (r["pool_mm"] - reft["pool_mm"]) / reft["pool_mm"], 2)

json.dump(dict(mesh=mesh, tstep=tstep),
          open(f"{OUT}/convergence.json", "w"), indent=2)
print("\nwritten:", f"{OUT}/convergence.json")
