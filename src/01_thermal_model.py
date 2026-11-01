#!/usr/bin/env python3
"""
Transient thermal model of a WAAM thin wall (Ti-6Al-4V).

2-D explicit finite-difference solution of the heat-conduction equation in the
longitudinal section (x = travel direction, z = build direction), with:
  * a moving Gaussian heat source,
  * layer-by-layer element activation,
  * convective and radiative loss from the side faces and the free top face,
  * latent heat of fusion by the effective-heat-capacity method.

Two interpass policies are simulated and compared:
  'controlled' - dwell until the peak wall temperature falls below 200 C
  'short'      - fixed 10 s dwell (no interpass temperature control)
"""
import numpy as np, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(ROOT, "results", "data")
os.makedirs(OUT, exist_ok=True)

# ---------------- Material: Ti-6Al-4V ----------------
RHO, CP, K = 4430.0, 610.0, 7.2          # kg/m3, J/kg.K, W/m.K
LATENT = 286e3                            # J/kg
T_SOL, T_LIQ, T_BETA = 1877.0, 1933.0, 1268.0   # K
EMISS, SIGMA, H_CONV, T_AMB = 0.4, 5.670e-8, 15.0, 298.0

# ---------------- Process (GTAW-WAAM, cf. Table 1) ----------------
CURRENT, VOLTAGE, ETA = 160.0, 11.0, 0.62
V_TRAV  = 4.5e-3        # m/s
R_SRC   = 3.0e-3        # m
W_WALL  = 6.0e-3        # m
H_LAYER = 1.5e-3        # m
L_WALL  = 100.0e-3      # m
N_LAYER = 8
T_INTER = 473.0         # K (200 C)
DWELL_MAX, DWELL_SHORT = 150.0, 10.0

Q_ARC   = ETA * VOLTAGE * CURRENT               # W
H_INPUT = Q_ARC / (V_TRAV * 1000.0)             # J/mm

# ---------------- Grid ----------------
DX = DZ = 0.5e-3
NX = int(round(L_WALL / DX))
NZ_SUB = int(round(6.0e-3 / DZ))
NZ = NZ_SUB + int(round(N_LAYER * H_LAYER / DZ))
ALPHA = K / (RHO * CP)
DT = 0.12 * DX**2 / ALPHA

x = (np.arange(NX) + 0.5) * DX
z = (np.arange(NZ) + 0.5) * DZ
rows_per_layer = int(round(H_LAYER / DZ))
mon_row = NZ_SUB + 2 * rows_per_layer + rows_per_layer // 2
mon_col = NX // 2


def eff_cp(Tf):
    c = np.full_like(Tf, CP)
    mushy = (Tf > T_SOL) & (Tf < T_LIQ)
    c[mushy] += LATENT / (T_LIQ - T_SOL)
    return c


def step(Tf, act, qvol):
    Tp = np.where(act, Tf, T_AMB)
    Tn = np.empty_like(Tp); Tn[:-1] = Tp[1:];  Tn[-1] = Tp[-1]
    Ts = np.empty_like(Tp); Ts[1:]  = Tp[:-1]; Ts[0]  = Tp[0]
    Te = np.empty_like(Tp); Te[:, :-1] = Tp[:, 1:];  Te[:, -1] = Tp[:, -1]
    Tw = np.empty_like(Tp); Tw[:, 1:]  = Tp[:, :-1]; Tw[:, 0]  = Tp[:, 0]

    an  = np.empty_like(act); an[:-1] = act[1:];  an[-1] = False
    as_ = np.empty_like(act); as_[1:] = act[:-1]; as_[0] = False
    ae  = np.empty_like(act); ae[:, :-1] = act[:, 1:];  ae[:, -1] = False
    aw  = np.empty_like(act); aw[:, 1:]  = act[:, :-1]; aw[:, 0]  = False

    Tn = np.where(an, Tn, Tp); Ts = np.where(as_, Ts, Tp)
    Te = np.where(ae, Te, Tp); Tw = np.where(aw, Tw, Tp)

    lap = (Te - 2 * Tp + Tw) / DX**2 + (Tn - 2 * Tp + Ts) / DZ**2
    dT = Tp - T_AMB
    rad = EMISS * SIGMA * (Tp**4 - T_AMB**4)
    q_surf = -(2.0 / W_WALL) * (H_CONV * dT + rad)
    top_face = act & ~an
    q_surf += np.where(top_face, -(1.0 / DZ) * (H_CONV * dT + rad), 0.0)

    dTdt = (K * lap + qvol + q_surf) / (RHO * eff_cp(Tp))
    return np.clip(np.where(act, Tp + DT * dTdt, T_AMB), T_AMB, 4000.0)


def run(mode):
    T = np.full((NZ, NX), T_AMB)
    act = np.zeros((NZ, NX), dtype=bool)
    act[:NZ_SUB, :] = True

    hist_t, hist_T = [], []
    peak, pool, t85, gr = [], [], [], []
    field = None
    _ = None
    tc = 0.0

    for layer in range(N_LAYER):
        r0 = NZ_SUB + layer * rows_per_layer
        r1 = r0 + rows_per_layer
        act[r0:r1, :] = True
        T[r0:r1, :] = T_AMB
        src = r1 - 1

        n_steps = int(np.ceil((L_WALL / V_TRAV) / DT))
        layer_peak = T_AMB
        pool_now, prev_front = 0.0, None
        band = slice(max(r0 - 6, 0), r1)
        zz_band = z[band]

        for n in range(n_steps):
            xs = V_TRAV * (n * DT)
            qvol = np.zeros_like(T)
            xx, zz = np.meshgrid(x, zz_band)
            d2 = (xx - xs)**2 + (zz - z[src])**2
            qvol[band] = (3.0 * Q_ARC / (np.pi * R_SRC**2 * W_WALL)
                          * np.exp(-3.0 * d2 / R_SRC**2))
            qvol[~act] = 0.0

            T = step(T, act, qvol)
            tc += DT
            layer_peak = max(layer_peak, T[src, mon_col])

            if n % 5 == 0:
                hist_t.append(tc); hist_T.append(T[mon_row, mon_col])
            if n == n_steps // 2:
                pool_now = float(np.sum(T[src] >= T_LIQ) * DX * 1000.0)
                if layer == 4:
                    field = T.copy()

            if 0.30 * n_steps < n < 0.75 * n_steps and n % 12 == 0:
                Ttop = T[src]
                liq = np.where(Ttop >= T_LIQ)[0]
                if liq.size > 2 and liq.min() > 2:
                    i = int(liq.min())
                    # linear interpolation of the T_LIQ crossing between i-1 and i
                    dT_cell = Ttop[i] - Ttop[i - 1]
                    frac = (T_LIQ - Ttop[i - 1]) / dT_cell if abs(dT_cell) > 1e-9 else 0.0
                    front = (i - 1 + np.clip(frac, 0.0, 1.0)) * DX
                    G = abs(Ttop[i + 1] - Ttop[i - 1]) / (2 * DX)
                    if prev_front is not None:
                        R = abs(front - prev_front) / (12 * DT)
                        if 1e-4 < R < 0.05 and G > 1e3:
                            gr.append((G, R))
                    prev_front = front

        peak.append(layer_peak - 273.15)
        pool.append(pool_now)

        # ---- interpass dwell ----
        t0, hi, t_hi, got = tc, False, None, False
        cap = DWELL_MAX if mode == "controlled" else DWELL_SHORT
        while tc - t0 < cap:
            T = step(T, act, np.zeros_like(T))
            tc += DT
            if int(tc / DT) % 5 == 0:
                hist_t.append(tc); hist_T.append(T[mon_row, mon_col])
            Tm = T[r1 - 1, mon_col]
            if not hi and Tm <= 1073.0:
                hi, t_hi = True, tc
            if hi and not got and Tm <= 773.0:
                t85.append(tc - t_hi); got = True
            if mode == "controlled" and T[act].max() < T_INTER:
                break
        if not got:
            t85.append(float("nan"))

    return dict(t=np.array(hist_t), T=np.array(hist_T), peak=np.array(peak),
                pool=np.array(pool), t85=np.array(t85),
                gr=np.array(gr) if gr else np.zeros((0, 2)),
                field=field, act=act.copy(), total=tc)


res = {m: run(m) for m in ("controlled", "short")}

# ---------------- Export ----------------
for m, r in res.items():
    k = np.linspace(0, len(r["t"]) - 1, min(1200, len(r["t"]))).astype(int)
    np.savetxt(f"{OUT}/cycle_{m}.dat", np.c_[r["t"][k], r["T"][k] - 273.15],
               header="t_s T_C", comments="")
    L = np.arange(1, N_LAYER + 1)
    np.savetxt(f"{OUT}/layer_{m}.dat", np.c_[L, r["peak"], r["pool"], r["t85"]],
               header="layer Tpeak_C pool_mm t85_s", comments="")
    if len(r["gr"]):
        np.savetxt(f"{OUT}/gr_{m}.dat", r["gr"], header="G_Kperm R_mpers", comments="")

np.save(f"{OUT}/field.npy", res["controlled"]["field"])
np.save(f"{OUT}/active.npy", res["controlled"]["act"])

g = res["controlled"]["gr"]
meta = dict(
    Q_ARC=Q_ARC, H_INPUT_J_per_mm=round(H_INPUT, 1), DT=DT, DX=DX,
    NX=NX, NZ=NZ, NZ_SUB=NZ_SUB, rows_per_layer=rows_per_layer,
    N_LAYER=N_LAYER, alpha=ALPHA, V_TRAV=V_TRAV,
    controlled=dict(peak_first=round(res["controlled"]["peak"][0], 1),
                    peak_last=round(res["controlled"]["peak"][-1], 1),
                    pool_first=res["controlled"]["pool"][0],
                    pool_last=res["controlled"]["pool"][-1],
                    total_s=round(res["controlled"]["total"], 1)),
    short=dict(peak_first=round(res["short"]["peak"][0], 1),
               peak_last=round(res["short"]["peak"][-1], 1),
               pool_first=res["short"]["pool"][0],
               pool_last=res["short"]["pool"][-1],
               total_s=round(res["short"]["total"], 1)),
    G_median=float(np.median(g[:, 0])) if len(g) else None,
    R_median=float(np.median(g[:, 1])) if len(g) else None,
    GR_median=float(np.median(g[:, 0] / g[:, 1])) if len(g) else None,
    n_gr=int(len(g)))
json.dump(meta, open(f"{OUT}/meta.json", "w"), indent=2)
print(json.dumps(meta, indent=2))
