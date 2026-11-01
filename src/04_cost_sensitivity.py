#!/usr/bin/env python3
"""
One-at-a-time sensitivity analysis of the techno-economic model of Section 7.

Each input is swept across a plausible range with all others held at baseline.
Reported: the single-unit (N=1) WAAM total, the margin against sand casting,
and the WAAM-versus-casting crossover mass with a standalone qualification
campaign charged.
"""
import json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(ROOT, "results", "data")
os.makedirs(OUT, exist_ok=True)

BASE = dict(mf=180.0, m_waam=230.0, m_pour=290.0, btf_b=5.0, phi=0.10,
            c_wire=25.0, c_ingot=12.0, mdot_dep=3.0, mdot_mrr=4.0,
            r_cell=90.0, r_cnc=75.0, t_mach_w=35.0, t_mach_c=30.0,
            c_fdy=2400.0, c_tool=18000.0, c_qual_w=14000.0, c_qual_c=1200.0)


def costs(p, mf=None):
    mf = p["mf"] if mf is None else mf
    k = mf / p["mf"]                      # geometric scaling of the baseline part
    waam_mat = p["m_waam"] * k * p["c_wire"]
    waam_dep = (p["m_waam"] * k / p["mdot_dep"]) * p["r_cell"]
    waam_cnc = p["t_mach_w"] * k * p["r_cnc"]
    W = waam_mat + waam_dep + waam_cnc
    net = mf + p["phi"] * (p["m_pour"] * k - mf)
    C = net * p["c_ingot"] + p["c_fdy"] * k + p["t_mach_c"] * k * p["r_cnc"]
    B = (mf * p["btf_b"] * p["c_ingot"]
         + (mf * (p["btf_b"] - 1) / p["mdot_mrr"]) * p["r_cnc"])
    return W, C, B


def totals(p, mf=None):
    W, C, B = costs(p, mf)
    return (W + p["c_qual_w"],
            C + p["c_tool"] + p["c_qual_c"],
            B + p["c_qual_c"])


def crossover(p):
    """Finished mass at which WAAM total equals sand-casting total, N = 1."""
    lo, hi = 10.0, 3000.0
    f = lambda m: totals(p, m)[0] - totals(p, m)[1]
    if f(lo) * f(hi) > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


bW, bC, bB = totals(BASE)
bX = crossover(BASE)
print(f"BASELINE  WAAM {bW:,.0f} | cast {bC:,.0f} | billet {bB:,.0f} "
      f"| margin {bW-bC:+,.0f} | crossover {bX:.0f} kg\n")

SWEEPS = [
    ("WAAM qualification campaign", "c_qual_w", 8000.0, 25000.0, "USD"),
    ("WAAM deposition rate",        "mdot_dep",     2.0,     4.0, "kg/h"),
    ("Wire cost",                   "c_wire",      20.0,    30.0, "USD/kg"),
    ("WAAM cell rate",              "r_cell",      70.0,   110.0, "USD/h"),
    ("Casting tooling",             "c_tool",   12000.0, 24000.0, "USD"),
]

rows = []
print(f"{'parameter':30} {'low':>10} {'high':>10} "
      f"{'WAAM lo':>10} {'WAAM hi':>10} {'margin lo':>10} {'margin hi':>10} "
      f"{'xover lo':>9} {'xover hi':>9}")
print("-" * 118)
for name, key, lo, hi, unit in SWEEPS:
    out = {}
    for tag, v in (("lo", lo), ("hi", hi)):
        p = dict(BASE); p[key] = v
        W, C, B = totals(p)
        out[tag] = dict(val=v, waam=W, cast=C, margin=W - C, xover=crossover(p))
    # Ranked by influence on the MARGIN, not on the WAAM total: casting
    # tooling moves the comparison without touching the WAAM side at all.
    rows.append(dict(name=name, key=key, unit=unit,
                     lo=out["lo"], hi=out["hi"],
                     swing=abs(out["hi"]["margin"] - out["lo"]["margin"])))
    print(f"{name:30} {lo:>10,.0f} {hi:>10,.0f} "
          f"{out['lo']['waam']:>10,.0f} {out['hi']['waam']:>10,.0f} "
          f"{out['lo']['margin']:>+10,.0f} {out['hi']['margin']:>+10,.0f} "
          f"{out['lo']['xover'] or float('nan'):>9.0f} "
          f"{out['hi']['xover'] or float('nan'):>9.0f}")

rows.sort(key=lambda r: r["swing"], reverse=True)
print("\nranked by influence on the single-unit cost margin:")
for r in rows:
    print(f"  {r['name']:30} swing {r['swing']:>8,.0f} USD")

# the threshold that decides the headline conclusion
flip = bC - costs(BASE)[0]
print(f"\nThe N=1 verdict reverses when the qualification campaign falls below "
      f"{flip:,.0f} USD (baseline assumption 14,000).")

json.dump(dict(baseline=dict(waam=bW, cast=bC, billet=bB,
                             margin=bW - bC, crossover=bX),
               sweeps=rows, flip_threshold=flip),
          open(f"{OUT}/sensitivity.json", "w"), indent=2, default=float)
print("written:", f"{OUT}/sensitivity.json")
