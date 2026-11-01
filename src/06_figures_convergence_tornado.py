#!/usr/bin/env python3
"""Convergence and sensitivity figures, in the report's palette."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "results", "data")
F = os.path.join(ROOT, "results", "figures")
os.makedirs(F, exist_ok=True)
STEEL, RUST, TEAL, GRAPH = "#2A3D4D", "#A54A2A", "#3E6C6A", "#5A626A"
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"], "font.size": 9,
    "axes.edgecolor": STEEL, "axes.labelcolor": STEEL, "axes.labelsize": 9.5,
    "xtick.color": GRAPH, "ytick.color": GRAPH, "axes.linewidth": 0.8,
    "legend.frameon": True, "legend.edgecolor": "#C8CCD1", "legend.fontsize": 8,
    "grid.color": "#C8CCD1", "grid.linewidth": 0.5,
})

# ---------------------------------------------------------------- convergence
cv = json.load(open(f"{D}/convergence.json"))
df = json.load(open(f"{D}/convergence_diff.json"))

fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(7.0, 2.6))

dx = [r["dx_mm"] for r in cv["mesh"]]
a1.plot(dx, [r["peak_C"] for r in cv["mesh"]], "o-", color=STEEL, lw=1.4, ms=4.5)
a1.set_xlabel("Cell size (mm)"); a1.set_ylabel("Peak temperature ($^\\circ$C)")
a1.invert_xaxis(); a1.grid(alpha=0.6)
a1.set_title("Mesh refinement", fontsize=9, color=STEEL)

dt = [r["dt_ms"] for r in cv["tstep"]]
a2.plot(dt, [r["peak_C"] for r in cv["tstep"]], "s-", color=TEAL, lw=1.4, ms=4.5)
a2.set_xlabel("Time step (ms)"); a2.set_ylabel("Peak temperature ($^\\circ$C)")
a2.invert_xaxis(); a2.grid(alpha=0.6)
a2.set_title("Time-step refinement", fontsize=9, color=STEEL)

lab, sht, ctl = [], [], []
for d in (0.5, 0.25):
    lab.append(f"{d:g} mm")
    sht.append(next(r["growth_pc"] for r in df if r["dx_mm"] == d and r["mode"] == "short"))
    ctl.append(next(r["growth_pc"] for r in df if r["dx_mm"] == d and r["mode"] == "controlled"))
i = np.arange(len(lab)); w = 0.35
a3.bar(i - w/2, sht, w, color=RUST, label="no control")
a3.bar(i + w/2, ctl, w, color=STEEL, label="controlled")
a3.set_xticks(i); a3.set_xticklabels(lab)
a3.set_xlabel("Cell size"); a3.set_ylabel("Pool growth, L1$\\rightarrow$L4 (\\%)")
a3.grid(alpha=0.6, axis="y"); a3.legend(loc="upper left", fontsize=7.5)
a3.set_title("Differential result", fontsize=9, color=STEEL)

fig.tight_layout(); [fig.savefig(f"{F}/fig_convergence.%s" % _e, dpi=200) for _e in ("pdf", "png")]; plt.close(fig)

# ---------------------------------------------------------------- tornado
sn = json.load(open(f"{D}/sensitivity.json"))
base, sweeps = sn["baseline"], sn["sweeps"]
order = sorted(sweeps, key=lambda s: abs(s["hi"]["margin"] - s["lo"]["margin"]))
names = [s["name"] for s in order]
los = [s["lo"]["margin"] for s in order]
his = [s["hi"]["margin"] for s in order]
b = base["margin"]

fig, ax = plt.subplots(figsize=(6.5, 2.9))
for i, (lo, hi) in enumerate(zip(los, his)):
    left, right = min(lo, hi), max(lo, hi)
    ax.barh(i, left - b, left=b, height=0.55, color=RUST, alpha=0.9)
    ax.barh(i, right - b, left=b, height=0.55, color=STEEL, alpha=0.9)
ax.axvline(b, color=GRAPH, lw=1.2)
ax.axvline(0, color=TEAL, lw=1.5, ls="--")
ax.set_yticks(np.arange(len(names))); ax.set_yticklabels(names, fontsize=8.5)
ax.set_xlabel("WAAM total minus sand-casting total at $N=1$ (USD)")
ax.grid(alpha=0.5, axis="x")
ax.text(-200, len(names) - 0.55, "verdict flips here", color=TEAL, fontsize=8,
        va="top", ha="right")
ax.set_xlim(-5000, 16000)
fig.tight_layout(); [fig.savefig(f"{F}/fig_tornado.%s" % _e, dpi=200) for _e in ("pdf", "png")]; plt.close(fig)

print("written fig_convergence.pdf, fig_tornado.pdf")
