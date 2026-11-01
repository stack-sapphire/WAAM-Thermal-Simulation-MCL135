#!/usr/bin/env python3
"""Generate the thermal-model figures as PDFs in the report's palette."""
import numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "results", "data")
F = os.path.join(ROOT, "results", "figures")
os.makedirs(F, exist_ok=True)
M = json.load(open(f"{D}/meta.json"))

STEEL, RUST, TEAL, GRAPH = "#2A3D4D", "#A54A2A", "#3E6C6A", "#5A626A"
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"], "font.size": 9,
    "axes.edgecolor": STEEL, "axes.labelcolor": STEEL, "axes.labelsize": 9.5,
    "xtick.color": GRAPH, "ytick.color": GRAPH, "axes.linewidth": 0.8,
    "legend.frameon": True, "legend.edgecolor": "#C8CCD1", "legend.fontsize": 8,
    "grid.color": "#C8CCD1", "grid.linewidth": 0.5,
})

# ---------------------------------------------------------------- Fig A
c = np.loadtxt(f"{D}/cycle_controlled.dat", skiprows=1)
fig, ax = plt.subplots(figsize=(6.4, 3.1))
ax.plot(c[:, 0], c[:, 1], color=STEEL, lw=1.4,
        label="Computed temperature at the monitor point")
ax.axhline(1660, color=RUST, ls="--", lw=1.0, label="Liquidus $T_L$ = 1660 $^\\circ$C")
ax.axhline(995,  color=TEAL, ls="-.", lw=1.0, label="$\\beta$-transus = 995 $^\\circ$C")
ax.axhline(200,  color=GRAPH, ls=":", lw=1.0, label="Interpass target = 200 $^\\circ$C")
ax.set_xlabel("Time from start of build, $t$ (s)")
ax.set_ylabel("Temperature (" + "$^\\circ$C)")
ax.set_xlim(320, c[-1, 0]); ax.set_ylim(0, 3100)
ax.grid(alpha=0.6); ax.legend(loc="upper right", ncol=1)
fig.tight_layout(); [fig.savefig(f"{F}/fig_cycle.%s" % _e, dpi=200) for _e in ("pdf", "png")]; plt.close(fig)

# ---------------------------------------------------------------- Fig B
lc = np.loadtxt(f"{D}/layer_controlled.dat", skiprows=1)
ls = np.loadtxt(f"{D}/layer_short.dat", skiprows=1)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.75))

a1.plot(ls[:, 0], ls[:, 1], "o-", color=RUST, lw=1.4, ms=4,
        label="10 s dwell (no control)")
a1.plot(lc[:, 0], lc[:, 1], "s-", color=STEEL, lw=1.4, ms=4,
        label="Dwell to 200 $^\\circ$C")
a1.set_xlabel("Layer number"); a1.set_ylabel("Peak temperature (" + "$^\\circ$C)")
a1.xaxis.set_major_locator(MaxNLocator(integer=True))
a1.grid(alpha=0.6); a1.legend(loc="lower right")

a2.plot(ls[:, 0], ls[:, 2], "o-", color=RUST, lw=1.4, ms=4)
a2.plot(lc[:, 0], lc[:, 2], "s-", color=STEEL, lw=1.4, ms=4)
a2.set_xlabel("Layer number"); a2.set_ylabel("Melt-pool length (mm)")
a2.xaxis.set_major_locator(MaxNLocator(integer=True))
a2.grid(alpha=0.6)
fig.tight_layout(); [fig.savefig(f"{F}/fig_accum.%s" % _e, dpi=200) for _e in ("pdf", "png")]; plt.close(fig)

# ---------------------------------------------------------------- Fig C
T = np.load(f"{D}/field.npy") - 273.15
act = np.load(f"{D}/active.npy")
Tm = np.ma.masked_where(~act, T)
NZ_SUB, rpl, DX = M["NZ_SUB"], M["rows_per_layer"], M["DX"]
nz, nx = T.shape
ext = [0, nx * DX * 1000, -NZ_SUB * DX * 1000, (nz - NZ_SUB) * DX * 1000]
z_top = (5 * rpl) * DX * 1000 + 1.0   # crop just above the layer being built

fig, ax = plt.subplots(figsize=(6.6, 2.5))
im = ax.imshow(Tm, origin="lower", extent=ext, aspect="auto",
               cmap="inferno", vmin=25, vmax=2600)
cs = ax.contour(np.linspace(ext[0], ext[1], nx),
                np.linspace(ext[2], ext[3], nz), T,
                levels=[200, 995, 1660], colors=["#FFFFFF", "#9ad7d3", "#ffffff"],
                linewidths=[0.6, 0.8, 1.1])
ax.clabel(cs, fmt="%d", fontsize=6)
ax.axhline(0, color="w", lw=0.8)
ax.set_ylim(ext[2], z_top)
ax.set_xlabel("Travel direction, $x$ (mm)")
ax.set_ylabel("Build direction, $z$ (mm)")
cb = fig.colorbar(im, ax=ax, pad=0.015)
cb.set_label("Temperature (" + "$^\\circ$C)", color=STEEL, fontsize=9)
cb.outline.set_edgecolor(STEEL)
fig.tight_layout(); [fig.savefig(f"{F}/fig_field.%s" % _e, dpi=300) for _e in ("pdf", "png")]; plt.close(fig)

# ---------------------------------------------------------------- Fig D
g = np.loadtxt(f"{D}/gr_controlled.dat", skiprows=1)
G, R = g[:, 0], g[:, 1]
fig, ax = plt.subplots(figsize=(4.6, 3.3))
ax.loglog(R, G, "o", ms=4.5, mfc=STEEL, mec="white", mew=0.4, alpha=0.85,
          label="Computed front conditions")
Rl = np.logspace(-4, -1, 50)
ax.loglog(Rl, 1e5 * Rl**0.0 * (Rl / Rl), color=RUST, ls="--", lw=0)   # placeholder
# columnar / equiaxed transition: G/R = 1e6 K.s/m^2  (Kobryn-Semiatin form)
ax.loglog(Rl, 1e6 * Rl, color=RUST, ls="--", lw=1.3,
          label="CET boundary, $G/R = 10^{6}$ K$\\cdot$s/m$^2$")
ax.fill_between(Rl, 1e6 * Rl, 1e8, color=STEEL, alpha=0.07)
ax.text(2.2e-3, 8.0e4, "columnar", color=STEEL, fontsize=8.5, style="italic")
ax.text(2.2e-3, 6.0e2, "equiaxed", color=GRAPH, fontsize=8.5, style="italic")
ax.set_xlabel("Solidification rate, $R$ (m/s)")
ax.set_ylabel("Thermal gradient, $G$ (K/m)")
ax.set_xlim(5e-4, 2e-2); ax.set_ylim(1e2, 1e6)
ax.grid(alpha=0.5, which="both"); ax.legend(loc="lower right")
fig.tight_layout(); [fig.savefig(f"{F}/fig_gr.%s" % _e, dpi=200) for _e in ("pdf", "png")]; plt.close(fig)

print("figures written to", F)
for f in sorted(os.listdir(F)):
    print(" ", f, os.path.getsize(os.path.join(F, f)), "bytes")
