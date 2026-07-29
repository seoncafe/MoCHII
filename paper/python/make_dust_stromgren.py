"""H II region size versus dust treatment (MoCHII dusty Stromgren gate).

Overlays the radial neutral-fraction profile x_HI(r) for three dust
treatments of the same Stromgren configuration (uniform n_H = 100 sphere,
Q_H = 1e49 s^-1, T_star = 4e4 K, case B, fixed T_e = 1e4 K, level-6
octree), and marks each effective ionized radius
R_eff = (3 V_ion / 4 pi)^(1/3), V_ion = sum_cell x_HII V_cell:

  (a) dust-free  (ion_add_dust = .false.)          reused from
      tests/g1_stromgren/g1_stromgren_rates.h5
  (b) full MW dust (ion_add_dust = .true., D03 EUV extinction)
      tests/d_dusty/d_dust_full_rates.h5
  (c) laursen09_live (dust density tracks the computed x_HI, so dust
      vanishes in ionized gas)   tests/d_dusty/d_dust_l09_rates.h5

All three are the same uniform level-6 grid (64^3, box half 4 pc, cell
0.125 pc).  Cells outside the density sphere (center r > 4 pc, n_H = 0)
are masked, since their x_HI is not physically meaningful.

Output: figures/dust_stromgren.png ; run:  python3 make_dust_stromgren.py
"""

import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "figures")
ROOT = os.path.join(HERE, "..", "..")
OUT = os.path.join(FIGDIR, "dust_stromgren.pdf")

DCELL = 0.125            # pc, uniform level-6 cell size (box half 4 pc, 64^3)
RSPHERE = 4.0            # pc, density-sphere radius (center-based cut)

CASES = [
    ("dust-free", os.path.join(ROOT, "tests", "g1_stromgren",
                               "g1_stromgren_rates.h5"), "#000000", "-"),
    ("full MW dust", os.path.join(ROOT, "tests", "d_dusty",
                                  "d_dust_full_rates.h5"), "#d62728", "-"),
    (r"\texttt{laursen09\_live}", os.path.join(ROOT, "tests", "d_dusty",
                                  "d_dust_l09_rates.h5"), "#1f77b4", "--"),
]


def reff_from(xhi, r):
    m = r <= RSPHERE
    vion = np.sum((1.0 - xhi)[m]) * DCELL ** 3
    return (3.0 * vion / (4.0 * np.pi)) ** (1.0 / 3.0)


def radial_median(q, r, edges):
    out = np.full(len(edges) - 1, np.nan)
    for i in range(len(edges) - 1):
        m = (r >= edges[i]) & (r < edges[i + 1])
        if m.any():
            out[i] = np.median(q[m])
    return out


# ---------------------------------------------------------------- read
edges = np.linspace(0.0, RSPHERE, 65)
mid = 0.5 * (edges[:-1] + edges[1:])

profiles = []
r0 = None
for name, path, col, ls in CASES:
    with h5py.File(path, "r") as f:
        xhi = f["x_HI/data"][:]
        xyz = f["LeafXYZ/data"][:]
    r = np.sqrt((xyz ** 2).sum(axis=0))
    inside = r <= RSPHERE
    xhi_m = radial_median(np.where(inside, xhi, np.nan)[inside], r[inside],
                          edges)
    reff = reff_from(xhi, r)
    if r0 is None:
        r0 = reff
    profiles.append((name, col, ls, xhi_m, reff))
    print("%-16s R_eff = %.4f pc   R/R0 = %.4f" % (name, reff, reff / r0))

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 11, "axes.linewidth": 0.9})
# Paper figures: axis, tick, and title fonts sized to ~85% of the body text.
_S = 1.63
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "figure.titlesize": _b * 1.2 * _S,
})
fig, ax = plt.subplots(figsize=(7.2, 5.0))

for name, col, ls, xhi_m, reff in profiles:
    rr = reff / r0
    lab = (name + r',  $R_{\rm eff}=' + ("%.2f" % reff) +
           r'$~pc  ($R/R_0=' + ("%.3f" % rr) + r'$)')
    ax.plot(mid, xhi_m, ls, color=col, lw=2.0, label=lab)
    ax.axvline(reff, color=col, ls=":", lw=1.2, alpha=0.8)

ax.set_yscale("log")
ax.set_ylim(3e-5, 2.0)
ax.set_xlim(0.0, RSPHERE)
ax.set_xlabel(r'$r$~[pc]')
ax.set_ylabel(r'neutral fraction $x_{\rm HI}$')
ax.legend(loc="lower right", frameon=False, fontsize=9.5)
ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.4)
ax.set_title(r'H\,II region size versus dust treatment '
             r'($n_{\rm H}=100$~cm$^{-3}$, $Q_{\rm H}=10^{49}$~s$^{-1}$)',
             fontsize=11.5)

fig.tight_layout()
fig.savefig(OUT, dpi=180, bbox_inches="tight")
print("wrote", OUT)
