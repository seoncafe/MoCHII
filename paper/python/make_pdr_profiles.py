"""Radial profiles of the MoCHII PDR (photodissociation-region) test.

Reads the converged level-5 octree rates file
``tests/pdr/pdr_L5_rates.h5`` and plots, versus radius r [pc] from the
central source at the box center:

  top panel  -- ionization fractions x_HI (neutral H) and x_CII (C II),
                showing the fully ionized H II interior, the sharp
                ionization front, and x_CII -> 1 in the PDR shell;
  bottom     -- electron temperature T_e [K] and electron density n_e
                [cm^-3], showing the PDR shell (n_e ~ A_C n_H set by
                carbon, T_e ~ 265 K with a gradient rather than pinned
                at the 100 K floor).

Output: figures/pdr_profiles.png

Run:  python3 make_pdr_profiles.py
"""

import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "figures")
RUN = os.path.join(HERE, "..", "..", "tests", "pdr", "pdr_L5_rates.h5")
OUT = os.path.join(FIGDIR, "pdr_profiles.pdf")


def rd(f, name):
    return f[name + "/data"][:]


# ---------------------------------------------------------------- read
with h5py.File(RUN, "r") as f:
    xyz = rd(f, "LeafXYZ")          # (3, nleaf), code units = pc here
    ne = rd(f, "n_e")
    te = rd(f, "T_e")
    xhi = rd(f, "x_HI")
    xc = rd(f, "x_c_stages")        # (nleaf, 4): C I, C II, C III, C IV
    if xc.shape[0] != xyz.shape[1]:
        xc = xc.T

r = np.sqrt((xyz ** 2).sum(axis=0))     # pc (source at box center)
xcii = xc[:, 1]                          # C II fraction

# ------------------------------------------------ radial median binning
rmax = 4.0
nbin = 48
edges = np.linspace(0.0, rmax, nbin + 1)
mid = 0.5 * (edges[:-1] + edges[1:])


def radial_stat(q, fn):
    out = np.full(nbin, np.nan)
    for i in range(nbin):
        m = (r >= edges[i]) & (r < edges[i + 1])
        if m.any():
            out[i] = fn(q[m])
    return out


xhi_med = radial_stat(xhi, np.median)
xcii_med = radial_stat(xcii, np.median)
te_med = radial_stat(te, np.median)
ne_med = radial_stat(ne, np.median)
# 5-95 percentile band for T_e (octree scatter within a shell)
te_lo = radial_stat(te, lambda a: np.percentile(a, 5))
te_hi = radial_stat(te, lambda a: np.percentile(a, 95))

# ------------------------------------------------ figure
plt.rcParams.update({"font.size": 11, "axes.linewidth": 0.9})
# Paper figures: axis, tick, and title fonts sized to ~85% of the body text.
_S = 1.58
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "figure.titlesize": _b * 1.2 * _S,
})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.0, 7.2), sharex=True)

col_hi = "#1f77b4"
col_cii = "#d62728"
col_te = "#c1440e"
col_ne = "#1f6f4a"

# --- panel 1: ionization fractions -----------------------------------
ax1.plot(mid, xhi_med, "-", color=col_hi, lw=2.0,
         label=r'$x_{\rm HI}$ (neutral H)')
ax1.plot(mid, xcii_med, "--", color=col_cii, lw=2.0,
         label=r'$x_{\rm CII}$ (C\,II)')
ax1.set_yscale("log")
ax1.set_ylim(3e-5, 3.0)
ax1.set_ylabel(r'ionization fraction')
ax1.legend(loc="lower right", frameon=False, fontsize=10)
ax1.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
ax1.text(0.03, 0.90,
         r'H\,II interior $\;\rightarrow\;$ I-front $\;\rightarrow\;$ PDR shell',
         transform=ax1.transAxes, fontsize=10, va="top")

# --- panel 2: T_e (left) and n_e (right) -----------------------------
ax2.fill_between(mid, te_lo, te_hi, color=col_te, alpha=0.15, lw=0)
l_te, = ax2.plot(mid, te_med, "-", color=col_te, lw=2.0,
                 label=r'$T_{\rm e}$')
ax2.set_yscale("log")
ax2.set_ylabel(r'$T_{\rm e}$~[K]', color=col_te)
ax2.tick_params(axis="y", colors=col_te)
ax2.set_ylim(1e2, 2e4)
ax2.set_xlabel(r'$r$~[pc]')
ax2.set_xlim(0.0, rmax)
ax2.grid(True, which="major", ls=":", lw=0.5, alpha=0.5)

ax2r = ax2.twinx()
l_ne, = ax2r.plot(mid, ne_med, "--", color=col_ne, lw=2.0,
                  label=r'$n_{\rm e}$')
ax2r.set_yscale("log")
ax2r.set_ylabel(r'$n_{\rm e}$~[cm$^{-3}$]', color=col_ne)
ax2r.tick_params(axis="y", colors=col_ne)
ax2r.set_ylim(1e-3, 3e2)

# annotate the carbon-set PDR electron floor
ax2r.axhline(2.2e-4 * 100.0, color=col_ne, ls=":", lw=1.0, alpha=0.7)
ax2r.text(3.9, 2.2e-4 * 100.0 * 1.3, r'$A_{\rm C}\,n_{\rm H}$',
          color=col_ne, fontsize=9, ha="right", va="bottom")

ax2.legend([l_te, l_ne], [r'$T_{\rm e}$', r'$n_{\rm e}$'],
           loc="upper right", frameon=False, fontsize=10)

fig.suptitle(r'Photodissociation-region model: radial profiles '
             r'($n_{\rm H}=100$~cm$^{-3}$ sphere)', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.98))
fig.savefig(OUT, dpi=180, bbox_inches="tight")
print("wrote", OUT)
