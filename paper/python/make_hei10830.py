#!/usr/bin/env python3
"""He I 2^3S metastable / 10830 A diagnostic maps for the paper.

Reads the high-density He I metastable example
(examples/hei_metastable/hei_high_density.in: an n_H = 1e4 cm^-3, R = 0.15 pc
sphere at a fixed T_e = 1e4 K, par%hei_metastable) through the HeiMetaData
reader in tools/python/mochii_output.py, and makes the two observable maps of
the 2^3S reservoir:

  (a) the 2^3S column density N(2^3S) [cm^-2] (line-of-sight integral of the
      metastable density);
  (b) the 10830 A line-center optical depth tau_0 of the blended J'=1,2
      doublet, with the tau_0 = 1 contour -- the core is optically THICK, so
      the optically-thin diagnostic (its assumption) breaks down there, which
      the reader also flags at run time.

The 10830 opacity uses the fine-structure constants (vacuum wavelength,
oscillator strength, Einstein A) from data/atomic/hei_metastable.txt, the same
file the Fortran diagnostic reads.

Output: figures/hei10830.pdf   (vector; the imshow rasters live in the PDF)

Data (existing example output, no MoCHII re-run):
  ../../examples/hei_metastable/hei_high_density_heimeta.h5

Run:  python3 make_hei10830.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
_HERE_AX = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE_AX)
from axis_style import square_ticks   # noqa: E402


plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"

# Axis, tick, and title fonts sized to ~85% of the body text once the two-panel
# map row is scaled down by \includegraphics (0.95\textwidth) in the paper.
_S = 1.35
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "figure.titlesize": _b * 1.2 * _S,
})

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "python"))
from mochii_output import HeiMetaData             # noqa: E402

FIGDIR = os.path.join(HERE, "..", "figures")
ATOM = os.path.join(HERE, "..", "..", "data", "atomic")
RUN = os.path.join(HERE, "..", "..", "examples", "hei_metastable",
                   "hei_high_density_heimeta.h5")
OUT = os.path.join(FIGDIR, "hei10830.pdf")

NPIX = 256


def main():
    hm = HeiMetaData(RUN, atomic_dir=ATOM)

    col, ext = hm.column_map("n_2s3", npix=NPIX)          # cm^-2
    tau, _ = hm.tau10830_map(component="doublet", npix=NPIX, warn=True)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))

    show_col = np.log10(np.where(col > 0, col, np.nan))
    im0 = axes[0].imshow(show_col.T, origin="lower", extent=ext,
                         cmap="inferno", aspect="equal")
    axes[0].set_title(r"$2^3$S column density")
    cb0 = fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
    cb0.set_label(r"$\log_{10} N(2^3{\rm S})\ [\mathrm{cm^{-2}}]$")

    show_tau = np.log10(np.where(tau > 0, tau, np.nan))
    im1 = axes[1].imshow(show_tau.T, origin="lower", extent=ext,
                         cmap="inferno", aspect="equal")
    axes[1].set_title(r"10830~\AA\ line-center $\tau_0$")
    cb1 = fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
    cb1.set_label(r"$\log_{10} \tau_0(\mathrm{10830\ doublet})$")
    # the tau_0 = 1 contour: inside it the optically-thin diagnostic breaks down
    cs = axes[1].contour(show_tau.T, levels=[0.0], colors="cyan",
                         linewidths=1.3, extent=ext)
    axes[1].clabel(cs, fmt={0.0: r"$\tau_0=1$"}, fontsize=_b * _S * 0.8)

    for ax in axes:
        ax.set_xlabel(r"$x$~[pc]")
        ax.set_ylabel(r"$y$~[pc]")
        square_ticks(ax)

    fig.tight_layout()
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)

    print("wrote", OUT)
    print("  max N(2^3S)      = %.3e cm^-2" % np.nanmax(col))
    print("  max tau_0(10830) = %.3e  (doublet; optically THICK core)"
          % np.nanmax(tau))
    print("  pixels with tau_0 > 1: %d of %d"
          % (int(np.count_nonzero(tau > 1.0)), tau.size))
    print("  A_31 (file) = %.4e s^-1;  10830 components (vac A, f, A s^-1):"
          % hm.a31)
    for row in hm.lines10830:
        print("    %.3f  %.4f  %.3e" % (row[0], row[1], row[2]))


if __name__ == "__main__":
    main()
