#!/usr/bin/env python3
"""He I 2^3S metastable / 10833 A diagnostic profiles for the paper.

Reads the high-density He I metastable example
(examples/hei_metastable/hei_high_density.in: an n_H = 1e4 cm^-3, R = 0.15 pc
sphere at a fixed T_e = 1e4 K, par%hei_metastable + par%emis_output) and shows
the two sides of the same 2^3S - 2^3P transition as radial profiles:

  (a) absorption: the 10833 A line-center optical depth tau_0 of the blended
      J'=1,2 doublet against projected radius (impact parameter), with the
      2^3S column density N(2^3S) on the right axis.  Both are line-of-sight
      integrals, so the projected radius is their natural abscissa; at the
      fixed T_e of this run the Doppler width is constant and the two are
      proportional, so the right axis is the conversion of the left one.
      The tau_0 = 1 level is marked: inside it the optically-thin diagnostic
      (its own assumption) breaks down, which the reader also flags at run
      time.
  (b) emission: the 10833 A emissivity against (spherical) radius, a local
      quantity, taken from the code's own Porter et al. (2012/2013) He I block
      in the emissivity file.  Note the Porter emissivity is a
      collisional-radiative total, so it is read as written and nothing is
      summed onto it.

The opacity uses the fine-structure constants (vacuum wavelengths,
oscillator strength, Einstein A) from data/atomic/hei_metastable.txt, the same
file the Fortran diagnostic reads.

Output: figures/hei10833.pdf   (vector)

Data (existing example output, no MoCHII re-run):
  ../../examples/hei_metastable/hei_high_density_heimeta.h5
  ../../examples/hei_metastable/hei_high_density_emis.h5

Run:  python3 make_hei10833.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"

# Axis, tick, and title fonts sized to ~85% of the body text once the two-panel
# row is scaled down by \includegraphics (0.95\textwidth) in the paper.
_S = 1.35
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "legend.fontsize":  _b * _S * 0.85,
    "figure.titlesize": _b * 1.2 * _S,
})

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "python"))
from mochii_output import HeiMetaData, EmisData      # noqa: E402

FIGDIR = os.path.join(HERE, "..", "figures")
ATOM = os.path.join(HERE, "..", "..", "data", "atomic")
RUNDIR = os.path.join(HERE, "..", "..", "examples", "hei_metastable")
META = os.path.join(RUNDIR, "hei_high_density_heimeta.h5")
EMIS = os.path.join(RUNDIR, "hei_high_density_emis.h5")
OUT = os.path.join(FIGDIR, "hei10833.pdf")

NPIX = 256
RMAX_PC = 0.16             # box half-width of the run (sphere R = 0.15 pc)
WL10833 = 10833.22         # He I triplet line, vacuum [A] (10830.25 air)


def azimuthal_profile(img, extent, nbin=40):
    """Mean of a map over annuli of projected radius (NaN pixels dropped)."""
    nx, ny = img.shape
    x = np.linspace(extent[0], extent[1], nx)
    y = np.linspace(extent[2], extent[3], ny)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    r = np.sqrt(xx ** 2 + yy ** 2).ravel()
    v = np.asarray(img, float).ravel()
    good = np.isfinite(v) & (v > 0.0)
    r, v = r[good], v[good]
    edges = np.linspace(0.0, r.max(), nbin + 1)
    idx = np.digitize(r, edges) - 1
    rc, vc = [], []
    for i in range(nbin):
        m = idx == i
        if m.sum():
            rc.append(0.5 * (edges[i] + edges[i + 1]))
            vc.append(v[m].mean())
    return np.asarray(rc), np.asarray(vc)


def radial_emissivity(em, xyz, nbin=40):
    """Volume-weighted mean emissivity in spherical shells [erg/s/cm^3]."""
    r = np.sqrt((xyz ** 2).sum(axis=0))
    edges = np.linspace(0.0, r.max(), nbin + 1)
    idx = np.digitize(r, edges) - 1
    rc, ec = [], []
    for i in range(nbin):
        m = idx == i
        if m.sum():
            rc.append(0.5 * (edges[i] + edges[i + 1]))
            ec.append(em[m].mean())
    return np.asarray(rc), np.asarray(ec)


def main():
    hm = HeiMetaData(META, atomic_dir=ATOM)
    col, ext = hm.column_map("n_2s3", npix=NPIX)              # cm^-2
    tau, _ = hm.tau10833_map(component="doublet", npix=NPIX, warn=True)

    b_tau, p_tau = azimuthal_profile(tau, ext)
    b_col, p_col = azimuthal_profile(col, ext)

    ed = EmisData(EMIS)
    k = ed.find_line("hei", WL10833)
    wl_used = ed.wl["hei"][k]
    r_em, p_em = radial_emissivity(ed.emis["hei"][k], ed.xyz)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))

    #--- (a) absorption: tau_0 (left) and N(2^3S) (right)
    ax = axes[0]
    col_tau, col_col = "#1f77b4", "#d62728"
    ax.plot(b_tau, p_tau, color=col_tau, lw=1.8)
    ax.set_yscale("log")
    ax.set_xlabel(r"projected radius $b$~[pc]")
    ax.set_ylabel(r"$\tau_0(\mathrm{10833\ doublet})$", color=col_tau)
    ax.tick_params(axis="y", colors=col_tau)
    ax.axhline(1.0, color="0.4", ls=":", lw=1.2)
    ax.text(0.02, 1.35, r"$\tau_0 = 1$", color="0.35",
            transform=ax.get_yaxis_transform(), fontsize=_b * _S * 0.85)
    ax.set_title(r"absorption: 10833~\AA\ opacity")

    axr = ax.twinx()
    axr.plot(b_col, p_col, color=col_col, lw=1.8, ls="--")
    axr.set_yscale("log")
    axr.set_ylabel(r"$N(2^3{\rm S})$~[cm$^{-2}$]", color=col_col)
    axr.tick_params(axis="y", colors=col_col)
    # the two are proportional at fixed T_e; lock the axes so the curves
    # overlie and the right axis reads as the conversion of the left one
    ratio = np.median(p_col) / np.median(p_tau)
    axr.set_ylim(*(np.asarray(ax.get_ylim()) * ratio))

    #--- (b) emission: the Porter He I 10833 emissivity
    ax = axes[1]
    ax.plot(r_em, p_em, color="#2ca02c", lw=1.8)
    ax.set_yscale("log")
    ax.set_xlabel(r"radius $r$~[pc]")
    ax.set_ylabel(r"$j(10833)$~[erg\,s$^{-1}$\,cm$^{-3}$]")
    ax.set_title(r"emission: 10833~\AA\ emissivity")
    # the exterior is empty by construction (n_H = 0 outside R = 0.15 pc), so
    # show the decades that carry the He+ zone rather than the numerical floor
    peak = np.nanmax(p_em)
    ax.set_ylim(peak*1.0e-5, peak*3.0)

    for ax in axes:
        ax.grid(alpha=0.3, lw=0.4)
        ax.set_xlim(0.0, RMAX_PC)

    fig.tight_layout()
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)

    print("wrote", OUT)
    print("  max N(2^3S)      = %.3e cm^-2" % np.nanmax(col))
    print("  max tau_0(10833) = %.3e  (doublet; optically THICK core)"
          % np.nanmax(tau))
    print("  central tau_0 (profile) = %.3e at b = %.4f pc" % (p_tau[0], b_tau[0]))
    print("  tau_0 = 1 crossed at b = %.4f pc"
          % np.interp(0.0, np.log10(p_tau[::-1]), b_tau[::-1]))
    print("  He I line used: %.2f A (nearest %.1f)" % (wl_used, WL10833))
    print("  j(10833): centre %.3e, edge %.3e erg/s/cm^3" % (p_em[0], p_em[-1]))
    print("  L(10833) = %.4e erg/s" % ed.line_luminosity("hei", WL10833))


if __name__ == "__main__":
    main()
