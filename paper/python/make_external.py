#!/usr/bin/env python3
"""A cloud photoionized by an EXTERNAL isotropic radiation field.

The composable-sources example examples/external_field/external_field.in has no
internal star: an isotropic background of band-integrated mean intensity J eats
into an n_H = 1.7 cm^-3, R = 1.8 pc cloud from OUTSIDE, entering the box faces
(par%ext_geometry = 'rec') with a 4e4 K blackbody shape.  The one-sided surface
flux pi*J keeps a skin ionized against recombination while the core stays
neutral -- the opposite topology of a central-source Stromgren sphere.

The figure has two panels:

  (a) a planar slice of x_HI through the cloud center (z = 0): an ionized skin
      (x_HI small) wrapping a neutral core (x_HI -> 1);
  (b) shell-averaged radial profiles of x_HI and x_HII, showing x_HII -> 1 at
      the illuminated surface and the neutral core at small radius.

Output: figures/external_field.pdf   (vector; the slice raster lives in the PDF)

Data (existing example output, no MoCHII re-run):
  ../../examples/external_field/external_field_rates.h5

Run:  python3 make_external.py
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
# figure is scaled down by \includegraphics (0.95\textwidth) in the paper.
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
from mochii_output import RatesData             # noqa: E402

FIGDIR = os.path.join(HERE, "..", "figures")
RUN = os.path.join(HERE, "..", "..", "examples", "external_field",
                   "external_field_rates.h5")
OUT = os.path.join(FIGDIR, "external_field.pdf")

RCLOUD = 1.8                     # par%rmax [pc]; nH = 0 outside
NPIX = 400


def main():
    d = RatesData(RUN)
    xyz = d.xyz                                   # (3, nleaf), code units = pc
    r = np.sqrt((xyz ** 2).sum(axis=0))
    x_HI = d.fields["x_HI"]
    x_HII = 1.0 - x_HI                            # hydrogen: only H0 and H+
    incloud = r <= RCLOUD

    # --- panel (a): x_HI slice through the center, cloud only -------------
    # slice a copy of x_HI masked to the cloud so the vacuum outside is blank.
    d.fields["x_HI_cloud"] = np.where(incloud, x_HI, np.nan)
    img, extent = d.slice_field("x_HI_cloud", axis="z", coord=0.0, npix=NPIX)

    fig, (axm, axp) = plt.subplots(1, 2, figsize=(11.0, 4.6))

    im = axm.imshow(np.log10(np.where(img > 0, img, np.nan)).T, origin="lower",
                    extent=extent, cmap="viridis", aspect="equal",
                    vmin=-4.0, vmax=0.0)
    axm.set_xlabel(r"$x$~[pc]")
    axm.set_ylabel(r"$y$~[pc]")
    axm.set_title(r"$x_{\rm HI}$ slice ($z=0$): ionized skin, neutral core")
    cb = fig.colorbar(im, ax=axm, fraction=0.046, pad=0.04)
    cb.set_label(r"$\log_{10} x_{\rm HI}$")

    # --- panel (b): shell-averaged radial profiles -----------------------
    nbin = 40
    edges = np.linspace(0.0, RCLOUD, nbin + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])

    def shell(q):
        out = np.full(nbin, np.nan)
        for i in range(nbin):
            m = incloud & (r >= edges[i]) & (r < edges[i + 1])
            if m.any():
                out[i] = np.median(q[m])
        return out

    axp.plot(mid, shell(x_HII), "-", color="#d62728", lw=2.0,
             label=r'$x_{\rm HII}$ (ionized)')
    axp.plot(mid, shell(x_HI), "--", color="#1f77b4", lw=2.0,
             label=r'$x_{\rm HI}$ (neutral)')
    axp.set_yscale("log")
    square_ticks(axm)
    axp.set_ylim(1e-4, 1.5)
    axp.set_xlim(0.0, RCLOUD)
    axp.set_xlabel(r"$r$~[pc]")
    axp.set_ylabel(r"ionization fraction")
    axp.set_title(r"radial structure (field enters from $r=R$)")
    axp.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
    axp.legend(loc="center left", frameon=False, fontsize=_b * _S * 0.9)

    fig.tight_layout()
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)

    # report the topology numbers
    core = incloud & (r < 0.5)
    skin = incloud & (r > RCLOUD - 0.3)
    print("wrote", OUT)
    print("  cloud R = %.1f pc, nleaf(in cloud) = %d" % (RCLOUD, incloud.sum()))
    print("  core   (r<0.5 pc):  median x_HI  = %.3f  (neutral)"
          % np.median(x_HI[core]))
    print("  skin (r>%.1f pc):    median x_HII = %.3f  (ionized)"
          % (RCLOUD - 0.3, np.median(x_HII[skin])))
    print("  -> external field ionizes the surface and leaves a neutral core "
          "(opposite of a central source)")


if __name__ == "__main__":
    main()
