import sys
#!/usr/bin/env python3
"""Peel-off synthetic imaging of the dusty, FUV-illuminated Stromgren model.

Reads the peel-off image file of the dust-scattering + FUV run
(tests/peel/peel_scat.in): a level-6 octree Stromgren sphere with the FUV
band extension, MW dust absorption AND scattering in the band, and observer
peel-off toward the default +z observer.  The image carries direct and
dust-scattered surface-brightness [erg/s/cm^2/sr] maps in an ionizing and an
FUV channel; the pixel scale (DXIM [deg], DIST [pc]) sets the physical axes.

The figure is a three-panel row of the FUV channel, the channel that carries
the escaping light (in this converged, radiation-bounded H II region the
ionizing band is absorbed before it can scatter out, so the ionizing
scattered image is negligible -- physically correct for an embedded nebula):

  (a) direct FUV    -- the attenuated stellar image (a compact source);
  (b) scattered FUV -- the dust-scattered halo filling the dusty sphere;
  (c) total FUV     -- direct + scattered, the synthetic observation.

All three share one log surface-brightness scale, so the panels are directly
comparable: the scattered halo is fainter in surface brightness than the
direct source yet spatially extended, and their sum is what a telescope
would record.

Output: figures/peel.pdf   (vector; the imshow rasters live inside the PDF)

Data (existing peel output, no MoCHII re-run):
  ../../tests/peel/peel_scat_image.h5   direct/scatt x ion/fuv images
  ../../tests/peel/peel_scat_rates.h5   (only for the DIST_CM code-unit scale)

Run:  python3 make_peel.py
"""

import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
_HERE_AX = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE_AX)
from axis_style import square_ticks   # noqa: E402


plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"

# Axis, tick, and title fonts sized to ~85% of the body text once the wide
# three-panel row is scaled down by \includegraphics (\textwidth) in the paper.
_S = 1.6
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "figure.titlesize": _b * 1.2 * _S,
})

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "figures")
PEEL = os.path.join(HERE, "..", "..", "tests", "peel")
IMG = os.path.join(PEEL, "peel_scat_image.h5")
OUT = os.path.join(FIGDIR, "peel.pdf")

ZOOM_PC = 3.0                          # half-width of the plotted field [pc]
DYN = 4.0                              # log10 dynamic range shown
ILABEL = r"$\log_{10} I\ [\mathrm{erg\,s^{-1}\,cm^{-2}\,sr^{-1}}]$"


def read_image():
    """Return the four channel maps and the physical pixel scale [pc]."""
    with h5py.File(IMG, "r") as f:
        chan = {k: f[k + "/data"][:] for k in
                ("direct_ion", "scatt_ion", "direct_fuv", "scatt_fuv")}
        at = {k: float(np.asarray(v).ravel()[0])
              for k, v in f["direct_fuv"].attrs.items()
              if k in ("DXIM", "DIST")}
    npix = chan["direct_fuv"].shape[0]
    pc_per_pix = np.deg2rad(at["DXIM"]) * at["DIST"]      # DIST in code units = pc
    half = 0.5 * npix * pc_per_pix
    extent = (-half, half, -half, half)
    return chan, extent


def logmap(img):
    """log10 of a non-negative surface-brightness image; non-positive -> NaN."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log10(np.where(img > 0.0, img, np.nan))


def panel(ax, img, extent, title, cmap, vmin, vmax):
    # floor the empty (non-positive) pixels to the color-bar bottom so every
    # panel has the same black background -- the direct panel is then a point
    # source on black, not a blank white frame.
    shown = np.where(np.isfinite(img), img, vmin)
    im = ax.imshow(shown.T, origin="lower", extent=extent, cmap=cmap,
                   aspect="equal", vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlabel(r"$x$~[pc]")
    ax.set_ylabel(r"$y$~[pc]")
    ax.set_xlim(-ZOOM_PC, ZOOM_PC)
    ax.set_ylim(-ZOOM_PC, ZOOM_PC)
    square_ticks(ax)
    return im


def main():
    chan, extent = read_image()
    direct = chan["direct_fuv"]
    scatt = chan["scatt_fuv"]
    total = direct + scatt

    ld, ls, lt = logmap(direct), logmap(scatt), logmap(total)
    vmax = float(np.nanmax(lt))
    vmin = vmax - DYN

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.5))
    panel(axes[0], ld, extent, r"direct FUV (stellar)", "inferno", vmin, vmax)
    panel(axes[1], ls, extent, r"scattered FUV (dust halo)", "inferno",
          vmin, vmax)
    im = panel(axes[2], lt, extent, r"total FUV (direct $+$ scattered)",
               "inferno", vmin, vmax)

    # one shared colorbar to the right of the row
    fig.subplots_adjust(left=0.05, right=0.90, wspace=0.28)
    cax = fig.add_axes([0.915, 0.17, 0.013, 0.66])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(ILABEL)

    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)

    # report the numbers behind the panels
    sr_max = float(np.nanmax(ls))
    dr_max = float(np.nanmax(ld))
    scfrac = float(np.count_nonzero(scatt > 0) / scatt.size)
    print("wrote", OUT)
    print("  field of view (image) = +/- %.2f pc; plotted +/- %.1f pc"
          % (0.5 * (extent[1] - extent[0]), ZOOM_PC))
    print("  direct FUV peak  log10 I = %.2f  (compact stellar image)" % dr_max)
    print("  scatt  FUV peak  log10 I = %.2f  (halo covers %.1f%% of pixels)"
          % (sr_max, 100.0 * scfrac))
    print("  color scale: log10 I in [%.2f, %.2f]" % (vmin, vmax))
    # the ionizing scattered channel is negligible (band absorbed before escape)
    si = float(np.nanmax(chan["scatt_ion"]))
    print("  (ionizing scattered peak I = %.2e erg/s/cm^2/sr -- negligible: "
          "the band is absorbed before it can scatter out)" % si)


if __name__ == "__main__":
    main()
