#!/usr/bin/env python3
"""Regenerate the dust-emission surface-brightness maps for the paper.

Reads the SEDust dust-band emissivity files (emis_dust [erg/s/cm^3] per band,
with wl_dust in Angstrom) and projects each band, flux-conservingly, to a
surface-brightness map S [erg/s/cm^2] with the same leaf deposit the paper's
line maps use.  Wavelengths are labeled in micron; all text is TeX.

Outputs (PDF, vector):
  dustemis_maps.pdf : two dust-emission bands (8 and 160 um) of the standard
                      run (PAHs everywhere).
  pahlive_maps.pdf  : the 8 um band, standard run vs the x_HII-weighted PAH
                      survival run -- the central peak turns into a front ring.

Data (existing outputs, no MoCHII re-run):
  ../../tests/d_dusty/dustemis_smoke_dustemis.h5   (PAHs everywhere)
  ../../tests/d_dusty/pahlive_smoke_dustemis.h5    (live PAH survival)
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"

# Axis, tick, and title fonts sized to ~85% of the body text once the figure
# is scaled down by \includegraphics (0.95\textwidth of a 10 pt document).
_S = 1.35
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "figure.titlesize": _b * 1.2 * _S,
})

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..",
                                "tools", "python"))
from mochii_output import read_sections, _LeafProjectable   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "figures")
DDIR = os.path.join(HERE, "..", "..", "tests", "d_dusty")
DUSTEMIS = os.path.join(DDIR, "dustemis_smoke_dustemis.h5")
PAHLIVE = os.path.join(DDIR, "pahlive_smoke_dustemis.h5")

SLABEL = r"$\log_{10} S\ [\mathrm{erg\,s^{-1}\,cm^{-2}}]$"
AXLABEL = r"(code units)"


class _DustEmis(_LeafProjectable):
    """Minimal projectable over the emis_dust blocks of a *_dustemis file."""

    def __init__(self, fname):
        sec = read_sections(fname)
        xyz = sec["LeafXYZ"]["data"]
        if xyz.shape[0] != 3:
            xyz = xyz.T
        self.xyz = np.asarray(xyz, float)             # (3, nleaf), code units
        self.size = np.asarray(sec["LeafSize"]["data"], float).ravel()
        self.dist_cm = float(sec["LeafXYZ"]["attrs"].get("DIST_CM", 1.0))
        self.emis = np.asarray(sec["emis_dust"]["data"], float)  # (nleaf, nband)
        if self.emis.shape[0] != self.size.size:
            self.emis = self.emis.T
        self.wl_um = np.asarray(sec["wl_dust"]["data"], float).ravel() / 1.0e4
        self.vol_cm3 = (self.size * self.dist_cm) ** 3

    def band_map(self, iband, npix=64):
        """Surface brightness S [erg/s/cm^2] of band iband, log10, + extent."""
        lum = self.emis[:, iband] * self.vol_cm3          # erg/s per leaf
        img, extent, apix = self._deposit(lum, "z", npix, None)
        S = img / apix                                    # erg/s/cm^2
        with np.errstate(divide="ignore"):
            logS = np.log10(np.where(S > 0, S, np.nan))
        return logS, extent


def _panel(ax, logS, extent, title, cmap="magma", vmin=None, vmax=None):
    im = ax.imshow(logS.T, origin="lower", extent=extent, cmap=cmap,
                   aspect="equal", vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlabel(AXLABEL)
    ax.set_ylabel(AXLABEL)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(SLABEL)


def make_dustemis():
    d = _DustEmis(DUSTEMIS)
    # bands nearest 8 and 160 um
    b8 = int(np.argmin(np.abs(d.wl_um - 8.0)))
    b160 = int(np.argmin(np.abs(d.wl_um - 160.0)))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, ib in zip(axes, (b8, b160)):
        logS, extent = d.band_map(ib)
        wl = d.wl_um[ib]
        title = r"dust $%g\,\mu$m" % wl
        _panel(ax, logS, extent, title)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "dustemis_maps.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out, "bands (um):", d.wl_um[b8], d.wl_um[b160])


def make_pahlive():
    ref = _DustEmis(DUSTEMIS)     # PAHs everywhere
    live = _DustEmis(PAHLIVE)     # x_HII-weighted PAH survival
    b8r = int(np.argmin(np.abs(ref.wl_um - 8.0)))
    b8l = int(np.argmin(np.abs(live.wl_um - 8.0)))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    logS_ref, extent = ref.band_map(b8r)
    logS_live, _ = live.band_map(b8l)
    # both panels are the same 8 um band, so they share one brightness scale:
    # the comparison is between the two morphologies, not between two scales.
    both = np.concatenate([logS_ref[np.isfinite(logS_ref)].ravel(),
                           logS_live[np.isfinite(logS_live)].ravel()])
    vmin, vmax = np.percentile(both, 1.0), both.max()
    _panel(axes[0], logS_ref, extent, r"$8\,\mu$m: PAHs everywhere",
           vmin=vmin, vmax=vmax)
    _panel(axes[1], logS_live, extent,
           r"$8\,\mu$m: PAHs destroyed in ionized gas",
           vmin=vmin, vmax=vmax)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "pahlive_maps.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    make_dustemis()
    make_pahlive()
