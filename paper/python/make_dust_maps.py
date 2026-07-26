#!/usr/bin/env python3
"""Regenerate the dust-emission surface-brightness maps for the paper.

Reads the SEDust dust-band emissivity files (emis_dust [erg/s/cm^3/um] per
band, with wl_dust in Angstrom) and projects each band, flux-conservingly,
with the same leaf deposit the paper's line maps use.

Units: emis_dust is the power radiated per unit volume and per unit
wavelength into the full 4 pi sr.  The deposit therefore gives L_lambda per
pixel; dividing by the pixel area and by 4 pi turns it into the specific
intensity I_lambda [erg/s/cm^2/sr/um], which is what the panels show (the
nebula is optically thin in the infrared and the emission is isotropic, so
I_lambda is just the emissivity integrated along the line of sight).
Wavelengths are labeled in micron; all text is TeX.

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
_HERE_AX = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE_AX)
from axis_style import square_ticks   # noqa: E402


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

SLABEL = (r"$\log_{10} I_\lambda\ "
          r"[\mathrm{erg\,s^{-1}\,cm^{-2}\,sr^{-1}}\,\mu\mathrm{m}^{-1}]$")
# LeafXYZ is in code length units; these runs set par%distance_unit = 'pc',
# so DIST_CM = 3.0857e18 and one code unit is exactly 1 pc.
XLABEL = r"$x$ [pc]"
YLABEL = r"$y$ [pc]"
RLABEL = r"$r$ [pc]"


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
        """Specific intensity I_lambda of band iband, log10, + extent.

        [erg/s/cm^2/sr/um]: the deposited L_lambda per pixel divided by the
        pixel area and by the 4 pi sr the emissivity is integrated over.
        """
        lum = self.emis[:, iband] * self.vol_cm3          # erg/s/um per leaf
        img, extent, apix = self._deposit(lum, "z", npix, None)
        S = img / apix / (4.0 * np.pi)                    # erg/s/cm^2/sr/um
        with np.errstate(divide="ignore"):
            logS = np.log10(np.where(S > 0, S, np.nan))
        return logS, extent


def _panel(ax, logS, extent, title, cmap="magma", vmin=None, vmax=None):
    im = ax.imshow(logS.T, origin="lower", extent=extent, cmap=cmap,
                   aspect="equal", vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    square_ticks(ax)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(SLABEL)


def radial_profile(logS, extent, nbin=28):
    """Azimuthally averaged surface brightness against projected radius.

    The average is taken in linear S over the pixels of each annulus (the
    physically meaningful mean), and returned as log10 S at the annulus
    centers, so the profile is the radial cut of the map above.
    """
    ny, nx = logS.shape
    x = np.linspace(extent[0], extent[1], nx)
    y = np.linspace(extent[2], extent[3], ny)
    # logS is indexed [ix, iy] (the maps are drawn with logS.T)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    r = np.sqrt(xx ** 2 + yy ** 2).ravel()
    S = (10.0 ** logS).ravel()
    good = np.isfinite(S)
    r, S = r[good], S[good]
    edges = np.linspace(0.0, r.max(), nbin + 1)
    idx = np.digitize(r, edges) - 1
    rc, Sc = [], []
    for i in range(nbin):
        m = idx == i
        if m.sum():
            rc.append(0.5 * (edges[i] + edges[i + 1]))
            Sc.append(S[m].mean())
    rc, Sc = np.asarray(rc), np.asarray(Sc)
    with np.errstate(divide="ignore"):
        return rc, np.log10(np.where(Sc > 0, Sc, np.nan))


def pah_scale(ref, live):
    """The common 8 um brightness scale of the PAH-survival comparison.

    Returned so the 8 um panel of the two-band figure carries the same
    scale, making the two figures directly comparable.
    """
    b8r = int(np.argmin(np.abs(ref.wl_um - 8.0)))
    b8l = int(np.argmin(np.abs(live.wl_um - 8.0)))
    a, _ = ref.band_map(b8r)
    b, _ = live.band_map(b8l)
    both = np.concatenate([a[np.isfinite(a)].ravel(),
                           b[np.isfinite(b)].ravel()])
    return np.percentile(both, 1.0), both.max()


def make_dustemis(d, vlim8=None):
    # bands nearest 8 and 160 um
    b8 = int(np.argmin(np.abs(d.wl_um - 8.0)))
    b160 = int(np.argmin(np.abs(d.wl_um - 160.0)))
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    prof = []
    for ax, ib in zip(axes[:2], (b8, b160)):
        logS, extent = d.band_map(ib)
        wl = d.wl_um[ib]
        title = r"dust $%g\,\mu$m" % wl
        vmin, vmax = vlim8 if (ib == b8 and vlim8 is not None) else (None, None)
        _panel(ax, logS, extent, title, vmin=vmin, vmax=vmax)
        prof.append((radial_profile(logS, extent), wl))
    ax = axes[2]
    for (rc, lS), wl in prof:
        ax.plot(rc, lS, lw=1.6, label=r"$%g\,\mu$m" % wl)
    ax.set_xlabel(RLABEL)
    ax.set_ylabel(SLABEL)
    ax.set_title("azimuthally averaged profile")
    ax.legend(frameon=False)
    ax.grid(alpha=0.3, lw=0.4)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "dustemis_maps.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out, "bands (um):", d.wl_um[b8], d.wl_um[b160])


def make_pahlive(ref, live, vlim8):
    b8r = int(np.argmin(np.abs(ref.wl_um - 8.0)))
    b8l = int(np.argmin(np.abs(live.wl_um - 8.0)))
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    logS_ref, extent = ref.band_map(b8r)
    logS_live, _ = live.band_map(b8l)
    # both panels are the same 8 um band, so they share one brightness scale:
    # the comparison is between the two morphologies, not between two scales.
    vmin, vmax = vlim8
    _panel(axes[0], logS_ref, extent, r"$8\,\mu$m: PAHs everywhere",
           vmin=vmin, vmax=vmax)
    _panel(axes[1], logS_live, extent,
           r"$8\,\mu$m: PAHs destroyed in ionized gas",
           vmin=vmin, vmax=vmax)
    ax = axes[2]
    for logS, lab in ((logS_ref, "PAHs everywhere"),
                      (logS_live, "PAHs destroyed in ionized gas")):
        rc, lS = radial_profile(logS, extent)
        ax.plot(rc, lS, lw=1.6, label=lab)
    ax.set_xlabel(RLABEL)
    ax.set_ylabel(SLABEL)
    ax.set_title(r"azimuthally averaged profile ($8\,\mu$m)")
    ax.legend(frameon=False, fontsize="small")
    ax.grid(alpha=0.3, lw=0.4)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "pahlive_maps.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    ref = _DustEmis(DUSTEMIS)     # PAHs everywhere
    live = _DustEmis(PAHLIVE)     # x_HII-weighted PAH survival
    vlim8 = pah_scale(ref, live)
    make_dustemis(ref, vlim8)
    make_pahlive(ref, live, vlim8)
