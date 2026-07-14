"""Readers for MoCHII output files and 2D map making.

Covers the HDF5/FITS section files written by MoCHII (the '<base>_rates'
and '<base>_emis' files share the same layout: one section for each
EXTNAME, keywords as attributes) and the text outputs (_lines.txt,
_nebcont.txt, _dustir.txt, _dustsed.txt).

Notebook / script use — import and plot directly:

    import sys; sys.path.insert(0, ".../MoCHII/tools/python")
    from mochii_output import EmisData
    d = EmisData("smoke_emis_emis.h5")
    d.info()                               # blocks, fields, line count
    d.line_list()                          # every line: (block, wl[A])
    img, ext = d.line_map("o_3", 5007.0, axis="z", npix=512)
    img, ext = d.line_map("h", 4861.0)     # SH95 H I lines: block "h"
    img, ext = d.field_map("T_e", weight="EM")
    fig, ax = d.plot_line("o_3", 5007.0, log=True)      # ready-made figure
    fig, ax = d.plot_line("h", 6563.0, unit="intensity", photons=True)
    fig, ax = d.plot_field("T_e", weight="EM")

Maps are flux-conserving: each leaf's luminosity (emissivity x volume)
is deposited onto the pixel grid with exact area overlap, then divided
by the pixel area.  unit='flux' gives surface brightness [erg/s/cm^2];
unit='intensity' gives specific intensity [erg/s/cm^2/sr] (= flux/4pi,
isotropic emission, optically thin).  Field maps are line-of-sight
weighted averages (weight = 'EM' n_e^2, 'ne', 'nH', or 'V').

Command line:
    python3 mochii_output.py <base>_emis.h5 --line o_3 5007 --npix 512
    python3 mochii_output.py <base>_emis.h5 --field T_e --weight EM
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from leaf_field import leaf_slice, leaf_slice_nn, plot_slice, leaf_box_mask

HC_ERG_A = 1.9864458571489287e-08     # h*c [erg Angstrom]


# ---------------------------------------------------------------------------
# low-level: read every section of a MoCHII HDF5/FITS file
# ---------------------------------------------------------------------------
def read_sections(fname):
    """Return {section_name: {'data': ndarray, 'attrs': {key: value}}}.

    Works for both backends: HDF5 (groups '/<EXTNAME>/data' with the
    keywords as group attributes) and FITS (image HDUs with EXTNAME).
    Array axis order is as stored; use the nleaf axis to orient.
    """
    out = {}
    if fname.endswith((".h5", ".hdf5")):
        import h5py
        with h5py.File(fname, "r") as f:
            for name, grp in f.items():
                attrs = {}
                for k, v in grp.attrs.items():
                    a = np.asarray(v).ravel()
                    attrs[k] = a[0] if a.size == 1 else a
                data = grp["data"][:] if "data" in grp else None
                out[name] = {"data": data, "attrs": attrs}
    else:
        from astropy.io import fits
        with fits.open(fname) as hdul:
            for hdu in hdul:
                if hdu.data is None:
                    continue
                name = hdu.header.get("EXTNAME", f"hdu{hdul.index(hdu)}")
                out[name] = {"data": np.asarray(hdu.data),
                             "attrs": dict(hdu.header)}
    return out


def read_lines_table(fname):
    """Parse '<base>_lines.txt' -> list of dicts (elem, stage, wl, L, ratio)."""
    rows = []
    with open(fname) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            t = ln.split()
            rows.append(dict(elem=t[0], stage=int(t[1]), wl=float(t[2]),
                             L=float(t[3]), ratio=float(t[4]),
                             label=t[5] if len(t) > 5 else ""))
    return rows


def read_spectrum(fname):
    """2-column text spectra (_nebcont/_dustir/_dustsed): (x, y) arrays."""
    return np.loadtxt(fname, unpack=True)


def _section_names(fname):
    """List the section / EXTNAME names of a MoCHII file without loading the
    array data (a cheap way to tell an emissivity file from a rates file)."""
    if fname.endswith((".h5", ".hdf5")):
        import h5py
        with h5py.File(fname, "r") as f:
            return list(f.keys())
    from astropy.io import fits
    with fits.open(fname) as hdul:
        return [h.header.get("EXTNAME", "") for h in hdul]


# ---------------------------------------------------------------------------
# planar slices / cutouts over the leaf cells (shared by the file readers)
# ---------------------------------------------------------------------------
class _LeafSliceable:
    """Slice and cutout helpers over the cubic leaf cells.

    A subclass must set ``self.xyz`` (3, nleaf; code units) and ``self.fields``
    (name -> 1D leaf array), and may set ``self.size`` (the full cell width;
    None when unknown, which makes the slice fall back to a nearest-leaf
    lookup).  A slice is the cell value AT a plane (from leaf_field.leaf_slice),
    distinct from the line-of-sight projection maps.
    """

    # usetex-safe color-bar labels for the common scalar fields
    _CBAR = {
        "T_e":    r"$T_e$ [K]",
        "n_e":    r"$n_e$ [cm$^{-3}$]",
        "n_H":    r"$n_H$ [cm$^{-3}$]",
        "x_HI":   r"$x_{\rm HI}$",
        "T_dust": r"$T_{\rm dust}$ [K]",
    }

    @property
    def _half(self):
        """Leaf half-size array (0.5 * full width), or None when unknown."""
        size = getattr(self, "size", None)
        if size is None:
            return None
        return 0.5 * np.asarray(size, float).ravel()

    def _cbar_label(self, name):
        return self._CBAR.get(name, name.replace("_", r"\_"))

    def slice_field(self, name, axis="z", coord=0.0, npix=512, bounds=None):
        """True planar slice of a scalar leaf field: the cell value at the plane.

        Returns (img, extent).  When leaf half-sizes are available the exact
        containing-cell slice (leaf_slice) is used; otherwise a nearest-leaf
        slice (leaf_slice_nn) is used and a one-time note is printed.  Pass
        ``bounds`` = (umin, umax, vmin, vmax) with a ``coord`` inside the box
        to render a cutout.
        """
        value = np.asarray(self.fields[name], float)
        if value.ndim != 1:
            raise ValueError(f"'{name}' is not a scalar leaf field")
        half = self._half
        if half is not None:
            return leaf_slice(self.xyz[0], self.xyz[1], self.xyz[2], half,
                              value, axis=axis, coord=coord, bounds=bounds,
                              npix=npix)
        if not getattr(self, "_nn_noted", False):
            print("note: leaf half-sizes are absent; using a nearest-leaf "
                  "slice (leaf_slice_nn)")
            self._nn_noted = True
        return leaf_slice_nn(self.xyz[0], self.xyz[1], self.xyz[2], value,
                             axis=axis, coord=coord, bounds=bounds, npix=npix)

    def plot_field_slice(self, name, axis="z", coord=0.0, npix=512,
                         bounds=None, log=False, ax=None, cmap="inferno"):
        """Slice a scalar leaf field and draw it; returns (fig, ax)."""
        img, extent = self.slice_field(name, axis=axis, coord=coord,
                                       npix=npix, bounds=bounds)
        title = f"{name} slice {axis}={coord:g}".replace("_", r"\_")
        return plot_slice(img, extent, ax=ax, log=log, cmap=cmap,
                          cbar_label=self._cbar_label(name), title=title,
                          axis=axis)

    def cutout_mask(self, box):
        """Boolean leaf mask, True where a leaf overlaps the axis-aligned box
        (xmin, xmax, ymin, ymax, zmin, zmax).

        A cutout FOR PLOTTING is obtained instead by passing ``bounds`` to
        slice_field / plot_field_slice with a ``coord`` inside the box.
        """
        half = self._half
        h = half if half is not None else 0.0
        return leaf_box_mask(self.xyz[0], self.xyz[1], self.xyz[2], h, box)


# ---------------------------------------------------------------------------
# the emissivity file
# ---------------------------------------------------------------------------
class EmisData(_LeafSliceable):
    """Leaf emissivities + state from '<base>_emis.h5' (par%emis_output)."""

    def __init__(self, fname):
        sec = read_sections(fname)
        self.fname = fname

        xyz = sec["LeafXYZ"]["data"]
        if xyz.shape[0] != 3:
            xyz = xyz.T
        self.xyz = np.asarray(xyz, dtype=float)      # (3, nleaf), code units
        self.nleaf = self.xyz.shape[1]
        self.size = np.asarray(sec["LeafSize"]["data"], float).ravel()
        self.dist_cm = float(sec["LeafXYZ"]["attrs"].get("DIST_CM", 1.0))

        self.fields = {}      # state arrays (nleaf,)
        self.emis = {}        # block -> (nl, nleaf) [erg/s/cm^3]
        self.wl = {}          # block -> (nl,) [A]
        for name, s in sec.items():
            if name in ("LeafXYZ", "LeafSize"):
                continue
            data = s["data"]
            if name.startswith("wl_"):
                self.wl[name[3:]] = np.asarray(data, float).ravel()
            elif name.startswith("emis_"):
                em = np.asarray(data, float)
                if em.ndim == 1:
                    em = em[None, :]
                if em.shape[0] == self.nleaf:
                    em = em.T
                self.emis[name[5:]] = em
            else:
                arr = np.asarray(data, float)
                if arr.ndim == 1 and arr.size == self.nleaf:
                    self.fields[name] = arr
                else:                       # x_<el>_stages (nstage, nleaf)
                    if arr.shape[0] == self.nleaf:
                        arr = arr.T
                    self.fields[name] = arr
        self.vol_cm3 = (self.size * self.dist_cm)**3

    def info(self):
        """Print a summary (nice in a notebook)."""
        nlines = sum(len(w) for w in self.wl.values())
        print(f"{self.fname}: {self.nleaf} leaves, "
              f"{len(self.emis)} emissivity blocks, {nlines} lines")
        print("blocks:", ", ".join(sorted(self.emis)))
        print("fields:", ", ".join(sorted(self.fields)))

    def __repr__(self):
        return (f"EmisData('{self.fname}', nleaf={self.nleaf}, "
                f"blocks={len(self.emis)})")

    # -- line access --------------------------------------------------------
    def line_list(self):
        """Every stored line as (block, wavelength [A])."""
        return [(b, w) for b in sorted(self.emis) for w in self.wl[b]]

    def find_line(self, block, wl):
        """Index of the line nearest wl [A] in a block ('o_3', 'h', ...)."""
        if block not in self.emis:
            raise KeyError(f"no emissivity block '{block}' "
                           f"(available: {sorted(self.emis)})")
        k = int(np.argmin(np.abs(self.wl[block] - wl)))
        if abs(self.wl[block][k] - wl) > 0.01 * wl + 2.0:
            raise ValueError(f"nearest line in '{block}' is "
                             f"{self.wl[block][k]:.1f} A, far from {wl}")
        return k

    def line_luminosity(self, block, wl, photons=False):
        """L_line [erg/s], or the photon rate [photons/s] with photons=True."""
        k = self.find_line(block, wl)
        L = float((self.emis[block][k] * self.vol_cm3).sum())
        if photons:
            L /= HC_ERG_A / self.wl[block][k]
        return L

    # -- projection ---------------------------------------------------------
    def _axes(self, axis):
        ia = "xyz".index(axis)
        iu, iv = [i for i in range(3) if i != ia]
        return ia, iu, iv

    def _grid(self, iu, iv, npix):
        h = 0.5 * self.size
        lo = min((self.xyz[iu] - h).min(), (self.xyz[iv] - h).min())
        hi = max((self.xyz[iu] + h).max(), (self.xyz[iv] + h).max())
        edges = np.linspace(lo, hi, npix + 1)
        return edges, (lo, hi, lo, hi)

    def _deposit(self, lum, axis, npix, slab):
        """Flux-conserving deposit of leaf luminosities [erg/s] onto the
        image grid; returns (sum image [erg/s], extent, pixel area cm^2)."""
        ia, iu, iv = self._axes(axis)
        edges, extent = self._grid(iu, iv, npix)
        dpix = edges[1] - edges[0]
        img = np.zeros((npix, npix))

        u = self.xyz[iu];  v = self.xyz[iv];  h = 0.5 * self.size
        sel = lum != 0.0
        if slab is not None:
            w = self.xyz[ia]
            sel &= (w + h > slab[0]) & (w - h < slab[1])
        idx = np.nonzero(sel)[0]
        lo = edges[0]
        for il in idx:
            L = lum[il]
            if slab is not None:      # partial slab overlap along the LOS
                ov = (min(self.xyz[ia][il] + h[il], slab[1])
                      - max(self.xyz[ia][il] - h[il], slab[0]))
                L = L * ov / (2.0 * h[il])
            i0 = max(int((u[il] - h[il] - lo) / dpix), 0)
            i1 = min(int((u[il] + h[il] - lo) / dpix) + 1, npix)
            j0 = max(int((v[il] - h[il] - lo) / dpix), 0)
            j1 = min(int((v[il] + h[il] - lo) / dpix) + 1, npix)
            if i1 <= i0 or j1 <= j0:
                continue
            # exact 1D overlaps of the leaf footprint with each pixel
            eu = edges[i0:i1 + 1]
            ou = (np.minimum(eu[1:], u[il] + h[il])
                  - np.maximum(eu[:-1], u[il] - h[il])).clip(min=0.0)
            ev = edges[j0:j1 + 1]
            ov2 = (np.minimum(ev[1:], v[il] + h[il])
                   - np.maximum(ev[:-1], v[il] - h[il])).clip(min=0.0)
            frac = np.outer(ou, ov2) / (2.0 * h[il])**2
            img[i0:i1, j0:j1] += L * frac
        apix_cm2 = (dpix * self.dist_cm)**2
        return img, extent, apix_cm2

    def line_map(self, block, wl, axis="z", npix=256, slab=None,
                 unit="flux", photons=False):
        """Projected map of one line.

        block : emissivity block ('h' for the SH95 H I lines, else
                '<el>_<stage>', e.g. 'o_3'); wl : wavelength [A] (nearest
        line is used); axis : projection axis; slab : (lo, hi) cut along
        the line of sight in code units.
        unit : 'flux'      -> surface brightness [erg/s/cm^2]
                              (sums to L_line when multiplied by the
                              pixel area);
               'intensity' -> specific intensity [erg/s/cm^2/sr]
                              (= flux / 4 pi; the emission is isotropic
                              and the nebula is optically thin to it).
        photons : divide by the photon energy hc/lambda of THIS line ->
                  [photons/s/cm^2] or [photons/s/cm^2/sr].
        Returns (image, extent).
        """
        k = self.find_line(block, wl)
        lum = self.emis[block][k] * self.vol_cm3
        img, extent, apix = self._deposit(lum, axis, npix, slab)
        img = img / apix
        if unit == "intensity":
            img = img / (4.0 * np.pi)
        elif unit != "flux":
            raise ValueError("unit must be 'flux' or 'intensity'")
        if photons:
            img = img / (HC_ERG_A / self.wl[block][k])
        return img, extent

    def field_map(self, field, axis="z", npix=256, weight="EM", slab=None):
        """Line-of-sight weighted-average map of a state field.

        field : 'T_e', 'n_e', 'x_HI', ... (see .fields); weight : 'EM'
        (n_e^2 V), 'ne' (n_e V), 'nH' (n_H V), or 'V'.  Returns
        (image, extent); empty pixels are NaN.
        """
        f = self.fields[field]
        if f.ndim != 1:
            raise ValueError(f"'{field}' is not a scalar field")
        if weight == "EM":
            w = self.fields["n_e"]**2 * self.vol_cm3
        elif weight == "ne":
            w = self.fields["n_e"] * self.vol_cm3
        elif weight == "nH":
            w = self.fields["n_H"] * self.vol_cm3
        elif weight == "V":
            w = self.vol_cm3.copy()
        else:
            raise ValueError("weight must be 'EM', 'ne', 'nH', or 'V'")
        num, extent, _ = self._deposit(f * w, axis, npix, slab)
        den, _, _ = self._deposit(w, axis, npix, slab)
        with np.errstate(invalid="ignore", divide="ignore"):
            img = num / den
        img[den <= 0.0] = np.nan
        return img, extent


    # -- ready-made figures (notebook use) -----------------------------------
    def _show(self, img, extent, title, cbar_label, log, ax):
        import matplotlib.pyplot as plt
        show = img
        if log:
            with np.errstate(divide="ignore", invalid="ignore"):
                show = np.log10(np.where(img > 0, img, np.nan))
            cbar_label = r"$\log_{10}$ " + cbar_label
        if ax is None:
            fig, ax = plt.subplots(figsize=(5.4, 4.4))
        else:
            fig = ax.figure
        im = ax.imshow(show.T, origin="lower", extent=extent, cmap="inferno")
        fig.colorbar(im, ax=ax, label=cbar_label)
        ax.set_xlabel("(code units)")
        ax.set_ylabel("(code units)")
        ax.set_title(title)
        return fig, ax

    def plot_line(self, block, wl, axis="z", npix=256, slab=None,
                  unit="flux", photons=False, log=True, ax=None):
        """Line map as a ready-made matplotlib figure; returns (fig, ax).

        Same options as line_map plus log (log10 color scale, default on)
        and ax (draw into an existing Axes, e.g. a subplot grid).
        """
        img, ext = self.line_map(block, wl, axis=axis, npix=npix,
                                 slab=slab, unit=unit, photons=photons)
        k = self.find_line(block, wl)
        title = block.replace("_", " ") + f" {self.wl[block][k]:.0f}" + r" \AA"
        eu = "photons" if photons else "erg"
        if unit == "intensity":
            lbl = r"$I$ [" + eu + r" s$^{-1}$ cm$^{-2}$ sr$^{-1}$]"
        else:
            lbl = r"$S$ [" + eu + r" s$^{-1}$ cm$^{-2}$]"
        return self._show(img, ext, title, lbl, log, ax)

    def plot_field(self, field, axis="z", npix=256, weight="EM",
                   slab=None, log=False, ax=None):
        """Field map as a ready-made matplotlib figure; returns (fig, ax)."""
        img, ext = self.field_map(field, axis=axis, npix=npix,
                                  weight=weight, slab=slab)
        lbl = field.replace("_", r"\_")
        title = lbl + f" ({weight}-weighted)"
        return self._show(img, ext, title, lbl, log, ax)

    # -- planar slices (value at a plane, not a projection) -----------------
    def slice_line(self, block, wl, axis="z", coord=0.0, npix=512,
                   bounds=None):
        """True planar slice of one line's leaf emissivity [erg/s/cm^3].

        value = the leaf emissivity of the line nearest ``wl`` [A] in
        ``block`` ('o_3', 'h', ...); returns (img, extent).  The emissivity
        file always carries leaf sizes, so the exact containing-cell slice is
        used.
        """
        k = self.find_line(block, wl)
        value = self.emis[block][k]
        return leaf_slice(self.xyz[0], self.xyz[1], self.xyz[2], self._half,
                          value, axis=axis, coord=coord, bounds=bounds,
                          npix=npix)

    def plot_line_slice(self, block, wl, axis="z", coord=0.0, npix=512,
                        bounds=None, log=True, ax=None, cmap="inferno"):
        """Slice one line's emissivity and draw it; returns (fig, ax)."""
        img, ext = self.slice_line(block, wl, axis=axis, coord=coord,
                                    npix=npix, bounds=bounds)
        k = self.find_line(block, wl)
        title = (block.replace("_", " ") + f" {self.wl[block][k]:.0f}"
                 + r" \AA{} slice " + f"{axis}={coord:g}")
        cbar = r"$\epsilon$ [erg s$^{-1}$ cm$^{-3}$]"
        return plot_slice(img, ext, ax=ax, log=log, cmap=cmap,
                          cbar_label=cbar, title=title, axis=axis)


# ---------------------------------------------------------------------------
# the rates / state file
# ---------------------------------------------------------------------------
class RatesData(_LeafSliceable):
    """Leaf state + rates from '<base>_rates.h5' (or any section file that
    carries LeafXYZ plus scalar leaf fields).

    Inherits the slice / cutout helpers of _LeafSliceable.  ``self.size`` is
    the 'LeafSize' block (full cell width) when the file has one, else None --
    in which case the slice falls back to a nearest-leaf lookup.
    """

    # sections that are not scalar leaf fields
    _SKIP = ("LeafXYZ", "LeafSize", "E_bin", "dE_bin")

    def __init__(self, fname):
        sec = read_sections(fname)
        self.fname = fname

        xyz = sec["LeafXYZ"]["data"]
        if xyz.shape[0] != 3:
            xyz = xyz.T
        self.xyz = np.asarray(xyz, dtype=float)      # (3, nleaf), code units
        self.nleaf = self.xyz.shape[1]
        self.dist_cm = float(sec["LeafXYZ"]["attrs"].get("DIST_CM", 1.0))

        if "LeafSize" in sec and sec["LeafSize"]["data"] is not None:
            self.size = np.asarray(sec["LeafSize"]["data"], float).ravel()
        else:
            self.size = None

        self.fields = {}      # name -> 1D leaf field
        for name, s in sec.items():
            if name in self._SKIP or name.startswith("NNU_"):
                continue
            data = s["data"]
            if data is None:
                continue
            arr = np.asarray(data, float)
            if arr.ndim == 1 and arr.size == self.nleaf:
                self.fields[name] = arr
            # 2D blocks (J_nu, x_<el>_stages) are not scalar fields -> skipped

    def info(self):
        """Print a summary (nice in a notebook)."""
        has = "yes" if self.size is not None else "no"
        print(f"{self.fname}: {self.nleaf} leaves, LeafSize present: {has}")
        print("fields:", ", ".join(sorted(self.fields)))

    def __repr__(self):
        return (f"RatesData('{self.fname}', nleaf={self.nleaf}, "
                f"leafsize={self.size is not None}, "
                f"fields={len(self.fields)})")


# ---------------------------------------------------------------------------
# command line: quick-look maps
# ---------------------------------------------------------------------------
def _main():
    import argparse
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ap = argparse.ArgumentParser(
        description="MoCHII quick-look maps (projection) and slices")
    ap.add_argument("infile", help="<base>_emis.h5 or <base>_rates.h5")
    ap.add_argument("--line", nargs=2, metavar=("BLOCK", "WL"),
                    help="projected line map: block ('o_3', 'h', ...) + "
                         "wavelength [A] (emis file)")
    ap.add_argument("--unit", default="flux", choices=["flux", "intensity"],
                    help="line map unit: flux [erg/s/cm^2] or intensity "
                         "[erg/s/cm^2/sr] (default flux)")
    ap.add_argument("--photons", action="store_true",
                    help="photon units instead of erg (divide by hc/lambda)")
    ap.add_argument("--field", help="projected field map: T_e, n_e, x_HI, ...")
    ap.add_argument("--weight", default="EM",
                    help="field-map weight: EM, ne, nH, V (default EM)")
    ap.add_argument("--slice", metavar="FIELD",
                    help="planar slice of a scalar leaf field "
                         "(T_e, n_e, x_HI, ...) at --axis = --coord")
    ap.add_argument("--coord", type=float, default=0.0,
                    help="plane position along --axis [code units] "
                         "(default 0)")
    ap.add_argument("--axis", default="z", choices=["x", "y", "z"])
    ap.add_argument("--npix", type=int, default=256)
    ap.add_argument("--log", action="store_true", help="log10 color scale")
    ap.add_argument("--out", default=None, help="output PNG")
    args = ap.parse_args()

    is_emis = any(n.startswith("emis_") for n in _section_names(args.infile))

    if args.slice:
        d = EmisData(args.infile) if is_emis else RatesData(args.infile)
        fig, ax = d.plot_field_slice(args.slice, axis=args.axis,
                                     coord=args.coord, npix=args.npix,
                                     log=args.log)
        base = f"slice_{args.slice}_{args.axis}{args.coord:g}"
    elif args.line:
        d = EmisData(args.infile)
        block, wl = args.line[0], float(args.line[1])
        fig, ax = d.plot_line(block, wl, axis=args.axis, npix=args.npix,
                              unit=args.unit, photons=args.photons,
                              log=args.log)
        k = d.find_line(block, wl)
        base = f"map_{block}_{int(round(d.wl[block][k]))}"
    elif args.field:
        d = EmisData(args.infile)
        fig, ax = d.plot_field(args.field, axis=args.axis, npix=args.npix,
                               weight=args.weight, log=args.log)
        base = f"map_{args.field}_{args.weight}"
    else:
        d = EmisData(args.infile) if is_emis else RatesData(args.infile)
        d.info()
        if is_emis:
            print("stored lines:")
            for b, w in d.line_list():
                print(f"  {b:8s} {w:12.2f} A")
        return

    fig.tight_layout()
    out = args.out or base + ".png"
    fig.savefig(out, dpi=200)
    print("written:", out)


if __name__ == "__main__":
    _main()
