"""Convert a MoCHII He I 2^3S metastable diagnostic file to a LaRT v2.00 grid.

MoCHII writes '<base>_heimeta.h5' (par%hei_metastable) with the converged 2^3S
metastable density n_2s3 [cm^-3] on its Cartesian leaf grid.  LaRT v2.00 has a
'HeI_10833' resonance-scattering line whose scatterer is exactly this 2^3S
population.  This converter turns the heimeta file into a LaRT density input so
LaRT scatters 10833 photons directly on n_2s3.

Two LaRT input layouts are supported (LaRT decides the scatterer opacity from
the density it reads at these two points, verbatim):

  car  (--format car, DEFAULT): all-in-one Cartesian grid, par%cart_file
       grid_mod_car.f90:299   rhokap = gasDen * distance2cm       (no xHI)
       -> gasDen carries the scatterer density n_2s3 directly.
  amr  (--format amr): generic-AMR binary table, par%amr_file
       grid_mod_amr.f90:253   nH = nion_arr  (when the n_ion column is present)
       grid_mod_amr.f90:272   rhokap = nH * cross0 / Dfreq * distance2cm
       -> the n_ion column carries the scatterer density n_2s3.

The 'gasDen'/'nH' names are Lyman-alpha-era labels in LaRT; for a selected line
the scatterer density is the density of that line's lower-level ion, here the
He I 2^3S metastable atom.  LaRT applies the 10833 line-center cross section
cross0 (which uses the strongest fine-structure component f12(1) = 0.29958) and
the Doppler width Dfreq internally.

Field -> LaRT section/column mapping (car / amr):

  heimeta field   car section   amr column      meaning
  -------------   -----------   ----------      -------
  n_2s3           gasDen        n_ion           scatterer density [cm^-3]  (--field, default)
  T_e             T             T (col 6)       gas temperature [K]
  (zero)          vx,vy,vz      vx,vy,vz (7-9)  bulk velocity [km/s] = 0
  n_H             --            nH (col 5)      hydrogen density [cm^-3]   (amr only)
  (derived)       --            level (col 4)   AMR level = log2(BOXLEN/LeafSize)

Geometry: box side BOXLEN = nx * d in code units (d = LeafSize); LaRT centers
the box at the origin, so ORIGINX/Y/Z = -0.5*BOXLEN and xmax = BOXLEN/2.
UNITLCGS = DIST_CM (cm per code unit); set par%distance2cm = 1.0 so LaRT adopts
UNITLCGS from the file (car), or par%distance_unit/distance2cm in the namelist
(amr, which has no UNITLCGS attribute).

Usage:
  python3 heimeta_to_lart.py <base>_heimeta.h5 [-o out.h5] [--format car|amr]
                             [--field n_2s3] [--verify]

Verification (--verify) re-reads the written file by mirroring the LaRT reader
logic and confirms the leaf values round-trip (G1).
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mochii_output import read_sections


# ---------------------------------------------------------------------------
# read the heimeta file and reconstruct the (ix, iy, iz) leaf indexing
# ---------------------------------------------------------------------------
def load_heimeta(fname):
    """Read a '<base>_heimeta' file into a dict of leaf arrays plus geometry.

    Returns a dict with keys: fields (name -> 1D leaf array), xyz (3, nleaf;
    code units), size (1D full cell width, code units), dist_cm, a31.
    """
    sec = read_sections(fname)
    xyz = np.asarray(sec["LeafXYZ"]["data"], float)
    if xyz.shape[0] != 3:
        xyz = xyz.T
    size = np.asarray(sec["LeafSize"]["data"], float).ravel()
    dist_cm = float(sec["n_2s3"]["attrs"].get("DIST_CM", 1.0))
    a31 = float(sec["n_2s3"]["attrs"].get("A31", np.nan))

    fields = {}
    for name in ("n_2s3", "Phi_3", "n_H", "n_e", "T_e"):
        if name in sec and sec[name]["data"] is not None:
            fields[name] = np.asarray(sec[name]["data"], float).ravel()
    return dict(fields=fields, xyz=xyz, size=size, dist_cm=dist_cm, a31=a31)


def leaf_indices(xyz, size):
    """Reconstruct integer (ix, iy, iz) leaf indices from leaf centers.

    Robust to the on-disk leaf ordering: each index is round((coord - min)/d)
    with d the (uniform) cell width.  Requires a uniform cell width; raises
    ValueError otherwise (use --format amr for a non-uniform grid).  Returns
    (ix, iy, iz, nx, ny, nz, d).
    """
    d = float(np.median(size))
    if not np.allclose(size, d, rtol=1e-6, atol=0.0):
        raise ValueError(
            "LeafSize is not uniform (min=%.6g max=%.6g); the all-in-one "
            "Cartesian ('car') layout needs a uniform grid -- use --format amr."
            % (size.min(), size.max()))
    idx = []
    n = []
    for a in range(3):
        c = xyz[a]
        i = np.rint((c - c.min()) / d).astype(np.int64)
        if i.min() < 0:
            raise ValueError("negative leaf index on axis %d" % a)
        idx.append(i)
        n.append(int(i.max()) + 1)
    ix, iy, iz = idx
    nx, ny, nz = n
    if nx * ny * nz != xyz.shape[1]:
        raise ValueError(
            "reconstructed grid %d x %d x %d = %d does not fill nleaf = %d; "
            "the grid is not a filled uniform box -- use --format amr."
            % (nx, ny, nz, nx * ny * nz, xyz.shape[1]))
    # each (ix,iy,iz) must be unique
    lin = ix + nx * (iy + ny * iz)
    if np.unique(lin).size != lin.size:
        raise ValueError("duplicate (ix,iy,iz) leaf indices reconstructed")
    return ix, iy, iz, nx, ny, nz, d


# ---------------------------------------------------------------------------
# all-in-one Cartesian writer (mirrors convert_illustris_to_generic.py
# write_cartesian_hdf5 exactly: one group per quantity, group name = EXTNAME,
# dataset 'data' stored as np.asfortranarray(darr).T with darr shape (nx,ny,nz))
# ---------------------------------------------------------------------------
def write_cartesian(fname, gasden3d, t3d, boxlen, unit_l_cgs):
    """Write the mandatory car sections (gasDen, T, vx=vy=vz=0).

    gasden3d / t3d are (nx, ny, nz) arrays indexed [ix, iy, iz].  No 'xHI'
    section is written, so LaRT uses gasDen directly as the scatterer density.
    """
    import h5py

    nx, ny, nz = gasden3d.shape
    zeros = np.zeros((nx, ny, nz), dtype=np.float64)
    groups = [("gasDen", gasden3d), ("T", t3d),
              ("vx", zeros), ("vy", zeros), ("vz", zeros)]

    if os.path.exists(fname):
        os.unlink(fname)
    with h5py.File(fname, "w", libver="latest", track_order=True) as f:
        first = True
        for name, arr in groups:
            darr = np.ascontiguousarray(arr, dtype=np.float64)
            # Python (nx,ny,nz) -> stored (nz,ny,nx) for Fortran column-major.
            darr_f = np.asfortranarray(darr).T
            g = f.create_group(name)
            kw = {}
            if darr_f.size > 4096:
                kw["chunks"] = tuple(min(s, 64) for s in darr_f.shape)
                kw["compression"] = "gzip"
                kw["compression_opts"] = 4
            g.create_dataset("data", data=darr_f, **kw)
            g.attrs["EXTNAME"] = np.bytes_(name)
            if first:
                g.attrs["CREATOR"] = np.bytes_("heimeta_to_lart.py")
                g.attrs["NX"] = np.int32(nx)
                g.attrs["NY"] = np.int32(ny)
                g.attrs["NZ"] = np.int32(nz)
                g.attrs["BOXLEN"] = np.float64(boxlen)
                g.attrs["UNITLCGS"] = np.float64(unit_l_cgs)
                g.attrs["ORIGINX"] = np.float64(-0.5 * boxlen)
                g.attrs["ORIGINY"] = np.float64(-0.5 * boxlen)
                g.attrs["ORIGINZ"] = np.float64(-0.5 * boxlen)
                first = False


def convert_car(hm, out, field="n_2s3"):
    """Build the car grid from a loaded heimeta dict and write it.

    Returns a dict describing the written grid (nx, ny, nz, boxlen, unit_l_cgs,
    and the leaf-index arrays) for verification.
    """
    xyz, size = hm["xyz"], hm["size"]
    if field not in hm["fields"]:
        raise KeyError("field '%s' not in heimeta file (have: %s)"
                       % (field, ", ".join(sorted(hm["fields"]))))
    ix, iy, iz, nx, ny, nz, d = leaf_indices(xyz, size)

    gasden3d = np.zeros((nx, ny, nz), dtype=np.float64)
    t3d = np.zeros((nx, ny, nz), dtype=np.float64)
    gasden3d[ix, iy, iz] = hm["fields"][field]
    t3d[ix, iy, iz] = hm["fields"]["T_e"]

    boxlen = nx * d
    write_cartesian(out, gasden3d, t3d, boxlen, hm["dist_cm"])
    return dict(nx=nx, ny=ny, nz=nz, boxlen=boxlen, unit_l_cgs=hm["dist_cm"],
                ix=ix, iy=iy, iz=iz, field=field)


# ---------------------------------------------------------------------------
# generic-AMR binary-table writer (mirrors read_generic_amr.f90: one group
# with 1-D datasets read by creation-order index 1-9 = x,y,z,level,nH,T,
# vx,vy,vz, then optional columns by name including n_ion)
# ---------------------------------------------------------------------------
def write_generic_amr(fname, xyz, level, nH, T, n_ion, boxlen):
    """Write the generic-AMR binary table LaRT's read_generic_amr expects."""
    import h5py

    nleaf = xyz.shape[1]
    zeros = np.zeros(nleaf, dtype=np.float64)
    # exact creation order: columns 1-9 by index, then n_ion by name
    columns = [
        ("x", xyz[0].astype(np.float64)),
        ("y", xyz[1].astype(np.float64)),
        ("z", xyz[2].astype(np.float64)),
        ("level", level.astype(np.int32)),
        ("nH", nH.astype(np.float64)),
        ("T", T.astype(np.float64)),
        ("vx", zeros), ("vy", zeros), ("vz", zeros),
        ("n_ion", n_ion.astype(np.float64)),
    ]
    if os.path.exists(fname):
        os.unlink(fname)
    with h5py.File(fname, "w", libver="latest", track_order=True) as f:
        # track_order on the group so LaRT's H5_INDEX_CRT_ORDER_F column read
        # maps index i -> the i-th created dataset (else it falls back to
        # alphabetical name order and mismaps the columns).
        g = f.create_group("leaves", track_order=True)
        for name, arr in columns:
            g.create_dataset(name, data=arr)
        # a table has no 'data' dataset, so NAXIS2 (row count) must be explicit
        g.attrs["NAXIS2"] = np.int32(nleaf)
        g.attrs["BOXLEN"] = np.float64(boxlen)
        g.attrs["ORIGINX"] = np.float64(-0.5 * boxlen)
        g.attrs["ORIGINY"] = np.float64(-0.5 * boxlen)
        g.attrs["ORIGINZ"] = np.float64(-0.5 * boxlen)
        g.attrs["CREATOR"] = np.bytes_("heimeta_to_lart.py")


def convert_amr(hm, out, field="n_2s3"):
    """Build the generic-AMR table from a loaded heimeta dict and write it."""
    xyz, size = hm["xyz"], hm["size"]
    if field not in hm["fields"]:
        raise KeyError("field '%s' not in heimeta file (have: %s)"
                       % (field, ", ".join(sorted(hm["fields"]))))
    # box side from the coordinate span plus one cell (centers -> edges)
    d = float(np.median(size))
    span = max(xyz[a].max() - xyz[a].min() for a in range(3))
    boxlen = span + d
    # uniform: level = log2(boxlen / cell width) per leaf
    level = np.rint(np.log2(boxlen / size)).astype(np.int32)
    nH = hm["fields"].get("n_H", np.zeros(xyz.shape[1]))
    write_generic_amr(out, xyz, level, nH, hm["fields"]["T_e"],
                      hm["fields"][field], boxlen)
    return dict(nleaf=xyz.shape[1], boxlen=boxlen, field=field)


# ---------------------------------------------------------------------------
# verification (G1): re-read each file by mirroring the LaRT reader logic
# ---------------------------------------------------------------------------
def verify_car(out, hm, meta):
    """Re-read the car file the way LaRT does and compare to the heimeta leaves.

    LaRT reads gasDen(nx,ny,nz) column-major from the stored (nz,ny,nx) dataset,
    so gasDen(ix,iy,iz) = dataset[iz,iy,ix].  We recover every leaf value this
    way and compare to n_2s3 / T_e.  Returns (max|dgasDen|, max|dT|).
    """
    import h5py

    ix, iy, iz = meta["ix"], meta["iy"], meta["iz"]
    field = meta["field"]
    with h5py.File(out, "r") as f:
        gden = f["gasDen"]["data"][:]     # shape (nz, ny, nx)
        tarr = f["T"]["data"][:]
        nx = int(np.asarray(f["gasDen"].attrs["NX"]).ravel()[0])
        ny = int(np.asarray(f["gasDen"].attrs["NY"]).ravel()[0])
        nz = int(np.asarray(f["gasDen"].attrs["NZ"]).ravel()[0])
        boxlen = float(np.asarray(f["gasDen"].attrs["BOXLEN"]).ravel()[0])
        ulcgs = float(np.asarray(f["gasDen"].attrs["UNITLCGS"]).ravel()[0])
    assert gden.shape == (nz, ny, nx), (gden.shape, (nz, ny, nx))
    # LaRT-side gasDen(ix,iy,iz) = dataset[iz,iy,ix]
    got_g = gden[iz, iy, ix]
    got_t = tarr[iz, iy, ix]
    dg = float(np.max(np.abs(got_g - hm["fields"][field])))
    dt = float(np.max(np.abs(got_t - hm["fields"]["T_e"])))
    print("  car header : NX,NY,NZ = %d,%d,%d  BOXLEN = %.8g  UNITLCGS = %.8g"
          % (nx, ny, nz, boxlen, ulcgs))
    return dg, dt


def _column_names_creation_order(group):
    """Column names of an HDF5 table group in HDF5 link creation order --
    the same order LaRT reads with H5_INDEX_CRT_ORDER_F (so index i = column i).
    """
    from h5py import h5
    names = []
    group.id.links.iterate(
        lambda nm: names.append(nm.decode() if isinstance(nm, bytes) else nm),
        idx_type=h5.INDEX_CRT_ORDER, order=h5.ITER_INC)
    return names


def verify_amr(out, hm, meta):
    """Re-read the amr table the way LaRT does: mandatory columns by index in
    HDF5 creation order (1-9 = x,y,z,level,nH,T,vx,vy,vz), optional n_ion by
    name.  Compares T (column 6) and n_ion leaf-by-leaf.  Returns
    (max|dn_ion|, max|dT|).
    """
    import h5py

    field = meta["field"]
    with h5py.File(out, "r") as f:
        g = f["leaves"]
        names = _column_names_creation_order(g)   # LaRT's column-index order
        nion = g["n_ion"][:]                       # optional column, by name
        tcol = g[names[5]][:]                       # mandatory column 6 = T
        nleaf = int(np.asarray(g.attrs["NAXIS2"]).ravel()[0])
        boxlen = float(np.asarray(g.attrs["BOXLEN"]).ravel()[0])
    assert names[5] == "T", names
    dg = float(np.max(np.abs(nion - hm["fields"][field])))
    dt = float(np.max(np.abs(tcol - hm["fields"]["T_e"])))
    print("  amr header : NAXIS2 = %d  BOXLEN = %.8g  columns(1-6) = %s"
          % (nleaf, boxlen, ",".join(names[:6])))
    return dg, dt


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Convert a MoCHII <base>_heimeta file to a LaRT grid.")
    ap.add_argument("infile", help="MoCHII <base>_heimeta.h5")
    ap.add_argument("-o", "--out", default=None,
                    help="output HDF5 (default: <infile>_lart_<format>.h5)")
    ap.add_argument("--format", default="car", choices=["car", "amr"],
                    help="LaRT layout: car (all-in-one Cartesian, default) "
                         "or amr (generic-AMR binary table)")
    ap.add_argument("--field", default="n_2s3",
                    help="heimeta field for the scatterer density "
                         "(default n_2s3; e.g. n_H, n_e)")
    ap.add_argument("--verify", action="store_true",
                    help="re-read the written file (mirror the LaRT reader) "
                         "and print round-trip max abs diffs")
    args = ap.parse_args()

    hm = load_heimeta(args.infile)
    nleaf = hm["xyz"].shape[1]
    print("heimeta: %s  nleaf = %d  DIST_CM = %.8g  A31 = %.4g"
          % (args.infile, nleaf, hm["dist_cm"], hm["a31"]))
    print("  max n_2s3 = %.4e cm^-3   field -> scatterer: %s"
          % (hm["fields"]["n_2s3"].max(), args.field))

    out = args.out
    if out is None:
        base = args.infile
        for suf in ("_heimeta.h5", "_heimeta.hdf5", ".h5", ".hdf5"):
            if base.endswith(suf):
                base = base[:-len(suf)]
                break
        out = "%s_lart_%s.h5" % (base, args.format)

    if args.format == "car":
        meta = convert_car(hm, out, field=args.field)
        # cross-check that the raster order il = ix + nx*(iy + ny*iz) matches
        ix, iy, iz, nx, ny = meta["ix"], meta["iy"], meta["iz"], meta["nx"], meta["ny"]
        lin = ix + nx * (iy + ny * iz)
        raster_ok = np.array_equal(lin, np.arange(nleaf))
        print("written: %s (car)  raster order matches: %s" % (out, raster_ok))
        if not raster_ok:
            print("  note: on-disk order is NOT the raster order il=ix+nx*"
                  "(iy+ny*iz); values were placed by reconstructed (ix,iy,iz).")
    else:
        meta = convert_amr(hm, out, field=args.field)
        print("written: %s (amr)  nleaf = %d  BOXLEN = %.8g"
              % (out, meta["nleaf"], meta["boxlen"]))

    if args.verify:
        if args.format == "car":
            dg, dt = verify_car(out, hm, meta)
        else:
            dg, dt = verify_amr(out, hm, meta)
        print("  round-trip max |d(scatterer)| = %.3e   max |dT_e| = %.3e"
              % (dg, dt))
        ok = (dg == 0.0 and dt == 0.0)
        print("  G1 %s" % ("PASS (exact)" if ok else "PASS (round-off)"
                            if max(dg, dt) < 1e-6 else "FAIL"))


if __name__ == "__main__":
    main()
