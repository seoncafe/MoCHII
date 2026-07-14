"""Backend-agnostic rasterizer over cubic AMR leaves.

A "leaf field" is a scalar defined on the leaf cells of an octree (or of a
uniform Cartesian grid, which is a degenerate octree): a set of cubic cells
described by their centers (cx, cy, cz), a half-size ``half`` (so the cell
spans [c-half, c+half] on each axis), and one scalar ``value`` for each cell.
All four are 1D numpy arrays of length nleaf (``half`` may also be a single
scalar shared by every cell, e.g. for a uniform grid).

This module rasterizes such a field two ways:

* ``leaf_slice`` -- a true PLANAR SLICE.  For a plane at ``axis = coord`` it
  returns, at every pixel, the value of the leaf cube that CONTAINS that
  pixel.  It is not a projection: values are painted by assignment, not
  accumulated along the line of sight.  Cells tile space without overlap, so
  each covered pixel receives its containing cell's value; pixels covered by
  no cell stay NaN.  Contrast this with the flux-conserving deposit in
  mochii_output.py, which sums leaf luminosities into a projected map.
* ``leaf_slice_nn`` -- a nearest-leaf fallback used when the half-sizes are
  unknown.  Each pixel takes the value of the nearest leaf center.

Conventions match tools/python/mochii_output.py so the two agree:

* Plane axes: ``ia = "xyz".index(axis); iu, iv = [i for i in range(3) if
  i != ia]``.  So axis='z' gives in-plane axes (u, v) = (x, y); axis='y'
  gives (x, z); axis='x' gives (y, z).
* Images are indexed ``img[iu_pixel, iv_pixel]`` and are meant to be shown
  with ``imshow(img.T, origin='lower', extent=(umin, umax, vmin, vmax))``.

Only numpy is required.  ``leaf_slice_nn`` imports scipy lazily.
"""

import numpy as np


# ---------------------------------------------------------------------------
def plane_axes(axis):
    """Return (ia, iu, iv): the plane-normal axis index and the two in-plane
    axis indices, following the mochii_output.py convention.

    axis is one of 'x', 'y', 'z'.  Example: axis='z' -> (2, 0, 1).
    """
    ia = "xyz".index(axis)
    iu, iv = [i for i in range(3) if i != ia]
    return ia, iu, iv


def _as_half(half, shape):
    """Broadcast a scalar or array half-size to a float array of `shape`."""
    half = np.asarray(half, dtype=float)
    if half.ndim == 0:
        half = np.full(shape, float(half))
    return half


# ---------------------------------------------------------------------------
def leaf_slice(cx, cy, cz, half, value, axis="z", coord=0.0,
               bounds=None, npix=512):
    """True planar slice of a leaf field: the actual cell value at axis=coord.

    For each pixel at (u, v, coord) the returned image holds the value of the
    leaf cube that contains it.  Leaves straddling the plane
    (``abs(coord - c_axis) <= half``) are selected, and each such leaf's
    in-plane footprint [cu-h, cu+h] x [cv-h, cv+h] is painted onto the pixel
    grid by ASSIGNMENT (not area-weighted accumulation), because leaves tile
    space without overlap.

    Parameters
    ----------
    cx, cy, cz : 1D arrays of leaf-center coordinates [code units].
    half       : leaf half-size (scalar or 1D array); a cell spans
                 [c-half, c+half] on each axis.
    value      : 1D array of the scalar to slice.
    axis       : plane-normal axis 'x' | 'y' | 'z'.
    coord      : position of the plane along `axis` [code units].
    bounds     : optional (umin, umax, vmin, vmax) cutout window.  Default is
                 the full in-plane extent of ALL leaves
                 (umin=min(cu-h), umax=max(cu+h), vmin=min(cv-h),
                 vmax=max(cv+h)) -- a stable frame that does not jump between
                 slices.
    npix       : image side length in pixels.

    Returns
    -------
    (img, extent) : img has shape (npix, npix) indexed [u, v], initialized to
    NaN (uncovered pixels stay NaN); extent = (umin, umax, vmin, vmax).
    """
    cx = np.asarray(cx, float)
    cy = np.asarray(cy, float)
    cz = np.asarray(cz, float)
    value = np.asarray(value, float)
    half = _as_half(half, cx.shape)

    centers = (cx, cy, cz)
    ia, iu, iv = plane_axes(axis)
    cu, cv, cw = centers[iu], centers[iv], centers[ia]

    if bounds is None:
        umin = float((cu - half).min())
        umax = float((cu + half).max())
        vmin = float((cv - half).min())
        vmax = float((cv + half).max())
    else:
        umin, umax, vmin, vmax = (float(b) for b in bounds)
    extent = (umin, umax, vmin, vmax)

    du = (umax - umin) / npix
    dv = (vmax - vmin) / npix
    img = np.full((npix, npix), np.nan)

    sel = np.abs(coord - cw) <= half        # leaves straddling the plane
    for il in np.nonzero(sel)[0]:
        h = half[il]
        i0 = max(int((cu[il] - h - umin) / du), 0)
        i1 = min(int((cu[il] + h - umin) / du) + 1, npix)
        j0 = max(int((cv[il] - h - vmin) / dv), 0)
        j1 = min(int((cv[il] + h - vmin) / dv) + 1, npix)
        if i1 <= i0 or j1 <= j0:
            continue
        img[i0:i1, j0:j1] = value[il]
    return img, extent


# ---------------------------------------------------------------------------
def leaf_slice_nn(cx, cy, cz, value, axis="z", coord=0.0, bounds=None,
                  npix=512):
    """Nearest-leaf planar slice, for when half-sizes are unknown.

    Builds a KD-tree on the leaf centers and, for each pixel center at
    (u, v, coord), takes the value of the nearest leaf.

    Parameters mirror ``leaf_slice`` (no ``half``).  ``bounds`` defaults to the
    full center extent of the two in-plane axes
    (min/max of cu and of cv).  Returns (img, extent), img indexed [u, v].
    """
    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise ImportError("leaf_slice_nn requires scipy "
                          "(scipy.spatial.cKDTree); please install scipy") \
            from exc
    cx = np.asarray(cx, float)
    cy = np.asarray(cy, float)
    cz = np.asarray(cz, float)
    value = np.asarray(value, float)

    centers = (cx, cy, cz)
    ia, iu, iv = plane_axes(axis)
    cu, cv = centers[iu], centers[iv]

    if bounds is None:
        umin, umax = float(cu.min()), float(cu.max())
        vmin, vmax = float(cv.min()), float(cv.max())
    else:
        umin, umax, vmin, vmax = (float(b) for b in bounds)
    extent = (umin, umax, vmin, vmax)

    du = (umax - umin) / npix
    dv = (vmax - vmin) / npix
    us = umin + (np.arange(npix) + 0.5) * du
    vs = vmin + (np.arange(npix) + 0.5) * dv
    uu, vv = np.meshgrid(us, vs, indexing="ij")     # (npix, npix) indexed [u, v]

    pts = np.empty((npix * npix, 3))
    pts[:, ia] = coord
    pts[:, iu] = uu.ravel()
    pts[:, iv] = vv.ravel()

    tree = cKDTree(np.column_stack([cx, cy, cz]))
    _, jj = tree.query(pts)
    img = value[jj].reshape(npix, npix)
    return img, extent


# ---------------------------------------------------------------------------
def leaf_box_mask(cx, cy, cz, half, box):
    """Boolean mask: True where a leaf cube overlaps the axis-aligned box.

    box = (xmin, xmax, ymin, ymax, zmin, zmax); ``half`` may be a scalar or a
    1D array.  Overlap is the standard AABB test (touching counts as inside).
    """
    cx = np.asarray(cx, float)
    cy = np.asarray(cy, float)
    cz = np.asarray(cz, float)
    half = _as_half(half, cx.shape)
    xmin, xmax, ymin, ymax, zmin, zmax = box
    return ((cx + half >= xmin) & (cx - half <= xmax) &
            (cy + half >= ymin) & (cy - half <= ymax) &
            (cz + half >= zmin) & (cz - half <= zmax))


# ---------------------------------------------------------------------------
def plot_slice(img, extent, ax=None, log=False, cmap="inferno",
               cbar_label="", title="", axis="z"):
    """Show a slice image with the mochii_output.py display convention.

    Draws ``imshow(img.T, origin='lower', extent=extent, cmap=cmap)``.  When
    ``log`` is set, non-positive values are masked and log10 is shown, with the
    color-bar label prefixed by a log10 tag.  NaN pixels render blank.  x/y
    axis labels are the in-plane axis names implied by ``axis``.

    Returns (fig, ax).  A new figure is created when ax is None.
    """
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

    im = ax.imshow(show.T, origin="lower", extent=extent, cmap=cmap)
    fig.colorbar(im, ax=ax, label=cbar_label)

    names = "xyz"
    _, iu, iv = plane_axes(axis)
    ax.set_xlabel(f"{names[iu]} (code units)")
    ax.set_ylabel(f"{names[iv]} (code units)")
    if title:
        ax.set_title(title)
    return fig, ax


# ---------------------------------------------------------------------------
def _smoke():
    """Build a synthetic 2-level leaf set and print a slice summary."""
    cx, cy, cz, half, val = [], [], [], [], []

    # level-1 cells (half = 0.5) filling the box [-1, 1]^3, except the
    # (+,+,+) octant which is refined into level-2 cells (half = 0.25).
    for b in range(8):
        ix, iy, iz = b & 1, (b >> 1) & 1, (b >> 2) & 1
        px, py, pz = (2 * ix - 1) * 0.5, (2 * iy - 1) * 0.5, (2 * iz - 1) * 0.5
        if ix == 1 and iy == 1 and iz == 1:
            for c in range(8):
                jx, jy, jz = c & 1, (c >> 1) & 1, (c >> 2) & 1
                cx.append(px + (2 * jx - 1) * 0.25)
                cy.append(py + (2 * jy - 1) * 0.25)
                cz.append(pz + (2 * jz - 1) * 0.25)
                half.append(0.25)
                val.append(2.0)             # level-2 value
        else:
            cx.append(px); cy.append(py); cz.append(pz)
            half.append(0.5)
            val.append(1.0)                 # level-1 value

    cx = np.array(cx); cy = np.array(cy); cz = np.array(cz)
    half = np.array(half); val = np.array(val)

    img, extent = leaf_slice(cx, cy, cz, half, val, axis="z",
                             coord=0.6, npix=128)
    finite = np.isfinite(img)
    print(f"leaf_slice smoke: {cx.size} leaves, img shape {img.shape}, "
          f"extent {extent}, finite fraction {finite.mean():.3f}, "
          f"value range {np.nanmin(img):g}-{np.nanmax(img):g}")


if __name__ == "__main__":
    _smoke()
