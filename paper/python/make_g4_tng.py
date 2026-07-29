import sys
#!/usr/bin/env python3
"""G4b demo maps: slice through the source plane of the TNG cutout.

Nearest-leaf-center sampling on a uniform pixel grid (adequate for maps;
level boundaries are approximate).  Panels: nH, x_HI, Te, cell size
(shows the I-front re-refinement).
"""
import numpy as np
import os as _os
# Figures are written to paper/figures/ regardless of the working directory.
FIGDIR = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "figures")
import h5py
from scipy.spatial import cKDTree

Z0 = 4.6875          # source plane [kpc]
SRC = (-4.6875, -15.9375)
BOX = 60.0
NPIX = 512

DDIR = "../../tests/g4_tng/"
f = h5py.File(DDIR + "g4_tng_rates.h5", "r")
xyz = f["LeafXYZ"]["data"][:]          # (3, nleaf)
xhi = f["x_HI"]["data"][:]
te = f["T_e"]["data"][:]
ne = f["n_e"]["data"][:]

# nH per leaf: reconstruct from n_e closure is unreliable in neutral gas;
# reload from the input grid by sampling its own tree (uniform L5 file).
from astropy.io import fits
d0 = fits.open(DDIR + "tng_uniform.fits")[1].data
tree0 = cKDTree(np.stack([d0["x"], d0["y"], d0["z"]], axis=1))

tree = cKDTree(xyz.T)
g = np.linspace(-BOX/2, BOX/2, NPIX)
X, Y = np.meshgrid(g, g, indexing="ij")
pts = np.stack([X.ravel(), Y.ravel(), np.full(X.size, Z0)], axis=1)
_, idx = tree.query(pts)
_, idx0 = tree0.query(pts)

nH_map = d0["nH"][idx0].reshape(NPIX, NPIX)
xhi_map = xhi[idx].reshape(NPIX, NPIX)
te_map = te[idx].reshape(NPIX, NPIX)

# cell-size panel: draw the actual leaf cells cut by the slice plane, colored
# by discrete refinement level, in the style of the LaRT AMR_grid slice_plot.
# The leaf half-size comes from the file's own LeafSize block (grid units =
# kpc here), so no proxy is needed.  A leaf intersects the plane z = Z0 when
# cz - h <= Z0 < cz + h; the half-open test picks exactly one leaf per column,
# so the selected footprints tile the box (verified: area ratio 1.0000).
half_all = f["LeafSize"]["data"][:].ravel() / 2.0
in_slice = (xyz[2] - half_all <= Z0) & (Z0 < xyz[2] + half_all)
hc = half_all[in_slice]
xc = xyz[0][in_slice]
yc = xyz[1][in_slice]
cell_level = np.rint(np.log2(BOX / (2.0 * hc))).astype(int)   # 5, 6, 7
LEVELS = np.array([5, 6, 7])
LEVEL_SIZE = BOX / 2.0 ** LEVELS                              # kpc: 1.88, 0.94, 0.47

cell_polys = np.empty((in_slice.sum(), 4, 2))
cell_polys[:, 0] = np.column_stack([xc - hc, yc - hc])
cell_polys[:, 1] = np.column_stack([xc + hc, yc - hc])
cell_polys[:, 2] = np.column_stack([xc + hc, yc + hc])
cell_polys[:, 3] = np.column_stack([xc - hc, yc + hc])

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
_HERE_AX = _os.path.dirname(_os.path.abspath(__file__))
sys.path.insert(0, _HERE_AX)
from axis_style import square_ticks   # noqa: E402

# Paper figures: axis, tick, and title fonts sized to ~85% of the body text.
_S = 2.02
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "figure.titlesize": _b * 1.2 * _S,
})
from matplotlib.colors import LogNorm, BoundaryNorm
from matplotlib.collections import PolyCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, axes = plt.subplots(2, 2, figsize=(12, 11), sharex=True, sharey=True)
ext = [-BOX/2, BOX/2, -BOX/2, BOX/2]


def cbar(im, ax, label):
    """Colorbar of exactly the drawn image height.

    The panels have a fixed data aspect, so the axes box is taller than the
    square image and a plain colorbar overshoots it top and bottom.  The
    axes divider tracks the aspect-adjusted box instead.
    """
    cax = make_axes_locatable(ax).append_axes("right", size="4%", pad=0.06)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label(label)
    return cb

im = axes[0, 0].imshow(nH_map.T, origin="lower", extent=ext,
                       norm=LogNorm(vmin=1e-5, vmax=1.0), cmap="viridis",
                       rasterized=True)
cbar(im, axes[0, 0], r"$n_{\rm H}$ [cm$^{-3}$]")
axes[0, 0].set_title("gas density (TNG cutout)")

im = axes[0, 1].imshow(xhi_map.T, origin="lower", extent=ext,
                       norm=LogNorm(vmin=1e-5, vmax=1.0), cmap="magma_r",
                       rasterized=True)
cbar(im, axes[0, 1], r"$x_{\rm HI}$")
axes[0, 1].set_title(r"neutral fraction (converged)")

im = axes[1, 0].imshow(te_map.T, origin="lower", extent=ext,
                       vmin=3e3, vmax=3e4, cmap="inferno",
                       rasterized=True)
cbar(im, axes[1, 0], r"$T_e$ [K]")
axes[1, 0].set_title(r"electron temperature")

# discrete color per refinement level, with the cell boundaries drawn
cmap_lvl = plt.get_cmap("cividis_r", len(LEVELS))
norm_lvl = BoundaryNorm(np.append(LEVELS - 0.5, LEVELS[-1] + 0.5), len(LEVELS))
pc = PolyCollection(cell_polys, cmap=cmap_lvl, norm=norm_lvl,
                    edgecolors="0.15", linewidths=0.08, closed=True,
                    rasterized=True)
pc.set_array(cell_level)
axes[1, 1].add_collection(pc)
axes[1, 1].set_xlim(-BOX/2, BOX/2)
axes[1, 1].set_ylim(-BOX/2, BOX/2)
cb = cbar(pc, axes[1, 1], r"cell size [kpc]")
cb.set_ticks(LEVELS)
cb.set_ticklabels([f"{s:.2f}" for s in LEVEL_SIZE])
axes[1, 1].set_title("cell size (I-front re-refinement)")

for ax in axes.flat:
    ax.plot(*SRC, "w*", ms=12, mec="k")
    ax.set_xlabel(r"$x$ [kpc]")
    ax.set_ylabel(r"$y$ [kpc]")
    square_ticks(ax)
fig.suptitle(r"MoCHII on an Illustris-TNG cutout: $Q_{\rm H}=10^{53}$"
             r" s$^{-1}$ source, $z = 4.7$ kpc plane", y=0.995)
fig.tight_layout()
fig.savefig(_os.path.join(FIGDIR, "g4_tng_maps.pdf"), dpi=200, bbox_inches="tight")
print("wrote g4_tng_maps.pdf")

# summary numbers
gas = ne > 0
w = ne[gas]
print(f"nleaf = {xyz.shape[1]}")
print(f"volume-mean x_HI (gas cells) = {xhi[gas].mean():.4f}")
print(f"ne-weighted mean Te = {np.sum(te[gas]*w)/np.sum(w):.0f} K")
