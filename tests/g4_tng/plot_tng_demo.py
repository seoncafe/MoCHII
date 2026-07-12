#!/usr/bin/env python3
"""G4b demo maps: slice through the source plane of the TNG cutout.

Nearest-leaf-center sampling on a uniform pixel grid (adequate for maps;
level boundaries are approximate).  Panels: nH, x_HI, Te, cell size
(shows the I-front re-refinement).
"""
import numpy as np
import h5py
from scipy.spatial import cKDTree

Z0 = 4.6875          # source plane [kpc]
SRC = (-4.6875, -15.9375)
BOX = 60.0
NPIX = 512

f = h5py.File("g4_tng_rates.h5", "r")
xyz = f["LeafXYZ"]["data"][:]          # (3, nleaf)
xhi = f["x_HI"]["data"][:]
te = f["T_e"]["data"][:]
ne = f["n_e"]["data"][:]

# nH per leaf: reconstruct from n_e closure is unreliable in neutral gas;
# reload from the input grid by sampling its own tree (uniform L5 file).
from astropy.io import fits
d0 = fits.open("tng_uniform.fits")[1].data
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
# cell size from nearest-neighbor distance among leaf centers: use the
# distance to the 2nd neighbor as ~cell size proxy; better: infer level
# by alignment as in the other checks
half = np.zeros(xyz.shape[1])
for L in (7, 6, 5):
    step = BOX/2**L
    on = np.all(np.abs(((xyz + BOX/2)/step - 0.5)
                       - np.round((xyz + BOX/2)/step - 0.5)) < 1e-6,
                axis=0) & (half == 0)
    half[on] = step/2
size_map = (2*half[idx]).reshape(NPIX, NPIX)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

fig, axes = plt.subplots(2, 2, figsize=(12, 11), sharex=True, sharey=True)
ext = [-BOX/2, BOX/2, -BOX/2, BOX/2]

im = axes[0, 0].imshow(nH_map.T, origin="lower", extent=ext,
                       norm=LogNorm(vmin=1e-5, vmax=1.0), cmap="viridis")
plt.colorbar(im, ax=axes[0, 0], label=r"$n_{\rm H}$ [cm$^{-3}$]")
axes[0, 0].set_title("gas density (TNG cutout)")

im = axes[0, 1].imshow(xhi_map.T, origin="lower", extent=ext,
                       norm=LogNorm(vmin=1e-5, vmax=1.0), cmap="magma_r")
plt.colorbar(im, ax=axes[0, 1], label=r"$x_{\rm HI}$")
axes[0, 1].set_title(r"neutral fraction (converged)")

im = axes[1, 0].imshow(te_map.T, origin="lower", extent=ext,
                       vmin=3e3, vmax=3e4, cmap="inferno")
plt.colorbar(im, ax=axes[1, 0], label=r"$T_e$ [K]")
axes[1, 0].set_title(r"electron temperature")

im = axes[1, 1].imshow(size_map.T, origin="lower", extent=ext,
                       cmap="cividis_r")
plt.colorbar(im, ax=axes[1, 1], label=r"cell size [kpc]")
axes[1, 1].set_title("cell size (I-front re-refinement)")

for ax in axes.flat:
    ax.plot(*SRC, "w*", ms=12, mec="k")
    ax.set_xlabel(r"$x$ [kpc]")
    ax.set_ylabel(r"$y$ [kpc]")
fig.suptitle(r"MoCHII on an Illustris-TNG cutout: $Q_{\rm H}=10^{53}$"
             r" s$^{-1}$ source, $z = 4.7$ kpc plane", y=0.995)
fig.tight_layout()
fig.savefig("g4_tng_maps.png", dpi=130)
print("wrote g4_tng_maps.png")

# summary numbers
gas = ne > 0
w = ne[gas]
print(f"nleaf = {xyz.shape[1]}")
print(f"volume-mean x_HI (gas cells) = {xhi[gas].mean():.4f}")
print(f"ne-weighted mean Te = {np.sum(te[gas]*w)/np.sum(w):.0f} K")
