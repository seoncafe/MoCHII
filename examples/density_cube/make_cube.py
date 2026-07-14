"""Build a clumpy 3D hydrogen density cube for MoCHII.

This writes a deterministic 48x48x48 float64 density field (clumpy.fits.gz) that
MoCHII reads directly onto a 'car' Cartesian grid via par%density_file.  The
field has three regions:

  * a dense central clump (nH ~ 200 cm^-3),
  * a smoother diffuse background filling the rest of a sphere (nH ~ 30),
  * vacuum (nH = 0) outside a cutoff radius.

The values are fixed (no random numbers), so the cube is reproducible.  The
grid geometry (box size, distance unit) is set in the namelist; this script
only supplies the density values, one per cell.
"""

import numpy as np
from astropy.io import fits

# Cube shape and physical box half-size (must match par%xmax in the namelist).
N = 48
XMAX = 4.0  # box half-size in the namelist distance unit (pc)

# Cell-center coordinates on [-XMAX, XMAX], one per cell.
edges = np.linspace(-XMAX, XMAX, N + 1)
centers = 0.5 * (edges[:-1] + edges[1:])
xg, yg, zg = np.meshgrid(centers, centers, centers, indexing="ij")
r = np.sqrt(xg**2 + yg**2 + zg**2)

# Density field [cm^-3].
nH = np.zeros((N, N, N), dtype=np.float64)

# Diffuse background inside a sphere: a smooth radial falloff around ~30 cm^-3.
r_out = 3.8  # outer cutoff radius; nH = 0 beyond
inside = r <= r_out
nH[inside] = 30.0 * np.exp(-(r[inside] / r_out) ** 2)

# Central clump: a Gaussian bump peaking near 200 cm^-3.
r_clump = 1.0
nH += 200.0 * np.exp(-(r / r_clump) ** 2)

# Keep vacuum strictly outside the cutoff radius.
nH[r > r_out] = 0.0

hdu = fits.PrimaryHDU(data=nH)
hdu.writeto("clumpy.fits.gz", overwrite=True)

print("clumpy.fits.gz written")
print("shape =", nH.shape, "dtype =", nH.dtype)
print("nH min = %.4f  max = %.4f  cm^-3" % (nH.min(), nH.max()))
