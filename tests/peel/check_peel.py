"""Peel-off gate (dust-free Stromgren): the image-integrated direct flux
must equal the independent reference sum_b L_b <e^{-tau_b}>/(4 pi d^2),
with tau_b integrated radially through the converged x_HI/x_HeI/x_HeII
profiles from the rates file (VFKY96 cross sections, same as the code).

The reference uses the angle-averaged escape through the sphere, which
for a central point source is the radial tau in every direction — exact
for this geometry.

Run:  python3 check_peel.py peel_free
"""

import sys
import numpy as np
import h5py

BASE = sys.argv[1] if len(sys.argv) > 1 else "peel_free"

MB2CM2 = 1.0e-18
ETH_HI, ETH_HEI, ETH_HEII = 13.598, 24.587, 54.418


def vfky96(E, Eth, E0, s0, ya, P, yw, y0, y1):
    E = np.asarray(E, float)
    x = E / E0 - y0
    z = np.sqrt(x * x + y1 * y1)
    Q = 5.5 - 0.5 * P
    Fy = ((x - 1.0)**2 + yw**2) * z**(-Q) * (1.0 + np.sqrt(z / ya))**(-P)
    return np.where(E >= Eth, s0 * Fy * MB2CM2, 0.0)


def sig_HI(E):
    return vfky96(E, ETH_HI, 4.298e-1, 5.475e4, 3.288e1, 2.963, 0, 0, 0)


def sig_HeI(E):
    return vfky96(E, ETH_HEI, 1.361e1, 9.492e2, 1.469, 3.188,
                  2.039, 4.434e-1, 2.136)


def sig_HeII(E):
    return vfky96(E, ETH_HEII, 1.720, 1.369e4, 3.288e1, 2.963, 0, 0, 0)


# --- converged state -------------------------------------------------------
with h5py.File(f"{BASE}_rates.h5", "r") as f:
    xyz = f["LeafXYZ/data"][:]
    xhi = f["x_HI/data"][:]
    xhe1 = f["x_HeI/data"][:]
    xhe2 = f["x_HeII/data"][:]
    Eb = f["E_bin/data"][:]
    dist_cm = f["Gamma_HI"].attrs["DIST_CM"]
    dist_cm = float(np.asarray(dist_cm).ravel()[0])
r = np.sqrt((xyz**2).sum(axis=0))

# radial profiles (uniform sphere, nH = 100 within R = 4 pc)
NH, RSPH, YHE = 100.0, 4.0, 0.1
nr = 400
redge = np.linspace(0.0, RSPH, nr + 1)
rmid = 0.5 * (redge[1:] + redge[:-1])
xhi_r = np.ones(nr)
xhe1_r = np.ones(nr)
xhe2_r = np.zeros(nr)
for i in range(nr):
    m = (r >= redge[i]) & (r < redge[i + 1])
    if m.any():
        xhi_r[i] = xhi[m].mean()
        xhe1_r[i] = xhe1[m].mean()
        xhe2_r[i] = xhe2[m].mean()
dr_cm = (redge[1] - redge[0]) * dist_cm

# tau_b along a radius; L_b: Planck(4e4 K) over the stored bins,
# normalized to the band luminosity of the input.
LBAND, TSTAR = 3.177837e38, 4.0e4


def planck_w(E, T):
    x = E * 1.602176634e-12 / (1.380649e-16 * T)
    return E**3 / np.expm1(x)


wb = planck_w(Eb, TSTAR)
Lb = wb / wb.sum() * LBAND
tau_b = np.array([
    (NH * (xhi_r * sig_HI(E) + YHE * (xhe1_r * sig_HeI(E)
     + xhe2_r * sig_HeII(E))) * dr_cm).sum() for E in Eb])
d_code = None

# --- image -----------------------------------------------------------------
with h5py.File(f"{BASE}_image.h5", "r") as f:
    img = f["direct_ion/data"][:]
    at = {k: float(np.asarray(v).ravel()[0]) for k, v in
          f["direct_ion"].attrs.items() if k in
          ("DXIM", "DIST")}
    names = list(f.keys())
sr_pix = np.deg2rad(at["DXIM"])**2
d_cm = at["DIST"] * dist_cm

F_img = img.sum() * sr_pix
F_ref = (Lb * np.exp(-tau_b)).sum() / (4.0 * np.pi * d_cm**2)
print(f"blocks in image file: {names}")
print(f"tau at band edges: {tau_b[0]:.3f} (13.6 eV)  {tau_b[-1]:.4f} (100 eV)")
print(f"F(direct image)  = {F_img:.5e} erg/s/cm^2")
print(f"F(reference)     = {F_ref:.5e} erg/s/cm^2")
print(f"ratio - 1        = {F_img/F_ref-1:+.4%}")
ny, nx = img.shape
ipk = np.unravel_index(np.argmax(img), img.shape)
print(f"peak pixel at {ipk} (center = ({nx//2}, {ny//2})), "
      f"peak fraction of total = {img.max()/img.sum():.4f}")
