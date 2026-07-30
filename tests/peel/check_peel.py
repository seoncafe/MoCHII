"""Peel-off gate (optically thin analytic case, peel_thin.in): a G0
single pass through a uniform sphere (nH = 100, x_HI = 1e-4 fixed, He
fully stripped, no dust) leaves only H I opacity, so the direct image
flux has the exact reference
    F = sum_b L_b e^{-tau_b} / (4 pi d^2),
    tau_b = nH x_HI sigma_HI(E_b) R        (tau ~ 0.7 at threshold)
with L_b the Planck bin luminosities on the stored bin grid.  The point
source must land in exactly the central pixel and the scattered image
must be identically zero.

sigma_HI comes from tools/python/mochii_rates.py, which mirrors
src/photo_xsec.f90 -- the exact hydrogenic expression, not the VFKY96 fit this
gate used to carry.  The fit is 0.775% high at the H I threshold, which is
where tau is largest here, so the reference was not the attenuation the run
transported.  tests/rates/check_python_rates.py holds the module to the
Fortran.

Gate record (2026-07-30, the stored 4e7-packet run): +0.0195%.  The earlier
record of +0.017% predates the stored output; with the reference carrying the
fitted H I cross section this same run read +0.1507%, nearly all of it the fit's
0.775% error at the threshold where tau is largest.

NOTE the converged-Stromgren variant (peel_free.in) is NOT a
quantitative gate: the neutral shell absorbs every ionizing photon
(tau = 57-3100 across the band), so image and reference are both
exponentially tiny and the comparison amplifies the reference's radial
discretization (physically correct — an embedded H II region emits no
ionizing photons — but ill-conditioned as a test).

Run:  python3 check_peel.py [peel_thin]
"""

import os
import sys
import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "tools", "python"))
from mochii_rates import sigma_HI

BASE = sys.argv[1] if len(sys.argv) > 1 else "peel_thin"

NH, XHI, RSPH = 100.0, 1.0e-4, 4.0          # sphere of the test grid
LBAND, TSTAR = 3.177837e38, 4.0e4


with h5py.File(f"{BASE}_rates.h5", "r") as f:
    Eb = f["E_bin/data"][:]
    dEb = f["dE_bin/data"][:]
    dist_cm = float(np.asarray(f["Gamma_HI"].attrs["DIST_CM"]).ravel()[0])
with h5py.File(f"{BASE}_image.h5", "r") as f:
    img = f["direct_ion/data"][:]
    sca = f["scatt_ion/data"][:]
    at = {k: float(np.asarray(v).ravel()[0]) for k, v in
          f["direct_ion"].attrs.items() if k in ("DXIM", "DIST")}
sr_pix = np.deg2rad(at["DXIM"])**2
d_cm = at["DIST"] * dist_cm

tau = NH * XHI * sigma_HI(Eb) * RSPH * dist_cm


def planck_w(E, T):
    x = E * 1.602176634e-12 / (1.380649e-16 * T)
    return E**3 / np.expm1(x)


Lb = planck_w(Eb, TSTAR) * dEb
Lb = Lb / Lb.sum() * LBAND
F_ref = (Lb * np.exp(-tau)).sum() / (4.0 * np.pi * d_cm**2)
F_img = img.sum() * sr_pix

print(f"tau(13.6 eV) = {tau[0]:.4f},  tau(100 eV) = {tau[-1]:.5f}")
print(f"F(direct image) = {F_img:.6e} erg/s/cm^2")
print(f"F(analytic)     = {F_ref:.6e} erg/s/cm^2")
print(f"ratio - 1       = {F_img/F_ref-1:+.4%}   (gate: |.| < 0.5%)")
ipk = np.unravel_index(np.argmax(img), img.shape)
nx = img.shape[0]
print(f"peak pixel {ipk} (center ({nx//2}, {nx//2})), "
      f"peak fraction = {img.max()/img.sum():.6f}")
print(f"scattered image total = {sca.sum():.3e} (must be 0 without dust)")
ok = abs(F_img/F_ref - 1) < 0.005 and img.max()/img.sum() > 0.999 \
     and sca.sum() == 0.0
print("GATE:", "PASS" if ok else "FAIL")
#--- The gate must fail the PROCESS, not just print: a CI runner (and a
#--- human skimming a log) otherwise reads a FAIL as a successful run.
sys.exit(0 if ok else 1)

