#!/usr/bin/env python3
"""External-field peel-off gate.

An isotropic external radiation field of band-integrated mean intensity
J = ext_intensity enters the box faces and is Lambert-emitted (cosine-weighted
about the inward normal).  The DIRECT peel toward a +z observer weights each
entry point by max(k_obs . n_in, 0)/pi, so only the far (-z) face illuminates
the observer.  The unattenuated image is therefore the background J filling the
box's projected solid angle Omega = A_proj / d^2:

    F(direct, unattenuated) = J * Omega = J * A_proj / d^2,
    A_proj = (2*xmax)*(2*ymax)          (the +z observer sees the z-face)

which is the tau-independent analytic gate.  Attenuation lowers the tau-walked
`direct` image by ~exp(-tau).

Two runs:
  peel_ext_thin (tau_ion ~ 7e-3): the band-integrated `direct_ion` image
      integrates to J*A_proj/d^2 to within the small (~0.3-1%) attenuation;
      scattered image is identically 0 (no dust).
  peel_ext_tau  (chord tau ~ 2, peel_bins + save_direc0): (a) the
      unattenuated `direc0_cube` still integrates to J*A_proj/d^2 (tau-free);
      (b) the central-pixel bin-by-bin ratio direct_cube/direc0_cube = exp(-tau)
      recovers the chord optical depth tau(E) = nH*sigma_HI(E)*L, L = 2*zmax,
      whose wavelength dependence must follow the H I cross section.

Run:  python3 check_peel_ext.py [peel_ext_thin | peel_ext_tau | both]
"""
import os
import sys
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
PC2CM = 3.0856776e18
EV2ERG = 1.602176634e-12
KB = 1.380649e-16
SIG0_HI = 6.30e-18          # exact hydrogenic threshold cross section [cm^2]
ETH_HI = 13.598             # eV


def parse_in(fname):
    vals = {}
    with open(fname) as fh:
        for ln in fh:
            s = ln.split("!", 1)[0].strip()
            if not s.lower().startswith("par%") or "=" not in s:
                continue
            key, rhs = s.split("=", 1)
            key = key.strip().lower().replace("par%", "")
            try:
                vals[key] = float(rhs.strip().split()[0])
            except ValueError:
                vals[key] = rhs.strip().split()[0].strip("'\"")
    return vals


def sig_HI_exact(E):
    """Exact ground-state hydrogenic H I cross section (Osterbrock & Ferland
    2006 Eq. 2.4), the same formula the code now uses."""
    E = np.asarray(E, float)
    out = np.zeros_like(E)
    m = E >= ETH_HI
    eps = np.sqrt(np.maximum(E[m] / ETH_HI - 1.0, 0.0))
    eps = np.where(eps < 1e-30, 1e-30, eps)
    num = np.exp(4.0 - 4.0 * np.arctan(eps) / eps)
    den = 1.0 - np.exp(-2.0 * np.pi / eps)
    out[m] = SIG0_HI * (ETH_HI / E[m])**4 * num / den
    return out


def sig_HI_pow(E):
    """Simple E^-3 power law near threshold."""
    E = np.asarray(E, float)
    return np.where(E >= ETH_HI, SIG0_HI * (E / ETH_HI)**(-3.0), 0.0)


def load_rates(base):
    with h5py.File(f"{base}_rates.h5", "r") as f:
        Eb = f["E_bin/data"][:].ravel()
        dEb = f["dE_bin/data"][:].ravel()
        dist_cm = float(np.asarray(f["Gamma_HI"].attrs["DIST_CM"]).ravel()[0])
    return Eb, dEb, dist_cm


def block_attrs(f, name):
    return {k: float(np.asarray(v).ravel()[0])
            for k, v in f[name].attrs.items() if k in ("DXIM", "DIST")}


def orient_cube(a, nbin):
    """Return the cube with the bin axis first: (nbin, s1, s2)."""
    ax = [i for i, n in enumerate(a.shape) if n == nbin]
    if not ax:
        raise ValueError(f"no cube axis of length {nbin} in shape {a.shape}")
    return np.moveaxis(a, ax[0], 0)


def gate_thin(base="peel_ext_thin"):
    inf = parse_in(os.path.join(HERE, f"{base}.in"))
    J = inf["ext_intensity"]
    xmax, ymax = inf["xmax"], inf["ymax"]
    Eb, dEb, dist_cm = load_rates(os.path.join(HERE, base))
    with h5py.File(os.path.join(HERE, f"{base}_image.h5"), "r") as f:
        img = f["direct_ion/data"][:]
        sca = f["scatt_ion/data"][:]
        at = block_attrs(f, "direct_ion")
    sr_pix = np.deg2rad(at["DXIM"])**2
    d_cm = at["DIST"] * dist_cm
    A_proj = (2.0 * xmax) * (2.0 * ymax) * dist_cm**2      # z-face, cm^2
    F_ref = J * A_proj / d_cm**2
    F_img = img.sum() * sr_pix
    print(f"\n=== {base}  (optically thin analytic gate) ===")
    print(f"  J (ext_intensity)          = {J:.4e} erg/s/cm^2/sr")
    print(f"  A_proj = (2xmax)(2ymax)    = {A_proj:.4e} cm^2")
    print(f"  d (observer)               = {d_cm:.4e} cm")
    print(f"  F(direct image)            = {F_img:.6e} erg/s/cm^2")
    print(f"  F(analytic J*A_proj/d^2)   = {F_ref:.6e} erg/s/cm^2")
    print(f"  ratio - 1                  = {F_img/F_ref-1:+.4%}  "
          f"(attenuation lowers it ~0.3-1%; gate |.| < 2%)")
    print(f"  scattered image total      = {sca.sum():.3e} (must be 0)")
    ok = abs(F_img / F_ref - 1) < 0.02 and sca.sum() == 0.0
    print(f"  --> {'PASS' if ok else 'FAIL'}")
    return ok


def gate_tau(base="peel_ext_tau", nblk=10, nlow=6):
    inf = parse_in(os.path.join(HERE, f"{base}.in"))
    J = inf["ext_intensity"]
    xmax, ymax, zmax = inf["xmax"], inf["ymax"], inf["zmax"]
    nH = inf["nh_const"]      # parse_in lowercases all keys
    eion_min = inf.get("eion_min", 13.598)
    Eb, dEb, dist_cm = load_rates(os.path.join(HERE, base))
    nbin = Eb.size
    with h5py.File(os.path.join(HERE, f"{base}_image.h5"), "r") as f:
        dcube = orient_cube(f["direct_cube/data"][:], nbin)
        d0cube = orient_cube(f["direc0_cube/data"][:], nbin)
        scube = orient_cube(f["scatt_cube/data"][:], nbin)
        at = block_attrs(f, "direct_cube")
    sr_pix = np.deg2rad(at["DXIM"])**2
    d_cm = at["DIST"] * dist_cm
    L = 2.0 * zmax * dist_cm                              # central chord [cm]
    A_proj = (2.0 * xmax) * (2.0 * ymax) * dist_cm**2
    F_ref = J * A_proj / d_cm**2

    print(f"\n=== {base}  (attenuation morphology gate) ===")
    # (a) unattenuated direc0 flux gate (tau-independent)
    F0 = d0cube.sum() * sr_pix
    print(f"  (a) unattenuated flux gate:")
    print(f"      F(direc0 cube)         = {F0:.6e} erg/s/cm^2")
    print(f"      F(J*A_proj/d^2)        = {F_ref:.6e} erg/s/cm^2")
    print(f"      ratio - 1              = {F0/F_ref-1:+.4%}  (gate |.| < 2%)")
    ok_a = abs(F0 / F_ref - 1) < 0.02

    # (b) central-block chord tau spectrum
    s1, s2 = dcube.shape[1], dcube.shape[2]
    c1, c2 = s1 // 2, s2 // 2
    sl1 = slice(c1 - nblk, c1 + nblk + 1)
    sl2 = slice(c2 - nblk, c2 + nblk + 1)
    db = dcube[:, sl1, sl2].sum(axis=(1, 2))
    d0 = d0cube[:, sl1, sl2].sum(axis=(1, 2))
    ion = (Eb >= eion_min) & (d0 > 0) & (db > 0)
    idx = np.where(ion)[0][:nlow]
    tau_meas = -np.log(db[idx] / d0[idx])
    E = Eb[idx]
    tau_exact = nH * sig_HI_exact(E) * L
    tau_pow = nH * sig_HI_pow(E) * L
    print(f"  (b) central-chord tau (block {2*nblk+1}x{2*nblk+1}, L = "
          f"{L/PC2CM:.3f} pc):")
    print(f"      {'E[eV]':>8} {'tau_meas':>9} {'tau_exact':>10} "
          f"{'tau_E^-3':>9} {'meas/exact':>11}")
    for e, tm, te, tp in zip(E, tau_meas, tau_exact, tau_pow):
        print(f"      {e:8.3f} {tm:9.4f} {te:10.4f} {tp:9.4f} "
              f"{tm/te:11.4f}")
    # first-bin absolute check vs nH*sig0*L and vs exact
    tau1_sig0 = nH * SIG0_HI * L
    print(f"      tau(bin1) measured     = {tau_meas[0]:.4f}")
    print(f"      nH*6.3e-18*L (sigma_0) = {tau1_sig0:.4f}  "
          f"(ratio {tau_meas[0]/tau1_sig0:.4f}, expect <1: center>threshold)")
    print(f"      nH*sigma_exact(E1)*L   = {tau_exact[0]:.4f}  "
          f"(ratio {tau_meas[0]/tau_exact[0]:.4f}; gate within 10%)")
    # wavelength dependence: measured ratio vs exact-hydrogenic ratio
    r_meas = tau_meas / tau_meas[0]
    r_exact = tau_exact / tau_exact[0]
    r_pow = tau_pow / tau_pow[0]
    dev_exact = np.abs(r_meas / r_exact - 1.0)
    dev_pow = np.abs(r_meas / r_pow - 1.0)
    print(f"      tau(E)/tau(E1) vs exact-hydrogenic: max dev "
          f"{dev_exact.max():.2%} (gate < 5%)")
    print(f"      tau(E)/tau(E1) vs E^-3 law:         max dev "
          f"{dev_pow.max():.2%}")
    ok_b1 = abs(tau_meas[0] / tau_exact[0] - 1.0) < 0.10
    ok_b2 = dev_exact.max() < 0.05
    ok_b3 = scube.sum() == 0.0
    print(f"      scattered cube total   = {scube.sum():.3e} (must be 0)")
    ok = ok_a and ok_b1 and ok_b2 and ok_b3
    print(f"  --> {'PASS' if ok else 'FAIL'} "
          f"(a={ok_a}, tau1={ok_b1}, spectrum={ok_b2}, scatt0={ok_b3})")
    return ok


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    passed = True
    if which in ("both", "peel_ext_thin"):
        passed &= gate_thin()
    if which in ("both", "peel_ext_tau"):
        passed &= gate_tau()
    print(f"\nOVERALL: {'PASS' if passed else 'FAIL'}")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
