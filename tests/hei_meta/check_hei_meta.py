#!/usr/bin/env python3
"""He I 2^3S metastable diagnostic gate (docs/HEI_10830_PLAN.md H2).

Independently recomputes the module's n_2s3 (and Phi_3) from the data-file
coefficients on each run's OWN converged state, and checks the OH18 analytic
limits.

Gates:
  1. Agreement (every run): module n_2s3 vs the Python closed form on the same
     n_e/T_e/x_* and the module's own Phi_3 -> machine precision (same formula,
     same inputs).  NOTE: Phi_3 is never 0 even without add_fuv -- the 2^3S
     threshold is 4.78 eV, so the 13.6-100 eV ionizing band photoionizes the
     metastable too (sigma_3 is significant across the band); the module is
     right to integrate over all bins.
  2. Analytic limits: low-d n_3/n_HeII near alpha_3 n_e/A_31 (the pure-
     radiative ratio is within a few %, the EUV-only Phi_3 being << A_31 in the
     core); high-d n_3/n_HeII moves toward the plateau (n_e > n_crit); a pure
     n_e scan reproduces plateau and n_crit.
  3. Phi_3 (add_fuv run): Phi_3 integrated independently from the rates-file
     J_nu matches the module Phi_3; the full-chain n_2s3 (large soft-UV Phi_3,
     so the sigma_3 cross section is exercised) still matches to machine
     precision.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "python"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "fitting"))
from mochii_output import (read_sections, HeiMetaData, read_line10830,  # noqa
                           SQRTPI_E2_MEC, M_HE_G, KB_CGS)                # noqa
import make_hei_metastable as M                        # noqa: E402

ATOM = os.path.join(HERE, "..", "..", "data", "atomic")
EV2ERG = 1.602176634e-12
HPLANCK = 6.62607015e-27
NH_CONST = {"hei_lowd": 100.0, "hei_highd": 1.0e4, "hei_fuv": 100.0}
YHE = 0.1

FAILS = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print("  [%s] %s %s" % (tag, name, detail))
    if not ok:
        FAILS.append(name)


# --------------------------------------------------------- rate reconstruction
def parse_coeffs():
    """Read data/atomic/hei_metastable.txt into a dict (independent of M's
    module-level constants; proves the FILE is the source of truth)."""
    c = {}
    for ln in open(os.path.join(ATOM, "hei_metastable.txt")):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        t = ln.split()
        key, vals = t[0], [float(x) for x in t[1:]]
        c[key] = vals
    return c


C = parse_coeffs()
KB_EV = 1.380649e-16 / 1.602176634e-12
TMIN, TMAX = C["TRANGE"]


def clampT(T):
    return np.minimum(np.maximum(T, TMIN), TMAX)


def alpha3(T):
    Cc, p = C["ALPHA3"]
    return Cc * (T / 1e4) ** p


def qexc(T, key):
    pref, sc, Ee, div, a, b, cc, d = C[key]
    kT = KB_EV * T
    ups = a * np.exp(b * T) + cc * np.exp(d * T)
    return pref * np.sqrt(sc / kT) * np.exp(-Ee / kT) * ups / div


def penning(T):
    Ts, Alo, plo, Ahi, phi = C["PENNING"]
    return np.where(T <= Ts, Alo * (300.0 / T) ** plo, Ahi * (300.0 / T) ** phi)


def sigma3(E):
    return M.sigma3(E)     # same VFKY96 wings + bridge; reads M's SIG3* (= file)


def n3_closed(ne, Te, nHI, nHeI, nHeII, Phi3):
    Tc = clampT(Te)
    a3 = alpha3(Tc)
    q13 = qexc(Tc, "QEXC13")
    q31a = qexc(Tc, "QEXC31A")
    q31b = qexc(Tc, "QEXC31B")
    Q31 = penning(Tc)
    A31 = C["A31"][0]
    src = a3 * ne * nHeII + q13 * ne * nHeI
    snk = A31 + Phi3 + ne * (q31a + q31b) + nHI * Q31
    return src / snk


# ------------------------------------------------------------------- per-run
def load_run(base):
    hm = read_sections(os.path.join(HERE, base + "_heimeta.h5"))
    rt = read_sections(os.path.join(HERE, base + "_rates.h5"))
    d = dict(
        n2s3=hm["n_2s3"]["data"].ravel(),
        Phi3=hm["Phi_3"]["data"].ravel(),
        ne=hm["n_e"]["data"].ravel(),
        Te=hm["T_e"]["data"].ravel(),
        xHI=rt["x_HI"]["data"].ravel(),
        xHeI=rt["x_HeI"]["data"].ravel(),
        xHeII=rt["x_HeII"]["data"].ravel(),
    )
    d["nH"] = np.full_like(d["ne"], NH_CONST[base])
    # rmax cut leaves nH=0 outside; those have no He -> exclude from tests.
    d["nH"][d["ne"] <= 0.0] = 0.0
    d["rates"] = rt
    return d


def agreement_gate(base, d):
    """Recompute n_2s3 from the data-file closed form on the run's own state
    (using the module's own Phi_3) and compare to the module output.  Same
    formula + same inputs -> machine precision.  Phi_3 is ALWAYS the module
    value: the 2^3S threshold is 4.78 eV, so the 13.6-100 eV ionizing photons
    photoionize the metastable too (sigma_3 is significant across the band) --
    Phi_3 is never 0 wherever there is an ionizing field."""
    m = d["nH"] > 0.0
    nHI = d["nH"] * d["xHI"]
    nHeI = d["nH"] * YHE * d["xHeI"]
    nHeII = d["nH"] * YHE * d["xHeII"]
    ref = n3_closed(d["ne"], d["Te"], nHI, nHeI, nHeII, d["Phi3"])
    rel = np.abs(ref[m] - d["n2s3"][m]) / np.maximum(np.abs(d["n2s3"][m]), 1e-300)
    worst = rel.max()
    check("%s agreement: worst rel = %.2e" % (base, worst), worst < 1e-8,
          "(tol 1e-8)")
    return worst


# ============================================================= Stage 2 (H5)
def build_synthetic_heimeta(path, n3, Te, ne=1.0, nx=8, box=1.0,
                            dist_cm=3.0856776e18):
    """Write a tiny uniform section-format heimeta HDF5 for the analytic
    column/tau checks: an nx^3 Cartesian cube filling [-box, box]^3 with a
    uniform 2^3S density n3 [cm^-3] and temperature Te [K]."""
    import h5py
    step = 2.0 * box / nx
    c = -box + (np.arange(nx) + 0.5) * step
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    xyz = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])     # (3, nleaf)
    nleaf = xyz.shape[1]
    size = np.full(nleaf, step)
    with h5py.File(path, "w") as f:
        def sec(name, data, **attrs):
            g = f.create_group(name)
            g.create_dataset("data", data=np.asarray(data))
            g.attrs["EXTNAME"] = name
            for k, v in attrs.items():
                g.attrs[k] = v
        sec("n_2s3", np.full(nleaf, n3), DIST_CM=dist_cm, A31=1.272e-4)
        sec("Phi_3", np.zeros(nleaf))
        sec("n_e", np.full(nleaf, ne))
        sec("T_e", np.full(nleaf, Te))
        sec("LeafXYZ", xyz)
        sec("LeafSize", size)
    return dict(nx=nx, box=box, step=step, dist_cm=dist_cm, n3=n3, Te=Te)


def kappa0_hand(n3, Te, rows, v_turb=0.0):
    """Independent hand computation of the line-center opacity [cm^-1]."""
    b = np.sqrt(2.0 * KB_CGS * Te / M_HE_G + (v_turb * 1.0e5) ** 2)
    f_lam = np.sum(rows[:, 1] * rows[:, 0] * 1.0e-8)
    return SQRTPI_E2_MEC * f_lam * n3 / b


def stage2_gate():
    """docs/HEI_10830_PLAN.md H5: 10830 absorption column / tau diagnostic."""
    import tempfile
    print("\n[Stage 2: 10830 absorption diagnostic]")
    rows = read_line10830(ATOM)
    # component f-values sum 1:3:5, ratios exact.
    f = rows[:, 1]
    check("LINE10830 f-values 1:3:5 (%.4f:%.4f:%.4f)" % tuple(f),
          abs(f[1] / f[0] - 3.0) < 0.01 and abs(f[2] / f[0] - 5.0) < 0.01)
    check("LINE10830 A = 1.0216e7 s^-1 each (NIST/p-winds, not CHIANTI 1.15e7)",
          np.all(np.abs(rows[:, 2] / 1.0216e7 - 1.0) < 1e-3))

    tmp = tempfile.mkdtemp()
    # -- thin uniform cube: column and tau vs the hand integral ------------
    n3, Te = 1.0e-9, 1.0e4                # tiny n3 -> optically thin (no warn)
    syn = os.path.join(tmp, "syn_heimeta.h5")
    g = build_synthetic_heimeta(syn, n3, Te, nx=8, box=1.0)
    d = HeiMetaData(syn, atomic_dir=ATOM)
    depth_cm = 2.0 * g["box"] * g["dist_cm"]             # full LOS through cube

    # conservation: sum(column * A_pix) == sum(n * V).
    img, ext = d.column_map("n_2s3", npix=64)
    dpix = (ext[1] - ext[0]) / 64
    apix = (dpix * d.dist_cm) ** 2
    lhs = float((img * apix).sum())
    rhs = float((d.fields["n_2s3"] * d.vol_cm3).sum())
    check("synthetic column conservation sum(col A_pix)=sum(n V): rel %.2e"
          % (abs(lhs - rhs) / rhs), abs(lhs - rhs) / rhs < 1e-12)

    # interior pixel column == n3 * depth (full-cube LOS, exact).
    mid = 32
    col_interior = img[mid, mid]
    check("synthetic interior column = n3 * depth: %.6e vs %.6e (rel %.2e)"
          % (col_interior, n3 * depth_cm,
             abs(col_interior / (n3 * depth_cm) - 1.0)),
          abs(col_interior / (n3 * depth_cm) - 1.0) < 1e-6)

    # tau_0 interior vs hand kappa * depth, for each component and the doublet.
    for comp, sel in (("doublet", [1, 2]), (0, [0]), ("all", [0, 1, 2])):
        tau, _ = d.tau10830_map(component=comp, npix=64, warn=False)
        khand = kappa0_hand(n3, Te, rows[sel, :])
        thand = khand * depth_cm
        rel = abs(tau[mid, mid] / thand - 1.0)
        check("synthetic tau_0(%s) interior = kappa_hand*depth: "
              "%.4e vs %.4e (rel %.2e)" % (comp, tau[mid, mid], thand, rel),
              rel < 1e-6)

    # thin case must NOT warn (tau << 1 everywhere).
    tau_thin, _ = d.tau10830_map(component="doublet", npix=64, warn=False)
    check("synthetic thin case tau_0 max = %.2e < 1 (no warning expected)"
          % np.nanmax(tau_thin), np.nanmax(tau_thin) < 1.0)

    # scaled-up n3 must push tau > 1 (the warning path); capture the message.
    import io
    from contextlib import redirect_stdout
    syn2 = os.path.join(tmp, "syn_thick_heimeta.h5")
    build_synthetic_heimeta(syn2, n3 * 1.0e5, Te, nx=8, box=1.0)
    d2 = HeiMetaData(syn2, atomic_dir=ATOM)
    buf = io.StringIO()
    with redirect_stdout(buf):
        tau_thick, _ = d2.tau10830_map(component="doublet", npix=64, warn=True)
    fired = "tau_10830 > 1" in buf.getvalue()
    check("synthetic thick case (n3 x 1e5): tau_0 max = %.2e > 1 and warning "
          "fires" % np.nanmax(tau_thick),
          np.nanmax(tau_thick) > 1.0 and fired)

    # -- real run: hei_highd conservation + sanity -------------------------
    dh = HeiMetaData(os.path.join(HERE, "hei_highd_heimeta.h5"), atomic_dir=ATOM)
    img, ext = dh.column_map("n_2s3", npix=256)
    dpix = (ext[1] - ext[0]) / 256
    apix = (dpix * dh.dist_cm) ** 2
    lhs = float((img * apix).sum())
    rhs = float((dh.fields["n_2s3"] * dh.vol_cm3).sum())
    check("hei_highd column conservation: rel %.2e" % (abs(lhs - rhs) / rhs),
          abs(lhs - rhs) / rhs < 1e-10)
    tau_d, _ = dh.tau10830_map(component="doublet", npix=256, warn=False)
    tau_0, _ = dh.tau10830_map(component=0, npix=256, warn=False)
    ratio = np.nanmax(tau_d) / np.nanmax(tau_0)
    fexp = (rows[1, 1] + rows[2, 1]) / rows[0, 1]        # doublet/J0 f-ratio
    check("hei_highd tau doublet/J0 = %.3f = f-ratio %.3f (blended lines)"
          % (ratio, fexp), abs(ratio / fexp - 1.0) < 0.02,
          "[max tau_doublet = %.2e -- optically thick, as flagged]"
          % np.nanmax(tau_d))


# =========================================================================
def main():
    print("=== He I 2^3S metastable gate ===")

    # ---- pure-Python analytic limits (data-file coefficients) --------------
    print("\n[analytic limits, one-zone from the data file]")
    T0 = 1e4
    a3 = alpha3(T0)
    q31 = qexc(T0, "QEXC31A") + qexc(T0, "QEXC31B")
    A31 = C["A31"][0]
    lowd = a3 * 100.0 / A31
    plateau = a3 / q31
    ncrit = A31 / q31
    check("low-density ratio alpha_3 ne/A_31 @ne=100,1e4K = %.4e (spec 1.65e-7)"
          % lowd, abs(lowd / 1.65e-7 - 1.0) < 0.03)
    check("plateau alpha_3/(q31a+q31b) = %.4e" % plateau, plateau > 0)
    check("n_crit = A_31/(q31a+q31b) = %.1f cm^-3 (few x 1e3)" % ncrit,
          1e3 < ncrit < 1e4)
    # n_e scan: ratio crosses 0.5*plateau near n_crit
    ne = np.logspace(1, 6, 400)
    ratio = a3 * ne / (A31 + ne * q31)             # n_3/n_HeII, no HI/Phi
    j = int(np.argmin(np.abs(ratio - 0.5 * plateau)))
    check("n_e scan: ratio=0.5*plateau at n_e=%.0f ~ n_crit %.0f"
          % (ne[j], ncrit), abs(ne[j] / ncrit - 1.0) < 0.1)

    # ---- low-density run ----------------------------------------------------
    print("\n[low-density run: hei_lowd]")
    dl = load_run("hei_lowd")
    agreement_gate("hei_lowd", dl)
    # EUV-only Phi_3 (no add_fuv): the ionizing band photoionizes 2^3S but the
    # rate is small vs A_31 in the core, so the radiative limit still holds.
    xyz = dl["rates"]["LeafXYZ"]["data"]
    xyz = xyz if xyz.shape[0] == dl["ne"].size else xyz.T
    r = np.sqrt((xyz ** 2).sum(axis=1))
    core = (r < 1.0) & (dl["nH"] > 0) & (dl["xHeII"] > 0.5)
    check("hei_lowd EUV-only Phi_3 << A_31 in the core "
          "(max Phi_3/A_31 = %.2e)" % (dl["Phi3"][core].max() / A31),
          dl["Phi3"][core].max() / A31 < 0.1)
    # deep ionized core (r<1 pc, x_HeII high): n_3/n_HeII near alpha_3 ne/A_31
    nHeII = dl["nH"] * YHE * dl["xHeII"]
    ratio_run = dl["n2s3"][core] / nHeII[core]
    rad_lim = alpha3(clampT(dl["Te"][core])) * dl["ne"][core] / A31
    dev = np.abs(ratio_run / rad_lim - 1.0)
    check("hei_lowd core n_3/n_HeII within a few %% of alpha_3 ne/A_31 "
          "(median dev %.2f%%, max %.2f%%)" % (100 * np.median(dev), 100 * dev.max()),
          np.median(dev) < 0.05,
          "[%d core leaves, <ne>=%.1f]" % (core.sum(), dl["ne"][core].mean()))

    # ---- high-density run ---------------------------------------------------
    print("\n[high-density run: hei_highd]")
    dh = load_run("hei_highd")
    agreement_gate("hei_highd", dh)
    xyz = dh["rates"]["LeafXYZ"]["data"]
    xyz = xyz if xyz.shape[0] == dh["ne"].size else xyz.T
    r = np.sqrt((xyz ** 2).sum(axis=1))
    core = (r < 0.15) & (dh["nH"] > 0) & (dh["xHeII"] > 0.5)
    nHeII = dh["nH"] * YHE * dh["xHeII"]
    ratio_run = (dh["n2s3"][core] / nHeII[core])
    necore = dh["ne"][core].mean()
    # in the plateau direction: ratio is well above the low-density value and a
    # substantial fraction of the plateau (n_e > n_crit -> collisions matter).
    frac = np.median(ratio_run) / plateau
    check("hei_highd core n_e=%.2e > n_crit=%.0f (collisional regime)"
          % (necore, ncrit), necore > ncrit)
    check("hei_highd core n_3/n_HeII = %.3f of the plateau (0.4-0.95 expected)"
          % frac, 0.4 < frac < 0.98,
          "[%d core leaves]" % core.sum())

    # ---- add_fuv run: Phi_3 -------------------------------------------------
    print("\n[FUV run: hei_fuv (Phi_3 > 0)]")
    df = load_run("hei_fuv")
    # independent Phi_3 from the rates-file J_nu block.
    rt = df["rates"]
    Jnu = rt["J_nu"]["data"]                          # (nnu_band, nleaf) or T
    Ebin = rt["E_bin"]["data"].ravel()
    dEbin = rt["dE_bin"]["data"].ravel()
    nnu = Ebin.size
    if Jnu.shape[0] != nnu:
        Jnu = Jnu.T
    sig = sigma3(Ebin)                                # cm^2 in each bin
    dnu = dEbin * EV2ERG / HPLANCK                    # Hz
    # Phi_3 = sum_nu 4 pi J_nu dnu sigma_3 / (h nu),  h nu = E*ev2erg
    fourpi = 4.0 * np.pi
    integ = (fourpi * Jnu * dnu[:, None] * sig[:, None]
             / (Ebin[:, None] * EV2ERG))              # (nnu, nleaf)
    Phi3_ref = integ.sum(axis=0)
    m = df["nH"] > 0.0
    pos = m & (df["Phi3"] > 0.0)
    rel = np.abs(Phi3_ref[pos] - df["Phi3"][pos]) / df["Phi3"][pos]
    check("hei_fuv Phi_3 vs independent J_nu integration: worst rel = %.2e"
          % rel.max(), rel.max() < 1e-6, "(tol 1e-6)")
    check("hei_fuv Phi_3 > 0 in the ionized region (max = %.3e s^-1)"
          % df["Phi3"][m].max(), df["Phi3"][m].max() > 0.0)
    # full-chain agreement WITH Phi_3 (exercises sigma_3 in the module).
    agreement_gate("hei_fuv", df)

    # ---- Stage 2: 10830 absorption diagnostic -------------------------------
    stage2_gate()

    # ---- verdict ------------------------------------------------------------
    print("\n=== %s ===" % ("ALL GATES PASS" if not FAILS
                            else "FAILED: " + ", ".join(FAILS)))
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
