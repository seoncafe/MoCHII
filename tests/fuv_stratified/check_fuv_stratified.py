#!/usr/bin/env python3
"""Gate for the band-stratified launch of the FUV extension (par%fuv_sampling).

  check_fuv_stratified.py make <template.in> <workdir> <repo>
      write the inputs run_check.sh runs (base.in and the fixed-state set)
  check_fuv_stratified.py check <workdir>
      evaluate gates (a) and (b) on the outputs

(a) Launched luminosity: in every 'separate' run the ionizing stratum launches
    L(eion_min-eion_max), the FUV stratum L(efuv_min-eion_min) and the two sum
    to the band luminosity, to rounding (the ' ION: launched L' lines); the
    ionizing part equals par%luminosity (the fixed quantity of a Planck source)
    and L_FUV/L_ion equals the ratio of the two Planck integrals computed here
    by quadrature, to the accuracy of the sampler's piecewise-linear CDF
    (1e-6; source_cdf_tol = 1e-8 is a local tolerance).
(b) One fixed state (restart_file, gas_niter = 0), transported with
    'proportional', 'separate' and add_fuv = .false. at four seeds each, for a
    40 kK and a 20 kK blackbody: the volume integrals of the photoionization
    rate n_HI Gamma_HI, the 2^3S photoionization rate Phi_3, the FUV field G0
    (from the FUV bins of J_nu) and the grain heating agree between the modes
    within their seed scatter; the leaf-level Monte Carlo noise of Gamma_HI in
    'separate' is that of the add_fuv = .false. run and 'proportional''s is
    larger by 1/sqrt(ionizing share); and with the Sobol launch 'separate''s
    Gamma_HI equals the add_fuv = .false. run's to rounding (same Sobol
    points, the FUV packets score no hydrogen ionization), with the served
    (master-slave) schedule giving the same numbers as the static one.  "To
    rounding" there means to the 1e-8 by which the band sampler's conditional
    CDF and a sampler built on [eion_min, eion_max] alone differ (two
    piecewise-linear fits of one Planck function); the gate allows 1e-6.
"""
import glob, os, re, sys
import numpy as np

H_PLANCK = 6.62607015e-27
EV2ERG = 1.602176634e-12
KB_EV = 8.617333262145e-5

SEEDS = [1234, 4321, 2468, 8642]
TSTARS = {"t40": 4.0e4, "t20": 2.0e4}
MODES = {
    "prop": ["add_fuv=.true.", "fuv_sampling='proportional'"],
    "sep":  ["add_fuv=.true.", "fuv_sampling='separate'"],
    "off":  ["add_fuv=.false."],
}


def edit_namelist(template, out, repo, items):
    """copy template, replacing or adding par%key = value items"""
    over = dict(kv.split("=", 1) for kv in items)
    lines, seen = [], set()
    for ln in open(template).read().split("\n"):
        m = re.match(r"\s*par%(\w+)\s*=", ln)
        if m and m.group(1) in over:
            seen.add(m.group(1))
            lines.append(f" par%{m.group(1)} = {over[m.group(1)]}")
            continue
        lines.append(ln.replace("'../../data/", f"'{repo}/data/"))
    end = max(i for i, l in enumerate(lines) if l.strip() == "/")
    for k, v in over.items():
        if k not in seen:
            lines.insert(end, f" par%{k} = {v}")
            end += 1
    open(out, "w").write("\n".join(lines))


def make(template, work, repo):
    os.makedirs(work, exist_ok=True)
    edit_namelist(template, f"{work}/base.in", repo, ["out_file='base.h5'"])
    names = []
    fixed = ["gas_niter=0", f"restart_file='{work}/base_rates.h5'"]
    for tk, T in TSTARS.items():
        for mode, items in MODES.items():
            for s in SEEDS:
                name = f"{tk}_{mode}_s{s}"
                edit_namelist(template, f"{work}/{name}.in", repo,
                              fixed + items + [f"tstar={T}", f"iseed={s}",
                                               f"out_file='{name}.h5'"])
                names.append(name)
    # Sobol: the ionizing stratum is the add_fuv = .false. launch; the served
    # schedule transports the same packets as the static one.  Dust scattering
    # is off here: its draws stay on the Mersenne Twister of each rank, so
    # with it on a run is fixed by its launch only up to the scattering.
    for name, items in [("sobol_sep", MODES["sep"]), ("sobol_off", MODES["off"]),
                        ("sobol_sep_served", MODES["sep"] + ["use_master_slave=.true.",
                                                            "num_send_at_once=5000"])]:
        edit_namelist(template, f"{work}/{name}.in", repo,
                      fixed + items + ["launch_sequence='sobol'", "no_photons=1048576",
                                       "no_photons_fuv=262144", "ion_dust_scatter=.false.",
                                       f"out_file='{name}.h5'"])
        names.append(name)
    open(f"{work}/runs.txt", "w").write("\n".join(names) + "\n")


def read_run(work, name):
    import h5py
    f = h5py.File(f"{work}/{name}_rates.h5", "r")
    d = {}
    def ds(key):
        obj = f[key]
        if isinstance(obj, h5py.Dataset):
            return np.asarray(obj[()])
        for sub in obj.values():
            if isinstance(sub, h5py.Dataset):
                return np.asarray(sub[()])
        raise KeyError(key)
    size = ds("LeafSize").ravel()
    vol = size**3
    xHI = ds("x_HI").ravel()
    d["vol"] = vol
    d["xHI"] = xHI
    d["ne"] = ds("n_e").ravel()
    d["G"] = ds("Gamma_HI").ravel()
    d["Hd"] = ds("Heat_dust").ravel()
    jnu = ds("J_nu")
    E = ds("E_bin").ravel()
    dE = ds("dE_bin").ravel()
    if jnu.shape[0] != E.size:
        jnu = jnu.T
    fuv = E < 13.598
    dnu = dE*EV2ERG/H_PLANCK
    d["G0"] = 4*np.pi*(jnu[fuv, :]*dnu[fuv, None]).sum(axis=0)/1.6e-3
    hm = h5py.File(f"{work}/{name}_heimeta.h5", "r")
    obj = hm["Phi_3"]
    d["Phi3"] = np.asarray(obj[()] if isinstance(obj, h5py.Dataset)
                           else list(obj.values())[0][()]).ravel()
    return d


def launched(work, name):
    out = {}
    for ln in open(f"{work}/{name}.log"):
        m = re.search(r"launched L, (ionizing stratum|FUV stratum|both strata|whole band)"
                      r"\s*=\s*(\S+),.*=\s*(\S+)\s*$", ln)
        if m:
            out[m.group(1)] = (float(m.group(2)), float(m.group(3)))
        m = re.search(r"L\(eion_min-eion_max\) =\s*(\S+),\s*L\(efuv_min-eion_min\) =\s*(\S+),", ln)
        if m:
            out["Lion"], out["Lfuv"] = float(m.group(1)), float(m.group(2))
    return out


def planck_ratio(T, e1, e2, e3):
    from scipy.integrate import quad
    kT = KB_EV*T
    f = lambda e: e**3/np.expm1(e/kT)
    lo = quad(f, e1, e2, epsabs=0, epsrel=1e-13, limit=500)[0]
    hi = quad(f, e2, e3, epsabs=0, epsrel=1e-13, limit=500)[0]
    return lo/hi


def check(work):
    ok = True
    print("(a) launched luminosity of the strata ('separate')")
    worst = 0.0
    for tk, T in TSTARS.items():
        for s in SEEDS:
            L = launched(work, f"{tk}_sep_s{s}")
            for k in ("ionizing stratum", "FUV stratum", "both strata"):
                worst = max(worst, abs(L[k][1]))
        L = launched(work, f"{tk}_sep_s{SEEDS[0]}")
        r_code = L["Lfuv"]/L["Lion"]
        r_quad = planck_ratio(T, 4.78, 13.598, 100.0)
        dev_lum = L["Lion"]/3.178e38 - 1
        dev_r = r_code/r_quad - 1
        print(f"    {tk}: L_ion = {L['Lion']:.12e} (par%luminosity x {1+dev_lum:.15f}), "
              f"L_FUV/L_ion = {r_code:.10f}, quadrature {r_quad:.10f}, rel {dev_r:+.2e}")
        # the Planck sampler is a piecewise-linear CDF refined to
        # source_cdf_tol = 1e-8 locally; its band integrals carry ~1e-7
        ok &= abs(dev_lum) < 1e-12 and abs(dev_r) < 1e-6
    print(f"    worst |launched/L - 1| over strata and seeds = {worst:.2e}")
    ok &= worst < 1e-12

    print("(b) field estimators on one state, four seeds per mode")
    qn = {"G": "n_HI Gamma_HI dV", "Phi3": "Phi_3 dV", "G0": "G0 dV", "Hd": "Heat_dust dV"}
    for tk in TSTARS:
        runs = {m: [read_run(work, f"{tk}_{m}_s{s}") for s in SEEDS] for m in MODES}
        base = runs["off"][0]
        vol, xHI = base["vol"], base["xHI"]
        ionized = (xHI < 0.1) & (base["ne"] > 0.0)      # ionized gas leaves
        share = None
        L = launched(work, f"{tk}_prop_s{SEEDS[0]}")
        share = L["Lion"]/(L["Lion"] + L["Lfuv"])
        print(f"  {tk}: ionizing share of the one CDF = {share:.4f}")
        for q, label in qn.items():
            vals = {}
            for m in MODES:
                if m == "off" and q in ("G0",):
                    continue
                w = xHI if q == "G" else 1.0
                I = np.array([np.sum(r[q]*w*vol) for r in runs[m]])
                vals[m] = (I.mean(), I.std(ddof=1)/np.sqrt(len(I)))
            s = "    " + f"{label:18s}"
            for m, (mu, se) in vals.items():
                s += f"  {m} {mu:.6e} +- {se/mu:.1e}"
            mp, sp = vals["prop"]; ms, ss = vals["sep"]
            z = (ms - mp)/np.hypot(sp, ss)
            s += f"  sep/prop-1 = {ms/mp-1:+.2e} ({z:+.1f} sigma)"
            print(s)
            ok &= abs(z) < 4.0
        # leaf-level noise of Gamma_HI (ionized leaves) and of G0, Phi_3
        for q in ("G", "Phi3", "G0"):
            noise = {}
            for m in MODES:
                if m == "off" and q == "G0":
                    continue
                A = np.array([r[q][ionized] for r in runs[m]])
                mu = A.mean(axis=0)
                good = mu > 0
                noise[m] = np.sqrt(np.mean((A[:, good].std(axis=0, ddof=1)/mu[good])**2))
            s = f"    leaf rms noise {q:5s}:" + "".join(f"  {m} {v:.3e}" for m, v in noise.items())
            if q == "G":
                s += (f"   prop/sep = {noise['prop']/noise['sep']:.3f} "
                      f"(1/sqrt(share) = {1/np.sqrt(share):.3f}), sep/off = "
                      f"{noise['sep']/noise['off']:.3f}")
                ok &= abs(noise["sep"]/noise["off"] - 1) < 0.1
                ok &= abs(noise["prop"]/noise["sep"]*np.sqrt(share) - 1) < 0.15
            print(s)

    print("(b') Sobol: the ionizing stratum is the add_fuv = .false. launch")
    a, b, c = (read_run(work, n) for n in ("sobol_sep", "sobol_off", "sobol_sep_served"))
    # ionized gas leaves: beyond the front Gamma_HI falls to 1e-100 and its
    # relative change is tau times the relative change of sigma(E)
    m = (a["xHI"] < 0.1) & (a["ne"] > 0.0)
    dGo = np.max(np.abs(a["G"][m]/b["G"][m] - 1))
    print(f"    max |Gamma_HI(separate)/Gamma_HI(add_fuv off) - 1| = {dGo:.2e} "
          f"over {m.sum()} ionized leaves")
    # not rounding: the conditional CDF of the band sampler (knots from 4.78
    # eV) and the sampler built on [13.598, 100] eV alone are two
    # piecewise-linear fits of one Planck function, 1e-8 apart
    ok &= dGo < 1e-6
    for q in ("G", "Phi3", "G0", "Hd"):
        mm = np.abs(a[q]) > 0
        dsv = np.max(np.abs(c[q][mm]/a[q][mm] - 1)) if mm.any() else 0.0
        print(f"    served vs static schedule, max rel {q:5s} = {dsv:.2e}")
        ok &= dsv < 1e-9
    print("RESULT:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    if sys.argv[1] == "make":
        make(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        sys.exit(0 if check(sys.argv[2]) else 1)
