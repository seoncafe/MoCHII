#!/usr/bin/env python3
"""PN150 (Lexington/Meudon "Paris" planetary nebula): MoCHII against Cloudy
c25.00 and the published median of Jin, Kewley & Sutherland (2022, ApJ 927,
37) Table 4.

Inputs (paths relative to this directory):
  pn150_converged_lines.txt, pn150_converged_rates.h5, pn150_converged.log
  ../../../tests/cloudy_c25/pn150_c25_{abs,rel}.lin, pn150_c25.{rad,ovr},
  pn150_c25.ele_<el>, pn150_c25.out

Outputs (next to this script):
  compare_pn150.txt   line table, condition table, structure numbers
  pn150_profiles.pdf  T_e and ion-fraction radial profiles, MoCHII spherical
                      averages over Cloudy

    python3 compare_l6.py --base pn150_hbeta_residual/<run>
    (copy of ../compare_pn150.py writing compare_<run>.txt and
    <run>_profiles.pdf in pn150_hbeta_residual/)

MoCHII wavelengths are vacuum values; Cloudy labels optical/IR lines by
their air wavelengths and UV lines (below 2000 A) by vacuum ones.  Every
published blend is the sum of the MoCHII lines of one ion inside the
listed vacuum windows; the number of lines summed is printed beside it.
"""
import argparse
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # m13_restoration
CLOUDY = os.path.join(HERE, "..", "..", "..", "tests", "cloudy_c25")
PC = 3.0856775814913673e18
RIN_CM = 1.0e17

# ------------------------------------------------------------------ lines
# (row label, Jin Table 4 median, MoCHII (element, stage, [vacuum windows A]),
#  Cloudy labels as in pn150_c25_rel.lin)
# A window is (lo, hi); a single line is given as its vacuum wavelength and
# matched within max(0.6 A, 3e-4 lambda).
ROWS = [
    ("He I 5876",          0.099, ("he", 1, [5877.27]),         ["Blnd 5875.66A"]),
    ("He II 4686",         0.324, ("he", 2, [4687.00]),         ["He 2 4685.68A"]),
    ("C II] 2325+",        0.216, ("c", 2, [(2322.0, 2330.0)]), ["blnd 2326.00A"]),
    ("C II 1335",          0.119, ("c", 2, [(1334.0, 1336.5)]), ["blnd 1335.00A"]),
    ("C III] 1907+1909",   1.89,  ("c", 3, [(1906.0, 1910.0)]), ["blnd 1909.00A"]),
    ("C IV 1549+",         2.25,  ("c", 4, [(1547.0, 1552.0)]), ["blnd 1549.00A"]),
    ("[N I] 5200+5198",    0.020, ("n", 1, [(5198.5, 5202.5)]), ["blnd 5199.00A"]),
    ("[N II] 6584+6548",   1.18,  ("n", 2, [6585.27, 6549.86]), ["n  2 6583.45A", "n  2 6548.05A"]),
    ("[N II] 5755",        0.018, ("n", 2, [5756.19]),          ["n  2 5754.59A"]),
    ("N III] 1750 (Jin: N II] 1749+)", 0.106,
                                  ("n", 3, [(1745.0, 1755.0)]), ["blnd 1750.00A"]),
    ("[N III] 57.3um",     0.130, ("n", 3, [573394.5]),         ["n  3 57.3238m"]),
    ("N IV] 1487+",        0.199, ("n", 4, [(1482.0, 1488.0)]), ["blnd 1486.00A"]),
    ("N V 1240+",          0.147, ("n", 5, [(1237.0, 1244.0)]), ["blnd 1240.00A"]),
    ("[O I] 63.1um",       0.028, ("o", 1, [631851.6]),         ["O  1 63.1679m"]),
    ("[O I] 6300+6363",    0.154, ("o", 1, [6302.05, 6365.54]), ["O  1 6300.30A", "O  1 6363.78A"]),
    ("[O II] 3726+3729",   2.28,  ("o", 2, [3727.09, 3729.84]), ["blnd 3727.00A"]),
    ("[O III] 51.8um",     1.39,  ("o", 3, [518134.7]),         ["o  3 51.8004m"]),
    ("[O III] 88.3um",     0.274, ("o", 3, [883392.2]),         ["o  3 88.3323m"]),
    ("[O III] 5007+4959",  22.00, ("o", 3, [5008.24, 4960.29]), ["o  3 5006.84A", "o  3 4958.91A"]),
    ("[O III] 4363",       0.168, ("o", 3, [4364.45]),          ["O  3 4363.21A"]),
    ("[O IV] 25.9um",      4.09,  ("o", 4, [258933.2]),         ["o  4 25.8863m"]),
    ("O IV] 1403+",        0.183, ("o", 4, [(1396.0, 1408.5)]), ["blnd 1402.00A"]),
    ("O V] 1218+",         0.193, ("o", 5, [(1212.5, 1219.5)]), ["blnd 1218.00A"]),
    ("[Ne II] 12.8um",     0.042, ("ne", 2, [128139.4]),        ["ne 2 12.8101m"]),
    ("[Ne III] 15.5um",    2.02,  ("ne", 3, [155545.2]),        ["ne 3 15.5509m"]),
    ("[Ne III] 3869+3968", 2.63,  ("ne", 3, [3869.85, 3968.58]), ["ne 3 3868.76A", "ne 3 3967.47A"]),
    ("[Ne IV] 2423+",      0.723, ("ne", 4, [(2421.5, 2426.5)]), ["blnd 2424.00A"]),
    ("[Ne V] 3426+3346",   0.591, ("ne", 5, [3426.85, 3346.79]), ["ne 5 3425.88A", "ne 5 3345.82A"]),
    ("[Ne V] 24.2um",      0.936, ("ne", 5, [243190.7]),        ["ne 5 24.3108m"]),
    ("Mg II 2798+",        2.240, ("mg", 2, [(2795.0, 2805.0)]), ["blnd 2798.00A"]),
    ("[Mg IV] 4.49um",     0.121, ("mg", 4, [44883.3]),         ["mg 4 4.48711m"]),
    ("Si II] 2335+",       0.160, ("si", 2, [(2328.0, 2352.0)]), ["blnd 2335.00A"]),
    ("Si III] 1892+",      0.325, ("si", 3, [(1882.0, 1893.0)]), ["blnd 1888.00A"]),
    ("Si IV 1397+",        0.172, ("si", 4, [(1393.0, 1403.5)]), ["blnd 1397.00A"]),
    ("[S II] 6716+6731",   0.359, ("s", 2, [6718.29, 6732.67]), ["Blnd 6720.00A"]),
    ("[S II] 4069+4076",   0.064, ("s", 2, [4069.75, 4077.50]), ["S  2 4068.60A", "S  2 4076.35A"]),
    ("[S III] 18.7um",     0.527, ("s", 3, [187055.7]),         ["S  3 18.7078m"]),
    ("[S III] 33.5um",     0.210, ("s", 3, [335008.4]),         ["S  3 33.4704m"]),
    ("[S III] 9532+9069",  2.02,  ("s", 3, [9532.25, 9070.05]), ["S  3 9530.62A", "S  3 9068.62A"]),
    ("[S IV] 10.5um",      2.25,  ("s", 4, [105105.0]),         ["s  4 10.5076m"]),
]
JIN_HBETA = 2.68e35          # Jin Table 4 median L(Hbeta) [erg/s]

# ------------------------------------------------------------ conditions
# Each input item of the two runs: (item, MoCHII, Cloudy c25.00, matched?)
CONDITIONS = [
    ("source", "blackbody 1.5e5 K, R = 1e10 cm", "black body 150000 K, radius 10", "yes"),
    ("bolometric L", "3.607336e37 erg/s (4 pi R^2 sigma T^4)", "same R and T", "yes"),
    ("transported band", "4.78-200 eV (FUV from 4.78 eV; band L 3.464554e37 in 13.598-200)",
     "whole continuum", "no: 1.4e-4 of the H-ionizing L (2.6e-5 of the photons) above 200 eV "
     "is absent; refused above Ne VII 207.27 eV"),
    ("Q(H)", "5.417e47 s^-1 minus the 2.6e-5 above 200 eV", "5.417e47 s^-1 (same blackbody)", "yes, to 2.6e-5"),
    ("n_H", "3000 cm^-3, uniform", "hden 3.4771213, constant density", "yes"),
    ("inner radius", "1e17 cm (0.0324078 pc)", "radius 17", "yes"),
    ("geometry", "3D Cartesian octree, level 7, dx = 2.5e-3 pc", "1D sphere, 187 zones", "no: front sampled by ~2 cells (10-90% of x(H0))"),
    ("outer boundary", "gas to 0.16 pc, te_min 3000 K, all gas counted",
     "stops at 4000 K (default), 0.13516 pc, x(H0) = 0.967", "no: MoCHII counts the neutral gas beyond"),
    ("He, C, N, O, Ne, Mg, Si, S", "0.1, 2.999e-4, 1e-4, 5.998e-4, 1.500e-4, 2.999e-5, 2.999e-5, 1.500e-5",
     "pn_paris abund line (same numbers)", "yes"),
    ("other elements", "none (Cl, Ar, Fe, Ca at 0)", "ism.ini + element off: all but H He C N O Ne Mg Si S off", "yes"),
    ("metal stages", "I-VII for every element", "all stages", "no for stages above VII (Ne VIII .. reached only above 207 eV, outside the band)"),
    ("dust", "none", "no grains command", "yes"),
    ("CMB", "none", "not included (Cloudy warns)", "yes"),
    ("diffuse field", "transported (H, He I, He II cascades, metal free-bound, free-free)", "outward-only, Cloudy's own", "method differs"),
    ("H I lines", "SH95 case B (on the spot Lyman lines)", "H I 4861 total from the model atom", "method differs"),
    ("line intensities", "sum over the leaves, no line transfer", "intrinsic (save line list default), 1 iteration", "not checked for the resonance lines"),
    ("age", "steady state", "age 1e4 yr (a check on equilibrium time scales)", "yes"),
    ("convergence", "drift test, conv_crit vol, gas_tol 5e-3", "single iteration, no iterate", "differs"),
]


def read_mochii_lines(fn):
    lhb, rows = None, []
    for ln in open(fn):
        if ln.startswith("#"):
            if "L(Hbeta)" in ln and "=" in ln:
                lhb = float(ln.split("=")[1].split()[0])
            continue
        t = ln.split()
        if len(t) >= 5:
            rows.append((t[0].lower(), int(t[1]), float(t[2]), float(t[3])))
    return lhb, rows


def mochii_blend(rows, el, stage, wins):
    tot, n = 0.0, 0
    for w in wins:
        if isinstance(w, tuple):
            lo, hi = w
            hit = [r for r in rows if r[0] == el and r[1] == stage
                   and lo <= r[2] <= hi]
        else:
            tol = max(0.6, 3.0e-4 * w)
            hit = [r for r in rows if r[0] == el and r[1] == stage
                   and abs(r[2] - w) < tol]
            if len(hit) > 1:     # keep the nearest one for a single line
                hit = [min(hit, key=lambda r: abs(r[2] - w))]
        tot += sum(r[3] for r in hit)
        n += len(hit)
    return tot, n


def read_cloudy_lin(fn):
    out = {}
    for ln in open(fn):
        p = ln.rstrip("\n").split("\t")
        if len(p) == 2 and not ln.startswith("#"):
            try:
                out[p[0].strip()] = float(p[1])
            except ValueError:
                pass
    return out


def read_cloudy_ele(fn):
    with open(fn) as f:
        head = f.readline().lstrip("#").split()
    d = np.loadtxt(fn, comments="#")
    return head, d


def cloudy_averages(fn):
    """Volume Te(NeNp) and Te from the Averaged Quantities block."""
    lines = open(fn).read().splitlines()
    for i, ln in enumerate(lines):
        if "Averaged Quantities" in ln:
            hdr = lines[i + 1]
            vol = [l for l in lines[i + 2:i + 6] if l.strip().startswith("Volume:")][0]
            vals = [float(x) for x in vol.split(":")[1].split()]
            return vals[0], vals[2]      # Te, Te(NeNp)
    return None, None


def mochii_state(fn):
    import h5py
    st = {}
    with h5py.File(fn, "r") as f:
        dist = float(f["Gamma_HI"].attrs["DIST_CM"][0])
        xyz = np.asarray(f["LeafXYZ"]["data"][:])
        if xyz.shape[0] != 3:
            xyz = xyz.T
        st["r"] = np.sqrt((xyz ** 2).sum(axis=0)) * dist
        st["size"] = f["LeafSize"]["data"][:].ravel() * dist
        for k in ("T_e", "n_e", "x_HI", "x_HII", "x_HeI", "x_HeII", "x_HeIII"):
            st[k] = f[k]["data"][:].ravel()
        st["metal"] = {}
        for k in f.keys():
            if k.startswith("x_") and k.count("_") == 2 and k.split("_")[1] not in ("H", "He"):
                st["metal"][k] = np.asarray(f[k]["data"][:])
        st["keys"] = list(f.keys())
    return st


def shell_profile(r, y, w, edges):
    idx = np.digitize(r, edges)
    prof = np.full(len(edges) - 1, np.nan)
    for i in range(1, len(edges)):
        m = idx == i
        if m.any() and w[m].sum() > 0:
            prof[i - 1] = (y[m] * w[m]).sum() / w[m].sum()
    return prof


def front_radius(rc, prof, level, rising):
    """First radius where a shell profile crosses `level`."""
    ok = np.isfinite(prof)
    rc, prof = rc[ok], prof[ok]
    for i in range(1, len(prof)):
        a, b = prof[i - 1], prof[i]
        if (rising and a < level <= b) or (not rising and a > level >= b):
            return rc[i - 1] + (level - a) * (rc[i] - rc[i - 1]) / (b - a)
    return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="pn150_converged")
    args = ap.parse_args()
    base = os.path.join(HERE, args.base)
    out = []
    pr = out.append

    lhb, mrows = read_mochii_lines(base + "_lines.txt")
    crel = read_cloudy_lin(os.path.join(CLOUDY, "pn150_c25_rel.lin"))
    cabs = read_cloudy_lin(os.path.join(CLOUDY, "pn150_c25_abs.lin"))
    c_hb = cabs["H  1 4861.32A"]
    c_cab = cabs["Ca B 4861.32A"]

    # ------------------------------------------------------- conditions
    pr("Input conditions, MoCHII against Cloudy c25.00")
    for it, a, b, m in CONDITIONS:
        pr("  %-26s MoCHII: %s" % (it, a))
        pr("  %-26s Cloudy: %s" % ("", b))
        pr("  %-26s matched: %s" % ("", m))
    pr("")

    # ------------------------------------------------------- line table
    pr("PN150 line table: L/L(Hbeta), Hbeta = H I 4861 total in all three")
    pr("(Cloudy's own deck normalizes to Ca B 4861; Ca B / H I = %.5f, so a"
       % (c_cab / c_hb))
    pr(" Ca B-normalized Cloudy value is %.4f x the one below)" % (c_hb / c_cab))
    pr("")
    fmt = "%-32s %9s %9s %9s  %7s %7s %7s  %s"
    pr(fmt % ("line", "MoCHII", "C25", "Jin med", "Mo/C25", "Mo/Jin", "C25/Jin", "n(Mo)"))
    pr(fmt % ("L(Hbeta) [1e35 erg/s]", "%.4f" % (lhb / 1e35), "%.4f" % (c_hb / 1e35),
              "%.2f" % (JIN_HBETA / 1e35), "%.3f" % (lhb / c_hb),
              "%.3f" % (lhb / JIN_HBETA), "%.3f" % (c_hb / JIN_HBETA), ""))
    table = []
    for label, med, (el, st, wins), clabs in ROWS:
        mo, n = mochii_blend(mrows, el, st, wins)
        mo /= lhb
        cl = sum(crel[k] for k in clabs)
        table.append((label, mo, cl, med))
        pr(fmt % (label, "%.4g" % mo, "%.4g" % cl, "%.4g" % med,
                  "%.3f" % (mo / cl) if cl > 0 else "-",
                  "%.3f" % (mo / med), "%.3f" % (cl / med), n))
    r = np.array([t[1] / t[2] for t in table])
    pr("")
    pr("MoCHII/C25 over the %d rows: median %.3f, within 10%% %d, within 25%% %d"
       % (len(r), np.median(r), (abs(r - 1) <= 0.10).sum(), (abs(r - 1) <= 0.25).sum()))
    pr("")
    tab = {t[0]: t for t in table}

    def c_over_jin(lab):
        return tab[lab][2] / tab[lab][3]
    pr("Known c25.00 offsets from the Jin median, measured on this Cloudy run")
    pr("(H I 4861 normalization):")
    for lab in ("[S II] 6716+6731", "[S II] 4069+4076", "[S III] 9532+9069",
                "He II 4686", "[Ne III] 3869+3968", "[O III] 5007+4959"):
        pr("  %-22s C25/Jin %.3f   MoCHII/Jin %.3f" % (lab, c_over_jin(lab),
                                                     tab[lab][1] / tab[lab][3]))
    pr("  %-22s C25/Jin %.3f   MoCHII/Jin %.3f" % ("L(Hbeta)", c_hb / JIN_HBETA,
                                                 lhb / JIN_HBETA))

    # ------------------------------------------------------- structure
    st = mochii_state(base + "_rates.h5")
    gas = st["n_e"] > 0
    vol = st["size"] ** 3
    r = st["r"][gas] / PC
    v = vol[gas]
    edges = np.linspace(RIN_CM / PC, r.max(), 121)
    rc = 0.5 * (edges[1:] + edges[:-1])
    xh0 = shell_profile(r, st["x_HI"][gas], v, edges)
    xhe3 = shell_profile(r, st["x_HeIII"][gas], v, edges)
    te = shell_profile(r, st["T_e"][gas], v, edges)
    ne, xp = st["n_e"][gas], st["x_HII"][gas]
    w = ne * xp * v
    tnpne = (st["T_e"][gas] * w).sum() / w.sum()
    c_te, c_tnpne = cloudy_averages(os.path.join(CLOUDY, "pn150_c25.out"))

    ch_h = read_cloudy_ele(os.path.join(CLOUDY, "pn150_c25.ele_h"))[1]
    ch_he = read_cloudy_ele(os.path.join(CLOUDY, "pn150_c25.ele_he"))[1]
    crh = (RIN_CM + ch_h[:, 0]) / PC
    crhe = (RIN_CM + ch_he[:, 0]) / PC

    def cfront(rr, y, level, rising):
        for i in range(1, len(y)):
            a, b = y[i - 1], y[i]
            if (rising and a < level <= b) or (not rising and a > level >= b):
                return rr[i - 1] + (level - a) * (rr[i] - rr[i - 1]) / (b - a)
        return np.nan

    pr("")
    pr("Structure (MoCHII: volume-weighted shell averages, 120 shells)")
    pr("%-34s %10s %10s %8s" % ("quantity", "MoCHII", "C25", "Mo/C25"))
    for lab, lev, rising, my, cy, cr in [
            ("H front x(H0)=0.1 [pc]", 0.1, True, xh0, ch_h[:, 1], crh),
            ("H front x(H0)=0.5 [pc]", 0.5, True, xh0, ch_h[:, 1], crh),
            ("H front x(H0)=0.9 [pc]", 0.9, True, xh0, ch_h[:, 1], crh),
            ("He III front x(He2+)=0.9 [pc]", 0.9, False, xhe3, ch_he[:, 3], crhe),
            ("He III front x(He2+)=0.5 [pc]", 0.5, False, xhe3, ch_he[:, 3], crhe),
            ("He III front x(He2+)=0.1 [pc]", 0.1, False, xhe3, ch_he[:, 3], crhe)]:
        a = front_radius(rc, my, lev, rising)
        b = cfront(cr, cy, lev, rising)
        pr("%-34s %10.5f %10.5f %8.4f" % (lab, a, b, a / b))
    pr("%-34s %10.1f %10.1f %8.4f" % ("<T[Np Ne]> [K]", tnpne, c_tnpne, tnpne / c_tnpne))
    inner = np.argsort(r)[:200]
    ovr = np.loadtxt(os.path.join(CLOUDY, "pn150_c25.ovr"), comments="#", usecols=(0, 1))
    pr("%-34s %10.1f %10.1f %8.4f" % ("T_e at the inner edge [K] (200 leaves)",
                                       st["T_e"][gas][inner].mean(), ovr[0, 1],
                                       st["T_e"][gas][inner].mean() / ovr[0, 1]))
    wl = ne * v
    he_ratio = (ne * st["x_HeII"][gas] * v).sum() / (ne * xp * v).sum()
    he3_ratio = (ne * st["x_HeIII"][gas] * v).sum() / (ne * xp * v).sum()
    pr("%-34s %10.4f" % ("<He+>/<H+> (n_e weighted, x n_H/n_He)", he_ratio))
    pr("%-34s %10.4f" % ("<He2+>/<H+>", he3_ratio))

    # ------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["text.usetex"] = True
    roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"]
    elements = [("h", "H", 2), ("he", "He", 3), ("c", "C", 7), ("n", "N", 7),
                ("o", "O", 7), ("ne", "Ne", 7), ("mg", "Mg", 7), ("si", "Si", 7),
                ("s", "S", 7)]
    fig, axes = plt.subplots(5, 2, figsize=(10, 15), sharex=True)
    axes = axes.ravel()
    ax = axes[0]
    ax.plot(rc, te, "k-", label=r"\textsc{MoCHII}")
    ax.plot((RIN_CM + ovr[:, 0]) / PC, ovr[:, 1], "r--", label=r"Cloudy c25.00")
    ax.set_ylabel(r"$T_e$ [K]")
    ax.legend(fontsize=8)
    colors = plt.cm.viridis(np.linspace(0, 0.95, 8))
    for k, (el, name, nst) in enumerate(elements):
        ax = axes[k + 1]
        hd, cd = read_cloudy_ele(os.path.join(CLOUDY, "pn150_c25.ele_" + el))
        crr = (RIN_CM + cd[:, 0]) / PC
        if el == "h":
            mys = [st["x_HI"], st["x_HII"]]
        elif el == "he":
            mys = [st["x_HeI"], st["x_HeII"], st["x_HeIII"]]
        else:
            key = [kk for kk in st["metal"] if kk.lower() == "x_%s_stages" % el]
            mys = []
            if key:
                arr = st["metal"][key[0]]
                if arr.shape[0] != st["r"].size:
                    arr = arr.T
                mys = [arr[:, i] for i in range(arr.shape[1])]
        for i in range(min(nst, len(mys), cd.shape[1] - 1)):
            p = shell_profile(r, mys[i][gas], v, edges)
            lab = name + r"\,\textsc{" + roman[i].lower() + r"}"
            ax.plot(rc, p, "-", color=colors[i], label=lab)
            ax.plot(crr, cd[:, i + 1], "--", color=colors[i])
        ax.set_ylim(-0.02, 1.02)
        ax.set_ylabel(name + r" stage fraction")
        ax.legend(fontsize=7, ncol=2)
    for ax in axes[-2:]:
        ax.set_xlabel(r"$r$ [pc]")
    axes[0].set_title(r"PN150: \textsc{MoCHII} (solid, shell average) and Cloudy c25.00 (dashed)",
                      fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "pn150_hbeta_residual", args.base.split("/")[-1] + "_profiles.pdf"))

    txt = "\n".join(out) + "\n"
    open(os.path.join(HERE, "pn150_hbeta_residual", "compare_" + args.base.split("/")[-1] + ".txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
