"""Generate the two Lexington benchmark tables of Wood, Mathis & Ercolano
(2004, MNRAS 348, 1337) with a MoCHII column appended.

  wood_hii40.tex   -- their Table 1 (HII40, 40 kK blackbody)
  wood_hii20.tex   -- their Table 2 (HII20, 20 kK blackbody)

Both are \\input by mochii.tex.

Sources, and nothing else:

* Columns GF ... WME are transcribed from the LaTeX source of the paper
  (arXiv:astro-ph/0311584, file kwood.tex, the two `table*` environments).
  A dash in the published table becomes \\nodata.  The GF ... BE columns
  agree entry by entry with Tables 4 and 5 of Ercolano et al. (2003), which
  the paper already uses for the published benchmark columns.
  One correction is applied and flagged in the table note: the HII40 H-beta
  row is headed 10^36 erg/s in the published table, while the tabulated
  values (2.01-2.10) are the 10^37 erg/s values of Ercolano et al. (2003).

* The MoCHII line column is computed here from the converged line tables of
  the continuous-energy production runs,

      ../../results/continuous_energy/production/hii40_continuous_lines.txt
      ../../results/continuous_energy/production/hii20_continuous_lines.txt

  (2^25 packets, level-6 octree, conv_crit='vol'), by summing exactly the
  components that make up each published blend.  MoCHII lists vacuum
  wavelengths.

* The two MoCHII temperatures are computed here from the rates files of the
  same runs,

      ../../results/continuous_energy/production/hii40_continuous_rates.h5
      ../../results/continuous_energy/production/hii20_continuous_rates.h5

  as sum(T_e n_p n_e) / sum(n_p n_e) over the gas cells, which is the
  published <T[Np Ne]> definition.  The level-6 grid is uniform, so every
  cell carries the same volume and no volume weight is needed; n_H is 100 by
  construction wherever there is gas, the same closure check_g2_hii.py uses.
  Quantities MoCHII does not produce for these runs (the Balmer jump, the
  inner-edge temperature, the outer radius, and the published helium ratio)
  are left as \\nodata and explained in the table note.

* The C25 column is Cloudy c25.00, run for this paper from

      ../../tests/cloudy_c25/hii40_c25.in
      ../../tests/cloudy_c25/hii20_c25.in

  whose physics blocks are the unmodified c25.00 test-suite inputs
  tsuite/auto/hii_paris.in (HII40) and tsuite/auto/hii_coolstar.in (HII20);
  only output commands were appended.  Line luminosities are read from the
  `save line list ... absolute` files hii{40,20}_c25_abs.lin, driven by the
  list ../../tests/cloudy_c25/wood_lines.dat, and are divided here by the
  Cloudy "H  1 4861.32A" luminosity of the same file.  The averaged
  quantities come from hii{40,20}_c25.out (the "Volume:" row of the last
  "Averaged Quantities" block for <T[Np Ne]>, and the "Log10 Mean Ionisation
  (over volume)" block for He+/H+) and from hii{40,20}_c25.ovr (first zone
  for T_inner, last zone depth + R_in for R_out).  The Cloudy blends are
  spelled out in CLOUDY_BLEND below.

Run:  python3 make_wood_tables.py
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
CLOUDY = os.path.join(HERE, "..", "..", "tests", "cloudy_c25")

CODES = ["GF", "C25", "HN", "DP", "TK", "PH", "RS", "RR", "BE", "WME"]

# inner radius of both benchmarks [cm]; Cloudy reports depth, not radius
RIN_CM = 3.0e18

# ---------------------------------------------------------------- MoCHII
# element, ionization stage, and the MoCHII vacuum wavelengths [A] that make
# up each published blend.
BLEND = {
    # the published label of the He I line is its air wavelength, 5876; the
    # MoCHII table lists the vacuum value of the same line
    "He I 5876":        ("he", 1, [5877.27]),
    "C II 2325+":       ("c",  2, [2322.70, 2324.27, 2325.41, 2326.12,
                                   2327.70, 2328.84]),
    "C III 1907+1909":  ("c",  3, [1906.68, 1908.73]),
    "N II 122":         ("n",  2, [1218026.80]),
    "N II 6584+6548":   ("n",  2, [6585.27, 6549.86]),
    "N II 5755":        ("n",  2, [5756.19]),
    "N III 57.3":       ("n",  3, [573394.50]),
    "O I 6300+6363":    ("o",  1, [6302.05, 6365.54]),
    "O II 7320+7330":   ("o",  2, [7321.13, 7322.20, 7331.75, 7332.83]),
    "O II 3726+3729":   ("o",  2, [3727.09, 3729.84]),
    "O III 52+88":      ("o",  3, [518134.72, 883392.23]),
    "O III 5007+4959":  ("o",  3, [5008.24, 4960.29]),
    "O III 4363":       ("o",  3, [4364.45]),
    "Ne II 12.8":       ("ne", 2, [128139.42]),
    "Ne III 15.5":      ("ne", 3, [155545.19]),
    "Ne III 3869+3968": ("ne", 3, [3869.85, 3968.58]),
    "S II 6716+6731":   ("s",  2, [6718.29, 6732.67]),
    "S II 4068+4076":   ("s",  2, [4069.75, 4077.50]),
    "S III 18.7":       ("s",  3, [187055.74]),
    "S III 9532+9069":  ("s",  3, [9532.25, 9070.05]),
}

# ---------------------------------------------------------------- Cloudy
# The Cloudy c25.00 components summed for each published blend.  Labels are
# exactly as they appear in the *_c25_abs.lin files (label field, tab, value).
# "blnd"/"Blnd" entries are Cloudy multiplet blends; their own components are
# fixed in c25.00 data/blends.ini and are listed here so the two codes can be
# checked against each other component by component:
#   blnd 2326 = C II]  2323.50 2324.69 2325.40 2326.93 2328.12 (air)
#   blnd 1909 = C III] 1906.68 1908.73 (vacuum)
#   blnd 3727 = [O II] 3726.03 3728.81 (air)
#   blnd 7325 = [O II] 7318.92 7319.99 7329.67 7330.73 (air)
#   Blnd 6720 = [S II] 6716.44 6730.82 (air)
#   Blnd 5875.66A is the He I 5876 (3 3D - 2 3P) multiplet total.
# [O I] is summed from the two individual lines rather than Cloudy's
# blnd 6300, because that blend also carries 6391.73, which the published
# 6300+6363 row does not include.
CLOUDY_HBETA = "H  1 4861.32A"
CLOUDY_BLEND = {
    "He I 5876":        ["Blnd 5875.66A"],
    "C II 2325+":       ["blnd 2326.00A"],
    "C III 1907+1909":  ["blnd 1909.00A"],
    "N II 122":         ["n  2 121.769m"],
    "N II 6584+6548":   ["n  2 6583.45A", "n  2 6548.05A"],
    "N II 5755":        ["n  2 5754.59A"],
    "N III 57.3":       ["n  3 57.3238m"],
    "O I 6300+6363":    ["O  1 6300.30A", "O  1 6363.78A"],
    "O II 7320+7330":   ["blnd 7325.00A"],
    "O II 3726+3729":   ["blnd 3727.00A"],
    "O III 52+88":      ["o  3 51.8004m", "o  3 88.3323m"],
    "O III 5007+4959":  ["o  3 5006.84A", "o  3 4958.91A"],
    "O III 4363":       ["O  3 4363.21A"],
    "Ne II 12.8":       ["ne 2 12.8101m"],
    "Ne III 15.5":      ["ne 3 15.5509m"],
    "Ne III 3869+3968": ["ne 3 3868.76A", "ne 3 3967.47A"],
    "S II 6716+6731":   ["Blnd 6720.00A"],
    "S II 4068+4076":   ["S  2 4068.60A", "S  2 4076.35A"],
    "S III 18.7":       ["S  3 18.7078m"],
    "S III 9532+9069":  ["S  3 9530.62A", "S  3 9068.62A"],
}


def read_cloudy_lin(fn):
    """Absolute line luminosities [erg/s] of a `save line list` file."""
    out = {}
    for ln in open(fn):
        if ln.startswith("#") or "\t" not in ln:
            continue
        lab, val = ln.split("\t")
        out[lab.rstrip()] = float(val)
    return out


def read_cloudy_out(fn):
    """<T[Np Ne]> and the volume-averaged He+/H+ of the last iteration."""
    txt = open(fn, errors="replace").read().splitlines()
    tnpne = he_over_h = None
    for i, ln in enumerate(txt):
        if "Averaged Quantities" in ln:
            for sub in txt[i:i + 6]:
                # columns: Te  Te(Ne)  Te(NeNp)  ...
                if sub.startswith(" Volume:"):
                    tnpne = float(sub.split()[3])
        if "Log10 Mean Ionisation (over volume)" in ln:
            # fixed width: 11-char element label, then 7-char log10 fields
            def logs(row):
                return [float(row[11 + 7 * j: 18 + 7 * j])
                        for j in range(3) if row[11 + 7 * j: 18 + 7 * j].strip()
                        and "(" not in row[11 + 7 * j: 18 + 7 * j]]
            hyd = logs(txt[i])          # H0, H+
            hel = logs(txt[i + 1])      # He0, He+, He++
            he_over_h = 10.0 ** (hel[1] - hyd[1])
    if tnpne is None or he_over_h is None:
        raise SystemExit("could not read averaged quantities from " + fn)
    return {"tnpne": tnpne, "he_over_h": he_over_h}


def read_cloudy_ovr(fn):
    """Inner-edge temperature and outer radius from the structure file."""
    rows = [ln.split("\t") for ln in open(fn) if not ln.startswith("#")]
    return {"tinner": float(rows[0][1]),
            "rout": RIN_CM + float(rows[-1][0])}


def read_cloudy(tag):
    cl = read_cloudy_out(os.path.join(CLOUDY, "%s_c25.out" % tag))
    cl.update(read_cloudy_ovr(os.path.join(CLOUDY, "%s_c25.ovr" % tag)))
    cl["lines"] = read_cloudy_lin(
        os.path.join(CLOUDY, "%s_c25_abs.lin" % tag))
    return cl


def cloudy_cell(label, fmt, key, tag, cl):
    """The Cloudy c25.00 entry of one printed row."""
    hb = cl["lines"][CLOUDY_HBETA]
    if r"T_{\mathrm{inner}}" in label:
        return "%.0f" % cl["tinner"]
    if r"R_{\mathrm{out}}" in label:
        unit = 1.0e19 if tag == "hii40" else 1.0e18
        return "%.2f" % (cl["rout"] / unit)
    if r"He^{+}" in label:
        return "%.3f" % cl["he_over_h"]
    if key is None:                       # Balmer jump: not extracted here
        return r"\nodata"
    if key == "TNPNE":
        return "%.0f" % cl["tnpne"]
    if key == "HBETA":
        return fmt % (hb / (1.0e37 if tag == "hii40" else 1.0e36))
    return fmt % (sum(cl["lines"][c] for c in CLOUDY_BLEND[key]) / hb)


# ---------------------------------------------------------------- MoCHII
def mochii_tnpne(fn):
    """<T[Np Ne]> = sum(T_e n_p n_e) / sum(n_p n_e) over the gas cells.

    The benchmark grid is uniform at level 6, so every cell has the same
    volume and the volume weight cancels; n_H = 100 wherever there is gas,
    the closure also used by the HII benchmark checks.
    """
    import h5py
    import numpy as np
    with h5py.File(fn, "r") as f:
        te = f["T_e"]["data"][:]
        ne = f["n_e"]["data"][:]
        xhi = f["x_HI"]["data"][:]
    nh = np.where(ne > 0, 100.0, 0.0)
    w = nh * (1.0 - xhi) * ne
    return float((te * w).sum() / w.sum())



def mochii_structure(fn, tag):
    """T_inner, the ionization-front radius, and <He+>/<H+> from a rates file.

    T_inner  : mean T_e of the cells at the illuminated inner edge.
    R_front  : radius where the shell-averaged x_HII falls through 0.5 --- the
               benchmark's R_out for a radiation-bounded model.  Returned as
               None when the H II region still fills the outermost gas shell
               (the front is then set by the grid, not predicted).
    he_ratio : Ercolano et al. (2003) Eq. 17,
               <He+>/<H+> = int n_e n(He+) dV / int n_e n_p dV x n_H/n_He,
               which for a uniform medium is the n_e-weighted x_HeII over the
               n_e-weighted x_HII.
    """
    import h5py
    import numpy as np
    with h5py.File(fn, "r") as f:
        dist = float(f["Gamma_HI"].attrs["DIST_CM"][0])
        xyz = np.asarray(f["LeafXYZ"]["data"][:])
        if xyz.shape[0] != 3:
            xyz = xyz.T
        size = f["LeafSize"]["data"][:].ravel()
        te = f["T_e"]["data"][:].ravel()
        ne = f["n_e"]["data"][:].ravel()
        xhi = f["x_HI"]["data"][:].ravel()
        xheii = f["x_HeII"]["data"][:].ravel()
    r = np.sqrt((xyz ** 2).sum(axis=0)) * dist
    vol = (size * dist) ** 3
    xhii = 1.0 - xhi
    gas = ne > 0
    # inner edge: the 50 cells closest to the source
    inner = np.argsort(r[gas])[:50]
    tinner = float(te[gas][inner].mean())
    # ionization front from the shell-averaged profile
    rg, xg = r[gas], xhii[gas]
    edges = np.linspace(0.0, rg.max(), 80)
    idx = np.digitize(rg, edges)
    prof, rc = [], []
    for i in range(1, len(edges)):
        m = idx == i
        if m.sum():
            prof.append(xg[m].mean())
            rc.append(0.5 * (edges[i] + edges[i - 1]))
    prof, rc = np.asarray(prof), np.asarray(rc)
    rfront = None
    if prof[-1] < 0.5:                     # front captured inside the medium
        rfront = float(np.interp(0.5, prof[::-1], rc[::-1]))
    # Ercolano Eq. 17
    he_ratio = float((ne * xheii * vol).sum() / (ne * xhii * vol).sum())
    return tinner, rfront, he_ratio


def read_lines(fn):
    lhb, rows = None, []
    for ln in open(fn):
        if ln.startswith("#"):
            if "L(Hbeta)" in ln and "=" in ln:
                lhb = float(ln.split("=")[1].split()[0])
            continue
        t = ln.split()
        if len(t) >= 5:
            rows.append((t[0], int(t[1]), float(t[2]), float(t[4])))
    return lhb, rows


def blend_value(rows, key):
    el, stage, wls = BLEND[key]
    tot = 0.0
    for w in wls:
        hit = [r for r in rows if r[0] == el and r[1] == stage
               and abs(r[2] - w) < max(0.5, 1e-5 * w)]
        if len(hit) != 1:
            raise SystemExit("%s: %d matches for %.2f A" % (key, len(hit), w))
        tot += hit[0][3]
    return tot


# ------------------------------------------------------- published columns
# Wood et al. (2004), Table 1: HII40.  Order as printed.
HII40 = [
    # (row label, format, [GF, HN, DP, TK, PH, RS, RR, BE, WME])
    (r"$L(\mathrm{H}\beta)/10^{37}\,\mathrm{erg\,s^{-1}}$", "%.2f",
     ["2.06", "2.02", "2.02", "2.10", "2.05", "2.07", "2.05", "2.02", "2.01"],
     "HBETA"),
    (r"\HeI{} 5876", "%.3f",
     ["0.119", "0.112", "0.113", "0.116", "0.118", "0.116", "-", "0.114",
      "0.114"], "He I 5876"),
    (r"C\,\textsc{ii}] 2325+", "%.3f",
     ["0.157", "0.141", "0.139", "0.110", "0.166", "0.096", "0.178", "0.148",
      "0.172"], "C II 2325+"),
    (r"C\,\textsc{iii}] 1907+1909", "%.3f",
     ["0.071", "0.076", "0.069", "0.091", "0.060", "0.066", "0.074", "0.041",
      "0.078"], "C III 1907+1909"),
    (r"$[$N\,\textsc{ii}$]$ 122\,$\mu$m", "%.3f",
     ["0.027", "0.037", "0.034", "-", "0.032", "0.035", "0.030", "0.036",
      "0.031"], "N II 122"),
    (r"$[$N\,\textsc{ii}$]$ 6584+6548", "%.3f",
     ["0.669", "0.817", "0.725", "0.69", "0.736", "0.723", "0.807", "0.852",
      "0.742"], "N II 6584+6548"),
    (r"$[$N\,\textsc{ii}$]$ 5755", "%.4f",
     ["0.0050", "0.0054", "0.0050", "-", "0.0064", "0.0050", "0.0068",
      "0.0061", "0.0057"], "N II 5755"),
    (r"$[$N\,\textsc{iii}$]$ 57.3\,$\mu$m", "%.3f",
     ["0.306", "0.261", "0.311", "-", "0.292", "0.273", "0.301", "0.223",
      "0.308"], "N III 57.3"),
    (r"$[$O\,\textsc{i}$]$ 6300+6363", "%.4f",
     ["0.0094", "0.0086", "0.0088", "0.012", "0.0059", "0.0070", "-",
      "0.0065", "0.011"], "O I 6300+6363"),
    (r"$[$O\,\textsc{ii}$]$ 7320+7330", "%.3f",
     ["0.029", "0.030", "0.031", "0.023", "0.032", "0.024", "0.036", "0.025",
      "0.033"], "O II 7320+7330"),
    (r"$[$O\,\textsc{ii}$]$ 3726+3729", "%.2f",
     ["1.94", "2.17", "2.12", "1.6", "2.19", "1.88", "2.26", "1.92", "2.23"],
     "O II 3726+3729"),
    (r"$[$O\,\textsc{iii}$]$ 52+88\,$\mu$m", "%.2f",
     ["2.35", "2.10", "2.26", "2.17", "2.34", "2.29", "2.34", "2.28", "2.42"],
     "O III 52+88"),
    (r"$[$O\,\textsc{iii}$]$ 5007+4959", "%.2f",
     ["2.21", "2.38", "2.20", "3.27", "1.93", "2.17", "2.08", "1.64", "2.34"],
     "O III 5007+4959"),
    (r"$[$O\,\textsc{iii}$]$ 4363", "%.4f",
     ["0.00235", "0.0046", "0.0041", "0.0070", "0.0032", "0.0040", "0.0035",
      "0.0022", "0.0044"], "O III 4363"),
    (r"$[$Ne\,\textsc{ii}$]$ 12.8\,$\mu$m", "%.3f",
     ["0.177", "0.195", "0.192", "-", "0.181", "0.217", "0.196", "0.212",
      "0.192"], "Ne II 12.8"),
    (r"$[$Ne\,\textsc{iii}$]$ 15.5\,$\mu$m", "%.3f",
     ["0.294", "0.264", "0.270", "0.35", "0.429", "0.350", "0.417", "0.267",
      "0.263"], "Ne III 15.5"),
    (r"$[$Ne\,\textsc{iii}$]$ 3869+3968", "%.3f",
     ["0.084", "0.087", "0.071", "0.092", "0.087", "0.083", "0.086", "0.053",
      "0.063"], "Ne III 3869+3968"),
    (r"$[$S\,\textsc{ii}$]$ 6716+6731", "%.3f",
     ["0.137", "0.166", "0.153", "0.315", "0.155", "0.133", "0.130", "0.141",
      "0.18"], "S II 6716+6731"),
    (r"$[$S\,\textsc{ii}$]$ 4068+4076", "%.4f",
     ["0.0093", "0.0090", "0.0100", "0.026", "0.0070", "0.005", "0.0060",
      "0.0060", "0.0094"], "S II 4068+4076"),
    (r"$[$S\,\textsc{iii}$]$ 18.7\,$\mu$m", "%.3f",
     ["0.627", "0.750", "0.726", "0.535", "0.556", "0.567", "0.580", "0.574",
      "0.623"], "S III 18.7"),
    (r"$[$S\,\textsc{iii}$]$ 9532+9069", "%.2f",
     ["1.13", "1.19", "1.16", "1.25", "1.23", "1.25", "1.28", "1.21", "1.44"],
     "S III 9532+9069"),
    (r"$10^{3}\times\Delta(\mathrm{BC}\,3645)$/\AA", None,
     ["4.88", "-", "4.95", "-", "5.00", "5.70", "-", "5.47", "4.96"], None),
    (r"$T_{\mathrm{inner}}$ / K", None,
     ["7719", "7668", "7663", "8318", "7440", "7644", "7399", "7370", "7743"],
     "TINNER"),
    (r"$\langle T[N_{\mathrm{p}}N_{\mathrm{e}}]\rangle$ / K", None,
     ["7940", "7936", "8082", "8199", "8030", "8022", "8060", "7720", "8183"],
     "TNPNE"),
    (r"$R_{\mathrm{out}}/10^{19}$\,cm", "%.2f",
     ["1.46", "1.46", "1.46", "1.45", "1.46", "1.47", "1.46", "1.46", "1.45"],
     "RFRONT"),
    (r"$\langle\mathrm{He^{+}}\rangle/\langle\mathrm{H^{+}}\rangle$", "%.3f",
     ["0.787", "0.727", "0.754", "0.77", "0.764", "0.804", "0.829", "0.715",
      "0.771"], "HEHRATIO"),
]

# Wood et al. (2004), Table 2: HII20.
HII20 = [
    (r"$L(\mathrm{H}\beta)/10^{36}\,\mathrm{erg\,s^{-1}}$", "%.2f",
     ["4.85", "4.85", "4.83", "4.90", "4.93", "5.04", "4.89", "4.97", "4.87"],
     "HBETA"),
    (r"\HeI{} 5876", "%.4f",
     ["0.0072", "0.008", "0.0073", "0.008", "0.0074", "0.0110", "-", "0.0069",
      "0.00684"], "He I 5876"),
    (r"C\,\textsc{ii}] 2325+", "%.3f",
     ["0.054", "0.047", "0.046", "0.040", "0.060", "0.038", "0.063", "0.042",
      "0.053"], "C II 2325+"),
    (r"$[$N\,\textsc{ii}$]$ 122\,$\mu$m", "%.3f",
     ["0.068", "-", "0.072", "0.007", "0.072", "0.071", "0.071", "0.071",
      "0.067"], "N II 122"),
    (r"$[$N\,\textsc{ii}$]$ 6584+6548", "%.3f",
     ["0.745", "0.786", "0.785", "0.925", "0.843", "0.803", "0.915", "0.846",
      "0.778"], "N II 6584+6548"),
    (r"$[$N\,\textsc{ii}$]$ 5755", "%.4f",
     ["0.0028", "0.0024", "0.0023", "0.0029", "0.0033", "0.0030", "0.0033",
      "0.0025", "0.0025"], "N II 5755"),
    (r"$[$N\,\textsc{iii}$]$ 57.3\,$\mu$m", "%.4f",
     ["0.0040", "0.0030", "0.0032", "0.0047", "0.0031", "0.0020", "0.0022",
      "0.0019", "0.0049"], "N III 57.3"),
    (r"$[$O\,\textsc{i}$]$ 6300+6363", "%.4f",
     ["0.0080", "0.0060", "0.0063", "0.0059", "0.0047", "0.0050", "-",
      "0.0088", "0.055"], "O I 6300+6363"),
    (r"$[$O\,\textsc{ii}$]$ 7320+7330", "%.4f",
     ["0.0087", "0.0085", "0.0089", "0.0037", "0.0103", "0.0080", "0.0100",
      "0.0064", "0.0089"], "O II 7320+7330"),
    (r"$[$O\,\textsc{ii}$]$ 3726+3729", "%.2f",
     ["1.01", "1.13", "1.10", "1.04", "1.22", "1.08", "1.17", "0.909",
      "1.12"], "O II 3726+3729"),
    (r"$[$O\,\textsc{iii}$]$ 52+88\,$\mu$m", "%.4f",
     ["0.0030", "0.0026", "0.0026", "0.0040", "0.0037", "0.0020", "0.0017",
      "0.0022", "0.0044"], "O III 52+88"),
    (r"$[$O\,\textsc{iii}$]$ 5007+4959", "%.4f",
     ["0.0021", "0.0016", "0.0015", "0.0024", "0.0014", "0.0010", "0.0010",
      "0.0011", "0.0026"], "O III 5007+4959"),
    (r"$[$Ne\,\textsc{ii}$]$ 12.8\,$\mu$m", "%.3f",
     ["0.264", "0.267", "0.276", "0.27", "0.271", "0.286", "0.290", "0.295",
      "0.289"], "Ne II 12.8"),
    (r"$[$S\,\textsc{ii}$]$ 6716+6731", "%.3f",
     ["0.499", "0.473", "0.459", "1.02", "0.555", "0.435", "0.492", "0.486",
      "0.521"], "S II 6716+6731"),
    (r"$[$S\,\textsc{ii}$]$ 4068+4076", "%.3f",
     ["0.022", "0.017", "0.020", "0.052", "0.017", "0.012", "0.015", "0.013",
      "0.017"], "S II 4068+4076"),
    (r"$[$S\,\textsc{iii}$]$ 18.7\,$\mu$m", "%.3f",
     ["0.445", "0.460", "0.441", "0.34", "0.365", "0.398", "0.374", "0.371",
      "0.381"], "S III 18.7"),
    (r"$[$S\,\textsc{iii}$]$ 9532+9069", "%.3f",
     ["0.501", "0.480", "0.465", "0.56", "0.549", "0.604", "0.551", "0.526",
      "0.602"], "S III 9532+9069"),
    (r"$10^{3}\times\Delta(\mathrm{BC}\,3645)$/\AA", None,
     ["5.54", "-", "5.62", "-", "5.57", "5.50", "-", "6.18", "5.52"], None),
    (r"$T_{\mathrm{inner}}$ / K", None,
     ["7224", "6815", "6789", "6607", "6742", "6900", "6708", "6562", "6852"],
     "TINNER"),
    (r"$\langle T[N_{\mathrm{p}}N_{\mathrm{e}}]\rangle$ / K", None,
     ["6680", "6650", "6626", "6662", "6749", "6663", "6679", "6402", "6692"],
     "TNPNE"),
    (r"$R_{\mathrm{out}}/10^{18}$\,cm", "%.2f",
     ["8.89", "8.88", "8.88", "8.7", "8.95", "9.01", "8.92", "8.89", "8.73"],
     "RFRONT"),
    (r"$\langle\mathrm{He^{+}}\rangle/\langle\mathrm{H^{+}}\rangle$", "%.3f",
     ["0.048", "0.051", "0.049", "0.048", "0.044", "0.077", "0.034", "0.041",
      "0.047"], "HEHRATIO"),
]

CAPTION = {
    "hii40": r"""Lexington \HII{} region benchmark HII40 (40~kK blackbody):
the published multi-code results of \citet[][their Table 1]{wood2004}, with a
\textsc{Cloudy} c25.00 calculation and the \MoCHII{} result appended.  Line
intensities are relative to H$\beta$.\label{tab:wood40}""",
    "hii20": r"""Lexington \HII{} region benchmark HII20 (20~kK blackbody):
the published multi-code results of \citet[][their Table 2]{wood2004}, with a
\textsc{Cloudy} c25.00 calculation and the \MoCHII{} result appended.  Line
intensities are relative to H$\beta$.\label{tab:wood20}""",
}

NOTE_COMMON = r"""Codes, in the notation of \citet{wood2004}: GF =
Ferland's \textsc{Cloudy}; HN = Netzer's \textsc{Ion}; DP = P\'equignot's
\textsc{Nebu}; TK = Kallman's \textsc{XStar}; PH = Harrington's code;
RS = Sutherland's \textsc{Mappings}; RR = Rubin's \textsc{Nebula};
BE = Ercolano's \textsc{Mocassin}; WME = the code of \citet{wood2004}
itself.  Columns GF--BE agree entry by entry with Tables 4 and 5 of
\citet{ercolano2003}, the compilation of the published columns.
C25 is \textsc{Cloudy} c25.00 \citep{cloudy2025} run for this paper and
provides a uniform one-dimensional reference alongside the published
multi-code spread.
Its physics input is the unmodified c25.00 test-suite file
\code{tsuite/auto/hii\_paris.in} (HII40) or \code{hii\_coolstar.in} (HII20)
--- \textsc{Cloudy}'s own maintained version of these benchmarks --- with
only output commands appended; line luminosities come from a \code{save line
list absolute} file and are divided by the \textsc{Cloudy} H\,\textsc{i}
4861.32\,\AA{} luminosity, and the \textsc{Cloudy} blends summed for each row
are listed in \code{tables/make\_wood\_tables.py}.  Unlike the published
columns, C25 runs to the stopping criterion of its own deck rather than to a
prescribed outer radius --- an electron fraction of $10^{-2}$ at HII40, the
1000~K temperature floor at HII20 --- which mainly affects the
low-ionization rows fed by the partly neutral edge.
The \MoCHII{} column uses continuous stellar photon energies sampled from
the Planck spectrum with a Sobol sequence (seed 20260728), and evaluates
photoionization cross sections at each sampled energy; the 32-bin ionizing
array is a diagnostic tally only.  It uses the explicit Monte Carlo diffuse
field including the \HeI{} channels (case~A, matching the transfer treatment
of the published codes), and Badnell 2023 total recombination with the
ground-level rate from the Milne relation applied to the same photoionization
cross sections the transport uses, at $2^{25}$ stellar packets with the
volume-integrated convergence criterion; the ionizing luminosity is set from
the benchmark $L({\rm BB})$ so that $Q({\rm H})$ matches the specification;
blends are the sum of the
\MoCHII{} components listed in \code{tables/make\_wood\_tables.py}, and
\MoCHII{} wavelengths are vacuum values.  The \MoCHII{} structural entries
are computed from the same run: $T_{\mathrm{inner}}$ is the mean $\Te$ of the
cells at the illuminated inner edge, $R_{\mathrm{out}}$ the radius at which
the shell-averaged $x_{\mathrm{HII}}$ falls through 0.5, and the helium ratio
follows Eq.~(17) of \citet{ercolano2003},
$\langle\mathrm{He^{+}}\rangle/\langle\mathrm{H^{+}}\rangle =
\int \nelec n(\mathrm{He^{+}})\,{\rm d}V / \int \nelec n_{\mathrm{p}}\,{\rm d}V
\times n_{\mathrm{H}}/n_{\mathrm{He}}$.  Where the ionized gas still fills the
outermost shell of the prescribed sphere, $R_{\mathrm{out}}$ would be set by the
grid rather than predicted, and the entry is left empty.  The Balmer jump is
not computed."""

NOTE = {
    "hii40": NOTE_COMMON + r"""  The H$\beta$ row is headed $10^{36}$
erg\,s$^{-1}$ in the published table; the values listed there are the
$10^{37}$~erg\,s$^{-1}$ values of \citet{ercolano2003}, and the heading is
corrected accordingly here.""",
    "hii20": NOTE_COMMON,
}


def emit(tag, table):
    production = os.path.join(REPO, "results", "continuous_energy", "production")
    lhb, rows = read_lines(os.path.join(production, "%s_continuous_lines.txt" % tag))
    cl = read_cloudy(tag)
    rates = os.path.join(production, "%s_continuous_rates.h5" % tag)
    tnpne = mochii_tnpne(rates)
    tinner, rfront, he_ratio = mochii_structure(rates, tag)
    hb_unit = 1.0e37 if tag == "hii40" else 1.0e36
    out = []
    out.append("%% Generated by tables/make_wood_tables.py -- do not edit.")
    out.append(r"\begin{deluxetable*}{l" + "c" * (len(CODES) + 1) + "}")
    out.append(r"\tablecaption{" + CAPTION[tag] + "}")
    out.append(r"\tablewidth{0pt}")
    out.append(r"\tabletypesize{\footnotesize}")
    out.append(r"\tablehead{\colhead{Quantity} & " +
               " & ".join(r"\colhead{%s}" % c for c in CODES) +
               r" & \colhead{\MoCHII{}}}")
    out.append(r"\startdata")
    for label, fmt, pub, key in table:
        # published columns; the C25 run goes next to GF, its 1990s ancestor
        cells = [c if c != "-" else r"\nodata" for c in pub]
        cells.insert(1, cloudy_cell(label, fmt, key, tag, cl))
        if key is None:
            moc = r"\nodata"
        elif key == "TNPNE":
            moc = "%.0f" % tnpne
        elif key == "TINNER":
            moc = "%.0f" % tinner
        elif key == "RFRONT":
            moc = (r"\nodata" if rfront is None
                   else "%.2f" % (rfront / (1.0e19 if tag == "hii40" else 1.0e18)))
        elif key == "HEHRATIO":
            moc = "%.3f" % he_ratio
        elif key == "HBETA":
            moc = fmt % (lhb / hb_unit)
        else:
            moc = fmt % blend_value(rows, key)
        out.append("%s & %s & %s \\\\" % (label, " & ".join(cells), moc))
    out.append(r"\enddata")
    out.append(r"\tablecomments{" + NOTE[tag] + "}")
    out.append(r"\end{deluxetable*}")
    path = os.path.join(HERE, "wood_%s.tex" % tag)
    open(path, "w").write("\n".join(out) + "\n")
    print("wrote", path)


emit("hii40", HII40)
emit("hii20", HII20)
