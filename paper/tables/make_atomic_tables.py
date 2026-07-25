"""Generate LaTeX tables of the MoCHII fitted atomic data.

Parses the provenance-headed coefficient files in ``data/atomic/`` and emits
two AASTeX deluxetable files that ``mochii.tex`` \\input's:

  tier1_cooling.tex  -- the collisional line-cooling coefficients
                        Lambda_ion(T) = T^{-1/2} sum_i A_i exp(-T_i/T),
                        one row for each (ion, term), from cooling_<ion>.txt.
  tier2_nlevel.tex   -- one summary row for each ion of the n-level
                        models (levels, fitted Upsilon transitions, worst
                        transition fit error), from nlevel_<ion>.txt.

No number is invented: every coefficient and error is read straight from the
files.  Run:  python3 make_atomic_tables.py
"""

import os
import re
import glob

HERE = os.path.dirname(os.path.abspath(__file__))
ATOMIC = os.path.join(HERE, "..", "..", "data", "atomic")

# element print order and roman numerals for the ionization stage
ELEM_ORDER = ["h", "c", "n", "o", "ne", "s", "ar", "mg", "fe", "si", "cl", "ca"]
ELEM_SYM = {"h": "H", "c": "C", "n": "N", "o": "O", "ne": "Ne", "s": "S",
            "ar": "Ar", "mg": "Mg", "fe": "Fe", "si": "Si", "cl": "Cl",
            "ca": "Ca"}
ROMAN = {1: "i", 2: "ii", 3: "iii", 4: "iv", 5: "v", 6: "vi"}


def ion_key(stem):
    """cooling_o_3 / nlevel_o_3 -> ('o', 3)."""
    m = re.match(r"(?:cooling|nlevel)_([a-z]+)_(\d+)", stem)
    return m.group(1), int(m.group(2))


def ion_label(el, stage):
    """LaTeX ion label, e.g. O\\,\\textsc{iii}."""
    return "%s\\,\\textsc{%s}" % (ELEM_SYM[el], ROMAN[stage])


def sort_ions(items):
    return sorted(items, key=lambda t: (ELEM_ORDER.index(t[0]), t[1]))


def sci(x, sig=4):
    """Format a float as $m.mmm\\times10^{e}$ with `sig` significant figures."""
    if x == 0.0:
        return "0"
    s = "%.*e" % (sig - 1, x)
    mant, exp = s.split("e")
    exp = int(exp)
    return r"$%s\times10^{%d}$" % (mant, exp)


# ----------------------------------------------------------------- parse
def parse_cooling(path):
    A, T, maxerr, signif = [], [], None, None
    with open(path) as f:
        lines = f.readlines()
    body = []
    for ln in lines:
        if ln.startswith("#"):
            m = re.search(r"max fit error:\s*([\d.]+)%\s*\(([\d.]+)%\s*where", ln)
            if m:
                maxerr = float(m.group(1))
                signif = float(m.group(2))
            elif re.search(r"max fit error:\s*([\d.]+)%", ln):
                maxerr = float(re.search(r"max fit error:\s*([\d.]+)%", ln).group(1))
        elif ln.strip():
            body.append(ln.strip())
    n = int(body[0])
    for row in body[1:1 + n]:
        a, t = row.split()
        A.append(float(a))
        T.append(float(t))
    return A, T, maxerr, signif


def parse_nlevel(path):
    nlev = nups = worst = None
    with open(path) as f:
        for ln in f:
            if ln.startswith("#"):
                m = re.search(r"worst transition fit error:\s*([\d.]+)%", ln)
                if m:
                    worst = float(m.group(1))
            elif ln.startswith("NLEV"):
                nlev = int(ln.split()[1])
            elif ln.startswith("NUPS"):
                nups = int(ln.split()[1])
    return nlev, nups, worst


def balanced_split(groups):
    """Split a list of (key, rows) into a left and a right half, keeping each
    group (ion) intact and balancing the total row count on the two sides."""
    total = sum(len(g[1]) for g in groups)
    best_k, best_diff, run = 1, None, 0
    for k in range(1, len(groups)):
        run += len(groups[k - 1][1])
        diff = abs(2 * run - total)
        if best_diff is None or diff < best_diff:
            best_diff, best_k = diff, k
    left = [r for g in groups[:best_k] for r in g[1]]
    right = [r for g in groups[best_k:] for r in g[1]]
    return left, right


def pair_rows(left, right, ncell):
    """Emit LaTeX rows placing the left and right cell lists side by side; a
    short side is padded with empty cells."""
    blank = [""] * ncell
    out = []
    n = max(len(left), len(right))
    for i in range(n):
        l = left[i] if i < len(left) else blank
        r = right[i] if i < len(right) else blank
        out.append(" & ".join(list(l) + [""] + list(r)) + " \\\\")
    return out


# ----------------------------------------------------------------- tier 1
cool = {}
for p in glob.glob(os.path.join(ATOMIC, "cooling_*.txt")):
    if "fit_parameters" in p:        # the collected listing, not an ion file
        continue
    el, st = ion_key(os.path.basename(p)[:-4])
    cool[(el, st)] = parse_cooling(p)

# blocks of 5-cell rows (Ion, N, A_i, T_i, err), one block for each ion
t1_groups = []
for (el, st) in sort_ions(cool.keys()):
    A, T, err, sig_err = cool[(el, st)]
    lab = ion_label(el, st)
    block = []
    for i in range(len(A)):
        if i == 0:
            block.append([lab, "%d" % len(A), sci(A[i]), sci(T[i]),
                          "%.2f\\%%" % err])
        else:
            block.append(["", "", sci(A[i]), sci(T[i]), ""])
    t1_groups.append(((el, st), block))

t1_left, t1_right = balanced_split(t1_groups)

# right-aligned header cell (to match the numeric r-columns) and a
# left-aligned header cell (to match the l Ion column); the vrule strut
# reproduces the row spacing that \colhead builds in.
HEADMACRO = (
    "\\providecommand{\\rhead}[1]{\\multicolumn{1}{r}{\\vrule depth 6pt "
    "height 12pt width 0pt #1}}\n"
    "\\providecommand{\\lhead}[1]{\\multicolumn{1}{l}{\\vrule depth 6pt "
    "height 12pt width 0pt #1}}")

head5 = ("\\colhead{Ion} & \\colhead{$N$} & \\colhead{$A_i$} & "
         "\\colhead{$T_i$ (K)} & \\colhead{max err}")

t1 = []
t1.append("% Auto-generated by paper/tables/make_atomic_tables.py -- do not edit by hand.")
t1.append(HEADMACRO)
t1.append("\\startlongtable")
t1.append("\\begin{deluxetable*}{lccccclcccc}")
t1.append("\\tablecaption{Collisional line-cooling coefficients for the "
          "hot thermal loop, $\\Lambda_{\\rm ion}(T)=T^{-1/2}\\sum_{i=1}^{N} "
          "A_i\\,e^{-T_i/T}$ in erg\\,cm$^3$\\,s$^{-1}$ (the volumetric "
          "cooling rate is $n_e n_{\\rm ion}\\Lambda_{\\rm ion}$), "
          "fitted to CHIANTI v11.0.2 over $10^3$--$10^5$~K.  $A_i$ is in "
          "erg\\,cm$^3$\\,s$^{-1}$\\,K$^{1/2}$.  Ions are listed in two "
          "side-by-side groups (left, then right).  The fitted quantity is the "
          "$n_e\\rightarrow0$ limit of the level populations, so $T_i$ are "
          "excitation temperatures.  The error column is the "
          "maximum fit error over the full range; for the few ions where it is "
          "large (Si\\,\\textsc{iii}, the Cl ions, Ca\\,\\textsc{v}, "
          "Ar\\,\\textsc{iii}, C\\,\\textsc{iii}) the error sits at the ends of "
          "the fitted range, where the cooling contribution is negligible, and "
          "the line diagnostics run through the exact $n$-level "
          "solve.\\label{tab:tier1}}")
t1.append("\\tablehead{%s & \\colhead{} & %s}" % (head5, head5))
t1.append("\\startdata")
t1.extend(pair_rows(t1_left, t1_right, 5))
t1.append("\\enddata")
t1.append("\\tablecomments{Coefficients are reproduced to four significant "
          "figures from the provenance-headed files in "
          "\\texttt{data/atomic/}.}")
t1.append("\\end{deluxetable*}")

with open(os.path.join(HERE, "tier1_cooling.tex"), "w") as f:
    f.write("\n".join(t1) + "\n")

# ----------------------------------------------------------------- tier 2
nlev = {}
for p in glob.glob(os.path.join(ATOMIC, "nlevel_*.txt")):
    stem = os.path.basename(p)[:-4]
    if stem.endswith("_full"):
        continue                    # optional expanded models, noted in caption
    el, st = ion_key(stem)
    nlev[(el, st)] = parse_nlevel(p)

t2cells = []
for (el, st) in sort_ions(nlev.keys()):
    nl, nu, w = nlev[(el, st)]
    t2cells.append([ion_label(el, st), "%d" % nl, "%d" % nu, "%.2f\\%%" % w])

# split into a left and right half (13 + 12), keeping element order
nleft = (len(t2cells) + 1) // 2
t2_left, t2_right = t2cells[:nleft], t2cells[nleft:]

head4 = ("\\colhead{Ion} & \\colhead{levels} & \\colhead{transitions} & "
         "\\colhead{worst err}")

t2 = []
t2.append("% Auto-generated by paper/tables/make_atomic_tables.py -- do not edit by hand.")
t2.append(HEADMACRO)
t2.append("\\startlongtable")
t2.append("\\begin{deluxetable*}{lcccclccc}")
t2.append("\\tablecaption{$n$-level diagnostic models: number of levels, "
          "number of fitted effective-collision-strength transitions, and the "
          "worst transition fit error.  The level energies and radiative "
          "$A$-values are from CHIANTI v11.0.2, and $\\Upsilon(T)$ is a "
          "low-order Chebyshev fit in Burgess--Tully scaled space over "
          "$10^3$--$10^5$~K; the full coefficient sets are in the "
          "provenance-headed files \\texttt{data/atomic/nlevel\\_<ion>.txt}."
          "\\label{tab:tier2}}")
t2.append("\\tablehead{%s & \\colhead{} & %s}" % (head4, head4))
t2.append("\\startdata")
t2.extend(pair_rows(t2_left, t2_right, 4))
t2.append("\\enddata")
t2.append("\\tablecomments{Expanded 16/34-level Fe\\,\\textsc{ii}/\\textsc{iii} "
          "models (\\texttt{nlevel\\_fe\\_\\{2,3\\}\\_full.txt}) are available "
          "for iron-line studies.  Ar\\,\\textsc{ii}, Cl\\,\\textsc{v}, and "
          "Ca\\,\\textsc{iii}/\\textsc{iv} have no CHIANTI~v11 level data and "
          "are tracked in the ionization balance without lines.}")
t2.append("\\end{deluxetable*}")

with open(os.path.join(HERE, "tier2_nlevel.tex"), "w") as f:
    f.write("\n".join(t2) + "\n")

print("tier1: %d ions, %d coefficient rows" %
      (len(cool), sum(len(v[0]) for v in cool.values())))
print("tier2: %d ions" % len(nlev))
print("wrote tier1_cooling.tex, tier2_nlevel.tex")


# ------------------------------------------------- tier 1: fit quality only
# Compact companion to the cooling fits: the number of terms and the maximum
# fit error of each ion.  The coefficients themselves are distributed under
# data/atomic/ (see the appendix footnote), so only the fit quality is shown.
err_cells = []
for (el, st) in sort_ions(cool.keys()):
    A, T, err, sig_err = cool[(el, st)]
    err_cells.append([ion_label(el, st), "%d" % len(A), "%.2f\\%%" % err,
                      "%.2f\\%%" % (sig_err if sig_err is not None else err)])

half = (len(err_cells) + 1) // 2
eleft, eright = err_cells[:half], err_cells[half:]

head3 = ("\\colhead{Ion} & \\colhead{$N$} & \\colhead{max err} & "
         "\\colhead{where it cools}")
te = []
te.append("% Auto-generated by paper/tables/make_atomic_tables.py -- do not edit by hand.")
# seven columns: Ion N err | gutter | Ion N err
te.append("\\begin{deluxetable}{lcccclccc}")
te.append("\\tablecaption{Quality of the collisional line-cooling fits: the "
          "number of exponential terms $N$, the maximum fit error of each ion "
          "over the full $10^3$--$10^5$~K range, and the maximum error "
          "restricted to the temperatures that carry the cooling "
          "($\\Lambda > 10^{-3}\\Lambda_{\\max}$).  The two agree for most ions; "
          "where they differ, the large full-range error sits in a temperature "
          "tail in which the ion radiates negligibly.  The fitted coefficients "
          "themselves are distributed with the code "
          "(Appendix~\\ref{app:atomic}).\\label{tab:tier1err}}")
te.append("\\tablehead{%s & \\colhead{} & %s}" % (head3, head3))
te.append("\\startdata")
te.extend(pair_rows(eleft, eright, 4))
te.append("\\enddata")
te.append("\\end{deluxetable}")

with open(os.path.join(HERE, "tier1_maxerr.tex"), "w") as f:
    f.write("\n".join(te) + "\n")

print("tier1 max-error table:", len(err_cells), "ions")
