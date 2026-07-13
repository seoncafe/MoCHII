#!/usr/bin/env python3
"""Element rate-data files for the MoCHII species registry (Stage G2b).

For each registry element (C, N, O, Ne, S) writes
data/atomic/element_<el>.txt with one block for each ionization transition
i <-> i+1:

  PHOTO — Verner, Ferland, Korista & Yakovlev (1996) VFKY96 outer-shell
          photoionization fit parameters, parsed directly from the authors'
          phfit2.f PH2 table (~/RT_Codes/Verner_Fortran_sub/phfit2.f);
          thresholds are the published E_th (NIST ionization potentials).
  CI    — Voronov (1997) collisional-ionization fit, parsed from the
          authors' cfit.f CF table (same directory).
  RR/DR — Badnell (2006) radiative recombination and Badnell-group
          dielectronic recombination total ground-level fits of the
          RECOMBINING ion (stage i+1), parsed from the CHIANTI v11
          .rrparams/.drparams files (verified identical to the EXHALE
          Badnell table for the shared ions).
  CXI/CXR — charge exchange with H+ / H0 (first stage only), from the
          verified transcription of Huang et al. (2023) Table 4
          (EXHALE docs/charge_exchange_table4.md; original sources
          KF96/S98/S99/L05/Z05).

The Fortran species registry parses these files; it never sees CHIANTI or
the Verner sources directly.  Run:  python3 make_element_data.py
"""
import os
import re

DATE = "2026-07-13"
VERNER_DIR = os.path.expanduser("~/RT_Codes/Verner_Fortran_sub")
DBASE = os.path.expanduser("~/RT_Codes/CHIANTI/dbase")
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "..", "data", "atomic")
# Badnell RR/DR tables: data/atomic/badnell_{rr,dr}.dat are the Strathclyde
# AMDPP TAMOC clist_K files (https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR|DR/,
# magic dates 20230511/20230512; verified md5-identical to the copies Cloudy
# c23.01 redistributes).  Primary RR+DR source: the recombining ion of
# transition i is X^i+ (N = Z-i electrons); the M=1 ground-level fit is used.
# Ions absent from these tables (low Ar/Fe/Cl/Ca stages — confirmed absent at
# every M in the original) fall back to the CHIANTI .rrparams/.drparams read.
BADNELL_RR = os.path.join(OUTDIR, "badnell_rr.dat")
BADNELL_DR = os.path.join(OUTDIR, "badnell_dr.dat")

# (element, Z, number of stages tracked; the highest stage is terminal)
ELEMENTS = [("c", 6, 4), ("n", 7, 3), ("o", 8, 3), ("ne", 10, 3),
            ("s", 16, 4), ("ar", 18, 4), ("mg", 12, 3), ("fe", 26, 4),
            ("si", 14, 5), ("cl", 17, 5), ("ca", 20, 5)]

# published VFKY96 thresholds E_th [eV] by (Z, n_electron)
ETH = {(6, 6): 11.260, (6, 5): 24.383, (6, 4): 47.888,
       (7, 7): 14.534, (7, 6): 29.601, (7, 5): 47.449,
       (8, 8): 13.618, (8, 7): 35.121, (8, 6): 54.936,
       (10, 10): 21.564, (10, 9): 40.963, (10, 8): 63.450,
       (16, 16): 10.360, (16, 15): 23.330, (16, 14): 34.830,
       (16, 13): 47.300,
       (18, 18): 15.760, (18, 17): 27.630, (18, 16): 40.740,
       (18, 15): 59.810,
       (12, 12): 7.646, (12, 11): 15.035,
       (26, 26): 7.902, (26, 25): 16.199, (26, 24): 30.651,
       (14, 14): 8.152, (14, 13): 16.346, (14, 12): 33.493,
       (14, 11): 45.142,
       (17, 17): 12.970, (17, 16): 23.810, (17, 15): 39.610,
       (17, 14): 53.470,
       (20, 20): 6.113, (20, 19): 11.872, (20, 18): 50.913,
       (20, 17): 67.273}

# charge exchange with H, first transition: (form, coefficients)
# forms: 1 a*(T/300)^b*exp(-c/T4); 2 a*T4^b + c*T4^d;
#        3 (a*T4^b + c*T4^d)*exp(-e/T); 4 exp(sum c_k lnT^k);
#        5 form4 * exp(-c6/T4); 6 a*t4^b*(1 + c*exp(d*t4))  (KF96);
#        7 a*t4^b*(1 + c*exp(d*t4))*exp(-e/t4)  (KF96 with barrier)
CX = {
    "c":  dict(CXI=(1, [1.31e-15, 0.213, 0.0, 0, 0, 0]),
               CXR=(1, [6.3e-17, 1.96, 17.0, 0, 0, 0])),
    "n":  dict(CXI=(4, [-35.4, 1.94, -0.154, -6.3e-3, -1.16e-3, 0]),
               CXR=(4, [-40.1, 6.4, -1.75, 0.18, -5.96e-3, 0])),
    "o":  dict(CXI=(2, [2.08e-9, 0.405, 1.11e-11, -0.458, 0, 0]),
               CXR=(3, [1.26e-9, 0.517, 4.25e-10, 6.69e-3, 227.0, 0])),
    "ne": dict(CXI=(0, [0]*6), CXR=(0, [0]*6)),
    "s":  dict(CXI=(4, [-50.0, 13.3, -2.77, 0.243, -7.24e-3, 0]),
               CXR=(5, [-50.14, 13.3, -2.77, 0.243, -7.24e-3, 3.76])),
    "ar": dict(CXI=(0, [0]*6), CXR=(0, [0]*6)),   # Ar+ + H negligible (KF96)
    # Mg (Huang Table 4 A1/A2, KF96): fast forward Mg + H+, endothermic
    # reverse with a 6.91e4-K barrier.
    "mg": dict(CXI=(6, [9.76e-12, 3.14, 55.54, -1.12, 0, 0]),
               CXR=(7, [2.95e-12, 3.28, 55.54, -1.12, 6.91, 0])),
    # Fe (Huang A5/A6): Fe + H+ fast constant (RV72); reverse with barrier.
    "fe": dict(CXI=(6, [4.0e-9, 0.0, 0.0, 0.0, 0, 0]),
               CXR=(1, [1.16e-9, 0.072, 6.61, 0, 0, 0])),
    # Si (Huang et al. 2023 A9/A10): Si0 + H+ -> Si+ ionization and
    # Si+ + H0 -> Si0 recombination at the first boundary (A10 = 2.371e-12
    # at 1e4 K).  The old code put the MOCASSIN chex(14,2) = Si2+ + H0 rate
    # here, one boundary too low; that rate now sits at transition 2 below.
    "si": dict(CXI=(1, [7.41e-11, 0.85, 0.0, 0, 0, 0]),
               CXR=(1, [4.71e-11, 0.95, 6.32, 0, 0, 0])),
    # Cl, Ca: no KF96/Huang transcription available — no CX.
    "cl": dict(CXI=(0, [0]*6), CXR=(0, [0]*6)),
    "ca": dict(CXI=(0, [0]*6), CXR=(0, [0]*6)),
}

# charge exchange with H, higher transitions (recombination direction only;
# the ionization direction is endothermic and negligible there).  Kingdon &
# Ferland (1996) fits k = a*1e-9 * t4^b * (1 + c*exp(d*t4)) [cm^3 s^-1],
# transcribed from the MOCASSIN update_mod.f90 chex table.  Fit validity
# ~5e3-5e4 K.
#
# Stage convention (corrected): MOCASSIN chex(Z,ion) is the reaction
# X^(ion+) + H0 -> X^((ion-1)+) (reactant charge = ion; its comment labels
# the PRODUCT ion).  MoCHII transition i recombines X^(i+) -> X^((i-1)+),
# so CX_HIGH[(el, i)] = chex(Z, i).  The earlier table used chex(Z, i+1)
# for C/N/O/Ne/S/Ar/Si (reading MOCASSIN's product-labeled comment as the
# reactant), which put every higher-transition rate one ionization stage
# too low — e.g. C2+ + H0 -> C+ was 3.27e-9 (the C3+ rate) instead of the
# correct 1.04e-12.  Mg/Fe were already right (Huang labels).
CX_HIGH = {
    ("c", 2):  (6, [1.67e-13, 2.79, 304.72, -4.07, 0, 0]),   # chex(6,2)  C2+ + H0 -> C+
    ("c", 3):  (6, [3.25e-9, 0.21, 0.19, -3.29, 0, 0]),      # chex(6,3)  C3+ + H0 -> C2+
    ("n", 2):  (6, [3.05e-10, 0.60, 2.65, -0.93, 0, 0]),     # chex(7,2)  N2+ + H0 -> N+
    ("o", 2):  (6, [1.04e-9, 0.27, 2.02, -5.92, 0, 0]),      # chex(8,2)  O2+ + H0 -> O+
    ("ne", 2): (6, [1.0e-14, 0.0, 0.0, 0.0, 0, 0]),          # chex(10,2) Ne2+ + H0 -> Ne+ (negligible)
    ("s", 2):  (6, [1.0e-14, 0.0, 0.0, 0.0, 0, 0]),          # chex(16,2) S2+ + H0 -> S+ (negligible)
    ("s", 3):  (6, [2.29e-9, 4.02e-2, 1.59, -6.06, 0, 0]),   # chex(16,3) S3+ + H0 -> S2+
    ("ar", 2): (6, [1.0e-14, 0.0, 0.0, 0.0, 0, 0]),          # chex(18,2) Ar2+ + H0 -> Ar+ (negligible)
    ("ar", 3): (6, [4.57e-9, 0.27, -0.18, -1.57, 0, 0]),     # chex(18,3) Ar3+ + H0 -> Ar2+
    # Mg/Fe: already correct (Huang reactant labels land on chex(Z,i)).
    ("mg", 2): (6, [8.58e-14, 2.49e-3, 0.0293, -4.33, 0, 0]),# chex(12,2) Mg2+ + H0 -> Mg+
    ("fe", 2): (6, [1.26e-9, 0.0772, -0.41, -7.31, 0, 0]),   # chex(26,2) Fe2+ + H0 -> Fe+
    ("fe", 3): (6, [3.42e-9, 0.51, -2.06, -8.99, 0, 0]),     # chex(26,3) Fe3+ + H0 -> Fe2+
    ("si", 2): (6, [1.23e-9, 0.24, 3.17, 4.18e-3, 0, 0]),    # chex(14,2) Si2+ + H0 -> Si+
    ("si", 3): (6, [4.90e-10, -8.74e-2, -0.36, -0.79, 0, 0]),# chex(14,3) Si3+ + H0 -> Si2+
    ("si", 4): (6, [7.58e-9, 0.37, 1.06, -4.09, 0, 0]),      # chex(14,4) Si4+ + H0 -> Si3+
}

def parse_badnell_rr():
    """(Z,N) -> [A,B,T0,T1,C,T2] for the M=1 ground level (C,T2 = 0 for the
    4-parameter form)."""
    out = {}
    for ln in open(BADNELL_RR):
        t = ln.split()
        if len(t) < 8 or not t[0].isdigit():
            continue
        Z, N, M = int(t[0]), int(t[1]), int(t[2])
        if M != 1:
            continue
        v = [float(x) for x in t[4:]]
        if len(v) == 4:
            v += [0.0, 0.0]
        out[(Z, N)] = v[:6]
    return out


def parse_badnell_dr():
    """(Z,N) -> ([c_i], [E_i]) for M=1.  The file lists all C rows, then a
    second header, then all E rows."""
    cdict, edict, sec = {}, {}, None
    for ln in open(BADNELL_DR):
        s = ln.split()
        if len(s) >= 5 and s[0] == "Z" and s[1] == "N":
            sec = "C" if s[4].startswith("C") else "E"
            continue
        if len(s) < 5 or not s[0].isdigit():
            continue
        Z, N, M = int(s[0]), int(s[1]), int(s[2])
        if M != 1:
            continue
        (cdict if sec == "C" else edict)[(Z, N)] = [float(x) for x in s[4:]]
    return {k: (cdict[k], edict[k]) for k in cdict
            if k in edict and len(cdict[k]) == len(edict[k])}


def parse_ph2():
    txt = open(os.path.join(VERNER_DIR, "phfit2.f")).read()
    txt = re.sub(r"\n     \S", " ", txt)
    out = {}
    for m in re.finditer(r"DATA \(PH2\(I,\s*(\d+),\s*(\d+)\),I=1,7\)\s*/([^/]+)/",
                         txt):
        nz, ne = int(m.group(1)), int(m.group(2))
        out[(nz, ne)] = [float(x) for x in m.group(3).split(",")]
    return out


MOCASSIN_DATA = os.path.expanduser(
    "~/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/data")
# subshell orbital quantum number by shell index (1s 2s 2p 3s 3p 3d 4s)
VY95_LEVEL = [0, 0, 1, 0, 1, 2, 0]
VY95_NTOT = [1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5,
             6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7]


def parse_ph1():
    """Verner & Yakovlev (1995) subshell parameters from the MOCASSIN
    data/ph1.dat (sequential; the read order replicates phInit): for
    each (nz, ne, shell) six values E_th, E_0, sigma_0, y_a, P, y_w.
    Returns {(nz, ne): (params of the OUTER shell, l of that shell)} —
    the outer-shell VY95 fit used for the elements absent from the
    VFKY96 PH2 table (P, Cl, K; MOCASSIN phFitEl does the same)."""
    toks = open(os.path.join(MOCASSIN_DATA, "ph1.dat")).read().split()
    p = 0
    out = {}
    for j in range(1, 31):              # nz
        for k in range(1, min(j, 30)+1):   # ne
            nt = VY95_NTOT[k-1]
            if j == k and k > 18:
                nt = 7
            if j == k+1 and j in (20, 21, 22, 25, 26):
                nt = 7
            rows = []
            for _ in range(nt):
                rows.append([float(x) for x in toks[p:p+6]])
                p += 6
            nout = nt                     # outer shell = last one read
            out[(j, k)] = (rows[nout-1], VY95_LEVEL[min(nout, 7)-1])
    return out


def parse_cf():
    txt = open(os.path.join(VERNER_DIR, "cfit.f")).read()
    txt = re.sub(r"\n     \S", " ", txt)
    out = {}
    for m in re.finditer(r"DATA\(CF\(I,\s*(\d+),\s*(\d+)\),I=1,5\)\s*/([^/]+)/",
                         txt):
        nz, ne = int(m.group(1)), int(m.group(2))
        out[(nz, ne)] = [float(x) for x in m.group(3).split(",")]
    return out


def read_rrparams(el, stage):
    """RR fit of the recombining ion: ('bad', [A,B,T0,T1,C,T2]) Badnell,
    or ('pl', [A, eta]) the type-3 power law alpha = A (T/1e4)^-eta
    (Shull & Van Steenberg style; used by CHIANTI for e.g. Fe II)."""
    path = os.path.join(DBASE, el, f"{el}_{stage}", f"{el}_{stage}.rrparams")
    with open(path) as fh:
        rrtype = int(fh.readline().split()[0])
        toks = fh.readline().split()
    if rrtype == 1:
        A, B, T0, T1 = (float(t) for t in toks[3:7])
        return ("bad", [A, B, T0, T1, 0.0, 1.0])
    if rrtype == 2:
        A, B, T0, T1, C, T2 = (float(t) for t in toks[3:9])
        return ("bad", [A, B, T0, T1, C, T2])
    if rrtype == 3:
        A, eta = (float(t) for t in toks[2:4])
        return ("pl", [A, eta])
    raise ValueError(f"{path}: unsupported rrparams type {rrtype}")


def read_drparams(el, stage):
    """DR fit of the recombining ion: ('t1', c_i, E_i) Badnell form,
    ('t2', [A, B, T0, T1]) Shull & Van Steenberg (1982) form, or None."""
    path = os.path.join(DBASE, el, f"{el}_{stage}", f"{el}_{stage}.drparams")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        drtype = int(fh.readline().split()[0])
        row1 = [float(t) for t in fh.readline().split()[2:]]
        if drtype == 1:
            row2 = [float(t) for t in fh.readline().split()[2:]]
            E, c = row1, row2
            n = sum(1 for v in c if v != 0.0)
            return ("t1", c[:max(n, 1)], E[:max(n, 1)])
        if drtype == 2:
            return ("t2", row1[:4])
    raise ValueError(f"{path}: unsupported drparams type {drtype}")


def main():
    ph2 = parse_ph2()
    ph1 = parse_ph1()
    cf = parse_cf()
    rr_bad = parse_badnell_rr()
    dr_bad = parse_badnell_dr()
    with open(os.path.join(DBASE, "VERSION")) as fh:
        chianti_version = fh.read().strip()

    for el, Z, nstage in ELEMENTS:
        path = os.path.join(OUTDIR, f"element_{el}.txt")
        with open(path, "w") as fh:
            fh.write(f"# MoCHII element rate data: {el.upper()} "
                     f"(Z={Z}, {nstage} stages, highest terminal)\n")
            fh.write("# PHOTO: E_th[eV] E0 sigma0[Mb] ya P yw y0 y1 "
                     "(VFKY96 outer shell, phfit2.f PH2)\n")
            fh.write("# CI: dE[eV] P A X K (Voronov 1997, cfit.f)\n")
            fh.write("# RR: A B T0 T1 C T2 (Badnell/Strathclyde TAMOC, "
                     "data/atomic/badnell_rr.dat M=1; else CHIANTI v"
                     f"{chianti_version} .rrparams)\n")
            fh.write("# DR: n, c_i [cm^3 s^-1 K^1.5], E_i [K] (Badnell/"
                     "Strathclyde TAMOC, data/atomic/badnell_dr.dat M=1; else "
                     "CHIANTI .drparams; n=0 if absent)\n")
            fh.write("# DR2: A B T0 T1 (Shull & Van Steenberg 1982 form, "
                     "alpha = A T^-3/2 e^-T0/T (1 + B e^-T1/T))\n")
            fh.write("# RR2: A eta (power law alpha = A (T/1e4)^-eta; "
                     "CHIANTI rrparams type 3)\n")
            fh.write("# CXI/CXR: form, 6 coefficients (charge exchange with "
                     "H+ / H0; Huang et al. 2023 Table 4 transcription;\n")
            fh.write("#   forms: 0 none; 1 a(T/300)^b exp(-c/T4); "
                     "2 aT4^b + cT4^d; 3 (aT4^b + cT4^d)exp(-e/T);\n")
            fh.write("#   4 exp(sum c_k lnT^k); 5 form4*exp(-c6/T4))\n")
            fh.write(f"# generated by tools/fitting/make_element_data.py on {DATE}\n")
            fh.write(f"ELEMENT {el} {Z} {nstage}\n")
            for i in range(1, nstage):
                ne_ion = Z - (i - 1)          # electrons of stage i
                fh.write(f"TRANSITION {i}\n")
                fh.write(f"ETH {ETH[(Z, ne_ion)]:.4f}\n")
                if (Z, ne_ion) in ph2:
                    p = ph2[(Z, ne_ion)]
                    fh.write("PHOTO " + " ".join(f"{v:.6e}" for v in
                             [ETH[(Z, ne_ion)]] + p) + "\n")
                else:
                    # element absent from the VFKY96 outer-shell table
                    # (P, Cl, K): the Verner & Yakovlev (1995) outer-shell
                    # fit, as in MOCASSIN phFitEl.
                    row, lsub = ph1[(Z, ne_ion)]
                    eth, e0, s0, ya, pp, yw = row
                    fh.write("PHOTO2 " + " ".join(f"{v:.6e}" for v in
                             [eth, e0, s0, ya, pp, yw,
                              float(lsub)]) + "\n")
                v = cf[(Z, ne_ion)]
                fh.write("CI " + " ".join(f"{x:.6e}" for x in
                         [v[0], v[1], v[2], v[3], v[4]]) + "\n")
                # recombining ion of transition i is X^i+ (N = Z-i electrons).
                Nrec = Z - i
                #--- RR: Badnell 2023 where tabulated, else CHIANTI.
                if (Z, Nrec) in rr_bad:
                    fh.write("RR " + " ".join(f"{x:.6e}" for x in rr_bad[(Z, Nrec)])
                             + "\n")
                else:
                    rr = read_rrparams(el, i + 1)
                    if rr[0] == "bad":
                        fh.write("RR " + " ".join(f"{x:.6e}" for x in rr[1]) + "\n")
                    else:   # power law alpha = A (T/1e4)^-eta
                        fh.write("RR2 " + " ".join(f"{x:.6e}" for x in rr[1]) + "\n")
                #--- DR: Badnell 2023 where tabulated, else CHIANTI.
                if (Z, Nrec) in dr_bad:
                    c, E = dr_bad[(Z, Nrec)]
                    fh.write(f"DR {len(c)} " +
                             " ".join(f"{x:.6e}" for x in c) + "  " +
                             " ".join(f"{x:.6e}" for x in E) + "\n")
                else:
                    dr = read_drparams(el, i + 1)
                    if dr is None:
                        fh.write("DR 0\n")
                    elif dr[0] == "t1":
                        _, c, E = dr
                        fh.write(f"DR {len(c)} " +
                                 " ".join(f"{x:.6e}" for x in c) + "  " +
                                 " ".join(f"{x:.6e}" for x in E) + "\n")
                    else:   # SVS82: alpha = A T^-3/2 e^-T0/T (1 + B e^-T1/T)
                        fh.write("DR2 " +
                                 " ".join(f"{x:.6e}" for x in dr[1]) + "\n")
                if i == 1:
                    for key in ("CXI", "CXR"):
                        form, coef = CX[el][key]
                        fh.write(f"{key} {form} " +
                                 " ".join(f"{x:.6e}" for x in coef) + "\n")
                elif (el, i) in CX_HIGH:
                    form, coef = CX_HIGH[(el, i)]
                    fh.write("CXI 0 0 0 0 0 0 0\n")
                    fh.write(f"CXR {form} " +
                             " ".join(f"{x:.6e}" for x in coef) + "\n")
                else:
                    fh.write("CXI 0 0 0 0 0 0 0\nCXR 0 0 0 0 0 0 0\n")
        print(f"wrote {os.path.relpath(path)}")


if __name__ == "__main__":
    main()
