#!/usr/bin/env python3
"""Write src/hydrogenic_recombination_data.inc: the recombination of the two
hydrogenic ions MoCHII carries, H II -> H I (Z = 1) and He III -> He II
(Z = 2), from ONE nl-resolved radiative cascade:

  alpha_B(T)        captures into n >= 2 (recomb_mod: the case-B rate the
                    ionization balance uses; alpha_A = alpha_B + alpha_1 with
                    alpha_1 the Milne rate on the transport cross section);
  f_2s(T, case)     the share of those captures that ends in 2s, case A and
                    case B (the H I and He II two-photon continua, and the
                    2^2S path of the He III cascade in the diffuse field);
  He II f_n(T)      case A only: the share that leaves by the Lyman line
                    np -> 1s, n = 2..30, with n >= 31 in one bin (the Lyman
                    path of the He III cascade when the Lyman lines escape).
                    In case B the Lyman path is Ly-alpha alone, 1 - f_2s.

What is computed
----------------
An nl-resolved radiative cascade in the zero-density limit (no collisional
l-mixing), for Z = 1 and Z = 2, from Storey & Hummer's hydrogenic routines
(Storey & Hummer 1991, Comput. Phys. Commun. 66, 129; the CPC library file
ABZS, `~/RT_Codes/Storey_Hummer/CPC_code/abzs_v1_0`), used unmodified:

  ACOFZ1   electric-dipole A-values for every nu,lu -> nl,ll,
  PIXZ1    photoionization cross sections sigma_nl(E) of every nl.

  alpha_nl(T)  from sigma_nl by the Milne relation,
      alpha_nl = [g_nl/g_+] sqrt(2/pi) c^-2 (m_e kT)^-3/2
                 Int_chi^inf E^2 sigma_nl(E) exp[-(E - chi)/kT] dE ,
      g_nl = 2(2l+1), g_+ = 1, Simpson in ln E over [chi, chi + 60 kT]
      (401 nodes).  The ionization energies carry the reduced mass.
  cascade      every level passes its effective recombination coefficient
      down in proportion to its A-values, top (n = NMAX) first:
          alpha_eff(nl) = alpha_nl + sum_{n'l'} alpha_eff(n'l') A(n'l'->nl)/A(n'l')
      Case A keeps every decay.  Case B deletes np -> 1s from both the
      branching and the total, so every Lyman-line photon is reabsorbed until
      the cascade leaves by another route.  The case-A Lyman-line rate of np
      is alpha_eff(np) A(np->1s)/A(np).

NMAX = 499, the largest n the routines' factorial arrays (FAL, ALO; 1000)
allow.  The captures above it are added as a tail: a power law fitted to the
n-summed alpha_n at n = 449 and 499 is summed from 500 to infinity (as the
integral from 499.5).  Each cascade output takes the tail in proportion to
its own increment between NMAX = 399 and 499, the share the highest levels
actually pass to it, so the case-A exits still sum to alpha_B (checked).
The tail's share of alpha_B is printed; for H I it is 0.03% at 1e4 K.

Validation (printed at the end): H I alpha_A and alpha_B against the SH95
alpha-tot of the CDS files e1a/e1b at n_e = 100 (0.1% at 1e4 K, where the
Badnell RR fit minus the Milne alpha_1 gave alpha_B 1.1% high); He II
alpha_B against e2b; H I case-B alpha_eff(2s) against Osterbrock & Ferland
(2006) Table 4.11; the cascade's case-A and case-B alpha_eff(H-beta) against
the SH95 emissivities (data/atomic/sh95_hi_case{A,B}.txt; the case-A value
tests the case-A branching directly); He II alpha_eff(2^2S) against SH95.

Why one cascade for everything
------------------------------
f_2s and the Lyman shares are fractions of the alpha_B of the same cascade,
so each closes one calculation, and case A and case B differ only by the one
rule that defines them.  They are relative to alpha_B in both cases: a
capture into 1s never feeds 2s or a Lyman line.

This replaces three older sources, none with a case A:
  alpha_B  Badnell (2023) total RR fit minus the Milne alpha_1: at 1e4 K
        H II 2.6142e-13 against SH95's 2.5850e-13 (+1.13%; +0.7% to +1.2%
        over 5e3-2e4 K), He III +0.35%.  Subtracting a 0.6%-high fit to the
        total from an exact alpha_1 amplified its error in alpha_B;
  H I f_2s  the Pengelly (1964) power law 0.8368e-13 T4^-0.723 (Osterbrock &
        Ferland 2006 Table 4.11, case B), reproduced here to 0.3% over
        1e4-2e4 K.  At 5000 K the table's 1.38e-13 is 4% above this cascade
        (1.326e-13), and so is it above SH95's He II value at 2e4 K taken to
        H I by hydrogenic scaling, alpha_H(T) = alpha_HeII(4T)/2 = 1.324e-13;
  He II f_2s  Storey & Hummer (1995) alpha_eff(2^2S) at n_e = 100 cm^-3, four
        anchors 5000-30000 K clamped outside -- held at the endpoint values
        beyond that range, 24% low at 3000 K and 62% high at 6e4 K -- and
        divided by the SH95 alpha-tot of the same density.  The zero-density
        cascade is below those anchors by 0.92% at 5000 K and 0.11% at 3e4 K;
        that difference is the collisional l-mixing SH95 includes at 100 cm^-3
        and this calculation leaves out.

Needs gfortran and numpy.  Run: python3 make_hydrogenic_recombination.py
"""
import gzip
import os
import subprocess
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.expanduser("~/RT_Codes/Storey_Hummer/CPC_code/abzs_v1_0")
SH95 = os.path.expanduser("~/RT_Codes/Storey_Hummer/data")
OUT = os.path.join(HERE, "..", "..", "src", "hydrogenic_recombination_data.inc")
ATOMIC = os.path.join(HERE, "..", "..", "data", "atomic")

NMAX = 499
NMAX_LOW = 399
NLYO = 30                     # Lyman lines kept one by one; n > NLYO binned
LOGT = np.round(np.arange(1.5, 7.0001, 0.05), 4)
# nuclear masses [u] (the routines form the reduced mass from them)
IONS = [(1, 1.007276467, "H I"), (2, 4.001506179, "He II")]

EV = 1.602176634e-12
H_ERG_S = 6.62607015e-27
C_CM = 2.99792458e10

DRIVER = r"""
      PROGRAM CASCADE
*     Hydrogenic nl-resolved recombination cascade, n_e -> 0.
*     Arguments: Z, AM (nuclear mass, u), NMAX, then log10 T values on stdin.
*     Output per T: log T, alpha_1s, alpha_B (n = 2..NMAX), alpha_2s(A),
*     alpha_2s(B), alpha_eff(4 -> 2) in case A and case B,
*     alpha_n(NMAX-50), alpha_n(NMAX), then the case-A Lyman-line
*     rates np -> 1s for n = 2..NLYO and their sum over n > NLYO
*     [cm^3 s^-1]
      IMPLICIT NONE
      INTEGER NMX,NQ
      PARAMETER (NMX=500,NQ=401)
      DOUBLE PRECISION Z,AM,ALF(NMX,0:NMX-1),AEF(NMX,0:NMX-1,2)
      DOUBLE PRECISION AC(NMX+1,2),ANC,ATOT(0:NMX-1,2),PHE(NQ),WQ(NQ)
      DOUBLE PRECISION PC(NMX+1,NQ),PCP(NMX+1,2,NQ),PNC(NQ)
      DOUBLE PRECISION ATAB(NMX,0:NMX,2)
      DOUBLE PRECISION T,KT,CHI,E,PRE,RYDEV,HART,PI,ME,C,KB,EV,S,UMAX,DU
      DOUBLE PRECISION AB,H42(2),W,REM,LT,AN1,AN2,LY(NMX),LYR
      INTEGER NLYO
      PARAMETER (NLYO=30)
      INTEGER NMAX,N,L,NP,LP,IC,IE,IOPT,IOS
      CHARACTER*32 ARG
      CALL GETARG(1,ARG)
      READ(ARG,*) Z
      CALL GETARG(2,ARG)
      READ(ARG,*) AM
      CALL GETARG(3,ARG)
      READ(ARG,*) NMAX
      PI=3.14159265358979324D0
      ME=9.1093837015D-28
      C=2.99792458D10
      KB=1.380649D-16
      EV=1.602176634D-12
      HART=27.211386245988D0
*     the same reduced-mass factor ACOFZ1 and PIXZ1 apply
      REM=1.D0/(1.D0+1.D0/(1822.889D0*AM))
      RYDEV=0.5D0*HART*REM
      IOPT=2
      N=NMAX
      CALL ACOFZ1(Z,AM,N,1,AC,ANC,NMX+1,IOPT)
      IF(IOPT.LT.0) STOP 'ACOFZ1 initialisation failed'
      IOPT=2
      N=NMAX
      CALL PIXZ1(Z,AM,N,1,PHE,PC,PCP,PNC,NQ,NMX+1,IOPT)
      IF(IOPT.LT.0) STOP 'PIXZ1 initialisation failed'
    1 READ(5,*,IOSTAT=IOS) LT
      IF(IOS.NE.0) STOP
      T=10.D0**LT
      KT=KB*T
      PRE=SQRT(2.D0/PI)/(C*C)/(ME*KT)**1.5D0
      DO 20 N=1,NMAX
         CHI=Z*Z*RYDEV/DBLE(N*N)*EV
         UMAX=LOG((CHI+60.D0*KT)/CHI)
         DU=UMAX/(NQ-1)
         DO 10 IE=1,NQ
            E=CHI*EXP((IE-1)*DU)
            PHE(IE)=(E-CHI)/EV/HART
            WQ(IE)=2.D0*DU/3.D0
            IF(MOD(IE,2).EQ.0) WQ(IE)=4.D0*DU/3.D0
            IF(IE.EQ.1.OR.IE.EQ.NQ) WQ(IE)=DU/3.D0
   10    CONTINUE
         CALL PIXZ1(Z,AM,NQ,N,PHE,PC,PCP,PNC,NQ,NMX+1,0)
         DO 20 L=0,N-1
            S=0.D0
            DO 15 IE=1,NQ
               E=CHI*EXP((IE-1)*DU)
               S=S+WQ(IE)*E**3*PC(L+1,IE)*EXP(-(E-CHI)/KT)
   15       CONTINUE
            ALF(N,L)=PRE*2.D0*(2*L+1)*S
   20 CONTINUE
      DO 30 N=1,NMAX
      DO 30 L=0,N-1
         AEF(N,L,1)=ALF(N,L)
         AEF(N,L,2)=ALF(N,L)
   30 CONTINUE
      H42(1)=0.D0
      DO 29 N=1,NMAX
   29 LY(N)=0.D0
      H42(2)=0.D0
*     IC = 1 case A, IC = 2 case B (np -> 1s deleted)
      DO 60 NP=NMAX,2,-1
         DO 40 N=1,NP-1
            CALL ACOFZ1(Z,AM,NP,N,AC,ANC,NMX+1,0)
*           AC(L,1): NP,L-1 -> N,L ;  AC(L,2): NP,L -> N,L-1
            DO 35 L=1,N
               ATAB(N,L,1)=AC(L,1)
               ATAB(N,L,2)=AC(L,2)
   35       CONTINUE
   40    CONTINUE
         DO 50 LP=0,NP-1
         DO 50 IC=1,2
            ATOT(LP,IC)=0.D0
            DO 45 N=1,NP-1
               IF(LP+1.LE.N-1) ATOT(LP,IC)=ATOT(LP,IC)+ATAB(N,LP+1,1)
               IF(LP.GE.1.AND.LP-1.LE.N-1) THEN
                  IF(.NOT.(IC.EQ.2.AND.N.EQ.1))
     :               ATOT(LP,IC)=ATOT(LP,IC)+ATAB(N,LP,2)
               ENDIF
   45       CONTINUE
   50    CONTINUE
         DO 58 LP=0,NP-1
         DO 58 IC=1,2
            IF(ATOT(LP,IC).LE.0.D0) GOTO 58
            DO 55 N=1,NP-1
               IF(LP+1.LE.N-1) THEN
                  W=AEF(NP,LP,IC)*ATAB(N,LP+1,1)/ATOT(LP,IC)
                  AEF(N,LP+1,IC)=AEF(N,LP+1,IC)+W
                  IF(NP.EQ.4.AND.N.EQ.2) H42(IC)=H42(IC)+W
               ENDIF
               IF(LP.GE.1.AND.LP-1.LE.N-1) THEN
                IF(.NOT.(IC.EQ.2.AND.N.EQ.1)) THEN
                  W=AEF(NP,LP,IC)*ATAB(N,LP,2)/ATOT(LP,IC)
                  AEF(N,LP-1,IC)=AEF(N,LP-1,IC)+W
                  IF(IC.EQ.1.AND.N.EQ.1) LY(NP)=LY(NP)+W
                  IF(NP.EQ.4.AND.N.EQ.2) H42(IC)=H42(IC)+W
                ENDIF
               ENDIF
   55       CONTINUE
   58    CONTINUE
   60 CONTINUE
      AB=0.D0
      DO 70 N=2,NMAX
      DO 70 L=0,N-1
         AB=AB+ALF(N,L)
   70 CONTINUE
      AN1=0.D0
      AN2=0.D0
      DO 75 L=0,NMAX-51
   75 AN1=AN1+ALF(NMAX-50,L)
      DO 76 L=0,NMAX-1
   76 AN2=AN2+ALF(NMAX,L)
      LYR=0.D0
      DO 77 N=NLYO+1,NMAX
   77 LYR=LYR+LY(N)
      WRITE(6,'(F8.4,999ES16.8)') LT,ALF(1,0),AB,AEF(2,0,1),AEF(2,0,2),
     :     H42(1),H42(2),AN1,AN2,(LY(N),N=2,NLYO),LYR
      GOTO 1
      END
"""


def library_source():
    """ACOFZ1 and PIXZ1 from the CPC distribution file, columns 1-72 (the
    card sequence field ABZSnnnn in 73-80 dropped), everything from the
    first subroutine on.  The routines themselves are not edited."""
    lines = [ln.rstrip("\n")[:72] for ln in open(SRC)]
    i0 = next(i for i, ln in enumerate(lines)
              if ln.strip().startswith("SUBROUTINE ACOFZ1"))
    # back up over the comment block that belongs to ACOFZ1
    while i0 > 0 and lines[i0 - 1].startswith("*"):
        i0 -= 1
    body = lines[i0:]
    # stop before the input deck appended to the distribution file
    iend = max(i for i, ln in enumerate(body) if ln.strip() == "END")
    return "\n".join(body[:iend + 1]) + "\n"


def build(tmp):
    with open(os.path.join(tmp, "radz1_lib.f"), "w") as fh:
        fh.write(library_source())
    with open(os.path.join(tmp, "cascade.f"), "w") as fh:
        fh.write(DRIVER)
    exe = os.path.join(tmp, "cascade.x")
    subprocess.run(["gfortran", "-O2", "-std=legacy", "-o", exe,
                    os.path.join(tmp, "cascade.f"),
                    os.path.join(tmp, "radz1_lib.f")], check=True)
    return exe


def run(exe, z, am, nmax, logt):
    inp = "\n".join(f"{t:.4f}" for t in logt) + "\n"
    out = subprocess.run([exe, str(float(z)), repr(am), str(nmax)],
                         input=inp, capture_output=True, text=True,
                         check=True).stdout
    return np.array([[float(v) for v in ln.split()]
                     for ln in out.splitlines() if ln.strip()])


def tail(a_lo, a_hi, n_lo, n_hi):
    """Sum over n > n_hi of a power law through (n_lo, a_lo), (n_hi, a_hi),
    taken as the integral from n_hi + 1/2."""
    p = np.log(a_lo / a_hi) / np.log(n_hi / n_lo)
    if np.any(p <= 1.0):
        raise RuntimeError("alpha_n tail exponent <= 1; raise NMAX")
    return a_hi * n_hi**p * (n_hi + 0.5)**(1.0 - p) / (p - 1.0), p


def sh95_line(fname, label):
    """alpha_eff of a SH95 H I line at n_e = 100: 4 pi j/(n_e n_p) / h nu,
    on the table's temperature grid."""
    lines = open(os.path.join(ATOMIC, fname)).read().splitlines()
    temps = np.array(next(ln for ln in lines if ln.startswith("T ")).split()[1:],
                     dtype=float)
    i = next(i for i, ln in enumerate(lines) if ln.startswith("LINE " + label))
    wl = float(lines[i].split()[4])
    vals = np.array([float(lines[i + 1 + k].split()[0])
                     for k in range(len(temps))])
    return temps, vals / (H_ERG_S * C_CM / (wl * 1e-8))


def sh95_alpha_tot(fname):
    """SH95 alpha-tot [cm^3 s^-1] at n_e = 100 on the CDS file's T grid (the
    block after the emissivities; see make_sh95_lines.py)."""
    toks = gzip.open(os.path.join(SH95, fname), "rt").read().split()
    nt, nd = int(toks[0]), int(toks[1])
    p, temps = 2, []
    for ia in range(nt):
        for ib in range(nd):
            if ib == 0:
                temps.append(float(toks[p + 2]))
            ncut = int(toks[p + 5])
            p += 6 + ncut * (ncut - 1) // 2
    a = np.array(toks[p:p + nt * nd], dtype=float).reshape(nt, nd)
    return np.array(temps), a[:, 0]


def fortran_data(name, values, per_line=4):
    """One DATA statement, continued; values in E notation."""
    v = [f"{x:.10e}".replace("e", "E") + "_wp" for x in values]
    lines = [", ".join(v[i:i + per_line]) for i in range(0, len(v), per_line)]
    out = f"  data {name} / &\n"
    out += ", &\n".join("       " + ln for ln in lines) + " /\n"
    return out


def main():
    rows = {}
    with tempfile.TemporaryDirectory() as tmp:
        exe = build(tmp)
        for z, am, name in IONS:
            hi = run(exe, z, am, NMAX, LOGT)
            lo = run(exe, z, am, NMAX_LOW, LOGT)
            rows[z] = (hi, lo)

    blocks = {}
    for z, am, name in IONS:
        hi, lo = rows[z]
        tb, p = tail(hi[:, 7], hi[:, 8], NMAX - 50, NMAX)
        d_b = hi[:, 2] - lo[:, 2]
        aB = hi[:, 2] + tb

        def with_tail(col):
            return hi[:, col] + tb * (hi[:, col] - lo[:, col]) / d_b

        a2sA, a2sB = with_tail(3), with_tail(4)
        ly = np.array([with_tail(9 + k) for k in range(NLYO)])   # n = 2..NLYO+1
        # ly[0..NLYO-2] are n = 2..NLYO, ly[NLYO-1] the n > NLYO bin
        closeA = (ly.sum(axis=0) + a2sA) / aB - 1.0
        blocks[z] = dict(logT=hi[:, 0], a1s=hi[:, 1], aB=aB, a2sA=a2sA,
                         a2sB=a2sB, ly=ly, closeA=closeA, tailB=tb / aB, p=p,
                         h42A=hi[:, 5], h42B=hi[:, 6])

    nt = len(LOGT)
    with open(OUT, "w") as fh:
        fh.write(
            "! MoCHII: hydrogenic recombination, H II -> H I (Z = 1) and\n"
            "! He III -> He II (Z = 2), from one nl-resolved radiative cascade in\n"
            "! the zero-density limit on the Storey & Hummer (1991, CPC 66, 129)\n"
            f"! routines ACOFZ1/PIXZ1, n <= {NMAX} plus a power-law tail.\n"
            "! GENERATED by tools/fitting/make_hydrogenic_recombination.py --\n"
            "! do not edit; rerun the script.  Included by\n"
            "! hydrogenic_recombination_mod.\n"
            "!   hyd_logT        log10 T [K]\n"
            "!   hyd_logaB       log10 alpha_B [cm^3 s^-1], captures into n >= 2\n"
            "!   hyd_f2s         alpha_eff(2s)/alpha_B, case A (:,1,z), B (:,2,z)\n"
            "!   heii_lyman_fA   He II, case A: share of alpha_B leaving by\n"
            "!                   np -> 1s, n = 2..30 (rows 1..29), n >= 31 (row 30)\n")
        for z, am, name in IONS:
            b = blocks[z]
            fh.write(f"! {name}: tail share of alpha_B {b['tailB'].min():.1e}"
                     f"-{b['tailB'].max():.1e}; case-A exits close on alpha_B "
                     f"to {np.abs(b['closeA']).max():.1e}\n")
        fh.write(f"  integer, parameter :: HYD_NT = {nt}\n")
        fh.write(f"  integer, parameter :: HEII_NLY = {NLYO}\n")
        fh.write("  real(kind=wp), save :: hyd_logT(HYD_NT)\n")
        fh.write("  real(kind=wp), save :: hyd_logaB(HYD_NT, 2)\n")
        fh.write("  real(kind=wp), save :: hyd_f2s(HYD_NT, 2, 2)\n")
        fh.write("  real(kind=wp), save :: heii_lyman_fA(HEII_NLY, HYD_NT)\n")
        fh.write(fortran_data("hyd_logT", LOGT))
        for z, am, name in IONS:
            b = blocks[z]
            fh.write(fortran_data(f"hyd_logaB(:, {z})", np.log10(b["aB"])))
            fh.write(fortran_data(f"hyd_f2s(:, 1, {z})", b["a2sA"] / b["aB"]))
            fh.write(fortran_data(f"hyd_f2s(:, 2, {z})", b["a2sB"] / b["aB"]))
        he = blocks[2]
        for it in range(nt):
            fh.write(fortran_data(f"heii_lyman_fA(:, {it + 1})",
                                  he["ly"][:, it] / he["aB"][it]))
    print(f"wrote {os.path.normpath(OUT)}")

    # ---------------------------------------------------------------- checks
    def at(b, key, T):
        return 10.0 ** np.interp(np.log10(T), b["logT"], np.log10(b[key]))

    h, he = blocks[1], blocks[2]
    for fn, key, label in (("e1b.d.gz", None, "H I alpha_B"),
                           ("e1a.d.gz", None, "H I alpha_A"),
                           ("e2b.d.gz", None, "He II alpha_B")):
        ts, vs = sh95_alpha_tot(fn)
        b = he if fn.startswith("e2") else h
        print(f"\n{label} vs SH95 alpha-tot (n_e = 100)")
        for T, v in zip(ts, vs):
            m = at(b, "aB", T)
            if label == "H I alpha_A":
                m = m + at(b, "a1s", T)
            print(f"  {T:8.0f} K  {m:.4e}  SH95 {v:.4e}  ratio {m/v:.4f}")
    print("\nH I alpha_2s case B vs Osterbrock & Ferland (2006) Table 4.11")
    for T, of in ((5e3, 1.38e-13), (1e4, 0.838e-13), (1.5e4, 0.625e-13),
                  (2e4, 0.506e-13)):
        v = at(h, "a2sB", T)
        print(f"  {T:7.0f} K  {v:.4e}  O&F {of:.3e}  ratio {v/of:.4f}")
    print("\nf_2s (A, B) and He II case-A Lyman shares at 1e4 K")
    for b, name in ((h, "H I"), (he, "He II")):
        print(f"  {name}: f_A {at(b, 'a2sA', 1e4)/at(b, 'aB', 1e4):.4f}  "
              f"f_B {at(b, 'a2sB', 1e4)/at(b, 'aB', 1e4):.4f}")
    i4 = int(np.argmin(np.abs(he["logT"] - 4.0)))
    fl = he["ly"][:, i4] / he["aB"][i4]
    print("  He II case A: Ly-alpha %.4f, Ly-beta %.4f, Ly-gamma %.4f, "
          "n = 5..30 %.4f, n >= 31 %.2e, sum + f_2s - 1 = %.1e"
          % (fl[0], fl[1], fl[2], fl[3:NLYO - 1].sum(), fl[NLYO - 1],
             fl.sum() + he["a2sA"][i4] / he["aB"][i4] - 1.0))
    print("\nH I alpha_eff(H-beta) vs SH95 n_e = 100 (tail not included)")
    for case, key in (("A", "h42A"), ("B", "h42B")):
        ts, vs = sh95_line(f"sh95_hi_case{case}.txt", "Hbeta")
        for T, v in zip(ts, vs):
            if 3e3 <= T <= 3e4:
                m = 10.0 ** np.interp(np.log10(T), h["logT"], np.log10(h[key]))
                print(f"  case {case} {T:7.0f} K  {m:.4e}  SH95 {v:.4e}  "
                      f"ratio {m/v:.4f}")
    print("\nHe II alpha_eff(2^2S) case B vs SH95 (n_e = 100)")
    sh_2s = {5e3: 6.161e-13, 1e4: 4.091e-13, 1.5e4: 3.189e-13, 3e4: 2.035e-13}
    for T in sh_2s:
        a = at(he, "a2sB", T)
        print(f"  {T:7.0f} K  {a:.4e}  SH95 {sh_2s[T]:.4e}  ratio {a/sh_2s[T]:.4f}")
    print("\ntail shares of alpha_B (max over T) and case-A closure:")
    for z, am, name in IONS:
        b = blocks[z]
        print(f"  {name}: {b['tailB'].max():.2e} at log T "
              f"{b['logT'][np.argmax(b['tailB'])]:.2f} (1e4 K: "
              f"{np.interp(4.0, b['logT'], b['tailB']):.2e}); power-law index "
              f"{b['p'].min():.2f}-{b['p'].max():.2f}; closure "
              f"{np.abs(b['closeA']).max():.1e}")


if __name__ == "__main__":
    main()
