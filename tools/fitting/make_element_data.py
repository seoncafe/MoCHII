#!/usr/bin/env python3
"""Element rate-data files for the MoCHII species registry (Stage G2b).

For each registry element (C, N, O, Ne, S, Ar, Mg, Fe, Si, Cl, Ca) writes
data/atomic/element_<el>.txt with one block for each ionization transition
i <-> i+1:

  PHOTO — Verner, Ferland, Korista & Yakovlev (1996) VFKY96 outer-shell
          photoionization fit parameters, parsed directly from the authors'
          phfit2.f PH2 table (~/RT_Codes/Verner_Fortran_sub/phfit2.f);
          ETH and occupied-valence edges use the physical ionization
          potentials of the identified NIST ASD snapshot; inner PH1 edges
          and Voronov fit parameters remain independent.
  PHINNER / PHSHELL — the rest of phfit2's prescription.  The PH2 fit is
          the outer-shell cross section BELOW the first inner-shell edge
          E_inn only; at E >= E_inn phfit2 switches to the sum over every
          occupied subshell of the Verner & Yakovlev (1995) fits (PH1).
          PHINNER gives E_inn (phfit2's einn: the PH1 threshold of shell
          NINN(ne); 0 for the non-OP elements P, Cl, K, Sc-Zn except Fe,
          whose cross section is the PH1 sum at every energy).  PHSHELL
          gives, for each occupied subshell, l E_th E0 sigma0 ya P yw,
          parsed from the same phfit2.f PH1 table (version 2, 25 Mar 1996,
          with the improved inner-shell energies of low ionization stages).
          Without them the PH2 fit was extrapolated through the Mg 2p,
          Ca 3p, Fe 3p and Cl 3s edges inside the transported band, which
          misses the inner-shell cross section entirely (Mg II at 85 eV:
          0.050 Mb from PH2 against 6.68 Mb from the shell sum).
  PHTABLE — n, then n rows "eps[eV] sigma[Mb]": a tabulated R-matrix
          photoionization cross section of the ground term against the
          photoelectron energy eps = E - ETH, which replaces the PHOTO fit
          below its last node E_join = ETH + eps_n (Ne III: Nahar 2019
          BPRM; ground_term_photoionization_table.py, where the choices and
          their validity are written).  The last node is the PHOTO value at
          E_join, so the cross section is continuous there.
  CI    — Voronov (1997) collisional-ionization fit, parsed from the
          authors' cfit.f CF table (same directory).
  RR/DR — Badnell (2006) radiative recombination and Badnell-group
          dielectronic recombination total ground-level fits of the
          RECOMBINING ion (stage i+1), primarily from local Badnell M=1
          tables, with CHIANTI v11 .rrparams/.drparams for absent records.
  CXI/CXR — charge exchange with H+ / H0.  Transition 1 both directions,
          from the verified transcription of Huang et al. (2023) Table 4
          (EXHALE docs/charge_exchange_table4.md; original sources
          KF96/S98/S99/Z05), except nitrogen, which takes Kingdon &
          Ferland (1996) Tables 1 and 3 (see the CX table below); higher
          transitions the recombination direction only (CX_HIGH: KF96, and
          Barragan et al. 2006 for N2+ and O2+).
  CXHEAT — dE_ion dE_rec g_split: the kinetic energy [eV] one event of each
          carried charge-exchange direction gives the gas (the CX_HEAT table
          below, with the product states and their sources).
  CXHEATTAB — n, then n pairs (T [K], dE_ion [eV]): the ionization defect
          of a reaction whose product channels branch with temperature (Mg +
          H+, Allan et al. 1988, extended below 5000 K by the threshold
          Boltzmann factor, BRANCH_T_EXTEND); the Fortran interpolates
          linearly in T and holds the end values; the CXHEAT row then
          carries the 1e4 K value.

The Fortran species registry parses these files; it never sees CHIANTI or
the Verner sources directly.  Run:  python3 make_element_data.py
"""
import math
import os
import re
import argparse
import io
from pathlib import Path
from atomic_output import atomic_write_text
from zhang_fe_recombination import DR_FITS as FE_DR_S1, RR_FITS as FE_RR, verify_final_recombination_source
from charge_exchange_expansion import (kf96_existing_rates, kf96_added_rates,
    kf96_recombination_total_rates, KF99_ADDED_HEAT, KF99_RECOMBINATION_MEAN_HEAT)
from kf96_rate_table import kf96_charge_transfer_table
from ground_term_photoionization_table import ground_term_photoionization_table
from charge_exchange_capture import langevin_recombination_rates, avoided_crossing_mean_heat
LANGEVIN_RECOMBINATION = langevin_recombination_rates()

from datetime import datetime
from zoneinfo import ZoneInfo

DATE = datetime.now(ZoneInfo("Asia/Seoul")).date().isoformat()
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

# (element, Z, number of stages tracked; the highest stage is terminal).
# Stages I-VII for every element (77 stages, 66 transitions; MOCASSIN v2's
# nstage = 7 scope, md/ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md Phase 3):
# every terminal ionization potential then lies above 108 eV (C VII is bare).
ELEMENTS = [("c", 6, 7), ("n", 7, 7), ("o", 8, 7), ("ne", 10, 7),
            ("s", 16, 7), ("ar", 18, 7), ("mg", 12, 7), ("fe", 26, 7),
            ("si", 14, 7), ("cl", 17, 7), ("ca", 20, 7)]

# Physical ionization energies for energy accounting and valence edges.
# Cloudy c25.00 atmdat_adfa.cpp::getEthresh uses the same NIST snapshot for
# the outermost occupied shell; all inner PH1 edges and CI fit energies remain.
from ionization_potentials import read_ionization_potentials
ETH = read_ionization_potentials()

# charge exchange with H, first transition: (form, coefficients)
# forms: 1 a*(T/300)^b*exp(-c/T4); 2 a*T4^b + c*T4^d;
#        3 (a*T4^b + c*T4^d)*exp(-e/T); 4 exp(sum c_k lnT^k);
#        5 form4 * exp(-c6/T4); 6 a*t4^b*(1 + c*exp(d*t4))  (KF96);
#        7 a*t4^b*(1 + c*exp(d*t4))*exp(-e/t4)  (KF96 with barrier);
#        8 tabulated rate coefficient (a CXRTAB row follows), log-log in T
#        9 resolved by product term (CXSTATE rows follow; both directions):
#          the forward rate is the sum over the product terms of the ion,
#          the reverse the detailed balance of each term's channel weighted
#          by that term's actual population, from the n-level atom
CX = {
    "c":  dict(CXI=(1, [1.31e-15, 0.213, 0.0, 0, 0, 0]),
               CXR=(1, [6.3e-17, 1.96, 17.0, 0, 0, 0])),
    # N (2026-09-24): Kingdon & Ferland (1996, ApJS 106, 205), read from the
    # published scan: Table 1 N+ + H -> N + H+, a = 1.01e-3 (1e-9 cm^3 s^-1),
    # b = -0.29, c = -0.92, d = -8.38, 1e2-5e4 K (Butler & Dalgarno 1979);
    # Table 3 N + H+ -> N+ + H, a = 4.55e-3 with the same b, c, d and
    # dE/k = 1.086e4 K.  The two a differ by 4.505, the ground-term weight
    # ratio g(N+) g(H) / (g(N) g(H+)) = 9*2/(4*1) = 4.5: one calculation and
    # its detailed-balance reverse.  Huang et al. (2023) Table 4 prints for
    # this pair two lnT polynomials (citing Lin et al. 2005) under which the
    # endothermic N + H+ FALLS with T, 1.6e-15 (1e3 K) -> 8.9e-20 (1e4 K) ->
    # 3.3e-23 (3e4 K), 2e7 below KF96 at 1e4 K, and N+ + H is 70-86x below;
    # Cloudy c25.00's own transcription of Lin et al. (2005) gives
    # 1.68e-12 / 1.26e-12 at 1e4 K, with KF96.  Lin et al. (2005) itself was
    # not obtained.
    "n":  dict(CXI=(7, [4.55e-12, -0.29, -0.92, -8.38, 1.086, 0]),
               CXR=(6, [1.01e-12, -0.29, -0.92, -8.38, 0, 0])),
    # O: the two directions are EXCHANGED relative to Huang et al. (2023)
    # Table 4, whose rows read "O + H+" and "O+ + H".  Four registry elements
    # have IP(X) > IP(H) -- N, O, Ne, Ar -- but Ne and Ar carry no CX data
    # (form 0) and nitrogen now takes KF96, whose Boltzmann factor sits on
    # N + H+ as it must, so oxygen is the Huang row whose Boltzmann-factor
    # placement had to be checked; it is live because O0 + H+ is
    # near-resonant (13.6181 against 13.5984 eV).  That makes
    # O0 + H+ -> O+ + H0 the ENDOTHERMIC direction, by 227.7 K, and the one
    # that must carry the Boltzmann factor.  Huang's table puts exp(-227/T) on
    # the "O+ + H" row instead: right magnitude, wrong side.  Read as written,
    # the pair violates detailed balance by 1.47x at 8000 K; exchanged, it
    # matches Cloudy c25's independent TableCurve fits
    # (atmdat_char_tran.cpp:112-176, CharExcIonOf/CharExcRecTo[H][O]) to 1-3%
    # over 4000-20000 K, and Cloudy's own pair satisfies detailed balance to 7%.
    # MOCASSIN agrees on the direction too: its chex(8,1) has no Boltzmann
    # factor and its comment labels the PRODUCT ion, i.e. O+ + H0 -> O0 + H+
    # (update_mod.f90:1671).  Same product-versus-reactant labeling that the
    # CX_HIGH note below records for the higher transitions.
    "o":  dict(CXI=(3, [1.26e-9, 0.517, 4.25e-10, 6.69e-3, 227.0, 0]),
               CXR=(2, [2.08e-9, 0.405, 1.11e-11, -0.458, 0, 0])),
    "ne": dict(CXI=(0, [0]*6), CXR=(0, [0]*6)),
    # S (2026-09-26): both directions resolved by the S+ term (form 9) from
    # Zhao et al. (2005) Table III (ZHAO2005 below).  Huang et al. (2023)
    # Table 4 carried the fit of Zhao's TOTAL forward rate and, for S+ + H,
    # "estimate[d] the reverse rates assuming microscopic balance" (their
    # sec. 2.4.2) from that total: a reverse that presumes the S+ 2D and 2P
    # terms populated in LTE.  At nebular densities S+ sits in 4S, and the
    # reverse from 4S is the detailed balance of the 4S channel alone,
    # 2.3e-17 at 1e4 K against the 3.7e-13 of the total -- 1.6e4 x.
    "s":  dict(CXI=(9, [0]*6), CXR=(9, [0]*6)),
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

# Zhao, Stancil, Gu, Liebermann, Funke, Buenker & Kimura (2005, PRA 71,
# 062713; read) Table III, S(3P) + H+ -> S+(term) + H, quantum-mechanical
# molecular-orbital close coupling, rate coefficients [cm^3 s^-1] by product
# term; their "Total" column is the sum.  Terms: the S+ levels of
# nlevel_s_2.txt they comprise (1: 4S3/2; 2, 3: 2D; 4, 5: 2P) and the
# CHIANTI term label, for the observed term energy.
ZHAO2005 = dict(
    T=[20., 40., 60., 80., 100., 200., 400., 600., 800., 1000., 2000.,
       4000., 6000., 8000., 1.0e4, 2.0e4, 4.0e4, 6.0e4, 8.0e4, 1.0e5, 2.0e5,
       4.0e5, 6.0e5, 8.0e5, 1.0e6],
    terms=[("4S", [1], "3s2 3p3 4S",
            [3.10e-16, 3.90e-16, 4.09e-16, 4.21e-16, 4.12e-16, 3.71e-16,
             3.28e-16, 3.07e-16, 2.96e-16, 2.93e-16, 3.11e-16, 4.05e-16,
             5.35e-16, 6.97e-16, 8.90e-16, 8.10e-15, 7.94e-13, 4.34e-12,
             1.03e-11, 1.77e-11, 6.30e-11, 1.59e-10, 2.48e-10, 3.40e-10,
             4.39e-10]),
           ("2D", [2, 3], "3s2 3p3 2D",
            [1.87e-15, 4.33e-15, 7.32e-15, 1.02e-14, 1.28e-14, 2.21e-14,
             3.03e-14, 3.24e-14, 3.33e-14, 3.43e-14, 4.48e-14, 1.23e-13,
             3.22e-13, 6.64e-13, 1.19e-12, 8.65e-12, 5.95e-11, 1.52e-10,
             2.69e-10, 4.02e-10, 1.13e-09, 2.53e-09, 3.76e-09, 4.83e-09,
             5.81e-09]),
           ("2P", [4, 5], "3s2 3p3 2P",
            [3.16e-13, 8.51e-13, 1.20e-12, 1.45e-12, 1.66e-12, 2.44e-12,
             3.28e-12, 3.65e-12, 3.85e-12, 3.98e-12, 4.39e-12, 5.34e-12,
             7.52e-12, 1.11e-11, 1.63e-11, 6.92e-11, 3.13e-10, 6.97e-10,
             1.17e-09, 1.68e-09, 4.40e-09, 9.34e-09, 1.35e-08, 1.70e-08,
             2.01e-08])],
    # the reactant S(3P): the term's statistical weight, the fine-structure
    # levels (<= 0.07 eV) being unresolved in the calculation
    g_react=9.0,
    # printed totals, for the transcription check in cx_state_rows
    total=[3.18e-13, 8.56e-13, 1.21e-12, 1.46e-12, 1.67e-12, 2.47e-12,
           3.31e-12, 3.68e-12, 3.88e-12, 4.01e-12, 4.44e-12, 5.46e-12,
           7.84e-12, 1.18e-11, 1.75e-11, 7.79e-11, 3.74e-10, 8.53e-10,
           1.45e-09, 2.10e-09, 5.59e-09, 1.20e-08, 1.75e-08, 2.22e-08,
           2.64e-08])
CX_STATE = {("s", 1): ("s_2", ZHAO2005)}

# Transition-1 reverse rates that are the detailed balance of the TOTAL
# forward rate.  For Mg+ + H that is stated: Kingdon & Ferland (1996) give
# only Mg + H+, and Huang et al. (2023, sec. 2.4.2) "estimate the reverse
# rates assuming microscopic balance".  The Si+ + H and Fe+ + H rows of their
# Table 4 (Glover & Jappsen 2007; Rutherford & Vroom 1972) carry the whole
# endothermicity as the Boltzmann factor (6.32 and 6.61 T4 against 5.45 and
# 5.70 eV) and have the same form; how they were derived is not stated.  Where the forward reaction ends mostly in excited terms
# of the ion (Mg+ 3p and 4s, Allan et al. 1988; Si+ 4P and Fe+ z4F/z4D by
# excitation energies consistent with KF99 mean heat values), that reverse
# presumes those terms populated
# in LTE, and at nebular densities the ion sits in its ground term: an UPPER
# BOUND on the reverse from the ground term, which needs the forward rate
# into that term and is not available.  Sulfur is resolved by term (form 9).
CX_LTE_UPPER = {"mg": "Huang et al. 2023 from the Kingdon & Ferland 1996 forward",
                "si": "form of Huang et al. 2023 Table 4",
                "fe": "form of Huang et al. 2023 Table 4"}


def cx_state_rows(el, it, eth):
    """The CXSTATE rows of a form-9 transition: header (number of terms, of
    temperatures, reactant weight, dIP = IP(H) - IP(X) [eV], the energy the
    forward reaction into the product ground term releases), the T grid,
    and one row per product term: label, its first and last level in the
    n-level atom, g, observed term energy [eV], rate at each T."""
    ion, tab = CX_STATE[(el, it)]
    for i, tot in enumerate(tab["total"]):
        s = sum(t[3][i] for t in tab["terms"])
        if abs(s/tot - 1.0) > 0.01:
            raise ValueError(f"{el} CX_STATE: terms sum to {s:.3e} at "
                             f"T = {tab['T'][i]}, printed total {tot:.3e}")
    rows = [f"CXSTATE {len(tab['terms'])} {len(tab['T'])} "
            f"{tab['g_react']:.1f} {IP_H - eth:.6f}",
            "CXSTATET " + " ".join(f"{x:.6e}" for x in tab["T"])]
    for label, levs, term, k in tab["terms"]:
        g = sum(2.0*j + 1.0 for j in term_J(ion, levs))
        E = 0.0 if levs == [1] else chianti_term_energy(ion, levs, term)
        rows.append(f"CXSTATEK {label} {levs[0]} {levs[-1]} {g:.1f} {E:.6f} " +
                    " ".join(f"{x:.6e}" for x in k))
    return rows


def term_J(ion, levels):
    """J of the given CHIANTI levels of ion."""
    el = ion.split("_")[0]
    path = os.path.join(DBASE, el, ion, f"{ion}.elvlc")
    from chianti_cooling import read_elvlc
    parsed = read_elvlc(path)
    return [(parsed[i][0]-1.0)/2.0 for i in levels]

# charge exchange with H, higher transitions (recombination direction only;
# the ionization direction is endothermic and negligible there).  Kingdon &
# Ferland (1996) fits k = a*1e-9 * t4^b * (1 + c*exp(d*t4)) [cm^3 s^-1],
# transcribed from the MOCASSIN update_mod.f90 chex table.  Each KF96 row
# has its own published temperature interval; the form-8 CXRTAB written below
# evaluates the fit over that interval and holds the endpoint rates outside.
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
    # N2+ and O2+: Barragan et al. (2006) close-coupling rates (BARRAGAN
    # below), replacing the KF96 fits (3.05e-10 0.60 2.65 -0.93 and
    # 1.04e-9 0.27 2.02 -5.92), which are 0.64x (N2+) and 1.22x (O2+) the
    # published k4/k1 at 1e4 K and 1.56x (O2+) at 5000 K.
    ("n", 2):  (8, "k4"),                                     # N2+(2P) + H0 -> N+ + H+
    ("o", 2):  (8, "k1"),                                     # O2+(3P) + H0 -> O+ + H+
    ("ne", 2): (6, [1.0e-14, 0.0, 0.0, 0.0, 0, 0]),          # chex(10,2) Ne2+ + H0 -> Ne+ (negligible)
    ("s", 2):  (6, [1.0e-14, 0.0, 0.0, 0.0, 0, 0]),          # chex(16,2) S2+ + H0 -> S+ (negligible)
    ("s", 3):  (6, [2.29e-9, 4.02e-2, 1.59, -6.06, 0, 0]),   # chex(16,3) S3+ + H0 -> S2+
    ("ar", 2): (6, [1.0e-14, 0.0, 0.0, 0.0, 0, 0]),          # chex(18,2) Ar2+ + H0 -> Ar+ (negligible)
    ("ar", 3): (6, [4.57e-9, 0.27, -0.18, -1.57, 0, 0]),     # chex(18,3) Ar3+ + H0 -> Ar2+
    # Mg/Fe: already correct (Huang reactant labels land on chex(Z,i)).
    ("mg", 2): (6, [8.58e-14, 2.49e-3, 0.0293, -4.33, 0, 0]),# chex(12,2) Mg2+ + H0 -> Mg+
    ("fe", 2): (6, [1.26e-9, 0.0772, -0.41, -7.31, 0, 0]),   # chex(26,2) Fe2+ + H0 -> Fe+
    ("fe", 3): (6, [3.42e-9, 0.51, -2.06, -8.99, 0, 0]),     # chex(26,3) Fe3+ + H0 -> Fe2+
    ("si", 2): (6, [6.77e-9, 0.0736, -0.43, -0.11, 0, 0]),   # Si2+ + H0 -> Si+ (updated KF96 database, Ferland et al. 2013)
    ("si", 3): (6, [4.90e-10, -8.74e-2, -0.36, -0.79, 0, 0]),# chex(14,3) Si3+ + H0 -> Si2+
    ("si", 4): (6, [7.58e-9, 0.37, 1.06, -4.09, 0, 0]),      # chex(14,4) Si4+ + H0 -> Si3+
}

# Added recombination directions: published KF96 total fits. Their constant
# mean heat is supplied independently by KF99; no excited-product reverse is
# inferred. Keep production stage counts unchanged until full atomic review.
for _key, _reaction in {**kf96_added_rates(), **kf96_recombination_total_rates()}.items():
    _a,_b,_c,_d = _reaction.coefficients
    CX_HIGH[_key] = (6, [_a*1e-9,_b,_c,_d,0,0])

# Barragan, Errea, Mendez, Rabadan & Riera (2006), ApJ 636, 544: electron
# capture O2+(2s2 2p2 3P) + H(1s) -> O+ + H+ (their reaction 1, Table 1 column
# k1) and N2+(2s2 2p 2Po) + H(1s) -> N+ + H+ (reaction 4, Table 2 column k4),
# molecular close-coupling on multireference-CI wavefunctions, 1e2-1e5 K;
# units 1e-9 cm^3 s^-1, transcribed from the journal PDF (checked
# 2026-09-24).  Both ions in their GROUND terms, which is what the registry
# stage fraction counts.  The same calculation is what Cloudy c25.00 fits
# (atmdat_char_tran.cpp:336-362) and CHIANTI v11 adopts.  The authors note
# that radiative electron capture, not included, may matter at the lowest
# temperatures; the table is evaluated clamped to its range.
BARRAGAN = dict(
    T=[1.0e2, 2.0e2, 3.0e2, 4.0e2, 5.0e2, 1.0e3, 1.5e3, 2.0e3, 2.5e3, 3.0e3,
       4.0e3, 5.0e3, 6.0e3, 7.0e3, 8.0e3, 9.0e3, 1.0e4, 2.0e4, 3.0e4, 4.0e4,
       5.0e4, 6.0e4, 7.0e4, 8.0e4, 9.0e4, 1.0e5],
    k1=[0.53362, 0.52153, 0.49852, 0.48016, 0.46613, 0.43607, 0.43975,
        0.45679, 0.47931, 0.50418, 0.55659, 0.60995, 0.66288, 0.71464,
        0.76481, 0.81319, 0.85973, 1.24359, 1.53854, 1.78217, 1.98734,
        2.16168, 2.31143, 2.44169, 2.55653, 2.65907],
    k4=[0.46312, 0.59312, 0.65966, 0.70483, 0.73830, 0.82782, 0.86587,
        0.88640, 0.89963, 0.90935, 0.92378, 0.93484, 0.94395, 0.95174,
        0.95857, 0.96469, 0.97027, 1.01334, 1.06250, 1.13596, 1.22873,
        1.33023, 1.43343, 1.53502, 1.63390, 1.73000])


# ---------------------------------------------------------------------------
# Charge-exchange heating: the kinetic energy each carried reaction gives the
# gas (CXHEAT rows).  One rule for every row:
#
#     dE_rec = [IP(X^it) - IP(H)] - E_rad     X^(it+1) + H0  -> X^it + H+
#     dE_ion = [IP(H) - IP(X^it)] - E_rad     X^it     + H+  -> X^(it+1) + H0
#
# with IP(X^it) the registry's own threshold ETH of transition it and E_rad the
# excitation of the product state that leaves the gas as radiation (an allowed
# or intercombination decay, a forbidden line below its critical density, or
# the photon of a radiative charge transfer).  A reaction into ground-term
# products has E_rad = 0; an endothermic one has dE < 0 and takes that energy
# from the gas.  Kingdon & Ferland (1996, ApJS 106, 205; KF96, read from the
# published scan) state the convention (Sec. 2: X^q + H -> X^(q-1) + H+ + dE,
# dE negative for endothermic reactions, Table 4) and that the recombined ion
# "can and usually will be in an excited state". KF96 Table 1 has no heat column.
# Cloudy c25.00 adds KF99 Table 1 recombination and Table 2 ionization heat in
# data/ctrecombdata.dat column 7 and data/ctiondata.dat column 8, respectively.
# The source/atmdat_char_tran.cpp reader stores recombination rows as
# CTRecombData[nelem][ion] = X^(ion+1) + H -> X^ion + H+; row mapping
# verified here by matching the a, b, c, d columns against KF96 Table 1 and
# this file's CX / CX_HIGH rows. Where a heat value equals
# dIP minus a level energy it suggests a candidate product; where it does not, it
# is used as distributed and the product is recorded as not identified.
# Level energies are CHIANTI observed energies (.elvlc column 6), averaged
# over the J levels of the term with weights 2J+1.
#
# Entries: ("level", element_ion, [level indices], term) -> E_rad = that term;
#          ("ground",) -> E_rad = 0;  ("radiative",) -> E_rad = dIP;
#          ("branch", [(fraction, spec), ...]) -> fraction-weighted E_rad;
#          ("defect", dE, note) -> dE as given (E_rad = dIP - dE).
# "split" > 0 marks a pair whose ionization rate includes radiative charge
# transfer; only its non-radiative part, the detailed-balance reverse of the
# recombination rate, heats (species_mod::metal_charge_exchange_heating), and
# the value is that reverse's statistical-weight ratio
# g(X^(it+1)) g(H0) / [g(X^it) g(H+)].
IP_H = 13.598434          # eV, NIST ASD
KB_EV = 8.617333262e-5    # Boltzmann constant [eV/K], CODATA 2018
HC_EV_CM = 1.239841984e-4  # eV per cm^-1 (CODATA 2018)

CX_HEAT = {
    # C + H+ -> C+ + H: the carried rate (Huang et al. 2023 Table 4, citing
    # Stancil et al. 1998, not read) is dominated by radiative charge
    # transfer; KF96 Sec. 3.1 (read) "assumed simple radiative charge
    # transfer for the reverse" of C+ + H.  Its non-radiative part, the
    # detailed-balance reverse of C+ + H -> C(3P) + H+, ends in C+ 2P and
    # releases the whole dIP; g(C+ 2P) g(H) / [g(C 3P) g(H+)] = 6*2/9.
    ("c", 1): dict(ion=("ground",), rec=("ground",), split=12.0/9.0,
                   note="C+H+ radiative except the detailed-balance reverse "
                        "of C+ + H (KF96 Sec. 3.1); C+ + H -> C(3P) + H+"),
    ("c", 2): dict(rec=("defect", 4.01,
                        "KF99 published mean heat (Cloudy c25 ctrecombdata); "
                        "dIP - 4.01 = 6.775 eV matches no C II term "
                        "(2s2p2 4P 5.33, 2D 9.29 eV): product not identified"),
                   note="C2+ + H -> C+ (KF96 Table 1, Butler et al. 1980)"),
    ("c", 3): dict(rec=("defect", 5.73,
                        "KF99 published mean heat 5.73 = 2.86 (q-1), the "
                        "Landau-Zener crossing-window value; product not "
                        "identified"),
                   note="C3+ + H -> C2+ (KF96 Table 1, Herrero et al. 1995b)"),
    # N: N 2D lies 2.38 eV up, above the 0.936 eV defect, so both directions
    # connect the ground terms by energetics; KF96 Tables 1/3 are one
    # calculation and its detailed-balance reverse.
    ("n", 1): dict(ion=("ground",), rec=("ground",),
                   note="N(4S) + H+ <-> N+(3P) + H, ground terms by "
                        "energetics (KF96 Tables 1 and 3)"),
    # N2+ + H: Barragan et al. (2006) Sec. 3 (read): their rate "agree[s]
    # with the values ... evaluated from the cross sections for populating
    # N+ (2s2p3 3Do) of Bienstock et al. (1986)"; it decays by the allowed
    # 1085 A multiplet. KF99 mean heat 4.56 eV is consistent with this channel;
    # mean heat alone does not establish a resolved product population.
    ("n", 2): dict(rec=("level", "n_2", [7, 8, 9], "2s.2p3 3D"),
                   note="N2+ + H -> N+(2s2p3 3Do) + H+ (Barragan et al. 2006)"),
    # O: ground terms both ways; the carried pair encodes exactly this dIP in
    # its 227 K Boltzmann factor.  Fine-structure branching of O(3P_J)
    # (Stancil et al. 1999, not read) is not resolved: |E_rad| <= 0.028 eV.
    ("o", 1): dict(ion=("ground",), rec=("ground",),
                   note="O(3P) + H+ <-> O+(4S) + H, ground terms "
                        "(Stancil et al. 1999 via Huang et al. 2023)"),
    # O2+ + H: Barragan et al. (2006) name no product, but place the avoided
    # crossing of the quartet subsystem at R0 = 4.57 a.u.; by KF96 eqs.
    # (3)-(4) that is an asymptotic defect of about 6.5 eV, the O+ 2s2p4 4P
    # channel (6.65 eV), consistent with KF99 mean heat 6.65 eV.
    # This consistency alone does not establish its branching. 4P decays by the allowed 833 A multiplet, above 13.6 eV: that
    # photon would photoionize hydrogen, a diffuse source not transported.
    ("o", 2): dict(rec=("level", "o_2", [6, 7, 8], "2s 2p4 4P"),
                   note="O2+ + H -> O+(2s2p4 4P) + H+ (Barragan et al. 2006; "
                        "KF99 published mean heat 6.65 eV)"),
    # 1e-14 cm^3 s^-1 rows: KF96 Sec. 3.1 "When a reaction has no
    # significant channels, it may proceed via radiative charge transfer.
    # The rates for these cases have all been set to 1 x 10^-14": the photon
    # carries dIP.  (Cloudy's table gives these rows -dIP, which its heating
    # sum books as cooling; that sign is not explained by KF96.)
    ("ne", 2): dict(rec=("radiative",),
                    note="Ne2+ + H -> Ne+ + H+ + photon (KF96 1e-14 row)"),
    # S + H+: Zhao et al. (2005, PRA 71, 062713) Table III (read), S(3P) +
    # H+ at 1e4 K: S+(4S) 8.90e-16, S+(2D) 1.19e-12, S+(2P) 1.63e-11 cm^3
    # s^-1.  2D and 2P radiate [S II] 6716/6731 and 4069/4076 below their
    # critical densities (~1e3 and ~1e6 cm^-3).  The 2P share is 0.978 at
    # 4000 K and 0.888 at 2e4 K, so E_rad moves by less than 0.05 eV over
    # 4e3-2e4 K; the 1e4 K branching is used.  The reverse, S+ + H (Huang's
    # row, exp(-3.76/T4)), starts from the S+ ground term.
    ("s", 1): dict(ion=("branch", [(8.90e-16, ("ground",)),
                                   (1.19e-12, ("level", "s_2", [2, 3],
                                               "3s2 3p3 2D")),
                                   (1.63e-11, ("level", "s_2", [4, 5],
                                               "3s2 3p3 2P"))]),
                   rec=("ground",),
                   note="S(3P) + H+ -> S+(4S, 2D, 2P) + H (Zhao et al. 2005 "
                        "Table III, 1e4 K); S+(4S) + H -> S(3P) + H+"),
    ("s", 2): dict(rec=("radiative",),
                   note="S2+ + H -> S+ + H+ + photon (KF96 1e-14 row)"),
    ("s", 3): dict(rec=("defect", 5.73,
                        "KF99 published mean heat 5.73 = 2.86 (q-1); product "
                        "not identified"),
                   note="S3+ + H -> S2+ (KF96 Table 1, Butler & Dalgarno "
                        "1980a)"),
    ("ar", 2): dict(rec=("radiative",),
                    note="Ar2+ + H -> Ar+ + H+ + photon (KF96 1e-14 row)"),
    ("ar", 3): dict(rec=("defect", 5.73,
                         "KF99 published mean heat 5.73 = 2.86 (q-1); product "
                         "not identified"),
                    note="Ar3+ + H -> Ar2+ (KF96 Table 1, Butler & Dalgarno "
                         "1980a)"),
    # Mg + H+ -> Mg+ + H: Allan, Clegg, Dickinson & Flower (1988, MNRAS 235,
    # 1245; read) Table 1 splits the rate between the 2 1Sigma channel, which
    # ends in Mg+ 3p 2P (exothermic, dIP - E(3p) = +1.52 eV), and the 4 1Sigma
    # channel, which ends in Mg+ 4s 2S, "closed for collision energies less
    # than 2.7 eV" (endothermic by E(4s) - dIP = 2.70 eV; the 4s level cascades
    # to 3p and both radiate).  At 1e4 K the 4s channel carries 66% of the
    # rate, so the reaction COOLS the gas by 1.26 eV per event on average;
    # the KF99 published mean heat (+1.52 eV) is the 3p channel alone.  The
    # printed "total" at 26000 K (7.954) is a misprint of 1.504 + 6.350 =
    # 7.854; only the two channel columns are used.  [T K, 2 1Sigma, 4 1Sigma]
    # in 1e-10 cm^3 s^-1:
    ("mg", 1): dict(ion=("branch_T", [
                        (5000., 0.2549, 0.07057), (6000., 0.3430, 0.1856),
                        (8000., 0.5055, 0.6125), (10000., 0.6419, 1.236),
                        (12000., 0.7582, 1.955), (14000., 0.8632, 2.696),
                        (15000., 0.9138, 3.060), (16000., 0.9641, 3.415),
                        (18000., 1.066, 4.093), (20000., 1.170, 4.723),
                        (22000., 1.278, 5.307), (24000., 1.390, 5.847),
                        (26000., 1.504, 6.350), (28000., 1.620, 6.819),
                        (30000., 1.736, 7.259), (32000., 1.850, 7.673)],
                        [("level", "mg_2", [2, 3], "3p 2P"),
                         ("level", "mg_2", [4], "4s 2S")]),
                    rec=("ground",),
                    note="Mg(1S) + H+ -> Mg+(3p 2P or 4s 2S) + H, branched by "
                         "Allan et al. 1988 Table 1 (CXHEATTAB); Mg+(3s) + H "
                         "-> Mg + H+"),
    # Mg2+ + H: Mg+ 3p (4.43 eV) lies above the 1.44 eV defect, so the
    # product is the ground term by energetics.
    ("mg", 2): dict(rec=("ground",),
                    note="Mg2+ + H -> Mg+(3s) + H+, ground by energetics"),
    # Fe + H+: Rutherford & Vroom (1972, JCP 57, 3091; read): "the state of
    # excitation of the products could not be determined" and the large
    # cross sections are attributed to near resonance.  DERIVED here: the
    # near-resonant Fe+ terms are 3d6(5D)4p z4F and z4D (their mean 5.545
    # eV, 0.15 eV below the defect); they decay by allowed lines to a4F/a4D,
    # whose [Fe II] lines radiate below their critical densities, so E_rad is
    # the z4F + z4D mean.  Not a determination of the product.
    # (the two terms weighted by their statistical weights, 20 and 28)
    ("fe", 1): dict(ion=("branch", [(20.0, ("level", "fe_2", [50, 51, 52, 53],
                                            "3d6(5D).4p 4D")),
                                    (28.0, ("level", "fe_2", [54, 55, 56, 57],
                                            "3d6(5D).4p 4F"))]),
                    rec=("ground",),
                    note="Fe + H+ -> Fe+(z4F, z4D; near-resonant, derived from "
                         "Rutherford & Vroom 1972) + H; Fe+(a6D) + H -> Fe + H+"),
    ("fe", 2): dict(rec=("ground",),
                    note="Fe2+ + H -> Fe+(a6D) + H+ (KF96 Table 1, Neufeld & "
                         "Dalgarno 1987; electronic defect 2.56 eV = the "
                         "ground-term dIP to 0.04 eV)"),
    ("fe", 3): dict(rec=("defect", 6.3,
                         "KF99 published mean heat 6.3 eV; product not "
                         "identified"),
                    note="Fe3+ + H -> Fe2+ (KF96 Table 1, Neufeld 1984)"),
    # Si + H+: Glover & Jappsen (2007) R42, which Huang et al. cite, takes
    # the rate from KF96; the KF99 published mean heat 0.12 eV = dIP - E(Si+
    # 3s3p2 4P), which decays by the Si II] 2335 A intercombination line.
    ("si", 1): dict(ion=("level", "si_2", [3, 4, 5], "3s.3p2 4P"),
                    rec=("ground",),
                    note="Si(3P) + H+ -> Si+(3s3p2 4P) + H (product-term "
                         "assignment inferred from Cloudy 0.12 eV heat); "
                         "Si+(2P) + H -> Si + H+"),
    ("si", 2): dict(rec=("ground",),
                    note="Si2+ + H -> Si+(2P) + H+, ground by energetics "
                         "(Si+ 4P at 5.33 eV > dIP)"),
    ("si", 3): dict(rec=("defect", 4.23,
                         "KF99 published mean heat 4.23 eV; dIP - 4.23 = 15.67 "
                         "eV matches no Si III term: product not identified"),
                    note="Si3+ + H -> Si2+ (KF96 Table 1, Butler & Dalgarno "
                         "1980a)"),
    ("si", 4): dict(rec=("level", "si_4", [6], "4s 2S"),
                    note="Si4+ + H -> Si3+(4s 2S) + H+ (product-term "
                         "assignment inferred from Cloudy 7.49 eV heat)"),
}


def chianti_term_energy(ion, levels, term):
    """(2J+1)-weighted mean observed energy [eV] of the given level indices of
    CHIANTI ion (e.g. 'o_2'), checking that each carries the expected term."""
    el = ion.split("_")[0]
    path = os.path.join(DBASE, el, ion, f"{ion}.elvlc")
    want_cfg, want_term = term.rsplit(" ", 1)
    from chianti_cooling import read_elvlc
    parsed = read_elvlc(path)
    rows = {}
    for line in Path(path).read_text().splitlines():
        if line.strip() in ('-1',''):
            break
        rows[int(line[:7])] = line
    sw, se = 0.0, 0.0
    for i in levels:
        line = rows[i]
        actual_term = str(int(line[42:47]))+line[47:52].strip()
        if actual_term != want_term:
            raise ValueError(f"{ion} level {i}: term {actual_term} is not {want_term}")
        have = line[7:37].replace('.', ' ').split()
        if any(c not in have for c in want_cfg.replace('.', ' ').split()):
            raise ValueError(f"{ion} level {i}: configuration {have} is not {want_cfg}")
        g, E, theoretical = parsed[i]
        if E <= 0.0 and i != 1:
            raise ValueError(f"{ion} level {i}: no observed energy")
        sw += g
        se += g*E
    return se/sw*HC_EV_CM


for _key, _energy in {**KF99_ADDED_HEAT, **KF99_RECOMBINATION_MEAN_HEAT}.items():
    CX_HEAT[_key] = dict(
        rec=('defect',_energy,'KF99 Table 1 mean kinetic heat'),
        note='KF96 total recombination rate; KF99 Table 1 mean heat; '
             'ground-state reactants, low-density radiative relaxation; '
             'product branching and diffuse cascade reabsorption not modeled')
    if _key in (('cl',3),('cl',4),('fe',4)):
        CX_HEAT[_key]['note'] += '; KF96 Table 2 Landau-Zener total rate estimate'
    if _key == ('fe',4):
        CX_HEAT[_key]['note'] += '; radiative charge transfer may dominate this estimated channel'
    if _key == ('n',4):
        CX_HEAT[_key]['note'] += '; rate from the updated KF96 database (Ferland et al. 2013)'

# User-adopted q>=5 unit-transfer capture approximation, distinct from KF96 fits.
# Form 1 with zero exponents is temperature independent, without fit clamping.
for _key, _rate in LANGEVIN_RECOMBINATION.items():
    CX_HIGH[_key]=(1,[_rate,0,0,0,0,0])
    CX_HEAT[_key]=dict(
        rec=('defect',avoided_crossing_mean_heat(_key[1]),'KF99 avoided crossing with exit charge: 2.86*(q-1) eV'),
        note='Ferland 1997 Langevin capture upper limit with reduced collision mass; '
             'unit electron-transfer probability assumed; no certified temperature interval; '
             'KF99 9.5-bohr avoided-crossing model, exit-charge mean kinetic heat 2.86*(q-1) eV; ground-state reactants and '
             'low-density radiative relaxation; reverse rates, product branching and '
             'cascade-photon reabsorption not determined')


def radiated_energy(spec, dip):
    """E_rad [eV] of one CX_HEAT entry for a reaction of energy defect dip."""
    kind = spec[0]
    if kind == "ground":
        return 0.0
    if kind == "radiative":
        return dip
    if kind == "level":
        return chianti_term_energy(spec[1], spec[2], spec[3])
    if kind == "branch":
        tot = sum(f for f, _ in spec[1])
        return sum(f*radiated_energy(s, dip) for f, s in spec[1])/tot
    if kind == "defect":
        return dip - spec[1]
    if kind == "branch_T":
        # the value at 1e4 K; the temperature table is branch_T_table
        T, E = branch_T_table(spec, dip)
        return E[T.index(1.0e4)]
    raise ValueError(f"unknown CX_HEAT entry {spec}")


# Below its first node the table is extended to these temperatures: the
# second channel is endothermic, so its share of the rate falls there as
# its threshold Boltzmann factor.  The ratio k_b/k_a of the first node is
# carried down as r(T) = r(T1) exp[-(E_th/k)(1/T - 1/T1)], with E_th the
# channel's endothermicity (for Mg + H+ -> Mg+ 4s: 2.70 eV, where Allan et
# al. state the channel closes), and the first channel held at its node
# value.  In the table itself the 4s rate rises from 5000 to 6000 K as an
# activation energy of 2.50 eV, close to that threshold, while the 3p rate
# rises by 1.35x; holding the 3p rate overstates the 4s share slightly
# below the node.  An extrapolation, not data: the 4s share is 21.7% at the
# 5000 K node, 5.5% at 4000 K and 0.4% at 3000 K, so the defect approaches
# the 3p channel's +1.52 eV.
BRANCH_T_EXTEND = (2000.0, 3000.0, 4000.0)


def branch_T_table(spec, dip):
    """(T, E_rad(T)) of a two-channel entry branched by a rate table,
    extended below the first node by BRANCH_T_EXTEND."""
    ea = radiated_energy(spec[2][0], dip)
    eb = radiated_energy(spec[2][1], dip)
    t1, ka1, kb1 = spec[1][0]
    e_th = eb - dip                        # endothermicity of channel b [eV]
    rows = [(t, 1.0, kb1/ka1*math.exp(-e_th/(KB_EV*t) + e_th/(KB_EV*t1)))
            for t in BRANCH_T_EXTEND if t < t1] + list(spec[1])
    T = [t for t, _, _ in rows]
    E = [(ka*ea + kb*eb)/(ka + kb) for _, ka, kb in rows]
    return T, E


def charge_exchange_heat_row(el, it, eth, has_ion, has_rec):
    """CXHEAT dE_ion dE_rec g_split [eV, eV, -] of transition it, and a
    comment line naming the product states.  Every carried direction must
    have an entry."""
    entry = CX_HEAT.get((el, it), {})
    for direction, carried in (("ion", has_ion), ("rec", has_rec)):
        if carried and direction not in entry:
            raise ValueError(f"{el} transition {it}: carried {direction} "
                             "charge exchange has no CX_HEAT entry")
        if not carried and direction in entry:
            raise ValueError(f"{el} transition {it}: CX_HEAT entry for a "
                             f"{direction} direction that carries no rate")
    dip_rec = eth - IP_H          # released by X^(it+1) + H -> X^it + H+
    dE_ion = dE_rec = 0.0
    if has_ion:
        dE_ion = -dip_rec - radiated_energy(entry["ion"], -dip_rec)
    if has_rec:
        dE_rec = dip_rec - radiated_energy(entry["rec"], dip_rec)
    table = None
    if has_ion and entry.get("ion", ("",))[0] == "branch_T":
        T, E = branch_T_table(entry["ion"], -dip_rec)
        table = [(t, -dip_rec - e) for t, e in zip(T, E)]
    split = entry.get("split", 0.0)
    if split > 0.0 and not (has_ion and has_rec and dE_ion > 0.0
                            and abs(dE_ion + dE_rec) < 1e-12):
        raise ValueError(f"{el} transition {it}: split needs an exothermic "
                         "ionization and its ground-channel reverse")
    return dE_ion, dE_rec, split, entry.get("note", ""), table


def parse_badnell_rr():
    """(Z,N) -> [A,B,T0,T1,C,T2] for the M=1 ground level (C,T2 = 0 for the
    4-parameter form)."""
    out = {}
    for ln in Path(BADNELL_RR).read_text().splitlines():
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
    for ln in Path(BADNELL_DR).read_text().splitlines():
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
    txt = Path(VERNER_DIR, "phfit2.f").read_text()
    txt = re.sub(r"\n     \S", " ", txt)
    out = {}
    for m in re.finditer(r"DATA \(PH2\(I,\s*(\d+),\s*(\d+)\),I=1,7\)\s*/([^/]+)/",
                         txt):
        nz, ne = int(m.group(1)), int(m.group(2))
        out[(nz, ne)] = [float(x) for x in m.group(3).split(",")]
    return out


# subshell orbital quantum number by shell index (1s 2s 2p 3s 3p 3d 4s)
VY95_LEVEL = [0, 0, 1, 0, 1, 2, 0]


def parse_ph1_phfit2():
    """Verner & Yakovlev (1995) subshell parameters from phfit2.f itself:
    {(nz, ne, shell): [E_th, E_0, sigma_0, y_a, P, y_w]}, plus the NINN and
    NTOT tables phfit2 uses to decide which shells a ground-state ion has
    and where its PH2 outer-shell fit stops being valid."""
    txt = Path(VERNER_DIR, "phfit2.f").read_text()
    txt = re.sub(r"\n     \S", " ", txt)
    ph1 = {}
    for m in re.finditer(r"DATA \(PH1\(I,\s*(\d+),\s*(\d+),\s*(\d+)\),"
                         r"I=1,6\)\s*/([^/]+)/", txt):
        key = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
        ph1[key] = [float(x) for x in m.group(4).split(",")]
    ninn = [int(x) for x in re.search(
        r"DATA \(NINN\(I\),I=1,30\)\s*/([^/]+)/", txt).group(1).split(",")]
    ntot = [int(x) for x in re.search(
        r"DATA \(NTOT\(I\),I=1,30\)\s*/([^/]+)/", txt).group(1).split(",")]
    return ph1, ninn, ntot


def phfit2_shells(ph1, ninn, ntot, nz, ne):
    """The inner-shell switch energy and the occupied subshells of ion
    (nz, ne), exactly as phfit2 selects them: returns (einn, rows) with rows
    [l, E_th, E_0, sigma_0, y_a, P, y_w] for shells 1..nout that carry data."""
    nout = ntot[ne-1]
    if nz == ne and nz > 18:
        nout = 7
    if nz == ne + 1 and nz in (20, 21, 22, 25, 26):
        nout = 7
    if nz in (15, 17, 19) or (nz > 20 and nz != 26):
        einn = 0.0
    elif ne < 3:
        einn = 1.0e30
    else:
        einn = ph1[(nz, ne, ninn[ne-1])][0]
    rows = []
    for shell in range(1, nout+1):
        p = ph1.get((nz, ne, shell))
        if p is None or p[2] <= 0.0:
            continue          # shell not occupied in the ground state
        row = [float(VY95_LEVEL[shell-1])] + p
        if shell == nout:
            row[1] = ETH[(nz,ne)]
        rows.append(row)
    return einn, rows


def parse_cf():
    txt = Path(VERNER_DIR, "cfit.f").read_text()
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


def generate(elements=None, output_dir=OUTDIR):
    if elements is None:
        elements = ELEMENTS
    ph2 = parse_ph2()
    ph1v, ninn, ntot = parse_ph1_phfit2()
    cf = parse_cf()
    rr_bad = parse_badnell_rr()
    dr_bad = parse_badnell_dr()
    with open(os.path.join(DBASE, "VERSION")) as fh:
        chianti_version = fh.read().strip()

    missing = []
    current_stages = {el: n for el, _, n in ELEMENTS}
    for el, Z, nstage in elements:
        if el == "fe" and nstage >= 4:
            verify_final_recombination_source()
        for stage in range(1, nstage):
            key = (Z, Z - stage + 1)
            if key not in ETH:
                missing.append(f"{el} stage {stage}: ETH")
            if key not in cf:
                missing.append(f"{el} stage {stage}: CI")
            if stage >= current_stages[el]:
                if (el, stage) not in CX_HIGH or (el, stage) not in CX_HEAT:
                    missing.append(f"{el} stage {stage}: unresolved added-boundary "
                                   "charge exchange and product energies")
    if missing:
        raise ValueError("Missing selected atomic sources: " + ", ".join(missing))

    outputs = []
    for el, Z, nstage in elements:
        path = os.path.join(output_dir, f"element_{el}.txt")
        fh = io.StringIO()
        fh.write(f"# MoCHII element rate data: {el.upper()} "
                 f"(Z={Z}, {nstage} stages, highest terminal)\n")
        fh.write("# ETH: physical ionization energy [eV], NIST ASD snapshot 2022-11-15\n")
        fh.write("# outer PHSHELL edge uses ETH; inner edges retain PH1; CI retains Voronov fit energy\n")
        fh.write("# PHOTO: E_th[eV] E0 sigma0[Mb] ya P yw y0 y1 "
                 "(VFKY96 outer shell, phfit2.f PH2; valid below E_inn)\n")
        fh.write("# PHINNER: E_inn[eV] (phfit2 einn; at E >= E_inn, or "
                 "always when E_inn = 0, sigma = sum of the PHSHELL fits)\n")
        fh.write("# PHSHELL: l E_th[eV] E0 sigma0[Mb] ya P yw "
                 "(Verner & Yakovlev 1995 subshell fit, phfit2.f PH1)\n")
        fh.write("# PHTABLE: n, then n rows eps[eV] sigma[Mb] (R-matrix ground-term "
                 "cross section against the photoelectron energy eps = E - ETH, linear "
                 "in eps; replaces PHOTO below E_join = ETH + eps_n, where it equals "
                 "PHOTO; tools/fitting/ground_term_photoionization_table.py)\n")
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
                 "H+ / H0; transition 1 from Huang et al. 2023 Table 4, "
                 "nitrogen from Kingdon & Ferland 1996;\n")
        fh.write("#   forms: 0 none; 1 a(T/300)^b exp(-c/T4); "
                 "2 aT4^b + cT4^d; 3 (aT4^b + cT4^d)exp(-e/T);\n")
        fh.write("#   4 exp(sum c_k lnT^k); 5 form4*exp(-c6/T4);\n")
        fh.write("#   6 a t4^b (1 + c exp(d t4)) (Kingdon & Ferland 1996; "
                 "first-transition t4 clamped to 0.5-5, except N: 0.01-5);\n")
        fh.write("#   7 form6*exp(-e/t4) (KF96 with an activation barrier; "
                 "only the form-6 part is clamped, exp(-e/t4) takes the "
                 "actual T);\n")
        fh.write("#   8 tabulated: CXRTAB n T_1..T_n [K] k_1..k_n [cm^3 s^-1], "
                 "log-log in T, T clamped to the table (Barragan et al. "
                 "2006 for N2+ and O2+; KF96 Table 1 for other higher rates);\n")
        fh.write("#   9 resolved by product term (sulfur, Zhao et al. 2005 "
                 "Table III): CXSTATE nterm nT g_react dIP[eV]; CXSTATET "
                 "T_1..T_n; CXSTATEK label lev_lo lev_hi g E[eV] k_1..k_n "
                 "(forward rate into that term).  Forward = sum over terms; "
                 "reverse = sum over terms of pop_j(T, n_e) k_j g_react/(2 g_j) "
                 "exp(-(dIP - E_j)/kT), pop_j from the n-level atom of the "
                 "product ion)\n")
        fh.write("# CXHEAT: dE_ion dE_rec g_split -- kinetic energy [eV] given "
                 "to the gas by one X^it + H+ -> X^(it+1) + H event (dE_ion) "
                 "and one X^(it+1) + H -> X^it + H+ event (dE_rec),\n")
        fh.write("#   dE = dIP - E_rad (E_rad: product excitation that "
                 "radiates, or the photon of radiative charge transfer; "
                 "< 0 endothermic); 0 where no rate is carried;\n")
        fh.write("#   g_split > 0: the ionization rate includes radiative "
                 "charge transfer and only its non-radiative part, the "
                 "detailed-balance reverse of the recombination rate with "
                 "this statistical-weight ratio, heats\n")
        fh.write(f"# generated by tools/fitting/make_element_data.py on {DATE}\n")
        fh.write(f"ELEMENT {el} {Z} {nstage}\n")
        for i in range(1, nstage):
            ne_ion = Z - (i - 1)          # electrons of stage i
            fh.write(f"TRANSITION {i}\n")
            fh.write(f"ETH {ETH[(Z, ne_ion)]:.12e}\n")
            einn, shells = phfit2_shells(ph1v, ninn, ntot, Z, ne_ion)
            if einn > 0.0:
                # the VFKY96 outer-shell fit, valid below E_inn
                p = ph2[(Z, ne_ion)]
                fh.write("PHOTO " + " ".join(f"{v:.12e}" for v in
                         [ETH[(Z, ne_ion)]] + p) + "\n")
            # element absent from the VFKY96 outer-shell table (P, Cl,
            # K): E_inn = 0 and the PH1 shell sum applies throughout.
            fh.write(f"PHINNER {einn:.6e}\n")
            for row in shells:
                fh.write("PHSHELL " + " ".join(f"{v:.12e}" for v in row)
                         + "\n")
            if einn > 0.0:
                tab = ground_term_photoionization_table(
                    el, i, ETH[(Z, ne_ion)], [ETH[(Z, ne_ion)]] + ph2[(Z, ne_ion)])
                if tab is not None:
                    eps, sig, notes = tab
                    fh.write(f"# PHTABLE: {notes['reference']}; source "
                             f"data/atomic/{notes['source']} (md5 {notes['md5']}); "
                             f"joins PHOTO at E_join = {notes['e_join']:.6f} eV, the first "
                             f"crossing above the last target threshold "
                             f"(ETH + {notes['e_last']:.3f} eV)\n")
                    fh.write(f"PHTABLE {notes['n']}\n")
                    for e, v in zip(eps, sig):
                        fh.write(f"{e:.10e} {v:.10e}\n")
            v = cf[(Z, ne_ion)]
            fh.write("CI " + " ".join(f"{x:.6e}" for x in
                     [v[0], v[1], v[2], v[3], v[4]]) + "\n")
            # recombining ion of transition i is X^i+ (N = Z-i electrons).
            Nrec = Z - i
            #--- RR: Badnell 2023 where tabulated, else CHIANTI.
            if el == "fe" and i in FE_RR:
                fh.write(f"# RR/DR: Zhang et al. 2026, Tables 6/7, ground-level zero-density; "
                         f"domain {20*i*i:g}--{1e7*i*i:g} K; analytic extrapolation "
                         "outside the domain is not validated.\n")
                fh.write("# DR S1 compact fits do not reproduce every author ADF09 node "
                         "within the reported 2%; see the measured source audit.\n")
                fh.write("RR " + " ".join(f"{x:.6e}" for x in FE_RR[i]) + "\n")
            elif (Z, Nrec) in rr_bad:
                fh.write("RR " + " ".join(f"{x:.6e}" for x in rr_bad[(Z, Nrec)])
                         + "\n")
            else:
                rr = read_rrparams(el, i + 1)
                if rr[0] == "bad":
                    fh.write("RR " + " ".join(f"{x:.6e}" for x in rr[1]) + "\n")
                else:   # power law alpha = A (T/1e4)^-eta
                    fh.write("RR2 " + " ".join(f"{x:.6e}" for x in rr[1]) + "\n")
            #--- DR: Badnell 2023 where tabulated, else CHIANTI.
            if el == "fe" and i in FE_DR_S1:
                c,E = FE_DR_S1[i]
                fh.write(f"DR {len(c)} " + " ".join(f"{x:.6e}" for x in c) +
                         "  " + " ".join(f"{x:.6e}" for x in E) + "\n")
            elif (Z, Nrec) in dr_bad:
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
                if el in CX_LTE_UPPER:
                    fh.write(f"# CXR: detailed balance of the total forward "
                             f"rate ({CX_LTE_UPPER[el]}), i.e. the ion's "
                             "excited terms in LTE -- an upper bound on the "
                             "reverse from the ground term at nebular "
                             "densities (no forward rate into that term is "
                             "available)\n")
                for key in ("CXI", "CXR"):
                    form, coef = CX[el][key]
                    fh.write(f"{key} {form} " +
                             " ".join(f"{x:.6e}" for x in coef) + "\n")
                if CX[el]["CXI"][0] == 9:
                    for row in cx_state_rows(el, i, ETH[(Z, ne_ion)]):
                        fh.write(row + "\n")
            elif (el, i) in CX_HIGH:
                form, coef = CX_HIGH[(el, i)]
                fh.write("CXI 0 0 0 0 0 0 0\n")
                if form == 8:
                    tab = BARRAGAN[coef]
                    fh.write("CXR 8 0 0 0 0 0 0\n")
                    fh.write(f"CXRTAB {len(BARRAGAN['T'])} " +
                             " ".join(f"{x:.6e}" for x in BARRAGAN["T"]) +
                             "  " + " ".join(f"{x*1e-9:.6e}" for x in tab)
                             + "\n")
                elif form == 6:
                    # KF96 Table 1 gives a different valid T interval for
                    # each reaction. Preserve CX_HIGH's published fit and
                    # check its coefficients before tabulating it. The
                    # existing Fortran form 8 interpolates log(k) in log(T).
                    reaction = {**kf96_existing_rates(), **kf96_added_rates(),
                                **kf96_recombination_total_rates()}.get((el, i))
                    if reaction is None:
                        raise ValueError(f"{el} transition {i}: no published "
                                         "KF96 temperature interval")
                    fit = (coef[0] * 1e9, *coef[1:4])
                    if any(not math.isclose(a, b, rel_tol=1e-12,
                                            abs_tol=1e-12)
                           for a, b in zip(fit, reaction.coefficients)):
                        raise ValueError(f"{el} transition {i}: CX_HIGH "
                                         "coefficients differ from KF96")
                    temperatures, rates = kf96_charge_transfer_table(reaction)
                    table = 2 if (el,i) in (("cl",3),("cl",4),("fe",4)) else 1
                    domain = "selected inside sampled temperatures" if table == 2 else "valid"
                    source = f"KF96 Table {table}"
                    if (el, i) in (("n", 4), ("si", 2)):
                        source = ("updated KF96 database (Ferland et al. 2013, "
                                  "Cloudy ctrecombdata.dat)")
                    fh.write(f"# CXR: {source}, {domain} from "
                             f"{reaction.min_temperature:g} to "
                             f"{reaction.max_temperature:g} K; held at the "
                             "endpoint rates outside this interval\n")
                    fh.write("CXR 8 0 0 0 0 0 0\n")
                    fh.write(f"CXRTAB {len(temperatures)} " +
                             " ".join(f"{x:.6e}" for x in temperatures) +
                             "  " + " ".join(f"{x:.6e}" for x in rates) +
                             "\n")
                else:
                    if (el,i) in LANGEVIN_RECOMBINATION:
                        fh.write('# CXR: user-adopted polarization capture upper limit; '
                                 'H(1s), reduced collision mass, unit-transfer probability; '
                                 'temperature-independent approximation, not a measured rate\n')
                    fh.write(f"CXR {form} " +
                             " ".join(f"{x:.6e}" for x in coef) + "\n")
            else:
                fh.write("CXI 0 0 0 0 0 0 0\nCXR 0 0 0 0 0 0 0\n")
            #--- charge-exchange heating of the reactions just written
            if i == 1:
                has_ion = CX[el]["CXI"][0] > 0
                has_rec = CX[el]["CXR"][0] > 0
            else:
                has_ion = False
                has_rec = (el, i) in CX_HIGH
            dE_ion, dE_rec, split, note, table = charge_exchange_heat_row(
                el, i, ETH[(Z, ne_ion)], has_ion, has_rec)
            if note:
                fh.write(f"# CXHEAT: {note}\n")
            fh.write(f"CXHEAT {dE_ion:.6f} {dE_rec:.6f} {split:.6f}\n")
            if table:
                #--- dE_ion(T) of a branched ionization; linear in T,
                #--- held at the end values outside the table
                fh.write(f"CXHEATTAB {len(table)} " + " ".join(
                    f"{t:.1f} {e:.6f}" for t, e in table) + "\n")
        outputs.append((path, fh.getvalue()))
    for path, contents in outputs:
        atomic_write_text(path, contents)
        print(f"wrote {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--element", action="append", choices=[row[0] for row in ELEMENTS],
                        help="Generate only this element; repeat to select several")
    parser.add_argument("--stages", action="append", metavar="ELEMENT=N",
                        help="Override the terminal stage for a selected element")
    parser.add_argument("--output-dir", default=OUTDIR,
                        help="Directory for generated element files")
    args = parser.parse_args(argv)
    roster = {el: (z, n) for el, z, n in ELEMENTS}
    selected = set(args.element or roster)
    stages = {}
    for item in args.stages or ():
        el, separator, value = item.partition("=")
        if not separator or el not in roster or el not in selected:
            parser.error(f"invalid or unselected stage specification: {item}")
        try:
            count = int(value)
        except ValueError:
            parser.error(f"invalid stage count: {item}")
        if not 2 <= count <= roster[el][0] + 1:
            parser.error(f"stage count outside 2..{roster[el][0] + 1}: {item}")
        stages[el] = count
    elements = [(el, z, stages.get(el, n)) for el, z, n in ELEMENTS if el in selected]
    generate(elements, args.output_dir)


if __name__ == "__main__":
    main()
