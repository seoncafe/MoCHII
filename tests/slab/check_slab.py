#!/usr/bin/env python3
"""Automatic pass/fail checker for the plane-parallel (xy-periodic) slab mode.

Runs the small tests/slab inputs (gas_niter=0, single band bin, pure-absorption
frozen neutral column) with mpirun and asserts analytic gates read from each
run's '<base>_slab_Imu.txt' boundary output:

  G-beam-normal   collimated beam from the top at theta=0: the transmitted
                  fraction (escape at the bottom face) = exp(-tau) = 0.8874,
                  the reflected fraction (escape at the top) = 0.
  G-beam-oblique  the SAME slab, beam at theta=60 (mu = cos60 = 0.5): the
                  transmission through the slab is exp(-tau/mu) = T_normal^(1/mu),
                  so ln(T_oblique)/ln(T_normal) = 1/mu = 2 (no hardcoded tau).
  G-dda-vs-shared the DDA and 'shared' car walks give the same transmission.
  G-nx1-vs-nx2    nx=1 and nx=2 periodic tiles give the same transmission
                  (the result is invariant under the tile size).
  G-budget        for every run (reflected + transmitted + absorbed)/L_in = 1.

Assumes ../../MoCHII.x is already built.  Prints PASS/FAIL per gate and exits
nonzero on any failure.  Set SLAB_NP to change the MPI rank count (default 8).
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXE  = os.path.join(HERE, "..", "..", "MoCHII.x")
NP   = os.environ.get("SLAB_NP", "8")

# recorded analytic normal-incidence transmission exp(-tau) for these inputs
T_NORMAL_ANALYTIC = 0.8874
MU_OBLIQUE        = 0.5          # cos(60 deg) in the theta=60 inputs


def run(inp):
    """Run one slab input with mpirun; abort the checker on a nonzero exit."""
    base = os.path.splitext(inp)[0]
    log  = os.path.join(HERE, base + ".log")
    with open(log, "w") as fh:
        rc = subprocess.call(["mpirun", "-np", NP, EXE, inp],
                             cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
    if rc != 0:
        print(f"  RUN FAILED: {inp} (exit {rc}); see {base}.log")
        sys.exit(2)
    return base


def read_budget(base):
    """Parse the '<base>_slab_Imu.txt' header into a dict of fractions."""
    path = os.path.join(HERE, base + "_slab_Imu.txt")
    frac = {}
    warn = None
    with open(path) as fh:
        for ln in fh:
            m = re.search(r"#\s*(L_\w+).*fraction\s*=\s*([-\d.Ee+]+)", ln)
            if m:
                frac[m.group(1)] = float(m.group(2))
            if "not separable" in ln:
                warn = "both"
    return frac, warn


def transmitted(base):
    """Transmitted fraction = escape at the bottom face (top-lit beam)."""
    frac, _ = read_budget(base)
    return frac["L_escape_bottom"]


def budget_closes(base, tol=1e-4):
    frac, _ = read_budget(base)
    s = frac["L_escape_top"] + frac["L_escape_bottom"] + frac["L_absorbed"]
    return abs(s - 1.0) < tol, s


def main():
    passed = True
    import math

    # ---- run the fast (gas_niter=0) inputs -----------------------------------
    b_norm = run("ga_beam_n1_dda.in")
    b_o60d = run("ga60_n1_dda.in")
    b_o60s = run("ga60_n1_sh.in")
    b_o60n2 = run("ga60_n2_dda.in")
    b_iso  = run("ga_iso_n1.in")

    T0  = transmitted(b_norm)
    T60d = transmitted(b_o60d)
    T60s = transmitted(b_o60s)
    T60n2 = transmitted(b_o60n2)

    # ---- G-beam-normal -------------------------------------------------------
    frac_n, _ = read_budget(b_norm)
    ok = (abs(T0 - T_NORMAL_ANALYTIC) < 2e-3) and (frac_n["L_escape_top"] < 1e-3)
    passed &= ok
    print(f"\n=== G-beam-normal (theta=0) ===")
    print(f"  transmitted = {T0:.5f}   analytic exp(-tau) = {T_NORMAL_ANALYTIC:.4f}"
          f"   (|diff| = {abs(T0-T_NORMAL_ANALYTIC):.2e})")
    print(f"  reflected   = {frac_n['L_escape_top']:.5f}  (expect 0, no scattering)")
    print(f"  --> {'PASS' if ok else 'FAIL'} (|T-0.8874|<2e-3, reflected<1e-3)")

    # ---- G-beam-oblique: ln(T60)/ln(T0) = 1/mu = 2 ---------------------------
    ratio = math.log(T60d) / math.log(T0)
    ok = abs(ratio - 1.0 / MU_OBLIQUE) < 2e-2
    passed &= ok
    print(f"\n=== G-beam-oblique (theta=60, mu=0.5) ===")
    print(f"  T_oblique = {T60d:.5f}   T_normal^(1/mu) = {T0**(1.0/MU_OBLIQUE):.5f}")
    print(f"  ln(T60)/ln(T0) = {ratio:.4f}   expect 1/mu = {1.0/MU_OBLIQUE:.4f}")
    print(f"  --> {'PASS' if ok else 'FAIL'} (|ratio-2|<0.02)")

    # ---- G-dda-vs-shared -----------------------------------------------------
    ok = abs(T60d - T60s) < 1e-4
    passed &= ok
    print(f"\n=== G-dda-vs-shared (car_walk) ===")
    print(f"  T(dda) = {T60d:.6f}   T(shared) = {T60s:.6f}   |diff| = {abs(T60d-T60s):.2e}")
    print(f"  --> {'PASS' if ok else 'FAIL'} (|diff|<1e-4)")

    # ---- G-nx1-vs-nx2 --------------------------------------------------------
    ok = abs(T60d - T60n2) < 1e-4
    passed &= ok
    print(f"\n=== G-nx1-vs-nx2 (periodic tile size invariance) ===")
    print(f"  T(nx=1) = {T60d:.6f}   T(nx=2) = {T60n2:.6f}   |diff| = {abs(T60d-T60n2):.2e}")
    print(f"  --> {'PASS' if ok else 'FAIL'} (|diff|<1e-4)")

    # ---- G-budget (all runs) -------------------------------------------------
    print(f"\n=== G-budget (reflected + transmitted + absorbed = L_in) ===")
    allok = True
    for base in (b_norm, b_o60d, b_o60s, b_o60n2, b_iso):
        bok, s = budget_closes(base)
        allok &= bok
        print(f"  {base:16s} sum = {s:.6f}   {'ok' if bok else 'FAIL'}")
    # both-face run: budget + the not-separable note
    b_both = run("ga_both.in")
    frac_both, warn = read_budget(b_both)
    s = frac_both["L_escape_top"] + frac_both["L_escape_bottom"] + frac_both["L_absorbed"]
    note_ok = warn == "both"
    allok &= (abs(s - 1.0) < 1e-4) and note_ok
    print(f"  {b_both:16s} sum = {s:.6f}   note='{warn}' (expect 'both')   "
          f"{'ok' if (abs(s-1.0)<1e-4 and note_ok) else 'FAIL'}")
    passed &= allok
    print(f"  --> {'PASS' if allok else 'FAIL'}")

    # ---- grazing-ray max-step guard must NOT fire on these runs --------------
    print(f"\n=== G-no-graze-guard (max-step guard silent on normal runs) ===")
    fired = []
    for base in (b_norm, b_o60d, b_o60s, b_o60n2, b_iso, b_both):
        log = os.path.join(HERE, base + ".log")
        with open(log) as fh:
            if any("max-step guard" in ln for ln in fh):
                fired.append(base)
    ok = not fired
    passed &= ok
    print(f"  guard fired in: {fired if fired else 'none'}")
    print(f"  --> {'PASS' if ok else 'FAIL'} (guard must stay silent)")

    print(f"\nOVERALL: {'PASS' if passed else 'FAIL'}")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
