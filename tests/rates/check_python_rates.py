#!/usr/bin/env python3
"""Gate: tools/python/mochii_rates.py must agree with the Fortran it mirrors.

The analytic gates and the paper's figure generators build their 1D references
from the Python module.  If it drifts from the Fortran, every deviation those
references report acquires a silent offset -- which is exactly the defect this
module was written to remove.  So the agreement is checked, not assumed.

Three groups are compared, each against the production objects through a small
driver:

  recombination   recomb_mod's alpha_A, alpha_B, alpha_1 for H II, He II and
                  He III at par%recomb_model = 'badnell_mao'.  Both sides
                  evaluate the same closed form, so they must agree to
                  roundoff;
  free-free       cooling_mod's gbar_ff(T, Z) for Z = 1..4.  The Python side
                  carries the nodes dumped from the Fortran table and repeats
                  its linear-in-log10 T interpolation, so this too is a
                  roundoff comparison -- it catches a stale table, a changed
                  node grid, or a changed u-quadrature in gaunt_ff_mean;
  cross sections  photo_xsec's sigma_HI, sigma_HeI, sigma_HeII.  Same closed
                  forms again.  The grid deliberately starts at each threshold,
                  where sigma_HI was 0.775% off while the gates used the VFKY96
                  fit.  It then steps away from threshold by more than the
                  10 ppm neighborhood in which the Python np.isclose branch
                  stands in for the Fortran's exact E == E_th test.

The temperature comparison covers 3e3-3e4 K, the production thermal range, plus
the wider 1e3-1e5 K the recombination fits are valid over; the gbar table is
checked over its full 1e2-1e6 K span.
"""
import os, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools", "python"))
import mochii_rates as R

#--- both sides evaluate the same closed form or the same table
TOL = 1.0e-10
NAMES = ["alphaA_HII", "alphaB_HII", "alpha1_HII",
         "alphaA_HeII", "alphaB_HeII", "alpha1_HeII",
         "alphaA_HeIII", "alphaB_HeIII", "alpha1_HeIII"]
GNAMES = ["gbar_ff(Z=1)", "gbar_ff(Z=2)", "gbar_ff(Z=3)", "gbar_ff(Z=4)"]
SNAMES = ["sigma_HI", "sigma_HeI", "sigma_HeII"]
CNAMES = ["betaA_HII", "betaB_HII", "beta_HeII", "beta_HeIII",
          "ci_HI", "ci_HeI", "ci_HeII"]

OBJS = ("define.o", "utility.o", "photo_xsec.o", "recomb_mod.o", "gaunt.o",
        "gaunt_vh14_mod.o", "nlevel_mod.o", "nlevel_cooling_mod.o",
        "cooling_mod.o")

DRIVER = r"""
program rates_dump
  use define
  use recomb_mod
  use photo_xsec, only : sigma_HI, sigma_HeI, sigma_HeII
  use cooling_mod, only : cooling_setup, gbar_ff, &
                         betaA_HII, betaB_HII, beta_HeII, beta_HeIII
  use mpi
  implicit none
  integer, parameter :: NT = 61, NG = 81, NE = 60
  integer, parameter :: NCH = 3
  real(kind=wp), parameter :: ETH(NCH) = &
     [13.598_wp, 24.587_wp, 54.416_wp]
  real(kind=wp) :: T, E
  integer :: i, ich, iz, ierr
  call MPI_INIT(ierr)
  par%recomb_model = 'badnell_mao'
  par%atomic_dir = '@ATOMIC@'
  par%data_dir   = '@DATA@'
  call cooling_setup()

  !--- recombination over 1e3-1e5 K
  do i = 1, NT
     T = 10.0_wp**(3.0_wp + 2.0_wp*real(i-1,wp)/real(NT-1,wp))
     write(11,'(10es24.16)') T, alphaA_HII(T), alphaB_HII(T), &
        alphaA_HII(T)-alphaB_HII(T), alphaA_HeII(T), alphaB_HeII(T), &
        alphaA_HeII(T)-alphaB_HeII(T), alphaA_HeIII(T), alphaB_HeIII(T), &
        alphaA_HeIII(T)-alphaB_HeIII(T)
  end do

  !--- gbar_ff over the full table span, on and off the nodes
  do i = 1, NG
     T = 10.0_wp**(2.0_wp + 4.0_wp*real(i-1,wp)/real(NG-1,wp))
     write(12,'(5es24.16)') T, (gbar_ff(T, iz), iz = 1, 4)
  end do

  !--- cross sections: from each threshold up to 1 keV
  do ich = 1, NCH
     do i = 1, NE
        if (i == 1) then
           E = ETH(ich)
        else
           E = ETH(ich)*(1000.0_wp/ETH(ich))**(real(i-1,wp)/real(NE-1,wp))
        end if
        write(13,'(i2,4es24.16)') ich, E, sigma_HI(E), sigma_HeI(E), sigma_HeII(E)
     end do
  end do

  !--- recombination cooling and collisional ionization over 3e3-5e4 K,
  !--- the span the thermal gates solve over
  do i = 1, NT
     T = 10.0_wp**(3.0_wp + 2.0_wp*real(i-1,wp)/real(NT-1,wp))
     !--- es26.16e3: ci_HI falls below 1e-100 at the cold end and a
     !--- two-digit exponent field drops the E
     write(14,'(8es26.16e3)') T, betaA_HII(T), betaB_HII(T), &
        beta_HeII(T, .false.), beta_HeIII(T, .false.), &
        ci_HI(T), ci_HeI(T), ci_HeII(T)
  end do

  call MPI_FINALIZE(ierr)
end program rates_dump
"""


def fortran_tables():
    """Compile and run the driver; return the three dumped tables."""
    d = tempfile.mkdtemp()
    src = os.path.join(d, "rates_dump.f90")
    exe = os.path.join(d, "rates_dump.x")
    text = (DRIVER.replace("@ATOMIC@", os.path.join(ROOT, "data", "atomic"))
                  .replace("@DATA@", os.path.join(ROOT, "data")))
    with open(src, "w") as fh:
        fh.write(text)
    objs = [os.path.join(ROOT, "src", o) for o in OBJS]
    cmd = [os.environ.get("FC", "mpiifort"), "-cpp", "-DMPI", "-ipo", "-O2",
           "-module", os.path.join(ROOT, "src"), "-I" + os.path.join(ROOT, "src"),
           "-o", exe, src] + objs
    subprocess.run(cmd, check=True, capture_output=True)
    subprocess.run([exe], check=True, capture_output=True, text=True, cwd=d)
    out = []
    for unit in (11, 12, 13, 14):
        with open(os.path.join(d, "fort.%d" % unit)) as fh:
            out.append(np.array([[float(x) for x in ln.split()]
                                 for ln in fh if ln.strip()]))
    return out


def report(label, names, fort, py):
    """Print one block of max relative deviations; return True if all pass."""
    ok = True
    print(f"  {label}")
    for j, nm in enumerate(names):
        f = fort[:, j]
        keep = f != 0.0
        dev = np.max(np.abs(py[keep, j] / f[keep] - 1.0)) if keep.any() else 0.0
        zero_ok = np.all(py[~keep, j] == 0.0)
        flag = "" if (dev <= TOL and zero_ok) else "   FAIL"
        if flag:
            ok = False
        print(f"    {nm:16s}  {dev:.3e}{flag}")
    return ok


def main():
    rec, gff, sig, cool = fortran_tables()

    T = rec[:, 0]
    py_rec = np.column_stack([R.alphaA_HII(T), R.alphaB_HII(T), R.alpha1_HII(T),
                              R.alphaA_HeII(T), R.alphaB_HeII(T), R.alpha1_HeII(T),
                              R.alphaA_HeIII(T), R.alphaB_HeIII(T), R.alpha1_HeIII(T)])

    Tg = gff[:, 0]
    py_gff = np.column_stack([R.gbar_ff(Tg, z) for z in (1, 2, 3, 4)])

    E = sig[:, 1]
    py_sig = np.column_stack([R.sigma_HI(E), R.sigma_HeI(E), R.sigma_HeII(E)])

    Tc = cool[:, 0]
    py_cool = np.column_stack([R.betaA_HII(Tc), R.betaB_HII(Tc),
                               R.beta_HeII(Tc), R.beta_HeIII(Tc),
                               R.ci_HI(Tc), R.ci_HeI(Tc), R.ci_HeII(Tc)])

    print("  max |python/fortran - 1|")
    ok = report("recombination, 1e3-1e5 K", NAMES, rec[:, 1:], py_rec)
    ok &= report("free-free Gaunt, 1e2-1e6 K", GNAMES, gff[:, 1:], py_gff)
    ok &= report("cross sections, threshold-1 keV", SNAMES, sig[:, 2:], py_sig)
    ok &= report("cooling and collisional ionization, 1e3-1e5 K",
                 CNAMES, cool[:, 1:], py_cool)
    print()
    print("GATE (Python rates match the Fortran): " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
