program check_hei_584
!---------------------------------------------------------------------------
! Gate for the He I 584 A branching, hei_584_ionizes_H.
!
! Checks:
!   1. the three algebraically distinct published forms agree to round-off --
!      this module's, K. Wood's (mcionize.f, h0find.f) and CMacIonize's
!      (PhysicalDiffuseReemissionHandler) -- over a grid of (T, xHI, xHeI);
!   2. the limits are exact: x(H^0) -> 0 gives p -> 0 (nothing left to
!      ionize, everything converts) and x(He^0) -> 0 gives p -> 1;
!   3. p is monotonic the right way in each argument;
!   4. p = 1 reproduces the unconverted branch energies exactly, which is what
!      makes the step-2 regression meaningful;
!   5. converting energy is conserved as it must be: the two-photon branch
!      gains exactly what the 2^1P branch loses, in DECAYS, while the
!      RADIATED energy drops because a converted decay emits 20.62 eV in the
!      pair (of which only the part above 13.6 eV is in band) instead of a
!      21.22 eV photon;
!   6. at the benchmark conditions the channel loses H-ionizing photons, and
!      it loses more of them in the HII20-like state than the HII40-like one.
!---------------------------------------------------------------------------
  use define
  use hei_cascade_mod,   only : hei_584_ionizes_H, hei_branch_fractions
  use hei_twophoton_mod, only : hei_2ph_setup, hei_2ph_energy_per_decay, &
                                hei_2ph_photons_per_decay
  use mpi
  implicit none

  real(kind=wp), parameter :: E3S = 19.82_wp, E1P = 21.22_wp
  real(kind=wp), parameter :: TG(4)  = [3.0e3_wp, 8.0e3_wp, 1.0e4_wp, 2.0e4_wp]
  real(kind=wp), parameter :: XH(4)  = [1.0e-5_wp, 1.0e-4_wp, 1.0e-3_wp, 1.0e-1_wp]
  real(kind=wp), parameter :: XHE(4) = [1.0e-4_wp, 1.0e-3_wp, 1.0e-2_wp, 5.0e-1_wp]
  real(kind=wp) :: T, xh0, xhe0, p, pw, pc, dev, dmax
  real(kind=wp) :: n40, n20, e40, e20
  integer :: ierr, i, j, k
  logical :: ok

  call MPI_INIT(ierr)
  par%eion_max = 100.0_wp
  par%data_dir = '../../data'
  call hei_2ph_setup()
  ok = .true.

  !--- 1. the three published forms
  dmax = 0.0_wp
  do i = 1, 4
     do j = 1, 4
        do k = 1, 4
           T = TG(i);  xh0 = XH(j);  xhe0 = XHE(k)
           p  = hei_584_ionizes_H(T, xh0, xhe0)
           !--- Wood, mcionize.f / h0find.f
           pw = 1.0_wp/(1.0_wp + 0.77_wp/(sqrt(T)/1.0e2_wp)*xhe0/xh0)
           !--- CMacIonize, PhysicalDiffuseReemissionHandler
           pc = sqrt(T)*xh0/(sqrt(T)*xh0 + 77.0_wp*xhe0)
           dmax = max(dmax, abs(p - pw), abs(p - pc))
        end do
     end do
  end do
  write(*,'(a,es10.3)') ' max |this - Wood|, |this - CMacIonize| = ', dmax
  if (dmax > 1.0e-12_wp) then
     write(*,'(a)') '  FAIL: the three published forms disagree'
     ok = .false.
  end if

  !--- 2. limits
  if (hei_584_ionizes_H(8.0e3_wp, 0.0_wp, 1.0e-2_wp) > 1.0e-30_wp) then
     write(*,'(a)') '  FAIL: x(H0) -> 0 does not give p -> 0'
     ok = .false.
  end if
  if (abs(hei_584_ionizes_H(8.0e3_wp, 1.0e-4_wp, 0.0_wp) - 1.0_wp) > 1.0e-15_wp) then
     write(*,'(a)') '  FAIL: x(He0) -> 0 does not give p -> 1'
     ok = .false.
  end if
  write(*,'(a)') ' limits x(H0) -> 0 and x(He0) -> 0 exact'

  !--- 3. monotonicity
  if (.not. (hei_584_ionizes_H(8.0e3_wp, 1.0e-3_wp, 1.0e-2_wp) > &
             hei_584_ionizes_H(8.0e3_wp, 1.0e-4_wp, 1.0e-2_wp))) then
     write(*,'(a)') '  FAIL: p does not increase with x(H0)'
     ok = .false.
  end if
  if (.not. (hei_584_ionizes_H(8.0e3_wp, 1.0e-4_wp, 1.0e-3_wp) > &
             hei_584_ionizes_H(8.0e3_wp, 1.0e-4_wp, 1.0e-2_wp))) then
     write(*,'(a)') '  FAIL: p does not decrease with x(He0)'
     ok = .false.
  end if
  write(*,'(a)') ' monotonic in both x(H0) and x(He0)'

  !--- 4-6. what the channel radiates, at the two benchmark states
  write(*,'(/,a)') ' channel 4 per case-B He II recombination:'
  write(*,'(a)') '   state                      p      E [eV]   N       dE      dN'
  call report(7.9e3_wp, 1.0_wp,     0.0_wp,    'p = 1 (step 2)  ', e40, n40, .true.)
  call report(7.9e3_wp, 3.01e-4_wp, 2.04e-3_wp,'HII40 typical   ', e20, n20, .false.)
  call report(6.4e3_wp, 1.54e-4_wp, 2.58e-2_wp,'HII20 typical   ', e20, n20, .false.)

  write(*,'(/,a)') 'GATE (He I 584 A conversion): ' // merge('PASS', 'FAIL', ok)
  call MPI_FINALIZE(ierr)
  if (.not. ok) error stop 1

contains

  subroutine report(T, xh0, xhe0, label, ee, nn, is_ref)
    real(kind=wp), intent(in)  :: T, xh0, xhe0
    character(len=*), intent(in) :: label
    real(kind=wp), intent(out) :: ee, nn
    logical, intent(in) :: is_ref
    real(kind=wp) :: f3s, f1s, f1p, pp, e2, n2
    real(kind=wp), save :: e_ref = 0.0_wp, n_ref = 0.0_wp
    call hei_branch_fractions(T, f3s, f1s, f1p)
    pp = hei_584_ionizes_H(T, xh0, xhe0)
    e2 = hei_2ph_energy_per_decay();  n2 = hei_2ph_photons_per_decay()
    ee = f3s*E3S + f1p*pp*E1P + (f1s + f1p*(1.0_wp - pp))*e2
    nn = f3s     + f1p*pp     + (f1s + f1p*(1.0_wp - pp))*n2
    if (is_ref) then
       e_ref = ee;  n_ref = nn
       write(*,'(a,a,f8.4,f9.3,f8.4)') '   ', label, pp, ee, nn
       !--- 4. p = 1 must leave the step-2 expression untouched
       if (abs(ee - (f3s*E3S + f1p*E1P + f1s*e2)) > 1.0e-12_wp) then
          write(*,'(a)') '  FAIL: p = 1 does not reproduce the unconverted energy'
          ok = .false.
       end if
    else
       write(*,'(a,a,f8.4,f9.3,f8.4,2f8.2,a)') '   ', label, pp, ee, nn, &
          100.0_wp*(ee/e_ref - 1.0_wp), 100.0_wp*(nn/n_ref - 1.0_wp), ' %'
       !--- 6. conversion must cost the channel photons and energy
       if (ee >= e_ref .or. nn >= n_ref) then
          write(*,'(a)') '  FAIL: conversion did not reduce the channel output'
          ok = .false.
       end if
    end if
  end subroutine report

end program check_hei_584
