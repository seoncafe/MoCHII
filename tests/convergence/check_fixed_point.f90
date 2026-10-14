program check_fixed_point
!---------------------------------------------------------------------------
! Gate: the par%conv_crit = 'vol' stopping test (fixed_point_distance_mod),
! on synthetic histories of the step d and the signed step s.  Noise is a
! fixed deterministic sequence, so the gate is reproducible.
!  (b) A geometric series d_n = d_0 rho^n with s = d (a front moving one
!      way) must not pass while its remaining drift d rho/(1 - rho) is
!      above tol, and must pass once it is below, from NWIN entries on.
!  (c) A slow decay (rho = 0.97) with 20% noise on d and a signed step of
!      half of d must never pass while its remaining drift is above tol --
!      the m9 HII40 restart, which the test of 2026-09-25 stopped after 4
!      iterations.
!  (d) A noise floor below tol, with s fluctuating about zero, must pass.
!  (e) The same floor with a net drift standing out of the scatter
!      (mean s = 4 sigma(s)/sqrt(NWIN)) must not.
!  (f) Fewer than NWIN entries never pass.
!---------------------------------------------------------------------------
  use define, only : wp
  use fixed_point_distance_mod
  implicit none
  real(kind=wp), parameter :: TOL = 5.0e-3_wp
  !--- zero-mean unit-variance noise, 16 values repeated
  real(kind=wp), parameter :: NZ(16) = [ 0.53_wp, -1.21_wp,  0.88_wp, -0.37_wp, &
                                          1.42_wp, -0.64_wp, -1.05_wp,  0.21_wp, &
                                         -0.12_wp,  1.17_wp, -1.46_wp,  0.71_wp, &
                                          0.33_wp, -0.93_wp,  1.26_wp, -0.73_wp]
  type(iteration_history) :: h
  real(kind=wp) :: d, e, z, drift
  integer :: i
  logical :: ok

  ok = .true.

  !--- (b) geometric, rho = 0.8: remaining = d*4, passes once d < tol/4
  call history_reset(h);  d = 1.0e-1_wp
  do i = 1, 60
     call history_push(h, d, d)
     e = d*0.8_wp/0.2_wp
     if (fixed_point_converged(h, TOL) .neqv. (i >= NWIN .and. e < TOL)) then
        write(*,'(a,i0,a,es10.3)') ' (b) FAIL at step ', i, ', e = ', e
        ok = .false.
     end if
     d = d*0.8_wp
  end do
  write(*,'(a)') ' (b) geometric series: passes exactly when d rho/(1-rho) < tol'

  !--- (c) slow decay under noise: never passes while e > tol
  call history_reset(h);  d = 1.2e-3_wp
  do i = 1, 40
     z = NZ(mod(i-1, 16) + 1)
     call history_push(h, d*(1.0_wp + 0.2_wp*z), 0.5_wp*d*(1.0_wp + 0.2_wp*NZ(mod(i+6, 16) + 1)))
     if (fixed_point_converged(h, TOL) .and. d*0.97_wp/0.03_wp > TOL) then
        write(*,'(a,i0)') ' (c) FAIL: noisy slow decay passed at step ', i
        ok = .false.
     end if
     d = d*0.97_wp
  end do
  write(*,'(a)') ' (c) noisy slow decay (rho = 0.97, s = d/2) does not pass'

  !--- (d) floor at 2e-3, s fluctuating about zero with sigma 2e-4
  call history_reset(h)
  do i = 1, NWIN
     call history_push(h, 2.0e-3_wp*(1.0_wp + 0.1_wp*NZ(i)), 2.0e-4_wp*NZ(i + 4))
  end do
  if (.not. fixed_point_converged(h, TOL)) then
     write(*,'(a)') ' (d) FAIL: noise floor below tol did not pass'
     ok = .false.
  else
     write(*,'(a)') ' (d) noise floor below tol, no net drift: passes'
  end if

  !--- (e) the same floor with a net drift of 4 sigma/sqrt(NWIN)
  call history_reset(h)
  drift = 4.0_wp*2.0e-4_wp/sqrt(real(NWIN, wp))
  do i = 1, NWIN
     call history_push(h, 2.0e-3_wp*(1.0_wp + 0.1_wp*NZ(i)), drift + 2.0e-4_wp*NZ(i + 4))
  end do
  if (fixed_point_converged(h, TOL)) then
     write(*,'(a)') ' (e) FAIL: floor with a resolved net drift passed'
     ok = .false.
  else
     write(*,'(a)') ' (e) floor with a resolved net drift: does not pass'
  end if
  !--- (e') the same history with the drift test off (par%conv_drift_test
  !--- = .false.) passes on d < tol alone; on explicitly passes as (e) does not
  if (.not. fixed_point_converged(h, TOL, drift_test=.false.)) then
     write(*,'(a)') " (e') FAIL: drift test off, floor below tol did not pass"
     ok = .false.
  else if (fixed_point_converged(h, TOL, drift_test=.true.)) then
     write(*,'(a)') " (e') FAIL: drift_test = .true. differs from the default"
     ok = .false.
  else
     write(*,'(a)') " (e') drift test off: passes on d < tol; explicit on: as (e)"
  end if

  !--- (f) a floor seen for NWIN - 1 iterations only
  call history_reset(h)
  do i = 1, NWIN - 1
     call history_push(h, 2.0e-3_wp*(1.0_wp + 0.1_wp*NZ(i)), 2.0e-4_wp*NZ(i + 4))
     if (fixed_point_converged(h, TOL)) then
        write(*,'(a,i0)') ' (f) FAIL: passed with ', i, ' entries'
        ok = .false.
     end if
  end do
  write(*,'(a,i0,a)') ' (f) fewer than ', NWIN, ' entries: does not pass'

  if (ok) then
     write(*,'(a)') ' GATE (vol stopping test): PASS'
  else
     write(*,'(a)') ' GATE (vol stopping test): FAIL'
     stop 1
  end if
end program check_fixed_point
