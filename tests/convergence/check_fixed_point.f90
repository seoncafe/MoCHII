program check_fixed_point
!---------------------------------------------------------------------------
! Gate: the par%conv_crit = 'vol' stopping test (fixed_point_distance_mod).
!  (a) Replayed on the recorded HII40 history (hii40_changes.txt, the
!      production input run for a fixed 80 iterations), the test must stop
!      at iteration 58 -- where the state equals the 80-iteration state --
!      and not at 37 or 47, where the step test without and with helium
!      stopped with the He I front 0.17 pc and ~0.04 pc short.
!  (b) A geometric series d_n = d_0 rho^n must not pass while its remaining
!      drift d rho/(1 - rho) is above tol, and must pass once it is below.
!  (c) A slow decay (rho = 0.97, every ratio below 1) must never pass on
!      its step, however small the step.
!  (d) A flat noise floor below tol (ratios above and below 1) must pass.
!---------------------------------------------------------------------------
  use define, only : wp
  use fixed_point_distance_mod
  implicit none
  real(kind=wp), parameter :: TOL = 5.0e-3_wp, TOL_TE = 1.0e-2_wp
  type(iteration_history) :: hH, hHe, h
  real(kind=wp) :: dH(200), dHe(200), dTe(200), d, e
  character(len=256) :: line
  integer :: n, i, it, stop_at, u, ios
  logical :: ok

  ok = .true.
  !--- (a)
  open(newunit=u, file='tests/convergence/hii40_changes.txt', status='old')
  n = 0
  do
     read(u,'(a)', iostat=ios) line
     if (ios /= 0) exit
     if (line(1:1) == '#') cycle
     n = n + 1
     read(line,*) it, dH(n), dHe(n), dTe(n)
  end do
  close(u)
  stop_at = 0
  do i = 1, n
     call history_push(hH, dH(i));  call history_push(hHe, dHe(i))
     if (fixed_point_converged(hH, TOL) .and. fixed_point_converged(hHe, TOL) &
         .and. dTe(i) < TOL_TE) then
        stop_at = i;  exit
     end if
  end do
  write(*,'(a,i0,a)') ' (a) recorded HII40 history: stops at iteration ', stop_at, &
     ' (expected 58)'
  if (stop_at /= 58) ok = .false.

  !--- (b) geometric, rho = 0.8: remaining = d*4, passes once d < tol/4
  call history_reset(h);  d = 1.0e-1_wp
  do i = 1, 60
     call history_push(h, d)
     e = d*0.8_wp/0.2_wp
     if (i >= 4) then
        if (fixed_point_converged(h, TOL) .neqv. (e < TOL)) then
           write(*,'(a,i0,a,es10.3)') ' (b) FAIL at step ', i, ', e = ', e
           ok = .false.
        end if
     end if
     d = d*0.8_wp
  end do
  write(*,'(a)') ' (b) geometric series: passes exactly when d rho/(1-rho) < tol'

  !--- (c) slow decay: never passes on a small step while e > tol
  call history_reset(h);  d = 4.0e-3_wp
  do i = 1, 30
     call history_push(h, d)
     if (i >= 4 .and. fixed_point_converged(h, TOL) .and. d*0.97_wp/0.03_wp > TOL) then
        write(*,'(a,i0)') ' (c) FAIL: slow decay passed at step ', i
        ok = .false.
     end if
     d = d*0.97_wp
  end do
  write(*,'(a)') ' (c) slow decay (rho = 0.97) does not pass on its step'

  !--- (d) noise floor at 2e-3, alternating up and down
  call history_reset(h)
  do i = 1, 6
     call history_push(h, 2.0e-3_wp*(1.0_wp + 0.1_wp*(-1.0_wp)**i))
  end do
  if (.not. fixed_point_converged(h, TOL)) then
     write(*,'(a)') ' (d) FAIL: noise floor below tol did not pass'
     ok = .false.
  else
     write(*,'(a)') ' (d) noise floor below tol passes'
  end if

  if (ok) then
     write(*,'(a)') ' GATE (vol stopping test): PASS'
  else
     write(*,'(a)') ' GATE (vol stopping test): FAIL'
     stop 1
  end if
end program check_fixed_point
