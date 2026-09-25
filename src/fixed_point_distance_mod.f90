module fixed_point_distance_mod
!---------------------------------------------------------------------------
! MoCHII: distance of the gas iteration from its fixed point, estimated from
! the history of one volume-integrated change measure (e.g. the volume-mean
! |delta x_HII| of successive iterations).
!
! If the iteration converges linearly, successive changes fall by a constant
! factor rho < 1, d_{n+1} = rho d_n, and what is still to come after
! iteration n sums to
!
!     e_n = d_n rho / (1 - rho).
!
! A test on the step d_n alone stops too early whenever rho is close to 1:
! the HII40 production input, run for 60 iterations, shows rho ~ 0.86-0.90,
! so the remaining drift is 6-9 times the last step, and the step test
! (gas_tol = 5e-3) stopped with the He I front 0.17 pc short of where it
! settles.
!
! rho is the median of the last three ratios d_k/d_(k-1), which rejects a
! single noisy step.  When the changes do not fall (rho >= RHO_FLOOR) the
! linear model does not apply -- early on because the iteration is still in
! its transient, late because the measure has reached its Monte Carlo
! floor -- and the estimate is reported as unavailable.
!
! STOPPING TEST (fixed_point_converged, par%conv_crit = 'vol', 2026-09-25).
! A measure has converged when
!   rho < RHO_FLOOR  and  e_n < tol          (still decaying: extrapolate), or
!   rho >= RHO_FLOOR and  d_n < tol  and one of the last three ratios >= 1
!                                            (the Monte Carlo floor: changes
!                                             that no longer fall monotonically).
! The second branch separates a floor, where successive steps fluctuate up
! and down, from a slow decay, where every ratio is below 1 but close to it
! and the drift still to come is large; a slow decay never passes on its
! step.  Calibrated on the HII40 production input run for 80 iterations
! (results/continuous_energy/convergence/): x_HII decays with rho ~ 0.85-0.87
! and reaches a ~1e-4 floor by iteration 60; the He0-He+ measure is ~5x
! larger at the same iteration, decays with rho ~ 0.85-0.94 and reaches a
! ~2e-3 floor near iteration 55.  With tol = 5e-3 the test stops at
! iteration 58, where the state equals the 80-iteration state (He I front
! within 0.0006 pc, <T[NpNe]> within 0.8 K); the step test it replaces
! stopped at 37, 0.17 pc short, and the step test with helium added would
! have stopped at 47.
!---------------------------------------------------------------------------
  use define, only : wp
  implicit none
  private

  public :: iteration_history, history_push, fixed_point_distance, &
            fixed_point_converged, history_reset

  integer, parameter :: NHIST = 512
  real(kind=wp), parameter :: RHO_FLOOR = 0.95_wp

  type :: iteration_history
     integer :: n = 0
     real(kind=wp) :: d(NHIST) = 0.0_wp
  end type iteration_history

contains

  subroutine history_push(h, d)
    type(iteration_history), intent(inout) :: h
    real(kind=wp),           intent(in)    :: d
    if (h%n < NHIST) then
       h%n = h%n + 1
    else
       h%d(1:NHIST-1) = h%d(2:NHIST)
    end if
    h%d(h%n) = d
  end subroutine history_push

  subroutine history_reset(h)
    type(iteration_history), intent(inout) :: h
    h%n = 0
  end subroutine history_reset

  !--- the stopping test of the header, for one measure
  logical function fixed_point_converged(h, tol) result(done)
    type(iteration_history), intent(in) :: h
    real(kind=wp),           intent(in) :: tol
    real(kind=wp) :: rho, dist, rmax
    logical :: available
    integer :: k
    done = .false.
    call fixed_point_distance(h, rho, dist, available)
    if (available) then
       done = dist < tol
    else if (h%n >= 4 .and. rho >= RHO_FLOOR) then
       rmax = 0.0_wp
       do k = 1, 3
          rmax = max(rmax, h%d(h%n-3+k)/h%d(h%n-4+k))
       end do
       done = rmax >= 1.0_wp .and. h%d(h%n) < tol
    end if
  end function fixed_point_converged

  !--- rho and e_n of the latest entry.  available = .false. with fewer than
  !--- four entries, a non-positive change, or rho >= RHO_FLOOR.
  subroutine fixed_point_distance(h, rho, dist, available)
    type(iteration_history), intent(in)  :: h
    real(kind=wp),           intent(out) :: rho, dist
    logical,                 intent(out) :: available
    real(kind=wp) :: r(3), t
    integer :: k
    rho = 0.0_wp;  dist = 0.0_wp;  available = .false.
    if (h%n < 4) return
    if (any(h%d(h%n-3:h%n) <= 0.0_wp)) return
    do k = 1, 3
       r(k) = h%d(h%n-3+k)/h%d(h%n-4+k)
    end do
    !--- median of three
    if (r(1) > r(2)) then; t = r(1); r(1) = r(2); r(2) = t; end if
    if (r(2) > r(3)) then; t = r(2); r(2) = r(3); r(3) = t; end if
    if (r(1) > r(2)) then; t = r(1); r(1) = r(2); r(2) = t; end if
    rho = r(2)
    if (rho >= RHO_FLOOR) return
    dist = h%d(h%n)*rho/(1.0_wp - rho)
    available = .true.
  end subroutine fixed_point_distance

end module fixed_point_distance_mod
