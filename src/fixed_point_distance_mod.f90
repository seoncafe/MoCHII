module fixed_point_distance_mod
!---------------------------------------------------------------------------
! MoCHII: distance of the gas iteration from its fixed point, estimated from
! the history of two volume-integrated measures of one fraction (e.g. x_HII)
! over successive iterations:
!
!   d_n = sum V |x_n - x_(n-1)| / sum V x_n     the step (L1), and
!   s_n = sum V (x_n - x_(n-1)) / sum V x_n     its signed counterpart, the
!                                               net change of the ionized
!                                               volume.
!
! If the iteration converges linearly, successive steps fall by a constant
! factor rho < 1, d_(n+1) = rho d_n, and what is still to come after
! iteration n sums to
!
!     e_n = d_n rho / (1 - rho).
!
! rho is the least-squares slope of ln d over the last NWIN iterations, and
! d_n in e_n is that fit's value at the latest iteration, so one noisy step
! moves neither.  When the steps do not fall (rho >= RHO_FLOOR) the linear
! model does not apply -- early on because the iteration is still in its
! transient, late because the step has reached its Monte Carlo floor -- and
! the estimate is reported as unavailable.
!
! The two cases are told apart by s.  On the floor the leaf fractions
! fluctuate about their fixed point, so the leaves' changes cancel in s: its
! mean over the window is zero within its scatter, and that scatter is
! smaller than d because d adds the same fluctuations without cancellation.
! A front that is still moving, however slowly, changes the ionized volume
! with one sign, iteration after iteration, and its mean over the window
! stands out of the scatter.
!
! STOPPING TEST (fixed_point_converged, par%conv_crit = 'vol').  With at
! least NWIN iterations of history, a measure has converged when
!   rho < RHO_FLOOR  and  e_n < tol          (still decaying: extrapolate), or
!   rho >= RHO_FLOOR and  d_n < tol  and  |mean s| < T_DRIFT sigma(s)/sqrt(NWIN)
!                                            (the Monte Carlo floor: no net
!                                             drift the window can resolve),
! with mean and sigma(s) the sample mean and standard deviation of the last
! NWIN values and T_DRIFT the two-sided 95% point of Student's t for NWIN - 1
! degrees of freedom.  The drift condition is the optional argument
! drift_test (on when absent), which the run passes as
! par%conv_drift_test (default .false. since 2026-10-18); without it the
! floor passes on d_n < tol alone.
!
! Why (2026-09-26).  The test of 2026-09-25 took rho as the median of the
! last three ratios d_k/d_(k-1) and read the floor from any of those ratios
! reaching 1.  A slow decay under Monte Carlo noise satisfies both: the m9
! restart of the HII40 benchmark (hydrogenic alpha_B, 1.1% lower than the
! state it restarted from) stopped after 4 iterations, while 16 further
! iterations moved the H front out by 0.011 pc and L(H beta) by +0.64% with
! the step decaying 1.25e-3 -> 8.05e-4 (rho ~ 0.97): a front drift that the
! old test could not see and the signed measure shows directly.  The test
! it replaced was in turn calibrated on the HII40 production input run for
! 80 iterations (results/continuous_energy/convergence/), where x_HII decays
! with rho ~ 0.85-0.87 to a ~1e-4 floor by iteration 60 and the He0-He+
! step, ~5x larger, with rho ~ 0.85-0.94 to a ~2e-3 floor near 55; the step
! test before THAT one (gas_tol = 5e-3 on d_n) stopped at 37, with the He I
! front 0.17 pc short.
!---------------------------------------------------------------------------
  use define, only : wp
  implicit none
  private

  public :: iteration_history, history_push, fixed_point_distance, &
            fixed_point_converged, history_reset, net_drift, NWIN
  public :: history_set, hist_dx, hist_dhe1, hist_dhe3

  integer, parameter :: NHIST = 512
  !--- iterations in the window of the fits and of the drift test
  integer, parameter :: NWIN = 8
  real(kind=wp), parameter :: RHO_FLOOR = 0.95_wp
  !--- Student's t, two-sided 95%, NWIN - 1 = 7 degrees of freedom
  real(kind=wp), parameter :: T_DRIFT = 2.365_wp

  type :: iteration_history
     integer :: n = 0
     real(kind=wp) :: d(NHIST) = 0.0_wp     ! step, L1
     real(kind=wp) :: s(NHIST) = 0.0_wp     ! signed step
  end type iteration_history

  !--- the run's histories of x_HII, He0-He+ and He+-He2+ (main.f90 pushes
  !--- them each iteration; the rates file carries them as ConvHist_*, and
  !--- par%restart_iteration_history restores them, previous_state_mod)
  type(iteration_history) :: hist_dx, hist_dhe1, hist_dhe3

contains

  subroutine history_push(h, d, s)
    type(iteration_history), intent(inout) :: h
    real(kind=wp),           intent(in)    :: d, s
    if (h%n < NHIST) then
       h%n = h%n + 1
    else
       h%d(1:NHIST-1) = h%d(2:NHIST)
       h%s(1:NHIST-1) = h%s(2:NHIST)
    end if
    h%d(h%n) = d
    h%s(h%n) = s
  end subroutine history_push

  !--- a history as written by a previous run (d(1:n), s(1:n) oldest first)
  subroutine history_set(h, d, s)
    type(iteration_history), intent(inout) :: h
    real(kind=wp),           intent(in)    :: d(:), s(:)
    integer :: n
    n = min(size(d), NHIST)
    h%n = n
    h%d(1:n) = d(size(d)-n+1:)
    h%s(1:n) = s(size(s)-n+1:)
  end subroutine history_set

  subroutine history_reset(h)
    type(iteration_history), intent(inout) :: h
    h%n = 0
  end subroutine history_reset

  !--- the stopping test of the header, for one measure
  logical function fixed_point_converged(h, tol, drift_test) result(done)
    type(iteration_history), intent(in) :: h
    real(kind=wp),           intent(in) :: tol
    logical, optional,       intent(in) :: drift_test
    real(kind=wp) :: rho, dist, smean, tstat
    logical :: available, test_drift
    test_drift = .true.
    if (present(drift_test)) test_drift = drift_test
    done = .false.
    if (h%n < NWIN) return
    call fixed_point_distance(h, rho, dist, available)
    if (available) then
       done = dist < tol
    else
       done = h%d(h%n) < tol
       if (test_drift) then
          call net_drift(h, smean, tstat)
          done = done .and. tstat < T_DRIFT
       end if
    end if
  end function fixed_point_converged

  !--- rho and e_n of the latest entry, from the fit of ln d over the last
  !--- min(n, NWIN) entries.  available = .false. with fewer than four
  !--- entries, a non-positive step in the window, or rho >= RHO_FLOOR
  !--- (rho is still returned in the last case).
  subroutine fixed_point_distance(h, rho, dist, available)
    type(iteration_history), intent(in)  :: h
    real(kind=wp),           intent(out) :: rho, dist
    logical,                 intent(out) :: available
    real(kind=wp) :: xm, ym, sxx, sxy, slope, dfit
    integer :: m, k, i0
    rho = 0.0_wp;  dist = 0.0_wp;  available = .false.
    m = min(h%n, NWIN)
    if (m < 4) return
    i0 = h%n - m
    if (any(h%d(i0+1:h%n) <= 0.0_wp)) return
    !--- least squares, ln d_k = a + b k, k = 1..m
    xm = 0.5_wp*real(m + 1, wp)
    ym = sum(log(h%d(i0+1:h%n)))/real(m, wp)
    sxx = 0.0_wp;  sxy = 0.0_wp
    do k = 1, m
       sxx = sxx + (real(k, wp) - xm)**2
       sxy = sxy + (real(k, wp) - xm)*(log(h%d(i0+k)) - ym)
    end do
    slope = sxy/sxx
    rho = exp(slope)
    if (rho >= RHO_FLOOR) return
    dfit = exp(ym + slope*(real(m, wp) - xm))
    dist = dfit*rho/(1.0_wp - rho)
    available = .true.
  end subroutine fixed_point_distance

  !--- mean of the signed step over the last min(n, NWIN) entries and its
  !--- Student statistic |mean| / (sigma/sqrt(m)); tstat = 0 with fewer than
  !--- two entries or a zero mean, huge() for a nonzero mean with no scatter.
  subroutine net_drift(h, smean, tstat)
    type(iteration_history), intent(in)  :: h
    real(kind=wp),           intent(out) :: smean, tstat
    real(kind=wp) :: var
    integer :: m, i0
    smean = 0.0_wp;  tstat = 0.0_wp
    m = min(h%n, NWIN)
    if (m < 2) return
    i0 = h%n - m
    smean = sum(h%s(i0+1:h%n))/real(m, wp)
    if (smean == 0.0_wp) return
    var = sum((h%s(i0+1:h%n) - smean)**2)/real(m - 1, wp)
    if (var > 0.0_wp) then
       tstat = abs(smean)/sqrt(var/real(m, wp))
    else
       tstat = huge(1.0_wp)
    end if
  end subroutine net_drift

end module fixed_point_distance_mod
