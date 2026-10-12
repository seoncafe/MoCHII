program check_adjacent_ion_fractions
  use define, only : wp, MAX_ION_STAGES
  use ion_stage_equilibrium_mod, only : adjacent_ion_fractions
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  real(kind=wp) :: up(MAX_ION_STAGES-1), down(MAX_ION_STAGES-1)
  real(kind=wp) :: f(MAX_ION_STAGES), expected(MAX_ION_STAGES), residual, scale
  integer :: status, i, ncase

  ncase = 0
  up = 2.0_wp; down = 1.0_wp
  do i = 1, MAX_ION_STAGES
     expected(i) = 2.0_wp**(i-1)
  end do
  expected = expected/sum(expected)
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 0 .and. maxval(abs(f-expected)) < 1.0e-14_wp, 'analytic positive-rate chain')
  call balance()

  ! Finite rates whose ratio and cumulative products overflow in ordinary
  ! arithmetic must still yield finite normalized stationary populations.
  up = 1.0e200_wp; down = 1.0e-200_wp
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 0 .and. f(MAX_ION_STAGES) == 1.0_wp, 'extreme ionization')
  call normalized()
  up = 1.0e-200_wp; down = 1.0e200_wp
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 0 .and. f(1) == 1.0_wp, 'extreme recombination')
  call normalized()

  ! Intermediate logarithms span an enormous range, but two populated
  ! stages have equal weight after the rates return to the initial ratio.
  up = 1.0e200_wp; down = 1.0e-200_wp
  up(2::2) = 1.0e-200_wp; down(2::2) = 1.0e200_wp
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 0 .and. maxval(abs(f(2::2)-1.0_wp/3.0_wp)) < 1.0e-14_wp, 'alternating extreme rates')
  call normalized()

  up = 1.0_wp; down = 0.0_wp
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 0 .and. f(MAX_ION_STAGES) == 1.0_wp, 'one-way ionization')
  up = 0.0_wp; down = 1.0_wp
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 0 .and. f(1) == 1.0_wp, 'one-way recombination')

  ! A unique connected interval with loss from both surrounding intervals.
  up = 2.0_wp; down = 1.0_wp
  down(2) = 0.0_wp; up(4) = 0.0_wp
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 0 .and. f(3) > 0.0_wp .and. f(4) > 0.0_wp, 'unique closed interval')
  call require(sum(f(1:2))+sum(f(5:)) == 0.0_wp, 'transient stages empty')
  call balance()

  up = 1.0_wp; down = 1.0_wp
  up(3) = 0.0_wp; down(3) = 0.0_wp
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 2, 'disconnected nonunique equilibrium')
  up = 0.0_wp; down = 0.0_wp
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 2, 'all zero rates are nonunique')
  up(1) = -1.0_wp
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 1, 'negative rate rejected')
  up(1) = ieee_value(1.0_wp, ieee_quiet_nan)
  call adjacent_ion_fractions(up, down, f, status)
  call require(status == 1, 'nonfinite rate rejected')
  call adjacent_ion_fractions([real(kind=wp)::], [real(kind=wp)::], f(1:1), status)
  call require(status == 0 .and. f(1) == 1.0_wp, 'single stage')
  write(*,'(a,i0,a)') 'PASS: ', ncase, ' adjacent-stage checks'
contains
  subroutine require(condition, name)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name
    ncase = ncase + 1
    if (.not. condition) then
       write(*,'(2a)') 'FAIL: ', name
       stop 1
    end if
  end subroutine require
  subroutine normalized()
    call require(all(ieee_is_finite(f)) .and. all(f >= 0.0_wp), 'finite nonnegative fractions')
    call require(abs(sum(f)-1.0_wp) < 1.0e-14_wp, 'unit population')
  end subroutine normalized
  subroutine balance()
    call normalized()
    residual = maxval(abs(up*f(:MAX_ION_STAGES-1)-down*f(2:)))
    scale = max(maxval(up*f(:MAX_ION_STAGES-1)), maxval(down*f(2:)), tiny(1.0_wp))
    call require(residual/scale < 1.0e-13_wp, 'stationary reaction balance')
  end subroutine balance
end program check_adjacent_ion_fractions
