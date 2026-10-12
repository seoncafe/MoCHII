program check_nonadjacent_ion_equilibrium
  use define, only: wp
  use ion_stage_stationary_mod, only: stationary_ion_fractions
  use ion_stage_equilibrium_mod, only: adjacent_ion_fractions
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_get_underflow_mode
  implicit none
  real(kind=wp) :: rates(5,5), fractions(5), adjacent(5)
  real(kind=wp) :: up(4), down(4), expected(3), balance, scale
  integer :: status, i
  logical :: gradual_underflow

  ! An Auger transition from stage one to stage three changes the stationary
  ! fractions from the adjacent-only result. The analytic ratio is 3:13:14.
  rates = 0.0_wp
  rates(2,1) = 2.0_wp
  rates(3,1) = 11.0_wp
  rates(1,2) = 3.0_wp
  rates(3,2) = 5.0_wp
  rates(2,3) = 7.0_wp
  call stationary_ion_fractions(rates(1:3,1:3), fractions(1:3), status)
  expected = [3.0_wp, 13.0_wp, 14.0_wp]/30.0_wp
  call require(status == 0, 'three-stage Auger status')
  call require(maxval(abs(fractions(1:3)-expected)) < 1.0e-14_wp, 'analytic Auger fractions')
  call require(abs(sum(fractions(1:3))-1.0_wp) < 1.0e-14_wp, 'Auger normalization')
  call require(abs(fractions(2)+2.0_wp*fractions(3)-41.0_wp/30.0_wp) &
       < 1.0e-14_wp, 'Auger mean charge')
  do i = 1, 3
     balance = sum(rates(i,1:3)*fractions(1:3)) &
          - sum(rates(1:3,i))*fractions(i)
     call require(abs(balance) < 1.0e-13_wp, 'Auger population balance')
  end do

  ! The general matrix must reduce to the established adjacent solution.
  up = [2.0_wp, 3.0_wp, 5.0_wp, 7.0_wp]
  down = [11.0_wp, 13.0_wp, 17.0_wp, 19.0_wp]
  call check_adjacent(up, down, 1.0e-13_wp, 'realistic adjacent rates')
  up = [1.0e50_wp, 1.0e-50_wp, 1.0e50_wp, 1.0e-50_wp]
  down = [1.0e-50_wp, 1.0e50_wp, 1.0e-50_wp, 1.0e50_wp]
  call check_adjacent(up, down, 1.0e-12_wp, 'large adjacent rate range')
  up = [1.0e100_wp, 1.0e-100_wp, 1.0e100_wp, 1.0e-100_wp]
  down = [1.0e-100_wp, 1.0e100_wp, 1.0e-100_wp, 1.0e100_wp]
  call check_adjacent(up, down, 1.0e-12_wp, 'extreme adjacent rate range')

  ! A representable small fraction must not be replaced by a population floor.
  rates = 0.0_wp
  rates(2,1) = 1.0e-200_wp
  rates(1,2) = 1.0_wp
  call stationary_ion_fractions(rates(1:2,1:2), fractions(1:2), status)
  call require(status == 0, 'small fraction status')
  call require(fractions(2) > 0.0_wp, 'small fraction retained')
  call require(abs(fractions(2)/1.0e-200_wp-1.0_wp) < 1.0e-12_wp, 'small fraction value')
  call ieee_get_underflow_mode(gradual_underflow)
  if (gradual_underflow) then
     rates(2,1) = 1.0e-310_wp
     call stationary_ion_fractions(rates(1:2,1:2), fractions(1:2), status)
     call require(status == 0, 'first subnormal fraction status')
     call require(fractions(2) > 0.0_wp, 'first subnormal fraction retained')
     call require(abs(fractions(2)/rates(2,1)-1.0_wp) < 1.0e-12_wp, 'first subnormal fraction value')
     rates(2,1) = 1.0e-320_wp
     call stationary_ion_fractions(rates(1:2,1:2), fractions(1:2), status)
     call require(status == 0, 'second subnormal fraction status')
     call require(fractions(2) > 0.0_wp, 'second subnormal fraction retained')
     call require(abs(fractions(2)/rates(2,1)-1.0_wp) < 1.0e-12_wp, 'second subnormal fraction value')
  else
     print *, 'SKIP subnormal fractions: gradual underflow disabled'
  end if

  ! A unique closed class can exceed the resolvable range of row-scaled
  ! elimination. Report numerical failure without manufacturing a population.
  rates = 0.0_wp
  rates(2,1) = 1.0e-300_wp
  rates(3,1) = 1.0e-300_wp
  rates(1,2) = 1.0e-300_wp
  rates(3,2) = 1.0e-300_wp
  rates(1,3) = 1.0e300_wp
  rates(2,3) = 1.0e300_wp
  call stationary_ion_fractions(rates(1:3,1:3), fractions(1:3), status)
  call require(status == 3, 'unresolved closed class status')
  call require(all(fractions(1:3) == 0.0_wp), 'unresolved closed class has no result')

  ! Transient stages feed one closed communicating class.
  rates = 0.0_wp
  rates(2,1) = 2.0_wp
  rates(3,2) = 3.0_wp
  rates(2,3) = 6.0_wp
  rates(3,4) = 5.0_wp
  call stationary_ion_fractions(rates(1:4,1:4), fractions(1:4), status)
  call require(status == 0, 'transient stages status')
  call require(all(fractions([1,4]) == 0.0_wp), 'transient stages vanish')
  call require(abs(fractions(2)-2.0_wp/3.0_wp) < 1.0e-14_wp, 'closed class stage two')
  call require(abs(fractions(3)-1.0_wp/3.0_wp) < 1.0e-14_wp, 'closed class stage three')

  rates = 0.0_wp
  rates(2,1) = 1.0_wp
  rates(3,2) = 1.0_wp
  call stationary_ion_fractions(rates(1:3,1:3), fractions(1:3), status)
  call require(status == 0, 'absorbing stage status')
  call require(all(fractions(1:3) == [0.0_wp,0.0_wp,1.0_wp]), 'absorbing stage')

  rates = 0.0_wp
  call stationary_ion_fractions(rates(1:1,1:1), fractions(1:1), status)
  call require(status == 0 .and. fractions(1) == 1.0_wp, 'single closed stage')
  call stationary_ion_fractions(rates(1:3,1:3), fractions(1:3), status)
  call require(status == 2, 'multiple isolated closed stages')
  rates(2,1) = 1.0_wp
  rates(1,2) = 1.0_wp
  call stationary_ion_fractions(rates(1:3,1:3), fractions(1:3), status)
  call require(status == 2, 'two closed classes')

  rates = 0.0_wp
  rates(2,1) = -1.0_wp
  call stationary_ion_fractions(rates(1:2,1:2), fractions(1:2), status)
  call require(status == 1, 'negative rate rejected')
  rates(2,1) = ieee_value(1.0_wp, ieee_quiet_nan)
  call stationary_ion_fractions(rates(1:2,1:2), fractions(1:2), status)
  call require(status == 1, 'nonfinite rate rejected')
  rates = 0.0_wp
  rates(1,1) = 1.0_wp
  call stationary_ion_fractions(rates(1:2,1:2), fractions(1:2), status)
  call require(status == 1, 'nonzero diagonal rejected')
  call stationary_ion_fractions(rates(1:2,1:2), fractions(1:3), status)
  call require(status == 1, 'mismatched dimensions rejected')
  rates = 0.0_wp
  rates(2,1) = huge(1.0_wp)
  rates(3,1) = huge(1.0_wp)
  call stationary_ion_fractions(rates(1:3,1:3), fractions(1:3), status)
  call require(status == 1, 'outgoing rate overflow rejected')
  print *, 'PASS nonadjacent stationary ion fractions'

contains
  subroutine check_adjacent(ionization, recombination, tolerance, label)
    real(kind=wp), intent(in) :: ionization(:), recombination(:), tolerance
    character(len=*), intent(in) :: label
    integer :: stage, result
    rates = 0.0_wp
    do stage = 1, 4
       rates(stage+1,stage) = ionization(stage)
       rates(stage,stage+1) = recombination(stage)
    end do
    call stationary_ion_fractions(rates, fractions, result)
    call require(result == 0, label//' stationary status')
    call adjacent_ion_fractions(ionization, recombination, adjacent, result)
    call require(result == 0, label//' adjacent status')
    call require(maxval(abs(fractions-adjacent)) < tolerance, label//' fractions')
    call require(abs(sum(fractions)-1.0_wp) < tolerance, label//' normalization')
    do stage = 1, 5
       balance = sum(rates(stage,:)*fractions) - sum(rates(:,stage))*fractions(stage)
       scale = max(sum(rates(stage,:)*fractions), sum(rates(:,stage))*fractions(stage))
       call require(abs(balance) < tolerance*max(scale,tiny(1.0_wp)), label//' balance')
    end do
  end subroutine check_adjacent

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (condition) return
    print *, 'FAIL ', label
    stop 1
  end subroutine require
end program check_nonadjacent_ion_equilibrium
