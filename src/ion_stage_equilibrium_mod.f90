module ion_stage_equilibrium_mod
! Adjacent ionization stages in stationary balance. Rates include all
! ionization and recombination processes in each direction, in s^-1.
! Nonadjacent transitions such as Auger ionization require a different
! balance equation and must not be represented by an adjacent rate here.
  use define, only : wp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: adjacent_ion_fractions
contains
  subroutine adjacent_ion_fractions(ionization, recombination, fractions, status)
    real(kind=wp), intent(in) :: ionization(:), recombination(:)
    real(kind=wp), intent(out) :: fractions(:)
    integer, intent(out) :: status
    real(kind=wp) :: log_population(size(fractions)), normalization
    integer :: ns, first, last, lower, upper, nclosed, i
    logical :: closed

    fractions = 0.0_wp
    status = 1                       ! invalid input
    ns = size(fractions)
    if (ns < 1) return
    if (size(ionization) /= ns-1 .or. size(recombination) /= ns-1) return
    if (.not. all(ieee_is_finite(ionization))) return
    if (.not. all(ieee_is_finite(recombination))) return
    if (any(ionization < 0.0_wp) .or. any(recombination < 0.0_wp)) return

    ! Bidirectional positive-rate links connect intervals of stages.
    ! Stationary populations can occupy only intervals with no loss to
    ! another interval. More than one such interval leaves their relative
    ! populations undetermined by stationary balance: report this rather
    ! than replacing exact zeros with an arbitrary rate floor.
    nclosed = 0
    lower = 1
    upper = ns
    first = 1
    do while (first <= ns)
       last = first
       do while (last < ns)
          if (ionization(last) == 0.0_wp .or. recombination(last) == 0.0_wp) exit
          last = last + 1
       end do
       closed = .true.
       if (first > 1) closed = recombination(first-1) == 0.0_wp
       if (last < ns) closed = closed .and. ionization(last) == 0.0_wp
       if (closed) then
          nclosed = nclosed + 1
          lower = first
          upper = last
       end if
       first = last + 1
    end do
    status = 2                       ! nonunique stationary populations
    if (nclosed /= 1) return

    ! Separate logarithms avoid overflow in an individual rate ratio and
    ! in a cumulative population product. Subtract the largest logarithm
    ! before exponentiation; negligible populations can underflow to zero.
    log_population(lower) = 0.0_wp
    do i = lower, upper-1
       log_population(i+1) = log_population(i) &
            + log(ionization(i)) - log(recombination(i))
    end do
    fractions(lower:upper) = exp(log_population(lower:upper) &
                                    - maxval(log_population(lower:upper)))
    normalization = sum(fractions(lower:upper))
    fractions(lower:upper) = fractions(lower:upper)/normalization
    status = 0
  end subroutine adjacent_ion_fractions
end module ion_stage_equilibrium_mod
