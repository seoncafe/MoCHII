module ion_stage_stationary_mod
! Stationary populations for directed ionization and recombination transitions.
! rates(destination, source) has units s^-1 and a zero diagonal. A unique
! closed communicating class supports the stationary population; states that
! feed that class are transient and have zero stationary population.
! Finite-precision elimination may reject a strongly connected rate network
! whose stationary vector is too ill-conditioned to resolve. Fractions below
! the smallest normal real number (tiny, 2.2e-308) are zero when the module
! is compiled with flush-to-zero (the production -fp-model fast=2), and below
! the smallest subnormal otherwise.
!
! Every normalization by a rate scale first multiplies both operands by a
! power of two (scale/exponent) and never divides by the scale itself: under
! -no-prec-div x/norm is evaluated as x*(1/norm), and 1/norm is subnormal,
! hence flushed to zero, for norm > 2**1022. The overflow tests on the total
! rates and on the back-substitution products compare binary exponents for
! the same reason. tests/ion_stages/check_nonadjacent_ion_equilibrium passes
! with the production flags and with -fp-model precise -no-ftz; its subnormal
! cases run only where gradual underflow is enabled.
  use define, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: stationary_ion_fractions

contains

  subroutine stationary_ion_fractions(rates, fractions, status)
    real(kind=wp), intent(in) :: rates(:,:)
    real(kind=wp), intent(out) :: fractions(:)
    integer, intent(out) :: status
    logical :: reach(size(fractions),size(fractions))
    logical :: member(size(fractions)), closed
    integer :: active(size(fractions)), n, i, j, k, first, nclosed, m, pivot, e
    real(kind=wp) :: outgoing(size(fractions)), norm, factor, value, term
    real(kind=wp) :: matrix(size(fractions),size(fractions)), rhs(size(fractions))
    real(kind=wp) :: tmp(size(fractions)), balance, flow

    fractions = 0.0_wp
    status = 1                         ! invalid rate matrix or dimensions
    n = size(fractions)
    if (n < 1) return
    if (size(rates,1) /= n .or. size(rates,2) /= n) return
    if (.not. all(ieee_is_finite(rates))) return
    if (any(rates < 0.0_wp)) return
    do i = 1, n
       if (rates(i,i) /= 0.0_wp) return
    end do
    outgoing = 0.0_wp
    do j = 1, n
       norm = maxval(rates(:,j))
       if (norm == 0.0_wp) cycle
       ! factor = sum/2**e lies in [0.5, n); the total rate factor*2**e is
       ! representable when exponent(factor) + e <= maxexponent.
       e = exponent(norm)
       factor = sum(scale(rates(:,j), -e))
       if (exponent(factor) > maxexponent(1.0_wp) - e) return
       outgoing(j) = scale(factor, e)
    end do

    ! Determine strongly connected classes from the original positive rates.
    ! Using exact zero here preserves a real but very weak transition.
    reach = .false.
    do i = 1, n
       reach(i,i) = .true.
       do j = 1, n
          if (rates(j,i) > 0.0_wp) reach(i,j) = .true.
       end do
    end do
    do k = 1, n
       do i = 1, n
          if (reach(i,k)) reach(i,:) = reach(i,:) .or. reach(k,:)
       end do
    end do
    nclosed = 0
    active = 0
    m = 0
    do first = 1, n
       member = reach(first,:) .and. reach(:,first)
       if (any(member(:first-1))) cycle
       closed = .true.
       do i = 1, n
          if (.not. member(i)) cycle
          do j = 1, n
             if (.not. member(j) .and. rates(j,i) > 0.0_wp) closed = .false.
          end do
       end do
       if (.not. closed) cycle
       nclosed = nclosed + 1
       if (nclosed > 1) exit
       m = count(member)
       j = 0
       do i = 1, n
          if (.not. member(i)) cycle
          j = j + 1
          active(j) = i
       end do
    end do
    status = 2                         ! nonunique stationary populations
    if (nclosed /= 1) return
    if (m == 1) then
       fractions(active(1)) = 1.0_wp
       status = 0
       return
    end if

    ! Each balance equation is divided by its own largest coefficient before
    ! elimination. This avoids a global rate scale hiding physically small
    ! rows. Partial pivoting still has finite-precision limitations for very
    ! ill-conditioned nonadjacent networks; report those as numerical failure.
    matrix = 0.0_wp
    rhs = 0.0_wp
    do i = 1, m
       do j = 1, m
          matrix(i,j) = rates(active(i),active(j))
       end do
       matrix(i,i) = -outgoing(active(i))
       norm = maxval(abs(matrix(i,1:m)))
       if (norm <= 0.0_wp) then
          status = 3
          return
       end if
       ! Exact equilibration (largest |coefficient| = 1 in every row): partial
       ! pivoting compares rows, so a power-of-two scale alone (largest
       ! coefficient anywhere in [0.5, 1)) changes the pivot sequence and
       ! loses the large-range adjacent networks of the test. Both operands
       ! are brought to [0.5, 1) first so the divisor is never above 2**1022.
       e = exponent(norm)
       matrix(i,1:m) = scale(matrix(i,1:m), -e)/scale(norm, -e)
    end do
    matrix(m,1:m) = 1.0_wp
    rhs(m) = 1.0_wp
    do k = 1, m-1
       pivot = k-1 + maxloc(abs(matrix(k:m,k)), dim=1)
       if (matrix(pivot,k) == 0.0_wp) then
          status = 3
          return
       end if
       if (pivot /= k) then
          tmp(1:m) = matrix(k,1:m)
          matrix(k,1:m) = matrix(pivot,1:m)
          matrix(pivot,1:m) = tmp(1:m)
          value = rhs(k); rhs(k) = rhs(pivot); rhs(pivot) = value
       end if
       do i = k+1, m
          factor = matrix(i,k)/matrix(k,k)
          matrix(i,k) = 0.0_wp
          matrix(i,k+1:m) = matrix(i,k+1:m) - factor*matrix(k,k+1:m)
          rhs(i) = rhs(i) - factor*rhs(k)
       end do
    end do
    if (.not. all(ieee_is_finite(matrix(1:m,1:m)))) then
       status = 3
       return
    end if
    if (.not. all(ieee_is_finite(rhs(1:m)))) then
       status = 3
       return
    end if
    if (matrix(m,m) == 0.0_wp) then
       status = 3
       return
    end if
    tmp = 0.0_wp
    do i = m, 1, -1
       value = rhs(i)
       do j = i+1, m
          ! |matrix(i,j)*tmp(j)| < 2**(exponent sum); reject before it can
          ! exceed huge.
          if (exponent(matrix(i,j)) + exponent(tmp(j)) > maxexponent(1.0_wp)) then
             status = 3
             return
          end if
          term = matrix(i,j)*tmp(j)
          if (term > 0.0_wp) then
             if (value < term-huge(1.0_wp)) then
                status = 3
                return
             end if
          else if (term < 0.0_wp) then
             if (value > huge(1.0_wp)+term) then
                status = 3
                return
             end if
          end if
          value = value-term
       end do
       if (matrix(i,i) == 0.0_wp) then
          status = 3
          return
       end if
       if (abs(matrix(i,i)) < 1.0_wp) then
          if (abs(value) > huge(1.0_wp)*abs(matrix(i,i))) then
             status = 3
             return
          end if
       end if
       tmp(i) = value/matrix(i,i)
    end do
    if (.not. all(ieee_is_finite(tmp(1:m)))) then
       status = 3
       return
    end if
    if (any(tmp(1:m) < 0.0_wp)) then
       status = 3
       return
    end if
    norm = maxval(tmp(1:m))
    if (norm <= 0.0_wp) then
       status = 3
       return
    end if
    tmp(1:m) = scale(tmp(1:m), -exponent(norm))
    factor = sum(tmp(1:m))
    fractions(active(1:m)) = tmp(1:m)/factor

    ! Check each original balance equation after local rate scaling. A failed
    ! residual is not repaired by clipping a negative or tiny population.
    do i = 1, m
       j = active(i)
       e = exponent(max(outgoing(j), maxval(rates(j,:))))
       flow = sum(scale(rates(j,:), -e)*fractions)
       value = scale(outgoing(j), -e)*fractions(j)
       balance = flow - value
       if (abs(balance) > 1.0e-10_wp*max(flow, value, tiny(1.0_wp))) then
          fractions = 0.0_wp
          status = 3
          return
       end if
    end do
    status = 0
  end subroutine stationary_ion_fractions
end module ion_stage_stationary_mod
