module ion_stage_stationary_mod
! Stationary populations for directed ionization and recombination transitions.
! rates(destination, source) has units s^-1 and a zero diagonal. A unique
! closed communicating class supports the stationary population; states that
! feed that class are transient and have zero stationary population.
! Finite-precision elimination may reject a strongly connected rate network
! whose stationary vector is too ill-conditioned to resolve. Fractions below
! the smallest representable real number can underflow to zero.
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
    integer :: active(size(fractions)), n, i, j, k, first, nclosed, m, pivot
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
       factor = sum(rates(:,j)/norm)
       if (norm > 1.0_wp) then
          if (factor > huge(1.0_wp)/norm) return
       end if
       outgoing(j) = norm*factor
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
       matrix(i,1:m) = matrix(i,1:m)/norm
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
          if (abs(matrix(i,j)) > 1.0_wp) then
             if (abs(tmp(j)) > huge(1.0_wp)/abs(matrix(i,j))) then
                status = 3
                return
             end if
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
    factor = sum(tmp(1:m)/norm)
    fractions(active(1:m)) = (tmp(1:m)/norm)/factor

    ! Check each original balance equation after local rate scaling. A failed
    ! residual is not repaired by clipping a negative or tiny population.
    do i = 1, m
       j = active(i)
       norm = max(outgoing(j), maxval(rates(j,:)))
       flow = sum((rates(j,:)/norm)*fractions)
       balance = flow - (outgoing(j)/norm)*fractions(j)
       if (abs(balance) > 1.0e-10_wp*max(flow, (outgoing(j)/norm)*fractions(j), &
                                        tiny(1.0_wp))) then
          fractions = 0.0_wp
          status = 3
          return
       end if
    end do
    status = 0
  end subroutine stationary_ion_fractions
end module ion_stage_stationary_mod
