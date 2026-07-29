program check_hei_branches
!---------------------------------------------------------------------------
! Gate for hei_cascade_mod: the branch fractions of the case-B He I cascade.
!
! Checks:
!   1. the three fractions are positive and sum to one, so the channel emits
!      exactly the recombinations alpha_B(He II) feeds it;
!   2. they reproduce the ratios of BSS99 Table 1 + JSM as tabulated in
!      K. Wood's ionize/probset.f, recomputed here independently;
!   3. the ordering and magnitude are physical -- 2^3S dominates, and the
!      fractions stay within a few percent of the low-density statistical
!      weights 3/4 : 1/12 : 1/6 rather than reordering them.  An earlier
!      version of the module misread nebcont_mod's singlet TOTAL as
!      alpha_eff(2^1S) and produced 0.746 : 0.224 : 0.031; this bound is what
!      catches that class of error;
!   4. the temperature trend has the sign the fits imply: the 2^1P exponent
!      is the steepest, so its fraction falls with T while 2^3S rises;
!   5. the fits are clamped, not extrapolated, outside 3e3-2e4 K.
!
! Deliberately NOT checked: that the three coefficients sum to alpha_B(He II).
! They do not -- 0.90 of it at 8000 K -- and that shortfall is a statement
! about two independent atomic data sets, reported by hei_cascade_setup and
! printed here, not a tolerance to enforce.
!---------------------------------------------------------------------------
  use define
  use hei_cascade_mod, only : hei_branch_fractions
  use recomb_mod,      only : alphaB_HeII
  use mpi
  implicit none

  !--- independent copy of the source coefficients (Wood, probset.f)
  real(kind=wp), parameter :: TCHK(6) = &
     [3.0e3_wp, 5.0e3_wp, 8.0e3_wp, 1.0e4_wp, 1.5e4_wp, 2.0e4_wp]
  real(kind=wp) :: T, f3s, f1s, f1p, r3, r1s, r1p, s, dev
  real(kind=wp) :: p3s(6), p1p(6)
  real(kind=wp) :: g3s, g1s, g1p
  integer :: ierr, i
  logical :: ok

  call MPI_INIT(ierr)
  ok = .true.

  write(*,'(a)') ' T [K]      2^3S      2^1S      2^1P       sum    ' // &
     'sum(alpha)/alphaB'
  do i = 1, 6
     T = TCHK(i)
     call hei_branch_fractions(T, f3s, f1s, f1p)
     r3  = 2.10e-13_wp*(T/1.0e4_wp)**(-0.381_wp)
     r1s = 2.06e-14_wp*(T/1.0e4_wp)**(-0.451_wp)
     r1p = 4.17e-14_wp*(T/1.0e4_wp)**(-0.695_wp)
     s   = r3 + r1s + r1p
     write(*,'(f7.0,4f10.5,f14.5)') T, f3s, f1s, f1p, f3s+f1s+f1p, &
        s/alphaB_HeII(T)
     p3s(i) = f3s;  p1p(i) = f1p

     !--- 1. closure of the fractions themselves
     if (min(f3s, f1s, f1p) <= 0.0_wp) then
        write(*,'(a)') '  FAIL: a branch fraction is not positive'
        ok = .false.
     end if
     if (abs(f3s + f1s + f1p - 1.0_wp) > 1.0e-12_wp) then
        write(*,'(a)') '  FAIL: fractions do not sum to one'
        ok = .false.
     end if

     !--- 2. against the source coefficients, recomputed here
     dev = max(abs(f3s - r3/s), abs(f1s - r1s/s), abs(f1p - r1p/s))
     if (dev > 1.0e-12_wp) then
        write(*,'(a,es10.3)') '  FAIL: fractions differ from BSS99 ratios by ', dev
        ok = .false.
     end if

     !--- 3. physical ordering and distance from the statistical weights
     if (.not. (f3s > f1p .and. f1p > f1s)) then
        write(*,'(a)') '  FAIL: branch ordering is not 2^3S > 2^1P > 2^1S'
        ok = .false.
     end if
     if (abs(f3s - 0.75_wp) > 0.05_wp .or. abs(f1s - 1.0_wp/12.0_wp) > 0.02_wp &
         .or. abs(f1p - 1.0_wp/6.0_wp) > 0.05_wp) then
        write(*,'(a)') '  FAIL: fractions depart from the statistical ' // &
           'weights by more than a correction'
        ok = .false.
     end if
  end do

  !--- 4. temperature trend
  write(*,'(/,a)') ' trend over 3e3-2e4 K:'
  write(*,'(a,f9.5,a,f9.5)') '   2^3S rises ', p3s(1), ' -> ', p3s(6)
  write(*,'(a,f9.5,a,f9.5)') '   2^1P falls ', p1p(1), ' -> ', p1p(6)
  do i = 2, 6
     if (p3s(i) <= p3s(i-1) .or. p1p(i) >= p1p(i-1)) then
        write(*,'(a)') '  FAIL: temperature trend has the wrong sign'
        ok = .false.
     end if
  end do

  !--- 5. clamping outside the fit range
  call hei_branch_fractions(1.0e2_wp, f3s, f1s, f1p)
  call hei_branch_fractions(3.0e3_wp, g3s, g1s, g1p)
  if (abs(f3s-g3s) + abs(f1s-g1s) + abs(f1p-g1p) > 1.0e-12_wp) then
     write(*,'(a)') '  FAIL: fits are extrapolated below 3e3 K instead of clamped'
     ok = .false.
  end if
  call hei_branch_fractions(1.0e6_wp, f3s, f1s, f1p)
  call hei_branch_fractions(2.0e4_wp, g3s, g1s, g1p)
  if (abs(f3s-g3s) + abs(f1s-g1s) + abs(f1p-g1p) > 1.0e-12_wp) then
     write(*,'(a)') '  FAIL: fits are extrapolated above 2e4 K instead of clamped'
     ok = .false.
  end if
  write(*,'(a)') ' clamped outside 3e3-2e4 K'

  write(*,'(/,a)') 'GATE (He I cascade branch fractions): ' // &
     merge('PASS', 'FAIL', ok)
  call MPI_FINALIZE(ierr)
  if (.not. ok) error stop 1

end program check_hei_branches
