module qmc_mod
!---------------------------------------------------------------------------
! MoCHII: random-access Owen-scrambled Sobol points for quasi-Monte Carlo
! photon launching (all source configurations + diffuse packets).
!
! Self-contained; no external libraries.  The launch of a stellar packet is
! a low-dimensional integral (frequency bin, direction, and, for a superset of
! sources, component and entry-surface variables); a scrambled Sobol net fills
! that cube more uniformly than independent draws, reducing
! launch noise while keeping every post-launch transport draw on the existing
! Mersenne Twister.  See docs/QUASI_RANDOM_LAUNCH.md.
!
! Sobol point by index (random access):
!   x_d(i) = XOR over the set bits k of i of V_d(k),
! with 32-bit direction numbers V_d(k).  Dimension 1 is van der Corput base 2
! (V(k) = 2^(32-k)); dimensions 2..12 use the Joe & Kuo (2008, new-joe-kuo-6)
! primitive polynomials and initial direction integers, expanded by the
! Joe-Kuo direction-number recurrence.  The recurrence and initial numbers
! reproduce scipy.stats.qmc.Sobol (which uses the same table); this is
! verified in tests/qmc.
!
! Owen scrambling by hashing (Burley 2020, "Practical Hash-based Owen
! Scrambling", JCGT): reverse the 32 bits, apply the Laine-Karras nested
! permutation keyed by a seed for each dimension, reverse back.  This preserves
! the net's equidistribution while producing independent randomized
! replicates.  The keys for each dimension are derived deterministically from
! par%qmc_seed, so a replicate is reproducible from the input seed and fixed
! across the nonlinear iterations (common random numbers).
!
! Author: Kwang-il Seon (KASI)
!---------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only : int64, real64
  implicit none
  private

  integer, parameter :: wp = real64

  !--- number of embedded Sobol dimensions and the number consumed by the
  !--- single-point launch (frequency, cos theta, phi).  QMC_MAXDIM covers the
  !--- superset stellar layout (7 dims: component, frequency, mu, phi, and up to
  !--- three surface coordinates) and the diffuse launch (9 dims: leaf, three
  !--- positions, mu, phi, channel, two energy variates); 12 is the largest
  !--- Joe & Kuo (2008) new-joe-kuo-6 dimension embedded here.
  integer, parameter, public :: QMC_MAXDIM     = 12
  integer, parameter, public :: QMC_NDIM_USED  = 3

  !--- two decorrelated scramble streams sharing the same Sobol net: the stellar
  !--- launch (stream 1) and the diffuse launch (stream 2).  Distinct keys for
  !--- each dimension make the two streams independent randomized replicates.
  integer, parameter, public :: QMC_STREAM_STELLAR = 1
  integer, parameter, public :: QMC_STREAM_DIFFUSE = 2
  integer, parameter :: NSTREAM = 2

  integer, parameter :: NBITS = 32
  integer(int64), parameter :: MASK32 = int(z'FFFFFFFF', int64)
  real(wp),       parameter :: TWO_M32 = 1.0_wp / 4294967296.0_wp   ! 2^-32

  !--- direction numbers V_d(k) (unsigned 32-bit values held in int64) and the
  !--- Owen scramble seeds, one for each dimension and stream.
  integer(int64), save :: Vdir(NBITS, QMC_MAXDIM)
  integer(int64), save :: dkey(QMC_MAXDIM, NSTREAM)
  logical,        save :: qmc_ready = .false.

  public :: qmc_setup, qmc_uniforms, qmc_uniforms_stream, qmc_uniforms_raw

contains

  !=========================================================================
  ! Precompute the direction numbers for dimensions 1..QMC_MAXDIM and the
  ! Owen scramble seeds (one for each dimension) from the input seed.  Cheap; every rank
  ! calls it once.
  !=========================================================================
  subroutine qmc_setup(seed)
    implicit none
    integer, intent(in) :: seed
    !--- Joe & Kuo (2008) new-joe-kuo-6 table for the Sobol dimensions after
    !--- van der Corput (which is dimension 1 here).  sdeg = polynomial degree
    !--- s; apoly = the interior-coefficient bit pattern a; minit = the s
    !--- initial direction integers m_1..m_s (odd, m_i < 2^i).
    integer, parameter :: sdeg(2:QMC_MAXDIM)  = [1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5]
    integer, parameter :: apoly(2:QMC_MAXDIM) = [0, 1, 1, 2, 1, 4, 2, 4, 7, 11, 13]
    integer :: minit(QMC_MAXDIM, 5)
    integer :: d, k, j, s, a
    integer(int64) :: v

    minit = 0
    minit(2, 1:1) = [1]
    minit(3, 1:2) = [1, 3]
    minit(4, 1:3) = [1, 3, 1]
    minit(5, 1:3) = [1, 1, 1]
    minit(6, 1:4) = [1, 1, 3, 3]
    minit(7, 1:4) = [1, 3, 5, 13]
    minit(8, 1:5) = [1, 1, 5, 5, 17]
    minit(9, 1:5) = [1, 1, 5, 5,  5]
    minit(10,1:5) = [1, 1, 7, 11, 19]
    minit(11,1:5) = [1, 1, 5, 1,  1]
    minit(12,1:5) = [1, 1, 1, 3, 11]

    !--- dimension 1: van der Corput base 2, V(k) = 2^(32-k).
    do k = 1, NBITS
       Vdir(k, 1) = ishft(1_int64, NBITS - k)
    end do

    !--- dimensions 2..12: Joe-Kuo direction-number recurrence on the scaled V.
    !---   V(k) = m_k << (32-k)                                  for k <= s
    !---   V(k) = V(k-s) XOR (V(k-s) >> s)
    !---          XOR [ sum_{j=1}^{s-1} a_j V(k-j) ]             for k >  s
    !--- with a_j = bit (s-1-j) of a.
    do d = 2, QMC_MAXDIM
       s = sdeg(d)
       a = apoly(d)
       do k = 1, s
          Vdir(k, d) = iand(ishft(int(minit(d, k), int64), NBITS - k), MASK32)
       end do
       do k = s + 1, NBITS
          v = ieor(Vdir(k-s, d), ishft(Vdir(k-s, d), -s))
          do j = 1, s - 1
             if (iand(ishft(int(a, int64), -(s-1-j)), 1_int64) == 1_int64) &
                v = ieor(v, Vdir(k-j, d))
          end do
          Vdir(k, d) = iand(v, MASK32)
       end do
    end do

    !--- Owen scramble seeds for each dimension and stream, from (qmc_seed, dimension,
    !--- stream).  Stream 1 (stellar) keeps the exact stage-1 expression, so the
    !--- stellar launch set is unchanged; stream 2 (diffuse) adds a distinct
    !--- 32-bit stream constant, decorrelating the two streams while both remain
    !--- balanced scrambles of the same net.
    do d = 1, QMC_MAXDIM
       dkey(d, QMC_STREAM_STELLAR) = hash32( iand( int(seed, int64) &
                 + int(d, int64) * int(z'9E3779B9', int64), MASK32 ) )
       dkey(d, QMC_STREAM_DIFFUSE) = hash32( iand( int(seed, int64) &
                 + int(d, int64) * int(z'9E3779B9', int64) &
                 + int(z'85EBCA6B', int64), MASK32 ) )
    end do

    qmc_ready = .true.
  end subroutine qmc_setup

  !=========================================================================
  ! Scrambled uniforms for the 0-based global index idx, dimensions 1..size(u).
  ! Each coordinate is strictly inside (0,1).
  !=========================================================================
  subroutine qmc_uniforms(idx, u)
    implicit none
    integer(int64), intent(in)  :: idx
    real(wp),       intent(out) :: u(:)
    call qmc_uniforms_stream(idx, u, QMC_STREAM_STELLAR)
  end subroutine qmc_uniforms

  !=========================================================================
  ! Scrambled uniforms from a chosen stream (QMC_STREAM_STELLAR /
  ! QMC_STREAM_DIFFUSE).  Same Sobol net, stream-specific scramble keys, so the
  ! two streams are decorrelated randomized replicates.
  !=========================================================================
  subroutine qmc_uniforms_stream(idx, u, stream)
    implicit none
    integer(int64), intent(in)  :: idx
    real(wp),       intent(out) :: u(:)
    integer,        intent(in)  :: stream
    integer :: d
    integer(int64) :: code
    do d = 1, size(u)
       code = owen_scramble( sobol_code(idx, d), dkey(d, stream) )
       u(d) = (real(code, wp) + 0.5_wp) * TWO_M32
    end do
  end subroutine qmc_uniforms_stream

  !=========================================================================
  ! Unscrambled uniforms (debug / cross-check against scipy).  Same mapping
  ! u = (code + 0.5) 2^-32 but with no Owen scramble.
  !=========================================================================
  subroutine qmc_uniforms_raw(idx, u)
    implicit none
    integer(int64), intent(in)  :: idx
    real(wp),       intent(out) :: u(:)
    integer :: d
    integer(int64) :: code
    do d = 1, size(u)
       code = sobol_code(idx, d)
       u(d) = (real(code, wp) + 0.5_wp) * TWO_M32
    end do
  end subroutine qmc_uniforms_raw

  !=========================================================================
  ! Unscrambled Sobol integer for index idx (0-based) and dimension d:
  !   x = XOR over set bits k of idx of V_d(k).
  ! Bits beyond 32 have no direction number and contribute nothing (so the
  ! index is used modulo 2^32; ample for any packet count in use).
  !=========================================================================
  pure function sobol_code(idx, d) result(code)
    implicit none
    integer(int64), intent(in) :: idx
    integer,        intent(in) :: d
    integer(int64) :: code, ii
    integer :: k
    code = 0_int64
    ii   = idx
    k    = 1
    do while (ii /= 0_int64 .and. k <= NBITS)
       if (iand(ii, 1_int64) == 1_int64) code = ieor(code, Vdir(k, d))
       ii = ishft(ii, -1)
       k  = k + 1
    end do
  end function sobol_code

  !=========================================================================
  ! Owen (nested-uniform) scramble of a 32-bit code by a keyed hash
  ! (Burley 2020): reverse bits, Laine-Karras permutation, reverse back.
  !=========================================================================
  pure function owen_scramble(x, seed) result(y)
    implicit none
    integer(int64), intent(in) :: x, seed
    integer(int64) :: y
    y = reverse32(x)
    y = laine_karras(y, seed)
    y = reverse32(y)
  end function owen_scramble

  pure function laine_karras(x0, seed) result(x)
    implicit none
    integer(int64), intent(in) :: x0, seed
    integer(int64) :: x
    x = iand(x0 + seed, MASK32)
    x = iand(ieor(x, iand(x * int(z'6C50B47C', int64), MASK32)), MASK32)
    x = iand(ieor(x, iand(x * int(z'B82F1E52', int64), MASK32)), MASK32)
    x = iand(ieor(x, iand(x * int(z'C7AFE638', int64), MASK32)), MASK32)
    x = iand(ieor(x, iand(x * int(z'8D22F6E6', int64), MASK32)), MASK32)
  end function laine_karras

  !--- 32-bit reversal.
  pure function reverse32(n0) result(n)
    implicit none
    integer(int64), intent(in) :: n0
    integer(int64) :: n
    n = iand(n0, MASK32)
    n = iand(ior(ishft(n, 16), ishft(n, -16)), MASK32)
    n = iand(ior(ishft(iand(n, int(z'00FF00FF', int64)),  8),  &
                 ishft(iand(n, int(z'FF00FF00', int64)), -8)), MASK32)
    n = iand(ior(ishft(iand(n, int(z'0F0F0F0F', int64)),  4),  &
                 ishft(iand(n, int(z'F0F0F0F0', int64)), -4)), MASK32)
    n = iand(ior(ishft(iand(n, int(z'33333333', int64)),  2),  &
                 ishft(iand(n, int(z'CCCCCCCC', int64)), -2)), MASK32)
    n = iand(ior(ishft(iand(n, int(z'55555555', int64)),  1),  &
                 ishft(iand(n, int(z'AAAAAAAA', int64)), -1)), MASK32)
  end function reverse32

  !--- lowbias32 integer avalanche hash (Wellons), for deriving scramble keys.
  pure function hash32(x0) result(x)
    implicit none
    integer(int64), intent(in) :: x0
    integer(int64) :: x
    x = iand(x0, MASK32)
    x = iand(ieor(x, ishft(x, -16)), MASK32)
    x = iand(x * int(z'7FEB352D', int64), MASK32)
    x = iand(ieor(x, ishft(x, -15)), MASK32)
    x = iand(x * int(z'846CA68B', int64), MASK32)
    x = iand(ieor(x, ishft(x, -16)), MASK32)
  end function hash32

end module qmc_mod
