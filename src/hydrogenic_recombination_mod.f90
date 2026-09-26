module hydrogenic_recombination_mod
!---------------------------------------------------------------------------
! MoCHII: recombination of the two hydrogenic ions, H II -> H I (z = 1) and
! He III -> He II (z = 2), from ONE nl-resolved radiative cascade in the
! zero-density limit on Storey & Hummer's (1991, CPC 66, 129) hydrogenic
! A-values and cross sections (the data are generated into
! hydrogenic_recombination_data.inc by
! tools/fitting/make_hydrogenic_recombination.py):
!
!   hydrogenic_alpha_B   captures into n >= 2 [cm^3 s^-1]; recomb_mod's case-B
!                        rate for these two ions (alpha_A adds the Milne
!                        alpha_1 on the transport cross section).
!   hydrogenic_f_2s      share of those captures that ends in 2s.
!   heii_lyman_shares    He II: share that leaves by each Lyman line np -> 1s.
!
! in the case of the recombination lines, lyman_line_case(): A when the
! Lyman-line photons escape (par%lyman_lines = 'escaping'), B when they are
! absorbed on the spot and re-emitted until the cascade leaves by another
! route (np -> 1s deleted).  In case B the only Lyman exit is Ly-alpha, and
! f_2s + f_Lya = 1; in case A every np -> 1s is an exit, and the shares
! still sum to 1 (the generator checks it).  All shares are relative to
! alpha_B in both cases: a capture into 1s never feeds 2s or a Lyman line.
!
! Validity: zero density.  Collisional l-mixing, which Storey & Hummer (1995)
! include, raises the He II case-B alpha_eff(2^2S) at n_e = 100 cm^-3 above
! this limit by 0.92% at 5000 K and 0.11% at 3e4 K; the totals are not
! affected by it (alpha_B within 0.1% of SH95's at 1e4 K for both ions).
! Tabulated for log T = 1.5-7.0 in steps of 0.05; outside, the end values
! are held.  alpha_B is interpolated by a Catmull-Rom cubic in log-log, so
! its logarithmic derivative, which the recombination cooling takes by
! differencing, is continuous; the shares linearly in log T.
!---------------------------------------------------------------------------
  use define, only : wp, par
  implicit none
  private

  public :: lyman_line_case
  public :: hydrogenic_alpha_B, hydrogenic_f_2s, hydrogenic_alpha_2s
  public :: heii_lyman_shares, heii_lyman_energy
  public :: HEII_NLY

  include 'hydrogenic_recombination_data.inc'

  !--- He II Ly-alpha [eV]; a line n -> 1 lies at E_LYA (1 - 1/n^2)/(3/4)
  real(kind=wp), parameter :: E_HEII_LYA = 40.813_wp

contains

  !=========================================================================
  ! The Baker & Menzel case of the Lyman lines: A when they escape
  ! (par%lyman_lines = 'escaping'), B when they are absorbed on the spot.
  !=========================================================================
  function lyman_line_case() result(c)
    character(len=1) :: c
    if (trim(par%lyman_lines) == 'escaping') then
       c = 'A'
    else
       c = 'B'
    end if
  end function lyman_line_case

  !--- grid position of T: node i and fraction w in [0,1] (clamped)
  pure subroutine grid_position(T, i, w)
    real(kind=wp), intent(in)  :: T
    integer,       intent(out) :: i
    real(kind=wp), intent(out) :: w
    real(kind=wp) :: x, dlt
    dlt = (hyd_logT(HYD_NT) - hyd_logT(1))/real(HYD_NT - 1, wp)
    x = (min(max(log10(T), hyd_logT(1)), hyd_logT(HYD_NT)) - hyd_logT(1))/dlt
    i = min(int(x) + 1, HYD_NT - 1)
    w = x - real(i - 1, wp)
  end subroutine grid_position

  !=========================================================================
  ! alpha_B [cm^3 s^-1] of nuclear charge z (1 or 2) at T [K].  Catmull-Rom
  ! cubic in log alpha vs log T, C1 across every node.
  !=========================================================================
  elemental real(kind=wp) function hydrogenic_alpha_B(z, T) result(a)
    integer,       intent(in) :: z
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: w, p0, p1, p2, p3, la
    integer :: i
    call grid_position(T, i, w)
    p1 = hyd_logaB(i, z);  p2 = hyd_logaB(i+1, z)
    !--- at the table ends the missing neighbor is the linear extrapolation,
    !--- so the slope is continuous there too
    if (i > 1) then
       p0 = hyd_logaB(i-1, z)
    else
       p0 = 2.0_wp*p1 - p2
    end if
    if (i < HYD_NT - 1) then
       p3 = hyd_logaB(i+2, z)
    else
       p3 = 2.0_wp*p2 - p1
    end if
    la = p1 + 0.5_wp*w*(p2 - p0 + w*(2.0_wp*p0 - 5.0_wp*p1 + 4.0_wp*p2 - p3 &
                                + w*(3.0_wp*(p1 - p2) + p3 - p0)))
    a = 10.0_wp**la
  end function hydrogenic_alpha_B

  !=========================================================================
  ! Share of the captures into n >= 2 that ends in 2s; lyman_case 'A' or 'B'.
  !=========================================================================
  pure real(kind=wp) function hydrogenic_f_2s(z, T, lyman_case) result(f)
    integer,          intent(in) :: z
    real(kind=wp),    intent(in) :: T
    character(len=*), intent(in) :: lyman_case
    real(kind=wp) :: w
    integer :: i, ic
    ic = 2
    if (lyman_case == 'A') ic = 1
    call grid_position(T, i, w)
    f = (1.0_wp - w)*hyd_f2s(i, ic, z) + w*hyd_f2s(i+1, ic, z)
  end function hydrogenic_f_2s

  !--- The cascade's own alpha_eff(2s) [cm^3 s^-1], for comparison with other
  !--- determinations; the code applies f_2s to recomb_mod's alpha_B.
  pure real(kind=wp) function hydrogenic_alpha_2s(z, T, lyman_case) result(a)
    integer,          intent(in) :: z
    real(kind=wp),    intent(in) :: T
    character(len=*), intent(in) :: lyman_case
    a = hydrogenic_f_2s(z, T, lyman_case)*hydrogenic_alpha_B(z, T)
  end function hydrogenic_alpha_2s

  !=========================================================================
  ! He II: share f(k) of the captures into n >= 2 that leaves by the Lyman
  ! line n = k + 1 -> 1, k = 1..HEII_NLY-1, and by all n >= HEII_NLY + 1
  ! together in k = HEII_NLY.  Case B: Ly-alpha alone, 1 - f_2s.
  !=========================================================================
  pure subroutine heii_lyman_shares(T, lyman_case, f)
    real(kind=wp),    intent(in)  :: T
    character(len=*), intent(in)  :: lyman_case
    real(kind=wp),    intent(out) :: f(HEII_NLY)
    real(kind=wp) :: w
    integer :: i
    f = 0.0_wp
    if (lyman_case == 'A') then
       call grid_position(T, i, w)
       f = (1.0_wp - w)*heii_lyman_fA(:, i) + w*heii_lyman_fA(:, i+1)
    else
       f(1) = 1.0_wp - hydrogenic_f_2s(2, T, 'B')
    end if
  end subroutine heii_lyman_shares

  !--- Energy [eV] of He II Lyman line k (n = k + 1 -> 1); the n >= 31 bin
  !--- (k = HEII_NLY) is put at n = 31, within 0.06 eV of every line in it.
  elemental real(kind=wp) function heii_lyman_energy(k) result(E)
    integer, intent(in) :: k
    real(kind=wp) :: n
    n = real(k + 1, wp)
    E = E_HEII_LYA*(1.0_wp - 1.0_wp/(n*n))/0.75_wp
  end function heii_lyman_energy

end module hydrogenic_recombination_mod
