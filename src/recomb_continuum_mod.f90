module recomb_continuum_mod
!---------------------------------------------------------------------------
! MoCHII: free-bound + free-free continuum coefficients (port of the MOCASSIN
! emission_mod::fb_ff tables and expressions).
!
! Per (n_e n_ion), in units of 1e-40 erg cm^3 s^-1 Hz^-1:
!  - BELOW each species' first ionization edge: the tabulated coefficients
!    data/gamma{HI,HeI,HeII}.dat (21 T x paired-edge nu points; log-log T
!    interpolation, linear-in-nu interpolation of log gamma between edge
!    pairs, as MOCASSIN fb_ff), plus the ranges the tables leave out (the
!    gaps around each edge, the stretch up to the ground edge, and below the
!    first node), which returned 0 until 2026-09-26; see gamma_tab.
!    VERIFIED (2026-07-12, tools check vs a hydrogenic Kramers Milne sum):
!    the tables are free-bound PLUS free-free combined — do not add ff again
!    inside the table range.  The Ercolano & Storey (2006) machine-readable
!    swap remains a drop-in at this interface.
!  - ABOVE the table range: Milne-relation free-bound from the ground
!    photoionization cross sections (photo_xsec: exact hydrogenic for H I
!    and He II, VFKY96 for He I) + free-free with the
!    Hummer getGauntFF Gaunt factor (Z=1 for H II/He II, Z=2 for He III).
!
! These coefficients have two consumers that must not drift apart: the
! written nebular continuum (nebcont_mod) and the He II n>=2 free-bound
! component of the diffuse ionizing field (heii_cascade_mod).  They live here
! so that both read one table through one interpolator.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: recomb_continuum_setup, recomb_continuum_is_ready
  public :: gamma_tab, gamma_milne, gamma_ff, gamma_tab_nodes
  public :: RYD_ERG, HCRYD_K, MILNE_C, NU_HI, NU_HEI, NU_HEII

  interface
     subroutine getGauntFF(z, log10Te, xlf, g, iflag)
        implicit none
        integer, intent(out) :: iflag
        real, intent(in)  :: log10Te, z
        real, dimension(:), intent(in)  :: xlf
        real, dimension(size(xlf)), intent(out) :: g
     end subroutine getGauntFF
  end interface

  !--- tabulated continuum coefficients (log10 of 1e-40 erg cm^3/s/Hz)
  type gamma_table
     integer :: nT = 0, nnu = 0
     real(kind=wp), allocatable :: tk(:), nu(:), lg(:,:)   ! (nnu, nT)
  end type gamma_table
  type(gamma_table) :: gtab(3)          ! 1=HI, 2=HeI, 3=HeII
  logical, save :: ready = .false.

  real(kind=wp), parameter :: RYD_ERG   = 2.1798723611035e-11_wp
  real(kind=wp), parameter :: HCRYD_K   = 157807.0_wp       ! Ryd/k [K]
  real(kind=wp), parameter :: MILNE_C   = 4.9874105e-6_wp   ! (h^2/2 pi me k)^1.5 Ryd^2 [cm K^1.5]
  real(kind=wp), parameter :: NU_HI     = 0.999466_wp       ! edges [Ryd]
  real(kind=wp), parameter :: NU_HEI    = 1.80804_wp
  real(kind=wp), parameter :: NU_HEII   = 4.001067_wp

contains

  !=========================================================================
  logical function recomb_continuum_is_ready() result(ok)
    ok = ready
  end function recomb_continuum_is_ready

  !=========================================================================
  ! Read the three tables once.  Idempotent, so any consumer may call it.
  !=========================================================================
  subroutine recomb_continuum_setup()
    use utility, only : fatal_error
    implicit none
    character(len=16) :: fn(3)
    integer :: k, i, unit, ios

    if (ready) return
    fn = [character(len=16) :: 'gammaHI.dat', 'gammaHeI.dat', 'gammaHeII.dat']
    do k = 1, 3
       open(newunit=unit, file=trim(par%data_dir)//'/'//trim(fn(k)), &
            status='old', iostat=ios)
       if (ios /= 0) call fatal_error('cannot open ' // &
          trim(par%data_dir)//'/'//trim(fn(k)) // ' (par%data_dir)')
       read(unit,*) gtab(k)%nT, gtab(k)%nnu
       allocate(gtab(k)%tk(gtab(k)%nT), gtab(k)%nu(gtab(k)%nnu), &
                gtab(k)%lg(gtab(k)%nnu, gtab(k)%nT))
       read(unit,*) gtab(k)%tk
       do i = 1, gtab(k)%nnu
          read(unit,*) gtab(k)%nu(i), gtab(k)%lg(i,:)
       end do
       close(unit)
       gtab(k)%lg = log10(max(gtab(k)%lg, 1.0e-30_wp))
    end do
    ready = .true.
  end subroutine recomb_continuum_setup

  !=========================================================================
  ! fb+ff coefficient of recombining species k (1 = H+, 2 = He+, 3 = He2+)
  ! at (nu, T) BELOW its ground edge (0 at and above it, where the caller
  ! takes the Milne ground continuum plus free-free).  Inside the table:
  ! linear log-gamma interpolation within an edge pair, log-log linear in T
  ! (clamped) -- the fb_ff scheme.  Three ranges the table does not cover,
  ! all of which used to return 0, i.e. no emission at all:
  !  - the gaps between pairs.  Each pair runs from just above one edge to
  !    just below the next, so the nodes (rounded to 1e-4 Ryd) leave a gap of
  !    about 2e-4 Ryd with the edge inside it: 10 in the H+ table, 7 in He+,
  !    1 in He2+, and five of the written-continuum frequencies fell in one.
  !    The edge is placed at the gap's midpoint, within 1e-4 Ryd of the true
  !    edge (the precision of the nodes themselves; for H at 1/121 Ryd it is
  !    3.5e-5 Ryd off), and each side extends its own pair.
  !  - from the last node to the ground edge: the last pair is extended
  !    (H+ 0.9994-0.99947, He+ 1.7999-1.8080, He2+ 3.9996-4.0011 Ryd).
  !  - below the first node.  He+ (table from 0.0911 Ryd, 1.0 micron): the
  !    levels whose edges lie there (n >= 4) are nearly hydrogenic, and a
  !    He+ capture into them weighs g_n/g_ion = 4n^2/2, the same 2n^2 as for
  !    H+, so the coefficient per ion is the H+ one; at the join the two
  !    tables agree to 0.4% (1e4 K) and 2.2% (5000 K).  Below the H+ table
  !    (0.0069 Ryd) and for He2+ below its table (0.0278 Ryd), free-free plus
  !    the captures into every level whose edge lies below h nu, summed in
  !    the Kramers approximation (bound-free Gaunt factor 1, the sum over n
  !    as an integral): with n1 the lowest such level and chi_n1 its binding
  !    energy, sum_{n >= n1} fb = [exp(chi_n1/kT) - 1] ff / g_ff
  !    (high_level_capture).  Just below the table's first node, 1e4 K, it
  !    is -0.3% (H+, 0.00689 against 0.0069 Ryd) and +3.6% (He2+, 0.02779
  !    against 0.0278 Ryd) off the table; free-free alone would be -6.6% and
  !    -30%.
  !=========================================================================
  recursive real(kind=wp) function gamma_tab(k, nu, T) result(g)
    implicit none
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: nu, T
    real(kind=wp), parameter :: NU_GROUND(3) = [NU_HI, NU_HEI, NU_HEII]
    integer :: i, ilo, n

    g = 0.0_wp
    n = gtab(k)%nnu
    if (nu >= NU_GROUND(k)) return
    if (nu < gtab(k)%nu(1)) then
       select case (k)
       case (1)
          g = high_level_capture(nu, T, 1)
       case (2)
          g = gamma_tab(1, nu, T)
       case (3)
          g = high_level_capture(nu, T, 2)
       end select
       return
    end if
    !--- the pair whose span, extended halfway across the gaps on either
    !--- side (and to the ground edge above the last node), contains nu
    ilo = n - 1
    do i = 1, n/2 - 1
       if (nu < 0.5_wp*(gtab(k)%nu(2*i) + gtab(k)%nu(2*i+1))) then
          ilo = 2*i - 1
          exit
       end if
    end do
    g = pair_value(k, ilo, nu, T)
  end function gamma_tab

  !--- Free-free plus Kramers free-bound into all levels n >= n1 of a
  !--- hydrogenic ion of charge Z, n1 the lowest level with Z^2/n^2 <= nu
  !--- [Ryd]: sum_n (2 chi_n/(n kT)) exp(chi_n/kT) taken as an integral over
  !--- u = chi_n/kT gives exp(chi_n1/kT) - 1 times the Kramers free-free,
  !--- i.e. that factor over the free-free Gaunt factor times gamma_ff.
  real(kind=wp) function high_level_capture(nu, T, Z) result(g)
    real(kind=wp), intent(in) :: nu, T
    integer,       intent(in) :: Z
    real(kind=wp) :: ff, gff, chi1
    integer :: n1
    ff  = gamma_ff(nu, T, Z)
    gff = ff/(6.8391e-38_wp*real(Z*Z,wp)*exp(-nu*HCRYD_K/T)/sqrt(T)*1.0e40_wp)
    n1  = ceiling(real(Z,wp)/sqrt(nu) - 1.0e-9_wp)
    chi1 = real(Z*Z,wp)/real(n1*n1,wp)
    g = ff*(1.0_wp + (exp(chi1*HCRYD_K/T) - 1.0_wp)/gff)
  end function high_level_capture

  !--- log gamma linear in nu through the pair (ilo, ilo+1), log-log in T
  !--- (clamped); used inside the pair and, extended, across a gap.
  real(kind=wp) function pair_value(k, ilo, nu, T) result(g)
    integer,       intent(in) :: k, ilo
    real(kind=wp), intent(in) :: nu, T
    integer :: iT
    real(kind=wp) :: w, lg1, lg2, lg
    if (T <= gtab(k)%tk(1)) then
       iT = 1;  w = 0.0_wp
    else if (T >= gtab(k)%tk(gtab(k)%nT)) then
       iT = gtab(k)%nT - 1;  w = 1.0_wp
    else
       do iT = 1, gtab(k)%nT - 1
          if (T < gtab(k)%tk(iT+1)) exit
       end do
       w = log10(T/gtab(k)%tk(iT)) / log10(gtab(k)%tk(iT+1)/gtab(k)%tk(iT))
    end if
    lg1 = gtab(k)%lg(ilo,  iT)*(1.0_wp - w) + gtab(k)%lg(ilo,  iT+1)*w
    lg2 = gtab(k)%lg(ilo+1,iT)*(1.0_wp - w) + gtab(k)%lg(ilo+1,iT+1)*w
    lg  = lg2 + (lg1 - lg2)*(gtab(k)%nu(ilo+1) - nu) &
          / (gtab(k)%nu(ilo+1) - gtab(k)%nu(ilo))
    g = 10.0_wp**lg
  end function pair_value

  !--- The frequencies [Ryd] where gamma_tab of table k is not smooth: the
  !--- table nodes and the edges gamma_tab places at the gap midpoints.  A
  !--- quadrature of it must break there.
  function gamma_tab_nodes(k) result(nu)
    integer, intent(in) :: k
    real(kind=wp), allocatable :: nu(:)
    integer :: i, n
    n = gtab(k)%nnu
    allocate(nu(n + n/2 - 1))
    nu(1:n) = gtab(k)%nu
    do i = 1, n/2 - 1
       nu(n+i) = 0.5_wp*(gtab(k)%nu(2*i) + gtab(k)%nu(2*i+1))
    end do
  end function gamma_tab_nodes

  !=========================================================================
  ! Milne free-bound to the ground level above the edge (fb_ff expression):
  ! gamma = 4 pi sigma statW hcRyd MILNE_C nu^3/T^1.5 exp(-(nu-nu0) Ryd/kT)
  ! in 1e-40 units.  statW = g_ground(recombined)/(2 g_ground(ion)).
  !=========================================================================
  real(kind=wp) function gamma_milne(nu, nu0, sigma_cm2, statW, T) result(g)
    implicit none
    real(kind=wp), intent(in) :: nu, nu0, sigma_cm2, statW, T
    g = 0.0_wp
    if (nu < nu0) return
    g = fourpi*sigma_cm2*statW*RYD_ERG*MILNE_C*nu**3/(T*sqrt(T)) &
        * exp(-(nu - nu0)*HCRYD_K/T) * 1.0e40_wp
  end function gamma_milne

  !=========================================================================
  ! Free-free coefficient [1e-40 erg cm^3/s/Hz]:
  ! gamma_ff = 6.8391e-38 Z^2 g_ff exp(-h nu/kT)/sqrt(T) * 1e40
  !=========================================================================
  real(kind=wp) function gamma_ff(nu, T, Z) result(g)
    implicit none
    real(kind=wp), intent(in) :: nu, T
    integer,       intent(in) :: Z
    real :: xlf(1), gf(1)
    integer :: iflag
    if (par%gaunt_vh14) then
       block
         use gaunt_vh14_mod, only : gaunt_vh14_setup, gauntff_vh14
         call gaunt_vh14_setup()
         gf(1) = real(gauntff_vh14(real(Z,wp), T, nu))
       end block
    else
       xlf(1) = real(log10(nu))
       call getGauntFF(real(log10(real(Z,wp))), real(log10(T)), xlf, gf, iflag)
    end if
    g = 6.8391e-38_wp*real(Z*Z,wp)*real(gf(1),wp) &
        * exp(-nu*HCRYD_K/T)/sqrt(T) * 1.0e40_wp
  end function gamma_ff

end module recomb_continuum_mod
