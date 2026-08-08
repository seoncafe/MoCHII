module recomb_continuum_mod
!---------------------------------------------------------------------------
! MoCHII: free-bound + free-free continuum coefficients (port of the MOCASSIN
! emission_mod::fb_ff tables and expressions).
!
! Per (n_e n_ion), in units of 1e-40 erg cm^3 s^-1 Hz^-1:
!  - BELOW each species' first ionization edge: the tabulated coefficients
!    data/gamma{HI,HeI,HeII}.dat (21 T x paired-edge nu points; log-log T
!    interpolation, linear-in-nu interpolation of log gamma between edge
!    pairs, exactly as MOCASSIN fb_ff).
!    VERIFIED (2026-07-12, tools check vs a hydrogenic Kramers Milne sum):
!    the tables are free-bound PLUS free-free combined — do not add ff again
!    inside the table range.  The Ercolano & Storey (2006) machine-readable
!    swap remains a drop-in at this interface.
!  - ABOVE the table range: Milne-relation free-bound from the ground
!    photoionization cross sections (photo_xsec VFKY96) + free-free with the
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
  public :: gamma_tab, gamma_milne, gamma_ff
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
  ! Tabulated fb+ff coefficient for species k at (nu, T); 0 outside the
  ! table nu range.  Linear log-gamma interpolation between edge pairs,
  ! log-log linear in T (clamped) — the fb_ff scheme.
  !=========================================================================
  real(kind=wp) function gamma_tab(k, nu, T) result(g)
    implicit none
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: nu, T
    integer :: iT, i, ilo
    real(kind=wp) :: w, lg1, lg2, lg

    g = 0.0_wp
    if (nu < gtab(k)%nu(1) .or. nu >= gtab(k)%nu(gtab(k)%nnu)) return
    !--- nu segment: pairs (2i-1, 2i) span the inter-edge intervals
    ilo = 0
    do i = 1, gtab(k)%nnu/2
       if (nu >= gtab(k)%nu(2*i-1) .and. nu < gtab(k)%nu(2*i)) then
          ilo = 2*i - 1
          exit
       end if
    end do
    if (ilo == 0) return
    !--- T bracket (clamped)
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
  end function gamma_tab

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
