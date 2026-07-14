module nebcont_mod
!---------------------------------------------------------------------------
! MoCHII: nebular continuum emission (port of MOCASSIN emission_mod::fb_ff
! + twoPhoton/HeI2photSub).
!
! Components, per (n_e n_ion), in units of 1e-40 erg cm^3 s^-1 Hz^-1:
!  - BELOW each species' first ionization edge: the tabulated continuum
!    coefficients data/gamma{HI,HeI,HeII}.dat (21 T x paired-edge nu
!    points; log-log T interpolation, linear-in-nu interpolation of
!    log gamma between edge pairs, exactly as MOCASSIN fb_ff).
!    VERIFIED (2026-07-12, tools check vs a hydrogenic Kramers Milne sum):
!    the tables are free-bound PLUS free-free combined — do not add ff
!    again inside the table range.  The Ercolano & Storey (2006)
!    machine-readable swap remains a drop-in at this interface.
!  - ABOVE the table range: Milne-relation free-bound from the ground
!    photoionization cross sections (photo_xsec VFKY96) + free-free with
!    the Hummer getGauntFF Gaunt factor (Z=1 for H II/He II, Z=2 for
!    He III), as in fb_ff.
!  - Two-photon: H I 2s (Nussbaumer & Schmutz 1984 spectral shape;
!    Pengelly 1964 alpha_eff(2s); Osterbrock q(2^2S-2^2P) collisional
!    suppression), He II 2s (NS84 Z-scaled; Storey & Hummer 1995
!    alpha_eff), He I 2^1S (Almog & Netzer 1989 A = 51.3 s^-1 shape from
!    data/HeI2phot.dat; Benjamin, Skillman & Smits 1999 alpha_eff).
!
! The MOCASSIN origin reallocates its logGamma scratch on every call; here
! everything is read and stored once at setup (nebcont_setup).
! Output (nebcont_write): '<base>_nebcont.txt' — the grid-integrated
! continuum L_nu on a log frequency grid, total and per component.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: nebcont_setup, nebcont_write

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

  !--- He I two-photon shape table (y, A(y)) from HeI2phot.dat
  integer :: nhe2q = 0
  real(kind=wp), allocatable :: he2q_y(:), he2q_A(:)

  !--- output frequency grid [Ryd]
  integer, parameter :: NNU_OUT = 400
  real(kind=wp) :: nu_out(NNU_OUT)

  real(kind=wp), parameter :: RYD_ERG   = 2.1798723611035e-11_wp
  real(kind=wp), parameter :: HCRYD_K   = 157807.0_wp       ! Ryd/k [K]
  real(kind=wp), parameter :: MILNE_C   = 4.9874105e-6_wp   ! (h^2/2 pi me k)^1.5 Ryd^2 [cm K^1.5]
  real(kind=wp), parameter :: NU_HI     = 0.999466_wp       ! edges [Ryd]
  real(kind=wp), parameter :: NU_HEI    = 1.80804_wp
  real(kind=wp), parameter :: NU_HEII   = 4.001067_wp

contains

  !=========================================================================
  subroutine nebcont_setup()
    use mpi
    implicit none
    character(len=16) :: fn(3)
    integer :: k, i, unit, ios, ierr

    fn = [character(len=16) :: 'gammaHI.dat', 'gammaHeI.dat', 'gammaHeII.dat']
    do k = 1, 3
       open(newunit=unit, file='../../data/'//trim(fn(k)), status='old', iostat=ios)
       if (ios /= 0) &
          open(newunit=unit, file='data/'//trim(fn(k)), status='old', iostat=ios)
       if (ios /= 0) then
          if (mpar%p_rank == 0) write(*,'(2a)') 'ERROR: cannot open ', trim(fn(k))
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
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

    !--- He I two-photon shape
    open(newunit=unit, file='../../data/HeI2phot.dat', status='old', iostat=ios)
    if (ios /= 0) open(newunit=unit, file='data/HeI2phot.dat', status='old', iostat=ios)
    if (ios == 0) then
       nhe2q = 0
       allocate(he2q_y(64), he2q_A(64))
       do
          read(unit,*,iostat=ios) he2q_y(nhe2q+1), he2q_A(nhe2q+1)
          if (ios /= 0) exit
          nhe2q = nhe2q + 1
          if (nhe2q >= 64) exit
       end do
       close(unit)
    end if

    !--- output grid 0.008 - 4.2 Ryd (lambda ~ 11.4 um - 217 A)
    do i = 1, NNU_OUT
       nu_out(i) = exp(log(0.008_wp) + (log(4.2_wp) - log(0.008_wp)) &
                   *real(i-1,wp)/real(NNU_OUT-1,wp))
    end do
    if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
       ' NEBC: continuum tables loaded (', nhe2q, ' He I 2q shape points)'
  end subroutine nebcont_setup

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

  !=========================================================================
  ! Two-photon coefficients (1e-40 erg cm^3/s/Hz), MOCASSIN formulas.
  !=========================================================================
  real(kind=wp) function gamma_2q_HI(nu, T, ne) result(g)
    implicit none
    real(kind=wp), intent(in) :: nu, T, ne
    real(kind=wp) :: y, Ay, gNu, a2s, q, fac
    real(kind=wp), parameter :: nu0 = 0.7496_wp
    g = 0.0_wp
    y = nu/nu0
    if (y >= 1.0_wp) return
    a2s = 0.8368_wp*(T*1.0e-4_wp)**(-0.723_wp)      ! [1e-13 cm^3/s], Pengelly 64
    Ay  = 202.0_wp*(y*(1.0_wp-y)*(1.0_wp-(4.0_wp*y*(1.0_wp-y))**0.8_wp) &
          + 0.88_wp*((y*(1.0_wp-y))**1.53_wp)*(4.0_wp*y*(1.0_wp-y))**0.8_wp)
    gNu = h_planck_cgs*y*Ay/8.2249_wp               ! [erg/Hz]
    if (T <= 1.0e4_wp) then
       q = 5.31e-4_wp
    else if (T >= 2.0e4_wp) then
       q = 4.71e-4_wp
    else
       fac = log10(4.71e-4_wp/5.31e-4_wp)/log10(2.0_wp)
       q = 10.0_wp**(log10(5.31e-4_wp) + fac*log10(T/1.0e4_wp))
    end if
    g = a2s*1.0e-13_wp*gNu/(1.0_wp + ne*q/8.23_wp) * 1.0e40_wp
  end function gamma_2q_HI

  real(kind=wp) function gamma_2q_HeII(nu, T) result(g)
    implicit none
    real(kind=wp), intent(in) :: nu, T
    real(kind=wp) :: y, Ay, gNu, a2s, fac
    real(kind=wp), parameter :: nu0 = 3.00_wp
    g = 0.0_wp
    y = nu/nu0
    if (y >= 1.0_wp) return
    !--- SH95 alpha_eff(2s) [1e-13 cm^3/s] interpolation (ne = 100 values)
    if (T <= 5000.0_wp) then
       a2s = 6.161_wp
    else if (T >= 30000.0_wp) then
       a2s = 2.035_wp
    else if (T <= 10000.0_wp) then
       fac = log10(4.091_wp/6.161_wp)/log10(2.0_wp)
       a2s = 10.0_wp**(log10(6.161_wp) + fac*log10(T/5000.0_wp))
    else if (T <= 15000.0_wp) then
       fac = log10(3.189_wp/4.091_wp)/log10(1.5_wp)
       a2s = 10.0_wp**(log10(4.091_wp) + fac*log10(T/10000.0_wp))
    else
       fac = log10(2.035_wp/3.189_wp)/log10(2.0_wp)
       a2s = 10.0_wp**(log10(3.189_wp) + fac*log10(T/15000.0_wp))
    end if
    Ay = 64.0_wp*0.9994667_wp*202.0_wp &
         *(y*(1.0_wp-y)*(1.0_wp-(4.0_wp*y*(1.0_wp-y))**0.8_wp) &
         + 0.88_wp*((y*(1.0_wp-y))**1.53_wp)*(4.0_wp*y*(1.0_wp-y))**0.8_wp)
    gNu = h_planck_cgs*y*Ay/(8.226_wp*64.0_wp)
    g = a2s*1.0e-13_wp*gNu * 1.0e40_wp
  end function gamma_2q_HeII

  real(kind=wp) function gamma_2q_HeI(nu, T) result(g)
    implicit none
    real(kind=wp), intent(in) :: nu, T
    real(kind=wp) :: y, Ay, a2s, w
    real(kind=wp), parameter :: nu0 = 1.514_wp
    integer :: j
    g = 0.0_wp
    y = nu/nu0
    if (y >= 1.0_wp .or. nhe2q < 2) return
    a2s = 6.23_wp*(T/1.0e4_wp)**(-0.827_wp)         ! [1e-14], BSS99
    do j = 1, nhe2q-1
       if (y <= he2q_y(j+1)) exit
    end do
    w  = (y - he2q_y(j))/(he2q_y(j+1) - he2q_y(j))
    Ay = he2q_A(j)*(1.0_wp - w) + he2q_A(j+1)*w
    if (Ay < 9.2163086e-3_wp) Ay = 0.0_wp
    !--- result in 1e-40 units: a2s [1e-14] x (h x 1e26 = 0.66262) x y Ay/A
    g = a2s*0.66262_wp*y*Ay/51.3_wp
  end function gamma_2q_HeI

  !=========================================================================
  ! Grid-integrated continuum spectrum: '<base>_nebcont.txt'.
  !=========================================================================
  subroutine nebcont_write()
    use octree_mod, only : amr_grid, leaf_half
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, &
                              gas_Te, gas_nleaf
    use photo_xsec,    only : sigma_HI, sigma_HeI, sigma_HeII
    use utility,       only : get_base_name
    implicit none
    real(kind=wp) :: Lfb(NNU_OUT), Lff(NNU_OUT), L2q(NNU_OUT)
    real(kind=wp) :: T, ne, nH, vol, nHII, nHeII_n, nHeIII, xHeIII
    real(kind=wp) :: nu, e_eV, gH, gHe1, gHe2, ff1, ff2
    character(len=192) :: outname
    integer :: il, i, unit

    if (mpar%p_rank /= 0) return
    Lfb = 0.0_wp;  Lff = 0.0_wp;  L2q = 0.0_wp

    do il = 1, gas_nleaf
       nH = gas_nH(il)
       if (nH <= 0.0_wp) cycle
       T  = gas_Te(il);  ne = gas_ne(il)
       vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
       xHeIII  = max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il))
       nHII    = nH*(1.0_wp - gas_xHI(il))
       nHeII_n = nH*par%He_abund*gas_xHeII(il)
       nHeIII  = nH*par%He_abund*xHeIII
       do i = 1, NNU_OUT
          nu = nu_out(i)
          e_eV = nu*13.605693_wp
          !--- table region: fb+ff combined; above: Milne fb + explicit ff
          gH   = gamma_tab(1, nu, T)
          gHe1 = gamma_tab(2, nu, T)
          gHe2 = gamma_tab(3, nu, T)
          if (nu >= NU_HI)   gH   = gamma_milne(nu, NU_HI,  sigma_HI(e_eV), &
                                    1.0_wp, T)
          if (nu >= NU_HEI)  gHe1 = gamma_milne(nu, NU_HEI, sigma_HeI(e_eV), &
                                    0.25_wp, T)
          if (nu >= NU_HEII) gHe2 = gamma_milne(nu, NU_HEII, sigma_HeII(e_eV), &
                                    1.0_wp, T)
          if (nu >= NU_HI) then
             ff1 = gamma_ff(nu, T, 1)
             gH = gH + ff1
             if (nu >= NU_HEI)  gHe1 = gHe1 + ff1
             if (nu >= NU_HEII) gHe2 = gHe2 + gamma_ff(nu, T, 2) - ff1
          end if
          Lfb(i) = Lfb(i) + ne*(nHII*gH + nHeII_n*gHe1 + nHeIII*gHe2) &
                   *1.0e-40_wp*vol
          !--- explicit ff component for the report (table region estimate)
          if (nu < NU_HI) then
             ff1 = gamma_ff(nu, T, 1)
             ff2 = gamma_ff(nu, T, 2)
             Lff(i) = Lff(i) + ne*((nHII + nHeII_n)*ff1 + nHeIII*ff2) &
                      *1.0e-40_wp*vol
          end if
          L2q(i) = L2q(i) + ne*( nHII*gamma_2q_HI(nu, T, ne) &
                   + nHeII_n*gamma_2q_HeI(nu, T) &
                   + nHeIII*gamma_2q_HeII(nu, T) )*1.0e-40_wp*vol
       end do
    end do

    outname = trim(get_base_name(par%out_file))//'_nebcont.txt'
    open(newunit=unit, file=trim(outname), status='replace')
    write(unit,'(a)') '# MoCHII nebular continuum (grid-integrated), '// &
       'fb_ff port (see nebcont_mod.f90 header for provenance)'
    write(unit,'(a)') '# col 1: nu [Ryd]; 2: lambda [um]; '// &
       '3: L_nu total = fb+ff(+Milne) + two-photon [erg/s/Hz]; '// &
       '4: ff estimate below the H edge; 5: two-photon'
    do i = 1, NNU_OUT
       write(unit,'(f10.5,f12.5,3es14.5)') nu_out(i), &
          0.0911267_wp/nu_out(i), Lfb(i) + L2q(i), Lff(i), L2q(i)
    end do
    close(unit)
    write(*,'(2a)') ' NEBC: nebular continuum written to: ', trim(outname)
  end subroutine nebcont_write

end module nebcont_mod
