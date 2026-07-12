module cooling_mod
!---------------------------------------------------------------------------
! MoCHII: gas cooling for the G2 thermal balance.
!
! Components (all per unit volume, erg cm^-3 s^-1, assembled in
! cooling_total):
!   - recombination cooling: Hui & Gnedin (1997) case A/B fits for H II
!     and He III (hydrogenic scaling); He II as k_B T alpha (their
!     recommendation; DR negligible at HII-region temperatures);
!   - free-free: Lambda_ff = 1.42554e-27 sqrt(T) Z^2 n_ion n_e gbar(T,Z),
!     with the thermally averaged Gaunt factor gbar(T,Z) obtained at setup
!     by integrating the ported MOCASSIN getGauntFF (Hummer 1988) over
!     u = h nu / kT with weight e^-u, tabulated on a log-T grid;
!   - collisional-ionization cooling: rate x ionization potential;
!   - line cooling from the Tier-1 fits (tools/fitting products):
!     Lambda(T) = T^-1/2 sum A_i exp(-T_i/T) per (n_e n_ion), loaded from
!     data/atomic/cooling_tier1_<ion>.txt (par%atomic_dir).  Stage G2a
!     loads H I only; the metal set enters with the registry (G2b/c).
!---------------------------------------------------------------------------
  use define
  use recomb_mod
  implicit none
  private

  public :: cooling_setup, cooling_total, tier1_eval
  public :: tier1_fit_type, tier1_load, cool_HI_fit

  interface
     subroutine getGauntFF(z, log10Te, xlf, g, iflag)
        implicit none
        integer, intent(out) :: iflag
        real, intent(in)  :: log10Te, z
        real, dimension(:), intent(in)  :: xlf
        real, dimension(size(xlf)), intent(out) :: g
     end subroutine getGauntFF
  end interface

  type tier1_fit_type
     integer :: n = 0
     real(kind=wp), allocatable :: A(:), Ti(:)
  end type tier1_fit_type

  type(tier1_fit_type) :: cool_HI_fit

  !--- thermally averaged free-free Gaunt factor gbar(T) for Z=1 and Z=2,
  !--- tabulated on a log-T grid at setup.
  integer, parameter :: NGT = 41
  real(kind=wp) :: gff_logT(NGT), gff_z1(NGT), gff_z2(NGT)

contains

  !=========================================================================
  subroutine cooling_setup()
    use mpi
    implicit none
    character(len=256) :: fname
    integer :: i

    !--- Tier-1 H I line cooling (Ly-alpha dominated).
    fname = trim(par%atomic_dir)//'/cooling_tier1_h_1.txt'
    call tier1_load(trim(fname), cool_HI_fit)

    !--- gbar_ff(T, Z) tables over log T = 2..6.
    do i = 1, NGT
       gff_logT(i) = 2.0_wp + 4.0_wp*real(i-1,wp)/real(NGT-1,wp)
       gff_z1(i) = gaunt_ff_mean(10.0_wp**gff_logT(i), 1.0_wp)
       gff_z2(i) = gaunt_ff_mean(10.0_wp**gff_logT(i), 2.0_wp)
    end do

    if (mpar%p_rank == 0) then
       write(*,'(2a)')      ' COOL: Tier-1 H I fit loaded from ', trim(fname)
       write(*,'(a,f7.4)')  ' COOL: gbar_ff(1e4 K, Z=1) = ', gbar_ff(1.0e4_wp, 1)
    end if
  end subroutine cooling_setup

  !=========================================================================
  ! Read a Tier-1 coefficient file: '#' comments, then nterm, then rows
  ! (A_i [erg cm^3 s^-1 K^1/2], T_i [K]).
  !=========================================================================
  subroutine tier1_load(fname, fit)
    use mpi
    implicit none
    character(len=*),     intent(in)  :: fname
    type(tier1_fit_type), intent(out) :: fit
    character(len=256) :: line
    integer :: unit, ios, i, ierr

    open(newunit=unit, file=fname, status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(2a)') 'ERROR: cannot open ', trim(fname)
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    do
       read(unit,'(a)') line
       line = adjustl(line)
       if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
    end do
    read(line,*) fit%n
    allocate(fit%A(fit%n), fit%Ti(fit%n))
    do i = 1, fit%n
       read(unit,*) fit%A(i), fit%Ti(i)
    end do
    close(unit)
  end subroutine tier1_load

  !=========================================================================
  elemental real(kind=wp) function tier1_eval(fit, T) result(lam)
    type(tier1_fit_type), intent(in) :: fit
    real(kind=wp),        intent(in) :: T
    integer :: i
    lam = 0.0_wp
    do i = 1, fit%n
       lam = lam + fit%A(i)*exp(-fit%Ti(i)/T)
    end do
    lam = lam/sqrt(T)
  end function tier1_eval

  !=========================================================================
  ! Hui & Gnedin (1997) recombination cooling [erg cm^3 s^-1].
  !=========================================================================
  elemental real(kind=wp) function betaA_HII(T) result(b)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*157807.0_wp/T
    b = 1.778e-29_wp*T*lam**1.965_wp/(1.0_wp + (lam/0.541_wp)**0.502_wp)**2.697_wp
  end function betaA_HII

  elemental real(kind=wp) function betaB_HII(T) result(b)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*157807.0_wp/T
    b = 3.435e-30_wp*T*lam**1.970_wp/(1.0_wp + (lam/2.250_wp)**0.376_wp)**3.720_wp
  end function betaB_HII

  elemental real(kind=wp) function beta_HeII(T, caseA) result(b)
    real(kind=wp), intent(in) :: T
    logical,       intent(in) :: caseA
    !--- power-law alpha: cooling ~ k_B T alpha (Hui & Gnedin).
    if (caseA) then
       b = kboltz_cgs*T*alphaA_HeII(T)
    else
       b = kboltz_cgs*T*alphaB_HeII(T)
    end if
  end function beta_HeII

  elemental real(kind=wp) function beta_HeIII(T, caseA) result(b)
    real(kind=wp), intent(in) :: T
    logical,       intent(in) :: caseA
    real(kind=wp) :: lam
    !--- hydrogenic Z=2: the H fit on lambda(T_TR = 631515 K), x Z.
    lam = 2.0_wp*631515.0_wp/T
    if (caseA) then
       b = 2.0_wp*1.778e-29_wp*T*lam**1.965_wp/(1.0_wp + (lam/0.541_wp)**0.502_wp)**2.697_wp
    else
       b = 2.0_wp*3.435e-30_wp*T*lam**1.970_wp/(1.0_wp + (lam/2.250_wp)**0.376_wp)**3.720_wp
    end if
  end function beta_HeIII

  !=========================================================================
  ! Thermally averaged free-free Gaunt factor: gbar = Int g_ff(u) e^-u du
  ! over u = h nu / kT (32-point midpoint on u in [1e-4, 20]), using the
  ! ported MOCASSIN getGauntFF (log10 nu in Rydberg units).
  !=========================================================================
  real(kind=wp) function gaunt_ff_mean(T, Zc) result(gm)
    implicit none
    real(kind=wp), intent(in) :: T, Zc
    integer, parameter :: NU = 64
    real, allocatable :: xlf(:), g(:)
    real(kind=wp) :: u(NU), du(NU), lnu_lo, lnu_hi, ryd_over_kt
    integer :: i, iflag

    allocate(xlf(NU), g(NU))
    !--- log-spaced u grid; weight e^-u du.
    lnu_lo = log(1.0e-4_wp);  lnu_hi = log(20.0_wp)
    do i = 1, NU
       u(i)  = exp(lnu_lo + (lnu_hi - lnu_lo)*(real(i,wp)-0.5_wp)/real(NU,wp))
       du(i) = u(i)*(lnu_hi - lnu_lo)/real(NU,wp)
    end do
    !--- h nu / Ryd = u * kT / Ryd
    ryd_over_kt = 13.605693_wp*ev2erg/(kboltz_cgs*T)
    do i = 1, NU
       xlf(i) = real(log10(u(i)/ryd_over_kt))
    end do
    call getGauntFF(real(log10(Zc)), real(log10(T)), xlf, g, iflag)
    gm = 0.0_wp
    do i = 1, NU
       gm = gm + real(g(i),wp)*exp(-u(i))*du(i)
    end do
    !--- normalize by Int e^-u du over the same grid (finite-range bias).
    gm = gm / (exp(-1.0e-4_wp) - exp(-20.0_wp))
    deallocate(xlf, g)
  end function gaunt_ff_mean

  !=========================================================================
  real(kind=wp) function gbar_ff(T, Z) result(gm)
    implicit none
    real(kind=wp), intent(in) :: T
    integer,       intent(in) :: Z
    real(kind=wp) :: lt, w
    integer :: i
    lt = min(max(log10(T), gff_logT(1)), gff_logT(NGT))
    i  = min(int((lt - gff_logT(1))/(gff_logT(2) - gff_logT(1))) + 1, NGT-1)
    w  = (lt - gff_logT(i))/(gff_logT(i+1) - gff_logT(i))
    if (Z == 2) then
       gm = gff_z2(i)*(1.0_wp - w) + gff_z2(i+1)*w
    else
       gm = gff_z1(i)*(1.0_wp - w) + gff_z1(i+1)*w
    end if
  end function gbar_ff

  !=========================================================================
  ! Total H/He cooling per unit volume [erg cm^-3 s^-1] at temperature T.
  ! Metal line cooling is added by the registry hook (G2b/c).
  !=========================================================================
  real(kind=wp) function cooling_total(T, nH, ne, xHI, xHeI, xHeII) result(cool)
    implicit none
    real(kind=wp), intent(in) :: T, nH, ne, xHI, xHeI, xHeII
    real(kind=wp) :: nHI, nHII, nHeI, nHeII_n, nHeIII, xHeIII
    logical :: caseA

    caseA  = trim(par%case_ab) == 'A'
    xHeIII = max(0.0_wp, 1.0_wp - xHeI - xHeII)
    nHI    = nH*xHI
    nHII   = nH*(1.0_wp - xHI)
    nHeI   = nH*par%He_abund*xHeI
    nHeII_n= nH*par%He_abund*xHeII
    nHeIII = nH*par%He_abund*xHeIII

    !--- recombination cooling
    if (caseA) then
       cool = ne*(nHII*betaA_HII(T) + nHeII_n*beta_HeII(T, .true.) &
              + nHeIII*beta_HeIII(T, .true.))
    else
       cool = ne*(nHII*betaB_HII(T) + nHeII_n*beta_HeII(T, .false.) &
              + nHeIII*beta_HeIII(T, .false.))
    end if

    !--- free-free (Z=1: H II + He II; Z=2: He III)
    cool = cool + 1.42554e-27_wp*sqrt(T)*ne* &
           ( (nHII + nHeII_n)*gbar_ff(T,1) + 4.0_wp*nHeIII*gbar_ff(T,2) )

    !--- collisional-ionization cooling
    cool = cool + ne*( nHI*ci_HI(T)*eth_HI + nHeI*ci_HeI(T)*eth_HeI &
           + nHeII_n*ci_HeII(T)*eth_HeII )*ev2erg

    !--- H I collisional-excitation line cooling (Tier-1 fit)
    cool = cool + ne*nHI*tier1_eval(cool_HI_fit, T)
  end function cooling_total

end module cooling_mod
