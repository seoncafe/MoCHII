module cooling_mod
!---------------------------------------------------------------------------
! MoCHII: gas cooling for the thermal balance.
!
! Components (all per unit volume, erg cm^-3 s^-1, assembled in
! cooling_total):
!   - recombination cooling: beta = k T [(3/2) alpha + d alpha/d ln T]
!     from the case-A/B alpha of recomb_mod for H II, He II and He III
!     (maxwellian_recombination_cooling; exact for a Maxwellian);
!   - free-free: Lambda_ff = FF_TOTAL sqrt(T) Z^2 n_ion n_e gbar(T,Z)
!     (free_free_mod: FF_TOTAL = 1.425557e-27 cgs, the frequency integral of
!     the spectral coefficient the written continuum and the transported
!     free-free packets use), with the thermally averaged Gaunt factor
!     gbar(T,Z) obtained at setup by integrating free_free_gaunt (Hummer
!     1988 within its stated range, the van Hoof et al. 2014 table outside
!     it) over u = h nu / kT with weight e^-u, tabulated on a log-T grid;
!   - collisional-ionization cooling: rate x ionization potential;
!   - line cooling from the cooling fits (tools/fitting products):
!     Lambda(T) = T^-1/2 sum A_i exp(-T_i/T) per (n_e n_ion), loaded from
!     data/atomic/cooling_<ion>.txt (par%atomic_dir).  H I is loaded here;
!     the metal set enters with the registry.
!---------------------------------------------------------------------------
  use define
  use recomb_mod
  use nlevel_cooling_mod, only : nlevel_cooling_add, nlevel_cooling_apply
  use free_free_mod,      only : free_free_gaunt, FF_TOTAL
  implicit none
  private

  public :: cooling_setup, cooling_total, cooling_eval
  public :: cooling_fit_type, cooling_load, cool_HI_fit
  public :: gbar_ff
  !--- The recombination-cooling coefficients are public so that the analytic
  !--- gates can build their 1D references from the same expressions the code
  !--- evaluates, and so tests/rates/check_python_rates.py can hold the Python
  !--- copies to them.  A reference that carries its own copy drifts silently
  !--- when one of these changes, which is how three such mismatches reached
  !--- the published gate numbers before they were found.
  public :: betaA_HII, betaB_HII, beta_HeII, beta_HeIII

  type cooling_fit_type
     integer :: n = 0
     real(kind=wp), allocatable :: A(:), Ti(:)
  end type cooling_fit_type

  type(cooling_fit_type) :: cool_HI_fit

  !--- ln T step of the centered derivative in the recombination cooling
  real(kind=wp), parameter :: DLNT_BETA = 1.0e-4_wp

  !--- index into the nlevel_cooling_mod (T, n_e) table for H I line cooling
  !--- (par%cooling_model = 'local_ne'); 0 = no table, keep the fit.
  integer :: hi_cool_idx = 0

  !--- thermally averaged free-free Gaunt factor gbar(T,Z) for net ion
  !--- charge Z = 1..MAX_ION_STAGES-1, tabulated on a log-T grid at setup.
  !--- 0.01 dex in T: linear interpolation on the former 0.1 dex grid left
  !--- up to 4.9e-4 between nodes (Z = 4 at 3000 K, beside the Hummer / van
  !--- Hoof switch at gamma^2 = 1e3; tests/freefree_emission (a))
  integer, parameter :: NGT = 401, NGZ = MAX_ION_STAGES - 1
  real(kind=wp) :: gff_logT(NGT), gff_z(NGT, NGZ)

contains

  !=========================================================================
  subroutine cooling_setup()
    use mpi
    implicit none
    character(len=PATH_LEN) :: fname
    integer :: i, iz

    !--- H I line cooling (Ly-alpha dominated).
    fname = trim(par%atomic_dir)//'/cooling_h_1.txt'
    call cooling_load(trim(fname), cool_HI_fit)

    !--- (T, n_e) suppression table for H I under the local-n_e cooling model.
    if (trim(par%cooling_model) == 'local_ne') hi_cool_idx = nlevel_cooling_add('h', 1)

    !--- gbar_ff(T, Z) tables over log T = 2..6, through the highest
    !--- represented ion charge.
    do i = 1, NGT
       gff_logT(i) = 2.0_wp + 4.0_wp*real(i-1,wp)/real(NGT-1,wp)
       do iz = 1, NGZ
          gff_z(i,iz) = gaunt_ff_mean(10.0_wp**gff_logT(i), real(iz,wp))
       end do
    end do

    if (mpar%p_rank == 0) then
       write(*,'(2a)')      ' COOL: H I cooling fit loaded from ', trim(fname)
       write(*,'(a,f7.4)')  ' COOL: gbar_ff(1e4 K, Z=1) = ', gbar_ff(1.0e4_wp, 1)
    end if
  end subroutine cooling_setup

  !=========================================================================
  ! Read a cooling coefficient file: '#' comments, then nterm, then rows
  ! (A_i [erg cm^3 s^-1 K^1/2], T_i [K]).
  !=========================================================================
  subroutine cooling_load(fname, fit)
    use mpi
        use utility, only : fatal_error
implicit none
    character(len=*),       intent(in)  :: fname
    type(cooling_fit_type), intent(out) :: fit
    character(len=256) :: line
    integer :: unit, ios, i, ierr

    open(newunit=unit, file=fname, status='old', iostat=ios)
    if (ios /= 0) then
       call fatal_error('cannot open '//trim(fname))
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
  end subroutine cooling_load

  !=========================================================================
  elemental real(kind=wp) function cooling_eval(fit, T) result(lam)
    type(cooling_fit_type), intent(in) :: fit
    real(kind=wp),          intent(in) :: T
    integer :: i
    lam = 0.0_wp
    do i = 1, fit%n
       lam = lam + fit%A(i)*exp(-fit%Ti(i)/T)
    end do
    !--- A fit with coefficients of both signs can go negative outside the
    !--- temperatures it was fitted on (e.g. a Mg III fit below 1000 K); a
    !--- collisional line cooling coefficient is never negative, so the sum
    !--- is held at zero there.  Inside the fit range it is positive.
    lam = max(lam, 0.0_wp)/sqrt(T)
  end function cooling_eval

  !=========================================================================
  ! Recombination cooling [erg cm^3 s^-1]: the kinetic energy the recombining
  ! electrons remove, beta = < E sigma_rec v >.  For a Maxwellian at T and any
  ! capture cross section sigma(E) -- every level, radiative and
  ! dielectronic alike -- differentiating alpha = < sigma v > with respect to
  ! T gives exactly
  !     beta(T) = k T [ (3/2) alpha(T) + d alpha / d ln T ],
  ! so beta is evaluated from the SAME alpha the ionization balance uses
  ! (recomb_mod), case A or B, and the two cannot disagree.  The derivative
  ! is a centered difference in ln T (step 1e-4: truncation ~1e-8 relative).
  ! This replaced (2026-09-24) the Hui & Gnedin (1997) fits for H II and
  ! He III, which fit a different alpha (Ferland et al. 1992; 2.3% below the
  ! identity on MoCHII's case-B alpha at 1e4 K), and k T alpha for He II, which
  ! is the mean kinetic energy only for alpha ~ T^-1/2 and overstated the
  ! He II recombination cooling 1.43x at 1e4 K (the mean energy removed is
  ! 0.70 kT there).
  !=========================================================================
  elemental real(kind=wp) function maxwellian_recombination_cooling(T, am, a0, ap) &
                                                                    result(b)
    real(kind=wp), intent(in) :: T, am, a0, ap     ! alpha at T e^-h, T, T e^+h
    b = kboltz_cgs*T*(1.5_wp*a0 + (ap - am)/(2.0_wp*DLNT_BETA))
  end function maxwellian_recombination_cooling

  elemental real(kind=wp) function betaA_HII(T) result(b)
    real(kind=wp), intent(in) :: T
    b = maxwellian_recombination_cooling(T, alphaA_HII(T*exp(-DLNT_BETA)), &
           alphaA_HII(T), alphaA_HII(T*exp(DLNT_BETA)))
  end function betaA_HII

  elemental real(kind=wp) function betaB_HII(T) result(b)
    real(kind=wp), intent(in) :: T
    b = maxwellian_recombination_cooling(T, alphaB_HII(T*exp(-DLNT_BETA)), &
           alphaB_HII(T), alphaB_HII(T*exp(DLNT_BETA)))
  end function betaB_HII

  elemental real(kind=wp) function beta_HeII(T, caseA) result(b)
    real(kind=wp), intent(in) :: T
    logical,       intent(in) :: caseA
    if (caseA) then
       b = maxwellian_recombination_cooling(T, alphaA_HeII(T*exp(-DLNT_BETA)), &
              alphaA_HeII(T), alphaA_HeII(T*exp(DLNT_BETA)))
    else
       b = maxwellian_recombination_cooling(T, alphaB_HeII(T*exp(-DLNT_BETA)), &
              alphaB_HeII(T), alphaB_HeII(T*exp(DLNT_BETA)))
    end if
  end function beta_HeII

  elemental real(kind=wp) function beta_HeIII(T, caseA) result(b)
    real(kind=wp), intent(in) :: T
    logical,       intent(in) :: caseA
    if (caseA) then
       b = maxwellian_recombination_cooling(T, alphaA_HeIII(T*exp(-DLNT_BETA)), &
              alphaA_HeIII(T), alphaA_HeIII(T*exp(DLNT_BETA)))
    else
       b = maxwellian_recombination_cooling(T, alphaB_HeIII(T*exp(-DLNT_BETA)), &
              alphaB_HeIII(T), alphaB_HeIII(T*exp(DLNT_BETA)))
    end if
  end function beta_HeIII

  !=========================================================================
  ! Thermally averaged free-free Gaunt factor: gbar = Int g_ff(u) e^-u du
  ! over u = h nu / kT, g_ff from free_free_mod::free_free_gaunt (Hummer
  ! 1988 inside its stated range 1e-3 <= gamma^2 <= 1e3, 1e-4 <= u <=
  ! 10^1.5, T >= 10^2.5 K; the van Hoof et al. 2014 table outside it), the
  ! Gaunt factor of the spectral coefficient as well, so the cooling and the
  ! emitted spectrum integrate one function.  Trapezoid rule in ln u over
  ! [1e-9, 300] (400 nodes): the integrand g u e^-u vanishes at both ends,
  ! so the rule converges exponentially.  It was a 64-point midpoint rule on
  ! [1e-4, 20] normalized by Int e^-u du, whose missing u < 1e-4 (where g is
  ! largest) and midpoint error left gbar 2.8e-4 below the integral of the
  ! spectrum (tests/freefree_emission (a)).
  !=========================================================================
  real(kind=wp) function gaunt_ff_mean(T, Zc) result(gm)
    implicit none
    real(kind=wp), intent(in) :: T, Zc
    integer, parameter :: NU = 400
    real(kind=wp), parameter :: U_LO = 1.0e-9_wp, U_HI = 300.0_wp
    real(kind=wp) :: lu, h, u, f
    integer :: i

    h  = (log(U_HI) - log(U_LO))/real(NU-1, wp)
    gm = 0.0_wp
    do i = 1, NU
       lu = log(U_LO) + h*real(i-1, wp)
       u  = exp(lu)
       f  = free_free_gaunt(Zc, T, u)*u*exp(-u)
       if (i == 1 .or. i == NU) f = 0.5_wp*f
       gm = gm + f*h
    end do
  end function gaunt_ff_mean

  !=========================================================================
  real(kind=wp) function gbar_ff(T, Z) result(gm)
    implicit none
    real(kind=wp), intent(in) :: T
    integer,       intent(in) :: Z
    real(kind=wp) :: lt, w
    integer :: i, iz
    iz = max(1, min(Z, NGZ))
    lt = min(max(log10(T), gff_logT(1)), gff_logT(NGT))
    i  = min(int((lt - gff_logT(1))/(gff_logT(2) - gff_logT(1))) + 1, NGT-1)
    w  = (lt - gff_logT(i))/(gff_logT(i+1) - gff_logT(i))
    gm = gff_z(i,iz)*(1.0_wp - w) + gff_z(i+1,iz)*w
  end function gbar_ff

  !=========================================================================
  ! Total H/He cooling per unit volume [erg cm^-3 s^-1] at temperature T.
  ! Metal line cooling is added by the registry hook.  x_HII and x_HeIII
  ! are passed beside x_HI and x_HeI, x_HeII (gas_state_mod): beyond the
  ! ionization fronts 1 - x_HI and 1 - x_HeI - x_HeII do not represent them.
  !=========================================================================
  real(kind=wp) function cooling_total(T, nH, ne, xHI, xHII, xHeI, xHeII, &
                                       xHeIII) result(cool)
    implicit none
    real(kind=wp), intent(in) :: T, nH, ne, xHI, xHII, xHeI, xHeII, xHeIII
    real(kind=wp) :: nHI, nHII, nHeI, nHeII_n, nHeIII
    logical :: caseA

    caseA  = trim(par%recombination_continuum) /= 'on_the_spot'
    nHI    = nH*xHI
    nHII   = nH*xHII
    nHeI   = nH*par%abund_He*xHeI
    nHeII_n= nH*par%abund_He*xHeII
    nHeIII = nH*par%abund_He*xHeIII

    !--- recombination cooling
    if (caseA) then
       cool = ne*(nHII*betaA_HII(T) + nHeII_n*beta_HeII(T, .true.) &
              + nHeIII*beta_HeIII(T, .true.))
    else
       cool = ne*(nHII*betaB_HII(T) + nHeII_n*beta_HeII(T, .false.) &
              + nHeIII*beta_HeIII(T, .false.))
    end if

    !--- free-free (Z=1: H II + He II; Z=2: He III)
    cool = cool + FF_TOTAL*sqrt(T)*ne* &
           ( (nHII + nHeII_n)*gbar_ff(T,1) + 4.0_wp*nHeIII*gbar_ff(T,2) )

    !--- collisional-ionization cooling
    cool = cool + ne*( nHI*ci_HI(T)*eth_HI + nHeI*ci_HeI(T)*eth_HeI &
           + nHeII_n*ci_HeII(T)*eth_HeII )*ev2erg

    !--- H I collisional-excitation line cooling (cooling fit, suppressed by
    !--- the level populations at n_e when par%cooling_model = 'local_ne';
    !--- hi_cool_idx = 0 otherwise -> the fit value returned unchanged).
    cool = cool + ne*nHI*nlevel_cooling_apply(hi_cool_idx, &
           cooling_eval(cool_HI_fit, T), T, ne)
  end function cooling_total

end module cooling_mod
