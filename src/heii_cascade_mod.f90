module heii_cascade_mod
!---------------------------------------------------------------------------
! MoCHII: the H-ionizing radiation a He III recombination emits on its way
! down through He II -- the fifth diffuse channel (par%heii_diffuse).
!
! Channel 3 already emits the He II GROUND continuum, alpha_1 = alpha_A -
! alpha_B, above 54.416 eV.  Everything else the recombination radiates used
! to be dropped, and unlike the He I cascade (channel 4) that omission is not
! small: for He II the n = 2 threshold sits at 13.604 eV, a hair ABOVE the
! hydrogen edge, so the whole Balmer continuum photoionizes hydrogen, and the
! pair from 2^2S spans 40.813 eV instead of He I's 20.62 eV.  Three paths:
!
!   free-bound   captures into n >= 2.  Emitted at the capture, not at the
!                n = 2 exit, so it carries its own coefficient rather than a
!                branch of alpha_B.  Taken from the same fb+ff tables the
!                written continuum uses (recomb_continuum_mod), minus the
!                free-free term, integrated from 13.598 eV to the He II
!                ground edge.
!   Ly-alpha     2^2P -> 1^2S, 40.813 eV in one photon.
!   two-photon   2^2S -> 1^2S, a pair sharing 40.813 eV; 88% of that energy
!                lands above the hydrogen edge (hydrogenic_twophoton_mod).
!
! Measured shares of the channel at n_e = 100: Ly-alpha 64-70%, two-photon
! 22-26%, free-bound 7.5-10% over 1e4-3e4 K.  The two-photon term alone --
! the only one of the three the written continuum carried -- is under a
! quarter of it.
!
! The n = 2 exit branching f_2s = alpha_eff(2^2S)/alpha_B comes from
! data/atomic/heii_caseB_n2.txt, both parts from Storey & Hummer (1995) so
! that the ratio closes one cascade calculation (see the file header and
! tools/fitting/make_heii_cascade.py).  The RATE that branching applies to is
! MoCHII's own alphaB_HeIII(T), the one the ionization balance uses, exactly
! as channel 4 applies its fractions to alphaB_HeII.
!
! Approximations, with their range:
!  - He II Ly-alpha is emitted at 40.813 eV and transported as a continuum
!    photon.  In reality it resonantly scatters off He II many times before
!    being absorbed; since it cannot photoionize He II (40.8 < 54.4 eV) the
!    scattering only diffuses it spatially before hydrogen or He I takes it,
!    and the He III zone is compact wherever it exists.  This is what
!    MOCASSIN does with its HeIILyman lines.
!  - Only Ly-alpha is emitted, not the higher He II Lyman lines.  Case B
!    routes every cascade through n = 2, so Ly-beta and above are converted
!    rather than emitted; MOCASSIN carries four of them because it works from
!    a Lyman-line emissivity set rather than from the n = 2 exit.
!  - Lines ending on n >= 2 are all below 13.604 eV (the series limit) and so
!    are correctly absent from the band.
!---------------------------------------------------------------------------
  use define
  use recomb_continuum_mod,   only : recomb_continuum_setup, gamma_tab, &
                                     gamma_ff, NU_HEII
  use hydrogenic_twophoton_mod, only : h2q_band_setup, h2q_band_sample_energy, &
                                     h2q_band_energy_per_decay, &
                                     h2q_band_photons_per_decay
  implicit none
  private

  public :: heii_cascade_setup, heii_cascade_is_ready
  public :: heii_branch_luminosity, heii_fb_sample_energy
  public :: heii_2q_sample_energy, heii_alpha_2s, heii_f_2s
  public :: E_HEII_LYA

  !--- He II 2 -> 1 transition energy [eV]; the two-photon pair shares it
  real(kind=wp), parameter :: E_HEII_LYA = 40.813_wp
  real(kind=wp), parameter :: RYD_HZ     = 3.2898419602508e15_wp   ! Hz/Ryd
  real(kind=wp), parameter :: RYD_EV     = 13.605693_wp            ! eV/Ryd

  !--- free-bound integration/sampling grid
  integer, parameter :: NT_FB  = 48        ! log T, 1e3 - 1e5 K
  integer, parameter :: NNU_FB = 512       ! nu from the H edge to the He II edge
  real(kind=wp), save :: fb_T(NT_FB) = 0.0_wp, fb_nu(NNU_FB) = 0.0_wp
  real(kind=wp), save :: fb_lum(NT_FB) = 0.0_wp        ! erg cm^3 s^-1
  real(kind=wp), save :: fb_cdf(NNU_FB, NT_FB) = 1.0_wp

  !--- SH95 branching table
  integer, parameter :: MAX_TAB = 64
  integer,       save :: n_aB = 0, n_a2s = 0
  real(kind=wp), save :: tab_aB_T(MAX_TAB), tab_aB(MAX_TAB)
  real(kind=wp), save :: tab_a2s_T(MAX_TAB), tab_a2s(MAX_TAB)

  logical, save :: ready = .false.

contains

  !=========================================================================
  logical function heii_cascade_is_ready() result(ok)
    ok = ready
  end function heii_cascade_is_ready

  !=========================================================================
  ! Log-log interpolation with endpoint clamping -- the rule the data file
  ! header states, so file and code agree by construction.
  !=========================================================================
  real(kind=wp) function loglog(T, xt, xv, n) result(y)
    real(kind=wp), intent(in) :: T, xt(:), xv(:)
    integer,       intent(in) :: n
    real(kind=wp) :: lt, w
    integer :: i
    if (T <= xt(1)) then
       y = xv(1);  return
    else if (T >= xt(n)) then
       y = xv(n);  return
    end if
    do i = 1, n-1
       if (T < xt(i+1)) exit
    end do
    lt = log10(T/xt(i));  w = lt/log10(xt(i+1)/xt(i))
    y = 10.0_wp**(log10(xv(i))*(1.0_wp - w) + log10(xv(i+1))*w)
  end function loglog

  !--- SH95 alpha_eff(2^2S) [cm^3 s^-1], n_e = 100.  Shared with the written
  !--- He II two-photon continuum so that the two cannot disagree.
  real(kind=wp) function heii_alpha_2s(T) result(a)
    real(kind=wp), intent(in) :: T
    a = loglog(T, tab_a2s_T, tab_a2s, n_a2s)
  end function heii_alpha_2s

  !--- Fraction of case-B recombinations leaving n = 2 through 2^2S.
  real(kind=wp) function heii_f_2s(T) result(f)
    real(kind=wp), intent(in) :: T
    f = heii_alpha_2s(T)/loglog(T, tab_aB_T, tab_aB, n_aB)
    f = min(max(f, 0.0_wp), 1.0_wp)
  end function heii_f_2s

  !=========================================================================
  subroutine heii_cascade_setup()
    use utility, only : fatal_error
    implicit none
    character(len=32) :: key
    real(kind=wp) :: nu_lo, nu_hi, T, g, seg, e_hi
    integer :: unit, ios, i, it

    call recomb_continuum_setup()

    !--- n = 2 exit branching
    open(newunit=unit, file=trim(par%atomic_dir)//'/heii_caseB_n2.txt', &
         status='old', iostat=ios)
    if (ios /= 0) call fatal_error('cannot open ' // &
       trim(par%atomic_dir)//'/heii_caseB_n2.txt (par%atomic_dir)')
    do
       read(unit,'(a)',iostat=ios) key
       if (ios /= 0) exit
       if (key(1:1) == '#' .or. len_trim(key) == 0) cycle
       if (index(key, 'ALPHA_B') == 1) then
          read(key(8:),*) n_aB
          if (n_aB > MAX_TAB) call fatal_error('heii_caseB_n2.txt ALPHA_B too long')
          do i = 1, n_aB
             read(unit,*) tab_aB_T(i), tab_aB(i)
          end do
       else if (index(key, 'ALPHA_2S') == 1) then
          read(key(9:),*) n_a2s
          if (n_a2s > MAX_TAB) call fatal_error('heii_caseB_n2.txt ALPHA_2S too long')
          do i = 1, n_a2s
             read(unit,*) tab_a2s_T(i), tab_a2s(i)
          end do
       end if
    end do
    close(unit)
    if (n_aB < 2 .or. n_a2s < 2) &
       call fatal_error('heii_caseB_n2.txt is missing a table')

    !--- Ly-alpha is a line at 40.813 eV: a band that does not reach it cannot
    !--- carry the largest of the three paths, and clamping would fabricate an
    !--- energy.  Refuse instead.  In practice this cannot bind on a nebula
    !--- that has any He III, which needs 54.416 eV photons to make.
    if (par%heii_diffuse .and. par%eion_max < E_HEII_LYA) &
       call fatal_error('par%heii_diffuse needs par%eion_max >= 40.813 eV ' // &
          '(the He II Ly-alpha energy)')

    !--- two-photon band: the pair spans E_HEII_LYA, the window runs from the
    !--- hydrogen edge to the band ceiling
    e_hi = min(E_HEII_LYA, par%eion_max)
    call h2q_band_setup(E_HEII_LYA, eth_HI, e_hi)

    !--- free-bound: (fb+ff table) - ff, from the hydrogen edge to the He II
    !--- ground edge, where channel 3 takes over.  Below the He II n = 2 edge
    !--- at 13.604 eV the sliver down to 13.598 eV is the n >= 3 tail, which
    !--- is ionizing too and costs nothing to keep.
    !--- stop at eth_HeII, the threshold channel 3 starts from, rather than at
    !--- the fb_ff edge constant NU_HEII (4.001067 Ryd = 54.437 eV): the two
    !--- differ by 0.02 eV and the boundary between the two channels must be
    !--- one number.  The table ends at 3.9996 Ryd anyway, and what lies
    !--- between is suppressed by exp(-40 eV/kT).
    nu_lo = eth_HI/RYD_EV
    nu_hi = min(eth_HeII, par%eion_max)/RYD_EV
    do it = 1, NT_FB
       fb_T(it) = 10.0_wp**(3.0_wp + 2.0_wp*real(it-1,wp)/real(NT_FB-1,wp))
    end do
    fb_lum = 0.0_wp;  fb_cdf = 0.0_wp
    do i = 1, NNU_FB
       fb_nu(i) = nu_lo + (nu_hi - nu_lo)*real(i-1,wp)/real(NNU_FB-1,wp)
    end do
    !--- Simpson panel by panel, not the trapezoid: gamma_tab interpolates
    !--- log-linearly in nu between the two table points that bracket this
    !--- band, so the integrand is an exponential of scale kT, and at the cold
    !--- end of the T grid one panel spans nearly a full e-fold.  Simpson puts
    !--- that error at 3e-4 where the trapezoid would leave 7%.
    if (nu_hi > nu_lo) then
       do it = 1, NT_FB
          T = fb_T(it)
          seg = 0.0_wp
          do i = 1, NNU_FB-1
             g = (fb_nu(i+1) - fb_nu(i))/6.0_wp &
                 *(fb_emis(fb_nu(i), T) &
                   + 4.0_wp*fb_emis(0.5_wp*(fb_nu(i) + fb_nu(i+1)), T) &
                   + fb_emis(fb_nu(i+1), T))
             seg = seg + g
             fb_cdf(i+1, it) = seg
          end do
          fb_lum(it) = seg*1.0e-40_wp*RYD_HZ
          if (seg > 0.0_wp) then
             fb_cdf(:, it) = fb_cdf(:, it)/seg
          else
             fb_cdf(:, it) = 1.0_wp
          end if
          fb_cdf(NNU_FB, it) = 1.0_wp
       end do
    else
       fb_cdf = 1.0_wp
    end if
    ready = .true.

    if (mpar%p_rank == 0) then
       write(*,'(a)') ' HE2CAS: He III case-B cascade channel built'
       write(*,'(a,f7.4,a,f7.4,a,f8.4,a)') &
          '         at 1e4 K: f(2^2S) = ', heii_f_2s(1.0e4_wp), &
          ', f(2^2P) = ', 1.0_wp - heii_f_2s(1.0e4_wp), &
          ', 2q in band = ', h2q_band_energy_per_decay(), ' eV/decay'
    end if
  end subroutine heii_cascade_setup

  !--- He II free-bound coefficient [1e-40 erg cm^3/s/Hz]: the tabulated
  !--- fb+ff of recomb_continuum_mod with its free-free term removed.
  real(kind=wp) function fb_emis(nu, T) result(g)
    real(kind=wp), intent(in) :: nu, T
    g = gamma_tab(3, nu, T) - gamma_ff(nu, T, 2)
    g = max(g, 0.0_wp)
  end function fb_emis

  !=========================================================================
  ! The three H-ionizing branch luminosities per (n_e n_HeIII) [erg cm^3/s].
  ! One routine for the emission table and for both launchers, so the branch
  ! CDF and the emitted luminosity cannot disagree.
  !=========================================================================
  subroutine heii_branch_luminosity(T, l_fb, l_lya, l_2q)
    use recomb_mod, only : alphaB_HeIII
    real(kind=wp), intent(in)  :: T
    real(kind=wp), intent(out) :: l_fb, l_lya, l_2q
    real(kind=wp) :: aB, f2s
    l_fb = 0.0_wp;  l_lya = 0.0_wp;  l_2q = 0.0_wp
    if (.not. ready) return
    aB  = alphaB_HeIII(T)
    f2s = heii_f_2s(T)
    if (fb_lum(1) > 0.0_wp) l_fb = loglog(T, fb_T, fb_lum, NT_FB)
    l_lya = aB*(1.0_wp - f2s)*E_HEII_LYA*ev2erg
    l_2q  = aB*f2s*h2q_band_energy_per_decay()*ev2erg
  end subroutine heii_branch_luminosity

  !--- Energy [eV] of one in-band two-photon packet.
  real(kind=wp) function heii_2q_sample_energy(u) result(E)
    real(kind=wp), intent(in) :: u
    E = h2q_band_sample_energy(u)
  end function heii_2q_sample_energy

  !=========================================================================
  ! Energy [eV] of one free-bound packet at temperature T from the uniform u.
  ! The CDF is built per tabulated temperature; the nearest one in log T is
  ! used, which over the 48-point grid is at most 2.5% in T.
  !=========================================================================
  real(kind=wp) function heii_fb_sample_energy(T, u) result(E)
    real(kind=wp), intent(in) :: T, u
    real(kind=wp) :: uu, nu
    integer :: it, lo, hi, mid
    it = nint(1.0_wp + (log10(min(max(T, fb_T(1)), fb_T(NT_FB))) - 3.0_wp) &
              *real(NT_FB-1,wp)/2.0_wp)
    it = min(max(it, 1), NT_FB)
    uu = min(max(u, 0.0_wp), 1.0_wp)
    lo = 1;  hi = NNU_FB
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       if (uu > fb_cdf(mid, it)) then
          lo = mid
       else
          hi = mid
       end if
    end do
    if (fb_cdf(hi, it) > fb_cdf(lo, it)) then
       nu = fb_nu(lo) + (fb_nu(hi) - fb_nu(lo)) &
            *(uu - fb_cdf(lo, it))/(fb_cdf(hi, it) - fb_cdf(lo, it))
    else
       nu = fb_nu(lo)
    end if
    E = nu*RYD_EV
  end function heii_fb_sample_energy

end module heii_cascade_mod
