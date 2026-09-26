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
!   Lyman lines  np -> 1^2S: in case B Ly-alpha alone (2^2P, 40.813 eV in
!                one photon); in case A, when the Lyman lines escape
!                (par%lyman_lines = 'escaping'), every np -> 1s up to the
!                series limit, 40.8-54.4 eV.
!   two-photon   2^2S -> 1^2S, a pair sharing 40.813 eV; 88% of that energy
!                lands above the hydrogen edge (hydrogenic_twophoton_mod).
!
! Measured shares of the channel in case B: Ly-alpha 64-70%, two-photon
! 22-26%, free-bound 7.5-10% over 1e4-3e4 K.  The two-photon term alone --
! the only one of the three the written continuum carried -- is under a
! quarter of it.
!
! The exit branching -- f_2s = alpha_eff(2^2S)/alpha_B and the Lyman-line
! shares -- comes from one zero-density nl-resolved cascade
! (hydrogenic_recombination_mod, Storey & Hummer 1991 hydrogenic data) in the
! case of the recombination lines (lyman_line_case), so the shares close one
! calculation and sum to 1.  Case A exists since 2026-09-26; the channel was
! case B whatever par%lyman_lines said.  The case-B f_2s used to be SH95's
! alpha_eff(2^2S) at n_e = 100 on four anchors,
! 5000-30000 K, held at the end values outside them, over the SH95
! alpha-tot: the two agree to 0.3% at 1e4 K, and the zero-density cascade
! lies 0.6% below at 5000 K, the collisional l-mixing SH95 include.  The
! RATE that branching applies to is MoCHII's own alphaB_HeIII(T), the one
! the ionization balance uses, exactly as channel 4 applies its fractions
! to alphaB_HeII.
!
! The same three paths, absorbed where they are made instead of transported
! (par%heii_on_the_spot, heii_recombination_photons_mod), are counted by
! heii_cascade_photons from the SAME branch luminosities, split into a part
! below the He0 edge (24.587 eV, which only H0 can absorb) and HEII_NGROUP
! groups above it (which H0 and He0 share by opacity).  The energy of each
! path closes on heii_branch_luminosity to rounding and its photon count on
! what the transported samplers imply to 2e-5 (gate tests/heii_on_the_spot),
! so the transported channel and the on-the-spot sources cannot count the
! cascade differently.
!
! Approximations, with their range:
!  - He II Ly-alpha is emitted at 40.813 eV and transported as a continuum
!    photon.  In reality it resonantly scatters off He II many times before
!    being absorbed; since it cannot photoionize He II (40.8 < 54.4 eV) the
!    scattering only diffuses it spatially before hydrogen or He I takes it,
!    and the He III zone is compact wherever it exists.  This is what
!    MOCASSIN does with its HeIILyman lines.  It holds only while the
!    continuum photon is itself absorbed near its origin: on a 1e5 K
!    blackbody in the HII40 geometry (results/heii_on_the_spot/) hydrogen
!    in the 2 pc He III zone gives tau = 0.01-0.04 per pc at 40.8 eV, so the
!    transported Ly-alpha leaves a zone that resonance scattering would hold
!    it in.
!  - In case B only Ly-alpha is emitted: every cascade is routed through
!    n = 2, so Ly-beta and above are converted rather than emitted.  In case A
!    the higher lines are emitted at their own energies, n = 2..30 one by one
!    and n >= 31 together at the n = 31 energy (within 0.06 eV).  The same
!    resonance-scattering caveat as for Ly-alpha applies to each.
!  - Lines ending on n >= 2 are all below 13.604 eV (the series limit) and so
!    are correctly absent from the band.
!---------------------------------------------------------------------------
  use define
  use recomb_continuum_mod,   only : recomb_continuum_setup, gamma_tab, &
                                     gamma_ff, gamma_tab_nodes, NU_HEII
  use hydrogenic_recombination_mod, only : hydrogenic_f_2s, heii_lyman_shares, &
                                          heii_lyman_energy, HEII_NLY, &
                                          lyman_line_case
  use hydrogenic_twophoton_mod, only : h2q_band_setup, h2q_band_sample_energy, &
                                     h2q_band_energy_per_decay, &
                                     h2q_band_photons_per_decay
  implicit none
  private

  public :: heii_cascade_setup, heii_cascade_is_ready, heii_cascade_transported
  public :: heii_branch_luminosity, heii_fb_sample_energy
  public :: heii_2q_sample_energy, heii_lyman_sample_energy, heii_f_2s
  public :: heii_cascade_photons, heii_cascade_photons_type
  public :: E_HEII_LYA, HEII_NGROUP

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

  logical, save :: ready = .false.

  !--- the cascade split at the He0 edge for on-the-spot absorption: part 0
  !--- lies below eth_HeI, parts 1..HEII_NGROUP are groups from eth_HeI to
  !--- the top of the path (part_bounds).  Energy fractions and photon-mean
  !--- energies [eV]: free-bound at each tabulated T, two-photon once.
  integer, parameter :: HEII_NGROUP = 8
  integer, parameter :: NPANEL_PART = 64   ! Simpson panels in each part
  real(kind=wp), save :: fb_efrac(0:HEII_NGROUP, NT_FB) = 0.0_wp
  real(kind=wp), save :: fb_emean(0:HEII_NGROUP, NT_FB) = 0.0_wp
  real(kind=wp), save :: q2_efrac(0:HEII_NGROUP) = 0.0_wp
  real(kind=wp), save :: q2_emean(0:HEII_NGROUP) = 0.0_wp

  !--- photons emitted per unit n_e n_HeIII [cm^3 s^-1] in each part of each
  !--- path (1 free-bound, 2 Lyman lines, 3 two-photon), and their mean
  !--- energy [eV]; the Lyman lines all lie above the He0 edge, line k in
  !--- part min(k, HEII_NGROUP)
  type :: heii_cascade_photons_type
     real(kind=wp) :: n(0:HEII_NGROUP, 3) = 0.0_wp
     real(kind=wp) :: e(0:HEII_NGROUP, 3) = 0.0_wp
  end type heii_cascade_photons_type

contains

  !=========================================================================
  logical function heii_cascade_is_ready() result(ok)
    ok = ready
  end function heii_cascade_is_ready

  !=========================================================================
  ! Whether the He III cascade is carried by channel 5 of the diffuse field:
  ! par%heii_diffuse with the transported continuum, or with the on-the-spot
  ! one, where it is then the only transported recombination radiation --
  ! the ground captures of H, He and He+ and the He I cascade stay on the
  ! spot, and the He II on-the-spot sources keep only the He2+ ground
  ! captures (heii_recombination_photons_mod::heii_ots_nodes).  That second
  ! form exists for a He III zone optically thin at 40.8 eV, where absorbing
  ! the cascade in place deposits its E - IP (up to ~27 eV a photon) where
  ! the photons in fact leave (thermal_mod's validity count).  With
  ! 'escaping' the switch has no effect.
  !=========================================================================
  logical function heii_cascade_transported() result(t)
    t = par%heii_diffuse .and. &
        (trim(par%recombination_continuum) == 'transported' .or. &
         trim(par%recombination_continuum) == 'on_the_spot')
  end function heii_cascade_transported

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

  !--- Share of the recombinations into n >= 2 leaving through 2^2S, in the
  !--- case of the Lyman lines; the written He II two-photon continuum takes
  !--- the same share (nebcont_mod).
  real(kind=wp) function heii_f_2s(T) result(f)
    real(kind=wp), intent(in) :: T
    f = hydrogenic_f_2s(2, T, lyman_line_case())
  end function heii_f_2s

  !--- The Lyman-line path per unit alpha_B: photons (sum of the shares) and
  !--- energy [eV] (sum of share x line energy), line by line in f, e.
  subroutine lyman_path(T, f, e)
    real(kind=wp), intent(in)  :: T
    real(kind=wp), intent(out) :: f(HEII_NLY), e(HEII_NLY)
    integer :: k
    call heii_lyman_shares(T, lyman_line_case(), f)
    do k = 1, HEII_NLY
       e(k) = heii_lyman_energy(k)
    end do
  end subroutine lyman_path

  !--- Energy [eV] of one Lyman-line packet at T from the uniform u: the line
  !--- is picked by its ENERGY share (equal-energy packets), so the photons
  !--- the packets represent are the shares themselves.
  real(kind=wp) function heii_lyman_sample_energy(T, u) result(E)
    real(kind=wp), intent(in) :: T, u
    real(kind=wp) :: fl(HEII_NLY), el(HEII_NLY), tot, acc
    integer :: k
    call lyman_path(T, fl, el)
    tot = sum(fl*el)
    acc = 0.0_wp
    E = el(1)
    do k = 1, HEII_NLY
       acc = acc + fl(k)*el(k)
       E = el(k)
       if (u*tot <= acc) exit
    end do
  end function heii_lyman_sample_energy

  !=========================================================================
  subroutine heii_cascade_setup()
    use utility, only : fatal_error
    implicit none
    character(len=32) :: key
    real(kind=wp) :: nu_lo, nu_hi, T, g, seg, e_hi, e_top, ph
    integer :: unit, ios, i, it

    call recomb_continuum_setup()

    !--- The top of the emitted spectrum.  Transported (channel 5), a photon
    !--- above the band ceiling par%eion_max cannot be followed, so the
    !--- tables stop there.  Absorbed on the spot (par%heii_on_the_spot) or
    !--- only written out (nebcont_mod), the photons
    !--- are made by the gas and the band plays no part: the tables then run
    !--- to the physical ends, 40.813 eV for the pair and the He II edge for
    !--- the free-bound.
    if (heii_cascade_transported()) then
       e_top = par%eion_max
    else
       e_top = eth_HeII
    end if

    !--- Ly-alpha is a line at 40.813 eV: a band that does not reach it cannot
    !--- carry the largest of the three paths, and clamping would fabricate an
    !--- energy.  Refuse instead.  In practice this cannot bind on a nebula
    !--- that has any He III, which needs 54.416 eV photons to make.  Only
    !--- the transported channel needs the band (heii_cascade_transported).
    if (heii_cascade_transported() .and. par%eion_max < E_HEII_LYA) &
       call fatal_error('par%heii_diffuse needs par%eion_max >= 40.813 eV ' // &
          '(the He II Ly-alpha energy)')
    !--- in case A the higher Lyman lines reach the series limit likewise
    if (heii_cascade_transported() .and. lyman_line_case() == 'A' &
        .and. par%eion_max < heii_lyman_energy(HEII_NLY)) &
       call fatal_error('par%heii_diffuse with par%lyman_lines = ''escaping'' ' // &
          'needs par%eion_max >= 54.36 eV (the He II Lyman series)')

    !--- two-photon band: the pair spans E_HEII_LYA, the window runs from the
    !--- hydrogen edge to the top of the spectrum
    e_hi = min(E_HEII_LYA, e_top)
    call h2q_band_setup(E_HEII_LYA, eth_HI, e_hi)
    call twophoton_parts(e_hi)

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
    nu_hi = min(eth_HeII, e_top)/RYD_EV
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
    !--- that error at 3e-4 where the trapezoid would leave 7%.  A panel is
    !--- broken at the table node it straddles (fb_simpson): the He II n = 2
    !--- edge sits at 0.9999 Ryd, inside the first cell above the hydrogen
    !--- edge, and a Simpson panel across that jump undercounted the
    !--- free-bound by 0.6-0.8% at 1e4 K (2026-09-25).
    if (nu_hi > nu_lo) then
       do it = 1, NT_FB
          T = fb_T(it)
          seg = 0.0_wp
          do i = 1, NNU_FB-1
             call fb_simpson(fb_nu(i), fb_nu(i+1), T, 1, g, ph)
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
    do it = 1, NT_FB
       call free_bound_parts(nu_lo, nu_hi, fb_T(it), fb_efrac(:, it), &
                             fb_emean(:, it))
    end do
    ready = .true.

    if (mpar%p_rank == 0) then
       write(*,'(3a)') ' HE2CAS: He III cascade channel built, case ', &
          lyman_line_case(), ' (par%lyman_lines)'
       write(*,'(a,f7.4,a,f7.4,a,f8.4,a)') &
          '         at 1e4 K: f(2^2S) = ', heii_f_2s(1.0e4_wp), &
          ', f(Lyman) = ', 1.0_wp - heii_f_2s(1.0e4_wp), &
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
  subroutine heii_branch_luminosity(T, l_fb, l_lyman, l_2q)
    use recomb_mod, only : alphaB_HeIII
    real(kind=wp), intent(in)  :: T
    real(kind=wp), intent(out) :: l_fb, l_lyman, l_2q
    real(kind=wp) :: aB, f2s, f(HEII_NLY), e(HEII_NLY)
    l_fb = 0.0_wp;  l_lyman = 0.0_wp;  l_2q = 0.0_wp
    if (.not. ready) return
    aB  = alphaB_HeIII(T)
    f2s = heii_f_2s(T)
    call lyman_path(T, f, e)
    if (fb_lum(1) > 0.0_wp) l_fb = loglog(T, fb_T, fb_lum, NT_FB)
    l_lyman = aB*sum(f*e)*ev2erg
    l_2q    = aB*f2s*h2q_band_energy_per_decay()*ev2erg
  end subroutine heii_branch_luminosity

  !=========================================================================
  ! The part boundaries [eV] of a path spanning [eth_HI, e_top]: part 0 from
  ! eth_HI to the He0 edge, parts 1..HEII_NGROUP from the edge to e_top with
  ! boundaries quadratic in the group index, so the first group is
  ! 1/HEII_NGROUP^2 of the span (0.47 eV for the free-bound): the free-bound
  ! falls as exp(-E/kT), 0.9 eV at 1e4 K, and the H0 share of a group is
  ! evaluated at one energy.  A path that stops below the edge is all
  ! part 0.
  !=========================================================================
  subroutine part_bounds(e_top, ea, eb)
    real(kind=wp), intent(in)  :: e_top
    real(kind=wp), intent(out) :: ea(0:HEII_NGROUP), eb(0:HEII_NGROUP)
    integer :: k
    ea = e_top;  eb = e_top
    ea(0) = eth_HI;  eb(0) = min(eth_HeI, e_top)
    if (e_top <= eth_HeI) return
    do k = 1, HEII_NGROUP
       ea(k) = eth_HeI + (e_top - eth_HeI)*(real(k-1,wp)/real(HEII_NGROUP,wp))**2
       eb(k) = eth_HeI + (e_top - eth_HeI)*(real(k,  wp)/real(HEII_NGROUP,wp))**2
    end do
  end subroutine part_bounds

  !=========================================================================
  ! Simpson on fb_emis over [a, b] [Ryd] in n panels, broken at the table's
  ! frequency nodes, across which gamma_tab jumps: en = Int fb_emis dnu,
  ! ph = Int fb_emis/nu dnu (energy and photons, up to constants).
  !=========================================================================
  subroutine fb_simpson(a, b, T, n, en, ph)
    real(kind=wp), intent(in)  :: a, b, T
    integer,       intent(in)  :: n
    real(kind=wp), intent(out) :: en, ph
    real(kind=wp), allocatable :: nodes(:)
    real(kind=wp) :: lo, hi, h, x(3), g(3)
    integer :: m, i, j
    en = 0.0_wp;  ph = 0.0_wp
    if (b <= a) return
    nodes = gamma_tab_nodes(3)
    lo = a
    do
       !--- the piece ends at the first node strictly inside (lo, b)
       hi = b
       do j = 1, size(nodes)
          if (nodes(j) > lo .and. nodes(j) < hi) hi = nodes(j)
       end do
       m = max(1, nint(real(n,wp)*(hi - lo)/(b - a)))
       h = (hi - lo)/real(m, wp)
       do i = 1, m
          x = [lo + h*real(i-1,wp), lo + h*(real(i,wp) - 0.5_wp), lo + h*real(i,wp)]
          !--- the left end of a piece is a node: take the limit from above
          if (i == 1) x(1) = lo*(1.0_wp + 1.0e-12_wp)
          if (i == m) x(3) = hi*(1.0_wp - 1.0e-12_wp)
          do j = 1, 3
             g(j) = fb_emis(x(j), T)
          end do
          en = en + h/6.0_wp*(g(1) + 4.0_wp*g(2) + g(3))
          ph = ph + h/6.0_wp*(g(1)/x(1) + 4.0_wp*g(2)/x(2) + g(3)/x(3))
       end do
       if (hi >= b) exit
       lo = hi
    end do
  end subroutine fb_simpson

  !--- Energy fraction and photon-mean energy [eV] of each part of the
  !--- free-bound spectrum at T (fb_simpson over each part).
  subroutine free_bound_parts(nu_lo, nu_hi, T, efrac, emean)
    real(kind=wp), intent(in)  :: nu_lo, nu_hi, T
    real(kind=wp), intent(out) :: efrac(0:HEII_NGROUP), emean(0:HEII_NGROUP)
    real(kind=wp) :: ea(0:HEII_NGROUP), eb(0:HEII_NGROUP), a, b
    real(kind=wp) :: en(0:HEII_NGROUP), ph(0:HEII_NGROUP)
    integer :: k
    call part_bounds(nu_hi*RYD_EV, ea, eb)
    en = 0.0_wp;  ph = 0.0_wp
    do k = 0, HEII_NGROUP
       a = max(ea(k)/RYD_EV, nu_lo);  b = eb(k)/RYD_EV
       emean(k) = 0.5_wp*(ea(k) + eb(k))
       if (b <= a) cycle
       call fb_simpson(a, b, T, NPANEL_PART, en(k), ph(k))
       if (ph(k) > 0.0_wp) emean(k) = en(k)/ph(k)*RYD_EV
    end do
    if (sum(en) > 0.0_wp) then
       efrac = en/sum(en)
    else
       efrac = 0.0_wp;  efrac(0) = 1.0_wp
    end if
  end subroutine free_bound_parts

  !--- The same for the two-photon pair over its window [eth_HI, e_hi], on
  !--- the Nussbaumer & Schmutz shape hydrogenic_twophoton_mod normalizes.
  subroutine twophoton_parts(e_hi)
    use hydrogenic_twophoton_mod, only : h2q_photons_per_y
    real(kind=wp), intent(in) :: e_hi
    real(kind=wp) :: ea(0:HEII_NGROUP), eb(0:HEII_NGROUP), a, b, h, y(3), p(3)
    real(kind=wp) :: en(0:HEII_NGROUP), ph(0:HEII_NGROUP)
    integer :: k, i, j
    call part_bounds(e_hi, ea, eb)
    en = 0.0_wp;  ph = 0.0_wp
    do k = 0, HEII_NGROUP
       a = ea(k)/E_HEII_LYA;  b = eb(k)/E_HEII_LYA
       q2_emean(k) = 0.5_wp*(ea(k) + eb(k))
       if (b <= a) cycle
       h = (b - a)/real(NPANEL_PART, wp)
       do i = 1, NPANEL_PART
          y = [a + h*real(i-1,wp), a + h*(real(i,wp) - 0.5_wp), a + h*real(i,wp)]
          do j = 1, 3
             p(j) = h2q_photons_per_y(y(j))
          end do
          ph(k) = ph(k) + h/6.0_wp*(p(1) + 4.0_wp*p(2) + p(3))
          en(k) = en(k) + h/6.0_wp*(p(1)*y(1) + 4.0_wp*p(2)*y(2) + p(3)*y(3)) &
                  *E_HEII_LYA
       end do
       if (ph(k) > 0.0_wp) q2_emean(k) = en(k)/ph(k)
    end do
    if (sum(en) > 0.0_wp) then
       q2_efrac = en/sum(en)
    else
       q2_efrac = 0.0_wp;  q2_efrac(0) = 1.0_wp
    end if
  end subroutine twophoton_parts

  !=========================================================================
  ! Photons emitted per unit n_e n_HeIII [cm^3 s^-1] by the three paths at T,
  ! part by part, with their mean energies.  The energy of each path is its
  ! heii_branch_luminosity term, distributed over the parts by the energy
  ! fractions (free-bound: interpolated linearly in log T between the
  ! tabulated temperatures and renormalized; the transported sampler takes
  ! the nearest one); each part's photons are its energy over its mean
  ! energy.  So sum_k n e = l/ev2erg path by path, to rounding.
  !=========================================================================
  type(heii_cascade_photons_type) function heii_cascade_photons(T) result(c)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: l_fb, l_lyman, l_2q, lt, w, ef(0:HEII_NGROUP), em(0:HEII_NGROUP)
    real(kind=wp) :: f(HEII_NLY), e(HEII_NLY), aB
    integer :: it, k, kp
    if (.not. ready) return
    call heii_branch_luminosity(T, l_fb, l_lyman, l_2q)
    !--- free-bound
    if (T <= fb_T(1)) then
       ef = fb_efrac(:, 1);  em = fb_emean(:, 1)
    else if (T >= fb_T(NT_FB)) then
       ef = fb_efrac(:, NT_FB);  em = fb_emean(:, NT_FB)
    else
       lt = 1.0_wp + (log10(T) - 3.0_wp)*real(NT_FB-1,wp)/2.0_wp
       it = min(max(int(lt), 1), NT_FB-1)
       w  = lt - real(it, wp)
       ef = (1.0_wp - w)*fb_efrac(:, it) + w*fb_efrac(:, it+1)
       em = (1.0_wp - w)*fb_emean(:, it) + w*fb_emean(:, it+1)
       if (sum(ef) > 0.0_wp) ef = ef/sum(ef)
    end if
    c%e(:, 1) = em
    c%n(:, 1) = l_fb/ev2erg*ef/em
    !--- Lyman lines, all above the He0 edge: line k in part min(k, NGROUP),
    !--- photons sum(f) aB, energy the l_lyman term to rounding
    call lyman_path(T, f, e)
    aB = 0.0_wp
    if (sum(f*e) > 0.0_wp) aB = l_lyman/(ev2erg*sum(f*e))
    c%e(:, 2) = E_HEII_LYA
    do k = 1, HEII_NLY
       kp = min(k, HEII_NGROUP)
       c%n(kp, 2) = c%n(kp, 2) + aB*f(k)
       c%e(kp, 2) = 0.0_wp
    end do
    do k = 1, HEII_NLY
       kp = min(k, HEII_NGROUP)
       if (c%n(kp, 2) > 0.0_wp) c%e(kp, 2) = c%e(kp, 2) + aB*f(k)*e(k)/c%n(kp, 2)
    end do
    do kp = 0, HEII_NGROUP
       if (c%n(kp, 2) <= 0.0_wp) c%e(kp, 2) = E_HEII_LYA
    end do
    !--- two-photon
    c%e(:, 3) = q2_emean
    c%n(:, 3) = l_2q/ev2erg*q2_efrac/q2_emean
  end function heii_cascade_photons

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
