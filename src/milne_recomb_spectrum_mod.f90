module milne_recomb_spectrum_mod
!---------------------------------------------------------------------------
! MoCHII: ground-state recombination continua from the Milne relation.
!
! Detailed balance between photoionization and radiative recombination fixes
! the free-bound spectrum of each ground continuum.  Per unit photon energy,
! the photon-number emissivity of X^(i+1) + e -> X^i(ground) + h nu is
!
!     n(E) dE  proportional to  E^2 sigma_pi(E) exp[-(E - E_th)/kT] dE ,
!
! evaluated with the SAME ground-state photoionization cross section the
! transport uses, so a packet is absorbed with the cross section it was
! emitted with.  Restoring the Milne prefactor turns the same integral into
! the ground-state recombination coefficient,
!
!     alpha_1(T) = (g_n / 2 g_+) (2/sqrt(pi)) (kT)^(-3/2) sqrt(2/m_e)
!                  / (m_e c^2)  Int E^2 sigma_pi(E) exp[-(E-E_th)/kT] dE ,
!
! so the spectrum and the emission rate are not independent assumptions.
! milne_setup evaluates that integral and reports it against the
! alpha_A - alpha_B that sets the packet count.  That alpha_1 is itself a fit
! to this same Milne relation on these same cross sections
! (tools/fitting/fit_alpha1_milne.py, recomb_mod), so the two agree to the
! fit accuracy (0.01%); the comparison is reported, never rescaled into
! agreement, so that a drift between the fit and the cross sections shows.
!
! SAMPLING WEIGHT.  MoCHII emits equal-energy packets and fixes their number
! from the channel's ENERGY luminosity.  The packet energies must therefore
! follow the energy-luminosity density E n(E), not the photon-number density
! n(E).  With that pairing the represented photon rate is exactly the
! physical one: for a channel emitting R recombinations per second,
!
!     L_ch = R <E>_n ,   N_packet = L_ch / L_packet ,
!     photons represented = N_packet L_packet <1/E>_{En}
!                         = R <E>_n (1/<E>_n) = R .
!
! Drawing from n(E) while counting packets from the energy luminosity would
! over-represent the photon rate by <E>_n <1/E>_n, which exceeds one for any
! spread distribution: 1.003 for H I at 10^4 K and 1.010 at 2x10^4 K.  The
! historical 'exponential' model carries exactly that bias and keeps it,
! because its purpose is to reproduce earlier runs.
!
! TEMPERATURE INTERPOLATION.  Tables are indexed by x = (E - E_th)/kT, so the
! shape varies only weakly with temperature.  The CDF is interpolated in
! log T first and inverted once.  Inverting two neighboring temperature rows
! and interpolating the resulting energies (as some codes do) does not sample
! any single distribution; taking the nearest row loses the interpolation
! altogether.
!
! BAND TRUNCATION.  The tables stop at par%eion_max, and the energy radiated
! per recombination that milne_energy_per_recomb returns is the IN-BAND part.
! Photons above the band are therefore dropped consistently from both the
! luminosity and the sampler instead of being clipped onto the last bin.  The
! dropped fraction is reported at setup.
!
! The tables are filled once and are read only afterwards, so packet launch
! stays deterministic and free of shared mutable state.
!
! METAL GROUND CONTINUA (diffuse_mod channel 6; since 2026-10-11).  Every
! metal transition X^(i+1) + e -> X^i(ground term) + h nu emits the same
! Milne spectrum on the outer-shell part of its transport cross section
! (species_mod::species_outer_shell_cross_section, where the choice and its
! validity are written), with the factor g_n/(2 g_+) of the two ground
! terms (data/atomic/ground_term_weights.txt, NIST ground terms checked
! against CHIANTI).  No separate fit sets their rate: the Milne integral
! itself is alpha_1 and the energy, so the emission rate, the spectrum
! and the absorption are one set.  Per n_e n(X^(i+1)) a transition radiates
! inside the band [E_lo, eion_max], E_lo = max(E_th, eion_min),
!     P(T) = Int_band E d alpha_1(E)  [erg cm^3 s^-1]
!          = Q(T) exp[-(E_lo - E_th)/kT],
! Q tabulated as ln Q linear in log T (the edge Boltzmann factor, which a
! neutral with E_th < eion_min carries, is applied at the leaf).  Over all
! E this is n_e n_+ [alpha_1 IP + beta_1], beta_1 = kT [(3/2) alpha_1 +
! d alpha_1/d ln T] the kinetic energy those captures take from the gas --
! which species_mod::metal_ionization_cooling already removes as part of
! the beta of the total RR + DR rate -- and IP alpha_1 the potential energy
! the photoionization that made the ion had spent; nothing is counted twice
! (gate tests/metal_recombination_continuum).  The part below eion_min (the
! neutrals' continua in the FUV band) and above eion_max is not
! transported.
! CAP.  The ground captures are a part of the recombinations the
! ionization balance uses (species_mod alpha_rec, Badnell RR + DR), so their
! rate is held to alpha_1 <= alpha_rec: where the Milne integral exceeds
! the total (Cl I at nebular temperatures, 3.3x at 1e4 K: VY95's
! Hartree-Dirac-Slater outer shell of neutral Cl against a Badnell RR with
! no DR) the channel is scaled by alpha_rec/alpha_1 (mcap, printed at
! setup).  Otherwise the gas would radiate more potential energy than its
! recombinations release.  S I and Si I exceed the radiative rate alone
! (2.2x, 1.2x at 1e4 K) but not RR + DR: VFKY96 averages the near-threshold
! autoionizing resonances, whose inverse is dielectronic capture.
! TABULATED CROSS SECTIONS (Ne III, since 2026-10-11).  Where the element
! file carries a PHTABLE, the outer-shell cross section is the tabulated
! R-matrix ground-term cross section below its last node (species_mod),
! with its autoionizing resonances resolved (1.4 meV nodes near threshold).
! A fixed x step would alias them, so for such a transition every integral
! of this module -- ln Q, alpha_1 and the moments -- steps on the table's
! nodes (Simpson on each piece between nodes, where the cross section is
! linear, pieces no longer than DX_PIECE).  Its packets are NOT drawn from
! the x-indexed CDF: the resonances sit at fixed E, so in x they move with
! T and a CDF interpolated between two T rows puts their weight into the
! wrong x intervals (measured: 30x too many packets 0.1 eV beside a
! resonance at 1e4 K).  Instead the band is cut into energy blocks fixed in
! E (BLOCK_NODES table nodes each, BLOCK_DE beyond the table), each block's
! energy integral is tabulated on the T rows and interpolated linearly in
! 1/T in its logarithm (exact for the Boltzmann factor, which carries the
! T dependence of a block), a block is chosen from these, and the energy
! inside it is placed by inverting the node-stepped integral at the leaf
! temperature itself.  The emitted spectrum then follows the resonances
! node by node, and emission, alpha_1 and absorption use one function.
! NOT COVERED: captures into excited levels of X^i emit E = IP - E_exc +
! KE, ionizing H when IP - E_exc > 13.598 eV (O II, O III ... into their
! low excited terms), and the cascades and dielectronic stabilizations that
! follow emit EUV lines (e.g. O III 374 A); none of them is in the diffuse
! field.  The gate bounds them by n_e n_+ (alpha_RR - alpha_1) (IP + kT).
!---------------------------------------------------------------------------
  use define
  use photo_xsec, only : sigma_HI, sigma_HeI, sigma_HeII
  implicit none
  private

  public :: milne_setup, milne_sample_energy, milne_energy_per_recomb
  public :: milne_is_ready
  public :: MILNE_HI, MILNE_HEI, MILNE_HEII
  public :: milne_metal_count, milne_metal_band_coefficient
  public :: milne_metal_sample_energy, milne_metal_alpha1
  public :: milne_metal_energy_integral, milne_metal_excited_share
  public :: milne_metal_ground_factor, milne_metal_capture_cap

  !--- ground continua, in the order diffuse_mod uses for its channels
  integer, parameter :: MILNE_HI   = 1        ! H II   + e -> H I  (1s)
  integer, parameter :: MILNE_HEI  = 2        ! He II  + e -> He I (1^1S)
  integer, parameter :: MILNE_HEII = 3        ! He III + e -> He II(1s)
  integer, parameter :: NCH = 3

  !--- table shape.  NX nodes are clustered towards the threshold, where the
  !--- emissivity peaks; XMAX = 60 puts the truncation at exp(-60) of the
  !--- near-threshold value, far below any other error in the channel.
  integer,       parameter :: NT   = 192
  integer,       parameter :: NX   = 512
  real(kind=wp), parameter :: TLO  = 1.0e2_wp, THI = 1.0e6_wp
  real(kind=wp), parameter :: XMAX = 60.0_wp

  real(kind=wp), save :: lgT(NT)                  ! log10 of the temperature nodes
  real(kind=wp), save :: xg(NX)                   ! x = (E - E_th)/kT nodes
  real(kind=wp), save :: cdf_energy(NX, NT, NCH)  ! energy-luminosity CDF (sampling)
  real(kind=wp), save :: xbar_number(NT, NCH)     ! <x> under the photon-number PDF
  real(kind=wp), save :: frac_in_band(NT, NCH)    ! photon fraction below eion_max
  real(kind=wp), save :: eth_ch(NCH)
  logical,       save :: ready = .false.

  !--- metal ground continua, one per transition in species_transition order
  integer,       save :: nmetal = 0
  real(kind=wp), save, allocatable :: mcdf(:,:,:)     ! (NX, NT, nmetal) energy CDF in x
  real(kind=wp), save, allocatable :: mlnq(:,:)       ! (NT, nmetal) ln Q [erg cm^3/s]
  real(kind=wp), save, allocatable :: meth(:), melo(:), mgfac(:), me1(:)
  real(kind=wp), save, allocatable :: mbreak(:)        ! einn [eV], 0 = none
  real(kind=wp), save, allocatable :: mcap(:,:)       ! (NT, nmetal) capture cap
  integer,       save, allocatable :: mie(:), mit(:)
  !--- transitions whose outer-shell cross section is a PHTABLE (species_mod)
  logical,       save, allocatable :: mtab(:)
  !--- longest piece [in x = (E - E0)/kT] of a node-stepped quadrature
  real(kind=wp), parameter :: DX_PIECE = 0.02_wp
  !--- energy blocks of a tabulated transition (module header): table nodes
  !--- per block, and the block width beyond the table [eV]
  integer,       parameter :: BLOCK_NODES = 32
  real(kind=wp), parameter :: BLOCK_DE = 0.25_wp
  type energy_blocks_type
     integer :: nb = 0
     real(kind=wp), allocatable :: e(:)        ! (0:nb) block edges [eV]
     real(kind=wp), allocatable :: lni(:,:)    ! (nb, NT) ln Int E f dx
  end type energy_blocks_type
  type(energy_blocks_type), save, allocatable :: mblk(:)

contains

  !=========================================================================
  logical function milne_is_ready() result(ok)
    ok = ready
  end function milne_is_ready

  !=========================================================================
  ! Ground-state photoionization cross section of channel ich at E [eV].
  ! These are the transport cross sections, not a private copy.
  !=========================================================================
  real(kind=wp) function milne_sigma(ich, E) result(sig)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: E
    select case (ich)
    case (MILNE_HI);   sig = sigma_HI(E)
    case (MILNE_HEI);  sig = sigma_HeI(E)
    case default;      sig = sigma_HeII(E)
    end select
  end function milne_sigma

  !=========================================================================
  ! Photon-number emissivity shape at x, up to a temperature-only factor:
  !     n(x) = E^2 sigma_pi(E) exp(-x),   E = E_th + x kT.
  !=========================================================================
  real(kind=wp) function milne_shape(ich, x, kT_eV) result(f)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: x, kT_eV
    real(kind=wp) :: E
    E = eth_ch(ich) + x*kT_eV
    f = E*E*milne_sigma(ich, E)*exp(-x)
  end function milne_shape

  !=========================================================================
  ! Photon-number shape of metal channel k at x, measured from the band's
  ! lower edge E_lo (= E_th unless E_th < eion_min):
  !     E^2 sigma_outer(E) exp(-x),   E = E_lo + x kT.
  !=========================================================================
  real(kind=wp) function milne_metal_shape(k, x, kT_eV, e0) result(f)
    use species_mod, only : species_outer_shell_cross_section
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: x, kT_eV, e0
    real(kind=wp) :: E
    E = e0 + x*kT_eV
    f = E*E*species_outer_shell_cross_section(mie(k), mit(k), E)*exp(-x)
  end function milne_metal_shape

  !--- Milne prefactor (g_n/2g_+) (2/sqrt(pi)) (kT)^(-3/2) sqrt(2/m_e)/(m_e c^2)
  !--- [erg^-3 cm s^-1 per unit integral in erg^3 cm^2]
  real(kind=wp) function milne_prefactor(gfac, T_K) result(c)
    real(kind=wp), intent(in) :: gfac, T_K
    c = gfac*(2.0_wp/sqrt(pi))*(kboltz_cgs*T_K)**(-1.5_wp) &
        *sqrt(2.0_wp/me_cgs)/(me_cgs*clight_cgs**2)
  end function milne_prefactor

  !=========================================================================
  ! Build the tables and report the alpha_1 cross-check.  Call once, after
  ! the band setup has fixed par%eion_max.  With par%use_metals the metal
  ! ground continua are built too (milne_metal_setup).
  !=========================================================================
  subroutine milne_setup()
    use recomb_mod, only : alphaA_HII, alphaB_HII, alphaA_HeII, alphaB_HeII, &
                           alphaA_HeIII, alphaB_HeIII
    implicit none
    !--- g_n/(2 g_+): H I(1s) 2/(2*1), He I(1^1S) 1/(2*2), He II(1s) 2/(2*1)
    real(kind=wp), parameter :: gfac(NCH) = [1.0_wp, 0.25_wp, 1.0_wp]
    character(len=6), parameter :: chname(NCH) = &
       [character(len=6) :: 'H I', 'He I', 'He II']
    real(kind=wp), parameter :: TCHK(3) = [5.0e3_wp, 1.0e4_wp, 2.0e4_wp]
    real(kind=wp) :: t, T_K, kT_eV, x_hi, xa, xb, xm, fa, fb, fm
    real(kind=wp) :: seg_n, seg_e, sum_n, sum_e, sum_xn, tail_n
    real(kind=wp) :: a1_milne, a1_rate, worst, fmin
    integer :: ich, it, ix

    eth_ch = [eth_HI, eth_HeI, eth_HeII]

    do ix = 1, NX
       t = real(ix-1, wp)/real(NX-1, wp)
       xg(ix) = XMAX*t*t                    ! quadratic clustering near x = 0
    end do
    do it = 1, NT
       lgT(it) = log10(TLO) + (log10(THI) - log10(TLO)) &
                 *real(it-1, wp)/real(NT-1, wp)
    end do

    do ich = 1, NCH
       do it = 1, NT
          T_K   = 10.0_wp**lgT(it)
          kT_eV = kboltz_cgs*T_K/ev2erg
          !--- band truncation in x; <= 0 means the edge is out of band
          x_hi = min(XMAX, (par%eion_max - eth_ch(ich))/kT_eV)
          if (x_hi <= 0.0_wp) then
             cdf_energy(:, it, ich) = 1.0_wp
             xbar_number(it, ich)   = 0.0_wp
             frac_in_band(it, ich)  = 0.0_wp
             cycle
          end if
          !--- cumulative Simpson over the truncated range
          sum_n = 0.0_wp;  sum_e = 0.0_wp;  sum_xn = 0.0_wp
          cdf_energy(1, it, ich) = 0.0_wp
          do ix = 1, NX-1
             xa = xg(ix)
             xb = min(xg(ix+1), x_hi)
             if (xb > xa) then
                xm = 0.5_wp*(xa + xb)
                fa = milne_shape(ich, xa, kT_eV)
                fm = milne_shape(ich, xm, kT_eV)
                fb = milne_shape(ich, xb, kT_eV)
                seg_n = (xb - xa)/6.0_wp*(fa + 4.0_wp*fm + fb)
                seg_e = (xb - xa)/6.0_wp &
                        *( fa*(eth_ch(ich) + xa*kT_eV) &
                         + 4.0_wp*fm*(eth_ch(ich) + xm*kT_eV) &
                         + fb*(eth_ch(ich) + xb*kT_eV) )
                sum_n  = sum_n  + seg_n
                sum_e  = sum_e  + seg_e
                sum_xn = sum_xn + (xb - xa)/6.0_wp &
                                  *(fa*xa + 4.0_wp*fm*xm + fb*xb)
             end if
             cdf_energy(ix+1, it, ich) = sum_e
          end do
          !--- the untruncated photon integral, for the in-band fraction
          tail_n = sum_n
          if (x_hi < XMAX) then
             xa = x_hi
             do while (xa < XMAX)
                xb = min(xa + 0.25_wp, XMAX)
                xm = 0.5_wp*(xa + xb)
                tail_n = tail_n + (xb - xa)/6.0_wp &
                         *( milne_shape(ich, xa, kT_eV) &
                          + 4.0_wp*milne_shape(ich, xm, kT_eV) &
                          + milne_shape(ich, xb, kT_eV) )
                xa = xb
             end do
          end if
          if (sum_e > 0.0_wp) then
             cdf_energy(:, it, ich) = cdf_energy(:, it, ich)/sum_e
          else
             cdf_energy(:, it, ich) = 1.0_wp
          end if
          cdf_energy(NX, it, ich) = 1.0_wp
          xbar_number(it, ich)  = sum_xn/max(sum_n, tinest)
          frac_in_band(it, ich) = sum_n/max(tail_n, tinest)
       end do
    end do
    ready = .true.

    !--- alpha_1 from the same spectral integral against the rate data that
    !--- sets the packet count.  Independent paths: report, do not rescale.
    if (mpar%p_rank == 0) then
       write(*,'(a)') ' MILNE: ground-continuum tables built (photon PDF ' // &
          'E^2 sigma(E) exp[-(E-Eth)/kT]; packets drawn from E n(E))'
       write(*,'(a)') ' MILNE: alpha_1 from the spectral integral vs ' // &
          'alpha_A - alpha_B used for the packet count'
       worst = 0.0_wp
       do ich = 1, NCH
          do it = 1, 3
             T_K = TCHK(it)
             a1_milne = milne_alpha1(ich, T_K, gfac(ich))
             select case (ich)
             case (MILNE_HI);  a1_rate = alphaA_HII(T_K)   - alphaB_HII(T_K)
             case (MILNE_HEI); a1_rate = alphaA_HeII(T_K)  - alphaB_HeII(T_K)
             case default;     a1_rate = alphaA_HeIII(T_K) - alphaB_HeIII(T_K)
             end select
             write(*,'(a,a6,a,f7.0,2(a,es11.4),a,f8.4)') &
                '        ', chname(ich), '  T =', T_K, &
                ' K   Milne ', a1_milne, '   rates ', a1_rate, &
                '   ratio ', a1_milne/max(a1_rate, tinest)
             worst = max(worst, abs(a1_milne/max(a1_rate, tinest) - 1.0_wp))
          end do
       end do
       write(*,'(a,f6.2,a)') ' MILNE: largest alpha_1 disagreement ', &
          100.0_wp*worst, '% (atomic-data difference, not rescaled)'
       !--- Only the temperatures the solver can reach matter: the table
       !--- itself spans 1e2-1e6 K, and at the top of that range kT alone
       !--- exceeds the band, which is not a property of any H II region.
       do ich = 1, NCH
          fmin = 1.0_wp
          do it = 1, NT
             T_K = 10.0_wp**lgT(it)
             if (T_K >= par%te_min .and. T_K <= par%te_max) &
                fmin = min(fmin, frac_in_band(it, ich))
          end do
          if (fmin < 0.999_wp) &
             write(*,'(3a,es10.3,a,f8.0,a)') ' WARNING: MILNE ', &
                trim(chname(ich)), ' loses photons above par%eion_max; ' // &
                'smallest in-band fraction = ', fmin, ' within [te_min,', &
                par%te_max, ']'
       end do
    end if
    if (par%use_metals) call milne_metal_setup()
  end subroutine milne_setup

  !=========================================================================
  ! alpha_1(T) [cm^3/s] implied by the tabulated spectrum, with the Milne
  ! prefactor restored.  Used only for the setup cross-check.
  !=========================================================================
  real(kind=wp) function milne_alpha1(ich, T_K, gfac) result(a1)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: T_K, gfac
    real(kind=wp) :: kT_eV, kT_erg, xa, xb, xm, acc
    kT_eV  = kboltz_cgs*T_K/ev2erg
    kT_erg = kboltz_cgs*T_K
    acc = 0.0_wp
    xa  = 0.0_wp
    do while (xa < XMAX)
       xb = min(xa + 0.05_wp, XMAX)
       xm = 0.5_wp*(xa + xb)
       acc = acc + (xb - xa)/6.0_wp &
             *( milne_shape(ich, xa, kT_eV) &
              + 4.0_wp*milne_shape(ich, xm, kT_eV) &
              + milne_shape(ich, xb, kT_eV) )
       xa = xb
    end do
    !--- acc is Int E^2 sigma exp(-x) dx with E in eV; convert to erg^3 cm^2
    acc = acc*kT_eV*ev2erg**3
    a1  = gfac*(2.0_wp/sqrt(pi))*kT_erg**(-1.5_wp) &
          *sqrt(2.0_wp/me_cgs)/(me_cgs*clight_cgs**2)*acc
  end function milne_alpha1

  !=========================================================================
  ! Temperature bracket: node it and weight w for log10(T), clamped.
  !=========================================================================
  subroutine milne_tbracket(T_K, it, w)
    real(kind=wp), intent(in)  :: T_K
    integer,       intent(out) :: it
    real(kind=wp), intent(out) :: w
    real(kind=wp) :: s
    s = (log10(max(T_K, tinest)) - lgT(1))/(lgT(2) - lgT(1))
    s = min(max(s, 0.0_wp), real(NT-1, wp))
    it = min(int(s) + 1, NT-1)
    w  = s - real(it-1, wp)
  end subroutine milne_tbracket

  !=========================================================================
  ! In-band energy radiated per ground-state recombination [eV].  Multiplying
  ! by n_e n_ion alpha_1 V gives the channel's in-band energy luminosity, so
  ! the packet count and the sampler below describe the same photons.
  !=========================================================================
  real(kind=wp) function milne_energy_per_recomb(ich, T_K) result(eper)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: T_K
    real(kind=wp) :: w, kT_eV, xbar, fkeep
    integer :: it
    call milne_tbracket(T_K, it, w)
    xbar  = (1.0_wp - w)*xbar_number(it, ich)  + w*xbar_number(it+1, ich)
    fkeep = (1.0_wp - w)*frac_in_band(it, ich) + w*frac_in_band(it+1, ich)
    kT_eV = kboltz_cgs*T_K/ev2erg
    eper  = fkeep*(eth_ch(ich) + xbar*kT_eV)
  end function milne_energy_per_recomb

  !=========================================================================
  ! Photon energy [eV] for one equal-energy packet of channel ich at T_K,
  ! from the uniform variate u.  The CDF is interpolated in log T and then
  ! inverted once; both the pseudorandom and the Sobol launch call this, so
  ! they sample the same distribution by construction.
  !=========================================================================
  real(kind=wp) function milne_sample_energy(ich, T_K, u) result(E)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: T_K, u
    real(kind=wp) :: w, kT_eV, uu, ca, cb, cm, x
    integer :: it, lo, hi, mid
    call milne_tbracket(T_K, it, w)
    kT_eV = kboltz_cgs*T_K/ev2erg
    uu    = min(max(u, 0.0_wp), 1.0_wp)
    !--- binary search on the temperature-interpolated CDF; the interpolant
    !--- of two monotone CDFs is monotone, so the bracket is well defined.
    lo = 1;  hi = NX
    ca = (1.0_wp - w)*cdf_energy(lo, it, ich) + w*cdf_energy(lo, it+1, ich)
    cb = (1.0_wp - w)*cdf_energy(hi, it, ich) + w*cdf_energy(hi, it+1, ich)
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       cm  = (1.0_wp - w)*cdf_energy(mid, it, ich) &
           + w*cdf_energy(mid, it+1, ich)
       if (uu > cm) then
          lo = mid;  ca = cm
       else
          hi = mid;  cb = cm
       end if
    end do
    if (cb > ca) then
       x = xg(lo) + (xg(hi) - xg(lo))*(uu - ca)/(cb - ca)
    else
       x = xg(lo)
    end if
    E = eth_ch(ich) + x*kT_eV
    E = min(E, par%eion_max*(1.0_wp - 1.0e-12_wp))
  end function milne_sample_energy

  !=========================================================================
  ! Metal ground continua: tables of every transition of the registry
  ! (diffuse_mod channel 6).  See the module header.
  !=========================================================================
  subroutine milne_metal_setup()
    use species_mod, only : species_transition_count, species_transition, &
                            elem_name, elem_eth, species_radiative_recombination, &
                            species_recombination_coefficient, &
                            species_outer_shell_break, species_photo_table_end
    use utility,     only : fatal_error
    implicit none
    real(kind=wp), parameter :: TCHK = 1.0e4_wp
    real(kind=wp) :: t, T_K, kT_eV, x_hi, xa, xb, xm, fa, fb, fm, sum_e
    real(kind=wp) :: g_lo, g_hi, e1_hi, a1, arr, ratio, worst, share, xbr
    integer :: k, it, ix, ie, ist, kworst, ib
    character(len=512) :: big, capped

    if (allocated(mcdf)) return
    nmetal = species_transition_count()
    if (nmetal <= 0) return
    allocate(mcdf(NX, NT, nmetal), mlnq(NT, nmetal), meth(nmetal), &
             melo(nmetal), mgfac(nmetal), me1(nmetal), mie(nmetal), mit(nmetal), &
             mbreak(nmetal), mcap(NT, nmetal), mtab(nmetal), mblk(nmetal))
    !--- x nodes and T grid: the ones milne_setup built for H/He
    do ix = 1, NX
       t = real(ix-1, wp)/real(NX-1, wp)
       xg(ix) = XMAX*t*t
    end do
    do it = 1, NT
       lgT(it) = log10(TLO) + (log10(THI) - log10(TLO)) &
                 *real(it-1, wp)/real(NT-1, wp)
    end do
    do k = 1, nmetal
       call species_transition(k, ie, ist)
       mie(k) = ie;  mit(k) = ist
       meth(k) = elem_eth(ie, ist)
       melo(k) = max(meth(k), par%eion_min)
       call ground_term_weight(trim(elem_name(ie)), ist,   g_lo, xa)
       call ground_term_weight(trim(elem_name(ie)), ist+1, g_hi, e1_hi)
       mgfac(k) = g_lo/(2.0_wp*g_hi)
       me1(k)   = e1_hi
       mbreak(k) = species_outer_shell_break(ie, ist)
       mtab(k)   = species_photo_table_end(ie, ist) > 0.0_wp
       if (mtab(k)) call energy_block_edges(k)
       do it = 1, NT
          T_K   = 10.0_wp**lgT(it)
          kT_eV = kboltz_cgs*T_K/ev2erg
          !--- the capture cap alpha_1 <= alpha_rec (module header)
          a1 = milne_metal_alpha1(k, T_K)
          mcap(it, k) = min(1.0_wp, species_recombination_coefficient(ie, ist, T_K) &
                                    /max(a1, tinest))
          x_hi  = min(XMAX, (par%eion_max - melo(k))/kT_eV)
          if (x_hi <= 0.0_wp) then
             mcdf(:, it, k) = 1.0_wp
             mlnq(it, k)    = -huge(1.0_wp)
             cycle
          end if
          sum_e = 0.0_wp
          mcdf(1, it, k) = 0.0_wp
          if (mtab(k)) then
             !--- tabulated cross section: the energy blocks, each integrated
             !--- on the table's nodes (module header); the x-indexed CDF is
             !--- not used for its packets
             do ib = 1, mblk(k)%nb
                xa = (mblk(k)%e(ib-1) - melo(k))/kT_eV
                xb = (mblk(k)%e(ib)   - melo(k))/kT_eV
                fm = table_moment_x(k, kT_eV, melo(k), xa, xb, 1)
                if (fm > 0.0_wp) then
                   mblk(k)%lni(ib, it) = log(fm)
                else
                   mblk(k)%lni(ib, it) = -huge(1.0_wp)
                end if
                sum_e = sum_e + fm
             end do
             mcdf(2:NX, it, k) = sum_e
          else
             xbr = -1.0_wp
             if (mbreak(k) > melo(k)) xbr = (mbreak(k) - melo(k))/kT_eV
             fb = milne_metal_shape(k, 0.0_wp, kT_eV, melo(k))
             do ix = 1, NX-1
                xa = xg(ix)
                xb = min(xg(ix+1), x_hi)
                if (xb > xa) then
                   fa = fb                ! the previous segment's upper end
                   if (xbr > xa .and. xbr < xb) then
                      !--- the two outer-shell fits meet at einn without
                      !--- joining: one Simpson piece on each side
                      sum_e = sum_e + energy_piece(k, kT_eV, xa, &
                              xbr*(1.0_wp - 1.0e-12_wp), fa) &
                            + energy_piece(k, kT_eV, xbr, xb, &
                              milne_metal_shape(k, xbr, kT_eV, melo(k)))
                      fb = milne_metal_shape(k, xb, kT_eV, melo(k))
                   else
                      xm = 0.5_wp*(xa + xb)
                      fm = milne_metal_shape(k, xm, kT_eV, melo(k))
                      fb = milne_metal_shape(k, xb, kT_eV, melo(k))
                      sum_e = sum_e + (xb - xa)/6.0_wp &
                              *( fa*(melo(k) + xa*kT_eV) &
                               + 4.0_wp*fm*(melo(k) + xm*kT_eV) &
                               + fb*(melo(k) + xb*kT_eV) )
                   end if
                end if
                mcdf(ix+1, it, k) = sum_e
             end do
          end if
          if (sum_e > 0.0_wp) then
             mcdf(:, it, k) = mcdf(:, it, k)/sum_e
             !--- Int E^3 sigma e^-x dx [eV^3 cm^2] * kT [eV] -> erg^4 cm^2
             mlnq(it, k) = log(mcap(it, k)*milne_prefactor(mgfac(k), T_K) &
                               *sum_e*kT_eV*ev2erg**4)
          else
             mcdf(:, it, k) = 1.0_wp
             mlnq(it, k)    = -huge(1.0_wp)
          end if
          mcdf(NX, it, k) = 1.0_wp
       end do
    end do

    if (mpar%p_rank == 0) then
       write(*,'(a,i0,a)') ' MILNE: metal ground continua built for ', nmetal, &
          ' transitions (outer-shell sigma, g_n/2g_+ from ground_term_weights.txt)'
       !--- alpha_1 cannot exceed the radiative recombination it is part of
       worst = 0.0_wp;  kworst = 1
       big = ''
       do k = 1, nmetal
          a1  = milne_metal_alpha1(k, TCHK)
          arr = species_radiative_recombination(mie(k), mit(k), TCHK)
          ratio = a1/max(arr, tinest)
          if (ratio > worst) then
             worst = ratio;  kworst = k
          end if
          share = milne_metal_excited_share(k, TCHK)
          if (share > 0.1_wp .and. len_trim(big) < 480) &
             write(big,'(a,1x,a,i0)') trim(big), trim(elem_name(mie(k))), mit(k)
       end do
       write(*,'(a,f7.4,3a,i0,a)') ' MILNE: largest alpha_1/alpha_RR at 1e4 K = ', &
          worst, ' (', trim(elem_name(mie(kworst))), ' transition ', mit(kworst), ')'
       if (worst > 1.0_wp) write(*,'(a)') ' MILNE: (alpha_1 above the ' // &
          'radiative rate alone; the cap below is against RR + DR)'
       capped = ''
       do k = 1, nmetal
          if (minval(mcap(:, k), mask = lgT >= log10(par%te_min) .and. &
              lgT <= log10(par%te_max)) < 1.0_wp .and. len_trim(capped) < 480) &
             write(capped,'(a,1x,a,i0,a,f6.3,a)') trim(capped), &
                trim(elem_name(mie(k))), mit(k), ' (x', milne_metal_capture_cap(k, TCHK), ' at 1e4 K)'
       end do
       if (len_trim(capped) > 0) write(*,'(2a)') ' MILNE: alpha_1 held to ' // &
          'alpha_rec (RR + DR) within [te_min, te_max] for:', trim(capped)
       if (len_trim(big) > 0) write(*,'(2a)') ' MILNE: ground-term channel ' // &
          'an upper limit (excited ion terms carry > 10% of the Milne integral ' // &
          'at 1e4 K) for transitions:', trim(big)
    end if
  end subroutine milne_metal_setup

  !=========================================================================
  ! Energy blocks of a tabulated transition k (module header): edges at
  ! E_lo, every BLOCK_NODES-th table node, the table's end, einn, and every
  ! BLOCK_DE beyond the table, up to eion_max.
  !=========================================================================
  subroutine energy_block_edges(k)
    use species_mod, only : species_photo_table_nodes, &
                            species_photo_table_node_energy, species_photo_table_end
    integer, intent(in) :: k
    real(kind=wp), allocatable :: tmp(:)
    real(kind=wp) :: etab, etop, e
    integer :: j1, j2, j, n
    etab = species_photo_table_end(mie(k), mit(k))
    etop = par%eion_max
    call species_photo_table_nodes(mie(k), mit(k), melo(k), min(etab, etop), j1, j2)
    allocate(tmp(max(j2 - j1 + 1, 0)/BLOCK_NODES + &
                 int(max(etop - etab, 0.0_wp)/BLOCK_DE) + 8))
    n = 0
    tmp(1) = melo(k)
    do j = j1 + BLOCK_NODES - 1, j2, BLOCK_NODES
       n = n + 1
       tmp(n+1) = species_photo_table_node_energy(mie(k), mit(k), j)
    end do
    if (min(etab, etop) > tmp(n+1)) then
       n = n + 1
       tmp(n+1) = min(etab, etop)
    end if
    e = tmp(n+1)
    do while (e < etop)
       e = min(e + BLOCK_DE, etop)
       if (mbreak(k) > tmp(n+1) .and. mbreak(k) < e) then
          n = n + 1
          tmp(n+1) = mbreak(k)
       end if
       n = n + 1
       tmp(n+1) = e
    end do
    mblk(k)%nb = n
    allocate(mblk(k)%e(0:n), mblk(k)%lni(n, NT))
    mblk(k)%e(0:n) = tmp(1:n+1)
    mblk(k)%lni = -huge(1.0_wp)           ! rows with the band edge out of band
    deallocate(tmp)
  end subroutine energy_block_edges

  !--- photon energy [eV] of a packet of tabulated transition k at T_K from u:
  !--- the block from the block integrals interpolated linearly in 1/T in
  !--- their logarithm, the energy inside it from the node-stepped integral
  !--- at T_K (module header)
  real(kind=wp) function energy_block_sample(k, T_K, u) result(E)
    real(kind=wp), intent(in) :: T_K, u
    integer,       intent(in) :: k
    real(kind=wp) :: w, kT_eV, uu, lmax, tot, cum, pb, frac, la, lb
    real(kind=wp), allocatable :: lnib(:)
    integer :: it, ib, nb
    nb = mblk(k)%nb
    call milne_tbracket(T_K, it, w)
    !--- the bracket's weight in 1/T (milne_tbracket's is in log T)
    w  = (1.0_wp/T_K - 10.0_wp**(-lgT(it))) &
         /(10.0_wp**(-lgT(it+1)) - 10.0_wp**(-lgT(it)))
    w  = min(max(w, 0.0_wp), 1.0_wp)
    kT_eV = kboltz_cgs*T_K/ev2erg
    allocate(lnib(nb))
    do ib = 1, nb
       la = mblk(k)%lni(ib, it);  lb = mblk(k)%lni(ib, it+1)
       if (la < -1.0e300_wp .or. lb < -1.0e300_wp) then
          lnib(ib) = -huge(1.0_wp)
       else
          lnib(ib) = (1.0_wp - w)*la + w*lb
       end if
    end do
    lmax = maxval(lnib)
    tot = 0.0_wp
    do ib = 1, nb
       if (lnib(ib) > -1.0e300_wp) tot = tot + exp(lnib(ib) - lmax)
    end do
    uu  = min(max(u, 0.0_wp), 1.0_wp)*tot
    cum = 0.0_wp
    do ib = 1, nb
       pb = 0.0_wp
       if (lnib(ib) > -1.0e300_wp) pb = exp(lnib(ib) - lmax)
       if (cum + pb >= uu .and. pb > 0.0_wp) exit
       cum = cum + pb
    end do
    ib = min(ib, nb)
    frac = 0.0_wp
    if (pb > 0.0_wp) frac = (uu - cum)/pb
    E = melo(k) + kT_eV*table_invert_x(k, kT_eV, melo(k), &
           (mblk(k)%e(ib-1) - melo(k))/kT_eV, (mblk(k)%e(ib) - melo(k))/kT_eV, frac)
    E = min(E, par%eion_max*(1.0_wp - 1.0e-12_wp))
    deallocate(lnib)
  end function energy_block_sample

  !=========================================================================
  ! Node-stepped quadrature of a PHTABLE transition k (module header):
  ! Int_{xa}^{xb} E^m f(x) dx, f = E^2 sigma_outer(E) e^-x, E = e0 + x kT.
  !=========================================================================
  real(kind=wp) function table_moment_x(k, kT_eV, e0, xa, xb, m) result(acc)
    integer,       intent(in) :: k, m
    real(kind=wp), intent(in) :: kT_eV, e0, xa, xb
    real(kind=wp) :: xhit
    call table_walk(k, kT_eV, e0, xa, xb, m, -1.0_wp, acc, xhit)
  end function table_moment_x

  !--- the x in [xa, xb] at which Int_{xa}^{x} E f dx reaches the fraction
  !--- frac of Int_{xa}^{xb} E f dx (the energy-luminosity density the
  !--- packets follow), on the same pieces as table_moment_x
  real(kind=wp) function table_invert_x(k, kT_eV, e0, xa, xb, frac) result(x)
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: kT_eV, e0, xa, xb, frac
    real(kind=wp) :: total, acc
    call table_walk(k, kT_eV, e0, xa, xb, 1, -1.0_wp, total, x)
    x = xa
    if (total <= 0.0_wp) return
    call table_walk(k, kT_eV, e0, xa, xb, 1, &
                    min(max(frac, 0.0_wp), 1.0_wp)*total, acc, x)
  end function table_invert_x

  !--- One walk over the pieces of [xa, xb]: consecutive PHTABLE nodes, the
  !--- table's end and einn (where the fits jump: the piece below ends at
  !--- einn (1 - 1e-12)), none longer than DX_PIECE; Simpson on each.  With
  !--- target < 0 it returns the whole integral in acc; otherwise it stops
  !--- in the piece where acc reaches target and returns in xhit the point
  !--- that leaves the right share of that piece below it, the density taken
  !--- linear across the piece (a piece is <= 4 meV below eps = 27 eV,
  !--- where the table resolves its resonances, and <= DX_PIECE kT
  !--- elsewhere).
  subroutine table_walk(k, kT_eV, e0, xa, xb, m, target, acc, xhit)
    use species_mod, only : species_photo_table_nodes, &
                            species_photo_table_node_energy, species_photo_table_end
    integer,       intent(in)  :: k, m
    real(kind=wp), intent(in)  :: kT_eV, e0, xa, xb, target
    real(kind=wp), intent(out) :: acc, xhit
    real(kind=wp) :: x, xn, xnode, xend, xbr, fa, fm, fb, piece, r, w, slope
    integer :: j, j1, j2
    acc = 0.0_wp;  xhit = xb
    if (xb <= xa) return
    call species_photo_table_nodes(mie(k), mit(k), e0 + xa*kT_eV, &
                                   e0 + xb*kT_eV, j1, j2)
    xend = (species_photo_table_end(mie(k), mit(k)) - e0)/kT_eV
    xbr  = -1.0_wp
    if (mbreak(k) > 0.0_wp) xbr = (mbreak(k) - e0)/kT_eV
    j  = j1
    x  = xa
    fa = gm(x)
    do while (x < xb)
       xn = min(x + DX_PIECE, xb)
       do while (j <= j2)
          xnode = (species_photo_table_node_energy(mie(k), mit(k), j) - e0)/kT_eV
          if (xnode > x) exit
          j = j + 1
       end do
       if (j <= j2) then
          if (xnode < xn) xn = xnode
       end if
       if (xend > x .and. xend < xn) xn = xend
       if (xbr > x .and. xbr <= xn) then
          !--- the fits jump at einn: this piece ends just below it
          xn = xbr*(1.0_wp - 1.0e-12_wp)
          fm = gm(0.5_wp*(x + xn));  fb = gm(xn)
          piece = (xn - x)/6.0_wp*(fa + 4.0_wp*fm + fb)
          if (target >= 0.0_wp .and. acc + piece >= target) exit
          acc = acc + piece
          x  = xbr
          fa = gm(x)
          cycle
       end if
       fm = gm(0.5_wp*(x + xn));  fb = gm(xn)
       piece = (xn - x)/6.0_wp*(fa + 4.0_wp*fm + fb)
       if (target >= 0.0_wp .and. acc + piece >= target) exit
       acc = acc + piece
       x  = xn
       fa = fb
    end do
    if (target < 0.0_wp .or. x >= xb) return
    !--- inside [x, xn]: the share r of the piece, density linear fa -> fb
    r = (target - acc)/max(piece, tinest)
    r = min(max(r, 0.0_wp), 1.0_wp)
    w = xn - x
    slope = (fb - fa)/w
    if (r <= 0.0_wp) then
       xhit = x
    else if (abs(slope)*w < 1.0e-8_wp*max(fa + fb, tinest)) then
       xhit = x + r*w
    else
       !--- fa t + slope t^2/2 = r (fa + fb) w/2, the root in [0, w]
       xhit = x + 2.0_wp*r*0.5_wp*(fa + fb)*w &
              /(fa + sqrt(max(fa*fa + 2.0_wp*slope*r*0.5_wp*(fa + fb)*w, 0.0_wp)))
    end if
    xhit = min(max(xhit, x), xn)
    acc  = target
  contains
    real(kind=wp) function gm(xx) result(v)
      real(kind=wp), intent(in) :: xx
      v = milne_metal_shape(k, xx, kT_eV, e0)*(e0 + xx*kT_eV)**m
    end function gm
  end subroutine table_walk

  !--- Simpson piece of Int (E_lo + x kT) shape(x) dx over [xa, xb], with the
  !--- shape at xa given (the side of a break it belongs to)
  real(kind=wp) function energy_piece(k, kT_eV, xa, xb, fa) result(r)
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: kT_eV, xa, xb, fa
    real(kind=wp) :: xm, fm, fb
    xm = 0.5_wp*(xa + xb)
    fm = milne_metal_shape(k, xm, kT_eV, melo(k))
    fb = milne_metal_shape(k, xb, kT_eV, melo(k))
    r = (xb - xa)/6.0_wp*( fa*(melo(k) + xa*kT_eV) &
        + 4.0_wp*fm*(melo(k) + xm*kT_eV) + fb*(melo(k) + xb*kT_eV) )
  end function energy_piece

  !--- ground-term weight g and first excited-term energy e1 [eV] of stage
  !--- ist of element el, from data/atomic/ground_term_weights.txt
  subroutine ground_term_weight(el, ist, g, e1)
    use utility, only : fatal_error
    character(len=*), intent(in)  :: el
    integer,          intent(in)  :: ist
    real(kind=wp),    intent(out) :: g, e1
    character(len=256) :: line
    character(len=16)  :: key, name, term
    integer :: unit, ios, z, s, gi
    real(kind=wp) :: e
    logical :: found
    found = .false.
    open(newunit=unit, file=trim(par%atomic_dir)//'/ground_term_weights.txt', &
         status='old', iostat=ios)
    if (ios /= 0) call fatal_error('cannot open '//trim(par%atomic_dir)// &
       '/ground_term_weights.txt (metal ground continua of the diffuse field)')
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
       read(line,*) key, name, z, s, gi, term, e
       if (trim(key) == 'GROUND' .and. trim(name) == trim(el) .and. s == ist) then
          g = real(gi, wp);  e1 = e;  found = .true.
          exit
       end if
    end do
    close(unit)
    if (.not. found) call fatal_error('ground_term_weights.txt has no stage of '//trim(el))
  end subroutine ground_term_weight

  integer function milne_metal_count() result(n)
    n = nmetal
  end function milne_metal_count

  real(kind=wp) function milne_metal_ground_factor(k) result(g)
    integer, intent(in) :: k
    g = mgfac(k)
  end function milne_metal_ground_factor

  !=========================================================================
  ! In-band energy radiated by the ground captures of metal transition k per
  ! n_e n_+ [erg cm^3 s^-1]: Q(T) exp[-(E_lo - E_th)/kT] from the table,
  ! ln Q by a Catmull-Rom cubic in log T (linear interpolation left 4e-4 on
  ! Fe II at 1e4 K, where sigma rises steeply above threshold).
  !=========================================================================
  real(kind=wp) function milne_metal_band_coefficient(k, T_K) result(c)
    use mathlib, only : catmull_rom_uniform
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: T_K
    real(kind=wp) :: w, lq
    integer :: it
    c = 0.0_wp
    call milne_tbracket(T_K, it, w)
    lq = catmull_rom_uniform(mlnq(:, k), it, w)
    if (lq < -1.0e300_wp) return
    c  = exp(lq - (melo(k) - meth(k))*ev2erg/(kboltz_cgs*T_K))
  end function milne_metal_band_coefficient

  !--- the capture cap of transition k at T (alpha_rec/alpha_1, at most 1),
  !--- interpolated as the table carries it
  real(kind=wp) function milne_metal_capture_cap(k, T_K) result(c)
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: T_K
    real(kind=wp) :: w
    integer :: it
    call milne_tbracket(T_K, it, w)
    c = (1.0_wp - w)*mcap(it, k) + w*mcap(it+1, k)
  end function milne_metal_capture_cap

  !=========================================================================
  ! Photon energy [eV] of one packet of metal channel k at T_K from u, as
  ! milne_sample_energy: CDF interpolated in log T, inverted once.
  !=========================================================================
  real(kind=wp) function milne_metal_sample_energy(k, T_K, u) result(E)
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: T_K, u
    real(kind=wp) :: w, kT_eV, uu, ca, cb, cm, x
    integer :: it, lo, hi, mid
    if (mtab(k)) then
       E = energy_block_sample(k, T_K, u)
       return
    end if
    call milne_tbracket(T_K, it, w)
    kT_eV = kboltz_cgs*T_K/ev2erg
    uu    = min(max(u, 0.0_wp), 1.0_wp)
    lo = 1;  hi = NX
    ca = (1.0_wp - w)*mcdf(lo, it, k) + w*mcdf(lo, it+1, k)
    cb = (1.0_wp - w)*mcdf(hi, it, k) + w*mcdf(hi, it+1, k)
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       cm  = (1.0_wp - w)*mcdf(mid, it, k) + w*mcdf(mid, it+1, k)
       if (uu > cm) then
          lo = mid;  ca = cm
       else
          hi = mid;  cb = cm
       end if
    end do
    if (cb > ca) then
       x = xg(lo) + (xg(hi) - xg(lo))*(uu - ca)/(cb - ca)
    else
       x = xg(lo)
    end if
    E = melo(k) + x*kT_eV
    E = min(E, par%eion_max*(1.0_wp - 1.0e-12_wp))
  end function milne_metal_sample_energy

  !=========================================================================
  ! Direct quadratures, independent of the tables (setup report and gate):
  ! alpha_1(T) [cm^3/s] of metal transition k over all E >= E_th, and the
  ! energy Int_{e_lo}^{e_hi} E d alpha_1 [erg cm^3/s] (Simpson in
  ! x = (E - E_th)/kT, step 0.02, cut at x = XMAX).
  !=========================================================================
  real(kind=wp) function milne_metal_alpha1(k, T_K) result(a1)
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: T_K
    a1 = milne_metal_moment(k, T_K, 0, meth(k), huge(1.0_wp))
  end function milne_metal_alpha1

  real(kind=wp) function milne_metal_energy_integral(k, T_K, e_lo, e_hi) result(p)
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: T_K, e_lo, e_hi
    p = milne_metal_moment(k, T_K, 1, e_lo, e_hi)
  end function milne_metal_energy_integral

  !--- share of the Milne photon integral above E_th + E_1, where the outer-
  !--- shell fit also carries photoionization into excited ion terms
  !--- (-1 when E_1 is not available)
  real(kind=wp) function milne_metal_excited_share(k, T_K) result(f)
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: T_K
    f = -1.0_wp
    if (me1(k) <= 0.0_wp) return
    f = milne_metal_moment(k, T_K, 0, meth(k) + me1(k), huge(1.0_wp)) &
        /max(milne_metal_alpha1(k, T_K), tinest)
  end function milne_metal_excited_share

  !--- Int_{max(e_lo,E_th)}^{min(e_hi, E_th + XMAX kT)} E^m d alpha_1 with the
  !--- Milne prefactor, m = 0 [cm^3/s] or 1 [erg cm^3/s]
  real(kind=wp) function milne_metal_moment(k, T_K, m, e_lo, e_hi) result(r)
    integer,       intent(in) :: k, m
    real(kind=wp), intent(in) :: T_K, e_lo, e_hi
    real(kind=wp), parameter :: DX = 0.02_wp
    real(kind=wp) :: kT_eV, xlo, xhi, xa, xb, xm, acc, xbr
    kT_eV = kboltz_cgs*T_K/ev2erg
    xlo = max(0.0_wp, (e_lo - meth(k))/kT_eV)
    xhi = XMAX
    if (e_hi < huge(1.0_wp)) xhi = min(XMAX, (e_hi - meth(k))/kT_eV)
    if (mtab(k)) then
       acc = 0.0_wp
       if (xhi > xlo) acc = table_moment_x(k, kT_eV, meth(k), xlo, xhi, m)
       r = milne_prefactor(mgfac(k), T_K)*acc*kT_eV*ev2erg**(3+m)
       return
    end if
    xbr = -1.0_wp
    if (mbreak(k) > 0.0_wp) xbr = (mbreak(k) - meth(k))/kT_eV
    acc = 0.0_wp
    xa = xlo
    do while (xa < xhi)
       xb = min(xa + DX, xhi)
       if (xbr > xa .and. xbr < xb) then
          !--- end the step at einn, approached from below, and restart there
          xb = xbr*(1.0_wp - 1.0e-12_wp)
          xm = 0.5_wp*(xa + xb)
          acc = acc + (xb - xa)/6.0_wp*( moment_integrand(xa) &
                + 4.0_wp*moment_integrand(xm) + moment_integrand(xb) )
          xa = xbr
          cycle
       end if
       xm = 0.5_wp*(xa + xb)
       acc = acc + (xb - xa)/6.0_wp*( moment_integrand(xa) &
             + 4.0_wp*moment_integrand(xm) + moment_integrand(xb) )
       xa = xb
    end do
    r = milne_prefactor(mgfac(k), T_K)*acc*kT_eV*ev2erg**(3+m)
  contains
    real(kind=wp) function moment_integrand(x) result(f)
      real(kind=wp), intent(in) :: x
      real(kind=wp) :: E
      E = meth(k) + x*kT_eV
      f = milne_metal_shape(k, x, kT_eV, meth(k))*E**m
    end function moment_integrand
  end function milne_metal_moment

end module milne_recomb_spectrum_mod
