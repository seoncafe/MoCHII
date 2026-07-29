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
! alpha_A - alpha_B that sets the packet count.  They come from unrelated
! data (Milne + Verner cross sections versus the fitted recombination
! coefficients) and are reported, never rescaled into agreement.
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
! log T first and inverted once.  Inverting two neighbouring temperature rows
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
!---------------------------------------------------------------------------
  use define
  use photo_xsec, only : sigma_HI, sigma_HeI, sigma_HeII
  implicit none
  private

  public :: milne_setup, milne_sample_energy, milne_energy_per_recomb
  public :: milne_is_ready
  public :: MILNE_HI, MILNE_HEI, MILNE_HEII

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
  ! Build the tables and report the alpha_1 cross-check.  Call once, after
  ! the band setup has fixed par%eion_max.
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

end module milne_recomb_spectrum_mod
