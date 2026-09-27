module diffuse_mod
!---------------------------------------------------------------------------
! MoCHII: explicit diffuse ionizing field.
!
! Case B absorbs ground-level recombination photons on the spot.
! The explicit field replaces that: the equilibrium runs with CASE A rates
! and every recombination to the ground level emits an ionizing packet,
! transported exactly like a stellar packet.  Three continua:
!     H II  + e -> H I(1s)   : E >= 13.598 eV      (channel 1)
!     He II + e -> He I(1^1S): E >= 24.587 eV      (channel 2)
!     He III+ e -> He II(1s) : E >= 54.416 eV      (channel 3)
! and two EXCITED-level cascades, which the ground continua leave out and
! which are ionizing because the decay energies exceed 13.598 eV:
!     He I  n=2 terms        : channel 4, par%hei_diffuse  (hei_cascade_mod)
!     He III cascade         : channel 5, par%heii_diffuse (heii_cascade_mod)
! Channels 1-4 run with par%recombination_continuum = 'transported'.
! Channel 5 runs with it too, and also alone with 'on_the_spot' when
! par%heii_diffuse is set (heii_cascade_transported): the He III cascade
! transported, every other recombination photon absorbed on the spot.
! with rates alpha_1 = alpha_A - alpha_B per ground channel.  The ground
! photon energies follow par%diffuse_energy_model: 'milne' (default) samples
! the Milne free-bound spectrum E^2 sigma_pi(E) exp[-(E-E_th)/kT] on the
! transport cross sections (milne_recomb_spectrum_mod); 'exponential' draws
! E = E_th + kT_e x, x ~ Exp(1); 'threshold' puts every photon at E_th.
!
! Each iteration the emission table is rebuilt from the current state:
! L_ch(leaf) = n_e n_recomb alpha_1(T) <E_ch> V_cm3, with <E_ch> the mean
! of the same energy distribution (ground_mean_energy).
! The diffuse packet count is n_diff = max(1, nint(L_diffuse / Lpacket)),
! Lpacket that of one ionizing source packet (the argument of diffuse_build:
! ion_band_mod::ionizing_packet_luminosity), and every diffuse packet carries
! L_diffuse / n_diff, so they carry exactly L_diffuse; sampling: leaf CDF over
! total leaf luminosity, then a channel CDF within the leaf.
!---------------------------------------------------------------------------
  use define
  use octree_mod,    only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHII, gas_xHeI, gas_xHeII, &
                            gas_xHeIII, gas_ne, gas_Te, gas_nleaf
  use recomb_mod
  use milne_recomb_spectrum_mod, only : milne_setup, milne_sample_energy, &
                                        milne_energy_per_recomb, milne_is_ready
  use hei_twophoton_mod, only : hei_2ph_sample_energy, &
                               hei_2ph_energy_per_decay
  use hei_recombination_photons_mod, only : hei_excited_decays, &
                                           hei_excited_decays_coupled, &
                                           HEI_E3S, HEI_E1P
  use heii_cascade_mod,  only : heii_branch_luminosity, heii_fb_sample_energy, &
                                heii_2q_sample_energy, heii_lyman_sample_energy, &
                                heii_cascade_transported
  implicit none
  private

  public :: diffuse_build, gen_diffuse_photon, gen_diffuse_photon_qmc, &
            diffuse_nphot, diffuse_lum
  public :: hei_branch_energies, hei_branch_energies_coupled

  integer, parameter :: NCH = 5                 ! diffuse emission channels
  real(kind=wp), allocatable :: dif_cdf(:)      ! leaf CDF (total energy)
  real(kind=wp), allocatable :: dif_ch(:,:)     ! (NCH, nleaf) channel luminosities
  !--- channel-4 branch weights (2^3S, 2^1P, 2^1S) of each leaf when the
  !--- 2^3S population is coupled (par%hei_metastable_coupling): they depend
  !--- on n_3 and Phi_3, so the launchers read what diffuse_build computed
  !--- instead of recomputing them.
  real(kind=wp), allocatable :: dif_hei(:,:)    ! (3, nleaf) [eV cm^-3 s^-1]
  real(kind=wp) :: diffuse_lum = 0.0_wp         ! total [erg/s]
  !--- luminosity carried by one diffuse packet [erg/s]: the diffuse
  !--- luminosity over the packet count diffuse_build set from the lpacket it
  !--- was given, so the packets carry exactly diffuse_lum (before
  !--- 2026-09-25 each carried lpacket itself, and the rounded count left the
  !--- transported total off by up to lpacket/2, or dropped a field worth
  !--- less than half a packet).
  real(kind=wp) :: diffuse_packet_lum = 0.0_wp
  integer(kind=int64) :: diffuse_nphot = 0

  !--- He I excited channel (channel 4, par%hei_diffuse): the 2^3S decay is a
  !--- single 19.82 eV photon and the 2^1P decay the 21.22 eV (584 A) line
  !--- (HEI_E3S, HEI_E1P); the exits of each recombination come from
  !--- hei_recombination_photons_mod, which the case-B on-the-spot sources
  !--- use too, with the branch fractions of hei_cascade_mod (Benjamin,
  !--- Skillman & Smits 1999 with J. S. Mathis).
  !--- The in-band energy radiated per 2^1S decay now comes from the Drake
  !--- et al. (1969) shape (hei_twophoton_mod), not from a fixed 0.56 photons
  !--- times the mean of a flat spectrum.  The photon count was right to
  !--- 0.6%, but the flat surrogate overestimated the mean in-band photon
  !--- energy by 6% and the radiated in-band energy by 6.5%.

contains

  !=========================================================================
  subroutine diffuse_build(lpacket)
    implicit none
    real(kind=wp), intent(in) :: lpacket
    real(kind=wp) :: T, ne, nH, vol, kT_eV, a1
    real(kind=wp) :: be3s, be1p, be1s
    real(kind=wp) :: lfb, llyman, l2q
    real(kind=wp) :: nHII, nHeII_n, nHeIII, xHeIII
    integer :: il, ic
    logical :: all_channels

    if (allocated(dif_cdf)) then
       if (size(dif_cdf) /= gas_nleaf) deallocate(dif_cdf, dif_ch)
    end if
    if (.not. allocated(dif_cdf)) then
       allocate(dif_cdf(gas_nleaf), dif_ch(NCH, gas_nleaf))
    end if
    if (par%hei_diffuse .and. par%hei_metastable_coupling) then
       if (allocated(dif_hei)) then
          if (size(dif_hei, 2) /= gas_nleaf) deallocate(dif_hei)
       end if
       if (.not. allocated(dif_hei)) allocate(dif_hei(3, gas_nleaf))
    end if

    !--- the ground continua and the He I cascade are transported only with
    !--- the transported continuum; with the on-the-spot one channel 5 may be
    !--- transported alone (heii_cascade_transported)
    all_channels = trim(par%recombination_continuum) == 'transported'
    do il = 1, gas_nleaf
       nH = gas_nH(il)
       dif_ch(:,il) = 0.0_wp
       if (allocated(dif_hei)) dif_hei(:,il) = 0.0_wp
       if (nH <= 0.0_wp) cycle
       T  = gas_Te(il)
       ne = gas_ne(il)
       vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
       kT_eV  = kboltz_cgs*T/ev2erg
       xHeIII = gas_xHeIII(il)
       nHII    = nH*gas_xHII(il)
       nHeII_n = nH*par%abund_He*gas_xHeII(il)
       nHeIII  = nH*par%abund_He*xHeIII
       !--- channel energy luminosities [erg/s]
       if (all_channels) then
          a1 = alphaA_HII(T)   - alphaB_HII(T)
          dif_ch(1,il) = ne*nHII   *a1*ground_mean_energy(1, eth_HI, kT_eV, T)*ev2erg*vol
          a1 = alphaA_HeII(T)  - alphaB_HeII(T)
          dif_ch(2,il) = ne*nHeII_n*a1*ground_mean_energy(2, eth_HeI, kT_eV, T)*ev2erg*vol
          a1 = alphaA_HeIII(T) - alphaB_HeIII(T)
          dif_ch(3,il) = ne*nHeIII *a1*ground_mean_energy(3, eth_HeII, kT_eV, T)*ev2erg*vol
          !--- channel 4: He I EXCITED-level recombination radiation
          !--- (par%hei_diffuse): all case-B He II recombinations cascade to
          !--- n = 2 and decay with the hei_cascade_mod branch fractions; only
          !--- the H-ionizing part carries band luminosity.
          if (par%hei_diffuse .and. par%hei_metastable_coupling) then
             !--- with the 2^3S population coupled, the branch energies are
             !--- per unit volume (hei_branch_energies_coupled)
             call hei_branch_energies_coupled(il, T, ne, nH*gas_xHI(il), &
                  nH*par%abund_He*gas_xHeI(il), nHeII_n, gas_xHI(il), &
                  gas_xHeI(il), be3s, be1p, be1s)
             dif_hei(:,il) = [be3s, be1p, be1s]
             dif_ch(4,il) = vol*ev2erg*(be3s + be1p + be1s)
          else if (par%hei_diffuse) then
             a1 = alphaB_HeII(T)
             call hei_branch_energies(T, gas_xHI(il), gas_xHeI(il), &
                                      be3s, be1p, be1s)
             dif_ch(4,il) = ne*nHeII_n*a1*vol*ev2erg*(be3s + be1p + be1s)
          end if
       end if
       !--- channel 5: He III cascade through He II (par%heii_diffuse), in
       !--- the case of the Lyman lines.  The free-bound term already carries
       !--- its own coefficient [erg cm^3/s]; the Lyman-line and 2^2S exits
       !--- are branches of alphaB_HeIII, which heii_branch_luminosity applies.
       if (heii_cascade_transported()) then
          call heii_branch_luminosity(T, lfb, llyman, l2q)
          dif_ch(5,il) = ne*nHeIII*vol*(lfb + llyman + l2q)
       end if
    end do

    !--- leaf CDF over total energy
    diffuse_lum = 0.0_wp
    do il = 1, gas_nleaf
       diffuse_lum = diffuse_lum + sum(dif_ch(:,il))
       dif_cdf(il) = diffuse_lum
    end do
    if (diffuse_lum > 0.0_wp) then
       dif_cdf = dif_cdf/diffuse_lum
       diffuse_nphot = max(1_int64, nint(diffuse_lum/lpacket, int64))
       diffuse_packet_lum = diffuse_lum/real(diffuse_nphot, wp)
    else
       diffuse_nphot = 0
       diffuse_packet_lum = lpacket
    end if
  end subroutine diffuse_build

  !=========================================================================
  subroutine gen_diffuse_photon(photon)
    use random,       only : rand_number
    use ion_band_mod, only : ion_bin_of
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp) :: u, cost, sint, phi, half, eph, kT_eV, w(NCH), ub
    real(kind=wp) :: be3s, be1p, be1s
    integer :: lo, hi, mid, il, ic, ch

    !--- leaf from the CDF (binary search)
    u = rand_number()
    lo = 1;  hi = gas_nleaf
    do while (lo < hi)
       mid = (lo + hi)/2
       if (dif_cdf(mid) < u) then
          lo = mid + 1
       else
          hi = mid
       end if
    end do
    il = lo
    half = leaf_half(il)

    !--- uniform position within the leaf
    photon%x = leaf_cx(il) + (2.0_wp*rand_number() - 1.0_wp)*half
    photon%y = leaf_cy(il) + (2.0_wp*rand_number() - 1.0_wp)*half
    photon%z = leaf_cz(il) + (2.0_wp*rand_number() - 1.0_wp)*half

    cost = 2.0_wp*rand_number() - 1.0_wp
    sint = sqrt(1.0_wp - cost*cost)
    phi  = twopi*rand_number()
    photon%kx = sint*cos(phi)
    photon%ky = sint*sin(phi)
    photon%kz = cost

    !--- channel; ground continua sample par%diffuse_energy_model
    !--- (ground_photon_energy); the two excited channels sample their path by energy luminosity and then the
    !--- energy within it -- channel 4 the He I decay branch (fixed line
    !--- energies, or the Drake et al. 1969 two-photon shape), channel 5 the
    !--- He III cascade path (free-bound CDF, the 40.813 eV line, or the
    !--- Nussbaumer & Schmutz two-photon shape).
    w  = dif_ch(:,il)
    ch = pick_channel(w, rand_number())
    if (ch == 1) then
       eph = eth_HI
    else if (ch == 2) then
       eph = eth_HeI
    else if (ch == 3) then
       eph = eth_HeII
    end if
    if (ch <= 3) then
       kT_eV = kboltz_cgs*gas_Te(il)/ev2erg
       eph = ground_photon_energy(ch, eph, kT_eV, gas_Te(il), rand_number())
    else if (ch == 4) then
       if (par%hei_metastable_coupling) then
          be3s = dif_hei(1,il);  be1p = dif_hei(2,il);  be1s = dif_hei(3,il)
       else
          call hei_branch_energies(gas_Te(il), gas_xHI(il), gas_xHeI(il), &
                                   be3s, be1p, be1s)
       end if
       ub = rand_number()*(be3s + be1p + be1s)
       if (ub <= be3s) then
          eph = HEI_E3S
       else if (ub <= be3s + be1p) then
          eph = HEI_E1P
       else
          eph = hei_2ph_sample_energy(rand_number())
       end if
    else
       !--- drawn into named variables first: argument evaluation order is not
       !--- specified, and the two draws are the u(8), u(9) of the Sobol
       !--- launcher in that order.
       ub = rand_number()
       u  = rand_number()
       eph = heii_cascade_energy(gas_Te(il), ub, u)
    end if

    !--- bin index by binary search of the band edges (the one source of truth
    !--- for BOTH the aligned and the legacy-log grid; recombination photons all
    !--- carry eph >= the channel threshold > eion_min, so they land in the
    !--- ionizing segment).  The old closed-form log formula assumed a pure
    !--- log grid and mis-binned under par%ion_align_edges.
    photon%inu = ion_bin_of(eph)

    photon%wgt     = 1.0_wp
    photon%Lpacket = diffuse_packet_lum
    !--- the packet carries the energy the emission channel actually sampled;
    !--- inu is only the diagnostic bin it lands in.
    photon%energy_eV = eph
    photon%nscatt  = 0
    photon%inside  = .true.
    photon%from_external = .false.      ! diffuse packets peel isotropically
    photon%icell_amr = il
  end subroutine gen_diffuse_photon

  !=========================================================================
  ! Quasi-random diffuse launch.  Nine launch uniforms with FIXED dimension
  ! semantics (md/QUASI_RANDOM_LAUNCH.md), from the diffuse scramble stream:
  !   u(1)   emitting leaf (dif_cdf inverse);
  !   u(2-4) uniform position within the leaf;
  !   u(5)   mu = 2u-1;                 u(6) phi = 2 pi u;
  !   u(7)   emission channel (dif_ch inverse over the NCH channels);
  !   u(8)   first energy variate: channels 1-3 the ground-continuum energy
  !          (ground_photon_energy: the Milne CDF inverse by default);
  !          channel 4 the He I decay-branch pick; channel 5 the He II
  !          cascade-path pick (free-bound / Ly-alpha / two-photon);
  !   u(9)   second energy variate: the He I two-photon energy from the
  !          tabulated Drake shape, or the channel-5 free-bound / two-photon
  !          energy (else unused).
  ! The fifth channel reuses dimensions 8 and 9 rather than adding any, since
  ! its path pick and its energy draw are the same two decisions channel 4
  ! makes.
  ! Every expression matches gen_diffuse_photon; the bin comes from ion_bin_of.
  !=========================================================================
  subroutine gen_diffuse_photon_qmc(photon, u)
    use ion_band_mod, only : ion_bin_of
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp),     intent(in)    :: u(:)
    real(kind=wp) :: uu, cost, sint, phi, half, eph, kT_eV, w(NCH), ub
    real(kind=wp) :: be3s, be1p, be1s
    integer :: lo, hi, mid, il, ch

    !--- leaf from the CDF (binary search)
    uu = u(1)
    lo = 1;  hi = gas_nleaf
    do while (lo < hi)
       mid = (lo + hi)/2
       if (dif_cdf(mid) < uu) then
          lo = mid + 1
       else
          hi = mid
       end if
    end do
    il = lo
    half = leaf_half(il)

    !--- uniform position within the leaf
    photon%x = leaf_cx(il) + (2.0_wp*u(2) - 1.0_wp)*half
    photon%y = leaf_cy(il) + (2.0_wp*u(3) - 1.0_wp)*half
    photon%z = leaf_cz(il) + (2.0_wp*u(4) - 1.0_wp)*half

    cost = 2.0_wp*u(5) - 1.0_wp
    sint = sqrt(1.0_wp - cost*cost)
    phi  = twopi*u(6)
    photon%kx = sint*cos(phi)
    photon%ky = sint*sin(phi)
    photon%kz = cost

    !--- channel; ground continua sample par%diffuse_energy_model
    !--- (ground_photon_energy); the He I excited channel (4) samples its decay branch by energy
    !--- luminosity (fixed line energies; the two-photon branch from the
    !--- tabulated Drake shape).
    w  = dif_ch(:,il)
    ch = pick_channel(w, u(7))
    if (ch == 1) then
       eph = eth_HI
    else if (ch == 2) then
       eph = eth_HeI
    else if (ch == 3) then
       eph = eth_HeII
    end if
    if (ch <= 3) then
       kT_eV = kboltz_cgs*gas_Te(il)/ev2erg
       eph = ground_photon_energy(ch, eph, kT_eV, gas_Te(il), u(8))
    else if (ch == 4) then
       if (par%hei_metastable_coupling) then
          be3s = dif_hei(1,il);  be1p = dif_hei(2,il);  be1s = dif_hei(3,il)
       else
          call hei_branch_energies(gas_Te(il), gas_xHI(il), gas_xHeI(il), &
                                   be3s, be1p, be1s)
       end if
       ub = u(8)*(be3s + be1p + be1s)
       if (ub <= be3s) then
          eph = HEI_E3S
       else if (ub <= be3s + be1p) then
          eph = HEI_E1P
       else
          eph = hei_2ph_sample_energy(u(9))
       end if
    else
       eph = heii_cascade_energy(gas_Te(il), u(8), u(9))
    end if

    photon%inu = ion_bin_of(eph)

    photon%wgt     = 1.0_wp
    photon%Lpacket = diffuse_packet_lum
    photon%energy_eV = eph
    photon%nscatt  = 0
    photon%inside  = .true.
    photon%from_external = .false.      ! diffuse packets peel isotropically
    photon%icell_amr = il
  end subroutine gen_diffuse_photon_qmc

  !--- In-band energy radiated per ground-state recombination [eV].  This is
  !--- what turns the recombination rate into the channel's energy luminosity,
  !--- so it MUST be the mean of the same distribution gen_diffuse_photon
  !--- samples; otherwise the emitted photon rate no longer equals
  !--- n_e n_ion alpha_1 V.
  !--- Energy radiated per He II recombination in each H-ionizing branch of
  !--- the He I cascade [eV], from the tabulated alpha_eff (hei_cascade_mod)
  !--- rather than the statistical 3/4 : 1/4 x 2/3 : 1/3 weights.  One routine
  !--- for the luminosity table and for both launchers, so the branch CDF and
  !--- the emitted luminosity cannot disagree.
  subroutine hei_branch_energies(T_K, xHI, xHeI, e3s, e1p, e1s)
    real(kind=wp), intent(in)  :: T_K, xHI, xHeI
    real(kind=wp), intent(out) :: e3s, e1p, e1s
    real(kind=wp) :: r3s, r1p, r2q
    call hei_excited_decays(T_K, xHI, xHeI, r3s, r1p, r2q)
    e3s = r3s*HEI_E3S
    e1p = r1p*HEI_E1P
    !--- A converted decay radiates the 2^1S pair, 20.62 eV, not 21.22 eV;
    !--- the remaining 0.60 eV leaves as 2^1P -> 2^1S at 2.06 micron, which is
    !--- non-ionizing and so correctly absent from the band.
    e1s = r2q*hei_2ph_energy_per_decay()
  end subroutine hei_branch_energies

  !--- The same three branch energies with the 2^3S population coupled
  !--- (par%hei_metastable_coupling), per unit volume [eV cm^-3 s^-1]; the
  !--- exits are hei_excited_decays_coupled's (hei_recombination_photons_mod,
  !--- where the physics is written).
  subroutine hei_branch_energies_coupled(il, T_K, ne, nHI, nHeI, nHeII, &
                                         xHI, xHeI, e3s, e1p, e1s)
    integer,       intent(in)  :: il
    real(kind=wp), intent(in)  :: T_K, ne, nHI, nHeI, nHeII, xHI, xHeI
    real(kind=wp), intent(out) :: e3s, e1p, e1s
    real(kind=wp) :: r3s, r1p, r2q
    call hei_excited_decays_coupled(il, T_K, ne, nHI, nHeI, nHeII, xHI, xHeI, &
                                    r3s, r1p, r2q)
    e3s = r3s*HEI_E3S
    e1p = r1p*HEI_E1P
    e1s = r2q*hei_2ph_energy_per_decay()
  end subroutine hei_branch_energies_coupled

  !--- Emission channel from the uniform u: the first channel whose cumulative
  !--- energy luminosity passes u, over the channels that carry any.  A
  !--- channel switched off has weight zero and is skipped rather than left to
  !--- a final else, so a rounding step past the cumulative sum cannot select
  !--- it and reach a sampler that was never set up.  The comparison is the
  !--- same expression the if-chain used, so the draw is unchanged.
  integer function pick_channel(w, u) result(ch)
    real(kind=wp), intent(in) :: w(NCH), u
    real(kind=wp) :: acc, tot
    integer :: i
    tot = sum(w)
    acc = 0.0_wp
    ch  = 1
    do i = 1, NCH
       if (w(i) <= 0.0_wp) cycle
       acc = acc + w(i)
       ch  = i
       if (u*tot <= acc) return
    end do
  end function pick_channel

  !--- Photon energy [eV] of one He III cascade packet (channel 5).  The path
  !--- is picked by ENERGY luminosity from the same three terms
  !--- heii_branch_luminosity puts in the emission table, so the branch CDF
  !--- and the channel luminosity cannot disagree.  One routine for both
  !--- launchers, for the same reason.
  real(kind=wp) function heii_cascade_energy(T_K, u_path, u_energy) result(eph)
    real(kind=wp), intent(in) :: T_K, u_path, u_energy
    real(kind=wp) :: lfb, llyman, l2q, ub
    call heii_branch_luminosity(T_K, lfb, llyman, l2q)
    ub = u_path*(lfb + llyman + l2q)
    if (ub <= lfb) then
       eph = heii_fb_sample_energy(T_K, u_energy)
    else if (ub <= lfb + llyman) then
       eph = heii_lyman_sample_energy(T_K, u_energy)
    else
       eph = heii_2q_sample_energy(u_energy)
    end if
  end function heii_cascade_energy

  real(kind=wp) function ground_mean_energy(ich, eth, kT_eV, T_K) result(mean_e)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: eth, kT_eV, T_K
    select case (trim(par%diffuse_energy_model))
    case ('milne')
       mean_e = milne_energy_per_recomb(ich, T_K)
    case ('exponential')
       mean_e = eth + kT_eV
    case default                          ! 'threshold'
       mean_e = eth
    end select
  end function ground_mean_energy

  !--- Photon energy [eV] of one ground-continuum packet from the uniform u.
  !--- One routine for both the pseudorandom and the Sobol launch, so the two
  !--- cannot drift apart.
  real(kind=wp) function ground_photon_energy(ich, eth, kT_eV, T_K, u) result(eph)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: eth, kT_eV, T_K, u
    select case (trim(par%diffuse_energy_model))
    case ('milne')
       eph = milne_sample_energy(ich, T_K, u)
    case ('exponential')
       eph = eth + kT_eV*(-log(max(u, tinest)))
    case default                          ! 'threshold'
       eph = eth
    end select
    eph = min(eph, par%eion_max*0.999_wp)
  end function ground_photon_energy

end module diffuse_mod
