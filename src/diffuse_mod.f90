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
! and two continua of every ion that the transported continuum carries too
! (since 2026-10-11; plan md/ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md
! Phase 4, both on because the photons exist and ionize: in MOCASSIN v2's
! diffuse field, and absent here before):
!     X^(i+1) + e -> X^i(ground), every metal transition: channel 6
!         (milne_recomb_spectrum_mod: the Milne spectrum on the outer-shell
!         transport cross section, rate and spectrum from the same integral)
!     free-free of H+, He+, He2+ and (par%metal_freefree) the metal ions
!         inside the band: channel 7 (recomb_continuum_mod::freefree_band_*,
!         the coefficient of the written continuum)
! and, under par%metal_recombination_photons = 'valence_fractions' only,
!     X^(i+1) + e -> X^i, RR + DR, X^i with >= 3 electrons: channel 8, one
!         photon at 0.75 and one at 0.25 of the ionization potential of X^i
!         (the Cloudy c25.00 approximation, kept for comparison, not a
!         physical model: no atomic data; with channel 6 a ground capture is
!         counted twice; see define.f90).  Photons outside
!         [eion_min, eion_max) are not launched and their luminosity is
!         printed.
! Channels 1-4 and 6-8 run with par%recombination_continuum = 'transported'.
! With 'escaping' none is transported (every recombination photon leaves).
! With 'on_the_spot' channels 6 and 7 are not transported and their photons
! are not absorbed on the spot either: they leave the nebula, as with
! 'escaping'.  The on-the-spot mode absorbs locally the H and He photons
! the case-B rates assume; the metal ground photons are absorbed mostly by
! H0 and He0, not by the ion that emits them, so a local metal case B would
! be wrong, and an on-the-spot deposition into H0/He0 was not built for an
! approximation kept for comparison.  Their size is the channel 6 + 7 share
! the transported run prints (results in
! md/ATOMIC_PHYSICS_RESTORATION_RESULTS_2026-10-11.md, Phase 4).
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
  use mpi
  use memory_mod,    only : create_shared_mem, destroy_shared_mem, get_window, &
                           shared_host_rank, shared_host_size, shared_host_comm
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
  use milne_recomb_spectrum_mod, only : milne_metal_count, &
                                        milne_metal_band_coefficient, &
                                        milne_metal_sample_energy
  use recomb_continuum_mod, only : freefree_band_coefficient, &
                                   freefree_band_sample_energy, FF_ZMAX
  implicit none
  private

  public :: diffuse_build, gen_diffuse_photon, gen_diffuse_photon_qmc, &
            diffuse_nphot, diffuse_lum
  public :: hei_branch_energies, hei_branch_energies_coupled
  public :: diffuse_channel_luminosity, metal_valence_luminosity_outside

  integer, parameter :: NCH = 8                 ! diffuse emission channels
  integer, parameter :: CH_METAL_GROUND = 6, CH_FREEFREE = 7, CH_METAL_VALENCE = 8
  !--- channel-8 photon energies as fractions of the ionization potential of
  !--- the recombined ion (Cloudy c25.00 cont_createpointers.cpp)
  real(kind=wp), parameter :: VALENCE_FRACTION(2) = [0.75_wp, 0.25_wp]
  !--- total luminosity of each channel in the last diffuse_build [erg/s]
  real(kind=wp) :: diffuse_channel_luminosity(NCH) = 0.0_wp
  !--- channel-8 luminosity below eion_min and at or above eion_max in the
  !--- last diffuse_build [erg/s], world rank 0 only (not launched)
  real(kind=wp) :: metal_valence_luminosity_outside(2) = 0.0_wp
  !--- the leaf tables are the same on every rank: node-shared (2026-10-12),
  !--- the node's ranks filling disjoint leaves
  real(kind=wp), pointer, contiguous :: dif_cdf(:) => null()    ! leaf CDF (total energy)
  real(kind=wp), pointer, contiguous :: dif_ch(:,:) => null()   ! (NCH, nleaf) channel luminosities
  !--- channel-4 branch weights (2^3S, 2^1P, 2^1S) of each leaf when the
  !--- 2^3S population is coupled (par%hei_metastable_coupling): they depend
  !--- on n_3 and Phi_3, so the launchers read what diffuse_build computed
  !--- instead of recomputing them.
  real(kind=wp), pointer, contiguous :: dif_hei(:,:) => null()  ! (3, nleaf) [eV cm^-3 s^-1]
  integer :: dif_win(3) = MPI_WIN_NULL     ! dif_cdf, dif_ch, dif_hei
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
    integer :: il, ic, ierr
    logical :: all_channels

    if (associated(dif_cdf)) then
       if (size(dif_cdf) /= gas_nleaf) then
          call destroy_shared_mem(dif_cdf);  call destroy_shared_mem(dif_ch)
       end if
    end if
    if (.not. associated(dif_cdf)) then
       call create_shared_mem(dif_cdf, [gas_nleaf]);  call get_window(dif_cdf, dif_win(1))
       call create_shared_mem(dif_ch, [NCH, gas_nleaf]);  call get_window(dif_ch, dif_win(2))
    end if
    if (par%hei_diffuse .and. par%hei_metastable_coupling) then
       if (associated(dif_hei)) then
          if (size(dif_hei, 2) /= gas_nleaf) call destroy_shared_mem(dif_hei)
       end if
       if (.not. associated(dif_hei)) then
          call create_shared_mem(dif_hei, [3, gas_nleaf])
          call get_window(dif_hei, dif_win(3))
       end if
    end if

    !--- the ground continua and the He I cascade are transported only with
    !--- the transported continuum; with the on-the-spot one channel 5 may be
    !--- transported alone (heii_cascade_transported)
    all_channels = trim(par%recombination_continuum) == 'transported'
    do il = shared_host_rank() + 1, gas_nleaf, shared_host_size()
       nH = gas_nH(il)
       dif_ch(:,il) = 0.0_wp
       if (associated(dif_hei)) dif_hei(:,il) = 0.0_wp
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
       !--- channels 6 and 7: the metal ground continua and the free-free
       !--- emission inside the band
       if (all_channels) then
          dif_ch(CH_METAL_GROUND,il) = ne*vol*sum(metal_ground_weights(il, T))
          dif_ch(CH_FREEFREE,il)     = ne*vol*sum(freefree_weights(il, T))
          if (trim(par%metal_recombination_photons) == 'valence_fractions') &
             dif_ch(CH_METAL_VALENCE,il) = ne*vol*sum(metal_valence_weights(il, T))
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

    call diffuse_shared_sync()

    !--- leaf CDF over total energy; every rank forms the same sums from the
    !--- shared table, the node's rank 0 stores the CDF
    diffuse_lum = 0.0_wp
    diffuse_channel_luminosity = 0.0_wp
    do il = 1, gas_nleaf
       diffuse_lum = diffuse_lum + sum(dif_ch(:,il))
       if (shared_host_rank() == 0) dif_cdf(il) = diffuse_lum
       diffuse_channel_luminosity = diffuse_channel_luminosity + dif_ch(:,il)
    end do
    if (diffuse_lum > 0.0_wp .and. shared_host_rank() == 0) dif_cdf = dif_cdf/diffuse_lum
    call diffuse_shared_sync()
    metal_valence_luminosity_outside = 0.0_wp
    if (all_channels .and. mpar%p_rank == 0 .and. &
        trim(par%metal_recombination_photons) == 'valence_fractions') then
       do il = 1, gas_nleaf
          if (gas_nH(il) <= 0.0_wp) cycle
          vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
          metal_valence_luminosity_outside = metal_valence_luminosity_outside &
             + gas_ne(il)*vol*metal_valence_outside(il, gas_Te(il))
       end do
    end if
    if (diffuse_lum > 0.0_wp) then
       diffuse_nphot = max(1_int64, nint(diffuse_lum/lpacket, int64))
       diffuse_packet_lum = diffuse_lum/real(diffuse_nphot, wp)
    else
       diffuse_nphot = 0
       diffuse_packet_lum = lpacket
    end if
  contains

    subroutine diffuse_shared_sync()
      integer :: k
      do k = 1, 3
         if (dif_win(k) /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, dif_win(k), ierr)
      end do
      call MPI_BARRIER(shared_host_comm(), ierr)
    end subroutine diffuse_shared_sync

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
       if (ch == 5) then
          eph = heii_cascade_energy(gas_Te(il), ub, u)
       else if (ch == CH_METAL_GROUND) then
          eph = metal_ground_energy(il, gas_Te(il), ub, u)
       else if (ch == CH_METAL_VALENCE) then
          eph = metal_valence_energy(il, gas_Te(il), ub, u)
       else
          eph = freefree_energy(il, gas_Te(il), ub, u)
       end if
    end if

    !--- bin index by binary search of the band edges (the one source of truth
    !--- for BOTH the aligned and the legacy-log grid; recombination photons all
    !--- carry eph >= the channel threshold > eion_min, so they land in the
    !--- ionizing segment).  The old closed-form log formula assumed a pure
    !--- log grid and mis-binned under par%ion_align_edges.
    photon%inu = ion_bin_of(eph)

    photon%wgt     = 1.0_wp
    photon%Lpacket = diffuse_packet_lum
    photon%source_channel = ch
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
  !          channel 6 the metal transition; channel 7 the ion charge;
  !          channel 8 the metal transition;
  !   u(9)   second energy variate: the He I two-photon energy from the
  !          tabulated Drake shape, the channel-5 free-bound / two-photon
  !          energy, the channel-6 Milne energy or the channel-7 free-free
  !          energy, or the channel-8 line (0.75 or 0.25 of the ionization
  !          potential) (else unused).
  ! Channels 5-8 reuse dimensions 8 and 9 rather than adding any, since
  ! their pick and their energy draw are the same two decisions channel 4
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
    else if (ch == 5) then
       eph = heii_cascade_energy(gas_Te(il), u(8), u(9))
    else if (ch == CH_METAL_GROUND) then
       eph = metal_ground_energy(il, gas_Te(il), u(8), u(9))
    else if (ch == CH_METAL_VALENCE) then
       eph = metal_valence_energy(il, gas_Te(il), u(8), u(9))
    else
       eph = freefree_energy(il, gas_Te(il), u(8), u(9))
    end if

    photon%inu = ion_bin_of(eph)

    photon%wgt     = 1.0_wp
    photon%Lpacket = diffuse_packet_lum
    photon%source_channel = ch
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

  !--- In-band luminosity of each metal ground continuum of leaf il per n_e
  !--- per unit volume [erg cm^-3 s^-1 per cm^-3]: n(X^(i+1)) times the
  !--- band coefficient of transition k.  One routine for the emission table
  !--- and both launchers, so the channel luminosity and the transition pick
  !--- cannot disagree; the stage fractions are the ones the packet opacity
  !--- uses (species_stage_cache_refresh, current gas state).
  function metal_ground_weights(il, T) result(w)
    use species_mod, only : species_transition, species_stage_fraction, elem_abund
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: w(MAX_METAL_TRANSITIONS)
    integer :: k, ie, it
    w = 0.0_wp
    if (.not. par%use_metals) return
    do k = 1, milne_metal_count()
       call species_transition(k, ie, it)
       w(k) = gas_nH(il)*elem_abund(ie)*species_stage_fraction(ie, it+1, il) &
              *milne_metal_band_coefficient(k, T)
    end do
  end function metal_ground_weights

  !--- In-band free-free luminosity of leaf il by net ion charge Z = 1..
  !--- FF_ZMAX, per n_e per unit volume: Z = 1 H+ and He+ (and the metal
  !--- stages II), Z = 2 He2+ (and stages III), Z >= 3 the metal stages
  !--- Z+1, the metals with par%metal_freefree (the switch of their cooling).
  function freefree_weights(il, T) result(w)
    use species_mod, only : n_elements, elem_nstage, elem_abund, &
                            species_stage_fraction
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: w(FF_ZMAX), nz(FF_ZMAX), nH
    integer :: ie, j, z
    nH = gas_nH(il)
    nz = 0.0_wp
    nz(1) = nH*(gas_xHII(il) + par%abund_He*gas_xHeII(il))
    nz(2) = nH*par%abund_He*gas_xHeIII(il)
    if (par%use_metals .and. par%metal_freefree) then
       do ie = 1, n_elements
          do j = 2, elem_nstage(ie)
             z = j - 1
             nz(z) = nz(z) + nH*elem_abund(ie)*species_stage_fraction(ie, j, il)
          end do
       end do
    end if
    do z = 1, FF_ZMAX
       w(z) = nz(z)*freefree_band_coefficient(z, T)
    end do
  end function freefree_weights

  !--- Photon energy [eV] of one channel-6 packet: the transition picked by
  !--- u_pick over metal_ground_weights, the energy from its Milne CDF.
  real(kind=wp) function metal_ground_energy(il, T_K, u_pick, u_energy) result(eph)
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T_K, u_pick, u_energy
    real(kind=wp) :: w(MAX_METAL_TRANSITIONS)
    integer :: k
    w = metal_ground_weights(il, T_K)
    k = pick_index(w(1:milne_metal_count()), u_pick)
    eph = milne_metal_sample_energy(k, T_K, u_energy)
  end function metal_ground_energy

  !--- Photon energy [eV] of one channel-7 packet: the charge picked by
  !--- u_pick over freefree_weights, the energy from its free-free CDF.
  real(kind=wp) function freefree_energy(il, T_K, u_pick, u_energy) result(eph)
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T_K, u_pick, u_energy
    integer :: z
    z = pick_index(freefree_weights(il, T_K), u_pick)
    eph = freefree_band_sample_energy(z, T_K, u_energy)
  end function freefree_energy

  !--- the first entry of w whose cumulative weight passes u (entries with
  !--- no weight skipped, as in pick_channel)
  integer function pick_index(w, u) result(k)
    real(kind=wp), intent(in) :: w(:), u
    real(kind=wp) :: acc, tot
    integer :: i
    tot = sum(w)
    acc = 0.0_wp
    k   = 1
    do i = 1, size(w)
       if (w(i) <= 0.0_wp) cycle
       acc = acc + w(i)
       k   = i
       if (u*tot <= acc) return
    end do
  end function pick_index

  !--- In-band channel-8 luminosity of each metal transition of leaf il per
  !--- n_e per unit volume [erg cm^-3 s^-1 per cm^-3]: n(X^(i+1)) times the
  !--- total (RR + DR) coefficient times the energy of the two photons, each
  !--- counted where it falls in [eion_min, eion_max).  X^i with fewer than
  !--- 3 electrons is left out, as in Cloudy (its iso-sequence ions).  One
  !--- routine for the emission table and both launchers.
  function metal_valence_weights(il, T) result(w)
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: w(MAX_METAL_TRANSITIONS), e(2), rate
    integer :: k, ie, it
    w = 0.0_wp
    do k = 1, metal_valence_count()
       if (.not. metal_valence_recombination(k, il, T, ie, it, rate, e)) cycle
       w(k) = rate*ev2erg*sum(e, mask=valence_line_in_band(e))
    end do
  end function metal_valence_weights

  !--- the same luminosity below eion_min (1) and at or above eion_max (2)
  function metal_valence_outside(il, T) result(lout)
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lout(2), e(2), rate
    integer :: k, ie, it
    lout = 0.0_wp
    do k = 1, metal_valence_count()
       if (.not. metal_valence_recombination(k, il, T, ie, it, rate, e)) cycle
       lout(1) = lout(1) + rate*ev2erg*sum(e, mask=e < par%eion_min)
       lout(2) = lout(2) + rate*ev2erg*sum(e, mask=e >= par%eion_max)
    end do
  end function metal_valence_outside

  integer function metal_valence_count() result(n)
    use species_mod, only : species_transition_count
    n = 0
    if (par%use_metals) n = species_transition_count()
  end function metal_valence_count

  !--- transition k: its recombination rate per n_e [s^-1 per cm^-3, i.e.
  !--- n(X^(i+1)) alpha_RR+DR] and its two photon energies [eV]; .false. when
  !--- the recombined ion has fewer than 3 electrons
  logical function metal_valence_recombination(k, il, T, ie, it, rate, e) result(ok)
    use species_mod, only : species_transition, species_stage_fraction, &
                            species_recombination_coefficient, elem_abund, &
                            elem_eth, elem_atomic_number
    integer,       intent(in)  :: k, il
    real(kind=wp), intent(in)  :: T
    integer,       intent(out) :: ie, it
    real(kind=wp), intent(out) :: rate, e(2)
    call species_transition(k, ie, it)
    rate = 0.0_wp;  e = 0.0_wp
    !--- stage it (X^(it-1)+) keeps Z - it + 1 electrons
    ok = elem_atomic_number(ie) - it + 1 >= 3
    if (.not. ok) return
    rate = gas_nH(il)*elem_abund(ie)*species_stage_fraction(ie, it+1, il) &
           *species_recombination_coefficient(ie, it, T)
    e = VALENCE_FRACTION*elem_eth(ie, it)
  end function metal_valence_recombination

  elemental logical function valence_line_in_band(e) result(inb)
    real(kind=wp), intent(in) :: e
    inb = e >= par%eion_min .and. e < par%eion_max
  end function valence_line_in_band

  !--- Photon energy [eV] of one channel-8 packet: the transition picked by
  !--- u_pick over metal_valence_weights, then of its in-band lines the one
  !--- picked by u_line over their energy luminosities (equal photon rates,
  !--- so in proportion to the energies).
  real(kind=wp) function metal_valence_energy(il, T_K, u_pick, u_line) result(eph)
    use utility, only : fatal_error
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T_K, u_pick, u_line
    real(kind=wp) :: w(MAX_METAL_TRANSITIONS), e(2), rate
    logical :: inb(2)
    integer :: k, ie, it
    w = metal_valence_weights(il, T_K)
    k = pick_index(w(1:metal_valence_count()), u_pick)
    if (.not. metal_valence_recombination(k, il, T_K, ie, it, rate, e)) &
       call fatal_error('channel 8 picked a transition without photons')
    inb = valence_line_in_band(e)
    if (inb(1) .and. inb(2)) then
       eph = merge(e(1), e(2), u_line*(e(1) + e(2)) <= e(1))
    else if (inb(1)) then
       eph = e(1)
    else
       eph = e(2)
    end if
  end function metal_valence_energy

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
