module diffuse_mod
!---------------------------------------------------------------------------
! MoCHII: explicit diffuse ionizing field.
!
! Case B absorbs ground-level recombination photons on the spot.
! The explicit field replaces that: the equilibrium runs with CASE A rates
! and every recombination to the ground level emits an ionizing packet,
! transported exactly like a stellar packet.  Three continua:
!     H II  + e -> H I(1s)   : E >= 13.598 eV
!     He II + e -> He I(1^1S): E >= 24.587 eV
!     He III+ e -> He II(1s) : E >= 54.416 eV
! with rates alpha_1 = alpha_A - alpha_B per channel and photon energies
! sampled as E = E_th + kT_e * x, x ~ Exp(1) (the near-threshold Milne
! shape to the accuracy of the band bins).
!
! Each iteration the emission table is rebuilt from the current state:
! L_ch(leaf) = n_e n_recomb alpha_1(T) <E_ch> V_cm3, <E_ch> = E_th + kT.
! Packets carry the same Lpacket as stellar packets, so the diffuse packet
! count is n_diff = L_diffuse / Lpacket (rounded); sampling: leaf CDF over
! total leaf luminosity, then a channel CDF within the leaf.
!---------------------------------------------------------------------------
  use define
  use octree_mod,    only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, &
                            gas_Te, gas_nleaf
  use recomb_mod
  use milne_recomb_spectrum_mod, only : milne_setup, milne_sample_energy, &
                                        milne_energy_per_recomb, milne_is_ready
  use hei_twophoton_mod, only : hei_2ph_sample_energy, &
                               hei_2ph_energy_per_decay
  use hei_cascade_mod,   only : hei_branch_fractions, hei_584_ionizes_H
  implicit none
  private

  public :: diffuse_build, gen_diffuse_photon, gen_diffuse_photon_qmc, &
            diffuse_nphot, diffuse_lum

  real(kind=wp), allocatable :: dif_cdf(:)      ! leaf CDF (total energy)
  real(kind=wp), allocatable :: dif_ch(:,:)     ! (4, nleaf) channel luminosities
  real(kind=wp) :: diffuse_lum = 0.0_wp         ! total [erg/s]
  integer(kind=int64) :: diffuse_nphot = 0

  !--- He I excited-channel branching (AGN3 low-density limit; channel 4,
  !--- par%hei_diffuse): 2^3S single 19.82-eV photon; 2^1P 584 A; 2^1S
  !--- two-photon with probability 0.56 of an H-ionizing photon.
  real(kind=wp), parameter :: HEI_F3S = 0.75_wp,  HEI_E3S = 19.82_wp
  real(kind=wp), parameter :: HEI_F1P = 0.25_wp*(2.0_wp/3.0_wp), &
                              HEI_E1P = 21.22_wp
  real(kind=wp), parameter :: HEI_F1S = 0.25_wp*(1.0_wp/3.0_wp), &
                              HEI_E2PH = 20.62_wp
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
    real(kind=wp) :: nHII, nHeII_n, nHeIII, xHeIII
    integer :: il, ic

    if (allocated(dif_cdf)) then
       if (size(dif_cdf) /= gas_nleaf) deallocate(dif_cdf, dif_ch)
    end if
    if (.not. allocated(dif_cdf)) then
       allocate(dif_cdf(gas_nleaf), dif_ch(4, gas_nleaf))
    end if

    do il = 1, gas_nleaf
       nH = gas_nH(il)
       if (nH <= 0.0_wp) then
          dif_ch(:,il) = 0.0_wp
          cycle
       end if
       T  = gas_Te(il)
       ne = gas_ne(il)
       vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
       kT_eV  = kboltz_cgs*T/ev2erg
       xHeIII = max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il))
       nHII    = nH*(1.0_wp - gas_xHI(il))
       nHeII_n = nH*par%He_abund*gas_xHeII(il)
       nHeIII  = nH*par%He_abund*xHeIII
       !--- channel energy luminosities [erg/s]
       a1 = alphaA_HII(T)   - alphaB_HII(T)
       dif_ch(1,il) = ne*nHII   *a1*ground_mean_energy(1, eth_HI, kT_eV, T)*ev2erg*vol
       a1 = alphaA_HeII(T)  - alphaB_HeII(T)
       dif_ch(2,il) = ne*nHeII_n*a1*ground_mean_energy(2, eth_HeI, kT_eV, T)*ev2erg*vol
       a1 = alphaA_HeIII(T) - alphaB_HeIII(T)
       dif_ch(3,il) = ne*nHeIII *a1*ground_mean_energy(3, eth_HeII, kT_eV, T)*ev2erg*vol
       !--- channel 4: He I EXCITED-level recombination radiation
       !--- (par%hei_diffuse): all case-B He II recombinations cascade to
       !--- n = 2 and decay with the AGN3 branching; only the H-ionizing
       !--- part carries band luminosity.
       dif_ch(4,il) = 0.0_wp
       if (par%hei_diffuse) then
          a1 = alphaB_HeII(T)
          call hei_branch_energies(T, gas_xHI(il), gas_xHeI(il), &
                                   be3s, be1p, be1s)
          dif_ch(4,il) = ne*nHeII_n*a1*vol*ev2erg*(be3s + be1p + be1s)
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
       diffuse_nphot = nint(diffuse_lum/lpacket, int64)
    else
       diffuse_nphot = 0
    end if
  end subroutine diffuse_build

  !=========================================================================
  subroutine gen_diffuse_photon(photon)
    use random,       only : rand_number
    use ion_band_mod, only : ion_Ltot, ion_bin_of
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp) :: u, cost, sint, phi, half, eph, kT_eV, w(4), ub
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

    !--- channel; ground continua sample E = E_th + kT x, x ~ Exp(1);
    !--- the He I excited channel (4) samples its decay branch by
    !--- energy luminosity (fixed line energies; the two-photon branch
    !--- draws from the Drake et al. 1969 shape, hei_2ph_sample_energy).
    w = dif_ch(:,il)
    u = rand_number()*sum(w)
    if (u <= w(1)) then
       ch = 1;  eph = eth_HI
    else if (u <= w(1) + w(2)) then
       ch = 2;  eph = eth_HeI
    else if (u <= w(1) + w(2) + w(3)) then
       ch = 3;  eph = eth_HeII
    else
       ch = 4
    end if
    if (ch <= 3) then
       kT_eV = kboltz_cgs*gas_Te(il)/ev2erg
       eph = ground_photon_energy(ch, eph, kT_eV, gas_Te(il), rand_number())
    else
       call hei_branch_energies(gas_Te(il), gas_xHI(il), gas_xHeI(il), &
                                be3s, be1p, be1s)
       ub = rand_number()*(be3s + be1p + be1s)
       if (ub <= be3s) then
          eph = HEI_E3S
       else if (ub <= be3s + be1p) then
          eph = HEI_E1P
       else
          eph = hei_2ph_sample_energy(rand_number())
       end if
    end if

    !--- bin index by binary search of the band edges (the one source of truth
    !--- for BOTH the aligned and the legacy-log grid; recombination photons all
    !--- carry eph >= the channel threshold > eion_min, so they land in the
    !--- ionizing segment).  The old closed-form log formula assumed a pure
    !--- log grid and mis-binned under par%ion_align_edges.
    photon%inu = ion_bin_of(eph)

    photon%wgt     = 1.0_wp
    photon%Lpacket = ion_Ltot/real(par%nphotons, wp)
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
  ! semantics (docs/QUASI_RANDOM_LAUNCH.md), from the diffuse scramble stream:
  !   u(1)   emitting leaf (dif_cdf inverse);
  !   u(2-4) uniform position within the leaf;
  !   u(5)   mu = 2u-1;                 u(6) phi = 2 pi u;
  !   u(7)   emission channel (dif_ch inverse over the 4 channels);
  !   u(8)   first energy variate: channels 1-3 the Exp(1) via -log(u);
  !          channel 4 the He I decay-branch pick;
  !   u(9)   second energy variate: the He I two-photon energy from the
  !          tabulated Drake shape (else unused).
  ! Every expression matches gen_diffuse_photon; the bin comes from ion_bin_of.
  !=========================================================================
  subroutine gen_diffuse_photon_qmc(photon, u)
    use ion_band_mod, only : ion_Ltot, ion_bin_of
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp),     intent(in)    :: u(:)
    real(kind=wp) :: uu, cost, sint, phi, half, eph, kT_eV, w(4), ub
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

    !--- channel; ground continua sample E = E_th + kT x, x ~ Exp(1);
    !--- the He I excited channel (4) samples its decay branch by energy
    !--- luminosity (fixed line energies; the two-photon branch from the
    !--- tabulated Drake shape).
    w = dif_ch(:,il)
    uu = u(7)*sum(w)
    if (uu <= w(1)) then
       ch = 1;  eph = eth_HI
    else if (uu <= w(1) + w(2)) then
       ch = 2;  eph = eth_HeI
    else if (uu <= w(1) + w(2) + w(3)) then
       ch = 3;  eph = eth_HeII
    else
       ch = 4
    end if
    if (ch <= 3) then
       kT_eV = kboltz_cgs*gas_Te(il)/ev2erg
       eph = ground_photon_energy(ch, eph, kT_eV, gas_Te(il), u(8))
    else
       call hei_branch_energies(gas_Te(il), gas_xHI(il), gas_xHeI(il), &
                                be3s, be1p, be1s)
       ub = u(8)*(be3s + be1p + be1s)
       if (ub <= be3s) then
          eph = HEI_E3S
       else if (ub <= be3s + be1p) then
          eph = HEI_E1P
       else
          eph = hei_2ph_sample_energy(u(9))
       end if
    end if

    photon%inu = ion_bin_of(eph)

    photon%wgt     = 1.0_wp
    photon%Lpacket = ion_Ltot/real(par%nphotons, wp)
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
    real(kind=wp) :: f3s, f1s, f1p, p584
    call hei_branch_fractions(T_K, f3s, f1s, f1p)
    !--- The trapped 584 A photon either photoionizes hydrogen or is
    !--- collisionally transferred 2^1P -> 2^1S and leaves as the two-photon
    !--- continuum instead.  Both outcomes are local, so the split folds into
    !--- the branch weights here and needs no extra random variate; the Sobol
    !--- dimension assignment is unchanged.
    p584 = hei_584_ionizes_H(T_K, xHI, xHeI)
    e3s = f3s*HEI_E3S
    e1p = f1p*p584*HEI_E1P
    !--- A converted decay radiates the 2^1S pair, 20.62 eV, not 21.22 eV;
    !--- the remaining 0.60 eV leaves as 2^1P -> 2^1S at 2.06 micron, which is
    !--- non-ionizing and so correctly absent from the band.
    e1s = (f1s + f1p*(1.0_wp - p584))*hei_2ph_energy_per_decay()
  end subroutine hei_branch_energies

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
