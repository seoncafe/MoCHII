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
  implicit none
  private

  public :: diffuse_build, gen_diffuse_photon, diffuse_nphot, diffuse_lum

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
                              HEI_P2PH = 0.56_wp, HEI_E2PH = 20.62_wp
  !--- mean H-ionizing energy of the 2^1S branch (flat in [13.6, 20.62])
  real(kind=wp), parameter :: HEI_E1S = 0.5_wp*(13.598_wp + HEI_E2PH)

contains

  !=========================================================================
  subroutine diffuse_build(lpacket)
    implicit none
    real(kind=wp), intent(in) :: lpacket
    real(kind=wp) :: T, ne, nH, vol, kT_eV, a1
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
       dif_ch(1,il) = ne*nHII   *a1*(eth_HI  + kT_eV)*ev2erg*vol
       a1 = alphaA_HeII(T)  - alphaB_HeII(T)
       dif_ch(2,il) = ne*nHeII_n*a1*(eth_HeI + kT_eV)*ev2erg*vol
       a1 = alphaA_HeIII(T) - alphaB_HeIII(T)
       dif_ch(3,il) = ne*nHeIII *a1*(eth_HeII+ kT_eV)*ev2erg*vol
       !--- channel 4: He I EXCITED-level recombination radiation
       !--- (par%hei_diffuse): all case-B He II recombinations cascade to
       !--- n = 2 and decay with the AGN3 branching; only the H-ionizing
       !--- part carries band luminosity.
       dif_ch(4,il) = 0.0_wp
       if (par%hei_diffuse) then
          a1 = alphaB_HeII(T)
          dif_ch(4,il) = ne*nHeII_n*a1*vol*ev2erg* &
             ( HEI_F3S*HEI_E3S + HEI_F1P*HEI_E1P &
             + HEI_F1S*HEI_P2PH*HEI_E1S )
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
    use ion_band_mod, only : ion_e, ion_de, ion_Ltot, nnu_band
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp) :: u, cost, sint, phi, half, eph, kT_eV, w(4), ub
    integer :: lo, hi, mid, il, ic, ch, inu

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
    !--- draws flat in [13.6, 20.62] eV).
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
       eph = eph + kT_eV*(-log(max(rand_number(), tinest)))
       eph = min(eph, par%eion_max*0.999_wp)
    else
       ub = rand_number()*( HEI_F3S*HEI_E3S + HEI_F1P*HEI_E1P &
                          + HEI_F1S*HEI_P2PH*HEI_E1S )
       if (ub <= HEI_F3S*HEI_E3S) then
          eph = HEI_E3S
       else if (ub <= HEI_F3S*HEI_E3S + HEI_F1P*HEI_E1P) then
          eph = HEI_E1P
       else
          eph = eth_HI + rand_number()*(HEI_E2PH - eth_HI)
       end if
    end if

    !--- bin index (log-uniform IONIZING segment: recombination photons
    !--- all carry eph >= the channel threshold > eion_min, so with the
    !--- FUV extension the first par%nnu_fuv bins are simply offset).
    block
      integer :: nfuv, nion
      nfuv = merge(par%nnu_fuv, 0, par%add_fuv)
      nion = par%nnu_ion
      inu  = nfuv + int(real(nion,wp)*log(eph/par%eion_min) &
             / log(par%eion_max/par%eion_min)) + 1
      photon%inu = min(max(inu, nfuv+1), nnu_band)
    end block

    photon%wgt     = 1.0_wp
    photon%Lpacket = ion_Ltot/real(par%nphotons, wp)
    photon%nscatt  = 0
    photon%inside  = .true.
    photon%from_external = .false.      ! diffuse packets peel isotropically
    photon%icell_amr = il
  end subroutine gen_diffuse_photon

end module diffuse_mod
