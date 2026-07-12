module diffuse_mod
!---------------------------------------------------------------------------
! MoCHII: explicit diffuse ionizing field (Stage G3).
!
! Case B (G1/G2) absorbs ground-level recombination photons on the spot.
! G3 replaces that: the equilibrium runs with CASE A rates and every
! recombination to the ground level emits an ionizing packet, transported
! exactly like a stellar packet (docs/PLAN.md section 2.8 — structurally
! the dust re-emission of lucy_mod).  Three continua:
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
  use octree_mod,    only : amr_grid
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, &
                            gas_Te, gas_nleaf
  use recomb_mod
  implicit none
  private

  public :: diffuse_build, gen_diffuse_photon, diffuse_nphot, diffuse_lum

  real(kind=wp), allocatable :: dif_cdf(:)      ! leaf CDF (total energy)
  real(kind=wp), allocatable :: dif_ch(:,:)     ! (3, nleaf) channel luminosities
  real(kind=wp) :: diffuse_lum = 0.0_wp         ! total [erg/s]
  integer(kind=int64) :: diffuse_nphot = 0

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
       allocate(dif_cdf(gas_nleaf), dif_ch(3, gas_nleaf))
    end if

    do il = 1, gas_nleaf
       nH = gas_nH(il)
       if (nH <= 0.0_wp) then
          dif_ch(:,il) = 0.0_wp
          cycle
       end if
       T  = gas_Te(il)
       ne = gas_ne(il)
       vol = (2.0_wp*amr_grid%ch(amr_grid%icell_of_leaf(il)) &
             * par%distance2cm)**3
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
    use ion_band_mod, only : ion_e, ion_de, ion_Ltot
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp) :: u, cost, sint, phi, half, eph, kT_eV, w(3)
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
    ic = amr_grid%icell_of_leaf(il)
    half = amr_grid%ch(ic)

    !--- uniform position within the leaf
    photon%x = amr_grid%cx(ic) + (2.0_wp*rand_number() - 1.0_wp)*half
    photon%y = amr_grid%cy(ic) + (2.0_wp*rand_number() - 1.0_wp)*half
    photon%z = amr_grid%cz(ic) + (2.0_wp*rand_number() - 1.0_wp)*half

    cost = 2.0_wp*rand_number() - 1.0_wp
    sint = sqrt(1.0_wp - cost*cost)
    phi  = twopi*rand_number()
    photon%kx = sint*cos(phi)
    photon%ky = sint*sin(phi)
    photon%kz = cost

    !--- channel, then photon energy E = E_th + kT x, x ~ Exp(1)
    w = dif_ch(:,il)
    u = rand_number()*sum(w)
    if (u <= w(1)) then
       ch = 1;  eph = eth_HI
    else if (u <= w(1) + w(2)) then
       ch = 2;  eph = eth_HeI
    else
       ch = 3;  eph = eth_HeII
    end if
    kT_eV = kboltz_cgs*gas_Te(il)/ev2erg
    eph = eph + kT_eV*(-log(max(rand_number(), tinest)))
    eph = min(eph, par%eion_max*0.999_wp)

    !--- bin index (log-uniform IONIZING segment: recombination photons
    !--- all carry eph >= the channel threshold > eion_min, so with the
    !--- FUV extension the first par%nnu_fuv bins are simply offset).
    block
      integer :: nfuv, nion
      nfuv = merge(par%nnu_fuv, 0, par%add_fuv)
      nion = par%nnu_ion - nfuv
      inu  = nfuv + int(real(nion,wp)*log(eph/par%eion_min) &
             / log(par%eion_max/par%eion_min)) + 1
      photon%inu = min(max(inu, nfuv+1), par%nnu_ion)
    end block

    photon%wgt     = 1.0_wp
    photon%Lpacket = ion_Ltot/real(par%nphotons, wp)
    photon%nscatt  = 0
    photon%inside  = .true.
    photon%icell_amr = il
  end subroutine gen_diffuse_photon

end module diffuse_mod
