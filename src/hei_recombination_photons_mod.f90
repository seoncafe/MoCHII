module hei_recombination_photons_mod
!---------------------------------------------------------------------------
! MoCHII: the hydrogen-ionizing photons of He+ -> He0 recombination, for the
! two ways the code can dispose of them:
!  - transported, as the He I channels of the diffuse field (diffuse_mod:
!    the ground continuum as channel 2, the excited cascade as channel 4,
!    par%hei_diffuse), or
!  - absorbed where they are made, the on-the-spot approximation
!    (par%recombination_continuum = 'on_the_spot'), with par%hei_on_the_spot.
! One routine gives the decay branches to both, so the transported field and
! the on-the-spot sources cannot count the cascade differently.
!
! EXCITED CASCADE.  Every case-B recombination (alpha_B) reaches n = 2 and
! leaves through one of three H-ionizing exits, all below the He0 edge
! (24.587 eV), so on the spot only hydrogen can take them:
!   2^3S -> 1^1S          19.82 eV, one photon;
!   2^1P -> 1^1S          21.22 eV (584 A), trapped; a fraction p_584
!                         (Wood, Mathis & Ercolano 2004, hei_584_ionizes_H)
!                         photoionizes H, the rest is transferred to 2^1S;
!   2^1S -> 1^1S          two photons, hei_2ph_photons_per_decay of them
!                         above 13.598 eV (0.556 per decay).
! Branch fractions f_3S : f_1S : f_1P from hei_cascade_mod.  With the 2^3S
! population coupled (par%hei_metastable_coupling) the triplet exit is what
! the 2^3S atoms radiate, n_3 A_31, and the collisional transfers feed the
! singlet exits (hei_excited_decays_coupled).
!
! GROUND CAPTURES.  In case B the He0 ground-capture photons (alpha_1 =
! alpha_A - alpha_B, E > 24.587 eV) are absorbed on the spot by H0 or He0 in
! proportion to their opacities (Osterbrock & Ferland 2006, sec. 2.4, the
! fraction y).  The share absorbed by He0 re-ionizes it (the case-B
! assumption); the share y absorbed by H0 ionizes hydrogen, so helium's
! effective recombination is alpha_B + y alpha_1 and hydrogen gains
! y alpha_1 n_e n(He+) ionizations.  y is averaged over the Milne
! ground-continuum spectrum E^2 sigma_HeI(E) exp[-(E - E_th)/kT] with a
! four-point Gauss-Laguerre rule in (E - E_th)/kT, on the transport cross
! sections (photo_xsec).
!
! HEAT.  Each photon absorbed by H0 releases E - 13.598 eV to the gas.  For
! the ground captures absorbed by H the capture itself removed the electron's
! kinetic energy, which the case-B recombination cooling (beta from alpha_B)
! does not count, and the photon carries it back: the two cancel exactly and
! the net gain is IP(He) - IP(H) = 10.989 eV per photon.  The share absorbed
! by He0 exchanges nothing, which is what the case-B cooling assumes.
!
! VALIDITY.  On-the-spot absorption by H and He only: dust and the metals,
! which also absorb 13.6-25 eV photons, are not given a share; nor is the
! He 2^3S population (n_3/n(He0) < 1e-5 in H II regions).
!---------------------------------------------------------------------------
  use define
  use mathlib, only : GL_X => GAUSS_LAGUERRE4_X, GL_W => GAUSS_LAGUERRE4_W
  implicit none
  private

  public :: hei_excited_decays, hei_excited_decays_coupled
  public :: hei_ground_nodes_type, hei_ground_nodes, hei_ground_share_H
  public :: hei_on_the_spot_sources, hei_on_the_spot_heating
  public :: hei_twophoton_decays
  public :: HEI_E3S, HEI_E1P

  !--- line energies of the two single-photon exits [eV]
  real(kind=wp), parameter :: HEI_E3S = 19.82_wp, HEI_E1P = 21.22_wp
  !--- IP(He) - IP(H) [eV] (NIST ASD 24.587387 - 13.598434)
  real(kind=wp), parameter :: DE_HE_H = 10.988953_wp


  !--- the Milne-weighted cross sections of the He0 ground continuum at T
  type :: hei_ground_nodes_type
     real(kind=wp) :: sH(4) = 0.0_wp, sHe(4) = 0.0_wp, w(4) = 0.0_wp
  end type hei_ground_nodes_type

contains

  !=========================================================================
  ! Exits of one case-B recombination, uncoupled 2^3S: r3s 19.82 eV photons,
  ! r1p 584 A photons that ionize H, r2q 2^1S decays (including the
  ! transferred 584 A).  Dimensionless, per recombination.
  !=========================================================================
  subroutine hei_excited_decays(T_K, xHI, xHeI, r3s, r1p, r2q)
    use hei_cascade_mod, only : hei_branch_fractions, hei_584_ionizes_H
    real(kind=wp), intent(in)  :: T_K, xHI, xHeI
    real(kind=wp), intent(out) :: r3s, r1p, r2q
    real(kind=wp) :: f3s, f1s, f1p, p584
    call hei_branch_fractions(T_K, f3s, f1s, f1p)
    !--- The trapped 584 A photon either photoionizes hydrogen or ends in
    !--- 2^1S (the 2.06 micron decay of 2^1P) and leaves as the two-photon
    !--- continuum instead.  Both outcomes are local.
    p584 = hei_584_ionizes_H(T_K, xHI, xHeI)
    r3s = f3s
    r1p = f1p*p584
    r2q = f1s + f1p*(1.0_wp - p584)
  end subroutine hei_excited_decays

  !=========================================================================
  ! The same exits with the 2^3S population coupled
  ! (par%hei_metastable_coupling), per unit volume [cm^-3 s^-1].  The
  ! 19.82 eV photons are what the 2^3S atoms actually radiate, n_3 A_31,
  ! not every triplet recombination: the rest of them leave 2^3S by
  ! photo- or electron-impact ionization, Penning ionization or
  ! super-elastic collisions (no photon) or by collisional transfer into
  ! the singlets.  The transfers add to the singlet exits: n_e n_3 q_31a to
  ! 2^1S, n_e n_3 q_31b to 2^1P, and the transfers to n^1L, n >= 3, are
  ! divided between 2^1S and 2^1P in the ratio f_1S : f_1P of the singlet
  ! recombination branches -- an approximation: it gives a collisionally
  ! populated n^1L the cascade of a recombining one.  In the limit n_e -> 0
  ! with Phi_3 = 0, n_HI = 0 and no collisional excitation from 1^1S the
  ! three reduce to n_e n_HeII alpha_B times hei_excited_decays (gate
  ! tests/hei_metastable_coupling).
  !=========================================================================
  subroutine hei_excited_decays_coupled(il, T_K, ne, nHI, nHeI, nHeII, &
                                        xHI, xHeI, r3s, r1p, r2q)
    use recomb_mod,         only : alphaB_HeII
    use hei_cascade_mod,    only : hei_584_ionizes_H
    use hei_metastable_mod, only : hei_triplet_rates_type, hei_triplet_rates, &
                                   hei_triplet_population, hei_phi3_rate, &
                                   hm_A31
    integer,       intent(in)  :: il
    real(kind=wp), intent(in)  :: T_K, ne, nHI, nHeI, nHeII, xHI, xHeI
    real(kind=wp), intent(out) :: r3s, r1p, r2q
    type(hei_triplet_rates_type) :: hr
    real(kind=wp) :: n3, rec, sing, n1s, n1p, p584
    hr   = hei_triplet_rates(T_K)
    n3   = hei_triplet_population(hr, ne, nHI, nHeI, nHeII, hei_phi3_rate(il))
    rec  = ne*nHeII*alphaB_HeII(T_K)
    sing = ne*n3*hr%q31s/(hr%f1s + hr%f1p)
    n1s  = rec*hr%f1s + ne*n3*hr%q31a + sing*hr%f1s
    n1p  = rec*hr%f1p + ne*n3*hr%q31b + sing*hr%f1p
    p584 = hei_584_ionizes_H(T_K, xHI, xHeI)
    r3s  = n3*hm_A31
    r1p  = n1p*p584
    r2q  = n1s + n1p*(1.0_wp - p584)
  end subroutine hei_excited_decays_coupled

  !=========================================================================
  ! H and He0 photoionization cross sections at the four Gauss-Laguerre
  ! nodes of the He0 Milne spectrum at T, with the spectral weights
  ! E^2 sigma_HeI(E) of each node.  T is fixed inside a cell solve, so this
  ! is evaluated once there.
  !=========================================================================
  type(hei_ground_nodes_type) function hei_ground_nodes(T) result(g)
    use photo_xsec, only : sigma_HI, sigma_HeI
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: kT, E
    integer :: k
    kT = kboltz_cgs*T/ev2erg
    do k = 1, 4
       E = eth_HeI + GL_X(k)*kT
       g%sH(k)  = sigma_HI(E)
       g%sHe(k) = sigma_HeI(E)
       g%w(k)   = GL_W(k)*E*E*g%sHe(k)
    end do
    g%w = g%w/sum(g%w)
  end function hei_ground_nodes

  !--- the share y of the He0 ground-capture photons that H0 absorbs
  real(kind=wp) function hei_ground_share_H(g, nHI, nHeI) result(y)
    type(hei_ground_nodes_type), intent(in) :: g
    real(kind=wp),               intent(in) :: nHI, nHeI
    real(kind=wp) :: tH, tHe
    integer :: k
    y = 0.0_wp
    do k = 1, 4
       tH  = max(nHI, 0.0_wp)*g%sH(k)
       tHe = max(nHeI, 0.0_wp)*g%sHe(k)
       if (tH + tHe > 0.0_wp) y = y + g%w(k)*tH/(tH + tHe)
    end do
  end function hei_ground_share_H

  !=========================================================================
  ! The on-the-spot hydrogen ionizations [cm^-3 s^-1] of one cell's He+
  ! recombinations, at the state (n_H, n_e, x_HI, x_HeI, x_HeII): s_gnd from
  ! the ground captures H0 takes (y alpha_1 n_e n_HeII), s_exc from the
  ! excited cascade.  alpha_1 = alpha_A - alpha_B of the recombination model
  ! in use.
  !=========================================================================
  subroutine hei_on_the_spot_sources(il, T, nH, ne, xHI, xHeI, xHeII, g, &
                                     y, s_gnd, s_exc)
    use recomb_mod,        only : alphaA_HeII, alphaB_HeII
    use hei_twophoton_mod, only : hei_2ph_photons_per_decay
    integer,                     intent(in)  :: il
    real(kind=wp),               intent(in)  :: T, nH, ne, xHI, xHeI, xHeII
    type(hei_ground_nodes_type), intent(in)  :: g
    real(kind=wp),               intent(out) :: y, s_gnd, s_exc
    real(kind=wp) :: r3s, r1p, r2q, aB, nHeII
    nHeII = par%abund_He*nH*xHeII
    y     = hei_ground_share_H(g, nH*xHI, par%abund_He*nH*xHeI)
    aB    = alphaB_HeII(T)
    s_gnd = y*(alphaA_HeII(T) - aB)*ne*nHeII
    call excited_exits(il, T, nH, ne, xHI, xHeI, xHeII, aB, r3s, r1p, r2q)
    s_exc = r3s + r1p + r2q*hei_2ph_photons_per_decay()
  end subroutine hei_on_the_spot_sources

  !=========================================================================
  ! Gas heating [erg cm^-3 s^-1] of the same absorptions: 10.989 eV for each
  ! ground-capture photon H0 takes (header), E - 13.598 eV for each excited
  ! photon, and for the two-photon exit the in-band energy per decay less
  ! 13.598 eV for each H-ionizing photon of the pair.
  !=========================================================================
  !--- f_exc: the fraction of the photons given to H0 that it can absorb
  !--- (solve_ion_cell's f_h_only; below 1 only where the helium cascades
  !--- outnumber hydrogen's recombinations, heii_recombination_photons_mod
  !--- "Absorption limit").
  real(kind=wp) function hei_on_the_spot_heating(il, T, nH, ne, xHI, xHeI, &
                                                 xHeII, g, f_exc) result(heat)
    use recomb_mod,        only : alphaA_HeII, alphaB_HeII
    use hei_twophoton_mod, only : hei_2ph_photons_per_decay, &
                                  hei_2ph_energy_per_decay
    integer,                     intent(in) :: il
    real(kind=wp),               intent(in) :: T, nH, ne, xHI, xHeI, xHeII
    type(hei_ground_nodes_type), intent(in) :: g
    real(kind=wp), optional,     intent(in) :: f_exc
    real(kind=wp) :: y, aB, r3s, r1p, r2q, nHeII, f
    nHeII = par%abund_He*nH*xHeII
    y  = hei_ground_share_H(g, nH*xHI, par%abund_He*nH*xHeI)
    aB = alphaB_HeII(T)
    call excited_exits(il, T, nH, ne, xHI, xHeI, xHeII, aB, r3s, r1p, r2q)
    f = 1.0_wp
    if (present(f_exc)) f = f_exc
    heat = f*( y*(alphaA_HeII(T) - aB)*ne*nHeII*DE_HE_H           &
             + r3s*(HEI_E3S - eth_HI) + r1p*(HEI_E1P - eth_HI)    &
             + r2q*(hei_2ph_energy_per_decay()                    &
                    - hei_2ph_photons_per_decay()*eth_HI) )*ev2erg
  end function hei_on_the_spot_heating

  !=========================================================================
  ! 2^1S -> 1^1S decays per unit volume [cm^-3 s^-1] (direct recombination,
  ! the 584 A photons transferred to 2^1S and, with the 2^3S population
  ! coupled, its collisional transfers): the rate the written He I
  ! two-photon continuum radiates at (nebcont_mod), from the same exits the
  ! diffuse channel 4 and the on-the-spot sources take.
  !=========================================================================
  real(kind=wp) function hei_twophoton_decays(il, T, nH, ne, xHI, xHeI, &
                                              xHeII) result(r2q)
    use recomb_mod, only : alphaB_HeII
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne, xHI, xHeI, xHeII
    real(kind=wp) :: r3s, r1p
    call excited_exits(il, T, nH, ne, xHI, xHeI, xHeII, alphaB_HeII(T), &
                       r3s, r1p, r2q)
  end function hei_twophoton_decays

  !--- the three exits per unit volume [cm^-3 s^-1], coupled or not
  subroutine excited_exits(il, T, nH, ne, xHI, xHeI, xHeII, aB, r3s, r1p, r2q)
    integer,       intent(in)  :: il
    real(kind=wp), intent(in)  :: T, nH, ne, xHI, xHeI, xHeII, aB
    real(kind=wp), intent(out) :: r3s, r1p, r2q
    real(kind=wp) :: rec, nHe
    nHe = par%abund_He*nH
    if (par%hei_metastable_coupling) then
       call hei_excited_decays_coupled(il, T, ne, nH*xHI, nHe*xHeI, &
            nHe*xHeII, xHI, xHeI, r3s, r1p, r2q)
    else
       rec = ne*nHe*xHeII*aB
       call hei_excited_decays(T, xHI, xHeI, r3s, r1p, r2q)
       r3s = rec*r3s;  r1p = rec*r1p;  r2q = rec*r2q
    end if
  end subroutine excited_exits

end module hei_recombination_photons_mod
