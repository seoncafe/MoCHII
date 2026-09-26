module heii_recombination_photons_mod
!---------------------------------------------------------------------------
! MoCHII: the ionizing photons of He2+ -> He+ recombination absorbed where
! they are made, the on-the-spot approximation
! (par%recombination_continuum = 'on_the_spot' with par%heii_on_the_spot).
! Transported instead, the same photons are diffuse channels 3 (the ground
! continuum) and 5 (the excited cascade, par%heii_diffuse, heii_cascade_mod).
! With par%heii_diffuse set in the on-the-spot mode, the excited cascade is
! transported as channel 5 and only the ground captures are absorbed here
! (heii_ots_nodes): the cascade reaches 40.8-54.4 eV and deposits up to
! ~27 eV a photon, and in a He III zone optically thin at those energies it
! is absorbed far from where it is made (tau ~ 2e-6 across the 0.24 pc core
! of the G2a sphere, which absorbing it in place pinned at te_max).
!
! Osterbrock & Ferland (2006, sec. 2.5) name the coupling: the He2+ zone is
! also an H+ zone, and the He II Ly-alpha (40.8 eV), the 2^2S two-photon
! continuum and the Balmer continuum (whose threshold, 13.604 eV, is the
! hydrogen edge) all ionize H0.  Their Table 2.6 gives 0.64 + 0.36 + 0.20
! H-ionizing photons per H recombination at n(He2+)/n(H+) = 0.15, 1e4 K.
!
! GROUND CAPTURES.  In case B the He2+ -> He+(1s) photons (alpha_1 =
! alpha_A - alpha_B of He III, E >= 54.416 eV, the Milne spectrum
! E^2 sigma_HeII(E) exp[-(E - E_th)/kT]) are absorbed on the spot by He+,
! He0 and H0 in proportion to their opacities at those energies.  The share
! He+ takes re-ionizes it, the case-B assumption, and changes nothing; the
! shares y_H and y_He that H0 and He0 take ionize them instead, so
!   He III's effective recombination is alpha_B + (y_H + y_He) alpha_1,
!   H0 gains  y_H  alpha_1 n_e n(He2+) ionizations,
!   He0 gains y_He alpha_1 n_e n(He2+) ionizations (He0 -> He+).
! The shares are averaged over the Milne spectrum with the four-point
! Gauss-Laguerre rule in (E - E_th)/kT, on the transport cross sections
! (photo_xsec), as the He I ground captures are
! (hei_recombination_photons_mod).
!
! EXCITED CASCADE.  heii_cascade_photons gives, per n_e n(He2+), the photons
! of the three paths channel 5 transports -- free-bound into n >= 2
! (13.598-54.416 eV), 2^2P Ly-alpha (40.813 eV) and the 2^2S pair
! (13.598-40.813 eV in band) -- from the same branch luminosities, split at
! the He0 edge.  None of them reaches 54.416 eV, so He+ absorbs none.
! Photons below 24.587 eV ionize H0 alone; above it, H0 and He0 share each
! part by opacity at the part's photon-mean energy.
!
! HEAT.  Each photon an absorber X takes gives the gas E - IP(X).  For the
! excited photons that is all: the kinetic energy the captures into n >= 2
! removed is in the case-B recombination cooling (beta from alpha_B of
! He III, cooling_mod).  The ground captures are not in it; each one
! removes its electron's kinetic energy KE and emits E = IP(He+) + KE, so
! absorption by X nets IP(He+) - IP(X) for every photon, whatever its
! energy:
!   H0   54.417765 - 13.598434 = 40.819331 eV,
!   He0  54.417765 - 24.587387 = 29.830378 eV,
!   He+  0 (the capture and the re-ionization exchange KE both ways, which
!        is what the case-B cooling assumes)
! (NIST ASD ionization energies).
!
! VALIDITY.  Absorption by H, He0 and He+ only: dust and the metals are not
! given a share, as for He I.  "On the spot" is taken literally: every
! photon is absorbed in the leaf that made it.  Osterbrock & Ferland note
! that Ly-alpha (resonance-scattered by He+) and the Balmer continuum
! (just above the H edge) have short paths there, but that most of the
! two-photon continuum escapes the He2+ zone into the He+ zone; the local
! absorption misplaces those photons, it does not lose them.
!
! ABSORPTION LIMIT.  In steady state H0 cannot absorb more photons than the
! recombinations (and charge exchanges) that make H0.  H0's share of the
! excited cascade does not fall with n(H0) -- below the He0 edge there is
! no other absorber, and above it He0 is as scarce as H0 in the He2+ zone --
! so solve_ion_cell takes it as an absolute source in hydrogen's balance,
! with the share evaluated at every trial x_HI
! (ion_balance_mod::hydrogen_fraction_with_cascade); the ground captures
! stay a rate per H0, since He+ competes for them.  VALIDITY: a local
! balance exists only while the cascades H0 and He0 are given stay below
! the recombinations that make them.  At He/H = 0.1 the H0 source reaches
! 0.92 of hydrogen's recombinations in a weakly illuminated fully He2+
! cell (0.71 volume mean in the 1e5 K blackbody run of
! results/heii_on_the_spot/), and the balance closes to 1e-11; at He/H =
! 0.15 (n(He2+)/n_p = 0.15, 1.2-1.8e4 K) it exceeds 1 and none exists --
! the on-the-spot approximation fails there, the cell settles at x_HI ~ 0,
! the H balance residual printed each iteration reports it, and the heat
! is limited to the photons the supply can absorb (f_H_exc).  Such runs
! (e.g. Type I planetary nebulae) need recombination_continuum =
! 'transported'; setup prints a NOTE when abund_He > 0.12.
!
! Against the transported channels (recombination_continuum = 'transported',
! heii_diffuse) the H front and L(Hbeta) agree, the He2+ zone is hotter;
! there the transported run is no reference either, since it carries
! Ly-alpha as a continuum photon that leaves the zone.
!---------------------------------------------------------------------------
  use define
  use mathlib, only : GL_X => GAUSS_LAGUERRE4_X, GL_W => GAUSS_LAGUERRE4_W
  use heii_cascade_mod, only : HEII_NGROUP
  implicit none
  private

  public :: heii_ots_nodes_type, heii_ots_nodes
  public :: heii_ground_shares, heii_on_the_spot_sources
  public :: heii_on_the_spot_heating, heii_cascade_to_H
  public :: DE_HEII_H, DE_HEII_HEI

  !--- IP(He+) - IP(H) and IP(He+) - IP(He) [eV] (NIST ASD 54.417765,
  !--- 13.598434, 24.587387)
  real(kind=wp), parameter :: DE_HEII_H   = 54.417765_wp - 13.598434_wp
  real(kind=wp), parameter :: DE_HEII_HEI = 54.417765_wp - 24.587387_wp

  !--- the excited-cascade parts, three paths x (HEII_NGROUP + 1) parts
  integer, parameter :: NPART = 3*(HEII_NGROUP + 1)

  !--- everything that depends on T alone, evaluated once per cell solve
  type :: heii_ots_nodes_type
     !--- He III ground-capture coefficient alpha_1 = alpha_A - alpha_B
     real(kind=wp) :: a1 = 0.0_wp
     !--- Milne-weighted cross sections of the He II ground continuum at the
     !--- four Gauss-Laguerre nodes: H0, He0, He+ and the photon weights
     real(kind=wp) :: gH(4) = 0.0_wp, gHe(4) = 0.0_wp, gHe2(4) = 0.0_wp
     real(kind=wp) :: gw(4) = 0.0_wp
     !--- excited cascade: photons per n_e n(He2+) [cm^3 s^-1], mean energy
     !--- [eV] and the H0 and He0 cross sections there, part by part
     real(kind=wp) :: n(NPART) = 0.0_wp, e(NPART) = 0.0_wp
     real(kind=wp) :: sH(NPART) = 0.0_wp, sHe(NPART) = 0.0_wp
  end type heii_ots_nodes_type

contains

  !=========================================================================
  type(heii_ots_nodes_type) function heii_ots_nodes(T) result(g)
    use photo_xsec,       only : sigma_HI, sigma_HeI, sigma_HeII
    use recomb_mod,       only : alphaA_HeIII, alphaB_HeIII
    use heii_cascade_mod, only : heii_cascade_photons_type, heii_cascade_photons, &
                                 heii_cascade_transported
    real(kind=wp), intent(in) :: T
    type(heii_cascade_photons_type) :: c
    real(kind=wp) :: kT, E
    integer :: k
    g%a1 = alphaA_HeIII(T) - alphaB_HeIII(T)
    kT = kboltz_cgs*T/ev2erg
    do k = 1, 4
       E = eth_HeII + GL_X(k)*kT
       g%gH(k)   = sigma_HI(E)
       g%gHe(k)  = sigma_HeI(E)
       g%gHe2(k) = sigma_HeII(E)
       g%gw(k)   = GL_W(k)*E*E*g%gHe2(k)
    end do
    g%gw = g%gw/sum(g%gw)
    c = heii_cascade_photons(T)
    g%n   = reshape(c%n, [NPART])
    !--- with par%heii_diffuse the cascade is channel 5 of the diffuse field
    !--- (heii_cascade_transported), so on the spot only the ground captures
    !--- remain: sources, heat and the absorption limit all go through g%n
    if (heii_cascade_transported()) g%n = 0.0_wp
    g%e   = reshape(c%e, [NPART])
    g%sH  = sigma_HI(g%e)
    g%sHe = sigma_HeI(g%e)
  end function heii_ots_nodes

  !--- the shares of the He II ground-capture photons that H0 and He0 absorb
  !--- (the rest, He+'s share, re-ionizes He+)
  subroutine heii_ground_shares(g, nHI, nHeI, nHeII, yH, yHe)
    type(heii_ots_nodes_type), intent(in)  :: g
    real(kind=wp),             intent(in)  :: nHI, nHeI, nHeII
    real(kind=wp),             intent(out) :: yH, yHe
    real(kind=wp) :: tH, tHe, tot
    integer :: k
    yH = 0.0_wp;  yHe = 0.0_wp
    do k = 1, 4
       tH  = max(nHI, 0.0_wp)*g%gH(k)
       tHe = max(nHeI, 0.0_wp)*g%gHe(k)
       tot = tH + tHe + max(nHeII, 0.0_wp)*g%gHe2(k)
       if (tot > 0.0_wp) then
          yH  = yH  + g%gw(k)*tH/tot
          yHe = yHe + g%gw(k)*tHe/tot
       end if
    end do
  end subroutine heii_ground_shares

  !=========================================================================
  ! The on-the-spot ionizations [cm^-3 s^-1] of one cell's He2+
  ! recombinations: s_H of H0 and s_He of He0, ground captures plus the
  ! excited cascade, with the ground-capture shares yH, yHe.
  !
  ! With s_H_exc present, s_H is the ground captures alone and s_H_exc the
  ! excited cascade H0 takes.  In the He2+ zone He0 is as scarce as H0, so
  ! H0 takes nearly all of that cascade (all of it below the He0 edge)
  ! whatever n(H0) is: as a rate per H0 it would diverge as x_HI -> 0, and
  ! where it outnumbers hydrogen's recombinations no positive x_HI
  ! balances it.  solve_ion_cell therefore takes it as an absolute source
  ! (header, "Absorption limit").  The ground captures stay a rate per H0:
  ! He+ competes for them, so H0's share falls with n(H0).
  !=========================================================================
  subroutine heii_on_the_spot_sources(g, ne, nHI, nHeI, nHeII, nHeIII, &
                                      yH, yHe, s_H, s_He, s_H_exc)
    type(heii_ots_nodes_type), intent(in)  :: g
    real(kind=wp),             intent(in)  :: ne, nHI, nHeI, nHeII, nHeIII
    real(kind=wp),             intent(out) :: yH, yHe, s_H, s_He
    real(kind=wp), optional,   intent(out) :: s_H_exc
    real(kind=wp) :: fH, fHe, s_exc
    integer :: k
    call heii_ground_shares(g, nHI, nHeI, nHeII, yH, yHe)
    s_H  = yH*g%a1
    s_He = yHe*g%a1
    s_exc = 0.0_wp
    do k = 1, NPART
       if (g%n(k) <= 0.0_wp) cycle
       call excited_shares(g, k, nHI, nHeI, fH, fHe)
       s_exc = s_exc + g%n(k)*fH
       s_He  = s_He  + g%n(k)*fHe
    end do
    if (present(s_H_exc)) then
       s_H_exc = s_exc*ne*max(nHeIII, 0.0_wp)
    else
       s_H = s_H + s_exc
    end if
    s_H  = s_H*ne*max(nHeIII, 0.0_wp)
    s_He = s_He*ne*max(nHeIII, 0.0_wp)
  end subroutine heii_on_the_spot_sources

  !=========================================================================
  ! H0's share of the excited cascade [cm^-3 s^-1] at n(H0) = nHI: the
  ! s_H_exc of heii_on_the_spot_sources alone, for the root find of
  ! hydrogen's balance (ion_balance_mod::hydrogen_fraction_with_cascade),
  ! which evaluates it at every trial x_HI.
  !=========================================================================
  real(kind=wp) function heii_cascade_to_H(g, ne, nHI, nHeI, nHeIII) result(s)
    type(heii_ots_nodes_type), intent(in) :: g
    real(kind=wp),             intent(in) :: ne, nHI, nHeI, nHeIII
    real(kind=wp) :: fH, fHe
    integer :: k
    s = 0.0_wp
    do k = 1, NPART
       if (g%n(k) <= 0.0_wp) cycle
       call excited_shares(g, k, nHI, nHeI, fH, fHe)
       s = s + g%n(k)*fH
    end do
    s = s*ne*max(nHeIII, 0.0_wp)
  end function heii_cascade_to_H

  !=========================================================================
  ! Gas heating [erg cm^-3 s^-1] of the same absorptions (header): the
  ! ground-capture photons IP(He+) - IP(X) each, the excited photons
  ! E - IP(X) at their part's mean energy.
  !=========================================================================
  !--- f_H_exc: the fraction of the photons given to H0 that it can absorb
  !--- (solve_ion_cell's f_h_only, "Absorption limit" in the header); 1 by
  !--- default.
  real(kind=wp) function heii_on_the_spot_heating(g, ne, nHI, nHeI, nHeII, &
                                                  nHeIII, f_H_exc) result(heat)
    type(heii_ots_nodes_type), intent(in) :: g
    real(kind=wp),             intent(in) :: ne, nHI, nHeI, nHeII, nHeIII
    real(kind=wp), optional,   intent(in) :: f_H_exc
    real(kind=wp) :: yH, yHe, fH, fHe, f
    integer :: k
    call heii_ground_shares(g, nHI, nHeI, nHeII, yH, yHe)
    f = 1.0_wp
    if (present(f_H_exc)) f = f_H_exc
    heat = g%a1*(f*yH*DE_HEII_H + yHe*DE_HEII_HEI)
    do k = 1, NPART
       if (g%n(k) <= 0.0_wp) cycle
       call excited_shares(g, k, nHI, nHeI, fH, fHe)
       heat = heat + g%n(k)*(f*fH*(g%e(k) - eth_HI) + fHe*(g%e(k) - eth_HeI))
    end do
    heat = heat*ne*max(nHeIII, 0.0_wp)*ev2erg
  end function heii_on_the_spot_heating

  !--- shares of excited-cascade part k taken by H0 and He0 (below the He0
  !--- edge sigma_HeI = 0, so H0 takes it all; with neither absorber the part
  !--- is H0's, and the absorption limit of solve_ion_cell decides what of
  !--- it is taken)
  subroutine excited_shares(g, k, nHI, nHeI, fH, fHe)
    type(heii_ots_nodes_type), intent(in)  :: g
    integer,                   intent(in)  :: k
    real(kind=wp),             intent(in)  :: nHI, nHeI
    real(kind=wp),             intent(out) :: fH, fHe
    real(kind=wp) :: tH, tHe
    tH  = max(nHI, 0.0_wp)*g%sH(k)
    tHe = max(nHeI, 0.0_wp)*g%sHe(k)
    fH = 1.0_wp;  fHe = 0.0_wp
    if (tH + tHe > 0.0_wp) then
       fH  = tH/(tH + tHe)
       fHe = tHe/(tH + tHe)
    end if
  end subroutine excited_shares

end module heii_recombination_photons_mod
