program check_heii_on_the_spot
!---------------------------------------------------------------------------
! Gate for the case-B on-the-spot He III recombination photons
! (par%heii_on_the_spot, heii_recombination_photons_mod):
!  (a) the shares y_H, y_He of the He II ground-capture photons: their
!      limits (only He+ -> 0 and 0; only H0 -> 1; only He0 -> y_He = 1) and
!      the four-point Gauss-Laguerre average against a direct quadrature of
!      the Milne spectrum E^2 sigma_HeII(E) exp[-(E - E_th)/kT] to 40 kT;
!  (b) the excited cascade counted by heii_cascade_photons closes on
!      heii_branch_luminosity path by path to rounding (energy), and its
!      photon count equals the one the transported channel's samplers imply,
!      L <1/E> over heii_fb_sample_energy / heii_2q_sample_energy;
!  (c) the H0/He0 shares of the excited photons (parts at their mean
!      energy) against a direct quadrature of the free-bound and two-photon
!      spectra with the share evaluated photon by photon;
!  (d) the source and heating limits: no He0 -> all to H; no H0 -> the part
!      below the He0 edge is still H0's (its share as n(H0) -> 0+, which the
!      absorption limit of solve_ion_cell then bounds), the rest to He0;
!  (e) the heating assembled here independently, from the parts and the
!      cross sections, against heii_on_the_spot_heating;
!  (f) solve_ion_cell in case B closes H, He0/He+ and He+/He2+ with these
!      terms (and the He I ones) at the state it returns;
!  (g) the absorption limit (below);
!  (h) par%heii_diffuse in the on-the-spot mode: heii_cascade_transported
!      is true only for 'on_the_spot' and 'transported' with the switch on;
!      with it on the nodes carry no cascade photons, so the sources and the
!      heat are the ground captures' alone, equal to those computed here
!      from the ground shares, and the ground-capture part is unchanged.
! The .false. path is checked at run level (bit identity against the
! pre-change binary; tests/README.md).
!---------------------------------------------------------------------------
  use define
  use recomb_mod
  use photo_xsec,           only : sigma_HI, sigma_HeI, sigma_HeII
  use recomb_continuum_mod, only : gamma_tab, gamma_ff
  use hydrogenic_twophoton_mod, only : h2q_photons_per_y, &
                                   h2q_band_photons_per_decay, &
                                   h2q_band_energy_per_decay
  use heii_cascade_mod
  use hei_twophoton_mod,    only : hei_2ph_setup
  use hei_cascade_mod,      only : hei_cascade_setup
  use hei_recombination_photons_mod, only : hei_ground_nodes, &
                                   hei_on_the_spot_sources
  use heii_recombination_photons_mod
  use ion_balance_mod,      only : solve_ion_cell
  use mpi
  implicit none
  integer, parameter :: NG1 = HEII_NGROUP + 1
  real(kind=wp), parameter :: TS(3) = [1.0e4_wp, 2.0e4_wp, 4.0e4_wp]
  real(kind=wp), parameter :: RYD_EV = 13.605693_wp
  real(kind=wp), parameter :: RATIO(3) = [0.01_wp, 0.1_wp, 1.0_wp]
  real(kind=wp), parameter :: GHS(4) = [3.0e-10_wp, 3.0e-12_wp, 3.0e-13_wp, &
                                        3.0e-14_wp]
  real(kind=wp), parameter :: HE2_TO_H(4) = [1.0e-2_wp, 3.0e-2_wp, 0.1_wp, 0.3_wp]
  type(heii_ots_nodes_type)       :: g
  type(heii_cascade_photons_type) :: c
  real(kind=wp) :: T, kT, E, dE, w, tH, tHe, tHe2, num1, num2, den, yH, yHe, f_abs
  real(kind=wp) :: worst, l_fb, l_lyman, l_2q, en, ph, u, x, fH, fHe, a, b, h
  real(kind=wp) :: ne, nHI, nHeI, nHeII, nHeIII, sH, sHe, heat, href, rec
  real(kind=wp) :: yref(2), tHref, nfine, sfine, spart, npart, lt
  real(kind=wp) :: nH, xHI, xHeI, xHeII, xHeIII, gH, gHe1, gHe2, nHe, lhs, rhs
  real(kind=wp) :: ys1, sg1, se1, aB, aB3, a13
  real(kind=wp) :: dens(3,4)
  integer :: it, k, i, j, ir, ierr, ig, NS, ia
  logical :: ok

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir   = '../../data/atomic'
  par%data_dir     = '../../data'
  par%use_metals   = .false.
  par%recombination_continuum = 'on_the_spot'
  par%abund_He     = 0.1_wp
  par%recomb_model = 'detailed_balance'
  par%ci_model     = 'voronov'
  par%eion_max     = 30.0_wp      ! below Ly-alpha: the on-the-spot tables must not care
  par%hei_metastable_coupling = .false.
  par%helium_charge_exchange  = .false.
  par%hei_on_the_spot         = .true.
  par%heii_on_the_spot        = .true.
  call hei_2ph_setup();  call hei_cascade_setup()
  call heii_cascade_setup()
  ok = .true.

  !--- (a) ground-capture shares -------------------------------------------
  !--- densities (n_H0, n_He0, n_He+) of the share tests
  dens(:,1) = [1.0_wp, 0.1_wp, 1.0_wp]
  dens(:,2) = [1.0e-3_wp, 1.0e-4_wp, 1.0_wp]
  dens(:,3) = [1.0_wp, 1.0_wp, 1.0e-2_wp]
  dens(:,4) = [1.0e-2_wp, 1.0_wp, 0.1_wp]
  worst = 0.0_wp
  write(*,'(/,a)') '(a) shares of the He II ground-capture photons (H0, He0)'
  do it = 1, 3
     T = TS(it);  kT = kboltz_cgs*T/ev2erg
     g = heii_ots_nodes(T)
     call heii_ground_shares(g, 0.0_wp, 0.0_wp, 1.0_wp, yH, yHe)
     if (yH /= 0.0_wp .or. yHe /= 0.0_wp) ok = .false.
     call heii_ground_shares(g, 1.0_wp, 0.0_wp, 0.0_wp, yH, yHe)
     if (abs(yH - 1.0_wp) > 1.0e-14_wp .or. yHe /= 0.0_wp) ok = .false.
     call heii_ground_shares(g, 0.0_wp, 1.0_wp, 0.0_wp, yH, yHe)
     if (abs(yHe - 1.0_wp) > 1.0e-14_wp .or. yH /= 0.0_wp) ok = .false.
     do ir = 1, 4
        num1 = 0.0_wp;  num2 = 0.0_wp;  den = 0.0_wp
        dE = 40.0_wp*kT/20000.0_wp
        do k = 0, 20000
           E = eth_HeII + k*dE
           w = E*E*sigma_HeII(E)*exp(-(E - eth_HeII)/kT)
           if (k == 0 .or. k == 20000) w = 0.5_wp*w
           tH = dens(1,ir)*sigma_HI(E);  tHe = dens(2,ir)*sigma_HeI(E)
           tHe2 = dens(3,ir)*sigma_HeII(E)
           num1 = num1 + w*tH/(tH + tHe + tHe2)
           num2 = num2 + w*tHe/(tH + tHe + tHe2)
           den  = den + w
        end do
        yref = [num1/den, num2/den]
        call heii_ground_shares(g, dens(1,ir), dens(2,ir), dens(3,ir), yH, yHe)
        worst = max(worst, abs(yH/yref(1) - 1.0_wp), abs(yHe/yref(2) - 1.0_wp))
        write(*,'(a,f8.0,a,3es9.1,a,2f10.6,a,2f10.6)') '   T =', T, &
           '  n(H0,He0,He+) =', dens(:,ir), '  y =', yH, yHe, '  direct', yref
     end do
  end do
  write(*,'(a,es10.2,a)') '   worst |y/direct - 1| = ', worst, '  (tol 1e-3)'
  if (worst > 1.0e-3_wp) ok = .false.

  !--- (b) closure on the branch luminosities and the samplers ---------------
  NS = 200000
  worst = 0.0_wp
  write(*,'(/,a)') '(b) excited cascade: energy closure on heii_branch_luminosity,'
  write(*,'(a)')   '    photons against the transported samplers L <1/E>'
  do it = 1, 3
     !--- a tabulated free-bound temperature (the sampler takes the nearest)
     lt = 3.0_wp + 2.0_wp*real(nint(1.0_wp + 23.5_wp*(log10(TS(it)) - 3.0_wp)) - 1, wp)/47.0_wp
     T = 10.0_wp**lt
     c = heii_cascade_photons(T)
     call heii_branch_luminosity(T, l_fb, l_lyman, l_2q)
     en = sum(c%n(:,1)*c%e(:,1))*ev2erg
     worst = max(worst, abs(en/l_fb - 1.0_wp))
     en = sum(c%n(:,2)*c%e(:,2))*ev2erg
     worst = max(worst, abs(en/l_lyman - 1.0_wp))
     en = sum(c%n(:,3)*c%e(:,3))*ev2erg
     worst = max(worst, abs(en/l_2q - 1.0_wp))
     !--- photons implied by the samplers: L/ev2erg x mean of 1/E
     num1 = 0.0_wp;  num2 = 0.0_wp
     do i = 1, NS
        u = (real(i,wp) - 0.5_wp)/real(NS,wp)
        num1 = num1 + 1.0_wp/heii_fb_sample_energy(T, u)
        num2 = num2 + 1.0_wp/heii_2q_sample_energy(u)
     end do
     num1 = l_fb/ev2erg*num1/real(NS,wp)
     num2 = l_2q/ev2erg*num2/real(NS,wp)
     write(*,'(a,f8.0,a,2f10.6,a,f10.6,a,f10.6)') '   T =', T, &
        '  photons fb, 2q / sampler implied:', sum(c%n(:,1))/num1, &
        sum(c%n(:,3))/num2, '   2q photons / (aB f2s N_2q):', &
        sum(c%n(:,3))/(alphaB_HeIII(T)*heii_f_2s(T)*h2q_band_photons_per_decay()), &
        '  Lya / (aB (1-f2s)):', c%n(1,2)/(alphaB_HeIII(T)*(1.0_wp - heii_f_2s(T)))
     if (abs(sum(c%n(:,1))/num1 - 1.0_wp) > 1.0e-3_wp) ok = .false.
     if (abs(sum(c%n(:,3))/num2 - 1.0_wp) > 1.0e-3_wp) ok = .false.
     if (abs(sum(c%n(:,3))/(alphaB_HeIII(T)*heii_f_2s(T) &
             *h2q_band_photons_per_decay()) - 1.0_wp) > 1.0e-6_wp) ok = .false.
     if (abs(c%n(1,2)/(alphaB_HeIII(T)*(1.0_wp - heii_f_2s(T))) - 1.0_wp) &
         > 1.0e-13_wp) ok = .false.
     !--- H-ionizing photons per H recombination at n(He2+)/n(H+) = 0.15,
     !--- the form of Osterbrock & Ferland (2006) Table 2.6
     write(*,'(a,3f7.3,a)') '      per H recombination at n(He2+)/n(H+) = 0.15: '// &
        'Ly-alpha, two-photon, free-bound =', &
        0.15_wp*[sum(c%n(:,2)), sum(c%n(:,3)), sum(c%n(:,1))]/alphaB_HII(T), &
        '  (O&F Table 2.6, 1e4 K: 0.64 0.36 0.20)'
  end do
  write(*,'(a,es10.2,a)') '   worst energy closure = ', worst, '  (tol 1e-13)'
  if (worst > 1.0e-13_wp) ok = .false.
  write(*,'(a,f9.5,a,f9.4,a)') '   two-photon window: N = ', &
     h2q_band_photons_per_decay(), ' photons,', h2q_band_energy_per_decay(), &
     ' eV per decay (full 13.598-40.813 eV although eion_max = 30 eV)'
  if (abs(h2q_band_energy_per_decay() - 36.09_wp) > 0.01_wp) ok = .false.

  !--- (c) excited-photon shares against photon-by-photon quadrature --------
  worst = 0.0_wp
  write(*,'(/,a)') '(c) H0 share of the excited photons above the He0 edge'
  do it = 1, 3
     lt = 3.0_wp + 2.0_wp*real(nint(1.0_wp + 23.5_wp*(log10(TS(it)) - 3.0_wp)) - 1, wp)/47.0_wp
     T = 10.0_wp**lt
     g = heii_ots_nodes(T);  c = heii_cascade_photons(T)
     do ir = 1, 3
        nHI = 1.0_wp;  nHeI = RATIO(ir)
        do j = 1, 3, 2        ! free-bound, two-photon
           !--- parts: sum n_k f_H(E_k) / sum n_k over the parts above the edge
           npart = 0.0_wp;  spart = 0.0_wp
           do k = 1, HEII_NGROUP
              i = (j - 1)*NG1 + k + 1
              tH = nHI*g%sH(i);  tHe = nHeI*g%sHe(i)
              npart = npart + g%n(i);  spart = spart + g%n(i)*tH/(tH + tHe)
           end do
           !--- direct: photon spectrum n(E) from the edge to the top
           if (j == 1) then
              a = eth_HeI;  b = eth_HeII
           else
              a = eth_HeI;  b = E_HEII_LYA
           end if
           nfine = 0.0_wp;  sfine = 0.0_wp
           h = (b - a)/40000.0_wp
           do i = 0, 40000
              E = a + h*real(i,wp)
              if (j == 1) then
                 x = max(gamma_tab(3, E/RYD_EV, T) - gamma_ff(E/RYD_EV, T, 2), 0.0_wp)/E
              else
                 x = h2q_photons_per_y(E/E_HEII_LYA)
              end if
              if (i == 0 .or. i == 40000) x = 0.5_wp*x
              tH = nHI*sigma_HI(E);  tHe = nHeI*sigma_HeI(E)
              nfine = nfine + x;  sfine = sfine + x*tH/(tH + tHe)
           end do
           worst = max(worst, abs((spart/npart)/(sfine/nfine) - 1.0_wp))
           write(*,'(a,f8.0,a,f5.2,a,a,a,f10.6,a,f10.6)') '   T =', T, &
              '  n(He0)/n(H0) =', nHeI, '  ', merge('free-bound', 'two-photon', &
              j == 1), '  parts', spart/npart, '  direct', sfine/nfine
        end do
     end do
  end do
  write(*,'(a,es10.2,a)') '   worst |parts/direct - 1| = ', worst, '  (tol 1e-3)'
  if (worst > 1.0e-3_wp) ok = .false.

  !--- (d) limits ---------------------------------------------------------
  write(*,'(/,a)') '(d) limits of the sources'
  T = 2.0e4_wp;  g = heii_ots_nodes(T)
  ne = 110.0_wp;  nHeIII = 9.0_wp
  call heii_on_the_spot_sources(g, ne, 1.0e-2_wp, 0.0_wp, 1.0_wp, nHeIII, &
                                yH, yHe, sH, sHe)
  rec = ne*nHeIII
  num1 = (sum(g%n) + yH*g%a1)*rec
  write(*,'(a,es12.5,a,es12.5,a,es9.2)') '   no He0: S_H =', sH, &
     '  all photons + y_H alpha_1 n_e n(He2+) =', num1, '  S_He =', sHe
  if (abs(sH/num1 - 1.0_wp) > 1.0e-13_wp .or. sHe /= 0.0_wp) ok = .false.
  call heii_on_the_spot_sources(g, ne, 0.0_wp, 1.0e-2_wp, 1.0_wp, nHeIII, &
                                yH, yHe, sH, sHe)
  num1 = 0.0_wp;  num2 = 0.0_wp
  do i = 1, size(g%n)
     if (g%e(i) >= eth_HeI) then
        num1 = num1 + g%n(i)
     else
        num2 = num2 + g%n(i)
     end if
  end do
  num1 = (num1 + yHe*g%a1)*rec
  num2 = num2*rec
  write(*,'(a,es12.5,a,es12.5,a,es12.5,a,es12.5)') '   no H0:  S_He =', sHe, &
     '  photons above 24.587 eV + y_He alpha_1 n_e n(He2+) =', num1, &
     '  S_H =', sH, '  photons below it =', num2
  if (abs(sHe/num1 - 1.0_wp) > 1.0e-13_wp .or. &
      abs(sH/num2 - 1.0_wp) > 1.0e-13_wp) ok = .false.
  call heii_on_the_spot_sources(g, ne, 1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, &
                                yH, yHe, sH, sHe)
  if (sH /= 0.0_wp .or. sHe /= 0.0_wp) ok = .false.

  !--- (e) heating assembled independently ---------------------------------
  worst = 0.0_wp
  write(*,'(/,a)') '(e) heating'
  do it = 1, 3
     T = TS(it);  g = heii_ots_nodes(T);  c = heii_cascade_photons(T)
     ne = 110.0_wp;  nHI = 1.0e-2_wp;  nHeI = 3.0e-3_wp;  nHeII = 1.0_wp
     nHeIII = 9.0_wp
     heat = heii_on_the_spot_heating(g, ne, nHI, nHeI, nHeII, nHeIII)
     call heii_ground_shares(g, nHI, nHeI, nHeII, yH, yHe)
     href = (alphaA_HeIII(T) - alphaB_HeIII(T)) &
            *(yH*(54.417765_wp - 13.598434_wp) + yHe*(54.417765_wp - 24.587387_wp))
     do j = 1, 3
        do k = 0, HEII_NGROUP
           if (c%n(k,j) <= 0.0_wp) cycle
           tH  = nHI*sigma_HI(c%e(k,j));  tHe = nHeI*sigma_HeI(c%e(k,j))
           href = href + c%n(k,j)*(tH*(c%e(k,j) - eth_HI) &
                                   + tHe*(c%e(k,j) - eth_HeI))/(tH + tHe)
        end do
     end do
     href = href*ne*nHeIII*ev2erg
     call heii_on_the_spot_sources(g, ne, nHI, nHeI, nHeII, nHeIII, yH, yHe, sH, sHe)
     worst = max(worst, abs(heat/href - 1.0_wp))
     write(*,'(a,f8.0,a,es12.5,a,f8.3,a,f8.4,a,f8.4)') '   T =', T, &
        '  heat =', heat, '  eV per ionization', heat/ev2erg/(sH + sHe), &
        '  y_H =', yH, '  y_He =', yHe
  end do
  write(*,'(a,es10.2,a)') '   worst relative = ', worst, '  (tol 1e-13)'
  if (worst > 1.0e-13_wp) ok = .false.

  !--- (f) balances closed by solve_ion_cell ------------------------------
  worst = 0.0_wp
  nH = 100.0_wp;  nHe = par%abund_He*nH
  write(*,'(/,a)') '(f) balances closed by solve_ion_cell, case B, He I and He II on the spot'
  do it = 1, 3
     T = TS(it)
     aB = alphaB_HeII(T);  aB3 = alphaB_HeIII(T);  a13 = alphaA_HeIII(T) - aB3
     g = heii_ots_nodes(T)
     do ig = 1, 4
        gH = GHS(ig)
        gHe1 = 0.3_wp*gH;  gHe2 = HE2_TO_H(ig)*gH
        xHI = 0.5_wp;  xHeI = 0.3_wp;  xHeII = 0.3_wp
        call solve_ion_cell(gH, gHe1, gHe2, nH, T, .false., xHI, xHeI, xHeII, &
                            ne, 1)
        xHeIII = max(0.0_wp, 1.0_wp - xHeI - xHeII)
        call hei_on_the_spot_sources(1, T, nH, ne, xHI, xHeI, xHeII, &
                                     hei_ground_nodes(T), ys1, sg1, se1)
        call heii_on_the_spot_sources(g, ne, nH*xHI, nHe*xHeI, nHe*xHeII, &
                                      nHe*xHeIII, yH, yHe, sH, sHe)
        lhs = nH*xHI*(gH + ci_HI(T)*ne) + sg1 + se1 + sH
        rhs = alphaB_HII(T)*ne*nH*(1.0_wp - xHI)
        worst = max(worst, abs(lhs/rhs - 1.0_wp))
        lhs = nHe*xHeI*(gHe1 + ci_HeI(T)*ne) + sHe
        rhs = (aB + ys1*(alphaA_HeII(T) - aB))*ne*nHe*xHeII
        worst = max(worst, abs(lhs/rhs - 1.0_wp))
        lhs = nHe*xHeII*(gHe2 + ci_HeII(T)*ne)
        rhs = (aB3 + (yH + yHe)*a13)*ne*nHe*xHeIII
        worst = max(worst, abs(lhs/rhs - 1.0_wp))
        write(*,'(a,f8.0,a,es9.2,a,3es10.3,a,f8.5,a,f8.5)') '   T =', T, &
           '  Gamma_H =', gH, '  x(H0,He0,He2+) =', xHI, xHeI, xHeIII, &
           '  S_H(He II)/(alpha_B n_e n_p) =', sH/(alphaB_HII(T)*ne*nH*(1.0_wp - xHI)), &
           '  y_H+y_He =', yH + yHe
     end do
  end do
  write(*,'(a,es10.2,a)') '   worst balance residual = ', worst, '  (tol 1e-8)'
  if (worst > 1.0e-8_wp) ok = .false.

  !--- (g) the absorption limit: He/H = 0.2 fully doubly ionized at 1.8e4 K,
  !--- where the cascades given to H0 outnumber hydrogen's recombinations.
  !--- solve_ion_cell must leave x_HI ~ 0, return the absorbed fraction f < 1,
  !--- and close hydrogen's balance with f times the source.  (A fully He2+
  !--- cell at He/H = 0.1 must keep f = 1.)
  write(*,'(/,a)') '(g) absorption limit (He/H = 0.2 and 0.1, 1.8e4 K, He2+ zone)'
  worst = 0.0_wp
  do ia = 1, 2
     par%abund_He = merge(0.2_wp, 0.1_wp, ia == 1)
     nHe = par%abund_He*nH
     T = 1.8e4_wp;  g = heii_ots_nodes(T)
     gH = 1.0e-11_wp;  gHe1 = gH;  gHe2 = 1.0e-9_wp
     xHI = 1.0e-3_wp;  xHeI = 1.0e-3_wp;  xHeII = 1.0e-2_wp
     call solve_ion_cell(gH, gHe1, gHe2, nH, T, .false., xHI, xHeI, xHeII, &
                         ne, 1, f_abs)
     xHeIII = max(0.0_wp, 1.0_wp - xHeI - xHeII)
     call hei_on_the_spot_sources(1, T, nH, ne, xHI, xHeI, xHeII, &
                                  hei_ground_nodes(T), ys1, sg1, se1)
     call heii_on_the_spot_sources(g, ne, nH*xHI, nHe*xHeI, nHe*xHeII, &
                                   nHe*xHeIII, yH, yHe, sH, sHe)
     lhs = nH*xHI*(gH + ci_HI(T)*ne) + f_abs*(sg1 + se1 + sH)
     rhs = alphaB_HII(T)*ne*nH*(1.0_wp - xHI)
     worst = max(worst, abs(lhs/rhs - 1.0_wp))
     write(*,'(a,f5.2,a,es10.3,a,f8.5,a,f8.5,a,es9.2)') '   He/H =', par%abund_He, &
        '  x_HI =', xHI, '  source/supply =', (sg1 + se1 + sH)/rhs, &
        '  f =', f_abs, '  |balance - 1| =', abs(lhs/rhs - 1.0_wp)
     write(*,'(a,3es11.3,a,es11.3)') '      x(He0, He+, He2+) =', xHeI, xHeII, &
        xHeIII, '  n_e/n_H =', ne/nH
     if (ia == 1 .and. .not. (xHI < 1.0e-9_wp .and. f_abs < 1.0_wp)) ok = .false.
     if (ia == 2 .and. f_abs /= 1.0_wp) ok = .false.
  end do
  par%abund_He = 0.1_wp
  write(*,'(a,es10.2,a)') '   worst balance residual = ', worst, '  (tol 1e-8)'
  if (worst > 1.0e-8_wp) ok = .false.

  !--- (h) the He III cascade transported alone from the on-the-spot mode
  write(*,'(/,a)') '(h) heii_diffuse with the on-the-spot continuum'
  par%recombination_continuum = 'on_the_spot';  par%heii_diffuse = .false.
  a = merge(1.0_wp, 0.0_wp, heii_cascade_transported())
  par%heii_diffuse = .true.
  b = merge(1.0_wp, 0.0_wp, heii_cascade_transported())
  par%recombination_continuum = 'escaping'
  h = merge(1.0_wp, 0.0_wp, heii_cascade_transported())
  par%recombination_continuum = 'transported'
  x = merge(1.0_wp, 0.0_wp, heii_cascade_transported())
  write(*,'(a,4f4.0)') '   transported? ots off, ots on, escaping on, transported on:', a, b, h, x
  if (.not. (a == 0.0_wp .and. b == 1.0_wp .and. h == 0.0_wp .and. x == 1.0_wp)) ok = .false.
  par%recombination_continuum = 'on_the_spot'
  T = 1.8e4_wp;  ne = 120.0_wp
  nHI = 1.0e-3_wp;  nHeI = 1.0e-4_wp;  nHeII = 2.0_wp;  nHeIII = 8.0_wp
  par%heii_diffuse = .false.
  g = heii_ots_nodes(T)
  call heii_on_the_spot_sources(g, ne, nHI, nHeI, nHeII, nHeIII, yH, yHe, a, b)
  href = heii_on_the_spot_heating(g, ne, nHI, nHeI, nHeII, nHeIII)
  par%heii_diffuse = .true.
  g = heii_ots_nodes(T)
  call heii_on_the_spot_sources(g, ne, nHI, nHeI, nHeII, nHeIII, yH, yHe, sH, sHe)
  heat = heii_on_the_spot_heating(g, ne, nHI, nHeI, nHeII, nHeIII)
  worst = max(abs(sH/(yH*g%a1*ne*nHeIII) - 1.0_wp), &
              abs(sHe/(yHe*g%a1*ne*nHeIII) - 1.0_wp), &
              abs(heat/(g%a1*ne*nHeIII*ev2erg*(yH*DE_HEII_H &
                  + yHe*DE_HEII_HEI)) - 1.0_wp))
  write(*,'(a,es10.2,a,i0)') '   ground-only sources and heat: worst ', worst, &
     ';  cascade photons left on the spot: ', count(g%n /= 0.0_wp)
  write(*,'(a,2f9.5)') '   H0 source, heat with / without the cascade on the spot:', &
     sH/a, heat/href
  if (worst > 1.0e-12_wp .or. any(g%n /= 0.0_wp)) ok = .false.
  if (.not. (sH < a .and. heat < href)) ok = .false.
  par%heii_diffuse = .false.

  if (ok) then
     write(*,'(/,a)') ' GATE (case-B on-the-spot He III photons): PASS'
  else
     write(*,'(/,a)') ' GATE (case-B on-the-spot He III photons): FAIL'
     call MPI_FINALIZE(ierr);  stop 1
  end if
  call MPI_FINALIZE(ierr)
end program check_heii_on_the_spot
