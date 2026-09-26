program check_hei_on_the_spot
!---------------------------------------------------------------------------
! Gate for the case-B on-the-spot He I recombination photons
! (par%hei_on_the_spot, hei_recombination_photons_mod), uncoupled 2^3S:
!  (a) the share y of the He0 ground-capture photons that H0 absorbs: its
!      limits (no He0 -> 1, no H0 -> 0) and the four-point Gauss-Laguerre
!      average against a direct quadrature of the Milne spectrum
!      E^2 sigma_HeI(E) exp[-(E - E_th)/kT] to 40 kT;
!  (b) the excited-cascade photon count per case-B recombination equals
!      f_3S + f_1P p_584 + (f_1S + f_1P (1 - p_584)) N_2q, assembled here from
!      hei_branch_fractions, hei_584_ionizes_H and the two-photon table;
!  (c) the heating equals the sum of (E - 13.598 eV) over the absorbed
!      photons plus 10.989 eV for each ground-capture photon H0 takes;
!  (d) solve_ion_cell in case B closes both balances with these terms at the
!      state it returns: n_H0 (Gamma_H + C_H n_e) + S = alpha_B n_e n_p and
!      n_He0 (Gamma_He + C_He n_e) = (alpha_B + y alpha_1) n_e n_He+.
! The .false. path is checked at run level (bit identity against the
! pre-change binary; tests/README.md).
!---------------------------------------------------------------------------
  use define
  use recomb_mod
  use photo_xsec,         only : sigma_HI, sigma_HeI
  use hei_cascade_mod,    only : hei_branch_fractions, hei_584_ionizes_H
  use hei_twophoton_mod,  only : hei_2ph_setup, hei_2ph_photons_per_decay, &
                                 hei_2ph_energy_per_decay
  use hei_recombination_photons_mod
  use ion_balance_mod,    only : solve_ion_cell
  use mpi
  implicit none
  real(kind=wp), parameter :: TS(3) = [5.0e3_wp, 1.0e4_wp, 2.0e4_wp]
  real(kind=wp), parameter :: RATIO(3) = [0.1_wp, 1.0_wp, 10.0_wp]
  real(kind=wp), parameter :: GHS(4) = [3.0e-10_wp, 3.0e-12_wp, 3.0e-13_wp, &
                                        3.0e-14_wp]
  type(hei_ground_nodes_type) :: g
  real(kind=wp) :: T, kT, y, yref, num, den, E, dE, w, tH, tHe, worst
  real(kind=wp) :: nH, ne, xHI, xHeI, xHeII, nHeII, s_gnd, s_exc, aB, a1
  real(kind=wp) :: f3s, f1s, f1p, p, cnt, heat, href, n2q, e2q, y_s
  real(kind=wp) :: gH, gHe1, gHe2, lhs, rhs, xHeIII, nHe
  integer :: it, ir, k, ierr, ig
  logical :: ok

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir   = '../../data/atomic'
  par%data_dir     = '../../data'
  par%use_metals   = .false.
  par%recombination_continuum = 'on_the_spot'
  par%abund_He     = 0.1_wp
  par%recomb_model = 'detailed_balance'
  par%ci_model     = 'voronov'
  par%hei_metastable_coupling = .false.
  par%helium_charge_exchange  = .false.
  par%hei_on_the_spot         = .true.
  !--- the He III photons are gated on their own (tests/heii_on_the_spot)
  par%heii_on_the_spot        = .false.
  call hei_2ph_setup()
  ok = .true.
  nH = 100.0_wp

  !--- (a)
  worst = 0.0_wp
  write(*,'(/,a)') '(a) share y of the He0 ground-capture photons taken by H0'
  do it = 1, 3
     T = TS(it);  kT = kboltz_cgs*T/ev2erg
     g = hei_ground_nodes(T)
     if (abs(hei_ground_share_H(g, 1.0_wp, 0.0_wp) - 1.0_wp) > 1.0e-14_wp .or. &
         hei_ground_share_H(g, 0.0_wp, 1.0_wp) /= 0.0_wp) ok = .false.
     do ir = 1, 3
        num = 0.0_wp;  den = 0.0_wp
        dE = 40.0_wp*kT/20000.0_wp
        do k = 0, 20000
           E = eth_HeI + k*dE
           w = E*E*sigma_HeI(E)*exp(-(E - eth_HeI)/kT)
           if (k == 0 .or. k == 20000) w = 0.5_wp*w
           tH = sigma_HI(E);  tHe = RATIO(ir)*sigma_HeI(E)
           num = num + w*tH/(tH + tHe);  den = den + w
        end do
        yref = num/den
        y = hei_ground_share_H(g, 1.0_wp, RATIO(ir))
        worst = max(worst, abs(y/yref - 1.0_wp))
        write(*,'(a,f8.0,a,f6.1,a,f10.6,a,f10.6)') '   T =', T, &
           '  n(He0)/n(H0) =', RATIO(ir), '  y =', y, '  direct', yref
     end do
  end do
  write(*,'(a,es10.2,a)') '   worst |y/direct - 1| = ', worst, '  (tol 1e-3)'
  if (worst > 1.0e-3_wp) ok = .false.

  !--- (b) and (c)
  worst = 0.0_wp
  n2q = hei_2ph_photons_per_decay();  e2q = hei_2ph_energy_per_decay()
  write(*,'(/,a,f8.5,a,f8.4,a)') '(b,c) N_2q = ', n2q, &
     ' H-ionizing photons and ', e2q, ' eV in band per 2^1S decay'
  do it = 1, 3
     T = TS(it);  g = hei_ground_nodes(T)
     ne = 1.1_wp*nH;  xHI = 1.0e-3_wp;  xHeI = 1.0e-2_wp;  xHeII = 0.99_wp
     nHeII = par%abund_He*nH*xHeII
     call hei_on_the_spot_sources(1, T, nH, ne, xHI, xHeI, xHeII, g, &
                                  y_s, s_gnd, s_exc)
     aB = alphaB_HeII(T);  a1 = alphaA_HeII(T) - aB
     call hei_branch_fractions(T, f3s, f1s, f1p)
     p = hei_584_ionizes_H(T, xHI, xHeI)
     cnt = f3s + f1p*p + (f1s + f1p*(1.0_wp - p))*n2q
     worst = max(worst, abs(s_exc/(ne*nHeII*aB)/cnt - 1.0_wp))
     y = hei_ground_share_H(g, nH*xHI, par%abund_He*nH*xHeI)
     worst = max(worst, abs(s_gnd/(y*a1*ne*nHeII) - 1.0_wp))
     heat = hei_on_the_spot_heating(1, T, nH, ne, xHI, xHeI, xHeII, g)
     href = ne*nHeII*( y*a1*(24.587387_wp - 13.598434_wp) &
          + aB*( f3s*(19.82_wp - eth_HI) + f1p*p*(21.22_wp - eth_HI) &
               + (f1s + f1p*(1.0_wp - p))*(e2q - n2q*eth_HI) ) )*ev2erg
     worst = max(worst, abs(heat/href - 1.0_wp))
     write(*,'(a,f8.0,a,f7.4,a,f7.4,a,f7.4,a,f8.3,a)') '   T =', T, &
        '  H-ionizing photons per recombination: excited', cnt, &
        '  ground y a1/aB', y*a1/aB, '  p_584', p, &
        '  mean heat per excited photon', &
        (href/(ne*nHeII*ev2erg) - y*a1*10.988953_wp)/(aB*cnt), ' eV'
  end do
  write(*,'(a,es10.2,a)') '   worst relative = ', worst, '  (tol 1e-12)'
  if (worst > 1.0e-12_wp) ok = .false.

  !--- (d)
  worst = 0.0_wp
  write(*,'(/,a)') '(d) balances closed by solve_ion_cell, case B'
  do it = 1, 3
     T = TS(it);  g = hei_ground_nodes(T)
     aB = alphaB_HeII(T);  a1 = alphaA_HeII(T) - aB
     do ig = 1, 4
        gH = GHS(ig);  gHe1 = 0.3_wp*GHS(ig);  gHe2 = 1.0e-3_wp*GHS(ig)
        xHI = 0.5_wp;  xHeI = 0.5_wp;  xHeII = 0.5_wp
        call solve_ion_cell(gH, gHe1, gHe2, nH, T, .false., xHI, xHeI, xHeII, &
                            ne, 1)
        call hei_on_the_spot_sources(1, T, nH, ne, xHI, xHeI, xHeII, g, &
                                     y_s, s_gnd, s_exc)
        lhs = nH*xHI*(gH + ci_HI(T)*ne) + s_gnd + s_exc
        rhs = alphaB_HII(T)*ne*nH*(1.0_wp - xHI)
        worst = max(worst, abs(lhs/rhs - 1.0_wp))
        nHe = par%abund_He*nH
        xHeIII = max(0.0_wp, 1.0_wp - xHeI - xHeII)
        lhs = nHe*xHeI*(gHe1 + ci_HeI(T)*ne)
        rhs = (aB + y_s*a1)*ne*nHe*xHeII
        worst = max(worst, abs(lhs/rhs - 1.0_wp))
        write(*,'(a,f8.0,a,es9.2,a,f9.6,a,f8.5,a,es9.2)') '   T =', T, &
           '  Gamma_H =', gH, '  x_HII =', 1.0_wp - xHI, '  y =', y_s, &
           '  S/(alpha_B n_e n_p) =', &
           (s_gnd + s_exc)/(alphaB_HII(T)*ne*nH*(1.0_wp - xHI))
     end do
  end do
  write(*,'(a,es10.2,a)') '   worst balance residual = ', worst, '  (tol 1e-8)'
  if (worst > 1.0e-8_wp) ok = .false.

  if (ok) then
     write(*,'(/,a)') ' GATE (case-B on-the-spot He I photons): PASS'
  else
     write(*,'(/,a)') ' GATE (case-B on-the-spot He I photons): FAIL'
     call MPI_FINALIZE(ierr);  stop 1
  end if
  call MPI_FINALIZE(ierr)
end program check_hei_on_the_spot
