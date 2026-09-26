program check_hei_metastable_coupling
!---------------------------------------------------------------------------
! Gate: the He I 2^3S population coupled into the ionization balance, the
! thermal balance and the He I diffuse channel (par%hei_metastable_coupling,
! hei_metastable_mod).
!  (a) the electron-impact rows reproduce Bray, Burgess, Fursa & Tully (2000,
!      A&AS 146, 481) at their own nodes: for the four rows that carry the
!      network (1^1S-2^3S, 2^3S-2^1S, 2^3S-2^1P, 2^3S-2^3P) the values below
!      are typed from the printed Tables 2 and 3, independently of the data
!      file; and the rate the module forms at a node is
!      beta Upsilon exp(-dE/kT)/(g_l sqrt T) of that value;
!  (b) alpha_3 = f_3S alpha_B(He II) (hei_cascade_mod x recomb_mod);
!  (c) the closed-form n_3 against a direct 2x2 linear solve of the
!      1^1S / 2^3S rate equations with n_1 + n_3 = n_HeI;
!  (d) the energy ledger closes: with n_3 steady, recombination E_3 plus the
!      photon energy Phi_3 absorbs equals the module's net gas heating plus
!      the photon energy the levels emit plus the ionization energy stored,
!      the photon side assembled here from the rows' own energies;
!  (e) the diffuse channel-4 weights reduce to the uncoupled ones in the
!      low-density limit (n_e -> 0, Phi_3 = 0, n_HI = 0, no collisional
!      excitation from 1^1S);
!  (f) solve_ion_cell with the coupling reaches a state in which the H and
!      He0/He+ balances, WITH the 2^3S terms evaluated at that state, close,
!      in partially neutral gas (x_HII 0.02-0.9) where the lagged Penning
!      term matters most, with and without helium charge exchange.
! The .false. path is checked at run level instead (bit identity against
! the pre-change binary; tests/README.md).
!---------------------------------------------------------------------------
  use define
  use recomb_mod
  use hei_cascade_mod,    only : hei_branch_fractions
  use hei_metastable_mod
  use ion_balance_mod,    only : solve_ion_cell
  use diffuse_mod,        only : hei_branch_energies, hei_branch_energies_coupled
  use hei_twophoton_mod,  only : hei_2ph_setup
  use mpi
  implicit none
  integer, parameter :: NN = 9
  real(kind=wp), parameter :: BRAY_LOGT(NN) = [3.75_wp, 4.00_wp, 4.25_wp, &
       4.50_wp, 4.75_wp, 5.00_wp, 5.25_wp, 5.50_wp, 5.75_wp]
  !--- Table 2, 1^1S - n^3S, n = 2
  real(kind=wp), parameter :: U13(NN) = [6.198e-2_wp, 6.458e-2_wp, 6.387e-2_wp, &
       6.157e-2_wp, 5.832e-2_wp, 5.320e-2_wp, 4.787e-2_wp, 4.018e-2_wp, 3.167e-2_wp]
  !--- Table 3, 2^3S - n^1S, n = 2
  real(kind=wp), parameter :: U31A(NN) = [2.389_wp, 2.456_wp, 2.275_wp, 1.916_wp, &
       1.496_wp, 1.111_wp, 8.003e-1_wp, 5.660e-1_wp, 3.944e-1_wp]
  !--- Table 3, 2^3S - n^1P, n = 2
  real(kind=wp), parameter :: U31B(NN) = [7.965e-1_wp, 9.579e-1_wp, 1.042_wp, &
       1.015_wp, 8.950e-1_wp, 7.265e-1_wp, 5.516e-1_wp, 3.948e-1_wp, 2.677e-1_wp]
  !--- Table 3, 2^3S - n^3P, n = 2
  real(kind=wp), parameter :: U33P(NN) = [1.508e1_wp, 2.580e1_wp, 4.185e1_wp, &
       6.455e1_wp, 9.526e1_wp, 1.351e2_wp, 1.842e2_wp, 2.401e2_wp, 2.976e2_wp]
  real(kind=wp), parameter :: kb_eV = kboltz_cgs/ev2erg
  type(hei_triplet_rates_type) :: r
  real(kind=wp) :: T, worst, ups, ref, q, qsum(6), ne, nHI, nHeI, nHeII, nH
  real(kind=wp) :: phi3, heat3, n3, n3d, a11, a12, a21, a22, b1, b2, det
  real(kind=wp) :: f3s, f1s, f1p, lhs, rhs, h, pem, stored, e3, D, x, h1, h2
  real(kind=wp) :: e3s, e1p, e1s, c3s, c1p, c1s, rec, gH, gHe1, gHe2, u3
  real(kind=wp) :: aHe2, eh, e1, xHeIII, cH, aH, cHe1, k21, k12, nHII
  integer :: i, j, ierr, it, ih, icx
  logical :: ok, found(4)
  real(kind=wp) :: Ts(5) = [3.0e3_wp, 8.0e3_wp, 1.0e4_wp, 1.5e4_wp, 2.5e4_wp]
  real(kind=wp) :: GHs(4) = [3.0e-10_wp, 3.0e-12_wp, 3.0e-13_wp, 3.0e-14_wp]

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir   = '../../data/atomic'
  par%data_dir     = '../../data'
  par%use_metals   = .false.
  par%recombination_continuum = 'transported'
  par%abund_He     = 0.1_wp
  par%ci_model     = 'voronov'
  par%recomb_model = 'detailed_balance'
  par%hei_metastable_coupling = .true.
  call hei_metastable_load()
  call hei_2ph_setup()                  ! the 2^1S two-photon shape, for (e)
  ok = .true.

  !--- (a) Bray nodes
  worst = 0.0_wp
  found = .false.
  do i = 1, hm_nrow
     do j = 1, NN
        T = 10.0_wp**BRAY_LOGT(j)
        ref = -1.0_wp
        if (hm_lower(i) == 1 .and. all(hm_upper(:,i) == [2, 3, 0])) then
           ref = U13(j);   found(1) = .true.
        else if (hm_lower(i) == 2 .and. all(hm_upper(:,i) == [2, 1, 0])) then
           ref = U31A(j);  found(2) = .true.
        else if (hm_lower(i) == 2 .and. all(hm_upper(:,i) == [2, 1, 1])) then
           ref = U31B(j);  found(3) = .true.
        else if (hm_lower(i) == 2 .and. all(hm_upper(:,i) == [2, 3, 1])) then
           ref = U33P(j);  found(4) = .true.
        end if
        if (ref > 0.0_wp) worst = max(worst, abs(hei_upsilon(i, T)/ref - 1.0_wp))
     end do
  end do
  write(*,'(a,es9.2,a,4l2)') ' (a) Bray Tables 2/3 at their nodes (4 network rows): worst ', &
     worst, '   rows found', found
  if (worst > 1.0e-12_wp .or. .not. all(found)) ok = .false.
  !--- the file's Maxwellian constant is h^2/[(2 pi m_e)^1.5 k^0.5]
  write(*,'(a,es9.2)') '     BETA against CODATA 2018: |ratio-1| = ', &
     abs(hm_beta/hm_beta_value() - 1.0_wp)
  if (abs(hm_beta/hm_beta_value() - 1.0_wp) > 1.0e-9_wp) ok = .false.
  !--- the rate the module sums at a node is the formula of that Upsilon
  worst = 0.0_wp
  do j = 1, NN
     T = 10.0_wp**BRAY_LOGT(j)
     r = hei_triplet_rates(T)
     qsum = 0.0_wp
     do i = 1, hm_nrow
        ups = hei_upsilon(i, T)
        q = hm_beta*ups*exp(-hm_dE(i)/(kb_eV*T))/(hm_gl(i)*sqrt(T))
        qsum(hm_role(i)) = qsum(hm_role(i)) + q
     end do
     worst = max(worst, abs(r%q13/qsum(ROLE_Q13) - 1.0_wp), &
                 abs(r%qfeed/qsum(ROLE_FEED) - 1.0_wp), &
                 abs(r%q31a/qsum(ROLE_Q31A) - 1.0_wp), &
                 abs(r%q31b/qsum(ROLE_Q31B) - 1.0_wp), &
                 abs(r%q31s/qsum(ROLE_SINGLET) - 1.0_wp), &
                 abs(r%q33/qsum(ROLE_TRIPLET) - 1.0_wp))
     !--- q_31g: detailed balance of the Q13 row, g_1/g_3 = 1/3
     worst = max(worst, abs(r%q31g*exp(-hm_E23S/(kb_eV*T))/(3.0_wp*r%q13) &
                 - 1.0_wp/9.0_wp)*9.0_wp)
  end do
  write(*,'(a,es9.2)') '     rate sums and q_31g detailed balance at the nodes: worst ', worst
  if (worst > 1.0e-12_wp) ok = .false.
  r = hei_triplet_rates(1.0e4_wp)
  write(*,'(a,6es11.4)') '     1e4 K: q13 qfeed q31a q31b q31s q31g =', &
     r%q13, r%qfeed, r%q31a, r%q31b, r%q31s, r%q31g
  write(*,'(a,3es11.4)') '            q33 C_3 Q_31 =', r%q33, r%c3, r%pen

  !--- (b) alpha_3
  worst = 0.0_wp
  do it = 1, 5
     r = hei_triplet_rates(Ts(it))
     call hei_branch_fractions(Ts(it), f3s, f1s, f1p)
     worst = max(worst, abs(r%alpha3/(f3s*alphaB_HeII(Ts(it))) - 1.0_wp))
  end do
  write(*,'(a,es9.2)') ' (b) alpha_3 = f_3S alpha_B(He II): worst ', worst
  if (worst > 1.0e-14_wp) ok = .false.
  !--- against the Oklopcic & Hirata (2018) fit the module used to carry
  do it = 1, 5
     r = hei_triplet_rates(Ts(it))
     write(*,'(a,f8.0,a,es11.4,a,f7.4)') '     T =', Ts(it), '  alpha_3 =', &
        r%alpha3, '  / OH18 2.10e-13 T4^-0.778 =', &
        r%alpha3/(2.10e-13_wp*(Ts(it)/1.0e4_wp)**(-0.778_wp))
  end do

  !--- (c) closed form vs the 2x2 solve, over T, n_e, n_HI, Phi_3
  worst = 0.0_wp
  do it = 1, 5
     r = hei_triplet_rates(Ts(it))
     do ih = 0, 3
        ne = 10.0_wp**(1 + 2*ih);  nHI = 0.3_wp*ne;  nHeI = 0.05_wp*ne
        nHeII = 0.08_wp*ne;  phi3 = 1.0e-5_wp*ih
        n3 = hei_triplet_population(r, ne, nHI, nHeI, nHeII, phi3)
        D  = hei_triplet_destruction(r, ne, nHI, phi3)
        !--- unknowns (n_1, n_3):  n_1 + n_3 = n_HeI,
        !---   ne (q13 + qfeed) n_1 - D n_3 = - ne alpha_3 n_HeII
        a11 = 1.0_wp;  a12 = 1.0_wp;  b1 = nHeI
        a21 = ne*(r%q13 + r%qfeed);  a22 = -D;  b2 = -ne*r%alpha3*nHeII
        det = a11*a22 - a12*a21
        n3d = (a11*b2 - a21*b1)/det
        !--- the closed form takes n_1 = n_HeI; the difference is the
        !--- fraction n_3/n_HeI, which bounds the tolerance
        worst = max(worst, abs(n3/n3d - 1.0_wp)/max(n3/nHeI, 1.0e-16_wp))
     end do
  end do
  write(*,'(a,es9.2)') ' (c) closed-form n_3 vs 2x2 solve, |ratio-1|/(n_3/n_HeI): worst ', worst
  if (worst > 1.1_wp) ok = .false.

  !--- (d) energy ledger
  worst = 0.0_wp
  e3 = hm_E23S*ev2erg
  do it = 1, 5
     r = hei_triplet_rates(Ts(it))
     do ih = 0, 3
        ne = 10.0_wp**(1 + 2*ih);  nHI = 0.3_wp*ne;  nHeI = 0.05_wp*ne
        nHeII = 0.08_wp*ne;  phi3 = 2.0e-5_wp*ih;  heat3 = phi3*3.0_wp*ev2erg
        n3 = hei_triplet_population(r, ne, nHI, nHeI, nHeII, phi3)
        h  = hei_triplet_net_heating(r, ne, nHI, nHeI, n3, heat3)
        !--- photon energy emitted, from the rows' own energies
        pem = n3*hm_A31*e3
        do i = 1, hm_nrow
           q = hm_beta*hei_upsilon(i, Ts(it)) &
               *exp(-hm_dE(i)/(kb_eV*Ts(it)))/(hm_gl(i)*sqrt(Ts(it)))
           select case (hm_role(i))
           case (ROLE_FEED)               ! the triplet cascade to 2^3S
              pem = pem + ne*nHeI*q*(hm_dE(i) - hm_E23S)*ev2erg
           case (ROLE_Q31A, ROLE_Q31B, ROLE_SINGLET)   ! the singlet decays
              pem = pem + ne*n3*q*(hm_dE(i) + hm_E23S)*ev2erg
           case (ROLE_TRIPLET)            ! the cascade back to 2^3S
              pem = pem + ne*n3*q*hm_dE(i)*ev2erg
           end select
        end do
        stored = n3*((ne*r%c3 + phi3)*IP_HE_EV + nHI*r%pen*IP_H_EV)*ev2erg
        lhs = ne*r%alpha3*nHeII*e3 + n3*(heat3 + phi3*(IP_HE_EV - hm_E23S)*ev2erg)
        rhs = h + pem + stored
        worst = max(worst, abs(lhs - rhs)/(abs(lhs) + abs(h) + pem + stored))
     end do
  end do
  write(*,'(a,es9.2)') ' (d) energy ledger closure: worst relative ', worst
  if (worst > 1.0e-12_wp) ok = .false.

  !--- (e) low-density limit of the channel-4 weights
  call hei_phi3_assign([0.0_wp], [0.0_wp])
  worst = 0.0_wp
  do it = 1, 5
     ne = 1.0e-10_wp;  nHeII = 1.0e-11_wp   ! departure O(n_e q/A_31)
     call hei_branch_energies_coupled(1, Ts(it), ne, 0.0_wp, 0.0_wp, nHeII, &
                                      0.0_wp, 0.0_wp, c3s, c1p, c1s)
     call hei_branch_energies(Ts(it), 0.0_wp, 0.0_wp, e3s, e1p, e1s)
     rec = ne*nHeII*alphaB_HeII(Ts(it))
     !--- x_HI = x_HeI = 0: p_584 = 1, every 584 A photon ionizes H
     if (.not. (e3s > 0.0_wp .and. e1p > 0.0_wp .and. e1s > 0.0_wp)) ok = .false.
     worst = max(worst, abs(c3s/(rec*e3s) - 1.0_wp), &
                 abs(c1p/(rec*e1p) - 1.0_wp), abs(c1s/(rec*e1s) - 1.0_wp))
     !--- and with a 584 A split (0 < p_584 < 1); n_HI still 0 in n_3
     call hei_branch_energies_coupled(1, Ts(it), ne, 0.0_wp, 0.0_wp, nHeII, &
                                      0.5_wp, 0.01_wp, c3s, c1p, c1s)
     call hei_branch_energies(Ts(it), 0.5_wp, 0.01_wp, e3s, e1p, e1s)
     worst = max(worst, abs(c3s/(rec*e3s) - 1.0_wp), &
                 abs(c1p/(rec*e1p) - 1.0_wp), abs(c1s/(rec*e1s) - 1.0_wp))
  end do
  write(*,'(a,es9.2)') ' (e) diffuse He I weights, low-density limit: worst ', worst
  if (worst > 1.0e-10_wp) ok = .false.
  !--- and what the coupling does at 1e4 K, n_e = 100 and 1e4, in the gas
  T = 1.0e4_wp
  do ih = 1, 2
     ne = 10.0_wp**(2*ih);  nHeII = 0.1_wp*ne
     call hei_branch_energies_coupled(1, T, ne, 1.0e-3_wp*ne, 1.0e-3_wp*ne, &
          nHeII, 1.0e-3_wp, 1.0e-2_wp, c3s, c1p, c1s)
     call hei_branch_energies(T, 1.0e-3_wp, 1.0e-2_wp, e3s, e1p, e1s)
     rec = ne*nHeII*alphaB_HeII(T)
     write(*,'(a,es8.1,a,3f9.5)') '     n_e =', ne, &
        '  coupled/uncoupled (19.82 eV, 584 A, 2^1S) =', c3s/(rec*e3s), &
        c1p/(rec*e1p), c1s/(rec*e1s)
  end do

  !--- (f) solve_ion_cell with the coupling, partially neutral gas
  worst = 0.0_wp
  nH = 100.0_wp;  T = 8.0e3_wp
  gHe1 = 5.0e-11_wp;  gHe2 = 1.0e-13_wp
  call hei_phi3_assign([3.0e-5_wp], [3.0e-5_wp*4.0_wp*ev2erg])
  r = hei_triplet_rates(T)
  aH = alphaA_HII(T);  aHe2 = alphaA_HeII(T);  cH = ci_HI(T);  cHe1 = ci_HeI(T)
  do icx = 0, 1
     par%helium_charge_exchange = icx == 1
     do ih = 1, 4
        gH = GHs(ih)
        x = 0.5_wp;  h1 = 0.5_wp;  h2 = 0.4_wp
        call solve_ion_cell(gH, gHe1, gHe2, nH, T, .true., x, h1, h2, ne, 1)
        xHeIII = max(0.0_wp, 1.0_wp - h1 - h2)
        nHI = nH*x;  nHII = nH*(1.0_wp - x)
        nHeI = par%abund_He*nH*h1;  nHeII = par%abund_He*nH*h2
        n3 = hei_triplet_population(r, ne, nHI, nHeI, nHeII, 3.0e-5_wp)
        k21 = 0.0_wp;  k12 = 0.0_wp
        if (icx == 1) then
           block
             use ion_balance_mod, only : cx_HeII_H, cx_HeI_Hp, cx_HeIII_H, cx_HeII_Hp
             k21 = cx_HeII_H(T);  k12 = cx_HeI_Hp(T)
             !--- hydrogen with He-H exchange and Penning
             eh = abs(nHI*(gH + cH*ne + par%abund_He*nH*(h2*k21 + xHeIII*cx_HeIII_H(T)) &
                        + r%pen*n3) &
                      - nHII*(aH*ne + par%abund_He*nH*(h1*k12 + h2*cx_HeII_Hp(T)))) &
                  / (nHI*(gH + cH*ne + r%pen*n3))
           end block
        else
           eh = abs(nHI*(gH + cH*ne + r%pen*n3) - nHII*aH*ne) &
                / (nHI*(gH + cH*ne + r%pen*n3))
        end if
        !--- He0 <-> He+: ionization from the ground and from 2^3S against
        !--- recombination
        e1 = abs(nHeI*(gHe1 + cHe1*ne + k12*nHII) + n3*(ne*r%c3 + 3.0e-5_wp) &
                 - nHeII*(aHe2*ne + k21*nHI)) / (nHeII*(aHe2*ne + k21*nHI))
        u3 = (ne*r%c3 + 3.0e-5_wp)/hei_triplet_destruction(r, ne, nHI, 3.0e-5_wp)
        write(*,'(a,l2,a,es9.2,a,f8.5,a,es10.3,a,es10.3,a,2es9.1)') &
           '  He-H cx', icx == 1, '  Gamma_H =', gH, '  x_HII =', 1.0_wp - x, &
           '  Penning/Gamma_H =', r%pen*n3/gH, '  u =', u3, '  resid', eh, e1
        worst = max(worst, eh, e1)
     end do
  end do
  write(*,'(a,es9.2)') ' (f) H and He balances with the 2^3S terms at the solved state: worst ', worst
  if (worst > 1.0e-8_wp) ok = .false.

  if (ok) then
     write(*,'(a)') 'PASS'
  else
     write(*,'(a)') 'FAIL'
     call MPI_FINALIZE(ierr)
     stop 1
  end if
  call MPI_FINALIZE(ierr)

contains

  !--- the Maxwellian constant, formed here from CODATA 2018 rather than read
  real(kind=wp) function hm_beta_value() result(b)
    b = 6.62607015e-27_wp**2/((2.0_wp*pi*9.1093837015e-28_wp)**1.5_wp &
        *sqrt(1.380649e-16_wp))
  end function hm_beta_value

end program check_hei_metastable_coupling
