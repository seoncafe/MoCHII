program check_charge_exchange_heating
!---------------------------------------------------------------------------
! Gate for charge-exchange heating (par%charge_exchange_heating).
!
! A charge-exchange reaction moves one electron between two binding
! energies.  What the gas gains is the ionization-potential difference less
! the excitation of the product state that leaves as radiation:
!
!   X^(it+1) + H0 -> X^it + H+       dE_rec = [IP(X^it) - IP(H)] - E_rad
!   X^it + H+ -> X^(it+1) + H0       dE_ion = [IP(H) - IP(X^it)] - E_rad
!
! dE < 0 for an endothermic reaction.  The defects are the CXHEAT rows of
! data/atomic/element_*.txt (tools/fitting/make_element_data.py, which names
! the product states and their sources); species_mod::
! metal_charge_exchange_heating and ion_balance_mod::
! helium_charge_exchange_heating turn them into heat.
!
! Checks:
!  (a) census: every carried reaction (45: 7 ionization + 38 recombination,
!      all 11 elements active on the 77-stage roster; with the argument
!      'langevin', par%langevin_charge_exchange adds X^q + H for q = 5, 6,
!      22 more) has an entry in the product-state table below, and no entry
!      names a reaction that is not carried;
!  (b) detailed balance: for a pair that is one channel and its
!      detailed-balance reverse -- N and O (ground terms both ways), the
!      non-radiative part of C + H+ against C+ + H, and the non-radiative
!      He+ + H against He + H+ -- the heat of the two directions cancels at
!      the populations where they balance, to TOL_DB relative;
!  (c) each file defect equals the ionization-potential difference minus the
!      product excitation typed below (transcribed independently from the
!      CHIANTI observed level energies, or the Kingdon & Ferland (1999, ApJ
!      516, L107) Table 1 mean energy defect where the product is not
!      identified), to TOL_DE eV;
!  (e) metal_charge_exchange_heating equals the sum of the reactions formed
!      here from species_fractions, species_cx_coefficients and
!      species_cx_kinetic_energy, to TOL_SUM relative -- the stage indices
!      and the radiative split of the code are the ones checked in (b);
!  (d) prints the heat of each reaction at 1e4 K and in two cells of the
!      converged Lexington HII40 run (interior and ionization front), against
!      the photoheating there.
!
! The pairs S, Mg, Si and Fe (transition 1) are NOT tested in (b), and why is
! the point: their exothermic direction ends in a product term that radiates
! (S+ 2P/2D, Mg+ 3p, Si+ 4P, Fe+ z4F/z4D), whereas the carried endothermic
! direction (Huang et al. 2023) is the detailed-balance reverse of the total
! rate from the ground term.  That pair is not one channel, so its heat does
! not cancel; the gate prints the imbalance.
!
! Build:  see tests/charge_exchange/run_check.sh
!---------------------------------------------------------------------------
  use define
  use species_mod,     only : species_setup, species_fractions, &
                              species_cx_coefficients, &
                              species_cx_kinetic_energy, &
                              metal_charge_exchange_heating, &
                              species_set_photoionization_rates, &
                              n_elements, elem_name, elem_nstage, elem_abund, &
                              elem_eth
  use ion_balance_mod, only : helium_charge_exchange_heating, &
                              nonradiative_HeII_H, cx_HeI_Hp, cx_HeII_H, &
                              cx_HeIII_H, cx_HeII_Hp
  use mpi
  implicit none

  real(kind=wp) :: worst_mg, e_th, t_ext, rr
  real(kind=wp), parameter :: KB_EV_GATE = 8.617333262e-5_wp  ! CODATA 2018 [eV/K]
  real(kind=wp), parameter :: IP_H     = 13.598434_wp      ! NIST ASD
  real(kind=wp), parameter :: IP_HE    = 24.587387_wp
  real(kind=wp), parameter :: HC_EV_CM = 1.239841984e-4_wp ! CODATA 2018
  real(kind=wp), parameter :: TOL_DB   = 1.0e-12_wp
  real(kind=wp), parameter :: TOL_DE   = 1.0e-6_wp
  real(kind=wp), parameter :: TOL_SUM  = 1.0e-13_wp
  integer,       parameter :: NCXI_EXPECT = 7, NCXR_EXPECT = 38, NLANGEVIN = 22
  !--- Langevin rows (q >= 5): KF99 Sec. 2, a mean avoided crossing at
  !--- 9.5 bohr gives the exit channel X^(q-1)+ + H+ a Coulomb energy
  !--- (q-1) E_h / 9.5 = 2.86 (q-1) eV.  KF99's text writes "2.86q eV" for
  !--- q > 4, but every avoided-crossing entry of its Table 1 is 2.86 (q-1)
  !--- (5.73 at q = 3, 8.60 at q = 4), and Cloudy c25.00 applies
  !--- 2.86*(double)ion to the recombining ion of charge ion+1
  !--- (atmdat_char_tran.cpp, ChargTranSumHeat), i.e. 2.86 (q-1).
  real(kind=wp), parameter :: LANGEVIN_DE_EV = 2.86_wp
  integer,       parameter :: MST = MAX_ION_STAGES

  !--- product-state table, one row per carried reaction:
  !--- element, transition, direction (1 ion, 2 rec), kind, value
  !---   kind 0 ground-term products          E_rad = 0
  !---   kind 1 radiative charge transfer     E_rad = dIP
  !---   kind 2 product term                  E_rad = value [cm^-1]
  !---   kind 3 KF99 Table 1 mean defect      dE = value [eV], product
  !---                                        not identified
  ! KF96 gives total rates only; the mean energy defect of a reaction whose
  ! product is not resolved is KF99 Table 1 (read from the publisher PDF,
  ! md/atomic_stage_expansion/references/Kingdon_Ferland_1999_ApJ_516_L107.pdf).
  ! Kind 3 checks the file against that table; it does not verify an
  ! unresolved product state.
  !--- Term energies: CHIANTI observed energies, (2J+1)-weighted over the J
  !--- levels (O II 2s2p4 4P, N II 2s2p3 3D, Mg II 3p 2P, Si II 3s3p2 4P,
  !--- Si IV 4s 2S); S + H+ weighted over S+ 4S/2D/2P with Zhao et al.
  !--- (2005) Table III at 1e4 K (8.90e-16, 1.19e-12, 1.63e-11); Fe + H+
  !--- the z4D (g = 20) and z4F (g = 28) mean, 44743.6587 and 44714.0716.
  !--- rows 30-45: the q = 3, 4 recombinations the 77-stage roster added,
  !--- all KF99 Table 1 means (C4+ 11.30, N3+ 6.40, N4+ 11.00, O3+ 5.00,
  !--- O4+ 8.47, Ne3+ 5.82, Ne4+ 8.60, Mg3+ 5.73, Mg4+ 8.60, S4+ 8.60,
  !--- Ar4+ 8.60, Fe4+ 10.00 -- printed "Fe12" twice, the second in the
  !--- q = 3 column -- Cl3+ 8.00, Cl4+ 9.00, Ca3+ 6.60, Ca4+ 9.90)
  integer, parameter :: NTAB = 45
  character(len=9), parameter :: KIND_NAME(0:3) = &
     ['ground   ', 'radiative', 'term     ', 'KF99 dE  ']
  character(len=2), parameter :: TEL(NTAB) = [character(len=2) :: &
     'c ','c ','c ','c ', 'n ','n ','n ', 'o ','o ','o ', 'ne', &
     's ','s ','s ','s ', 'ar','ar', 'mg','mg','mg', &
     'fe','fe','fe','fe', 'si','si','si','si','si', &
     'c ', 'n ','n ', 'o ','o ', 'ne','ne', 'mg','mg', 's ', 'ar', 'fe', &
     'cl','cl', 'ca','ca']
  integer, parameter :: TIT(NTAB) = [1,1,2,3, 1,1,2, 1,1,2, 2, &
     1,1,2,3, 2,3, 1,1,2, 1,1,2,3, 1,1,2,3,4, &
     4, 3,4, 3,4, 3,4, 3,4, 4, 4, 4, 3,4, 3,4]
  integer, parameter :: TDIR(NTAB) = [1,2,2,2, 1,2,2, 1,2,2, 2, &
     1,2,2,2, 2,2, 1,2,2, 1,2,2,2, 1,2,2,2,2, &
     2, 2,2, 2,2, 2,2, 2,2, 2, 2, 2, 2,2, 2,2]
  integer, parameter :: TKIND(NTAB) = [0,0,3,3, 0,0,2, 0,0,2, 1, &
     2,0,1,3, 1,3, 2,0,0, 2,0,0,3, 2,0,0,3,2, &
     3, 3,3, 3,3, 3,3, 3,3, 3, 3, 3, 3,3, 3,3]
  real(kind=wp), parameter :: S_BRANCH_CM = &
     (8.90e-16_wp*0.0_wp + 1.19e-12_wp*14872.0140_wp + 1.63e-11_wp*24555.9693_wp) &
     / (8.90e-16_wp + 1.19e-12_wp + 1.63e-11_wp)
  real(kind=wp), parameter :: FE_Z4_CM = (20.0_wp*44743.6587_wp + 28.0_wp*44714.0716_wp)/48.0_wp
  !--- Mg + H+ ends in Mg+ 3p 2P (the 2 1Sigma channel) or Mg+ 4s 2S (4 1Sigma),
  !--- Allan et al. (1988) Table 1 [1e-10 cm^3 s^-1], transcribed here
  !--- independently of make_element_data.py; the CXHEAT row holds 1e4 K
  integer, parameter :: NMG = 16
  real(kind=wp), parameter :: MG_T(NMG) = [5000., 6000., 8000., 10000., 12000., &
     14000., 15000., 16000., 18000., 20000., 22000., 24000., 26000., 28000., 30000., 32000.]
  real(kind=wp), parameter :: MG_K3P(NMG) = [0.2549, 0.3430, 0.5055, 0.6419, 0.7582, &
     0.8632, 0.9138, 0.9641, 1.066, 1.170, 1.278, 1.390, 1.504, 1.620, 1.736, 1.850]
  real(kind=wp), parameter :: MG_K4S(NMG) = [0.07057, 0.1856, 0.6125, 1.236, 1.955, &
     2.696, 3.060, 3.415, 4.093, 4.723, 5.307, 5.847, 6.350, 6.819, 7.259, 7.673]
  real(kind=wp), parameter :: MG_3P_CM = 35730.3657_wp, MG_4S_CM = 69804.898_wp
  real(kind=wp), parameter :: MG_ALLAN_1E4_CM = (0.6419_wp*MG_3P_CM + 1.236_wp*MG_4S_CM) &
                                                /(0.6419_wp + 1.236_wp)
  real(kind=wp), parameter :: TVAL(NTAB) = [0.0_wp, 0.0_wp, 4.01_wp, 5.73_wp, &
     0.0_wp, 0.0_wp, 92244.4867_wp, &
     0.0_wp, 0.0_wp, 119932.5500_wp, &
     0.0_wp, &
     S_BRANCH_CM, 0.0_wp, 0.0_wp, 5.73_wp, &
     0.0_wp, 5.73_wp, &
     MG_ALLAN_1E4_CM, 0.0_wp, 0.0_wp, &
     FE_Z4_CM, 0.0_wp, 0.0_wp, 6.3_wp, &
     43002.2100_wp, 0.0_wp, 0.0_wp, 4.23_wp, 193978.8910_wp, &
     11.30_wp, 6.40_wp, 11.00_wp, 5.00_wp, 8.47_wp, 5.82_wp, 8.60_wp, &
     5.73_wp, 8.60_wp, 8.60_wp, 8.60_wp, 10.00_wp, 8.00_wp, 9.00_wp, &
     6.60_wp, 9.90_wp]
  !--- detailed-balance classes of the transition-1 pairs: 1 = one channel
  !--- and its reverse, 2 = radiative ionization with its non-radiative part
  !--- the reverse of the recombination (g_split), 0 = not one channel.
  character(len=2), parameter :: PAIR_EL(7) = ['c ','n ','o ','s ','mg','fe','si']
  integer,          parameter :: PAIR_CLASS(7) = [2, 1, 1, 0, 0, 0, 0]

  real(kind=wp) :: k_i, k_r, de_i, de_r, gs, dip, erad, dexp, kT, frac(MST)
  real(kind=wp) :: T, nH, ne, nHI, nHII, hi, hr, net, gross, knr, s, hsum, hfun
  real(kind=wp) :: nHeI, nHeII, heat_ph, tot_cx
  integer :: ie, it, idir, j, ierr, ncxi, ncxr, found, k, ist, ncxr_required
  integer :: kind_found
  logical :: ok, carried, langevin
  real(kind=wp) :: value_found
  character(len=2) :: nm
  character(len=16) :: arg
  real(kind=wp), parameter :: TT(3) = [3.0e3_wp, 1.0e4_wp, 2.0e4_wp]

  call MPI_INIT(ierr);  mpar%p_rank = 0
  call get_command_argument(1, arg)
  langevin = trim(arg) == 'langevin'
  par%langevin_charge_exchange = langevin
  par%atomic_dir = 'data/atomic'
  par%use_metals = .true.
  par%charge_exchange_heating = .true.
  !--- the n-level (T, n_e) cooling tables play no part here; 'low_density'
  !--- keeps species_setup from building them for every stage
  par%cooling_model = 'low_density'
  par%abund_C  = 2.2e-4_wp;  par%abund_N  = 4.0e-5_wp
  par%abund_O  = 3.3e-4_wp;  par%abund_Ne = 5.0e-5_wp
  par%abund_S  = 9.0e-6_wp;  par%abund_Ar = 3.0e-6_wp
  par%abund_Mg = 3.0e-6_wp;  par%abund_Fe = 3.0e-6_wp
  par%abund_Si = 3.0e-6_wp;  par%abund_Cl = 3.0e-7_wp
  par%abund_Ca = 2.0e-8_wp
  call species_setup(1)
  ok = .true.

  !=======================================================================
  write(*,'(/,a)') '=== (a) census and (c) energy defects [eV] ==='
  write(*,'(a4,a4,a5,a12,a12,a12,a12,a10)') 'el', 'it', 'dir', 'k(1e4 K)', &
     'dIP', 'E_rad', 'dE file', 'kind'
  ncxi = 0;  ncxr = 0
  do ie = 1, n_elements
     nm = elem_name(ie)
     do it = 1, elem_nstage(ie) - 1
        call species_cx_coefficients(ie, it, 1.0e4_wp, k_i, k_r)
        call species_cx_kinetic_energy(ie, it, de_i, de_r, gs)
        do idir = 1, 2
           carried = merge(k_i, k_r, idir == 1) > 0.0_wp
           if (carried) then
              if (idir == 1) ncxi = ncxi + 1
              if (idir == 2) ncxr = ncxr + 1
           end if
           found = 0
           do j = 1, NTAB
              if (trim(TEL(j)) == trim(nm) .and. TIT(j) == it .and. &
                  TDIR(j) == idir) found = j
           end do
           if (found > 0) then
              kind_found = TKIND(found);  value_found = TVAL(found)
           end if
           if (langevin .and. it >= 5 .and. idir == 2) then
              found = NTAB + 1
              kind_found = 3;  value_found = LANGEVIN_DE_EV*real(it - 1, wp)
           end if
           if (carried .neqv. (found > 0)) then
              write(*,'(a,a,a,i0,a,i0,a,l2,a,l2)') '  FAIL (a): ', nm, ' it=', it, &
                 ' dir=', idir, ' carried=', carried, ' in table=', found > 0
              ok = .false.
              cycle
           end if
           if (found == 0) then
              !--- no rate: the file defect must be zero
              if (merge(de_i, de_r, idir == 1) /= 0.0_wp) then
                 write(*,'(a,a,i3,i3)') '  FAIL (a): defect without a rate ', nm, it, idir
                 ok = .false.
              end if
              cycle
           end if
           !--- energy released when the ELECTRON moves as written
           dip = elem_eth(ie, it) - IP_H
           if (idir == 1) dip = -dip
           select case (kind_found)
           case (0);  erad = 0.0_wp
           case (1);  erad = dip
           case (2);  erad = value_found*HC_EV_CM
           case (3);  erad = dip - value_found
           end select
           dexp = dip - erad
           if (abs(merge(de_i, de_r, idir == 1) - dexp) > TOL_DE) then
              write(*,'(a,a,i3,i3,2f12.6)') '  FAIL (c): ', nm, it, idir, &
                 merge(de_i, de_r, idir == 1), dexp
              ok = .false.
           end if
           write(*,'(a4,i4,a5,es12.3,f12.6,f12.6,f12.6,a10)') nm, it, &
              merge(' ion',' rec', idir == 1), merge(k_i, k_r, idir == 1), &
              dip, erad, merge(de_i, de_r, idir == 1), trim(KIND_NAME(kind_found))
        end do
     end do
  end do
  !--- (c') the temperature-branched Mg + H+ defect at Allan's nodes
  do ie = 1, n_elements
     if (trim(elem_name(ie)) /= 'mg') cycle
     worst_mg = 0.0_wp
     do j = 1, NMG
        call species_cx_kinetic_energy(ie, 1, de_i, de_r, gs, MG_T(j))
        dexp = -(elem_eth(ie, 1) - IP_H) - (MG_K3P(j)*MG_3P_CM + MG_K4S(j)*MG_4S_CM) &
               /(MG_K3P(j) + MG_K4S(j))*HC_EV_CM
        worst_mg = max(worst_mg, abs(de_i - dexp))
     end do
     write(*,'(a,es10.2,a)') " (c') Mg + H+ defect at Allan et al. (1988) nodes: worst ", &
        worst_mg, ' eV'
     if (worst_mg > TOL_DE) then
        write(*,'(a)') "  FAIL (c'): Mg + H+ temperature table"
        ok = .false.
     end if
     !--- (c'') below Allan's first node (5000 K) the 4s share is carried
     !--- down by its threshold Boltzmann factor, the 3p rate held
     !--- (make_element_data.py BRANCH_T_EXTEND): at 2000, 3000, 4000 K
     worst_mg = 0.0_wp
     e_th = MG_4S_CM*HC_EV_CM - (IP_H - elem_eth(ie, 1))
     do j = 2, 4
        t_ext = 1000.0_wp*j
        rr = MG_K4S(1)/MG_K3P(1)*exp(-e_th/(KB_EV_GATE*t_ext) + e_th/(KB_EV_GATE*MG_T(1)))
        call species_cx_kinetic_energy(ie, 1, de_i, de_r, gs, t_ext)
        dexp = -(elem_eth(ie, 1) - IP_H) - (MG_3P_CM + rr*MG_4S_CM)/(1.0_wp + rr)*HC_EV_CM
        worst_mg = max(worst_mg, abs(de_i - dexp))
     end do
     write(*,'(a,es10.2,a)') " (c'') Mg + H+ defect below 5000 K (threshold extension): worst ", &
        worst_mg, ' eV'
     if (worst_mg > TOL_DE) then
        write(*,'(a)') "  FAIL (c''): Mg + H+ extension below the Allan table"
        ok = .false.
     end if
  end do
  ncxr_required = NCXR_EXPECT
  if (langevin) ncxr_required = NCXR_EXPECT + NLANGEVIN
  write(*,'(a,i0,a,i0,a,i0,a,i0,a,l2)') ' carried: ', ncxi, ' ionization (expect ', &
     NCXI_EXPECT, '), ', ncxr, ' recombination (expect ', ncxr_required, &
     '); langevin_charge_exchange =', langevin
  if (ncxi /= NCXI_EXPECT .or. ncxr /= ncxr_required) then
     write(*,'(a)') '  FAIL (a): reaction census'
     ok = .false.
  end if

  !=======================================================================
  ! (b) detailed balance.  The two stages of transition 1 are given the
  ! populations at which the two directions balance, n(X+)/n(X) =
  ! k_ion n_HII / (k_rec n_HI) (for carbon, whose ionization rate is mostly
  ! radiative, the balance of the non-radiative channel, g exp(dE/kT)
  ! n_HII/n_HI), and a one-channel pair must then deposit nothing.  The
  ! rate coefficients and defects are the code's; that the code applies them
  ! to the right stages is check (e).
  write(*,'(/,a)') '=== (b) heat of each transition-1 pair where the two directions balance ==='
  write(*,'(a4,a10,a14,a14,a12,a8)') 'el', 'T [K]', 'heat ion', 'heat rec', &
     'net/gross', 'class'
  nHI = 50.0_wp;  nHII = 50.0_wp
  do ie = 1, n_elements
     nm = elem_name(ie)
     k = 0
     do j = 1, size(PAIR_EL)
        if (trim(PAIR_EL(j)) == trim(nm)) k = j
     end do
     if (k == 0) cycle
     do j = 1, size(TT)
        T = TT(j);  kT = kboltz_cgs*T
        call species_cx_coefficients(ie, 1, T, k_i, k_r)
        call species_cx_kinetic_energy(ie, 1, de_i, de_r, gs)
        knr = k_i
        frac = 0.0_wp;  frac(1) = 1.0_wp
        frac(2) = k_i*nHII/(k_r*nHI)
        if (gs > 0.0_wp) then
           knr = min(k_i, gs*k_r*exp(min(de_i*ev2erg/kT, 600.0_wp)))
           frac(2) = gs*exp(de_i*ev2erg/kT)*nHII/nHI
        end if
        hi = frac(1)*nHII*knr*de_i*ev2erg
        hr = frac(2)*nHI*k_r*de_r*ev2erg
        net = hi + hr;  gross = abs(hi) + abs(hr)
        write(*,'(a4,f10.0,es14.5,es14.5,es12.3,i8)') nm, T, hi, hr, &
           net/max(gross, tiny(1.0_wp)), PAIR_CLASS(k)
        if (PAIR_CLASS(k) > 0 .and. abs(net) > TOL_DB*gross) then
           write(*,'(a,a,f8.0)') '  FAIL (b): pair does not cancel: ', nm, T
           ok = .false.
        end if
        if (PAIR_CLASS(k) == 2 .and. knr >= k_i) &
           write(*,'(a)') '    (carbon: the non-radiative reverse reaches the whole '// &
              'carried C + H+ rate here; clipped)'
     end do
  end do
  write(*,'(a)') ' class 1: one channel and its reverse; 2: non-radiative part of a '// &
     'radiative ionization;'
  write(*,'(a)') ' class 0: exothermic direction ends in a radiating term -- '// &
     'not one channel, not required to cancel'
  !--- He+ + H <-> He + H+ (non-radiative channel): populations where the two
  !--- balance, n(He+) = k(He + H+) n_HII / (k_nr n_HI) per n(He).  The
  !--- He+ + H+ -> He++ + H term of the same function is not part of this
  !--- pair (its forward, He++ + H, gives the gas nothing) and is taken out.
  do j = 1, size(TT)
     T = TT(j)
     nHeI = 1.0_wp
     nHeII = cx_HeI_Hp(T)*nHII/(nonradiative_HeII_H(T)*nHI)
     net = helium_charge_exchange_heating(T, nHI, nHII, nHeI, nHeII) &
         + cx_HeII_Hp(T)*nHeII*nHII*(54.417765_wp - IP_H)*ev2erg
     gross = 2.0_wp*nonradiative_HeII_H(T)*nHeII*nHI*(IP_HE - IP_H)*ev2erg
     write(*,'(a4,f10.0,es14.5,a14,es12.3,i8)') 'he', T, gross/2.0_wp, '', &
        net/gross, 1
     if (abs(net) > TOL_DB*gross) then
        write(*,'(a,f8.0)') '  FAIL (b): He+ + H / He + H+ does not cancel at ', T
        ok = .false.
     end if
  end do

  !=======================================================================
  ! (e) and (d): the code's sum against the reactions formed here, and the
  ! heat of each reaction in two cells of the converged Lexington HII40 run
  ! (results/continuous_energy/m5_converged/hii40_converged_rates.h5, leaves
  ! 50652 at r = 3.50 pc and 21790 at r = 4.70 pc; n_H = 100).  Only C, N,
  ! O, Ne, S have photoionization rates there (the other elements are not
  ! in that run), so (d) prints those five; (e) sums all eleven.
  do ist = 1, 3
     select case (ist)
     case (1)
        T = 1.0e4_wp
        write(*,'(/,a)') '=== (d) heating rate coefficient k dE of each reaction '// &
           'at 1e4 K [erg cm^3 s^-1] ==='
        do ie = 1, n_elements
           do it = 1, elem_nstage(ie) - 1
              call species_cx_coefficients(ie, it, T, k_i, k_r)
              call species_cx_kinetic_energy(ie, it, de_i, de_r, gs)
              if (gs > 0.0_wp) k_i = min(k_i, gs*k_r*exp(de_i*ev2erg/(kboltz_cgs*T)))
              if (k_i > 0.0_wp) write(*,'(a4,i3,a,es12.4)') elem_name(ie), it, &
                 ' ion ', k_i*de_i*ev2erg
              if (k_r > 0.0_wp) write(*,'(a4,i3,a,es12.4)') elem_name(ie), it, &
                 ' rec ', k_r*de_r*ev2erg
           end do
        end do
        write(*,'(a,es12.4)') '  he    nonradiative He+ + H ', &
           nonradiative_HeII_H(T)*(IP_HE - IP_H)*ev2erg
        write(*,'(a,es12.4)') '  he    He + H+              ', &
           -cx_HeI_Hp(T)*(IP_HE - IP_H)*ev2erg
        write(*,'(a,es12.4,a)') '  he    radiative He+ + H    ', 0.0_wp, &
           ' (of k = '//trim(real_str(cx_HeII_H(T) - nonradiative_HeII_H(T)))//')'
        write(*,'(a,es12.4,a)') '  he    He++ + H             ', 0.0_wp, &
           ' (of k = '//trim(real_str(cx_HeIII_H(T)))//')'
        write(*,'(a,es12.4)') '  he    He+ + H+             ', &
           -cx_HeII_Hp(T)*(54.417765_wp - IP_H)*ev2erg
        cycle
     case (2)
        T = 7803.476189761532_wp;  nH = 100.0_wp
        ne = 109.97796719961642_wp
        nHI = nH*7.309034716569132e-4_wp;  nHII = nH - nHI
        nHeI = 0.1_wp*nH*6.340759600166193e-3_wp
        nHeII = 0.1_wp*nH*0.9936592403998338_wp
        heat_ph = nHI*3.168978905701486e-19_wp + nHeI*4.41464184160417e-20_wp
        call set_leaf_gamma( &
           [2.5341205379496697e-07_wp, 7.3821155869454786e-09_wp, 4.633155589014981e-12_wp], &
           [2.1238966515413195e-07_wp, 3.1612554742944816e-09_wp], &
           [1.5489381604815775e-07_wp, 1.0520477016259716e-09_wp], &
           [2.763562146113921e-08_wp, 2.4626475156533203e-10_wp], &
           [7.491640171373708e-07_wp, 1.1627456451703479e-08_wp, 6.109035151600747e-11_wp])
        write(*,'(/,a)') '=== (d) HII40 interior cell (r = 3.50 pc, T = 7803 K, '// &
           'x_HI = 7.3e-4) ==='
     case (3)
        T = 8795.943545868564_wp;  nH = 100.0_wp
        ne = 67.83091840073897_wp
        nHI = nH*0.32215948165262676_wp;  nHII = nH - nHI
        nHeI = 0.1_wp*nH*0.9999999681438861_wp
        nHeII = 0.1_wp*nH*3.185611394960725e-8_wp
        heat_ph = nHI*4.627229060275799e-22_wp + nHeI*5.2520772469504154e-33_wp
        call set_leaf_gamma( &
           [2.5814216453604854e-10_wp, 3.6075812511592276e-12_wp, 9.243074264901394e-23_wp], &
           [2.821545257317047e-10_wp, 2.2468400281141265e-22_wp], &
           [2.348986910372561e-10_wp, 4.01492402211405e-22_wp], &
           [7.785278460183454e-11_wp, 5.415582458776014e-22_wp], &
           [7.054861563008813e-10_wp, 3.875191084358128e-11_wp, 5.575647501877496e-23_wp])
        write(*,'(/,a)') '=== (d) HII40 ionization-front cell (r = 4.70 pc, '// &
           'T = 8796 K, x_HI = 0.32) ==='
     end select
     write(*,'(a4,a4,a5,a14,a12)') 'el', 'it', 'dir', 'heat', '/photoheat'
     hsum = 0.0_wp;  tot_cx = 0.0_wp
     do ie = 1, n_elements
        call species_fractions(ie, 1, T, ne, nHI, nHII, frac)
        do it = 1, elem_nstage(ie) - 1
           call species_cx_coefficients(ie, it, T, k_i, k_r, ne)
           call species_cx_kinetic_energy(ie, it, de_i, de_r, gs, T, ne)
           if (gs > 0.0_wp) &
              k_i = min(k_i, gs*k_r*exp(min(de_i*ev2erg/(kboltz_cgs*T), 600.0_wp)))
           hi = elem_abund(ie)*nH*frac(it)*nHII*k_i*de_i*ev2erg
           hr = elem_abund(ie)*nH*frac(it+1)*nHI*k_r*de_r*ev2erg
           hsum = hsum + hi + hr
           if (.not. lexington(elem_name(ie))) cycle
           tot_cx = tot_cx + hi + hr
           if (k_i > 0.0_wp) write(*,'(a4,i4,a5,es14.4,es12.3)') elem_name(ie), &
              it, ' ion', hi, hi/heat_ph
           if (k_r > 0.0_wp) write(*,'(a4,i4,a5,es14.4,es12.3)') elem_name(ie), &
              it, ' rec', hr, hr/heat_ph
        end do
     end do
     hi = helium_charge_exchange_heating(T, nHI, nHII, nHeI, nHeII)
     write(*,'(a4,a4,a5,es14.4,es12.3)') 'he', '', '', hi, hi/heat_ph
     write(*,'(a,es12.4,a,es12.4,a,es12.4)') ' metals (C N O Ne S) ', tot_cx, &
        '   + He ', tot_cx + hi, '   share of H+He photoheating ', &
        (tot_cx + hi)/heat_ph
     !--- (e) the code's sum over all eleven elements
     hfun = metal_charge_exchange_heating(1, T, nH, ne, nHI, nHII)
     write(*,'(a,es16.8,a,es16.8,a,es10.2)') ' (e) metal_charge_exchange_heating ', &
        hfun, '  reaction sum ', hsum, '  rel ', abs(hfun - hsum)/max(abs(hsum), tiny(1.0_wp))
     if (abs(hfun - hsum) > TOL_SUM*max(abs(hsum), tiny(1.0_wp))) then
        write(*,'(a)') '  FAIL (e): code sum differs from the reactions'
        ok = .false.
     end if
  end do

  if (ok) then
     write(*,'(/,a)') ' GATE (charge-exchange heating): PASS'
  else
     write(*,'(/,a)') ' GATE (charge-exchange heating): FAIL'
     call MPI_FINALIZE(ierr);  stop 1
  end if
  call MPI_FINALIZE(ierr)

contains

  logical function lexington(name)
    character(len=*), intent(in) :: name
    lexington = any(trim(name) == ['c ','n ','o ','ne','s '])
  end function lexington

  character(len=16) function real_str(x)
    real(kind=wp), intent(in) :: x
    write(real_str,'(es10.3)') x
  end function real_str

  !--- metal photoionization rates of the Lexington elements in the one leaf;
  !--- the other elements get none
  subroutine set_leaf_gamma(gc, gn, go, gne, gs_)
    real(kind=wp), intent(in) :: gc(3), gn(2), go(2), gne(2), gs_(3)
    real(kind=wp) :: gg(MST-1,1)
    integer :: je, nt
    do je = 1, n_elements
       nt = elem_nstage(je) - 1
       gg = 0.0_wp
       select case (trim(elem_name(je)))
       case ('c');   gg(1:3,1) = gc
       case ('n');   gg(1:2,1) = gn
       case ('o');   gg(1:2,1) = go
       case ('ne');  gg(1:2,1) = gne
       case ('s');   gg(1:3,1) = gs_
       end select
       call species_set_photoionization_rates(je, gg(1:nt,:))
    end do
  end subroutine set_leaf_gamma

end program check_charge_exchange_heating
