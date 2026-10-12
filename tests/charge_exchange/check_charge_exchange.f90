program check_charge_exchange
!---------------------------------------------------------------------------
! Gate for two-way charge exchange (par%charge_exchange = 'two_way'):
! every reaction
!
!     X^i + H^+  ->  X^(i+1) + H^0        (coefficient k_CXI, rate R_cx_in)
!     X^(i+1) + H^0  ->  X^i + H^+        (coefficient k_CXR, rate R_cx_out)
!
! moves one unit of charge between a metal and hydrogen, so its volumetric
! rate must be ONE number: whatever leaves hydrogen must arrive at the metal.
! MoCHII forms the two sides in two places — the cascade's own side in
! species_fractions (which re-evaluates the coefficient from the registry and
! runs its own stage chain) and hydrogen's side in species_metal_sums (the
! cached coefficients and the cached chain of the n_e fixed point).  This gate
! forms both, reaction by reaction, on a converged cell state and requires
! them to agree to double-precision round-off.  A stage-index error in either
! sum shows up here as an O(1) disagreement on the reactions it touches.
!
! Checks:
!   1. reaction census: with all 11 elements active on the 77-stage roster
!      (ELEMENT .. 7, stages I-VII) the registry carries exactly the
!      reactions listed in CX_EL/CX_IT/CX_DIR below -- 7 CXI + 38 CXR, or
!      7 + 60 with par%langevin_charge_exchange (run with the argument
!      'langevin'), which adds X^q + H for q = 5, 6 of every element.  The
!      list is the set of nonzero CXI/CXR rows of data/atomic/element_*.txt
!      read off the files and checked against the sources they cite (KF96
!      Tables 1-2, Barragan et al. 2006, Huang et al. 2023, Zhao et al. 2005,
!      Ferland et al. 1997); it is compared reaction by reaction, so a row
!      gained or lost anywhere fails, not only a changed count;
!   2. pairwise conservation: for every reaction,
!         k_CXI n(X^i)  n_HII   from hydrogen's side
!      == k_CXI n(X^i)  n_HII   from the cascade's side,   to TOL_PAIR;
!      likewise k_CXR n(X^(i+1)) n_HI;
!   3. the balance
!         n_HI (Gamma_H + C_H n_e + R_cx_out) = (alpha_H n_e + R_cx_in) n_HII
!      holds (a) on the state solve_ion_cell converged to, to TOL_STATE, and
!      (b) exactly, to TOL_BALANCE, when x_HI is formed here from the closed
!      form the elimination of n_HII = n_H (1 - x_HI) gives;
!   4. the denominator of that closed form: R_cx_in appears in the denominator
!      as well as the numerator, for the same reason alpha_H n_e does.
!      Dropping it there must break check 3b — the gate asserts that the
!      truncated denominator FAILS, so the algebra cannot silently regress to
!      the wrong form;
!   5. the recombination rates of transitions 2-4 taken from Kingdon &
!      Ferland (1996, ApJS 106, 205) Tables 1 and 2 reproduce the published
!      fit k = a 1e-9 t4^b [1 + c exp(d t4)], t4 = T/1e4 K, held at the end
!      of each published temperature interval (the files' rule), with the
!      parameters typed below from the paper and cross-read against Cloudy
!      c25.00 data/ctrecombdata.dat.  Two rows are not the printed KF96
!      fits: N4+ + H and Si2+ + H take the updated KF96 database distributed
!      with Cloudy (Ferland et al. 2013, RMxAA 49, 137, Sec. 2.1;
!      ctrecombdata.dat: N4+ 2.95, 0.55, -0.39, -1.07 on 1e3-1e6 K, Si2+
!      6.77, 0.0736, -0.43, -0.11 on 500-1e5 K), the author's choice of
!      2026-10-11, as the element files do; the other 27 rows equal both.  The files tabulate these rates (form 8), so the
!      tolerance TOL_KF96 is the log-log interpolation between their nodes;
!   6. with the argument 'langevin', the q = 5, 6 rows are the Ferland et
!      al. (1997, ApJ 481, L115) eq. (2) Langevin rate with the reduced
!      mass, k = 2 pi q e sqrt(alpha_H a0^3 / mu), alpha_H = 4.5 (atomic
!      units), independent of T; without it every q >= 5 rate is zero.
!
! TOL_STATE is looser than TOL_BALANCE because solve_ion_cell exits its damped
! n_e fixed point at a relative 1e-10 and takes the charge-exchange sums from
! the previous iterate: in a nearly neutral cell the metal stage fractions are
! sharply sensitive to n_HI/n_HII, so that lag leaves a residual near 1e-6
! there against 1e-11 in the ionized states.  TOL_BALANCE is the round-off
! criterion on the algebra itself, raised by the round-off FLOOR of the metric:
! the residual is formed with n_HII = n_H (1 - x_HI), whose relative precision
! is eps/(1 - x_HI), so it degrades as the gas becomes neutral.
!
! The cell states below carry no metal photoionization (the leaf Gamma blocks
! of a bare species_setup are zero), so the cascade is driven by collisional
! ionization, recombination and charge exchange.  That is the regime in which
! the charge-exchange terms are largest, and the identity under test holds at
! any state: it is a property of how the two sides are formed, not of the
! ionizing field.
!
! Build:  see tests/charge_exchange/run_check.sh
!---------------------------------------------------------------------------
  use define
  use species_mod,     only : species_setup, species_fractions, &
                              species_ne_prepare, species_metal_sums, &
                              species_cx_coefficients, &
                              n_elements, elem_name, elem_nstage, elem_abund
  use ion_balance_mod, only : solve_ion_cell
  use recomb_mod,      only : alphaA_HII, alphaB_HII, ci_HI
  use mpi
  implicit none

  !--- Shared capacities of the reaction-resolved species output.
  integer,       parameter :: MST = MAX_ION_STAGES, MEL = MAX_METAL_ELEMENTS
  !--- round-off criterion on the reaction rates and on the algebra.
  real(kind=wp), parameter :: TOL_PAIR    = 1.0e-14_wp
  real(kind=wp), parameter :: TOL_BALANCE = 1.0e-14_wp
  !--- the one-iterate lag of the two sums inside the fixed point, with margin
  !--- (measured worst 1.6e-6, in the cold neutral state).
  real(kind=wp), parameter :: TOL_STATE   = 1.0e-5_wp
  !--- census: the nonzero CXI/CXR rows of data/atomic/element_*.txt
  !--- (transition it joins stage it and it+1; CXR is X^it+ + H0).
  !--- Transition 1 (both directions where carried): Huang et al. (2023)
  !--- Table 4, nitrogen KF96, sulfur by product term (Zhao et al. 2005);
  !--- recombination of q = 2-4: KF96 Tables 1-2, N2+ and O2+ from Barragan
  !--- et al. (2006).  No rate for Ne+, Ar+, Cl+, Ca+ + H (KF96 Table 1 has
  !--- none), Cl2+ (KF96 Table 2: radiative 1e-14, not carried) and Ca2+.
  integer, parameter :: NCX = 45
  character(len=2), parameter :: CX_EL(NCX) = [character(len=2) :: &
     'c ','c ','c ','c ','c ', 'n ','n ','n ','n ','n ', &
     'o ','o ','o ','o ','o ', 'ne','ne','ne', 's ','s ','s ','s ','s ', &
     'ar','ar','ar', 'mg','mg','mg','mg','mg', 'fe','fe','fe','fe','fe', &
     'si','si','si','si','si', 'cl','cl', 'ca','ca']
  integer, parameter :: CX_IT(NCX) = [1,1,2,3,4, 1,1,2,3,4, 1,1,2,3,4, &
     2,3,4, 1,1,2,3,4, 2,3,4, 1,1,2,3,4, 1,1,2,3,4, 1,1,2,3,4, 3,4, 3,4]
  !--- 1 ionization (CXI), 2 recombination (CXR)
  integer, parameter :: CX_DIR(NCX) = [1,2,2,2,2, 1,2,2,2,2, 1,2,2,2,2, &
     2,2,2, 1,2,2,2,2, 2,2,2, 1,2,2,2,2, 1,2,2,2,2, 1,2,2,2,2, 2,2, 2,2]
  integer, parameter :: NCXI_EXPECT = 7, NCXR_EXPECT = 38
  !--- Langevin rows: X^q + H for q = 5, 6 of all 11 elements
  integer, parameter :: NLANGEVIN = 22
  !--- KF96 recombination fits [a (1e-9 cm^3 s^-1), b, c, d, T_low, T_high]
  !--- for transitions 2-4, Table 1 except Cl3+, Cl4+, Fe4+ (Table 2 totals)
  !--- and N4+, Si2+ (updated database distributed with Cloudy, see check 5)
  integer, parameter :: NKF = 29
  character(len=2), parameter :: KF_EL(NKF) = [character(len=2) :: &
     'c ','c ','c ', 'n ','n ', 'o ','o ', 'ne','ne','ne', 'mg','mg','mg', &
     'si','si','si', 's ','s ','s ', 'ar','ar','ar', 'cl','cl', 'ca','ca', &
     'fe','fe','fe']
  integer, parameter :: KF_IT(NKF) = [2,3,4, 3,4, 3,4, 2,3,4, 2,3,4, 2,3,4, &
     2,3,4, 2,3,4, 3,4, 3,4, 2,3,4]
  real(kind=wp), parameter :: KF_P(6,NKF) = reshape([ &
     1.67e-4_wp, 2.79_wp,     304.72_wp, -4.07_wp,    5.0e3_wp, 5.0e4_wp, & ! C2+
     3.25_wp,    0.21_wp,     0.19_wp,   -3.29_wp,    1.0e3_wp, 1.0e5_wp, & ! C3+
     332.46_wp, -0.11_wp,    -0.995_wp,  -1.58e-3_wp, 1.0e1_wp, 1.0e5_wp, & ! C4+
     4.54_wp,    0.57_wp,    -0.65_wp,   -0.89_wp,    1.0e1_wp, 1.0e5_wp, & ! N3+
     2.95_wp,    0.55_wp,    -0.39_wp,   -1.07_wp,    1.0e3_wp, 1.0e6_wp, & ! N4+ (Cloudy update)
     3.98_wp,    0.26_wp,     0.56_wp,   -2.62_wp,    1.0e3_wp, 5.0e4_wp, & ! O3+
     0.252_wp,   0.63_wp,     2.08_wp,   -4.16_wp,    1.0e3_wp, 3.0e4_wp, & ! O4+
     1.0e-5_wp,  0.0_wp,      0.0_wp,     0.0_wp,     5.0e3_wp, 5.0e4_wp, & ! Ne2+
     14.73_wp,   4.52e-2_wp, -0.84_wp,   -0.31_wp,    5.0e3_wp, 5.0e4_wp, & ! Ne3+
     6.47_wp,    0.54_wp,     3.59_wp,   -5.22_wp,    1.0e3_wp, 3.0e4_wp, & ! Ne4+
     8.58e-5_wp, 2.49e-3_wp,  2.93e-2_wp,-4.33_wp,    1.0e3_wp, 3.0e4_wp, & ! Mg2+
     6.49_wp,    0.53_wp,     2.82_wp,   -7.63_wp,    1.0e3_wp, 3.0e4_wp, & ! Mg3+
     6.36_wp,    0.55_wp,     3.86_wp,   -5.19_wp,    1.0e3_wp, 3.0e4_wp, & ! Mg4+
     6.77_wp,    7.36e-2_wp, -0.43_wp,   -0.11_wp,    5.0e2_wp, 1.0e5_wp, & ! Si2+ (Cloudy update)
     0.490_wp,  -8.74e-2_wp, -0.36_wp,   -0.79_wp,    1.0e3_wp, 3.0e4_wp, & ! Si3+
     7.58_wp,    0.37_wp,     1.06_wp,   -4.09_wp,    1.0e3_wp, 5.0e4_wp, & ! Si4+
     1.0e-5_wp,  0.0_wp,      0.0_wp,     0.0_wp,     1.0e3_wp, 3.0e4_wp, & ! S2+
     2.29_wp,    4.02e-2_wp,  1.59_wp,   -6.06_wp,    1.0e3_wp, 3.0e4_wp, & ! S3+
     6.44_wp,    0.13_wp,     2.69_wp,   -5.69_wp,    1.0e3_wp, 3.0e4_wp, & ! S4+
     1.0e-5_wp,  0.0_wp,      0.0_wp,     0.0_wp,     1.0e3_wp, 3.0e4_wp, & ! Ar2+
     4.57_wp,    0.27_wp,    -0.18_wp,   -1.57_wp,    1.0e3_wp, 3.0e4_wp, & ! Ar3+
     6.37_wp,    0.85_wp,     10.21_wp,  -6.22_wp,    1.0e3_wp, 3.0e4_wp, & ! Ar4+
     1.88_wp,    0.32_wp,     1.77_wp,   -5.70_wp,    1.0e3_wp, 3.0e4_wp, & ! Cl3+ (Table 2)
     7.27_wp,    0.29_wp,     1.04_wp,  -10.14_wp,    1.0e3_wp, 3.0e4_wp, & ! Cl4+ (Table 2)
     3.17e-2_wp, 2.12_wp,     12.06_wp,  -0.40_wp,    1.0e3_wp, 3.0e4_wp, & ! Ca3+
     2.68_wp,    0.69_wp,    -0.68_wp,   -4.47_wp,    1.0e3_wp, 3.0e4_wp, & ! Ca4+
     1.26_wp,    7.72e-2_wp, -0.41_wp,   -7.31_wp,    1.0e3_wp, 1.0e5_wp, & ! Fe2+
     3.42_wp,    0.51_wp,    -2.06_wp,   -8.99_wp,    1.0e3_wp, 1.0e5_wp, & ! Fe3+
     14.60_wp,   3.57e-2_wp, -0.92_wp,   -0.37_wp,    1.0e3_wp, 3.0e4_wp  & ! Fe4+ (Table 2)
     ], [6, NKF])
  integer, parameter :: NTKF = 5
  !--- inside the intervals, and outside them where the rate is held
  real(kind=wp), parameter :: TKF(NTKF) = [2.0e3_wp, 5.0e3_wp, 1.0e4_wp, &
                                           2.0e4_wp, 6.0e4_wp]
  !--- log-log interpolation of the files' form-8 tables: measured worst
  !--- 4.6e-6 over the 145 points (2026-10-11, N4+ and Si2+ updated)
  real(kind=wp), parameter :: TOL_KF96 = 1.0e-4_wp
  !--- Langevin rate: standard atomic weights [u] (CIAAW), H 1.00794.
  !--- Measured 5.4e-5 .. 5.6e-5 on all 22 rows (2026-10-11), the same on
  !--- every element, so a difference in the constants rather than the masses
  real(kind=wp), parameter :: TOL_LANGEVIN = 1.0e-3_wp
  real(kind=wp), parameter :: AMU_G = 1.66053906660e-24_wp, A0_CM = 0.529177210903e-8_wp
  real(kind=wp), parameter :: E_ESU = 4.80320471e-10_wp, ALPHA_H_AU = 4.5_wp
  character(len=2), parameter :: AW_EL(11) = [character(len=2) :: &
     'c ','n ','o ','ne','mg','si','s ','cl','ar','ca','fe']
  real(kind=wp), parameter :: AW(11) = [12.011_wp, 14.007_wp, 15.999_wp, &
     20.180_wp, 24.305_wp, 28.085_wp, 32.06_wp, 35.45_wp, 39.948_wp, &
     40.078_wp, 55.845_wp]
  real(kind=wp), parameter :: AW_H = 1.00794_wp
  !--- cell states: T [K], n_H [cm^-3], Gamma_H [s^-1].  Two ionized states
  !--- (one dense), one at the ionization front, one mostly neutral, and one
  !--- cold and dense, so every reaction is exercised with a stage
  !--- distribution that actually populates its partners.
  integer,       parameter :: NST = 5
  real(kind=wp), parameter :: T_ST(NST) = &
     [8.0e3_wp, 8.0e3_wp, 1.0e4_wp, 6.0e3_wp, 3.0e2_wp]
  real(kind=wp), parameter :: NH_ST(NST) = &
     [1.0e2_wp, 1.0e2_wp, 1.0e3_wp, 1.0e2_wp, 5.0e3_wp]
  real(kind=wp), parameter :: GH_ST(NST) = &
     [1.0e-9_wp, 1.0e-11_wp, 1.0e-8_wp, 1.0e-13_wp, 1.0e-16_wp]

  real(kind=wp) :: cx_in(MST,MEL), cx_out(MST,MEL)
  real(kind=wp) :: frac(MST)
  real(kind=wp) :: T, nH, gH, ne, xHI, xHeI, xHeII, nHI, nHII
  real(kind=wp) :: nem, r_in, r_out, k_cxi, k_cxr, aH, cH
  real(kind=wp) :: rate_H, rate_M, rel, worst_in, worst_out
  real(kind=wp) :: worst_in_all, worst_out_all
  real(kind=wp) :: resid, resid_form, resid_bad, xHI_form, xHI_bad
  real(kind=wp) :: worst_bal, worst_bad, worst_state, tol_form
  character(len=8) :: worst_in_id, worst_out_id
  character(len=16) :: arg
  integer :: ierr, ie, it, is, k, ncxi, ncxr, j, idir, ncxr_required
  logical :: ok, caseA, langevin, carried, listed
  real(kind=wp) :: kexp, t4, worst_kf, mu, worst_lv

  call MPI_INIT(ierr)
  mpar%p_rank = 0
  call get_command_argument(1, arg)
  langevin = trim(arg) == 'langevin'
  par%langevin_charge_exchange = langevin

  !--- registry with every element active, so every carried reaction is present.
  par%atomic_dir     = 'data/atomic'
  par%use_metals     = .true.
  par%metal_ne       = .true.
  par%charge_exchange = 'two_way'
  !--- the balance residual below is formed with the metal terms alone, so
  !--- the He-H exchange (gated in check_helium_charge_exchange) is held off
  par%helium_charge_exchange = .false.
  !--- and so is the Penning term of the He I 2^3S coupling
  !--- (gated in tests/hei_metastable_coupling)
  par%hei_metastable_coupling = .false.
  par%recombination_continuum = 'on_the_spot'
  par%abund_He       = 0.1_wp
  par%ci_model       = 'voronov'
  par%recomb_model   = 'detailed_balance'
  !--- the n-level (T, n_e) cooling tables play no part here and are built by
  !--- cooling_setup, which this gate does not call.
  par%cooling_model  = 'low_density'
  par%abund_C  = 2.2e-4_wp;  par%abund_N  = 4.0e-5_wp
  par%abund_O  = 3.3e-4_wp;  par%abund_Ne = 5.0e-5_wp
  par%abund_S  = 9.0e-6_wp;  par%abund_Ar = 3.0e-6_wp
  par%abund_Mg = 3.0e-6_wp;  par%abund_Fe = 3.0e-6_wp
  par%abund_Si = 3.0e-6_wp;  par%abund_Cl = 3.0e-7_wp
  par%abund_Ca = 2.0e-8_wp
  call species_setup(1)

  ok    = .true.
  caseA = trim(par%recombination_continuum) /= 'on_the_spot'

  !--- 1: the reaction census the two sums run over, reaction by reaction.
  ncxr_required = NCXR_EXPECT
  if (langevin) ncxr_required = NCXR_EXPECT + NLANGEVIN
  write(*,'(/,a,l2,a)') '=== reaction census (11 elements, 77 stages, T = 1e4 K, '// &
     'langevin_charge_exchange =', langevin, ') ==='
  ncxi = 0;  ncxr = 0
  do ie = 1, n_elements
     if (elem_nstage(ie) /= 7) then
        write(*,'(a,a,a,i0)') '  FAIL: ', trim(elem_name(ie)), ' carries stages ', &
           elem_nstage(ie)
        ok = .false.
     end if
     do it = 1, elem_nstage(ie)-1
        call species_cx_coefficients(ie, it, 1.0e4_wp, k_cxi, k_cxr)
        if (k_cxi > 0.0_wp) ncxi = ncxi + 1
        if (k_cxr > 0.0_wp) ncxr = ncxr + 1
        do idir = 1, 2
           carried = merge(k_cxi, k_cxr, idir == 1) > 0.0_wp
           listed = .false.
           do j = 1, NCX
              if (trim(CX_EL(j)) == trim(elem_name(ie)) .and. CX_IT(j) == it &
                  .and. CX_DIR(j) == idir) listed = .true.
           end do
           if (langevin .and. idir == 2 .and. it >= 5) listed = .true.
           if (carried .neqv. listed) then
              write(*,'(a,a,a,i0,a,a,a,l2,a,l2)') '  FAIL: ', trim(elem_name(ie)), &
                 ' transition ', it, merge(' CXI', ' CXR', idir == 1), &
                 ': carried', carried, ', expected', listed
              ok = .false.
           end if
        end do
     end do
  end do
  write(*,'(a,i3,a,i3,a)') ' CXI (H+ -> H0): ', ncxi, ' reactions (expected ', &
     NCXI_EXPECT, ')'
  write(*,'(a,i3,a,i3,a)') ' CXR (H0 -> H+): ', ncxr, ' reactions (expected ', &
     ncxr_required, ')'
  if (ncxi /= NCXI_EXPECT .or. ncxr /= ncxr_required) then
     write(*,'(a)') '  FAIL: the registry no longer carries the expected ' // &
        'set of charge-exchange reactions'
     ok = .false.
  end if

  !--- 5: KF96 recombination rates of transitions 2-4
  write(*,'(/,a)') '=== KF96 Tables 1-2 recombination rates, q = 2-4 ==='
  worst_kf = 0.0_wp
  do j = 1, NKF
     ie = element_index(KF_EL(j))
     do k = 1, NTKF
        t4 = min(max(TKF(k), KF_P(5,j)), KF_P(6,j))/1.0e4_wp
        kexp = KF_P(1,j)*1.0e-9_wp*t4**KF_P(2,j)*(1.0_wp + KF_P(3,j)*exp(KF_P(4,j)*t4))
        call species_cx_coefficients(ie, KF_IT(j), TKF(k), k_cxi, k_cxr)
        rel = abs(k_cxr/kexp - 1.0_wp)
        if (rel > worst_kf) worst_kf = rel
        if (rel > TOL_KF96) then
           write(*,'(a,a,i2,a,f8.0,a,es12.5,a,es12.5)') '  FAIL: ', KF_EL(j), &
              KF_IT(j), ' T =', TKF(k), ' code ', k_cxr, ' KF96 ', kexp
           ok = .false.
        end if
     end do
     call species_cx_coefficients(ie, KF_IT(j), 1.0e4_wp, k_cxi, k_cxr)
     write(*,'(a4,a,i1,a,es12.4)') KF_EL(j), ' q=', KF_IT(j), '  k(1e4 K) = ', k_cxr
  end do
  write(*,'(a,es10.3,a,es9.2,a,i0,a)') ' worst |k/k_KF96 - 1| = ', worst_kf, &
     '  (tolerance ', TOL_KF96, ', ', NKF*NTKF, ' points)'

  !--- 6: Langevin rows (q = 5, 6), or their absence
  write(*,'(/,a)') '=== Langevin rows, q = 5, 6 (Ferland et al. 1997 eq. 2) ==='
  worst_lv = 0.0_wp
  do ie = 1, n_elements
     j = 0
     do k = 1, size(AW_EL)
        if (trim(AW_EL(k)) == trim(elem_name(ie))) j = k
     end do
     mu = AW(j)*AW_H/(AW(j) + AW_H)*AMU_G
     do it = 5, elem_nstage(ie)-1
        kexp = 0.0_wp
        if (langevin) kexp = 2.0_wp*pi*real(it, wp)*E_ESU*sqrt(ALPHA_H_AU*A0_CM**3/mu)
        do k = 1, NTKF
           call species_cx_coefficients(ie, it, TKF(k), k_cxi, k_cxr)
           if (langevin) then
              rel = abs(k_cxr/kexp - 1.0_wp)
           else
              rel = abs(k_cxr) + abs(k_cxi)
           end if
           if (rel > worst_lv) worst_lv = rel
           if (rel > TOL_LANGEVIN .or. (.not. langevin .and. rel /= 0.0_wp)) then
              write(*,'(a,a,i2,a,f8.0,a,es12.5,a,es12.5)') '  FAIL: ', &
                 trim(elem_name(ie)), it, ' T =', TKF(k), ' code ', k_cxr, &
                 ' expected ', kexp
              ok = .false.
           end if
        end do
        if (langevin) write(*,'(a4,a,i1,a,es12.5,a,es12.5)') elem_name(ie), &
           ' q=', it, '  code ', k_cxr, '  Langevin ', kexp
     end do
  end do
  if (langevin) then
     write(*,'(a,es10.3,a,es9.2,a)') ' worst |k/k_Langevin - 1| = ', worst_lv, &
        '  (tolerance ', TOL_LANGEVIN, ', atomic weights)'
  else
     write(*,'(a,es10.3)') ' switched off: largest q >= 5 rate ', worst_lv
  end if

  !--- 2 + 3 + 4: on each converged cell state.
  write(*,'(/,a)') '=== pairwise conservation and hydrogen''s closed form ==='
  write(*,'(a6,a10,2a10,5a12)') 'state', 'T [K]', 'x_HI', 'x_HII', &
     'worst CXI', 'worst CXR', 'solved x_HI', 'closed form', 'no R_cx_in'
  worst_bal = 0.0_wp;  worst_state = 0.0_wp;  worst_bad = huge(1.0_wp)
  worst_in_all = 0.0_wp;  worst_out_all = 0.0_wp
  worst_in_id = '-';  worst_out_id = '-'
  do is = 1, NST
     T  = T_ST(is);  nH = NH_ST(is);  gH = GH_ST(is)
     !--- converge the cell: solve_ion_cell is a damped fixed point that
     !--- exits at a relative 1e-10 on n_e, so it is repeated (100 passes)
     !--- until the state stops moving and the one-iterate lag of the
     !--- charge-exchange sums has worked itself out.
     xHI = 0.5_wp;  xHeI = 0.5_wp;  xHeII = 0.5_wp
     do k = 1, 100
        call solve_ion_cell(gH, 0.0_wp, 0.0_wp, nH, T, caseA, &
                            xHI, xHeI, xHeII, ne, 1)
     end do
     nHI = nH*xHI;  nHII = nH*(1.0_wp - xHI)

     !--- hydrogen's side, reaction by reaction, from the cached chain.
     call species_ne_prepare(1, T)
     call species_metal_sums(nH, ne, nHI, nHII, nem, r_in, r_out, cx_in, cx_out)

     !--- the cascade's side, from species_fractions and the registry
     !--- coefficients it evaluates for itself.
     worst_in = 0.0_wp;  worst_out = 0.0_wp
     do ie = 1, n_elements
        call species_fractions(ie, 1, T, ne, nHI, nHII, frac)
        do it = 1, elem_nstage(ie)-1
           call species_cx_coefficients(ie, it, T, k_cxi, k_cxr, ne)
           !--- X^it + H+ -> X^(it+1) + H0
           if (k_cxi > 0.0_wp) then
              rate_M = k_cxi*elem_abund(ie)*nH*frac(it)*nHII
              rate_H = cx_in(it,ie)*nHII
              rel = abs(rate_H - rate_M)/max(abs(rate_M), tinest)
              worst_in = max(worst_in, rel)
              if (rel > worst_in_all) then
                 worst_in_all = rel
                 write(worst_in_id,'(a,a,i1)') trim(elem_name(ie)), ' ', it
              end if
              if (rel > TOL_PAIR) then
                 write(*,'(a,a,a,i1,a,es10.3,2(a,es14.7))') '  FAIL: CXI ', &
                    trim(elem_name(ie)), ' transition ', it, ' rel = ', rel, &
                    ', H side ', rate_H, ', metal side ', rate_M
                 ok = .false.
              end if
           end if
           !--- X^(it+1) + H0 -> X^it + H+
           if (k_cxr > 0.0_wp) then
              rate_M = k_cxr*elem_abund(ie)*nH*frac(it+1)*nHI
              rate_H = cx_out(it,ie)*nHI
              rel = abs(rate_H - rate_M)/max(abs(rate_M), tinest)
              worst_out = max(worst_out, rel)
              if (rel > worst_out_all) then
                 worst_out_all = rel
                 write(worst_out_id,'(a,a,i1)') trim(elem_name(ie)), ' ', it
              end if
              if (rel > TOL_PAIR) then
                 write(*,'(a,a,a,i1,a,es10.3,2(a,es14.7))') '  FAIL: CXR ', &
                    trim(elem_name(ie)), ' transition ', it, ' rel = ', rel, &
                    ', H side ', rate_H, ', metal side ', rate_M
                 ok = .false.
              end if
           end if
        end do
     end do

     !--- hydrogen's balance: on the state solve_ion_cell reached, on the
     !--- closed form evaluated here, and on the same form with R_cx_in
     !--- dropped from the denominator.
     if (caseA) then
        aH = alphaA_HII(T)
     else
        aH = alphaB_HII(T)
     end if
     cH = ci_HI(T)
     resid = balance_residual(xHI)
     worst_state = max(worst_state, resid)
     xHI_form = (aH*ne + r_in)/(gH + (cH + aH)*ne + r_in + r_out)
     resid_form = balance_residual(xHI_form)
     worst_bal = max(worst_bal, resid_form)
     xHI_bad = (aH*ne + r_in)/(gH + (cH + aH)*ne + r_out)
     resid_bad = balance_residual(xHI_bad)
     worst_bad = min(worst_bad, resid_bad)

     !--- round-off floor of the metric: n_HII = n_H (1 - x_HI) carries a
     !--- relative error eps/(1 - x_HI).
     tol_form = max(TOL_BALANCE, &
                    4.0_wp*epsilon(1.0_wp)/max(1.0_wp - xHI_form, tinest))

     write(*,'(i6,f10.1,2es10.2,5es12.3)') is, T, xHI, 1.0_wp - xHI, &
        worst_in, worst_out, resid, resid_form, resid_bad
     if (resid > TOL_STATE) then
        write(*,'(a,i2,a,es10.3)') '  FAIL: state ', is, &
           ': the solved state does not satisfy hydrogen''s balance, ' // &
           'residual = ', resid
        ok = .false.
     end if
     if (resid_form > tol_form) then
        write(*,'(a,i2,a,es10.3,a,es10.3)') '  FAIL: state ', is, &
           ': the closed form does not satisfy hydrogen''s balance, ' // &
           'residual = ', resid_form, ' against ', tol_form
        ok = .false.
     end if
     if (resid_bad <= tol_form) then
        write(*,'(a,i2,a,es10.3)') '  FAIL: state ', is, &
           ': dropping R_cx_in from the denominator is not detected, ' // &
           'residual = ', resid_bad
        ok = .false.
     end if
  end do
  write(*,'(a,a,a,es10.3,a,a,a,es10.3)') &
     ' worst reaction over all states: CXI ', trim(worst_in_id), ' at ', &
     worst_in_all, ',  CXR ', trim(worst_out_id), ' at ', worst_out_all

  write(*,'(a,es10.3,a,es10.3)') ' worst balance residual: solved state ', &
     worst_state, ', closed form ', worst_bal

  !--- 4: the denominator assertion.  Dropping R_cx_in from the denominator
  !--- must break the balance on every state; if it does not, the gate is
  !--- blind to the error it exists to catch.
  write(*,'(/,a)') '=== denominator assertion ==='
  write(*,'(a,es10.3)') ' smallest residual with R_cx_in dropped from the ' // &
     'denominator: ', worst_bad
  if (worst_bad <= TOL_BALANCE) then
     write(*,'(a)') '  FAIL: the truncated denominator satisfies the ' // &
        'balance, so this gate cannot detect it'
     ok = .false.
  end if

  write(*,'(/,a)') 'GATE (two-way charge exchange, pairwise conservation): ' &
     // merge('PASS', 'FAIL', ok)
  call MPI_FINALIZE(ierr)
  if (.not. ok) error stop 1

contains

  integer function element_index(name) result(ie_found)
    character(len=*), intent(in) :: name
    integer :: je
    ie_found = 0
    do je = 1, n_elements
       if (trim(elem_name(je)) == trim(name)) ie_found = je
    end do
    if (ie_found == 0) error stop 'element not in the registry'
  end function element_index

  !--- relative residual of  n_HI (Gamma_H + C_H n_e + R_cx_out)
  !--- = (alpha_H n_e + R_cx_in) n_HII  at the neutral fraction x, on the
  !--- state (T, nH, ne, r_in, r_out) of the enclosing loop.
  real(kind=wp) function balance_residual(x) result(eps)
    real(kind=wp), intent(in) :: x
    real(kind=wp) :: loss, gain
    loss = nH*x*(gH + cH*ne + r_out)
    gain = (aH*ne + r_in)*nH*(1.0_wp - x)
    eps  = abs(loss - gain)/max(loss, tinest)
  end function balance_residual

end program check_charge_exchange
