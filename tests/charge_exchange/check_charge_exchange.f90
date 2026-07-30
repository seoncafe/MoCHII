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
!   1. reaction census: the registry carries 29 reactions with all 11
!      elements active (7 CXI + 22 CXR), the number the two sums run over;
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
!      the wrong form.
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

  !--- MAX_ST, MAX_EL of species_mod (the shape of its reaction-resolved
  !--- output); lines_mod carries the same literal for frac.
  integer,       parameter :: MST = 6, MEL = 11
  !--- round-off criterion on the reaction rates and on the algebra.
  real(kind=wp), parameter :: TOL_PAIR    = 1.0e-14_wp
  real(kind=wp), parameter :: TOL_BALANCE = 1.0e-14_wp
  !--- the one-iterate lag of the two sums inside the fixed point, with margin
  !--- (measured worst 1.6e-6, in the cold neutral state).
  real(kind=wp), parameter :: TOL_STATE   = 1.0e-5_wp
  integer,       parameter :: NCXI_EXPECT = 7, NCXR_EXPECT = 22
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
  integer :: ierr, ie, it, is, k, ncxi, ncxr
  logical :: ok, caseA

  call MPI_INIT(ierr)
  mpar%p_rank = 0

  !--- registry with every element active, so all 29 reactions are present.
  par%atomic_dir     = 'data/atomic'
  par%use_metals     = .true.
  par%metal_ne       = .true.
  par%charge_exchange = 'two_way'
  par%case_ab        = 'B'
  par%He_abund       = 0.1_wp
  par%ci_model       = 'voronov'
  par%recomb_model   = 'badnell_milne'
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
  caseA = trim(par%case_ab) == 'A'

  !--- 1: the reaction census the two sums run over.
  ncxi = 0;  ncxr = 0
  do ie = 1, n_elements
     do it = 1, elem_nstage(ie)-1
        call species_cx_coefficients(ie, it, 1.0e4_wp, k_cxi, k_cxr)
        if (k_cxi > 0.0_wp) ncxi = ncxi + 1
        if (k_cxr > 0.0_wp) ncxr = ncxr + 1
     end do
  end do
  write(*,'(/,a)') '=== reaction census (11 elements, T = 1e4 K) ==='
  write(*,'(a,i3,a,i3,a)') ' CXI (H+ -> H0): ', ncxi, ' reactions (expected ', &
     NCXI_EXPECT, ')'
  write(*,'(a,i3,a,i3,a)') ' CXR (H0 -> H+): ', ncxr, ' reactions (expected ', &
     NCXR_EXPECT, ')'
  if (ncxi /= NCXI_EXPECT .or. ncxr /= NCXR_EXPECT) then
     write(*,'(a)') '  FAIL: the registry no longer carries the expected ' // &
        'set of charge-exchange reactions'
     ok = .false.
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
           call species_cx_coefficients(ie, it, T, k_cxi, k_cxr)
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
