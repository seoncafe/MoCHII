module ion_balance_mod
!---------------------------------------------------------------------------
! MoCHII: H/He photoionization equilibrium at fixed T_e.
!
! Per leaf, given the rate integrals Gamma_i from the ionizing-band J tally
! and T_e = par%te_fixed, solve
!     n_HI  (Gamma_HI  + n_e C_HI )  = alpha_[AB](H II)   n_e n_HII
!     n_HeI (Gamma_HeI + n_e C_HeI)  = alpha_[AB](He II)  n_e n_HeII
!     n_HeII(Gamma_HeII+ n_e C_HeII) = alpha_[AB](He III) n_e n_HeIII
! with the closure  n_e = nH [ x_HII + y_He (x_HeII + 2 x_HeIII) ]
! by damped fixed-point iteration on n_e.
! Case A/B via par%case_ab ('B' = on-the-spot, the default).
!
! The converged fractions are written to the shared gas-state arrays with
! under-relaxation par%ion_relax (1 = none; < 1 damps oscillation at sharp
! I-fronts, a known failure mode).
! Only the node-local rank 0 runs the solve (the inputs are ALLREDUCEd
! tallies, so every node's rank 0 is identical); it writes the shared
! arrays directly and broadcasts max|delta x_HII| to the other node ranks,
! which skip that work entirely.
!---------------------------------------------------------------------------
  use define
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, &
                            gas_Te, gas_nleaf
  use gas_rates_mod, only : gamma_HI, gamma_HeI, gamma_HeII, &
                            sec_dgamma_HI, sec_dgamma_HeI
  use recomb_mod
  implicit none
  private

  public :: gas_equilibrium_update, solve_ion_cell, charge_neutrality_ne
  public :: hydrogen_balance_residual, hydrogen_balance_residual_grid

contains

  !=========================================================================
  ! Single-cell H/He equilibrium at temperature T: damped fixed point on
  ! n_e.  x* are in/out (start = current state).
  ! Shared by gas_equilibrium_update (fixed Te) and thermal_mod (Te solve).
  !=========================================================================
  ! With par%metal_ne and the optional leaf index il, electrons from the
  ! metal cascade join the n_e closure (species_ne; beyond the I-front
  ! they are the only electrons, e.g. n_e ~ n_CII in the PDR zone).
  !
  ! With par%charge_exchange = 'two_way' the reactions
  !     X^i + H^+ -> X^(i+1) + H^0        (rate R_cx_in  per H^+)
  !     X^(i+1) + H^0 -> X^i + H^+        (rate R_cx_out per H^0)
  ! enter hydrogen's balance as well as the metal cascade's, so the charge
  ! they move is conserved in the reacting pair:
  !     n_HI (Gamma_H + C_H n_e + R_cx_out) = (alpha_H n_e + R_cx_in) n_HII.
  ! Eliminating n_HII = n_H (1 - x_HI) keeps the closed form; BOTH terms
  ! appear in the denominator, R_cx_in for the same reason alpha_H n_e does.
  subroutine solve_ion_cell(gH, gHe1, gHe2, nH, T, caseA, &
                            xHI, xHeI, xHeII, ne, il)
    use species_mod, only : species_ne_prepare, species_metal_sums, n_elements
    real(kind=wp), intent(in)    :: gH, gHe1, gHe2, nH, T
    logical,       intent(in)    :: caseA
    real(kind=wp), intent(inout) :: xHI, xHeI, xHeII
    real(kind=wp), intent(out)   :: ne
    integer,       intent(in), optional :: il
    real(kind=wp) :: aH, aHe2, aHe3, cH, cHe1, cHe2, yHe
    real(kind=wp) :: ne_old, xHeIII, r1, r2, denH
    real(kind=wp) :: nem, r_cx_in, r_cx_out
    logical :: with_metal_ne, with_metal_cx
    integer :: it

    !--- the two feedbacks are independent: counting the metal electrons is
    !--- par%metal_ne, conserving charge in the charge-exchange pair is
    !--- par%charge_exchange.  Either one needs the (leaf, T) rate cache.
    with_metal_ne = .false.
    with_metal_cx = .false.
    if (present(il) .and. par%use_metals .and. n_elements > 0) then
       with_metal_ne = par%metal_ne
       with_metal_cx = trim(par%charge_exchange) == 'two_way'
    end if
    !--- T is fixed inside the fixed point: evaluate the metal rate
    !--- coefficients once (species_metal_sums reduces each iteration
    !--- below to multiply-adds).
    if (with_metal_ne .or. with_metal_cx) call species_ne_prepare(il, T)

    yHe = par%He_abund
    if (caseA) then
       aH = alphaA_HII(T);  aHe2 = alphaA_HeII(T);  aHe3 = alphaA_HeIII(T)
    else
       aH = alphaB_HII(T);  aHe2 = alphaB_HeII(T);  aHe3 = alphaB_HeIII(T)
    end if
    cH = ci_HI(T);  cHe1 = ci_HeI(T);  cHe2 = ci_HeII(T)

    xHeIII = max(0.0_wp, 1.0_wp - xHeI - xHeII)
    ne = nH * ((1.0_wp - xHI) + yHe*(xHeII + 2.0_wp*xHeIII))
    ne = max(ne, 1.0e-12_wp*nH)
    nem = 0.0_wp;  r_cx_in = 0.0_wp;  r_cx_out = 0.0_wp
    do it = 1, 200
       ne_old = ne
       !--- the stage densities the two charge-exchange sums need depend on
       !--- n_HI and n_HII through the charge-exchange terms of the cascade
       !--- itself, so they are taken from the previous iterate of this fixed
       !--- point — the same lag the electron sum carries.
       !--- Domain of validity of that lag: where the charge-exchange terms
       !--- dominate hydrogen's balance the iteration still converges, but
       !--- slowly.  Measured at T = 1e4 K, n_H = 100 with the HII20
       !--- abundances, the exit test below is met after 124 passes at
       !--- x_HII = 0.18, 239 at 0.10 and 1285 at 0.021, against the 200-pass
       !--- cap — so a partially neutral cell (x_HII below ~0.1) can leave
       !--- this loop still moving, with x_HII high by ~8% and a few 1e-3
       !--- residual in the balance.  Repeated whole-cell solves converge it
       !--- (the outer gas iteration does exactly that), and
       !--- hydrogen_balance_residual reports what is left.
       if (with_metal_ne .or. with_metal_cx) then
          call species_metal_sums(nH, ne_old, nH*xHI, nH*(1.0_wp - xHI), &
                                  nem, r_cx_in, r_cx_out)
          if (.not. with_metal_cx) then
             r_cx_in  = 0.0_wp
             r_cx_out = 0.0_wp
          end if
       end if
       denH = gH + (cH + aH)*ne + r_cx_in + r_cx_out
       if (denH > 0.0_wp) then
          xHI = (aH*ne + r_cx_in) / denH
       else
          xHI = 1.0_wp
       end if
       !--- r1 = He+/He0, r2 = He++/He+.  Cap at 1e150 (physical values stay
       !--- below ~1e15 given the ne floor) so r1*r2 cannot overflow to Inf
       !--- and produce a 0*Inf = NaN in xHeII below.
       r1 = min((gHe1 + cHe1*ne) / (aHe2*ne), 1.0e150_wp)
       r2 = min((gHe2 + cHe2*ne) / (aHe3*ne), 1.0e150_wp)
       xHeI   = 1.0_wp / (1.0_wp + r1 + r1*r2)
       xHeII  = xHeI * r1
       xHeIII = xHeII * r2
       ne = nH * ((1.0_wp - xHI) + yHe*(xHeII + 2.0_wp*xHeIII))
       if (with_metal_ne) ne = ne + nem
       ne = max(0.5_wp*(ne + ne_old), 1.0e-12_wp*nH)
       if (abs(ne - ne_old) <= 1.0e-10_wp*ne) exit
    end do
  end subroutine solve_ion_cell

  !=========================================================================
  ! Electron density of a state whose H/He fractions are already fixed, from
  ! charge neutrality including the metal cascade (par%metal_ne):
  !     n_e = n_e(H,He) + sum_el n_H A_el sum_i (i-1) f_i(n_e)
  ! The metal stage fractions f_i depend on n_e themselves (recombination
  ! goes as n_e alpha), so this is a scalar fixed point, solved with the
  ! same damped iteration and tolerance solve_ion_cell uses for its own
  ! n_e closure.  A single evaluation of the metal term at the H/He-only
  ! electron density is NOT the solution: beyond the ionization front the
  ! H/He term vanishes, so that argument is the n_e -> 0 limit in which
  ! nothing can recombine and every low-IP metal comes out fully ionized.
  !
  ! Its purpose is the WRITTEN state: the under-relaxed H/He fractions are
  ! not the ones solve_ion_cell converged with, and the stage fractions
  ! written alongside (species_write) are evaluated at the written n_e and
  ! T_e — so n_e must satisfy charge neutrality at exactly those fractions
  ! or the output describes two different gases.
  !
  ! ne_seed is the electron density the cell solve converged to; starting
  ! there keeps the iteration on the branch the physics solve selected.
  ! T is the temperature the state is WRITTEN at (par%te_fixed, or the
  ! newly solved T_e), which is not in general the last trial temperature
  ! of the thermal bisection, so the (leaf, T) rate cache is refilled here
  ! rather than inherited.
  !=========================================================================
  real(kind=wp) function charge_neutrality_ne(il, T, nH, ne_HHe, xHI, ne_seed) &
                                             result(ne)
    use species_mod, only : species_ne_prepare, species_ne_cached
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne_HHe, xHI, ne_seed
    real(kind=wp) :: ne_old, nHI, nHII
    integer :: it

    call species_ne_prepare(il, T)
    nHI  = nH*xHI
    nHII = nH*(1.0_wp - xHI)
    ne   = max(ne_seed, 1.0e-12_wp*nH)
    do it = 1, 200
       ne_old = ne
       ne = ne_HHe + species_ne_cached(nH, ne_old, nHI, nHII)
       ne = max(0.5_wp*(ne + ne_old), 1.0e-12_wp*nH)
       if (abs(ne - ne_old) <= 1.0e-10_wp*ne) exit
    end do
  end function charge_neutrality_ne

  !=========================================================================
  ! Relative residual of hydrogen's ionization balance on the WRITTEN state
  ! of leaf il, including both charge-exchange terms:
  !
  !   eps = | n_HI (Gamma_H + C_H n_e + R_cx_out) - (alpha_H n_e + R_cx_in) n_HII |
  !         / [ n_HI (Gamma_H + C_H n_e + R_cx_out) ]
  !
  ! Inside solve_ion_cell the two sides of every charge-exchange pair are
  ! formed from the same cached coefficients and the same stage densities, so
  ! the pair balances by construction there.  What this measures is the state
  ! actually written out, which differs from the solved one because the H/He
  ! fractions are under-relaxed toward it (par%ion_relax) while the metal
  ! stage fractions written alongside are evaluated at the written n_e and T.
  !
  ! With par%charge_exchange = 'metal_only' the two R_cx terms are dropped
  ! here as well, so the residual always measures the equation the run
  ! actually solves.
  !
  ! Domain of validity: the state carries x_HI, so n_HII = n_H (1 - x_HI) has
  ! the relative precision eps/(1 - x_HI) and the residual cannot resolve
  ! anything below that.  In a fully neutral leaf 1 - x_HI rounds to zero and
  ! eps saturates at 1 whatever the physics does; such leaves are excluded
  ! from the reported maximum and mean by XHII_FLOOR below rather than being
  ! allowed to set them.
  !=========================================================================
  real(kind=wp) function hydrogen_balance_residual(il, T, nH, xHI, ne) result(eps)
    use species_mod, only : species_ne_prepare, species_metal_sums, n_elements
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, xHI, ne
    real(kind=wp) :: aH, cH, nHI, nHII, nem, r_cx_in, r_cx_out, loss, gain

    r_cx_in  = 0.0_wp
    r_cx_out = 0.0_wp
    if (par%use_metals .and. n_elements > 0 .and. &
        trim(par%charge_exchange) == 'two_way') then
       call species_ne_prepare(il, T)
       call species_metal_sums(nH, ne, nH*xHI, nH*(1.0_wp - xHI), &
                               nem, r_cx_in, r_cx_out)
    end if
    if (trim(par%case_ab) == 'A') then
       aH = alphaA_HII(T)
    else
       aH = alphaB_HII(T)
    end if
    cH   = ci_HI(T)
    nHI  = nH*xHI
    nHII = nH*(1.0_wp - xHI)
    loss = nHI*(gamma_HI(il) + sec_dgamma_HI(il) + cH*ne + r_cx_out)
    gain = (aH*ne + r_cx_in)*nHII
    eps  = abs(loss - gain) / max(loss, tinest)
  end function hydrogen_balance_residual

  !=========================================================================
  ! The same residual over the grid: its maximum over leaves and its
  ! volume-weighted mean, over the leaves where it resolves anything
  ! (x_HII above XHII_FLOOR, so n_HII = n_H (1 - x_HI) keeps at least six
  ! significant digits).  Diagnostic only — it is printed, not acted on.
  !=========================================================================
  subroutine hydrogen_balance_residual_grid(eps_max, eps_vol)
    use mpi
    use octree_mod, only : leaf_half
    implicit none
    real(kind=wp), intent(out) :: eps_max, eps_vol
    real(kind=wp), parameter :: XHII_FLOOR = 1.0e-10_wp
    real(kind=wp) :: eps, vol, sum_ev, sum_v
    integer :: il, ierr

    eps_max = 0.0_wp
    eps_vol = 0.0_wp
    !--- the gas state is node-shared and every node's rank 0 holds the same
    !--- values, so one rank per node evaluates and broadcasts.
    if (mpar%h_rank == 0) then
       sum_ev = 0.0_wp
       sum_v  = 0.0_wp
       do il = 1, gas_nleaf
          if (gas_nH(il) <= 0.0_wp) cycle
          if (1.0_wp - gas_xHI(il) < XHII_FLOOR) cycle
          eps = hydrogen_balance_residual(il, gas_Te(il), gas_nH(il), &
                                          gas_xHI(il), gas_ne(il))
          eps_max = max(eps_max, eps)
          vol = (2.0_wp*leaf_half(il))**3
          sum_ev = sum_ev + vol*eps
          sum_v  = sum_v  + vol
       end do
       eps_vol = sum_ev / max(sum_v, tinest)
    end if
    call MPI_BCAST(eps_max, 1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
    call MPI_BCAST(eps_vol, 1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
  end subroutine hydrogen_balance_residual_grid

  !=========================================================================
  ! Returns both convergence measures (par%conv_crit picks which gates the
  ! iteration): max_dx = max leaf |delta x_HII| (stalls at the front-cell
  ! Monte Carlo noise floor); dx_vol = sum V |delta x_HII| / sum V x_HII,
  ! the volume-weighted L1 change of the ionized volume — front jitter
  ! occupies little volume, so this measure converges when the global
  ! solution has.
  subroutine gas_equilibrium_update(max_dx, dx_vol)
    use mpi
    use octree_mod, only : amr_grid, leaf_half
    implicit none
    real(kind=wp), intent(out) :: max_dx, dx_vol

    real(kind=wp), allocatable :: xHI_new(:), xHeI_new(:), xHeII_new(:), ne_new(:)
    real(kind=wp) :: T, yHe, w, vol, sum_dxv, sum_xv
    real(kind=wp) :: nH, ne, xHI, xHeI, xHeII, xHeIII
    logical :: caseA
    integer :: il, ierr

    T   = par%te_fixed
    yHe = par%He_abund
    w   = par%ion_relax
    caseA = trim(par%case_ab) == 'A'

    !--- only the node-local rank 0 solves.  jt_ion is ALLREDUCEd over
    !--- MPI_COMM_WORLD, so every node's rank 0 sees the identical tally and
    !--- produces identical values; the other node ranks skip the (expensive)
    !--- solve for each leaf and receive the state through shared memory + the
    !--- convergence metrics by broadcast.
    if (mpar%h_rank == 0) then
    allocate(xHI_new(gas_nleaf), xHeI_new(gas_nleaf), xHeII_new(gas_nleaf), &
             ne_new(gas_nleaf))

    max_dx  = 0.0_wp
    sum_dxv = 0.0_wp
    sum_xv  = 0.0_wp
    do il = 1, gas_nleaf
       nH = gas_nH(il)
       if (nH <= 0.0_wp) then
          xHI_new(il) = gas_xHI(il);  xHeI_new(il) = gas_xHeI(il)
          xHeII_new(il) = gas_xHeII(il);  ne_new(il) = 0.0_wp
          cycle
       end if

       xHI = gas_xHI(il);  xHeI = gas_xHeI(il);  xHeII = gas_xHeII(il)
       call solve_ion_cell(gamma_HI(il) + sec_dgamma_HI(il), &
                           gamma_HeI(il) + sec_dgamma_HeI(il), gamma_HeII(il), &
                           nH, T, caseA, xHI, xHeI, xHeII, ne, il)

       !--- under-relaxation toward the solved state.
       xHI_new(il)   = (1.0_wp - w)*gas_xHI(il)   + w*xHI
       xHeI_new(il)  = (1.0_wp - w)*gas_xHeI(il)  + w*xHeI
       xHeII_new(il) = (1.0_wp - w)*gas_xHeII(il) + w*xHeII
       xHeIII = max(0.0_wp, 1.0_wp - xHeI_new(il) - xHeII_new(il))
       ne_new(il) = nH * ((1.0_wp - xHI_new(il)) &
                    + yHe*(xHeII_new(il) + 2.0_wp*xHeIII))
       !--- metal electrons in the written state (par%metal_ne), solved
       !--- self-consistently: the stage fractions depend on n_e, so the
       !--- closure has to be iterated at the under-relaxed H/He state, not
       !--- evaluated once at the H/He-only electron density.
       if (par%metal_ne .and. par%use_metals) then
          block
            use species_mod, only : n_elements
            if (n_elements > 0) ne_new(il) = &
               charge_neutrality_ne(il, T, nH, ne_new(il), xHI_new(il), ne)
          end block
       end if
       max_dx = max(max_dx, abs(xHI_new(il) - gas_xHI(il)))
       vol = (2.0_wp*leaf_half(il))**3
       sum_dxv = sum_dxv + vol*abs(xHI_new(il) - gas_xHI(il))
       sum_xv  = sum_xv  + vol*(1.0_wp - xHI_new(il))
    end do
    dx_vol = sum_dxv / max(sum_xv, tinest)

    !--- write the solved state straight into the node-shared arrays.
    do il = 1, gas_nleaf
       gas_xHI(il)   = xHI_new(il)
       gas_xHeI(il)  = xHeI_new(il)
       gas_xHeII(il) = xHeII_new(il)
       gas_ne(il)    = ne_new(il)
    end do
    deallocate(xHI_new, xHeI_new, xHeII_new, ne_new)
    end if   ! h_rank == 0

    !--- broadcast the metrics to the node's other ranks and barrier so the
    !--- shared-memory writes are visible before this returns.
    call MPI_BCAST(max_dx, 1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
    call MPI_BCAST(dx_vol, 1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine gas_equilibrium_update

end module ion_balance_mod
