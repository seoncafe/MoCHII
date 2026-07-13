module ion_balance_mod
!---------------------------------------------------------------------------
! MoCHII: H/He photoionization equilibrium at fixed T_e (Stage G1).
!
! Per leaf, given the rate integrals Gamma_i from the ionizing-band J tally
! and T_e = par%te_fixed, solve
!     n_HI  (Gamma_HI  + n_e C_HI )  = alpha_[AB](H II)   n_e n_HII
!     n_HeI (Gamma_HeI + n_e C_HeI)  = alpha_[AB](He II)  n_e n_HeII
!     n_HeII(Gamma_HeII+ n_e C_HeII) = alpha_[AB](He III) n_e n_HeIII
! with the closure  n_e = nH [ x_HII + y_He (x_HeII + 2 x_HeIII) ]
! by damped fixed-point iteration on n_e (docs/PLAN.md section 2.5).
! Case A/B via par%case_ab ('B' = on-the-spot, the G1 default).
!
! The converged fractions are written to the shared gas-state arrays with
! under-relaxation par%ion_relax (1 = none; < 1 damps oscillation at sharp
! I-fronts, the known failure mode flagged in docs/PLAN.md section 3).
! Every rank runs the identical solve (inputs are ALLREDUCEd tallies), so
! all ranks agree on max|delta x_HII|; only h_rank 0 writes the shared
! arrays, bracketed by node barriers.
!---------------------------------------------------------------------------
  use define
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, gas_nleaf
  use gas_rates_mod, only : gamma_HI, gamma_HeI, gamma_HeII
  use recomb_mod
  implicit none
  private

  public :: gas_equilibrium_update, solve_ion_cell

contains

  !=========================================================================
  ! Single-cell H/He equilibrium at temperature T: damped fixed point on
  ! n_e (docs/PLAN.md section 2.5).  x* are in/out (start = current state).
  ! Shared by gas_equilibrium_update (fixed Te) and thermal_mod (Te solve).
  !=========================================================================
  ! With par%metal_ne and the optional leaf index il, electrons from the
  ! metal cascade join the n_e closure (species_ne; beyond the I-front
  ! they are the only electrons, e.g. n_e ~ n_CII in the PDR zone).
  subroutine solve_ion_cell(gH, gHe1, gHe2, nH, T, caseA, &
                            xHI, xHeI, xHeII, ne, il)
    use species_mod, only : species_ne_prepare, species_ne_cached, n_elements
    real(kind=wp), intent(in)    :: gH, gHe1, gHe2, nH, T
    logical,       intent(in)    :: caseA
    real(kind=wp), intent(inout) :: xHI, xHeI, xHeII
    real(kind=wp), intent(out)   :: ne
    integer,       intent(in), optional :: il
    real(kind=wp) :: aH, aHe2, aHe3, cH, cHe1, cHe2, yHe
    real(kind=wp) :: ne_old, xHeIII, r1, r2, denH
    logical :: with_metal_ne
    integer :: it

    with_metal_ne = .false.
    if (present(il)) with_metal_ne = par%metal_ne .and. par%use_metals &
                                     .and. n_elements > 0
    !--- T is fixed inside the fixed point: evaluate the metal rate
    !--- coefficients once (species_ne_cached reduces each iteration
    !--- below to multiply-adds).
    if (with_metal_ne) call species_ne_prepare(il, T)

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
    do it = 1, 200
       ne_old = ne
       denH = gH + (cH + aH)*ne
       if (denH > 0.0_wp) then
          xHI = aH*ne / denH
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
       if (with_metal_ne) &
          ne = ne + species_ne_cached(nH, ne_old, nH*xHI, nH*(1.0_wp - xHI))
       ne = max(0.5_wp*(ne + ne_old), 1.0e-12_wp*nH)
       if (abs(ne - ne_old) <= 1.0e-10_wp*ne) exit
    end do
  end subroutine solve_ion_cell

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
       call solve_ion_cell(gamma_HI(il), gamma_HeI(il), gamma_HeII(il), &
                           nH, T, caseA, xHI, xHeI, xHeII, ne, il)

       !--- under-relaxation toward the solved state.
       xHI_new(il)   = (1.0_wp - w)*gas_xHI(il)   + w*xHI
       xHeI_new(il)  = (1.0_wp - w)*gas_xHeI(il)  + w*xHeI
       xHeII_new(il) = (1.0_wp - w)*gas_xHeII(il) + w*xHeII
       xHeIII = max(0.0_wp, 1.0_wp - xHeI_new(il) - xHeII_new(il))
       ne_new(il) = nH * ((1.0_wp - xHI_new(il)) &
                    + yHe*(xHeII_new(il) + 2.0_wp*xHeIII))
       !--- metal electrons in the written state (par%metal_ne).
       if (par%metal_ne .and. par%use_metals) then
          block
            use species_mod, only : species_ne, n_elements
            if (n_elements > 0) ne_new(il) = ne_new(il) &
               + species_ne(il, T, ne_new(il), nH*xHI_new(il), &
                            nH*(1.0_wp - xHI_new(il)))
          end block
       end if
       max_dx = max(max_dx, abs(xHI_new(il) - gas_xHI(il)))
       vol = (2.0_wp*leaf_half(il))**3
       sum_dxv = sum_dxv + vol*abs(xHI_new(il) - gas_xHI(il))
       sum_xv  = sum_xv  + vol*(1.0_wp - xHI_new(il))
    end do
    dx_vol = sum_dxv / max(sum_xv, tinest)

    !--- write back (h_rank 0 only; identical values on every rank).
    call MPI_BARRIER(mpar%hostcomm, ierr)
    if (mpar%h_rank == 0) then
       do il = 1, gas_nleaf
          gas_xHI(il)   = xHI_new(il)
          gas_xHeI(il)  = xHeI_new(il)
          gas_xHeII(il) = xHeII_new(il)
          gas_ne(il)    = ne_new(il)
       end do
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)
    deallocate(xHI_new, xHeI_new, xHeII_new, ne_new)
  end subroutine gas_equilibrium_update

end module ion_balance_mod
