module thermal_mod
!---------------------------------------------------------------------------
! MoCHII: thermal balance — solve T_e with the ionization state (Stage G2).
!
! Per leaf, find T_e such that photoheating balances cooling:
!     H(T) = n_HI H_HI + n_HeI H_HeI + n_HeII H_HeII     (rate integrals,
!            per particle, from gas_rates_mod — T-independent)
!     C(T) = cooling_total(T, ...)                        (cooling_mod)
! where the ionization fractions at each trial T come from solve_ion_cell
! (photo + collisional vs recombination at that T).  Since C rises and the
! ionization response is mild, the net function H - C is monotone through
! the root: bisection on log T in [par%te_min, par%te_max] (~40 iterations
! to 1e-4 relative).
!
! gas_thermal_update replaces gas_equilibrium_update when par%solve_te:
! it writes x, n_e AND T_e (under-relaxation on x as in G1) and returns
! max |delta x_HII| and max |delta T_e|/T_e.
!---------------------------------------------------------------------------
  use define
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, &
                            gas_Te, gas_nleaf
  use gas_rates_mod, only : gamma_HI, gamma_HeI, gamma_HeII, &
                            heat_HI, heat_HeI, heat_HeII
  use ion_balance_mod, only : solve_ion_cell
  use cooling_mod,     only : cooling_total
  implicit none
  private

  public :: gas_thermal_update

contains

  !=========================================================================
  ! Net heating - cooling [erg cm^-3 s^-1] at trial temperature T for one
  ! leaf; the ionization state at T is solved self-consistently.
  !=========================================================================
  subroutine net_rate(il, T, caseA, xHI, xHeI, xHeII, ne, net)
    integer,       intent(in)    :: il
    real(kind=wp), intent(in)    :: T
    logical,       intent(in)    :: caseA
    real(kind=wp), intent(inout) :: xHI, xHeI, xHeII
    real(kind=wp), intent(out)   :: ne, net
    real(kind=wp) :: nH, heat

    nH = gas_nH(il)
    call solve_ion_cell(gamma_HI(il), gamma_HeI(il), gamma_HeII(il), &
                        nH, T, caseA, xHI, xHeI, xHeII, ne, il)
    heat = nH*( xHI*heat_HI(il) &
           + par%He_abund*(xHeI*heat_HeI(il) + xHeII*heat_HeII(il)) )
    net = heat - cooling_total(T, nH, ne, xHI, xHeI, xHeII)
    !--- trace-metal line cooling (registry, G2b/c)
    if (par%use_metals) then
       block
         use species_mod, only : metal_cooling, metal_heating
         net = net - metal_cooling(il, T, nH, ne, nH*xHI, nH*(1.0_wp-xHI))
         !--- PDR: metal photoheating (par%metal_heat)
         if (par%metal_heat) &
            net = net + metal_heating(il, T, nH, ne, nH*xHI, nH*(1.0_wp-xHI))
       end block
    end if
    !--- PDR: grain photoelectric heating + grain recombination cooling
    !--- (Bakes & Tielens 1994 fits) driven by the local FUV field.
    if (par%grain_pe) then
       block
         use gas_rates_mod, only : g0_fuv
         real(kind=wp) :: g0, xpe, eps, beta, dscale
         g0 = g0_fuv(il)
         if (g0 > 0.0_wp) then
            dscale = par%pe_scale*par%DGR/1.0e-2_wp
            xpe = g0*sqrt(T)/max(ne, 1.0e-12_wp*nH)
            eps = 4.87e-2_wp/(1.0_wp + 4.0e-3_wp*xpe**0.73_wp) &
                + 3.65e-2_wp*(T/1.0e4_wp)**0.7_wp/(1.0_wp + 2.0e-4_wp*xpe)
            net = net + 1.3e-24_wp*eps*g0*nH*dscale
            beta = 0.74_wp/T**0.068_wp
            net = net - 4.65e-30_wp*T**0.94_wp*xpe**beta*ne*nH*dscale
         end if
       end block
    end if
  end subroutine net_rate

  !=========================================================================
  ! Returns cell-max AND volume-integrated convergence measures
  ! (par%conv_crit picks which pair gates the iteration):
  !   dx_vol  = sum V |delta x_HII| / sum V x_HII
  !   dte_vol = sum V n_e |delta T_e| / sum V n_e T_e
  ! The n_e weight removes the cells with vanishing heating (deep PDR /
  ! vacuum), whose te_min-pinned oscillation dominates max_dte.
  subroutine gas_thermal_update(max_dx, max_dte, dx_vol, dte_vol)
    use mpi
    use octree_mod, only : amr_grid
    implicit none
    real(kind=wp), intent(out) :: max_dx, max_dte, dx_vol, dte_vol

    real(kind=wp), allocatable :: xHI_new(:), xHeI_new(:), xHeII_new(:), &
                                  ne_new(:), te_new(:)
    real(kind=wp) :: nH, ne, xHI, xHeI, xHeII, xHeIII, w, vol
    real(kind=wp) :: sum_dxv, sum_xv, sum_dtev, sum_tev
    real(kind=wp) :: tlo, thi, tmid, net_lo, net_hi, net_mid, te
    logical :: caseA
    integer :: il, it, ierr

    caseA = trim(par%case_ab) == 'A'
    w     = par%ion_relax
    allocate(xHI_new(gas_nleaf), xHeI_new(gas_nleaf), xHeII_new(gas_nleaf), &
             ne_new(gas_nleaf), te_new(gas_nleaf))

    max_dx   = 0.0_wp
    max_dte  = 0.0_wp
    sum_dxv  = 0.0_wp;  sum_xv  = 0.0_wp
    sum_dtev = 0.0_wp;  sum_tev = 0.0_wp
    do il = 1, gas_nleaf
       nH = gas_nH(il)
       if (nH <= 0.0_wp) then
          xHI_new(il) = gas_xHI(il);  xHeI_new(il) = gas_xHeI(il)
          xHeII_new(il) = gas_xHeII(il);  ne_new(il) = 0.0_wp
          te_new(il) = gas_Te(il)
          cycle
       end if

       !--- bisection on log T; ionization re-solved at each trial T.
       xHI = gas_xHI(il);  xHeI = gas_xHeI(il);  xHeII = gas_xHeII(il)
       tlo = par%te_min;  thi = par%te_max
       call net_rate(il, tlo, caseA, xHI, xHeI, xHeII, ne, net_lo)
       call net_rate(il, thi, caseA, xHI, xHeI, xHeII, ne, net_hi)
       if (net_lo <= 0.0_wp) then
          te = tlo               ! cooling wins even at te_min (no heating)
       else if (net_hi >= 0.0_wp) then
          te = thi               ! heating wins even at te_max
       else
          do it = 1, 60
             tmid = sqrt(tlo*thi)
             call net_rate(il, tmid, caseA, xHI, xHeI, xHeII, ne, net_mid)
             if (net_mid > 0.0_wp) then
                tlo = tmid
             else
                thi = tmid
             end if
             if (thi/tlo - 1.0_wp < 1.0e-5_wp) exit
          end do
          te = sqrt(tlo*thi)
       end if
       !--- final consistent state at te.
       call net_rate(il, te, caseA, xHI, xHeI, xHeII, ne, net_mid)

       !--- under-relaxation on x (as in G1); te taken directly.
       xHI_new(il)   = (1.0_wp - w)*gas_xHI(il)   + w*xHI
       xHeI_new(il)  = (1.0_wp - w)*gas_xHeI(il)  + w*xHeI
       xHeII_new(il) = (1.0_wp - w)*gas_xHeII(il) + w*xHeII
       xHeIII = max(0.0_wp, 1.0_wp - xHeI_new(il) - xHeII_new(il))
       ne_new(il) = nH * ((1.0_wp - xHI_new(il)) &
                    + par%He_abund*(xHeII_new(il) + 2.0_wp*xHeIII))
       te_new(il) = te
       max_dx  = max(max_dx,  abs(xHI_new(il) - gas_xHI(il)))
       max_dte = max(max_dte, abs(te_new(il) - gas_Te(il))/gas_Te(il))
       vol = (2.0_wp*amr_grid%ch(amr_grid%icell_of_leaf(il)))**3
       sum_dxv  = sum_dxv  + vol*abs(xHI_new(il) - gas_xHI(il))
       sum_xv   = sum_xv   + vol*(1.0_wp - xHI_new(il))
       sum_dtev = sum_dtev + vol*ne_new(il)*abs(te_new(il) - gas_Te(il))
       sum_tev  = sum_tev  + vol*ne_new(il)*te_new(il)
    end do
    dx_vol  = sum_dxv  / max(sum_xv,  tinest)
    dte_vol = sum_dtev / max(sum_tev, tinest)

    call MPI_BARRIER(mpar%hostcomm, ierr)
    if (mpar%h_rank == 0) then
       do il = 1, gas_nleaf
          gas_xHI(il)   = xHI_new(il)
          gas_xHeI(il)  = xHeI_new(il)
          gas_xHeII(il) = xHeII_new(il)
          gas_ne(il)    = ne_new(il)
          gas_Te(il)    = te_new(il)
       end do
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)
    deallocate(xHI_new, xHeI_new, xHeII_new, ne_new, te_new)
  end subroutine gas_thermal_update

end module thermal_mod
