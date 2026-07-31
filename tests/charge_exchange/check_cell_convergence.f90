program check_cell_convergence
!---------------------------------------------------------------------------
! Gate G-M1f (docs/ION_BALANCE_CONVERGENCE_PLAN.md): with two-way charge
! exchange a partially neutral cell must reach its converged state inside one
! solve_ion_cell call, and the neutral fraction must stay physical.
!
! It is a real gate rather than a formality.  Substituting into hydrogen's
! balance with the charge-exchange terms lagged gives a map whose slope tends
! to 1 in partially neutral gas: the x_HII = 0.02 state below then left the
! solve 8% high, with a relative deviation of 1.8e-3 against the converged
! state.  Reinstating that substitution fails this gate on the last two states.
!---------------------------------------------------------------------------
  use define
  use species_mod,     only : species_setup
  use ion_balance_mod, only : solve_ion_cell
  use mpi
  implicit none
  !--- Gamma_HI chosen to span x_HII from a third-ionized cell down past the
  !--- 0.02 that the substitution map could not converge.
  integer, parameter :: NST = 5
  real(kind=wp), parameter :: GH(NST) = [3.0e-12_wp, 8.0e-13_wp, 1.5e-13_wp, &
                                         3.0e-14_wp, 1.0e-14_wp]
  real(kind=wp), parameter :: TOL = 1.0e-9_wp
  real(kind=wp) :: T, nH, x, h1, h2, ne, xr, r1, r2, ner, dev, worst
  integer :: is, k, ierr
  logical :: ok

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir      = 'data/atomic'
  par%use_metals      = .true.
  par%metal_ne        = .true.
  par%charge_exchange = 'two_way'
  par%case_ab         = 'B'
  par%He_abund        = 0.1_wp
  par%ci_model        = 'voronov'
  par%recomb_model    = 'badnell_milne'
  par%cooling_model   = 'low_density'
  par%abund_C  = 2.2e-4_wp;  par%abund_N  = 4.0e-5_wp
  par%abund_O  = 3.3e-4_wp;  par%abund_Ne = 5.0e-5_wp
  par%abund_S  = 9.0e-6_wp
  par%abund_Ar = 0.0_wp;  par%abund_Mg = 0.0_wp;  par%abund_Fe = 0.0_wp
  par%abund_Si = 0.0_wp;  par%abund_Cl = 0.0_wp;  par%abund_Ca = 0.0_wp
  call species_setup(1)
  nH = 100.0_wp;  T = 1.0e4_wp;  ok = .true.;  worst = 0.0_wp

  write(*,'(/,a)') '=== G-M1f: one solve_ion_cell call against the converged state ==='
  write(*,'(a10,a16,a16,a14)') 'x_HII', 'one call', 'converged', '|rel dev|'
  do is = 1, NST
     x = 0.5_wp;  h1 = 0.5_wp;  h2 = 0.5_wp
     call solve_ion_cell(GH(is), 0.0_wp, 0.0_wp, nH, T, .false., x, h1, h2, ne, 1)
     !--- the reference: repeated whole-cell solves, which converge whatever
     !--- the inner method is
     xr = 0.5_wp;  r1 = 0.5_wp;  r2 = 0.5_wp
     do k = 1, 400
        call solve_ion_cell(GH(is), 0.0_wp, 0.0_wp, nH, T, .false., xr, r1, r2, ner, 1)
     end do
     dev = abs(x - xr)/max(xr, tinest)
     worst = max(worst, dev)
     write(*,'(es10.2,es16.6,es16.6,es14.3)') 1.0_wp-xr, 1.0_wp-x, 1.0_wp-xr, dev
     if (dev > TOL) ok = .false.
     !--- the neutral fraction is bracketed in [0,1] by construction
     !--- (hydrogen_neutral_fraction): assert it rather than trust it.
     if (x < 0.0_wp .or. x > 1.0_wp) then
        write(*,'(a,es14.6)') '  FAIL: unphysical x_HI = ', x
        ok = .false.
     end if
  end do

  write(*,'(/,a,es10.3,a,es8.1)') ' worst relative deviation: ', worst, &
     '   tolerance ', TOL
  if (ok) then
     write(*,'(a)') ' GATE (G-M1f, cell converges inside one solve): PASS'
  else
     write(*,'(a)') ' GATE (G-M1f, cell converges inside one solve): FAIL'
     call MPI_FINALIZE(ierr);  stop 1
  end if
  call MPI_FINALIZE(ierr)
end program check_cell_convergence
