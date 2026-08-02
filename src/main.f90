program main
!---------------------------------------------------------------------------
! MoCHII — MOnte Carlo for H II regions (driver).
!
! Each source packet is transported by the analytic edge walk (the packet's
! zero-variance direct contribution to J_nu in every leaf it crosses); when
! there is no ionizing-band scattering that walk is the complete tally.
! The rate integrals Gamma/H then follow from the J tally, feeding the
! ionization equilibrium and thermal balance.
!---------------------------------------------------------------------------
  use define
  use setup_mod
  use grid_mod_amr
  use octree_mod,      only : amr_grid
  use ion_band_mod,    only : ion_setup, ion_Ltot
  use gas_opacity_mod, only : gas_opacity_setup, gas_opacity_fill
  use jtally_mod,      only : jtally_ion_setup, &
                              slab_tally_setup, slab_tally_reduce, &
                              slab_write_Imu, slab_tally_on
  use ion_score_mod,   only : ion_score_setup
  use raytrace_amr_mod,only : slab_walk_report
  use ionizing_field_mod, only : ionizing_field_setup, ionizing_field_and_rates
  use gas_rates_mod,   only : gas_rates_write, &
                              gamma_HI, gamma_HeI, gamma_HeII, &
                              heat_HI, heat_HeI, heat_HeII, &
                              run_converged, run_iters, run_final_dx, &
                              run_final_dte, run_field_consistent
  use ion_balance_mod, only : gas_equilibrium_update, &
                              hydrogen_balance_residual_grid
  use thermal_mod,     only : gas_thermal_update
  use cooling_mod,     only : cooling_setup
  use nlevel_cooling_mod, only : nlevel_cooling_report
  use species_mod,     only : species_setup, &
                              n_elements, elem_nstage, elem_abund, elem_eth
  use memory_mod,      only : destroy_shared_mem_all
  use utility
  use mpi
  implicit none

  type(grid_type)     :: grid
  real(kind=wp)       :: dtime, max_dx, max_dte, dx_vol, dte_vol
  real(kind=wp)       :: thermal_t0, thermal_dt
  real(kind=wp)       :: eps_hmax, eps_hvol
  logical             :: converged, use_vol
  integer             :: ierr, iter, niter
  character(len=32)   :: pass_label
  !--- active metal photoionization thresholds gathered for the
  !--- threshold-aligned band (par%ion_align_edges).
  real(kind=wp), allocatable :: metal_eth(:)
  integer             :: nmetal, ie, it

  call MPI_INIT(ierr)
  call time_stamp(dtime)

  !--- setup
  call read_input()
  call setup_procedure()
  call grid_create_amr(grid)
  !--- species BEFORE ion_setup: with par%ion_align_edges the band setup
  !--- queries the active metal photoionization thresholds.  species_setup
  !--- does not depend on the band (that use is inside the loops over leaves).
  !--- species also BEFORE opacity: with par%ion_metal_abs the opacity fill
  !--- consumes the registry stage fractions.
  if (par%use_metals) call species_setup(amr_grid%nleaf)
  !--- gather active metal photoionization thresholds for the aligned band.
  nmetal = 0
  if (par%use_metals) then
     do ie = 1, n_elements
        if (elem_abund(ie) > 0.0_wp) nmetal = nmetal + max(0, elem_nstage(ie)-1)
     end do
  end if
  allocate(metal_eth(max(nmetal,1)));  metal_eth = 0.0_wp
  if (nmetal > 0) then
     nmetal = 0
     do ie = 1, n_elements
        if (elem_abund(ie) <= 0.0_wp) cycle
        do it = 1, elem_nstage(ie) - 1
           nmetal = nmetal + 1
           metal_eth(nmetal) = elem_eth(ie, it)
        end do
     end do
  end if
  call ion_setup(metal_eth, nmetal)
  deallocate(metal_eth)
  call gas_opacity_setup()
  call jtally_ion_setup(amr_grid%nleaf)
  call ion_score_setup(amr_grid%nleaf)
  if (trim(par%source_geometry) == 'slab') then
     block
       logical :: ltop, lbot
       ltop = trim(par%slab_faces) == 'top'    .or. trim(par%slab_faces) == 'both'
       lbot = trim(par%slab_faces) == 'bottom' .or. trim(par%slab_faces) == 'both'
       call slab_tally_setup(par%slab_nmu, ion_Ltot, [ltop, lbot])
     end block
  end if
  !--- Milne ground-continuum tables: built once, after ion_setup has fixed
  !--- par%eion_max (the tables are truncated at the band edge) and after the
  !--- atomic data are loaded.  Read-only during transport.
  if (par%diffuse_field .and. trim(par%diffuse_energy_model) == 'milne') then
     block
       use milne_recomb_spectrum_mod, only : milne_setup
       call milne_setup()
     end block
  end if
  !--- He I 2^1S two-photon shape: needed whenever the excited-cascade
  !--- channel emits, independently of the ground-continuum model.
  if (par%diffuse_field .and. par%hei_diffuse) then
     block
       use hei_twophoton_mod, only : hei_2ph_setup
       use hei_cascade_mod,   only : hei_cascade_setup
       call hei_2ph_setup()
       call hei_cascade_setup()
     end block
  end if
  if (par%solve_te) call cooling_setup()
  !--- report the (T, n_e) cooling tables built by species_setup/cooling_setup
  !--- (par%cooling_model = 'local_ne'); a no-op otherwise.
  if (trim(par%cooling_model) == 'local_ne') call nlevel_cooling_report()

  !--- iteration: transport -> rates -> equilibrium -> opacity feedback,
  !--- until max |delta x_HII| < par%gas_tol.  gas_niter = 0 = single
  !--- pass (no equilibrium solve, fixed state).  The stellar tally is
  !--- recomputed from zero each iteration (opacity changed).
  niter = max(par%gas_niter, 1)
  call ionizing_field_setup()
  converged = .false.
  max_dx = 0.0_wp;  max_dte = 0.0_wp;  dx_vol = 0.0_wp;  dte_vol = 0.0_wp
  do iter = 1, niter
     write(pass_label,'(a,i4)') 'iteration ', iter
     call ionizing_field_and_rates(pass_label)
     if (par%gas_niter < 1) exit          ! rates only, no solve
     use_vol = trim(par%conv_crit) == 'vol'
     if (par%solve_te) then
        thermal_t0 = MPI_WTIME()
        call gas_thermal_update(max_dx, max_dte, dx_vol, dte_vol)
        thermal_dt = MPI_WTIME() - thermal_t0
        if (use_vol) then
           converged = dx_vol < par%gas_tol .and. dte_vol < par%gas_tol_te
        else
           converged = max_dx < par%gas_tol .and. max_dte < par%gas_tol_te
        end if
        if (mpar%p_rank == 0) then
           write(6,'(a,i4,2(a,es12.4))') &
              '     iteration ', iter, ': max |delta x_HII| = ', max_dx, &
              ',  max |delta Te|/Te = ', max_dte
           write(6,'(a,2(a,es12.4))') &
              '                 ', ' vol |delta x_HII| = ', dx_vol, &
              ',  vol |delta Te|/Te = ', dte_vol
           write(6,'(a,f10.3,a,i0,a)') &
              '                  thermal solve = ', thermal_dt, &
              ' s (', mpar%h_size, ' node-local MPI ranks)'
        end if
     else
        call gas_equilibrium_update(max_dx, dx_vol)
        if (use_vol) then
           converged = dx_vol < par%gas_tol
        else
           converged = max_dx < par%gas_tol
        end if
        if (mpar%p_rank == 0) write(6,'(a,i4,2(a,es12.4))') &
           '     iteration ', iter, ': max |delta x_HII| = ', max_dx, &
           ',  vol = ', dx_vol
     end if
     !--- residual of hydrogen's ionization balance on the written state,
     !--- including both charge-exchange terms: a reported diagnostic, it
     !--- does not gate convergence.
     call hydrogen_balance_residual_grid(eps_hmax, eps_hvol)
     if (mpar%p_rank == 0) write(6,'(a,2(a,es12.4))') &
        '                 ', ' max H balance residual = ', eps_hmax, &
        ',  vol = ', eps_hvol
     call gas_opacity_fill()
     !--- solution-driven I-front re-refinement (one event).
     if (par%refine_front .and. iter == par%refine_iter) then
        block
          use amr_refine_mod, only : amr_refine_front
          call amr_refine_front()
        end block
        converged = .false.        ! keep iterating on the new grid
     end if
     if (converged) then
        if (mpar%p_rank == 0) write(6,'(a,i4,a)') &
           '     converged after ', iter, ' iterations'
        exit
     end if
  end do

  !--- record the convergence state (written to the rates header) and warn
  !--- when the iteration hit the cap without converging.
  run_converged = converged
  run_iters     = min(iter, niter)
  run_final_dx  = max_dx
  run_final_dte = max_dte
  !--- Only rank 0 evaluates the (namelist-set) require_convergence flag and
  !--- issues the collective abort; a single-rank MPI_ABORT tears down the
  !--- whole communicator.  This avoids letting any other rank act on its own
  !--- copy of the flag.
  if (par%gas_niter >= 1 .and. .not. converged .and. mpar%p_rank == 0) then
     write(6,'(a,i0,2(a,es10.3),a)') &
        ' WARNING: gas iteration did NOT converge in ', niter, &
        ' iterations (max|dx_HII| = ', max_dx, ', max|dTe|/Te = ', max_dte, &
        '); the written state is not at equilibrium.'
     if (par%require_convergence) then
        write(6,'(a)') &
           ' ERROR: par%require_convergence is set; stopping before output.'
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if
  end if

  !--- Final consistency pass.  The loop above ends on a gas solve followed
  !--- by gas_opacity_fill, so the state it leaves behind has never been
  !--- transported through: the rates, the gas heating and heat_dust still
  !--- belong to the PREVIOUS iterate.  At convergence that mismatch is
  !--- bounded by par%gas_tol, but at the iteration cap it is not bounded at
  !--- all.  One more transport with no gas solve puts the written radiation
  !--- field and the written gas state on the same iterate; it costs 1/niter
  !--- of the run.  Skipped when par%ion_peel is set, because the imaging
  !--- pass below rebuilds the same tally from the same state, and when
  !--- par%gas_niter < 1, where no solve ever ran and the two already agree.
  !--- The flag starts .false. and is set only where the pairing is actually
  !--- established, so a path added later that forgets to rebuild the field
  !--- reports itself in the header instead of claiming consistency.
  if (par%gas_niter < 1) then
     run_field_consistent = .true.      ! no solve ran; state never changed
  else if (.not. par%ion_peel) then
     call ionizing_field_and_rates('final consistency pass')
     run_field_consistent = .true.
  end if

  !--- peel-off imaging pass: one extra transport of the CONVERGED state
  !--- with observer peeling (direct at emission for stellar and diffuse
  !--- packets; dust-scattered at every interaction via the hook).  The
  !--- band tally is rebuilt by this pass, so the written rates reflect
  !--- the same field the images carry.
  if (par%ion_peel) then
     call ionizing_field_and_rates('imaging pass (peel-off)', with_peel=.true.)
     !--- this pass transported the written state, so it doubles as the
     !--- consistency pass skipped above.
     run_field_consistent = .true.
  end if

  !--- output
  if (par%ion_add_dust) then
     block
       use dust_temp_mod, only : dust_temp_setup, dust_temp_compute
       use gas_rates_mod, only : heat_dust
       call dust_temp_setup()
       call dust_temp_compute(heat_dust)
     end block
     !--- transport the reradiated thermal infrared and close the dust
     !--- self-heating loop, so Heat_dust and T_dust carry the infrared the
     !--- grains reabsorb from each other before anything is written.
     if (par%dust_emis_transport) then
        block
          use dust_emis_transport_mod, only : dust_emis_transport_setup, &
                                              transport_dust_emission
          use gas_rates_mod, only : heat_dust
          call dust_emis_transport_setup()
          call transport_dust_emission(heat_dust)
        end block
     end if
  end if
  if (par%use_sec_ion) then
     block
       use gas_rates_mod, only : sec_dgamma_HI, sec_dgamma_HeI, &
            sec_heat_HI, heat_HI, heat_hard_HI
       if (mpar%p_rank == 0) write(*,'(a,5es16.8)') &
          ' SECION_TMP sum(dgHI,dgHeI,secHtHI,HtHI,hardHI)=', &
          sum(sec_dgamma_HI), sum(sec_dgamma_HeI), sum(sec_heat_HI), &
          sum(heat_HI), sum(heat_hard_HI)
     end block
  end if
  call gas_rates_write()
  if (slab_tally_on) then
     call slab_walk_report()
     call slab_tally_reduce()
     block
       use utility, only : get_base_name
       if (mpar%p_rank == 0) &
          call slab_write_Imu(trim(get_base_name(par%out_file))//'_slab_Imu.txt')
     end block
     if (mpar%p_rank == 0) write(*,'(2a)') ' ION: slab emergent I(mu) written to: ', &
        trim(get_base_name(par%out_file))//'_slab_Imu.txt'
  end if
  if (par%ion_add_dust) then
     block
       use dust_temp_mod, only : dust_ir_write
       use gas_rates_mod, only : heat_dust
       call dust_ir_write(heat_dust)
     end block
  end if
  if (par%ion_add_dust .and. par%dust_sed) then
     block
       use sedust_mod,    only : sedust_setup, sedust_compute_write
       use gas_rates_mod, only : heat_dust
       call sedust_setup()
       call sedust_compute_write(heat_dust)
     end block
  end if
  if (par%ion_peel) then
     block
       use ion_peel_mod, only : ion_peel_write
       call ion_peel_write()
     end block
  end if
  !--- The H I / He II (SH95) and He I (Porter) recombination lines need only
  !--- n_e, T_e and the He ionization state, so they are defined without the
  !--- metal registry and with a fixed temperature; lines_write loads its own
  !--- tables and its metal loop is empty when no element is active.  Gating
  !--- this on use_metals + solve_te made par%emis_output silently produce
  !--- nothing for an H/He run.
  block
    use lines_mod, only : lines_write
    call lines_write()
  end block
  if (par%hei_metastable) then
     block
       use hei_metastable_mod, only : hei_metastable_run
       call hei_metastable_run()
     end block
  end if
  if (par%solve_te) then
     block
       use nebcont_mod, only : nebcont_setup, nebcont_write
       call nebcont_setup()
       call nebcont_write()
     end block
  end if

  call time_stamp(dtime)
  if (mpar%p_rank == 0) then
     par%exetime = dtime/60.0_wp
     write(6,'(a,f8.3,a)') 'Total Execution Time         : ', par%exetime, ' mins'
     write(6,'(2a)')       ' >>> STOP  @ ', get_date_time()
  endif

  call destroy_shared_mem_all()
  call MPI_FINALIZE(ierr)
  stop

end program main
