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
  use ion_band_mod,    only : ion_setup, ion_Ltot, ion_eedge, ion_lum
  use gas_opacity_mod, only : gas_opacity_setup, gas_opacity_fill, &
                              dust_reference_extinction
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
                              run_final_dte, run_state_solved_in_field
  use ion_balance_mod, only : gas_equilibrium_update, &
                              hydrogen_balance_residual_grid
  use thermal_mod,     only : gas_thermal_update
  use fixed_point_distance_mod, only : iteration_history, history_push, &
                                       history_reset, fixed_point_distance, &
                                       fixed_point_converged, net_drift, NWIN, &
                                       hist_dx, hist_dhe1, hist_dhe3
  use cooling_mod,     only : cooling_setup
  use heii_cascade_mod, only : heii_cascade_transported
  use nlevel_cooling_mod, only : nlevel_cooling_report
  use species_mod,     only : species_setup, species_terminal_photoionization_check, &
                              n_elements, elem_nstage, elem_abund, elem_eth
  use memory_mod,      only : destroy_shared_mem_all
  use utility
  use mpi
  implicit none

  type(grid_type)     :: grid
  real(kind=wp)       :: dtime, max_dx, max_dte, dx_vol, dte_vol, dhe_vol(2)
  real(kind=wp)       :: drift_vol(3)
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
  !--- the run's grain optics BEFORE the grid: grid_create_amr multiplies
  !--- par%cext_dust into every leaf's grey rhokap, and par%cext_dust is the
  !--- reference extinction of those very grains (gas_opacity_mod).  Loading
  !--- the table needs nothing but par% and the files it names.
  call dust_reference_extinction()
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
  if (par%use_metals) call species_terminal_photoionization_check(par%eion_max, &
                                                                  ion_eedge, ion_lum)
  !--- a previous run's state, before the first opacity fill reads it
  if (len_trim(par%restart_file) > 0) then
     block
       use previous_state_mod, only : load_previous_state
       call load_previous_state(trim(par%restart_file))
     end block
  end if
  call gas_opacity_setup()
  call jtally_ion_setup(amr_grid%nleaf)
  call ion_score_setup(amr_grid%nleaf)
  if (len_trim(par%restart_file) > 0) then
     block
       use previous_state_mod, only : restore_hei_phi3
       call restore_hei_phi3()
     end block
  end if
  if (trim(par%source_geometry) == 'slab') then
     block
       logical :: ltop, lbot
       ltop = trim(par%slab_faces) == 'top'    .or. trim(par%slab_faces) == 'both'
       lbot = trim(par%slab_faces) == 'bottom' .or. trim(par%slab_faces) == 'both'
       call slab_tally_setup(par%slab_nmu, ion_Ltot, [ltop, lbot])
     end block
  end if
  !--- Milne ground-continuum tables (H, He and, with metals, every metal
  !--- transition: diffuse channels 1-3 and 6) and the band free-free tables
  !--- (channel 7): built once, after ion_setup has fixed the band (the
  !--- tables are truncated at its edges) and after the atomic data are
  !--- loaded.  Read-only during transport.  The metal and free-free channels
  !--- need them whatever par%diffuse_energy_model says about H and He.
  if (trim(par%recombination_continuum) == 'transported') then
     block
       use milne_recomb_spectrum_mod, only : milne_setup
       use recomb_continuum_mod,      only : freefree_band_setup
       call milne_setup()
       call freefree_band_setup()
     end block
  end if
  !--- He I 2^1S two-photon shape: needed whenever the excited-cascade
  !--- channel emits, independently of the ground-continuum model.
  !--- The same two tables serve the on-the-spot He I photons
  !--- (par%hei_on_the_spot).
  if ((trim(par%recombination_continuum) == 'transported' .and. par%hei_diffuse) .or. &
      (trim(par%recombination_continuum) == 'on_the_spot' .and. par%hei_on_the_spot)) then
     block
       use hei_twophoton_mod, only : hei_2ph_setup
       use hei_cascade_mod,   only : hei_cascade_setup
       call hei_2ph_setup()
       call hei_cascade_setup()
     end block
  end if
  !--- He III cascade through He II: the free-bound integral is built over
  !--- the transported band when channel 5 carries it
  !--- (heii_cascade_transported), so it needs par%eion_max fixed too.  The
  !--- on-the-spot He II photons (par%heii_on_the_spot) count the cascade
  !--- from the same tables, built then over the whole spectrum the gas
  !--- emits, whatever the band (heii_cascade_setup).
  if (heii_cascade_transported() .or. &
      (trim(par%recombination_continuum) == 'on_the_spot' .and. par%heii_on_the_spot)) then
     block
       use heii_cascade_mod, only : heii_cascade_setup
       call heii_cascade_setup()
     end block
  end if
  if (par%solve_te) call cooling_setup()
  !--- report the (T, n_e) cooling tables built by species_setup/cooling_setup
  !--- (par%cooling_model = 'local_ne'); a no-op otherwise.
  if (trim(par%cooling_model) == 'local_ne') call nlevel_cooling_report()

  !--- iteration: transport -> rates -> equilibrium -> opacity feedback,
  !--- until the stopping test of par%conv_crit passes ('vol': the drift
  !--- still to come in x_HII and helium below par%gas_tol, see
  !--- composition_converged; 'cell': the largest leaf step below
  !--- par%gas_tol).  gas_niter = 0 = single
  !--- pass (no equilibrium solve, fixed state).  The stellar tally is
  !--- recomputed from zero each iteration (opacity changed).
  !--- the grains of the gas-grain physics (par%grain_pe,
  !--- par%gas_grain_energy_exchange): size bins on the band's energy grid,
  !--- their charge ranges, and the T_dust table they read each iteration
  if (par%solve_te .and. (par%grain_pe .or. par%gas_grain_energy_exchange)) then
     block
       use grain_size_bins_mod,     only : grain_size_bins_setup
       use grain_photoelectric_mod, only : grain_photoelectric_setup
       use dust_temp_mod,           only : dust_temp_setup
       call grain_size_bins_setup()
       call grain_photoelectric_setup()
       call dust_temp_setup()
     end block
  end if
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
        call gas_thermal_update(max_dx, max_dte, dx_vol, dte_vol, dhe_vol, &
                                drift_vol)
        thermal_dt = MPI_WTIME() - thermal_t0
        call record_change_history()
        if (use_vol) then
           converged = composition_converged() .and. dte_vol < par%gas_tol_te
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
           write(6,'(a,3es12.4)') &
              '                  signed vol delta (x_HII, He0-He+, He+-He2+) =', &
              drift_vol
           write(6,'(a,f10.3,a,i0,a)') &
              '                  thermal solve = ', thermal_dt, &
              ' s (', mpar%nproc, ' MPI ranks, all nodes)'
        end if
     else
        call gas_equilibrium_update(max_dx, dx_vol, dhe_vol, drift_vol)
        call record_change_history()
        if (use_vol) then
           converged = composition_converged()
        else
           converged = max_dx < par%gas_tol
        end if
        if (mpar%p_rank == 0) then
           write(6,'(a,i4,2(a,es12.4))') &
              '     iteration ', iter, ': max |delta x_HII| = ', max_dx, &
              ',  vol = ', dx_vol
           write(6,'(a,3es12.4)') &
              '                  signed vol delta (x_HII, He0-He+, He+-He2+) =', &
              drift_vol
        end if
     end if
     !--- distance from the fixed point, H and He (fixed_point_distance_mod):
     !--- the step d, the decay ratio rho fitted over the last NWIN steps,
     !--- the drift still to come e = d rho/(1 - rho), and the mean signed
     !--- step s with its Student statistic t (the floor needs t < 2.365
     !--- with par%conv_drift_test, the default).
     if (mpar%p_rank == 0 .and. par%gas_niter >= 1) then
        call report_fixed_point_distance(' x_HII ', hist_dx)
        call report_fixed_point_distance(' He0-He+', hist_dhe1)
        if (dhe_vol(2) >= 0.0_wp) &
           call report_fixed_point_distance(' He+-He2+', hist_dhe3)
        if (use_vol .and. hist_dx%n < NWIN) write(6,'(a,i0,a)') &
           '                  (the stopping test needs ', NWIN, ' iterations)'
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
        call history_reset(hist_dx)      ! the old grid's steps say nothing
        call history_reset(hist_dhe1)    ! about the new one
        call history_reset(hist_dhe3)
     end if
     if (converged) then
        if (mpar%p_rank == 0) write(6,'(a,i4,a)') &
           '     converged after ', iter, ' iterations'
        exit
     end if
     !--- the iterate, for a resume with par%restart_file (rank 0 writes;
     !--- the header says it is neither converged nor solved in its field)
     if (par%rates_write_interval > 0 .and. iter < niter) then
        if (mod(iter, par%rates_write_interval) == 0) then
           run_converged = .false.
           run_iters     = iter
           run_final_dx  = max_dx
           run_final_dte = max_dte
           run_state_solved_in_field = .false.
           call gas_rates_write(iterate=.true.)
        end if
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

  !--- Final pass and final solve.  The loop above ends on a gas solve, so
  !--- the state it leaves has never been transported through.  One more
  !--- transport of that state (the imaging pass when par%ion_peel is set,
  !--- which rebuilds the same tally with observer peeling) gives the field
  !--- that is written; the gas is then solved once more in that field, so
  !--- the written state -- H, He, T_e, n_e and the metal stages that
  !--- species_write re-derives from the written rates -- is the solution
  !--- of the written field.  Without the final solve the written H/He
  !--- state and n_e belonged to the iterate before the written field while
  !--- the written metal stages belonged to the written field, and before
  !--- convergence the two disagreed (2026-10-12: channel 6 from the written
  !--- state 1.354e36 against 2.017e36 transported, after 2 iterations).
  !--- The field written is then the transport of the state one solve
  !--- earlier; FINALDX/FINALDTE record that final solve's step, i.e. how
  !--- far apart the two are.  Skipped when par%gas_niter < 1 (rates only:
  !--- the written state is the input state, never solved).
  if (par%ion_peel) then
     call ionizing_field_and_rates('imaging pass (peel-off)', with_peel=.true.)
  else if (par%gas_niter >= 1) then
     call ionizing_field_and_rates('final pass')
  end if
  if (par%gas_niter >= 1) then
     if (par%solve_te) then
        call gas_thermal_update(max_dx, max_dte, dx_vol, dte_vol, dhe_vol, &
                                drift_vol)
     else
        call gas_equilibrium_update(max_dx, dx_vol, dhe_vol, drift_vol)
        max_dte = 0.0_wp
     end if
     call hydrogen_balance_residual_grid(eps_hmax, eps_hvol)
     call gas_opacity_fill()
     run_final_dx  = max_dx
     run_final_dte = max_dte
     run_state_solved_in_field = .true.
     if (mpar%p_rank == 0) then
        write(6,'(a,2(a,es12.4))') '     final solve in the written field:', &
           ' max |delta x_HII| = ', max_dx, ',  max |delta Te|/Te = ', max_dte
        write(6,'(a,2(a,es12.4))') '                 ', &
           ' vol |delta x_HII| = ', dx_vol, ',  max H balance residual = ', eps_hmax
     end if
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

contains

  subroutine record_change_history()
    call history_push(hist_dx, dx_vol, drift_vol(1))
    call history_push(hist_dhe1, dhe_vol(1), drift_vol(2))
    if (dhe_vol(2) >= 0.0_wp) then
       call history_push(hist_dhe3, dhe_vol(2), drift_vol(3))
    else
       call history_reset(hist_dhe3)
    end if
  end subroutine record_change_history

  !--- par%conv_crit = 'vol': the ionization state has converged when the
  !--- drift still to come is below par%gas_tol for x_HII, for He0 <-> He+,
  !--- and for He+ <-> He2+ where He2+ is at least 1% of the ionized helium,
  !--- or each sits on its Monte Carlo floor with no net drift of the
  !--- ionized volume (par%conv_drift_test; fixed_point_distance_mod; the test on the step alone
  !--- stopped the HII40 benchmark with its He I front 0.17 pc short, and
  !--- a floor read from the step ratios alone stopped a restart after 4
  !--- iterations with its H front still moving).
  logical function composition_converged() result(done)
    done = fixed_point_converged(hist_dx, par%gas_tol, par%conv_drift_test) .and. &
           fixed_point_converged(hist_dhe1, par%gas_tol, par%conv_drift_test)
    if (dhe_vol(2) >= 0.0_wp) &
       done = done .and. fixed_point_converged(hist_dhe3, par%gas_tol, &
                                               par%conv_drift_test)
  end function composition_converged

  !--- one log line: the latest step d, the decay ratio rho, the drift still
  !--- to come e (or why e is not available), and the mean signed step s
  !--- over the window with its Student statistic t
  subroutine report_fixed_point_distance(label, h)
    character(len=*),        intent(in) :: label
    type(iteration_history), intent(in) :: h
    real(kind=wp) :: rho, dist, smean, tstat
    logical :: available
    character(len=40) :: tail
    call fixed_point_distance(h, rho, dist, available)
    call net_drift(h, smean, tstat)
    write(tail,'(a,es10.2,a,f7.2)') ' s = ', smean, ' t = ', min(tstat, 9999.0_wp)
    if (available) then
       write(6,'(a,a9,3(a,es11.3),a,a)') '                  fixed point', label, &
          ': d = ', h%d(h%n), ',  rho = ', rho, ',  e = ', dist, ', ', trim(tail)
    else if (h%n < 4) then
       write(6,'(a,a9,a,es11.3,a,a)') '                  fixed point', label, &
          ': d = ', h%d(h%n), ',  (rho needs 4 iterations), ', trim(tail)
    else
       write(6,'(a,a9,2(a,es11.3),a,a)') '                  fixed point', label, &
          ': d = ', h%d(h%n), ',  rho = ', rho, ',  (not decaying), ', trim(tail)
    end if
  end subroutine report_fixed_point_distance

end program main
