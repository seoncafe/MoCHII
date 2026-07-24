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
  use ion_band_mod,    only : ion_setup, gen_ion_photon, gen_ion_photon_qmc, &
                              ion_Ltot, ion_qmc_ndim
  use qmc_mod,         only : qmc_uniforms, qmc_uniforms_stream, &
                              QMC_MAXDIM, QMC_STREAM_STELLAR, QMC_STREAM_DIFFUSE
  use gas_opacity_mod, only : gas_opacity_setup, gas_opacity_fill
  use jtally_mod,      only : jtally_ion_setup, jtally_ion_reduce, jt_ion, &
                              slab_tally_setup, slab_tally_reduce, &
                              slab_write_Imu, slab_tally_on, slab_Iesc
  use raytrace_amr_mod,only : transport_ion_packet, slab_walk_report
  use gas_rates_mod,   only : gas_rates_compute, gas_rates_write, &
                              secion_apply, &
                              run_converged, run_iters, run_final_dx, &
                              run_final_dte
  use ion_balance_mod, only : gas_equilibrium_update
  use thermal_mod,     only : gas_thermal_update
  use cooling_mod,     only : cooling_setup
  use species_mod,     only : species_setup, species_gamma_compute, &
                              n_elements, elem_nstage, elem_abund, elem_eth
  use diffuse_mod,     only : diffuse_build, gen_diffuse_photon, &
                              gen_diffuse_photon_qmc, diffuse_nphot, diffuse_lum
  use memory_mod,      only : destroy_shared_mem_all
  use utility
  use mpi
  implicit none

  type(grid_type)     :: grid
  type(photon_type)   :: photon
  integer(kind=int64) :: ip, n_done, n_step
  real(kind=wp)       :: dtime, max_dx, max_dte, dx_vol, dte_vol
  real(kind=wp)       :: u_launch(QMC_MAXDIM)   ! stellar (1:nd_qmc) + diffuse (1:9) launch uniforms
  logical             :: converged, use_vol, use_sobol
  integer             :: ierr, iter, niter, nd_qmc
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
  !--- does not depend on the band (that use is inside the per-leaf loops).
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
  if (trim(par%source_geometry) == 'slab') then
     block
       logical :: ltop, lbot
       ltop = trim(par%slab_faces) == 'top'    .or. trim(par%slab_faces) == 'both'
       lbot = trim(par%slab_faces) == 'bottom' .or. trim(par%slab_faces) == 'both'
       call slab_tally_setup(par%slab_nmu, ion_Ltot, [ltop, lbot])
     end block
  end if
  if (par%solve_te) call cooling_setup()

  !--- iteration: transport -> rates -> equilibrium -> opacity feedback,
  !--- until max |delta x_HII| < par%gas_tol.  gas_niter = 0 = single
  !--- pass (no equilibrium solve, fixed state).  The stellar tally is
  !--- recomputed from zero each iteration (opacity changed).
  niter = max(par%gas_niter, 1)
  n_step = max(1_int64, int(par%nprint,int64)/mpar%nproc)
  !--- quasi-random launch (single point source): the scrambled Sobol point is
  !--- indexed by the GLOBAL photon number ip-1 (0-based), so the launch set is
  !--- independent of the MPI task count and fixed across iterations.
  use_sobol = trim(par%launch_sequence) == 'sobol'
  !--- stellar launch dimension for this configuration (3 for a single point
  !--- source, 7 for the fixed superset layout); the diffuse launch uses 9.
  nd_qmc = ion_qmc_ndim()
  converged = .false.
  max_dx = 0.0_wp;  max_dte = 0.0_wp;  dx_vol = 0.0_wp;  dte_vol = 0.0_wp
  do iter = 1, niter
     jt_ion(:,:) = 0.0_wp
     if (slab_tally_on) slab_Iesc(:,:) = 0.0_wp    ! keep the last iteration's escaping field
     call time_stamp(dtime)
     if (mpar%p_rank == 0) write(6,'(a,i4,a,f8.3,a)') &
        '---> iteration ', iter, ': ionizing-band transport...  @ ', &
        dtime/60.0_wp, ' mins'
     n_done = 0
     do ip = mpar%p_rank+1, par%nphotons, mpar%nproc
        if (use_sobol) then
           call qmc_uniforms(ip-1_int64, u_launch(1:nd_qmc))
           call gen_ion_photon_qmc(photon, u_launch(1:nd_qmc))
        else
           call gen_ion_photon(photon)
        end if
        photon%id      = ip
        photon%istream = QMC_STREAM_STELLAR
        call transport_ion_packet(photon)
        n_done = n_done + 1
        if (mpar%p_rank == 0 .and. mod(n_done, n_step) == 0) then
           call time_stamp(dtime)
           write(6,'(es14.3,a,f8.3,a)') real(n_done,wp)*mpar%nproc, &
              ' photons  @ ', dtime/60.0_wp, ' mins'
        endif
     end do
     !--- diffuse ground-recombination packets from the current state.
     if (par%diffuse_field) then
        call diffuse_build(ion_Ltot/real(par%nphotons, wp))
        if (mpar%p_rank == 0) write(6,'(a,es12.4,a,i12,a)') &
           '     diffuse field: L = ', diffuse_lum, ' erg/s, ', &
           diffuse_nphot, ' packets'
        !--- diffuse packets ride the SECOND (decorrelated) Sobol stream,
        !--- indexed by the global diffuse photon number (ip-1); diffuse_nphot
        !--- varies per iteration, and a Sobol prefix of any length is balanced.
        do ip = mpar%p_rank+1, diffuse_nphot, mpar%nproc
           if (use_sobol) then
              call qmc_uniforms_stream(ip-1_int64, u_launch(1:9), QMC_STREAM_DIFFUSE)
              call gen_diffuse_photon_qmc(photon, u_launch(1:9))
           else
              call gen_diffuse_photon(photon)
           end if
           photon%id      = ip
           photon%istream = QMC_STREAM_DIFFUSE
           call transport_ion_packet(photon)
        end do
     end if
     call jtally_ion_reduce()
     call gas_rates_compute()
     if (par%use_metals) call species_gamma_compute()
     if (par%use_sec_ion) call secion_apply()
     if (par%gas_niter < 1) exit          ! rates only, no solve
     use_vol = trim(par%conv_crit) == 'vol'
     if (par%solve_te) then
        call gas_thermal_update(max_dx, max_dte, dx_vol, dte_vol)
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

  !--- peel-off imaging pass: one extra transport of the CONVERGED state
  !--- with observer peeling (direct at emission for stellar and diffuse
  !--- packets; dust-scattered at every interaction via the hook).  The
  !--- band tally is rebuilt by this pass, so the written rates reflect
  !--- the same field the images carry.
  if (par%ion_peel) then
     block
       use ion_peel_mod, only : ion_peel_setup, ion_peel_direct, &
                                ion_peel_scatter
       call ion_peel_setup()
       if (par%ion_add_dust .and. par%ion_dust_scatter) &
          ion_peel_scatter_hook => ion_peel_scatter
       jt_ion(:,:) = 0.0_wp
       call time_stamp(dtime)
       if (mpar%p_rank == 0) write(6,'(a,f8.3,a)') &
          '---> imaging pass (peel-off)...  @ ', dtime/60.0_wp, ' mins'
       do ip = mpar%p_rank+1, par%nphotons, mpar%nproc
          if (use_sobol) then
             call qmc_uniforms(ip-1_int64, u_launch(1:nd_qmc))
             call gen_ion_photon_qmc(photon, u_launch(1:nd_qmc))
          else
             call gen_ion_photon(photon)
          end if
          photon%id      = ip
          photon%istream = QMC_STREAM_STELLAR
          call ion_peel_direct(photon)
          call transport_ion_packet(photon)
       end do
       if (par%diffuse_field) then
          call diffuse_build(ion_Ltot/real(par%nphotons, wp))
          do ip = mpar%p_rank+1, diffuse_nphot, mpar%nproc
             if (use_sobol) then
                call qmc_uniforms_stream(ip-1_int64, u_launch(1:9), QMC_STREAM_DIFFUSE)
                call gen_diffuse_photon_qmc(photon, u_launch(1:9))
             else
                call gen_diffuse_photon(photon)
             end if
             photon%id      = ip
             photon%istream = QMC_STREAM_DIFFUSE
             call ion_peel_direct(photon)
             call transport_ion_packet(photon)
          end do
       end if
       ion_peel_scatter_hook => null()
       call jtally_ion_reduce()
       call gas_rates_compute()
       if (par%use_metals) call species_gamma_compute()
       if (par%use_sec_ion) call secion_apply()
     end block
  end if

  !--- output
  if (par%ion_add_dust) then
     block
       use dust_temp_mod, only : dust_temp_setup, dust_temp_compute
       use gas_rates_mod, only : heat_dust
       call dust_temp_setup()
       call dust_temp_compute(heat_dust)
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
  if (par%use_metals .and. par%solve_te) then
     block
       use lines_mod, only : lines_write
       call lines_write()
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
