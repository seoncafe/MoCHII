module ionizing_field_mod
!---------------------------------------------------------------------------
! MoCHII: one pass of the ionizing band.
!
! Transport every stellar and diffuse packet through the CURRENT opacity,
! reduce the path-length tally into 4 pi J_nu dnu, and turn that into the
! photoionization and heating rate integrals.  The nonlinear iteration, the
! final consistency pass and the peel-off imaging pass all come through here,
! so the radiation field the written rates carry is built exactly as the
! field that drove the iteration -- same launch set, same diffuse rebuild,
! same reduction.
!
! With with_peel every emitted packet (stellar AND diffuse) also peels a
! direct contribution toward each observer, and -- when dust scattering is
! sampled -- every interaction peels a scattered contribution through the
! hook installed here for the duration of the pass.
!---------------------------------------------------------------------------
  use define
  use ion_band_mod,     only : gen_ion_photon, gen_ion_photon_qmc, &
                               ion_Ltot, ion_qmc_ndim
  use qmc_mod,          only : qmc_uniforms, qmc_uniforms_stream, &
                               QMC_MAXDIM, QMC_STREAM_STELLAR, QMC_STREAM_DIFFUSE
  use jtally_mod,       only : jtally_ion_reduce, jt_ion, &
                               slab_tally_on, slab_Iesc
  use ion_score_mod,    only : ion_score_reset, ion_score_reduce, &
                               ion_score_report_energy_closure, &
                               ion_score_report_metal_rates
  use raytrace_amr_mod, only : transport_ion_packet
  use gas_rates_mod,    only : gas_rates_compute, secion_apply
  use species_mod,      only : species_gamma_compute
  use diffuse_mod,      only : diffuse_build, gen_diffuse_photon, &
                               gen_diffuse_photon_qmc, diffuse_nphot, diffuse_lum
  use ion_peel_mod,     only : ion_peel_setup, ion_peel_direct, ion_peel_scatter
  use photon_schedule_mod, only : photon_schedule_begin, photon_schedule_next, &
                                  photon_schedule_end
  use utility,          only : time_stamp
  implicit none
  private

  public :: ionizing_field_setup, ionizing_field_and_rates

  !--- progress-report stride, launch sequence and stellar launch dimension;
  !--- fixed once by ionizing_field_setup and read-only during transport.
  integer(kind=int64) :: n_step    = 1_int64
  logical             :: use_sobol = .false.
  integer             :: nd_qmc    = 0

contains

  !=========================================================================
  ! Launch configuration shared by every pass.  Called after ion_setup has
  ! fixed the band (ion_qmc_ndim queries the source layout).
  !=========================================================================
  subroutine ionizing_field_setup()

    n_step = max(1_int64, int(par%nprint,int64)/mpar%nproc)
    !--- quasi-random launch (single point source): the scrambled Sobol point is
    !--- indexed by the GLOBAL photon number ip-1 (0-based), so the launch set is
    !--- independent of the MPI task count and fixed across iterations.
    use_sobol = trim(par%launch_sequence) == 'sobol'
    !--- stellar launch dimension for this configuration (3 for a single point
    !--- source, 7 for the fixed superset layout); the diffuse launch uses 9.
    nd_qmc = ion_qmc_ndim()
  end subroutine ionizing_field_setup

  !=========================================================================
  ! One transport pass and the rate integrals it feeds.  with_peel adds the
  ! observer peeling (direct at emission for stellar and diffuse packets;
  ! dust-scattered at every interaction via the hook) and suppresses the
  ! progress and diffuse-luminosity reports of the iterating passes.
  !=========================================================================
  subroutine ionizing_field_and_rates(label, with_peel)
    character(len=*), intent(in)  :: label
    logical, intent(in), optional :: with_peel

    type(photon_type)   :: photon
    integer(kind=int64) :: ip, n_done
    real(kind=wp)       :: u_launch(QMC_MAXDIM)   ! stellar (1:nd_qmc) + diffuse (1:9) launch uniforms
    real(kind=wp)       :: dtime
    logical             :: peel

    peel = .false.
    if (present(with_peel)) peel = with_peel

    if (peel) then
       call ion_peel_setup()
       if (par%ion_add_dust .and. par%ion_dust_scatter) &
          ion_peel_scatter_hook => ion_peel_scatter
    end if
    jt_ion(:,:) = 0.0_wp
    call ion_score_reset()
    !--- keep the last pass's escaping field: an iterating pass starts the
    !--- slab tally from zero, the imaging pass adds to what is already there.
    if (.not. peel .and. slab_tally_on) slab_Iesc(:,:) = 0.0_wp
    call time_stamp(dtime)
    if (mpar%p_rank == 0) then
       if (peel) then
          write(6,'(3a,f8.3,a)') &
             '---> ', trim(label), '...  @ ', dtime/60.0_wp, ' mins'
       else
          write(6,'(3a,f8.3,a)') &
             '---> ', trim(label), ': ionizing-band transport...  @ ', &
             dtime/60.0_wp, ' mins'
       end if
    end if
    n_done = 0
    !--- how the photon indices are split over the ranks (static cyclic or
    !--- served on demand) is par%use_master_slave; the body below is the
    !--- same either way.  An iterating pass also asks the dynamic master to
    !--- print the progress line, since it transports nothing itself and the
    !--- count below would never advance on rank 0.
    if (peel) then
       call photon_schedule_begin(par%nphotons)
    else
       call photon_schedule_begin(par%nphotons, n_step*int(mpar%nproc,int64))
    end if
    do while (photon_schedule_next(ip))
       if (use_sobol) then
          call qmc_uniforms(ip-1_int64, u_launch(1:nd_qmc))
          call gen_ion_photon_qmc(photon, u_launch(1:nd_qmc))
       else
          call gen_ion_photon(photon)
       end if
       photon%id      = ip
       photon%istream = QMC_STREAM_STELLAR
       if (peel) call ion_peel_direct(photon)
       call transport_ion_packet(photon)
       if (.not. peel) then
          n_done = n_done + 1
          if (mpar%p_rank == 0 .and. mod(n_done, n_step) == 0) then
             call time_stamp(dtime)
             write(6,'(es14.3,a,f8.3,a)') real(n_done,wp)*mpar%nproc, &
                ' photons  @ ', dtime/60.0_wp, ' mins'
          endif
       end if
    end do
    call photon_schedule_end()
    !--- diffuse ground-recombination packets from the current state.
    if (par%diffuse_field) then
       call diffuse_build(ion_Ltot/real(par%nphotons, wp))
       if (.not. peel .and. mpar%p_rank == 0) write(6,'(a,es12.4,a,i12,a)') &
          '     diffuse field: L = ', diffuse_lum, ' erg/s, ', &
          diffuse_nphot, ' packets'
       !--- diffuse packets ride the SECOND (decorrelated) Sobol stream,
       !--- indexed by the global diffuse photon number (ip-1); diffuse_nphot
       !--- varies per pass, and a Sobol prefix of any length is balanced.
       call photon_schedule_begin(diffuse_nphot)
       do while (photon_schedule_next(ip))
          if (use_sobol) then
             call qmc_uniforms_stream(ip-1_int64, u_launch(1:9), QMC_STREAM_DIFFUSE)
             call gen_diffuse_photon_qmc(photon, u_launch(1:9))
          else
             call gen_diffuse_photon(photon)
          end if
          photon%id      = ip
          photon%istream = QMC_STREAM_DIFFUSE
          if (peel) call ion_peel_direct(photon)
          call transport_ion_packet(photon)
       end do
       call photon_schedule_end()
    end if
    if (peel) ion_peel_scatter_hook => null()
    call jtally_ion_reduce()
    call ion_score_reduce()
    call gas_rates_compute()
    call ion_score_report_energy_closure()
    if (par%use_metals) call species_gamma_compute()
    call ion_score_report_metal_rates()
    if (par%use_sec_ion) call secion_apply()
  end subroutine ionizing_field_and_rates

end module ionizing_field_mod
