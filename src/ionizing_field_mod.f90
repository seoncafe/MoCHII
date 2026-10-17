module ionizing_field_mod
!---------------------------------------------------------------------------
! MoCHII: one pass of the ionizing band.
!
! Transport every stellar and diffuse packet through the CURRENT opacity,
! reduce the path-length tally into 4 pi J_nu dnu, and turn that into the
! photoionization and heating rate integrals.  The nonlinear iteration, the
! final pass (before the final gas solve) and the peel-off imaging pass all
! come through here, so the radiation field the written rates carry is built exactly as the
! field that drove the iteration -- same launch set, same diffuse rebuild,
! same reduction.
!
! With par%add_fuv and par%fuv_sampling = 'separate' the source packets are
! launched in two strata, the ionizing part of the band (par%nphotons packets
! on the stellar Sobol stream, the same indices as without the FUV band) and
! the FUV part (par%nphotons_fuv packets on their own stream); each carries its
! stratum's luminosity over its stratum's count (ion_band_mod).  Every
! estimator downstream sums Lpacket * wgt * path, so nothing past the launch
! needs to know which stratum a packet came from.
!
! With with_peel every emitted packet (stellar AND diffuse) also peels a
! direct contribution toward each observer, and -- when dust scattering is
! sampled -- every interaction peels a scattered contribution through the
! hook installed here for the duration of the pass.
!---------------------------------------------------------------------------
  use define
  use ion_band_mod,     only : gen_ion_photon, gen_ion_photon_qmc, &
                               ion_Ltot, ion_qmc_ndim, fuv_stratified, &
                               ion_L_ionizing, ion_L_fuv, &
                               band_packet_count, ionizing_packet_luminosity, &
                               LAUNCH_BAND_WHOLE, LAUNCH_BAND_IONIZING, &
                               LAUNCH_BAND_FUV
  use qmc_mod,          only : qmc_uniforms_stream, QMC_MAXDIM, &
                               QMC_STREAM_STELLAR, QMC_STREAM_DIFFUSE, &
                               QMC_STREAM_FUV
  use jtally_mod,       only : jtally_ion_reduce, jt_ion, jt_ion_on, &
                               slab_tally_on, slab_Iesc
  use ion_score_mod,    only : ion_score_reset, ion_score_reduce, &
                               ion_score_report_energy_closure, &
                               ion_score_report_ledger, &
                               ion_score_report_metal_rates
  use raytrace_amr_mod, only : transport_ion_packet
  use gas_rates_mod,    only : gas_rates_compute
  use species_mod,      only : species_gamma_compute
  use diffuse_mod,      only : diffuse_build, gen_diffuse_photon, &
                               gen_diffuse_photon_qmc, diffuse_nphot, diffuse_lum, &
                               diffuse_channel_luminosity, metal_valence_luminosity_outside
  use heii_cascade_mod, only : heii_cascade_transported
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
  !--- the launched-luminosity check of the band strata is printed once, on
  !--- the first pass that has an FUV extension.
  logical             :: strata_reported = .false.
  integer, parameter  :: LEDGER_WP = selected_real_kind(18, 300)

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
    real(kind=LEDGER_WP) :: launched(2)   ! launched Sum Lpacket*wgt: stratum 1 (or whole band), FUV
    logical             :: peel

    peel = .false.
    if (present(with_peel)) peel = with_peel

    if (peel) then
       call ion_peel_setup()
       if (par%ion_add_dust .and. par%ion_dust_scatter) &
          ion_peel_scatter_hook => ion_peel_scatter
    end if
    if (jt_ion_on) jt_ion(:,:) = 0.0_wp
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
    launched = 0.0_LEDGER_WP
    !--- how the photon indices are split over the ranks (static cyclic or
    !--- served on demand) is par%use_master_slave; the body below is the
    !--- same either way.  An iterating pass also asks the dynamic master to
    !--- print the progress line, since it transports nothing itself and the
    !--- count below would never advance on rank 0.
    if (fuv_stratified) then
       call launch_band_part(LAUNCH_BAND_IONIZING, QMC_STREAM_STELLAR, &
                             band_packet_count(LAUNCH_BAND_IONIZING), launched(1))
       call launch_band_part(LAUNCH_BAND_FUV, QMC_STREAM_FUV, &
                             band_packet_count(LAUNCH_BAND_FUV), launched(2))
    else
       call launch_band_part(LAUNCH_BAND_WHOLE, QMC_STREAM_STELLAR, &
                             band_packet_count(LAUNCH_BAND_WHOLE), launched(1))
    end if
    if (par%add_fuv .and. .not. strata_reported) call report_launched_luminosity(launched)
    !--- diffuse recombination packets from the current state: every channel
    !--- with the transported continuum, the He III cascade alone (channel 5)
    !--- with the on-the-spot one and par%heii_diffuse.
    if (trim(par%recombination_continuum) == 'transported' .or. &
        heii_cascade_transported()) then
       !--- the diffuse packets carry the luminosity of one IONIZING source
       !--- packet (ionizing_packet_luminosity: ion_Ltot/nphotons for the
       !--- whole band, L_ion/nphotons with the strata), so the recombination
       !--- field is sampled as finely as the field that drives it.
       call diffuse_build(ionizing_packet_luminosity())
       if (.not. peel .and. mpar%p_rank == 0) then
          write(6,'(a,es12.4,a,i12,a)') &
             '     diffuse field: L = ', diffuse_lum, ' erg/s, ', &
             diffuse_nphot, ' packets'
          !--- by channel: H I, He I, He II ground; He I, He III cascades;
          !--- metal ground continua; free-free; metal valence-fraction
          !--- photons (par%metal_recombination_photons)
          write(6,'(a,8es11.3)') '     diffuse L by channel 1-8 [erg/s]:', &
             diffuse_channel_luminosity
          if (trim(par%metal_recombination_photons) /= 'none') &
             write(6,'(a,2es11.3)') '     channel 8 not launched, below '// &
                'eion_min / at or above eion_max [erg/s]:', &
                metal_valence_luminosity_outside
       end if
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
    if (.not. peel) call ion_score_report_ledger()
    if (par%use_metals) call species_gamma_compute()
    call ion_score_report_metal_rates()

  contains

    !--- one launch loop over npk packets of band part iband, on the Sobol
    !--- stream istream (indices 0..npk-1) or the Mersenne Twister; adds the
    !--- launched luminosity sum Lpacket*wgt to Lsum.
    subroutine launch_band_part(iband, istream, npk, Lsum)
      integer,                  intent(in)    :: iband, istream
      integer(kind=int64),      intent(in)    :: npk
      real(kind=LEDGER_WP),     intent(inout) :: Lsum
      if (peel) then
         call photon_schedule_begin(npk)
      else
         call photon_schedule_begin(npk, n_step*int(mpar%nproc,int64))
      end if
      do while (photon_schedule_next(ip))
         if (use_sobol) then
            call qmc_uniforms_stream(ip-1_int64, u_launch(1:nd_qmc), istream)
            call gen_ion_photon_qmc(photon, u_launch(1:nd_qmc), iband)
         else
            call gen_ion_photon(photon, iband)
         end if
         photon%id      = ip
         photon%istream = istream
         photon%source_channel = merge(9, 0, photon%energy_eV < par%eion_min)
         Lsum = Lsum + real(photon%Lpacket*photon%wgt, LEDGER_WP)
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
    end subroutine launch_band_part

  end subroutine ionizing_field_and_rates

  !=========================================================================
  ! Launched luminosity of each band part against what it must carry:
  ! stratified, the two sums against L(eion_min-eion_max) and
  ! L(efuv_min-eion_min) and their total against the band luminosity
  ! ion_Ltot; whole band, the one sum against ion_Ltot.  By construction a
  ! stratum launches N packets of L/N, so the ratios are 1 to rounding
  ! (gate tests/fuv_stratified).  Printed once.
  !=========================================================================
  subroutine report_launched_luminosity(launched)
    use mpi
    real(kind=LEDGER_WP), intent(in) :: launched(2)
    real(kind=wp) :: tot(2)
    integer :: ierr
    tot = real(launched, wp)
    call MPI_ALLREDUCE(MPI_IN_PLACE, tot, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    strata_reported = .true.
    if (mpar%p_rank /= 0) return
    if (fuv_stratified) then
       write(*,'(a,es23.15,a,es10.2)') ' ION: launched L, ionizing stratum = ', &
          tot(1), ',  /L(eion_min-eion_max) - 1 = ', rel_dev(tot(1), ion_L_ionizing)
       write(*,'(a,es23.15,a,es10.2)') ' ION: launched L, FUV stratum      = ', &
          tot(2), ',  /L(efuv_min-eion_min) - 1 = ', rel_dev(tot(2), ion_L_fuv)
       write(*,'(a,es23.15,a,es10.2)') ' ION: launched L, both strata      = ', &
          sum(tot), ',  /L_band - 1               = ', rel_dev(sum(tot), ion_Ltot)
    else
       write(*,'(a,es23.15,a,es10.2)') ' ION: launched L, whole band       = ', &
          tot(1), ',  /L_band - 1               = ', rel_dev(tot(1), ion_Ltot)
    end if
  contains
    real(kind=wp) function rel_dev(a, b)
      real(kind=wp), intent(in) :: a, b
      rel_dev = 0.0_wp
      if (b > 0.0_wp) rel_dev = a/b - 1.0_wp
    end function rel_dev
  end subroutine report_launched_luminosity

end module ionizing_field_mod
