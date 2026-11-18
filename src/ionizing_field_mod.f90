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
!
! Under par%diffuse_launch = 'absorption_reemission' the diffuse emission of
! the pass (diffuse_build: the state emission L_c, the re-emission energy
! eps_s/max(psi_s, 1) of each absorber and the explicit launch X_c, from the
! frozen state and the photoionizations of the previous pass's field) is
! tabulated BEFORE the first source packet and before ion_score_reset (psi
! comes from the previous pass's estimators), since the source packets
! already re-emit at their gas absorptions; the
! remainder is then launched after the sources by the same diffuse loop.  The
! default mode keeps its order (sources, then diffuse_build and the diffuse
! loop).  diffuse_build draws no random numbers.
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
  use jtally_mod,       only : jtally_ion_reduce, jtally_ion_reset, jt_ion_on, &
                               slab_tally_on, slab_Iesc
  use ion_score_mod,    only : ion_score_reset, ion_score_reduce, &
                               ion_score_report_energy_closure, &
                               ion_score_report_ledger, &
                               ion_score_report_metal_rates, &
                               ion_score_report_reemission
  use raytrace_amr_mod, only : transport_ion_packet
  use gas_rates_mod,    only : gas_rates_compute
  use species_mod,      only : species_gamma_compute
  use diffuse_mod,      only : diffuse_build, gen_diffuse_photon, &
                               gen_diffuse_photon_qmc, diffuse_nphot, diffuse_lum, &
                               diffuse_channel_luminosity, metal_valence_luminosity_outside, &
                               diffuse_remainder_luminosity
  use heii_cascade_mod, only : heii_cascade_transported
  use ion_peel_mod,     only : ion_peel_setup, ion_peel_direct, ion_peel_scatter, &
                               ion_peel_reemission
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
    logical             :: peel, reemission

    peel = .false.
    if (present(with_peel)) peel = with_peel
    reemission = trim(par%diffuse_launch) == 'absorption_reemission'

    if (peel) then
       call ion_peel_setup()
       if (par%ion_add_dust .and. par%ion_dust_scatter) &
          ion_peel_scatter_hook => ion_peel_scatter
       ! Gas re-emission also peels in dust-free calculations.
       if (reemission) ion_peel_reemission_hook => ion_peel_reemission
    end if
    call jtally_ion_reset()
    ! Freeze absorber yields and explicit emission before the first packet,
    ! using the shared previous-pass field before ion_score_reset invalidates it.
    if (reemission) call diffuse_build(ionizing_packet_luminosity())
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
    if (reemission) call report_diffuse_emission()
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
       if (.not. reemission) call tabulate_diffuse_emission()
       !--- diffuse packets ride the SECOND (decorrelated) Sobol stream,
       !--- indexed by the global diffuse photon number (ip-1); diffuse_nphot
       !--- varies with each pass. Arbitrary prefixes are usable, but the
       !--- digital-net balance guarantee requires a power-of-two count.
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
    ion_peel_reemission_hook => null()
    call jtally_ion_reduce()
    call ion_score_reduce()
    call gas_rates_compute()
    if (trim(par%recombination_continuum) == 'transported' .and. par%hei_diffuse) then
       block
         use diffuse_mod, only : report_hei_584_absorption
         call report_hei_584_absorption()
       end block
    end if
    call ion_score_report_energy_closure()
    !--- every pass, the imaging pass included: with par%ion_peel it is the
    !--- last pass, and its ledger is the one written to the rates file.
    call ion_score_report_ledger(label)
    if (reemission) then
       call ion_score_report_reemission(diffuse_channel_luminosity, &
                                        diffuse_remainder_luminosity)
       call report_local_residuals()
    end if
    if (par%use_metals) call species_gamma_compute()
    call ion_score_report_metal_rates()

  contains

    !--- the diffuse emission of the pass from the current state
    !--- (diffuse_build) and its report.  The diffuse packets carry the
    !--- luminosity of one IONIZING source packet (ionizing_packet_luminosity:
    !--- ion_Ltot/nphotons for the whole band, L_ion/nphotons with the
    !--- strata), so the recombination field is sampled as finely as the
    !--- field that drives it; under absorption_reemission that holds for
    !--- the launched remainder Sum X.  The count is
    !--- n = min(par%nphotons, max(1, nint(L / Lpacket))) packets of L / n
    !--- (diffuse_build), L being the diffuse luminosity or Sum X.
    subroutine tabulate_diffuse_emission()
      call diffuse_build(ionizing_packet_luminosity())
      call report_diffuse_emission()
    end subroutine tabulate_diffuse_emission

    !--- the report of the last diffuse_build (rank 0, iterating passes)
    subroutine report_diffuse_emission()
      if (.not. peel .and. mpar%p_rank == 0) then
         if (reemission) then
            write(6,'(a,es12.4,a,es12.4,a,i12,a)') &
               '     diffuse field: state emission L = ', &
               sum(diffuse_channel_luminosity), ' erg/s, explicit remainder X = ', &
               diffuse_lum, ' erg/s in ', diffuse_nphot, ' packets'
         else
            write(6,'(a,es12.4,a,i12,a)') &
               '     diffuse field: L = ', diffuse_lum, ' erg/s, ', &
               diffuse_nphot, ' packets'
         end if
         !--- by channel: H I, He I, He II ground; He I, He III cascades;
         !--- metal ground continua; free-free; metal valence-fraction
         !--- photons (par%metal_recombination_photons)
         write(6,'(a,8es11.3)') '     diffuse L by channel 1-8 [erg/s]:', &
            diffuse_channel_luminosity
         if (reemission) then
            write(6,'(a,8es11.3)') &
               '     remainder X by channel 1-8 [erg/s]:', diffuse_remainder_luminosity
            call report_absorber_psi()
         end if
         if (trim(par%metal_recombination_photons) /= 'none') &
            write(6,'(a,2es11.3)') '     channel 8 not launched, below '// &
               'eion_min / at or above eion_max [erg/s]:', &
               metal_valence_luminosity_outside
      end if
    end subroutine report_diffuse_emission

    !--- emission of the pass against the state emission, by channel and
    !--- leaf (rank 0), plan Section 12.4: with L_c,i the state emission,
    !--- X_c,i the explicit launch (its expectation; the launch samples it)
    !--- and Re_c,i the sampled re-emission into channel c of absorbing leaf
    !--- i, R_c,i = X_c,i + Re_c,i.  By channel: Sum L_c, Sum X_c, Sum Re_c
    !--- and (Re_c + X_c)/L_c; over all (c, i): Sum |R - L| / Sum L, the
    !--- L-weighted rms of (R - L)/L over L > 0, and max |R - L| / Sum L.
    subroutine report_local_residuals()
      use diffuse_mod,   only : diffuse_leaf_emission, diffuse_explicit_leaf_emission
      use ion_score_mod, only : ion_score_reemission_leaf, NDIFFUSE_CHANNEL
      use gas_state_mod, only : gas_nleaf
      real(kind=wp), allocatable :: lem(:,:), xem(:,:), reem(:,:), gev(:)
      real(kind=wp) :: lsum, wl, d, sabs, srms, smax, ratio(NDIFFUSE_CHANNEL)
      integer :: il, c
      if (peel .or. mpar%p_rank /= 0) return
      allocate(lem(NDIFFUSE_CHANNEL,gas_nleaf), xem(NDIFFUSE_CHANNEL,gas_nleaf), &
               reem(NDIFFUSE_CHANNEL,gas_nleaf), gev(gas_nleaf))
      call diffuse_leaf_emission(lem)
      call diffuse_explicit_leaf_emission(xem)
      call ion_score_reemission_leaf(reem, gev)
      lsum = sum(lem)
      if (lsum <= 0.0_wp) return
      wl = 0.0_wp;  sabs = 0.0_wp;  srms = 0.0_wp;  smax = 0.0_wp
      do il = 1, gas_nleaf
         do c = 1, NDIFFUSE_CHANNEL
            d = xem(c,il) + reem(c,il) - lem(c,il)
            sabs = sabs + abs(d)
            smax = max(smax, abs(d))
            if (lem(c,il) > 0.0_wp) then
               wl = wl + lem(c,il)
               srms = srms + d*d/lem(c,il)
            end if
         end do
      end do
      do c = 1, NDIFFUSE_CHANNEL
         ratio(c) = 0.0_wp
         if (sum(lem(c,:)) > 0.0_wp) ratio(c) = (sum(xem(c,:)) + sum(reem(c,:)))/sum(lem(c,:))
      end do
      write(6,'(a,8es11.3)') ' REEMIT: (Re_c + X_c)/L_c by channel 1-8 :', ratio
      call report_count_residuals(lem, xem, reem, lsum)
      write(6,'(a,es11.3,a,es11.3,a,es11.3)') &
         ' REEMIT: local R - L over (channel, leaf): Sum|.|/Sum L = ', sabs/lsum, &
         ', L-weighted rms of ./L = ', sqrt(srms/wl), ', max|.|/Sum L = ', smax/lsum
    end subroutine report_local_residuals

    !--- count-based residual of the re-emission (rank 0, after the
    !--- reduction): for c in C(s), Re_c,i - b_c,i N_now,s,i with
    !--- b_c = f_s L_c / max(N_rec,s, N_abs,s) and N_now,s the photoionizations
    !--- of s in this pass (raw path sums, the pass's own state) -- zero in
    !--- expectation, the sampling noise of the re-emission alone --
    !--- as Sum|.|/Sum L, L-weighted rms of ./L_c and max|.|/Sum L over
    !--- (c, i), and the signed total over Sum L
    subroutine report_count_residuals(lem, xem, reem, lsum)
      use diffuse_mod,   only : diffuse_absorber_leaf, reemission_field_valid, &
                                NABSORBER, ABSORBER_CHANNELS
      use ion_score_mod, only : previous_photoionizations
      use gas_state_mod, only : gas_nleaf
      real(kind=wp), intent(in) :: lem(:,:), xem(:,:), reem(:,:), lsum
      real(kind=wp), allocatable :: psi(:,:), en(:,:), nrec(:,:), nabs(:,:), fs(:,:)
      real(kind=wp) :: nnow(NABSORBER), b, d, sabs, srms, smax, stot, wl
      integer :: il, is, k, c
      logical :: valid
      if (.not. reemission_field_valid) return
      allocate(psi(NABSORBER,gas_nleaf), en(NABSORBER,gas_nleaf), &
               nrec(NABSORBER,gas_nleaf), nabs(NABSORBER,gas_nleaf), fs(NABSORBER,gas_nleaf))
      call diffuse_absorber_leaf(psi, en, nrec, nabs, fs)
      sabs = 0.0_wp;  srms = 0.0_wp;  smax = 0.0_wp;  stot = 0.0_wp;  wl = 0.0_wp
      do il = 1, gas_nleaf
         call previous_photoionizations(il, nnow, valid)
         do is = 1, NABSORBER
            if (psi(is,il) < 0.0_wp) cycle
            do k = 1, 2
               c = ABSORBER_CHANNELS(k,is)
               if (c <= 0) cycle
               b = fs(is,il)*lem(c,il)/max(nrec(is,il), nabs(is,il))
               d = reem(c,il) - b*nnow(is)
               sabs = sabs + abs(d);  stot = stot + d;  smax = max(smax, abs(d))
               if (lem(c,il) > 0.0_wp) then
                  wl = wl + lem(c,il);  srms = srms + d*d/lem(c,il)
               end if
            end do
         end do
      end do
      if (wl <= 0.0_wp) return
      write(6,'(a,es11.3,a,es11.3,a,es11.3,a,es11.3)') &
         ' REEMIT: count residual Re - b N_now: Sum|.|/Sum L = ', sabs/lsum, &
         ', L-weighted rms of ./L = ', sqrt(srms/wl), ', max|.|/Sum L = ', smax/lsum, &
         ', signed Sum/Sum L = ', stot/lsum
    end subroutine report_count_residuals

    !--- the distribution of psi_s over the leaves where absorber s re-emits
    !--- (rank 0): minimum, median, maximum and the fraction with psi > 1,
    !--- for H0, He0, He+
    subroutine report_absorber_psi()
      use diffuse_mod,   only : diffuse_absorber_leaf, reemission_field_valid, NABSORBER, &
                                ABSORBER_CHANNELS, diffuse_leaf_emission
      use gas_state_mod, only : gas_nleaf
      integer, parameter :: NCHAN = 8
      character(len=4), parameter :: ABS_NAME(NABSORBER) = ['H0  ', 'He0 ', 'He+ ']
      real(kind=wp), allocatable :: psi(:,:), en(:,:), v(:), nrec(:,:), nabs(:,:), fs(:,:), lem(:,:)
      real(kind=wp) :: moved, ls, ls_tot
      integer :: is, n, il, k
      if (.not. reemission_field_valid) then
         write(6,'(a)') '     no valid previous field: re-emission off, '// &
            'every channel launched explicitly'
         return
      end if
      allocate(psi(NABSORBER,gas_nleaf), en(NABSORBER,gas_nleaf), &
               nrec(NABSORBER,gas_nleaf), nabs(NABSORBER,gas_nleaf), fs(NABSORBER,gas_nleaf), &
               lem(NCHAN,gas_nleaf))
      call diffuse_absorber_leaf(psi, en, nrec, nabs, fs)
      call diffuse_leaf_emission(lem)
      do is = 1, NABSORBER
         !--- the luminosity the cap eps_cap moves to the explicit launch,
         !--- Sum_c L_c (1 - f) min(psi, 1) over the leaves where s re-emits
         moved = 0.0_wp;  ls_tot = 0.0_wp
         do il = 1, gas_nleaf
            ls = 0.0_wp
            do k = 1, 2
               if (ABSORBER_CHANNELS(k,is) > 0) ls = ls + lem(ABSORBER_CHANNELS(k,is),il)
            end do
            ls_tot = ls_tot + ls
            if (psi(is,il) >= 0.0_wp) moved = moved + ls*(1.0_wp - fs(is,il))*min(psi(is,il), 1.0_wp)
         end do
         n = count(psi(is,:) >= 0.0_wp)
         write(6,'(3a,2es12.4,a,es11.3,a,i0,a,es11.3,a,es11.3)') '     ', ABS_NAME(is), &
            ' Sum N_rec, Sum N_abs [1/s] =', sum(nrec(is,:)), sum(nabs(is,:)), &
            ';  Sum eps_re N_rec/max(psi,1) =', sum(en(is,:)*nrec(is,:)), &
            ';  leaves with f < 1: ', count(psi(is,:) >= 0.0_wp .and. fs(is,:) < 1.0_wp), &
            ';  Sum L_s (1-f) min(psi,1), moved to the explicit launch by the cap [erg/s] =', moved, ' of Sum L_s =', ls_tot
      end do
      do is = 1, NABSORBER
         n = count(psi(is,:) >= 0.0_wp)
         if (n == 0) then
            write(6,'(3a)') '     psi(', trim(ABS_NAME(is)), '): no leaf re-emits'
            cycle
         end if
         v = pack(psi(is,:), psi(is,:) >= 0.0_wp)
         write(6,'(3a,i0,a,3es11.3,a,f9.6)') '     psi(', trim(ABS_NAME(is)), ') over ', &
            n, ' leaves: min/median/max =', minval(v), median_of(v), maxval(v), &
            ', fraction psi > 1: ', real(count(v > 1.0_wp), wp)/real(n, wp)
      end do
    end subroutine report_absorber_psi

    !--- the median of v (v is reordered): Hoare selection of the middle
    !--- element (the upper one for an even count)
    real(kind=wp) function median_of(v) result(m)
      real(kind=wp), intent(inout) :: v(:)
      integer :: lo, hi, i, j, k
      real(kind=wp) :: pivot, tmp
      k = size(v)/2 + 1
      lo = 1;  hi = size(v)
      do while (lo < hi)
         pivot = v((lo + hi)/2)
         i = lo;  j = hi
         do while (i <= j)
            do while (v(i) < pivot);  i = i + 1;  end do
            do while (v(j) > pivot);  j = j - 1;  end do
            if (i <= j) then
               tmp = v(i);  v(i) = v(j);  v(j) = tmp
               i = i + 1;  j = j - 1
            end if
         end do
         if (k <= j) then
            hi = j
         else if (k >= i) then
            lo = i
         else
            exit
         end if
      end do
      m = v(k)
    end function median_of

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
