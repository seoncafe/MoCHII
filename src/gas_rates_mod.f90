module gas_rates_mod
!---------------------------------------------------------------------------
! MoCHII: photoionization and photoheating rate integrals.
!
! Each transported packet scores its path against the cross sections at its own
! photon energy (ion_score_mod), so the rate integrals
!     Gamma_i = Int 4 pi J_nu sigma_i(E) dnu / (h nu)          [s^-1]
!     H_i     = Int 4 pi J_nu sigma_i(E) dnu (1 - E_th,i/E)    [erg/s]
! per particle of species i (i = HI, HeI, HeII) are accumulated continuously
! rather than reconstructed from a binned tally.  This module normalizes those
! reduced path sums by the leaf volume and hands them to the solver.  The
! ionizing-band tally jt_ion_sum(inu, il) supplies diagnostics and dust physics,
! and is written out
! with par%jnu_output as J_nu = jt_ion_sum / (4 pi V_leaf dnu distance2cm^2)
! [erg/s/cm^2/Hz/sr]; without it (and without the SEDust emission and the grain
! charge distribution of par%grain_pe, which also read it) the tally is not
! kept at all.
!
! gas_rates_write writes '<base>_rates.<ext>': Gamma / heating arrays,
! J_nu(nnu, nleaf) with par%jnu_output, the energy grid, and leaf centers.
!---------------------------------------------------------------------------
  use define
  use mpi
  use memory_mod,   only : shared_owner, allocate_shared, free_shared, publish_shared, shared_is_live, &
                          shared_host_rank, shared_host_size
  use octree_mod,   only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
  use jtally_mod,   only : jt_ion_sum
  use ion_band_mod, only : ion_e, ion_de, ion_dnu, nnu_band, source_luminosity
  implicit none
  private

  public :: gas_rates_release

  public :: gas_rates_compute, gas_rates_write
  public :: gamma_HI, gamma_HeI, gamma_HeII, heat_HI, heat_HeI, heat_HeII
  public :: heat_dust
  public :: run_converged, run_iters, run_final_dx, run_final_dte, run_stop_reason
  public :: run_state_solved_in_field
  public :: secondary_ionization

  !--- H/He photoionization rates and photoheating per leaf: the same on
  !--- every rank, so one node-shared copy (2026-10-12); the node's ranks
  !--- fill disjoint leaves from the reduced path estimators
  real(kind=wp), pointer, contiguous :: gamma_HI(:) => null(), &
       gamma_HeI(:) => null(), gamma_HeII(:) => null()
  real(kind=wp), pointer, contiguous :: heat_HI(:) => null(), &
       heat_HeI(:) => null(), heat_HeII(:) => null()
  type(shared_owner) :: rate_owner(6)

  !--- secondary ionization by fast photoelectrons (par%use_sec_ion;
  !--- photoelectron_degradation_mod): the photoelectron energy deposition
  !--- rate above E_SECONDARY [erg/s per absorber particle], sec_edep(k, a,
  !--- leaf) for energy node k and absorber a = H I, He I, He II.  It depends
  !--- on the radiation field only; the partition at the gas state is applied
  !--- by secondary_ionization, inside the cell solve.  Metals are not
  !--- partitioned: they are trace absorbers, and their photoelectron energy is
  !--- all heat.
  real(kind=wp), pointer, contiguous :: sec_edep(:,:,:) => null()   ! node-shared
  type(shared_owner) :: sec_owner(1)
  !--- EUV grain heating [erg s^-1 cm^-3]: the ionizing-band part of the
  !--- grain heating integral that the dust SED band misses;
  !--- consumed by the SEDust stage later.
  real(kind=wp), allocatable :: heat_dust(:)

  !--- convergence state of the gas iteration, set by main.f90 after the
  !--- loop (and for each iterate written under par%rates_write_interval)
  !--- and written to the rates-file header so an unconverged run is
  !--- self-documenting (the iteration can hit the cap without converging).
  logical       :: run_converged = .false.
  character(len=96) :: run_stop_reason = 'not converged (input state; no gas iteration)'
  integer       :: run_iters     = 0
  real(kind=wp) :: run_final_dx  = 0.0_wp
  real(kind=wp) :: run_final_dte = 0.0_wp
  !--- .true. when the written gas state was solved in the written field
  !--- (main.f90, final pass and final solve), so H, He, T_e, n_e and the
  !--- metal stages are the solution of the written rates; the field is the
  !--- transport of the state one solve earlier, run_final_dx/_dte apart.
  !--- .false. for rates-only runs (par%gas_niter < 1): the written state
  !--- is the input state.
  logical       :: run_state_solved_in_field = .false.

contains

  !=========================================================================
  subroutine gas_rates_compute()
    implicit none
    integer  :: nleaf

    call publish_shared(rate_owner)
    call publish_shared(sec_owner)
    nleaf = amr_grid%nleaf
    if (associated(gamma_HI)) then
       if (size(gamma_HI) /= nleaf) then                     ! re-refinement
          call free_shared(rate_owner(1), gamma_HI);  call free_shared(rate_owner(2), gamma_HeI)
          call free_shared(rate_owner(3), gamma_HeII)
          call free_shared(rate_owner(4), heat_HI);   call free_shared(rate_owner(5), heat_HeI)
          call free_shared(rate_owner(6), heat_HeII)
       end if
    end if
    if (.not. associated(gamma_HI)) then
       call allocate_shared(rate_owner(1), gamma_HI, [nleaf])
       call allocate_shared(rate_owner(2), gamma_HeI, [nleaf])
       call allocate_shared(rate_owner(3), gamma_HeII, [nleaf])
       call allocate_shared(rate_owner(4), heat_HI, [nleaf])
       call allocate_shared(rate_owner(5), heat_HeI, [nleaf])
       call allocate_shared(rate_owner(6), heat_HeII, [nleaf])
    end if
    !--- heat_dust stays private: the dust emission transport adds the
    !--- reabsorbed infrared to it on every rank
    if (allocated(heat_dust)) then
       if (size(heat_dust) /= nleaf) deallocate(heat_dust)
    end if
    if (.not. allocated(heat_dust)) allocate(heat_dust(nleaf))
    heat_dust = 0.0_wp

    !--- The already reduced packet-energy path estimators are the rates.
    !--- jt_ion_sum supplies spectral diagnostics and dust physics; H/He gas
    !--- photoionization and photoheating use the continuous estimators.
    block
      use ion_score_mod, only : ion_score_apply_hhe
      call ion_score_apply_hhe(gamma_HI, gamma_HeI, gamma_HeII, &
                               heat_HI, heat_HeI, heat_HeII)
    end block
    call publish_shared(rate_owner)
    if (par%ion_add_dust) then
       block
         use ion_score_mod, only : ion_score_apply_dust_heat
         call ion_score_apply_dust_heat(heat_dust)
       end block
    end if
    !--- photoelectron energy deposition above E_SECONDARY, by energy node and
    !--- absorber; the partition at the gas state is applied in the cell solve
    !--- (secondary_ionization).
    if (par%use_sec_ion) then
       block
         use ion_score_mod, only : ion_score_apply_sec
         use photoelectron_degradation_mod, only : N_EDEP
         if (associated(sec_edep)) then
            if (size(sec_edep, 3) /= nleaf) call free_shared(sec_owner(1), sec_edep)
         end if
         if (.not. associated(sec_edep)) then
            call allocate_shared(sec_owner(1), sec_edep, [N_EDEP, 3, nleaf])
         end if
         call ion_score_apply_sec(sec_edep)
         call publish_shared(sec_owner)
       end block
    end if
  end subroutine gas_rates_compute

  !=========================================================================
  ! Secondary ionization at the gas state of leaf il (par%use_sec_ion):
  !   x = n_e/(n_H + n_He) and the neutral densities n_HI, n_HeI of the
  !   CURRENT iterate of the cell solve, n_HeII the He II absorber density.
  ! Returns the secondary photoionization rates dgHI, dgHeI [s^-1 per target
  ! particle], added to Gamma_HI, Gamma_HeI in the balance, and the heat each
  ! absorber particle does NOT deposit, dhHI, dhHeI, dhHeII [erg/s], taken
  ! off heat_* in the thermal sum:
  !   dh_a  = Sum_k [1 - f_heat(k)] edep(k,a),
  !   dg_t  = Sum_k c_t(k) Sum_a n_a edep(k,a) / E_th,t .
  ! Energy is conserved node by node: n_HI c_HI + n_HeI c_HeI = f_ion(k)
  ! <= 1 - f_heat(k), the rest being excitation that escapes as lines.
  !=========================================================================
  subroutine secondary_ionization(il, x, nHI, nHeI, nHeII, dgHI, dgHeI, &
                                  dhHI, dhHeI, dhHeII)
    use photoelectron_degradation_mod, only : N_EDEP, photoelectron_partition
    implicit none
    integer,       intent(in)  :: il
    real(kind=wp), intent(in)  :: x, nHI, nHeI, nHeII
    real(kind=wp), intent(out) :: dgHI, dgHeI, dhHI, dhHeI, dhHeII
    real(kind=wp) :: fh(N_EDEP), cHI(N_EDEP), cHeI(N_EDEP), dep(N_EDEP)

    call photoelectron_partition(x, nHI, nHeI, fh, cHI, cHeI)
    dep  = nHI*sec_edep(:,1,il) + nHeI*sec_edep(:,2,il) + nHeII*sec_edep(:,3,il)
    dgHI   = sum(cHI *dep)/(eth_HI *ev2erg)
    dgHeI  = sum(cHeI*dep)/(eth_HeI*ev2erg)
    dhHI   = sum((1.0_wp - fh)*sec_edep(:,1,il))
    dhHeI  = sum((1.0_wp - fh)*sec_edep(:,2,il))
    dhHeII = sum((1.0_wp - fh)*sec_edep(:,3,il))
  end subroutine secondary_ionization

  !=========================================================================
  !--- iterate = .true. (main.f90, par%rates_write_interval) writes the
  !--- state of the gas iteration in progress to <base>_rates_iterate: it is
  !--- written to a temporary name and renamed over the previous iterate
  !--- only once complete.  The header then records CONVERGD F, STATESOL F
  !--- and the iteration count, so restart_file can resume from it.
  subroutine gas_rates_write(iterate)
    use iofile_mod
    use utility, only : get_base_name, abort_run
    implicit none
    logical, intent(in), optional :: iterate
    type(io_file_type) :: file
    character(len=PATH_LEN) :: filename, partname
    logical :: write_iterate
    integer :: cmdstat, exitstat
    real(kind=wp), allocatable :: leafxyz(:,:), leafsize(:), jnu(:,:)
    real(kind=wp) :: energy_emitted, energy_absorbed, energy_escaped, energy_lost
    integer :: status, il, inu, ic, nleaf, writer_ierr

    ! All ranks call this writer. Hold the other ranks until rank 0 has
    ! copied the shared frozen emission; the next diffuse_build overwrites
    ! it before its first collective. Keep default-mode execution unchanged.
    if (mpar%p_rank /= 0) then
       if (trim(par%diffuse_launch) == 'absorption_reemission') &
          call MPI_BARRIER(MPI_COMM_WORLD, writer_ierr)
       return
    end if
    nleaf = amr_grid%nleaf

    write_iterate = .false.
    if (present(iterate)) write_iterate = iterate
    status = 0
    if (write_iterate) then
       filename = trim(get_base_name(par%out_file))//'_rates_iterate'// &
                  trim(io_file_extension(par%file_format))
       partname = trim(get_base_name(par%out_file))//'_rates_iterate_part'// &
                  trim(io_file_extension(par%file_format))
    else
       filename = trim(get_base_name(par%out_file))//'_rates'// &
                  trim(io_file_extension(par%file_format))
       partname = filename
    end if
    call io_open_new(file, trim(partname), status)

    call io_append_image(file, gamma_HI, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Gamma_HI','H I photoionization rate [s^-1]',status)
    call io_put_keyword(file,'NNU_ION', par%nnu_ion, 'ionizing energy-bin count', status)
    call io_put_keyword(file,'NNU_BAND', nnu_band, 'total band bins (ionizing+FUV)', status)
    call io_put_keyword(file,'EION_MIN',par%eion_min,'band lower edge [eV]',    status)
    call io_put_keyword(file,'EION_MAX',par%eion_max,'band upper edge [eV]',    status)
    call io_put_keyword(file,'TOT_LUM', par%luminosity,'source luminosity scale [erg/s]',status)
    call io_put_keyword(file,'LUMREF',trim(par%luminosity_reference),'source normalization reference',status)
    call io_put_keyword(file,'DIST_CM', par%distance2cm,'distance unit (cm)',   status)
    call io_put_keyword(file,'nphotons',par%no_photons, 'number of photons',    status)
    if (par%add_fuv) then
       call io_put_keyword(file,'EFUV_MIN',par%efuv_min,'FUV extension lower edge [eV]', status)
       call io_put_keyword(file,'FUVSAMP', trim(par%fuv_sampling), &
          'FUV packets: separate stratum or one CDF', status)
       if (trim(par%fuv_sampling) == 'separate') &
          call io_put_keyword(file,'NPHOTFUV',par%no_photons_fuv, &
             'FUV-stratum packets per pass', status)
    end if
    call io_put_keyword(file,'CONVERGD',run_converged,trim(run_stop_reason),status)
    call io_put_keyword(file,'NITERDON',run_iters,   'gas iterations run',      status)
    call io_put_keyword(file,'FINALDX', run_final_dx, 'last-solve max|delta x_HII|', status)
    call io_put_keyword(file,'FINALDTE',run_final_dte,'last-solve max|delta Te|/Te', status)
    call io_put_keyword(file,'STATESOL', run_state_solved_in_field, &
       'written state solved in the written field', status)
    !--- provenance of the spectral model behind these rates: packet energies
    !--- were sampled continuously from each source's spectrum, so the spectral
    !--- bins below carry the diagnostic spectrum only.
    call io_put_keyword(file,'IONEMODE','continuous', &
       'ionizing packet energy sampling', status)
    call io_put_keyword(file,'CDFRTOL', par%source_cdf_tol, &
       'source sampler relative tolerance', status)
    call io_put_keyword(file,'NUBINUSE','diagnostic_only', &
       'role of ionizing spectral bins', status)
    call io_put_keyword(file,'ENRGSAMP',trim(par%launch_sequence), &
       'energy/source sampling sequence', status)
    block
      use ion_score_mod, only : ion_score_energy_totals
      call ion_score_energy_totals(energy_emitted, energy_absorbed, &
                                   energy_escaped, energy_lost)
    end block
    !--- the physical totals of the pass (ion_score_mod): L_ESC is the direct
    !--- transmission plus the dust-scattered flights that leave the grid,
    !--- L_LOST the weight of walks that ended without leaving the nebula
    !--- (normally 0), and L_ABS = L_EMIT - L_ESC - L_LOST.
    call io_put_keyword(file,'L_EMIT',energy_emitted, &
       'transported ionizing packet luminosity',status)
    call io_put_keyword(file,'L_ABS',energy_absorbed, &
       'absorbed: L_EMIT - L_ESC - L_LOST',status)
    call io_put_keyword(file,'L_ESC',energy_escaped, &
       'escaped: direct + scattered out of grid',status)
    call io_put_keyword(file,'L_LOST',energy_lost, &
       'lost in walks (step cap, slab re-entry)',status)
    if (trim(par%diffuse_launch) == 'absorption_reemission') then
       block
         use ion_score_mod, only : ion_score_reemission_totals, ion_score_interaction_energy_totals
         use diffuse_mod, only : diffuse_nphot, diffuse_nominal_packet_lum, reemission_field_valid
         real(kind=wp) :: re, leafem, netabs, src, launched, ext, trans, dep, roul, cesc, clost, gas, dust
         call ion_score_reemission_totals(re, leafem, netabs)
         call ion_score_interaction_energy_totals(src, launched, ext, trans, dep, roul, cesc, clost, gas, dust)
         call io_put_keyword(file,'DIFLAUNC',trim(par%diffuse_launch), &
              'diffuse emission position rule',status)
         call io_put_keyword(file,'L_REEMIT',re, &
              'gas re-emission [erg/s]',status)
         call io_put_keyword(file,'L_LEAFEM',leafem, &
              'explicit leaf emission [erg/s]',status)
         call io_put_keyword(file,'L_NETABS',netabs, &
              'L_ABS minus L_REEMIT [erg/s]',status)
         call io_put_keyword(file,'L_SRCEMT',src, &
              'launched sources [erg/s]',status)
         call io_put_keyword(file,'L_DIREMT',launched, &
              'sources plus explicit emission [erg/s]',status)
         call io_put_keyword(file,'L_EXTING',ext, &
              'extinguished direct light [erg/s]',status)
         call io_put_keyword(file,'L_TRANS',trans, &
              'direct transmitted light [erg/s]',status)
         call io_put_keyword(file,'L_DEPOS',dep, &
              'signed net event deposition [erg/s]',status)
         call io_put_keyword(file,'L_ROULET',roul, &
              'signed roulette correction [erg/s]',status)
         call io_put_keyword(file,'L_CESC',cesc, &
              'continued escaping light [erg/s]',status)
         call io_put_keyword(file,'L_CLOST',clost, &
              'continued lost light [erg/s]',status)
         call io_put_keyword(file,'L_GASEVT',gas, &
              'gas event absorption [erg/s]',status)
         call io_put_keyword(file,'L_DUSTEVT',dust, &
              'dust event absorption [erg/s]',status)
         call io_put_keyword(file,'NLEAFEM',diffuse_nphot,'explicit leaf packet count',status)
         call io_put_keyword(file,'LPACKREF',diffuse_nominal_packet_lum,'reference packet luminosity [erg/s]',status)
         call io_put_keyword(file,'REVALID',reemission_field_valid,'previous field used for re-emission',status)
       end block
    end if
    if (trim(par%launch_sequence) == 'sobol') &
       call io_put_keyword(file,'QMCSEED',par%qmc_seed,'Sobol scramble seed',status)
    call io_append_image(file, gamma_HeI, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Gamma_HeI','He I photoionization rate [s^-1]',status)
    call io_append_image(file, gamma_HeII, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Gamma_HeII','He II photoionization rate [s^-1]',status)
    call io_append_image(file, heat_HI, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Heat_HI','H I photoheating [erg/s per HI]',status)
    call io_append_image(file, heat_HeI, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Heat_HeI','He I photoheating [erg/s per HeI]',status)
    call io_append_image(file, heat_HeII, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Heat_HeII','He II photoheating [erg/s per HeII]',status)
    if (par%ion_add_dust) then
       call io_append_image(file, heat_dust, status, bitpix=-64)
       !--- with par%dust_emis_transport this carries the reabsorbed thermal
       !--- infrared as well as the transported EUV/FUV band.
       call io_put_keyword(file,'EXTNAME','Heat_dust','grain heating [erg/s/cm^3]',status)
       block
         use dust_temp_mod, only : t_dust
         if (allocated(t_dust)) then
            call io_append_image(file, t_dust, status, bitpix=-64)
            call io_put_keyword(file,'EXTNAME','T_dust', &
                 'equilibrium mixture dust temperature [K]',status)
         end if
       end block
    end if

    !--- J_nu(nnu, nleaf) [erg/s/cm^2/Hz/sr], with par%jnu_output only
    !--- (the tally is not kept otherwise; jtally_mod::jtally_ion_needed)
    if (par%jnu_output) then
       allocate(jnu(nnu_band, nleaf))
       do il = 1, nleaf
          do inu = 1, nnu_band
             jnu(inu,il) = jt_ion_sum(inu,il) / (fourpi * (2.0_wp*leaf_half(il))**3 &
                           * ion_dnu(inu) * par%distance2cm**2)
          end do
       end do
       call io_append_image(file, jnu, status, bitpix=-64)
       call io_put_keyword(file,'EXTNAME','J_nu','J(nu,leaf) [erg/s/cm^2/Hz/sr]',status)
       deallocate(jnu)
    end if

    !--- ionization structure (fixed initialization when gas_niter = 0).
    block
      use gas_state_mod, only : gas_nH, gas_xHI, gas_xHII, gas_xHeI, gas_xHeII, &
                                gas_xHeIII, gas_ne, gas_Te
      real(kind=wp), allocatable :: tmp(:)
      allocate(tmp(nleaf))
      tmp = gas_xHI(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','x_HI','n_HI/n_H per leaf',status)
      !--- x_HII beside x_HI, to its own relative precision (beyond the
      !--- ionization front 1 - x_HI does not represent it)
      tmp = gas_xHII(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','x_HII','n_HII/n_H per leaf',status)
      tmp = gas_xHeI(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','x_HeI','n_HeI/n_He per leaf',status)
      tmp = gas_xHeII(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','x_HeII','n_HeII/n_He per leaf',status)
      !--- x_HeIII to its own relative precision (1 - x_HeI - x_HeII is
      !--- rounding where He2+ is scarce)
      tmp = gas_xHeIII(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','x_HeIII','n_HeIII/n_He per leaf',status)
      !--- the hydrogen density, 0 in vacuum (the same block, name and
      !--- meaning as n_H in <base>_emis): restart_file tells vacuum from
      !--- gas with it (previous_state_mod)
      tmp = gas_nH(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','n_H','H density [cm^-3]',status)
      tmp = gas_ne(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','n_e','electron density [cm^-3]',status)
      tmp = gas_Te(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','T_e','electron temperature [K]',status)
      !--- He I 2^3S photoionization rate of the last transport pass: the
      !--- 2^3S population behind the He I diffuse channel of the next pass
      !--- reads it before that pass's own transport, so restart_file
      !--- restores it (previous_state_mod)
      if (par%hei_metastable_coupling) then
         block
           use hei_metastable_mod, only : hei_phi3_rate
           do il = 1, nleaf
              tmp(il) = hei_phi3_rate(il)
           end do
         end block
         call io_append_image(file, tmp, status, bitpix=-64)
         call io_put_keyword(file,'EXTNAME','Phi_HeI_23S', &
              'He I 2^3S photoionization rate [s^-1]',status)
      end if
      deallocate(tmp)
    end block

    !--- absorption ledger by leaf (par%absorption_ledger): photons/s absorbed
    !--- by H0, He0, He+ in each leaf from each source (columns: 0 stellar
    !--- band, 1-8 diffuse channels, 9 stellar band below eion_min)
    if (par%absorption_ledger) then
       block
         use ion_score_mod, only : ion_score_ledger_leaf, LEDGER_NSRC
         real(kind=wp), allocatable :: led(:,:)
         character(len=10), parameter :: lname(3) = ['Ledger_H0 ', 'Ledger_He0', 'Ledger_He+']
         integer :: ia, is
         allocate(led(nleaf, 0:LEDGER_NSRC-1))
         do ia = 1, 3
            do is = 0, LEDGER_NSRC-1
               call ion_score_ledger_leaf(ia, is, led(:,is))
            end do
            call io_append_image(file, led, status, bitpix=-64)
            call io_put_keyword(file,'EXTNAME',trim(lname(ia)), &
                 'photons/s absorbed, by leaf and source',status)
         end do
         deallocate(led)
       end block
    end if

    call io_append_image(file, ion_e, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','E_bin','bin centers [eV]',status)
    call io_append_image(file, ion_de, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','dE_bin','bin widths [eV]',status)

    allocate(leafxyz(nleaf,3))
    do il = 1, nleaf
       leafxyz(il,1) = leaf_cx(il)
       leafxyz(il,2) = leaf_cy(il)
       leafxyz(il,3) = leaf_cz(il)
    end do
    call io_append_image(file, leafxyz, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','LeafXYZ','leaf center x,y,z (code units)',status)
    deallocate(leafxyz)

    allocate(leafsize(nleaf))
    do il = 1, nleaf
       leafsize(il) = 2.0_wp * leaf_half(il)
    end do
    call io_append_image(file, leafsize, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','LeafSize','leaf full width (code units)',status)
    deallocate(leafsize)

    if (allocated(source_luminosity)) then
       call io_append_image(file, source_luminosity, status, bitpix=-64)
       call io_put_keyword(file,'EXTNAME','SourceLuminosity','internal source normalization',status)
       call io_put_keyword(file,'LUMREF',trim(par%luminosity_reference),'source normalization reference',status)
       call io_put_keyword(file,'COL1','Emin [eV]','full spectrum lower edge',status)
       call io_put_keyword(file,'COL2','Emax [eV]','-1 denotes infinity for Planck',status)
       call io_put_keyword(file,'COL3','Lfull [erg/s]','full source luminosity',status)
       call io_put_keyword(file,'COL4','Lion [erg/s]','transported ionizing luminosity',status)
       call io_put_keyword(file,'COL5','LFUV [erg/s]','transported FUV luminosity',status)
       call io_put_keyword(file,'COL6','Loutside [erg/s]','luminosity outside transport band',status)
       call io_put_keyword(file,'COL7','Ltransport [erg/s]','total transported source luminosity',status)
    end if

    !--- the stopping-test histories (fixed_point_distance_mod), oldest
    !--- iteration first, columns d and s; par%restart_iteration_history
    !--- reads them back
    block
      use fixed_point_distance_mod, only : iteration_history, hist_dx, &
                                           hist_dhe1, hist_dhe3
      call write_history(hist_dx,   'ConvHist_x_HII')
      call write_history(hist_dhe1, 'ConvHist_HeI_HeII')
      call write_history(hist_dhe3, 'ConvHist_HeII_HeIII')
    end block

    !--- converged metal stage fractions.
    if (par%use_metals) then
       block
         use species_mod, only : species_write
         call species_write(file, status)
       end block
    end if

    if (trim(par%diffuse_launch) == 'absorption_reemission') then
       block
         use ion_score_mod, only : ion_score_absorbed_power_leaf, ion_score_reemission_leaf, &
                                  ion_score_ledger_channels, LEDGER_NSRC, LEDGER_NABS
         use diffuse_mod, only : diffuse_leaf_emission, diffuse_explicit_leaf_emission, diffuse_absorber_leaf
         real(kind=wp), allocatable :: ch(:,:), re(:,:), event(:), psi(:,:), q(:,:), nr(:,:), na(:,:), f(:,:)
         real(kind=wp) :: le(0:LEDGER_NSRC-1), la(LEDGER_NABS,0:LEDGER_NSRC-1), &
                         ls(0:LEDGER_NSRC-1), ll(0:LEDGER_NSRC-1), ledger(LEDGER_NSRC,4)
         allocate(ch(8,nleaf), re(8,nleaf), event(nleaf), psi(3,nleaf), q(3,nleaf), &
                  nr(3,nleaf), na(3,nleaf), f(3,nleaf))
         call ion_score_absorbed_power_leaf(event)
         call io_append_image(file,event,status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','GasAbsorbedPower','measured path absorption [erg/s]',status)
         call diffuse_leaf_emission(ch)
         call io_append_image(file,transpose(ch),status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','DiffuseStateLuminosity','frozen L_c [erg/s], leaf then channel',status)
         call diffuse_explicit_leaf_emission(ch)
         call io_append_image(file,transpose(ch),status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','DiffuseExplicitLuminosity','explicit X_c [erg/s], leaf then channel',status)
         call ion_score_reemission_leaf(re,event)
         call io_append_image(file,transpose(re),status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','ReemittedLuminosity','measured Re_c [erg/s], leaf then channel',status)
         call io_append_image(file,event,status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','GasEventAbsorption','gas event absorption [erg/s]',status)
         call diffuse_absorber_leaf(psi,q,nr,na,f)
         call io_append_image(file,transpose(nr),status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','RecombinationRate','N_rec [1/s]; H0 He0 He+',status)
         call io_append_image(file,transpose(na),status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','PhotoionizationRatePrev','raw previous N_abs [1/s]; H0 He0 He+',status)
         call io_append_image(file,transpose(psi),status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','ReemissionPsi','N_abs/N_rec; -1 inactive; H0 He0 He+',status)
         call io_append_image(file,transpose(q),status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','ReemissionEnergy','eps_re/max(psi,1) [erg]; H0 He0 He+',status)
         call io_append_image(file,transpose(f),status,bitpix=-64)
         call io_put_keyword(file,'EXTNAME','ReemissionFraction','eps_re/eps; H0 He0 He+',status)
         if (par%absorption_ledger) then
            call ion_score_ledger_channels(le,la,ls,ll)
            ledger(:,1) = le; ledger(:,2) = sum(la,dim=1); ledger(:,3) = ls; ledger(:,4) = ll
            call io_append_image(file,ledger,status,bitpix=-64)
            call io_put_keyword(file,'EXTNAME','LedgerTotals','photons/s: emitted absorbed escaped lost',status)
         end if
       end block
    end if

    call io_close(file, status)
    if (status /= 0) call abort_run('failed to write rates file: '//trim(partname))
    if (write_iterate) then
       call execute_command_line('/bin/mv -f "'//trim(partname)//'" "'// &
            trim(filename)//'"', exitstat=exitstat, cmdstat=cmdstat)
       if (cmdstat /= 0 .or. exitstat /= 0) call abort_run( &
            'failed to rename '//trim(partname)//' to '//trim(filename))
       write(*,'(2a,i0)') ' ION: iterate written to: ', trim(filename)//', iteration ', run_iters
    else
       write(*,'(2a)') ' ION: rates written to: ', trim(filename)
    end if
    if (trim(par%diffuse_launch) == 'absorption_reemission') &
       call MPI_BARRIER(MPI_COMM_WORLD, writer_ierr)
  contains

    subroutine write_history(h, name)
      use fixed_point_distance_mod, only : iteration_history
      type(iteration_history), intent(in) :: h
      character(len=*),        intent(in) :: name
      real(kind=wp), allocatable :: tab(:,:)
      if (h%n < 1) return
      allocate(tab(h%n, 2))
      tab(:,1) = h%d(1:h%n);  tab(:,2) = h%s(1:h%n)
      call io_append_image(file, tab, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME',name, &
           'stopping-test history: step d, signed step s',status)
      deallocate(tab)
    end subroutine write_history

  end subroutine gas_rates_write

  subroutine gas_rates_release()
    call publish_shared(rate_owner)
    call publish_shared(sec_owner)
    if (shared_is_live(rate_owner(1))) call free_shared(rate_owner(1), gamma_HI)
    if (shared_is_live(rate_owner(2))) call free_shared(rate_owner(2), gamma_HeI)
    if (shared_is_live(rate_owner(3))) call free_shared(rate_owner(3), gamma_HeII)
    if (shared_is_live(rate_owner(4))) call free_shared(rate_owner(4), heat_HI)
    if (shared_is_live(rate_owner(5))) call free_shared(rate_owner(5), heat_HeI)
    if (shared_is_live(rate_owner(6))) call free_shared(rate_owner(6), heat_HeII)
    if (shared_is_live(sec_owner(1))) call free_shared(sec_owner(1), sec_edep)
  end subroutine gas_rates_release

end module gas_rates_mod
