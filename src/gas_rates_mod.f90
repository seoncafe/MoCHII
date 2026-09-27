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
! ionizing-band tally jt_ion(inu, il) remains a diagnostic spectrum, written out
! with par%jnu_output as J_nu = jt_ion / (4 pi V_leaf dnu distance2cm^2)
! [erg/s/cm^2/Hz/sr]; without it (and without the SEDust emission, which also
! reads it) the tally is not kept at all.
!
! gas_rates_write writes '<base>_rates.<ext>': Gamma / heating arrays,
! J_nu(nnu, nleaf) with par%jnu_output, the energy grid, and leaf centers.
!---------------------------------------------------------------------------
  use define
  use octree_mod,   only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
  use jtally_mod,   only : jt_ion
  use ion_band_mod, only : ion_e, ion_de, ion_dnu, nnu_band
  implicit none
  private

  public :: gas_rates_compute, gas_rates_write
  public :: gamma_HI, gamma_HeI, gamma_HeII, heat_HI, heat_HeI, heat_HeII
  public :: heat_dust, g0_fuv
  public :: run_converged, run_iters, run_final_dx, run_final_dte
  public :: run_field_consistent
  public :: secondary_ionization

  real(kind=wp), allocatable :: gamma_HI(:), gamma_HeI(:), gamma_HeII(:)
  real(kind=wp), allocatable :: heat_HI(:),  heat_HeI(:),  heat_HeII(:)

  !--- secondary ionization by fast photoelectrons (par%use_sec_ion;
  !--- photoelectron_degradation_mod): the photoelectron energy deposition
  !--- rate above E_SECONDARY [erg/s per absorber particle], sec_edep(k, a,
  !--- leaf) for energy node k and absorber a = H I, He I, He II.  It depends
  !--- on the radiation field only; the partition at the gas state is applied
  !--- by secondary_ionization, inside the cell solve.  Metals are not
  !--- partitioned: they are trace absorbers, and their photoelectron energy is
  !--- all heat.
  real(kind=wp), allocatable :: sec_edep(:,:,:)
  !--- EUV grain heating [erg s^-1 cm^-3]: the ionizing-band part of the
  !--- grain heating integral that the dust SED band misses;
  !--- consumed by the SEDust stage later.
  real(kind=wp), allocatable :: heat_dust(:)
  !--- local FUV field in Habing units, G0 = int_{FUV} 4 pi J_nu dnu
  !--- / 1.6e-3 erg/s/cm^2 (nonzero only with par%add_fuv); drives the
  !--- grain photoelectric heating (par%grain_pe).
  real(kind=wp), allocatable :: g0_fuv(:)

  !--- convergence state of the gas iteration, set by main.f90 after the
  !--- loop and written to the rates-file header so an unconverged run is
  !--- self-documenting (the iteration can hit the cap without converging).
  logical       :: run_converged = .false.
  integer       :: run_iters     = 0
  real(kind=wp) :: run_final_dx  = 0.0_wp
  real(kind=wp) :: run_final_dte = 0.0_wp
  !--- .true. when the written rates, gas heating and dust heating were
  !--- rebuilt by transporting the WRITTEN gas state, so that state and the
  !--- radiation field belong to the same iterate.  .false. means they are
  !--- one gas solve apart (rates-only runs, par%gas_niter < 1).
  logical       :: run_field_consistent = .false.

contains

  !=========================================================================
  subroutine gas_rates_compute()
    implicit none
    integer  :: nleaf

    nleaf = amr_grid%nleaf
    if (allocated(gamma_HI)) then
       if (size(gamma_HI) /= nleaf) &
          deallocate(gamma_HI, gamma_HeI, gamma_HeII, &
                     heat_HI, heat_HeI, heat_HeII, heat_dust, &
                     g0_fuv)                                    ! re-refinement
    end if
    if (.not. allocated(gamma_HI)) &
       allocate(gamma_HI(nleaf), gamma_HeI(nleaf), gamma_HeII(nleaf), &
                heat_HI(nleaf),  heat_HeI(nleaf),  heat_HeII(nleaf), &
                heat_dust(nleaf), g0_fuv(nleaf))
    gamma_HI = 0.0_wp;  gamma_HeI = 0.0_wp;  gamma_HeII = 0.0_wp
    heat_HI  = 0.0_wp;  heat_HeI  = 0.0_wp;  heat_HeII  = 0.0_wp
    heat_dust = 0.0_wp; g0_fuv = 0.0_wp

    !--- The already reduced packet-energy path estimators are the rates.
    !--- jt_ion stays populated only as a diagnostic spectrum and must not feed
    !--- solver rates.
    block
      use ion_score_mod, only : ion_score_apply_hhe
      call ion_score_apply_hhe(gamma_HI, gamma_HeI, gamma_HeII, &
                               heat_HI, heat_HeI, heat_HeII)
    end block
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
         if (allocated(sec_edep)) then
            if (size(sec_edep, 3) /= nleaf) deallocate(sec_edep)
         end if
         if (.not. allocated(sec_edep)) allocate(sec_edep(N_EDEP, 3, nleaf))
         call ion_score_apply_sec(sec_edep)
       end block
    end if
    if (par%add_fuv) then
       block
         use ion_score_mod, only : ion_score_apply_g0
         call ion_score_apply_g0(g0_fuv)
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
  subroutine gas_rates_write()
    use iofile_mod
    use utility, only : get_base_name, fatal_error
    implicit none
    type(io_file_type) :: file
    character(len=PATH_LEN) :: filename
    real(kind=wp), allocatable :: leafxyz(:,:), leafsize(:), jnu(:,:)
    real(kind=wp) :: energy_emitted, energy_absorbed, energy_escaped
    integer :: status, il, inu, ic, nleaf

    if (mpar%p_rank /= 0) return
    nleaf = amr_grid%nleaf

    status = 0
    filename = trim(get_base_name(par%out_file))//'_rates'// &
               trim(io_file_extension(par%file_format))
    call io_open_new(file, trim(filename), status)

    call io_append_image(file, gamma_HI, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Gamma_HI','H I photoionization rate [s^-1]',status)
    call io_put_keyword(file,'NNU_ION', par%nnu_ion, 'ionizing energy-bin count', status)
    call io_put_keyword(file,'NNU_BAND', nnu_band, 'total band bins (ionizing+FUV)', status)
    call io_put_keyword(file,'EION_MIN',par%eion_min,'band lower edge [eV]',    status)
    call io_put_keyword(file,'EION_MAX',par%eion_max,'band upper edge [eV]',    status)
    call io_put_keyword(file,'TOT_LUM', par%luminosity,'band luminosity [erg/s]',status)
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
    call io_put_keyword(file,'CONVERGD',run_converged,'gas iteration converged',status)
    call io_put_keyword(file,'NITERDON',run_iters,   'gas iterations run',      status)
    call io_put_keyword(file,'FINALDX', run_final_dx, 'final max|delta x_HII|',  status)
    call io_put_keyword(file,'FINALDTE',run_final_dte,'final max|delta Te|/Te',  status)
    call io_put_keyword(file,'FLDCONS', run_field_consistent, &
       'rates rebuilt by transporting the written state', status)
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
      call ion_score_energy_totals(energy_emitted, energy_absorbed, energy_escaped)
    end block
    call io_put_keyword(file,'L_EMIT',energy_emitted, &
       'transported ionizing packet luminosity',status)
    call io_put_keyword(file,'L_ABS',energy_absorbed, &
       'absorbed ionizing packet luminosity',status)
    call io_put_keyword(file,'L_ESC',energy_escaped, &
       'escaped ionizing packet luminosity',status)
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
             jnu(inu,il) = jt_ion(inu,il) / (fourpi * (2.0_wp*leaf_half(il))**3 &
                           * ion_dnu(inu) * par%distance2cm**2)
          end do
       end do
       call io_append_image(file, jnu, status, bitpix=-64)
       call io_put_keyword(file,'EXTNAME','J_nu','J(nu,leaf) [erg/s/cm^2/Hz/sr]',status)
       deallocate(jnu)
    end if

    !--- ionization structure (fixed initialization when gas_niter = 0).
    block
      use gas_state_mod, only : gas_xHI, gas_xHII, gas_xHeI, gas_xHeII, &
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
      tmp = gas_ne(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','n_e','electron density [cm^-3]',status)
      tmp = gas_Te(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','T_e','electron temperature [K]',status)
      deallocate(tmp)
    end block

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

    !--- converged metal stage fractions.
    if (par%use_metals) then
       block
         use species_mod, only : species_write
         call species_write(file, status)
       end block
    end if

    call io_close(file, status)
    if (status /= 0) call fatal_error('failed to write rates file: '//trim(filename))
    write(*,'(2a)') ' ION: rates written to: ', trim(filename)
  end subroutine gas_rates_write

end module gas_rates_mod
