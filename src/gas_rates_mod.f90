module gas_rates_mod
!---------------------------------------------------------------------------
! MoCHII: photoionization and photoheating rate integrals.
!
! From the reduced ionizing-band tally jt_ion(inu, il) = Sum(Lpacket*wgt*dl)
! [erg/s * code length], the mean intensity per leaf is
!     J_nu = jt_ion / (4 pi V_leaf dnu distance2cm^2)  [erg/s/cm^2/Hz/sr]
! and the rate integrals follow as bin sums:
!     Gamma_i = Sum_nu 4 pi J_nu sigma_i(E) dnu / (h nu)          [s^-1]
!     H_i     = Sum_nu 4 pi J_nu sigma_i(E) dnu (1 - E_th,i/E)    [erg/s]
! both per particle of species i (i = HI, HeI, HeII).
!
! gas_rates_write writes '<base>_rates.<ext>': Gamma / heating arrays,
! J_ion(nnu, nleaf), the energy grid, and leaf centers.
!---------------------------------------------------------------------------
  use define
  use octree_mod,   only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
  use jtally_mod,   only : jt_ion
  use ion_band_mod, only : ion_e, ion_de, ion_nu, ion_dnu, nnu_band
  use photo_xsec,   only : sigma_HI, sigma_HeI, sigma_HeII
  implicit none
  private

  public :: gas_rates_compute, gas_rates_write
  public :: gamma_HI, gamma_HeI, gamma_HeII, heat_HI, heat_HeI, heat_HeII
  public :: heat_dust, g0_fuv
  public :: run_converged, run_iters, run_final_dx, run_final_dte
  public :: secion_apply
  public :: sec_dgamma_HI, sec_dgamma_HeI
  public :: sec_heat_HI, sec_heat_HeI, sec_heat_HeII

  real(kind=wp), allocatable :: gamma_HI(:), gamma_HeI(:), gamma_HeII(:)
  real(kind=wp), allocatable :: heat_HI(:),  heat_HeI(:),  heat_HeII(:)

  !--- secondary ionization by fast photoelectrons (par%use_sec_ion;
  !--- Shull & van Steenberg 1985).  Only H I, He I and He II are treated
  !--- as absorbers (the ionizing-band gas absorbers); metals are trace
  !--- and are ignored as absorbers here, a deliberate approximation.
  !--- x-INDEPENDENT per-atom band integrals filled in gas_rates_compute
  !--- when the switch is on: the hard-part (E0 > 40 eV) excess-energy
  !--- heating, and the secondary-ionization potentials from each absorber.
  real(kind=wp), parameter :: E_SEC_ION = 40.0_wp
  real(kind=wp), allocatable :: heat_hard_HI(:), heat_hard_HeI(:), heat_hard_HeII(:)
  real(kind=wp), allocatable :: si_HI_HI(:),  si_HI_HeI(:),  si_HI_HeII(:)
  real(kind=wp), allocatable :: si_HeI_HI(:), si_HeI_HeI(:), si_HeI_HeII(:)
  !--- x-DEPENDENT effective rates filled by secion_apply, or set to the
  !--- trivial (switch-off) values at the end of every gas_rates_compute so
  !--- the call sites can use them unconditionally: additive secondary
  !--- ionization rate [s^-1/atom] for the balance, and effective heating
  !--- [erg/s/atom] replacing heat_* in the thermal sum.
  real(kind=wp), allocatable :: sec_dgamma_HI(:), sec_dgamma_HeI(:)
  real(kind=wp), allocatable :: sec_heat_HI(:), sec_heat_HeI(:), sec_heat_HeII(:)
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

contains

  !=========================================================================
  subroutine gas_rates_compute()
    implicit none
    integer  :: il, inu, nleaf
    real(kind=wp) :: vol, fac, hnu, sHI, sHeI, sHeII, fJ
    real(kind=wp) :: sig_HI(nnu_band), sig_HeI(nnu_band), &
                     sig_HeII(nnu_band)

    nleaf = amr_grid%nleaf
    if (allocated(gamma_HI)) then
       if (size(gamma_HI) /= nleaf) &
          deallocate(gamma_HI, gamma_HeI, gamma_HeII, &
                     heat_HI, heat_HeI, heat_HeII, heat_dust, &
                     g0_fuv, &
                     heat_hard_HI, heat_hard_HeI, heat_hard_HeII, &
                     si_HI_HI, si_HI_HeI, si_HI_HeII, &
                     si_HeI_HI, si_HeI_HeI, si_HeI_HeII, &
                     sec_dgamma_HI, sec_dgamma_HeI, &
                     sec_heat_HI, sec_heat_HeI, sec_heat_HeII)     ! re-refinement
    end if
    if (.not. allocated(gamma_HI)) &
       allocate(gamma_HI(nleaf), gamma_HeI(nleaf), gamma_HeII(nleaf), &
                heat_HI(nleaf),  heat_HeI(nleaf),  heat_HeII(nleaf), &
                heat_dust(nleaf), g0_fuv(nleaf), &
                heat_hard_HI(nleaf), heat_hard_HeI(nleaf), heat_hard_HeII(nleaf), &
                si_HI_HI(nleaf), si_HI_HeI(nleaf), si_HI_HeII(nleaf), &
                si_HeI_HI(nleaf), si_HeI_HeI(nleaf), si_HeI_HeII(nleaf), &
                sec_dgamma_HI(nleaf), sec_dgamma_HeI(nleaf), &
                sec_heat_HI(nleaf), sec_heat_HeI(nleaf), sec_heat_HeII(nleaf))
    gamma_HI = 0.0_wp;  gamma_HeI = 0.0_wp;  gamma_HeII = 0.0_wp
    heat_HI  = 0.0_wp;  heat_HeI  = 0.0_wp;  heat_HeII  = 0.0_wp
    heat_dust = 0.0_wp; g0_fuv = 0.0_wp
    heat_hard_HI = 0.0_wp;  heat_hard_HeI = 0.0_wp;  heat_hard_HeII = 0.0_wp
    si_HI_HI  = 0.0_wp;  si_HI_HeI  = 0.0_wp;  si_HI_HeII  = 0.0_wp
    si_HeI_HI = 0.0_wp;  si_HeI_HeI = 0.0_wp;  si_HeI_HeII = 0.0_wp

    !--- cross sections depend only on the (fixed) band, not on the leaf:
    !--- evaluate once per bin instead of nleaf times (identical values).
    do inu = 1, nnu_band
       sig_HI(inu)   = sigma_HI(ion_e(inu))
       sig_HeI(inu)  = sigma_HeI(ion_e(inu))
       sig_HeII(inu) = sigma_HeII(ion_e(inu))
    end do

    do il = 1, nleaf
       vol = (2.0_wp * leaf_half(il))**3
       fac = 1.0_wp / (vol * par%distance2cm**2)     ! -> 4 pi J_nu dnu per bin
       do inu = 1, nnu_band
          fJ  = jt_ion(inu, il) * fac                ! 4 pi J dnu [erg/s/cm^2]
          if (fJ <= 0.0_wp) cycle
          hnu = ion_e(inu) * ev2erg                  ! photon energy [erg]
          sHI   = sig_HI(inu)
          sHeI  = sig_HeI(inu)
          sHeII = sig_HeII(inu)
          gamma_HI(il)   = gamma_HI(il)   + fJ*sHI  /hnu
          gamma_HeI(il)  = gamma_HeI(il)  + fJ*sHeI /hnu
          gamma_HeII(il) = gamma_HeII(il) + fJ*sHeII/hnu
          heat_HI(il)   = heat_HI(il)   + fJ*sHI  *(1.0_wp - eth_HI  /ion_e(inu))
          heat_HeI(il)  = heat_HeI(il)  + fJ*sHeI *(1.0_wp - eth_HeI /ion_e(inu))
          heat_HeII(il) = heat_HeII(il) + fJ*sHeII*(1.0_wp - eth_HeII/ion_e(inu))
          !--- secondary-ionization band integrals (par%use_sec_ion): the
          !--- hard part (photoelectron energy E0 = h nu - E_th > 40 eV) of
          !--- the excess-energy heating, and the secondary H I / He I
          !--- ionization potentials from each absorber (h nu / E_th,sec
          !--- ionizations of the secondary species).  x-independent.
          if (par%use_sec_ion) then
             if (ion_e(inu) > eth_HI + E_SEC_ION) then
                heat_hard_HI(il) = heat_hard_HI(il) &
                   + fJ*sHI*(1.0_wp - eth_HI/ion_e(inu))
                si_HI_HI(il)  = si_HI_HI(il)  &
                   + (fJ*sHI/hnu)*(ion_e(inu) - eth_HI)/eth_HI
                si_HeI_HI(il) = si_HeI_HI(il) &
                   + (fJ*sHI/hnu)*(ion_e(inu) - eth_HI)/eth_HeI
             end if
             if (ion_e(inu) > eth_HeI + E_SEC_ION) then
                heat_hard_HeI(il) = heat_hard_HeI(il) &
                   + fJ*sHeI*(1.0_wp - eth_HeI/ion_e(inu))
                si_HI_HeI(il)  = si_HI_HeI(il)  &
                   + (fJ*sHeI/hnu)*(ion_e(inu) - eth_HeI)/eth_HI
                si_HeI_HeI(il) = si_HeI_HeI(il) &
                   + (fJ*sHeI/hnu)*(ion_e(inu) - eth_HeI)/eth_HeI
             end if
             if (ion_e(inu) > eth_HeII + E_SEC_ION) then
                heat_hard_HeII(il) = heat_hard_HeII(il) &
                   + fJ*sHeII*(1.0_wp - eth_HeII/ion_e(inu))
                si_HI_HeII(il)  = si_HI_HeII(il)  &
                   + (fJ*sHeII/hnu)*(ion_e(inu) - eth_HeII)/eth_HI
                si_HeI_HeII(il) = si_HeI_HeII(il) &
                   + (fJ*sHeII/hnu)*(ion_e(inu) - eth_HeII)/eth_HeI
             end if
          end if
          !--- EUV grain heating: 4 pi J dnu x kappa_abs,dust [per cm].
          if (par%ion_add_dust) then
             block
               use gas_opacity_mod, only : ion_dust_sabs
               heat_dust(il) = heat_dust(il) + fJ &
                  * amr_grid%rhokap(il)*ion_dust_sabs(inu)/par%distance2cm
             end block
          end if
          !--- Habing-unit FUV field (bins below the H threshold).
          if (par%add_fuv .and. ion_e(inu) < par%eion_min) &
             g0_fuv(il) = g0_fuv(il) + fJ/1.6e-3_wp
       end do
    end do

    !--- trivial (switch-off) values so the ion/thermal call sites can use
    !--- these unconditionally; secion_apply overwrites them when the
    !--- switch is on.  With use_sec_ion off, sec_dgamma_* = 0 and
    !--- sec_heat_* = heat_* exactly, so the balance is unchanged.
    sec_dgamma_HI = 0.0_wp;  sec_dgamma_HeI = 0.0_wp
    sec_heat_HI = heat_HI;  sec_heat_HeI = heat_HeI;  sec_heat_HeII = heat_HeII
  end subroutine gas_rates_compute

  !=========================================================================
  ! Shull & van Steenberg (1985, ApJ 298, 268) high-energy partition of a
  ! fast photoelectron's excess energy; x = ionized fraction of the H+He
  ! nuclei.  f_heat is the heat fraction, f_ion_HI / f_ion_HeI the
  ! fractions driving secondary H I / He I ionization.  (The SvS85 Ly-alpha
  ! excitation channel is assumed to escape as line radiation and is not
  ! put into any rate.)
  elemental real(kind=wp) function f_heat(x) result(f)
    real(kind=wp), intent(in) :: x
    f = 0.9971_wp*(1.0_wp - (1.0_wp - x**0.2663_wp)**1.3163_wp)
  end function f_heat

  elemental real(kind=wp) function f_ion_HI(x) result(f)
    real(kind=wp), intent(in) :: x
    f = 0.3908_wp*(1.0_wp - x**0.4092_wp)**1.7592_wp
  end function f_ion_HI

  elemental real(kind=wp) function f_ion_HeI(x) result(f)
    real(kind=wp), intent(in) :: x
    f = 0.0554_wp*(1.0_wp - x**0.4614_wp)**1.6660_wp
  end function f_ion_HeI

  !=========================================================================
  ! Apply the x-dependent SvS85 partition on the CURRENT (lagged) gas state
  ! (par%use_sec_ion; called from main.f90 after gas_rates_compute).  The
  ! secondary-ionization rate for the balance and the effective (reduced)
  ! photoheating for the thermal sum are built from the x-independent band
  ! integrals filled in gas_rates_compute.
  !=========================================================================
  subroutine secion_apply()
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII
    implicit none
    integer :: il, nleaf
    real(kind=wp) :: nH, yHe, xHI, xHeI, xHeII, x, fh, fiHI, fiHeI
    real(kind=wp) :: nHI, nHeI, nHeII, R_secHI, R_secHeI

    nleaf = amr_grid%nleaf
    yHe = par%He_abund
    do il = 1, nleaf
       nH   = gas_nH(il)
       xHI  = gas_xHI(il);  xHeI = gas_xHeI(il);  xHeII = gas_xHeII(il)
       !--- x = (n_HII + n_HeII + n_HeIII) / (n_H + n_He), clamped to [0,1].
       x = min(max(((1.0_wp - xHI) + yHe*(1.0_wp - xHeI))/(1.0_wp + yHe), &
                   0.0_wp), 1.0_wp)
       fh   = f_heat(x);  fiHI = f_ion_HI(x);  fiHeI = f_ion_HeI(x)
       nHI  = nH*xHI;  nHeI = nH*yHe*xHeI;  nHeII = nH*yHe*xHeII
       R_secHI  = fiHI *(nHI*si_HI_HI(il)  + nHeI*si_HI_HeI(il)  &
                         + nHeII*si_HI_HeII(il))
       R_secHeI = fiHeI*(nHI*si_HeI_HI(il) + nHeI*si_HeI_HeI(il) &
                         + nHeII*si_HeI_HeII(il))
       !--- per-atom; x -> 1 gives f_ion -> 0 so these -> 0 (the max() is
       !--- only a divide-by-zero guard for a fully depleted cell).
       sec_dgamma_HI(il)  = R_secHI  / max(nHI,  1.0e-99_wp)
       sec_dgamma_HeI(il) = R_secHeI / max(nHeI, 1.0e-99_wp)
       sec_heat_HI(il)   = heat_HI(il)   - (1.0_wp - fh)*heat_hard_HI(il)
       sec_heat_HeI(il)  = heat_HeI(il)  - (1.0_wp - fh)*heat_hard_HeI(il)
       sec_heat_HeII(il) = heat_HeII(il) - (1.0_wp - fh)*heat_hard_HeII(il)
    end do
  end subroutine secion_apply

  !=========================================================================
  subroutine gas_rates_write()
    use iofile_mod
    use utility, only : get_base_name
    implicit none
    type(io_file_type) :: file
    character(len=192) :: filename
    real(kind=wp), allocatable :: leafxyz(:,:), leafsize(:), jnu(:,:)
    integer :: status, il, inu, ic, nleaf

    if (mpar%p_rank /= 0) return
    nleaf = amr_grid%nleaf

    status = 0
    filename = trim(get_base_name(par%out_file))//'_rates'// &
               trim(io_file_extension(par%file_format))
    call io_open_new(file, trim(filename), status)

    call io_append_image(file, gamma_HI, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Gamma_HI','H I photoionization rate [s^-1]',status)
    call io_put_keyword(file,'NNU_ION', par%nnu_ion, 'ionizing frequency bins', status)
    call io_put_keyword(file,'NNU_BAND', nnu_band, 'total band bins (ionizing+FUV)', status)
    call io_put_keyword(file,'EION_MIN',par%eion_min,'band lower edge [eV]',    status)
    call io_put_keyword(file,'EION_MAX',par%eion_max,'band upper edge [eV]',    status)
    call io_put_keyword(file,'TOT_LUM', par%luminosity,'band luminosity [erg/s]',status)
    call io_put_keyword(file,'DIST_CM', par%distance2cm,'distance unit (cm)',   status)
    call io_put_keyword(file,'nphotons',par%no_photons, 'number of photons',    status)
    call io_put_keyword(file,'CONVERGD',run_converged,'gas iteration converged',status)
    call io_put_keyword(file,'NITERDON',run_iters,   'gas iterations run',      status)
    call io_put_keyword(file,'FINALDX', run_final_dx, 'final max|delta x_HII|',  status)
    call io_put_keyword(file,'FINALDTE',run_final_dte,'final max|delta Te|/Te',  status)
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
       call io_put_keyword(file,'EXTNAME','Heat_dust','EUV grain heating [erg/s/cm^3]',status)
       block
         use dust_temp_mod, only : t_dust
         if (allocated(t_dust)) then
            call io_append_image(file, t_dust, status, bitpix=-64)
            call io_put_keyword(file,'EXTNAME','T_dust', &
                 'equilibrium mixture dust temperature [K]',status)
         end if
       end block
    end if

    !--- J_nu(nnu, nleaf) [erg/s/cm^2/Hz/sr]
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

    !--- ionization structure (fixed initialization when gas_niter = 0).
    block
      use gas_state_mod, only : gas_xHI, gas_xHeI, gas_xHeII, gas_ne, gas_Te
      real(kind=wp), allocatable :: tmp(:)
      allocate(tmp(nleaf))
      tmp = gas_xHI(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','x_HI','n_HI/n_H per leaf',status)
      tmp = gas_xHeI(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','x_HeI','n_HeI/n_He per leaf',status)
      tmp = gas_xHeII(1:nleaf)
      call io_append_image(file, tmp, status, bitpix=-64)
      call io_put_keyword(file,'EXTNAME','x_HeII','n_HeII/n_He per leaf',status)
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
         call species_write(file)
       end block
    end if

    call io_close(file, status)
    write(*,'(2a)') ' ION: rates written to: ', trim(filename)
  end subroutine gas_rates_write

end module gas_rates_mod
