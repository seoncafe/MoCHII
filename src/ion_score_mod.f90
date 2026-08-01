module ion_score_mod
!---------------------------------------------------------------------------
! Common ionizing-packet path scorer.
!
! Every direct analytic segment and every sampled scattered segment supplies
! the same path quantity
!
!   path_lum = Lpacket * weight * effective_path_length
!
! [erg/s * code length].  Every rate coefficient is scored at the packet's own
! photon%energy_eV; those reduced arrays are authoritative for gas_rates_mod,
! while the jt_ion(inu,leaf) tally is only a diagnostic spectrum.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: ion_score_setup, ion_score_resize, ion_score_reset
  public :: score_ion_path, ion_score_reduce
  public :: ion_score_apply_hhe
  public :: ion_score_apply_dust_heat
  public :: ion_score_apply_sec
  public :: ion_score_apply_g0
  public :: ion_score_energy_packet, ion_score_energy_totals
  public :: ion_score_report_energy_closure, ion_score_report_metal_rates

  integer, parameter :: I_GHI = 1, I_GHEI = 2, I_GHEII = 3
  integer, parameter :: I_HHI = 4, I_HHEI = 5, I_HHEII = 6
  integer, parameter :: NHHE = 6
  !--- secondary-ionization (par%use_sec_ion) path-estimator layout: the nine
  !--- x-independent band integrals of gas_rates_mod -- three hard-part
  !--- (photoelectron energy E0 > E_SEC_ION) heating terms and the six
  !--- secondary H I / He I ionization potentials, one per (secondary,absorber)
  !--- pair.  E_SEC_ION matches the gas_rates_mod value.
  integer, parameter :: I_HH_HI = 1, I_HH_HEI = 2, I_HH_HEII = 3
  integer, parameter :: I_SI_HI_HI = 4, I_SI_HI_HEI = 5, I_SI_HI_HEII = 6
  integer, parameter :: I_SI_HEI_HI = 7, I_SI_HEI_HEI = 8, I_SI_HEI_HEII = 9
  integer, parameter :: NSEC = 9
  real(kind=wp), parameter :: E_SEC_ION = 40.0_wp
  integer, parameter :: ENERGY_ACC_WP = selected_real_kind(18,300)

  real(kind=wp), allocatable :: hhe_rate_path(:,:) ! (NHHE,nleaf), rank local then
                                                   ! reduced; Sum(path_lum * coeff(E))
                                                   ! for the three photoionization
                                                   ! rates and the three photoheating
                                                   ! integrals.
  real(kind=wp), allocatable :: dust_abs_path(:)   ! (nleaf), rank local then reduced;
                                                   ! Sum(path_lum * sigma_abs,dust(E))
                                                   ! per leaf, the EUV grain-heating
                                                   ! integral.
  real(kind=wp), allocatable :: sec_ion_path(:,:)  ! (NSEC,nleaf), rank local then
                                                   ! reduced; the nine x-independent
                                                   ! secondary-ionization band
                                                   ! integrals.
  real(kind=wp), allocatable :: fuv_field_path(:)  ! (nleaf), rank local then
                                                   ! reduced; Sum(path_lum) over
                                                   ! FUV segments (E < eion_min),
                                                   ! the Habing G0 integral.
  integer :: score_nleaf = 0
  integer(kind=int64) :: direct_segments = 0_int64
  integer(kind=int64) :: scattered_segments = 0_int64
  real(kind=wp) :: energy_emitted = 0.0_wp
  real(kind=wp) :: energy_absorbed = 0.0_wp
  real(kind=wp) :: energy_escaped = 0.0_wp
  real(kind=ENERGY_ACC_WP) :: energy_emitted_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_absorbed_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_escaped_acc = 0.0_ENERGY_ACC_WP

contains

  subroutine ion_score_setup(nleaf)
    implicit none
    integer, intent(in) :: nleaf

    if (allocated(hhe_rate_path)) deallocate(hhe_rate_path)
    if (allocated(dust_abs_path)) deallocate(dust_abs_path)
    if (allocated(sec_ion_path)) deallocate(sec_ion_path)
    if (allocated(fuv_field_path)) deallocate(fuv_field_path)
    score_nleaf = nleaf
    allocate(hhe_rate_path(NHHE,nleaf))
    hhe_rate_path = 0.0_wp

    !--- EUV grain-heating integral: allocated only when dust is on.  Left
    !--- unallocated otherwise, so score_ion_path, the reduction, and the reset
    !--- are no-ops on the dust-off path.
    if (par%ion_add_dust) then
       allocate(dust_abs_path(nleaf))
       dust_abs_path = 0.0_wp
    end if

    !--- secondary-ionization integrals: allocated only when use_sec_ion is on.
    if (par%use_sec_ion) then
       allocate(sec_ion_path(NSEC,nleaf))
       sec_ion_path = 0.0_wp
    end if

    !--- Habing-G0 integral: allocated only when the FUV band is on.
    if (par%add_fuv) then
       allocate(fuv_field_path(nleaf))
       fuv_field_path = 0.0_wp
    end if

    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_setup
         call species_rate_path_setup(nleaf)
       end block
    end if
  end subroutine ion_score_setup

  subroutine ion_score_resize(nleaf)
    integer, intent(in) :: nleaf
    call ion_score_setup(nleaf)
  end subroutine ion_score_resize

  subroutine ion_score_reset()
    if (allocated(hhe_rate_path)) hhe_rate_path = 0.0_wp
    if (allocated(dust_abs_path)) dust_abs_path = 0.0_wp
    if (allocated(sec_ion_path)) sec_ion_path = 0.0_wp
    if (allocated(fuv_field_path)) fuv_field_path = 0.0_wp
    direct_segments = 0_int64
    scattered_segments = 0_int64
    energy_emitted = 0.0_wp
    energy_absorbed = 0.0_wp
    energy_escaped = 0.0_wp
    energy_emitted_acc = 0.0_ENERGY_ACC_WP
    energy_absorbed_acc = 0.0_ENERGY_ACC_WP
    energy_escaped_acc = 0.0_ENERGY_ACC_WP
    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_reset
         call species_rate_path_reset()
       end block
    end if
  end subroutine ion_score_reset

  subroutine score_ion_path(photon, il, path_lum)
    use jtally_mod,       only : jt_ion
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR, &
                                 dust_ir_absorb_path
    implicit none
    type(photon_type), intent(in) :: photon
    integer,           intent(in) :: il
    real(kind=wp),     intent(in) :: path_lum
    integer :: inu
    real(kind=wp) :: energy, coeff(NHHE)

    if (path_lum == 0.0_wp) return
    !--- Dust thermal re-emission packets sit below every gas photoionization
    !--- threshold, so they carry no photoionization, no photoheating, no
    !--- secondary ionization and no Habing G0, and their energies lie outside
    !--- the ionizing-band bin grid jt_ion is defined on.  Their only
    !--- interaction is grain absorption, which is scored into the infrared
    !--- band's own estimator; every ionizing-band estimator is left untouched.
    if (transported_band == BAND_DUST_IR) then
       call dust_ir_absorb_path(il, path_lum*photon%ionphys%dust_abs)
       return
    end if
    inu = photon%inu
    jt_ion(inu,il) = jt_ion(inu,il) + path_lum
    !--- EUV grain-heating integral: accumulate the packet's dust absorption
    !--- cross section, evaluated at its exact energy.
    if (allocated(dust_abs_path)) &
       dust_abs_path(il) = dust_abs_path(il) + path_lum*photon%ionphys%dust_abs
    !--- the nine x-independent secondary-ionization band integrals, each
    !--- accumulated at the packet's own energy (hnu = E*ev2erg).
    if (allocated(sec_ion_path)) then
       block
         real(kind=wp) :: en, xHI, xHeI, xHeII, hnu
         en    = photon%ionphys%energy_eV
         xHI   = photon%ionphys%sigma_HI
         xHeI  = photon%ionphys%sigma_HeI
         xHeII = photon%ionphys%sigma_HeII
         hnu   = en*ev2erg
         if (en > eth_HI + E_SEC_ION) then
            sec_ion_path(I_HH_HI,il)     = sec_ion_path(I_HH_HI,il)     &
               + path_lum*xHI*(1.0_wp - eth_HI/en)
            sec_ion_path(I_SI_HI_HI,il)  = sec_ion_path(I_SI_HI_HI,il)  &
               + path_lum*(xHI/hnu)*(en - eth_HI)/eth_HI
            sec_ion_path(I_SI_HEI_HI,il) = sec_ion_path(I_SI_HEI_HI,il) &
               + path_lum*(xHI/hnu)*(en - eth_HI)/eth_HeI
         end if
         if (en > eth_HeI + E_SEC_ION) then
            sec_ion_path(I_HH_HEI,il)     = sec_ion_path(I_HH_HEI,il)     &
               + path_lum*xHeI*(1.0_wp - eth_HeI/en)
            sec_ion_path(I_SI_HI_HEI,il)  = sec_ion_path(I_SI_HI_HEI,il)  &
               + path_lum*(xHeI/hnu)*(en - eth_HeI)/eth_HI
            sec_ion_path(I_SI_HEI_HEI,il) = sec_ion_path(I_SI_HEI_HEI,il) &
               + path_lum*(xHeI/hnu)*(en - eth_HeI)/eth_HeI
         end if
         if (en > eth_HeII + E_SEC_ION) then
            sec_ion_path(I_HH_HEII,il)     = sec_ion_path(I_HH_HEII,il)     &
               + path_lum*xHeII*(1.0_wp - eth_HeII/en)
            sec_ion_path(I_SI_HI_HEII,il)  = sec_ion_path(I_SI_HI_HEII,il)  &
               + path_lum*(xHeII/hnu)*(en - eth_HeII)/eth_HI
            sec_ion_path(I_SI_HEI_HEII,il) = sec_ion_path(I_SI_HEI_HEII,il) &
               + path_lum*(xHeII/hnu)*(en - eth_HeII)/eth_HeI
         end if
       end block
    end if
    !--- Habing-G0 integral: accumulate path_lum for FUV segments, i.e. those
    !--- whose photon energy sits below the ionization floor.
    if (allocated(fuv_field_path)) then
       if (photon%ionphys%energy_eV < par%eion_min) &
          fuv_field_path(il) = fuv_field_path(il) + path_lum
    end if
    if (photon%nscatt > 0) then
       scattered_segments = scattered_segments + 1_int64
    else
       direct_segments = direct_segments + 1_int64
    end if
    energy = photon%ionphys%energy_eV
    coeff(I_GHI)   = photon%ionphys%sigma_HI/(energy*ev2erg)
    coeff(I_GHEI)  = photon%ionphys%sigma_HeI/(energy*ev2erg)
    coeff(I_GHEII) = photon%ionphys%sigma_HeII/(energy*ev2erg)
    coeff(I_HHI)   = photon%ionphys%sigma_HI * &
                     max(1.0_wp-eth_HI/energy,0.0_wp)
    coeff(I_HHEI)  = photon%ionphys%sigma_HeI * &
                     max(1.0_wp-eth_HeI/energy,0.0_wp)
    coeff(I_HHEII) = photon%ionphys%sigma_HeII * &
                     max(1.0_wp-eth_HeII/energy,0.0_wp)
    hhe_rate_path(:,il) = hhe_rate_path(:,il) + path_lum*coeff
    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_score
         call species_rate_path_score(photon%ionphys%energy_eV, &
              photon%ionphys%sigma_metal, photon%ionphys%n_metal, &
              il, path_lum)
       end block
    end if
  end subroutine score_ion_path

  subroutine ion_score_energy_packet(photon, tau_edge)
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR, &
                                 dust_ir_escape_packet
    type(photon_type), intent(in) :: photon
    real(kind=wp), intent(in) :: tau_edge
    real(kind=wp) :: luminosity, transmission
    !--- an infrared packet belongs to the dust thermal band's own ledger, not
    !--- to the ionizing band's.
    if (transported_band == BAND_DUST_IR) then
       call dust_ir_escape_packet(photon, tau_edge)
       return
    end if
    luminosity = photon%Lpacket*photon%wgt
    transmission = exp(-max(tau_edge,0.0_wp))
    energy_emitted_acc = energy_emitted_acc + real(luminosity, ENERGY_ACC_WP)
    energy_absorbed_acc = energy_absorbed_acc + real( &
       luminosity*(1.0_wp-transmission), ENERGY_ACC_WP)
    energy_escaped_acc = energy_escaped_acc + real( &
       luminosity*transmission, ENERGY_ACC_WP)
  end subroutine ion_score_energy_packet

  subroutine ion_score_reduce()
    use mpi
    implicit none
    integer :: ierr, nchunk, i0, n
    integer(kind=int64) :: counts(2)
    real(kind=wp) :: energy_totals(3)

    nchunk = max(1, 100000000/NHHE)
    i0 = 1
    do while (i0 <= score_nleaf)
       n = min(nchunk, score_nleaf-i0+1)
       call MPI_ALLREDUCE(MPI_IN_PLACE, hhe_rate_path(:,i0:i0+n-1), NHHE*n, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       i0 = i0 + n
    end do
    if (allocated(dust_abs_path)) &
       call MPI_ALLREDUCE(MPI_IN_PLACE, dust_abs_path, score_nleaf, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (allocated(sec_ion_path)) &
       call MPI_ALLREDUCE(MPI_IN_PLACE, sec_ion_path, NSEC*score_nleaf, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (allocated(fuv_field_path)) &
       call MPI_ALLREDUCE(MPI_IN_PLACE, fuv_field_path, score_nleaf, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    counts = [direct_segments, scattered_segments]
    call MPI_ALLREDUCE(MPI_IN_PLACE, counts, 2, MPI_INTEGER8, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    direct_segments = counts(1)
    scattered_segments = counts(2)
    energy_totals = real( &
       [energy_emitted_acc, energy_absorbed_acc, energy_escaped_acc], wp)
    call MPI_ALLREDUCE(MPI_IN_PLACE, energy_totals, 3, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    energy_emitted = energy_totals(1)
    energy_absorbed = energy_totals(2)
    energy_escaped = energy_totals(3)
    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_reduce
         call species_rate_path_reduce()
       end block
    end if
  end subroutine ion_score_reduce

  subroutine ion_score_apply_hhe(gHI, gHeI, gHeII, hHI, hHeI, hHeII)
    use octree_mod, only : leaf_half
    implicit none
    real(kind=wp), intent(out) :: gHI(:), gHeI(:), gHeII(:)
    real(kind=wp), intent(out) :: hHI(:), hHeI(:), hHeII(:)
    real(kind=wp) :: fac
    integer :: il

    if (.not. allocated(hhe_rate_path)) &
       error stop 'H/He path estimators are not allocated'
    do il = 1, score_nleaf
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       gHI(il)   = hhe_rate_path(I_GHI,il)*fac
       gHeI(il)  = hhe_rate_path(I_GHEI,il)*fac
       gHeII(il) = hhe_rate_path(I_GHEII,il)*fac
       hHI(il)   = hhe_rate_path(I_HHI,il)*fac
       hHeI(il)  = hhe_rate_path(I_HHEI,il)*fac
       hHeII(il) = hhe_rate_path(I_HHEII,il)*fac
    end do
  end subroutine ion_score_apply_hhe

  subroutine ion_score_apply_dust_heat(heat_dust)
    !--- EUV grain heating [erg/s/cm^3] per leaf.  dust_abs_path
    !--- = Sum(path_lum * sigma_abs,dust(E)) with the factor below is the
    !--- transport form of
    !---   heat_dust = 4 pi J_abs rhokap / distance2cm,
    !--- since jt_ion = Sum path_lum gives 4 pi J dnu after the same factor.
    use octree_mod, only : leaf_half, amr_grid
    implicit none
    real(kind=wp), intent(out) :: heat_dust(:)
    real(kind=wp) :: fac
    integer :: il

    if (.not. allocated(dust_abs_path)) &
       error stop 'dust-absorption path estimator is not allocated'
    do il = 1, score_nleaf
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       heat_dust(il) = dust_abs_path(il)*fac &
                       * amr_grid%rhokap(il)/par%distance2cm
    end do
  end subroutine ion_score_apply_dust_heat

  subroutine ion_score_apply_g0(g0_fuv)
    !--- Habing FUV field G0 [dimensionless] per leaf.  fuv_field_path
    !--- = Sum(path_lum) over the segments carried by photons below eion_min,
    !--- and the factor below turns it into 4 pi J_FUV / 1.6e-3 erg/s/cm^2.
    use octree_mod, only : leaf_half
    implicit none
    real(kind=wp), intent(out) :: g0_fuv(:)
    real(kind=wp) :: fac
    integer :: il
    if (.not. allocated(fuv_field_path)) &
       error stop 'FUV-field path estimator is not allocated'
    do il = 1, score_nleaf
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       g0_fuv(il) = fuv_field_path(il)*fac/1.6e-3_wp
    end do
  end subroutine ion_score_apply_g0

  subroutine ion_score_energy_totals(emitted, absorbed, escaped)
    real(kind=wp), intent(out) :: emitted, absorbed, escaped
    emitted = energy_emitted
    absorbed = energy_absorbed
    escaped = energy_escaped
  end subroutine ion_score_energy_totals

  subroutine ion_score_report_energy_closure()
    !--- Energy budget of the transported packets: emitted = absorbed + escaped
    !--- per packet by construction (each packet's luminosity is split by its
    !--- edge transmission), so the residual measures only the accumulation
    !--- error of the reduced sums.
    use mpi
    implicit none
    real(kind=wp) :: closure
    integer :: ierr

    closure = abs(energy_absorbed+energy_escaped-energy_emitted) / &
              max(energy_emitted,tiny(1.0_wp))
    if (mpar%p_rank == 0) then
       write(*,'(a,2(i0,1x))') &
          ' ION: continuous scored path segments direct/scattered = ', &
          direct_segments, scattered_segments
       write(*,'(a)') ' ION: continuous H/He direct rates: ACTIVE'
       write(*,'(a,3es14.6)') &
          ' ION: continuous energy emitted/absorbed/escaped = ', &
          energy_emitted, energy_absorbed, energy_escaped
       write(*,'(a,es11.3)') ' ION: continuous energy closure = ', closure
    end if
    if (closure > 5.0e-13_wp) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end subroutine ion_score_report_energy_closure


  subroutine ion_score_apply_sec(heat_hard_HI, heat_hard_HeI, heat_hard_HeII, &
       si_HI_HI, si_HI_HeI, si_HI_HeII, si_HeI_HI, si_HeI_HeI, si_HeI_HeII)
    !--- Secondary-ionization band integrals per leaf.  sec_ion_path
    !--- = Sum(path_lum * coeff(E)) with the factor below is the transport form
    !--- of fJ * coeff(E), since Sum(path_lum) times that factor is 4 pi J dnu.
    !--- gas_rates_mod's secion_apply then folds the gas state into these nine
    !--- arrays.
    use octree_mod, only : leaf_half
    implicit none
    real(kind=wp), intent(out) :: heat_hard_HI(:), heat_hard_HeI(:), heat_hard_HeII(:)
    real(kind=wp), intent(out) :: si_HI_HI(:),  si_HI_HeI(:),  si_HI_HeII(:)
    real(kind=wp), intent(out) :: si_HeI_HI(:), si_HeI_HeI(:), si_HeI_HeII(:)
    real(kind=wp) :: fac
    integer :: il

    if (.not. allocated(sec_ion_path)) &
       error stop 'secondary-ionization path estimators are not allocated'
    do il = 1, score_nleaf
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       heat_hard_HI(il)   = sec_ion_path(I_HH_HI,il)*fac
       heat_hard_HeI(il)  = sec_ion_path(I_HH_HEI,il)*fac
       heat_hard_HeII(il) = sec_ion_path(I_HH_HEII,il)*fac
       si_HI_HI(il)   = sec_ion_path(I_SI_HI_HI,il)*fac
       si_HI_HeI(il)  = sec_ion_path(I_SI_HI_HEI,il)*fac
       si_HI_HeII(il) = sec_ion_path(I_SI_HI_HEII,il)*fac
       si_HeI_HI(il)   = sec_ion_path(I_SI_HEI_HI,il)*fac
       si_HeI_HeI(il)  = sec_ion_path(I_SI_HEI_HEI,il)*fac
       si_HeI_HeII(il) = sec_ion_path(I_SI_HEI_HEII,il)*fac
    end do
  end subroutine ion_score_apply_sec


  subroutine ion_score_report_metal_rates()
    if (.not. par%use_metals) return
    if (mpar%p_rank == 0) &
       write(*,'(a)') ' ION: continuous metal direct rates: ACTIVE'
  end subroutine ion_score_report_metal_rates

end module ion_score_mod
