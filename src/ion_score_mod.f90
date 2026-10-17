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
  use memory_mod, only : shared_host_rank, shared_host_size
  use photoelectron_degradation_mod, only : N_EDEP, E_SECONDARY, &
                                            photoelectron_energy_weights
  implicit none
  private

  public :: ion_score_setup, ion_score_resize, ion_score_reset
  public :: score_ion_path, ion_score_reduce
  public :: ion_score_apply_hhe
  public :: ion_score_apply_dust_heat
  public :: ion_score_apply_sec
  public :: ion_score_energy_packet, ion_score_energy_totals
  public :: ion_score_report_energy_closure, ion_score_report_metal_rates
  public :: ion_score_report_ledger, ion_score_ledger_leaf
  public :: LEDGER_NSRC

  integer, parameter :: I_GHI = 1, I_GHEI = 2, I_GHEII = 3
  integer, parameter :: I_HHI = 4, I_HHEI = 5, I_HHEII = 6
  integer, parameter :: NHHE = 6
  !--- secondary-ionization (par%use_sec_ion) path-estimator layout: for each
  !--- absorber (H I, He I, He II) the photoelectron energy deposited above
  !--- E_SECONDARY, sigma (1 - E_th/E), split onto the N_EDEP energy nodes of
  !--- photoelectron_degradation_mod by the x-independent weights of E0, so
  !--- the state-dependent partition is applied afterwards exactly.  Row
  !--- (a-1)*N_EDEP + k is absorber a, node k.
  integer, parameter :: NSEC = 3*N_EDEP
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
                                                   ! reduced; the photoelectron
                                                   ! energy deposited above 30 eV,
                                                   ! by absorber and energy node.
  integer :: score_nleaf = 0
  integer(kind=int64) :: direct_segments = 0_int64
  integer(kind=int64) :: scattered_segments = 0_int64
  real(kind=wp) :: energy_emitted = 0.0_wp
  real(kind=wp) :: energy_absorbed = 0.0_wp
  real(kind=wp) :: energy_escaped = 0.0_wp
  real(kind=ENERGY_ACC_WP) :: energy_emitted_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_absorbed_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_escaped_acc = 0.0_ENERGY_ACC_WP
  !--- absorption ledger (par%absorption_ledger) [photons/s]: by source
  !--- (photon%source_channel 0..LEDGER_NSRC-1) the photons emitted and
  !--- escaping at the edge, and the photons absorbed by each absorber,
  !--- path_lum kappa_a / (E erg) on the scored segments
  integer, parameter :: LEDGER_NSRC = 10, LEDGER_NABS = 5
  character(len=6), parameter :: LEDGER_ABS_NAME(LEDGER_NABS) = &
     ['H0    ', 'He0   ', 'He+   ', 'metals', 'dust  ']
  real(kind=ENERGY_ACC_WP) :: ledger_emitted(0:LEDGER_NSRC-1) = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ledger_escaped(0:LEDGER_NSRC-1) = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ledger_absorbed(LEDGER_NABS,0:LEDGER_NSRC-1) = 0.0_ENERGY_ACC_WP
  !--- the same for H0, He0, He+ in every leaf [photons/s]:
  !--- (absorber, source, leaf), rank local then reduced, written to the
  !--- rates file as Ledger_H0 / Ledger_He0 / Ledger_He+ (nleaf, source)
  integer, parameter :: LEDGER_NLEAF_ABS = 3
  real(kind=wp), allocatable :: ledger_leaf(:,:,:)

contains

  subroutine ion_score_setup(nleaf)
    implicit none
    integer, intent(in) :: nleaf

    if (allocated(hhe_rate_path)) deallocate(hhe_rate_path)
    if (allocated(dust_abs_path)) deallocate(dust_abs_path)
    if (allocated(sec_ion_path)) deallocate(sec_ion_path)
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

    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_setup
         call species_rate_path_setup(nleaf)
       end block
    end if
    if (allocated(ledger_leaf)) deallocate(ledger_leaf)
    if (par%absorption_ledger) then
       allocate(ledger_leaf(LEDGER_NLEAF_ABS, 0:LEDGER_NSRC-1, nleaf))
       ledger_leaf = 0.0_wp
    end if
    !--- He I 2^3S photoionization rate Phi_3 and its heating (the output
    !--- diagnostic par%hei_metastable and the coupling
    !--- par%hei_metastable_coupling)
    if (par%hei_metastable .or. par%hei_metastable_coupling) then
       block
         use hei_metastable_mod, only : hei_metastable_path_setup
         call hei_metastable_path_setup(nleaf)
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
    direct_segments = 0_int64
    scattered_segments = 0_int64
    energy_emitted = 0.0_wp
    energy_absorbed = 0.0_wp
    energy_escaped = 0.0_wp
    energy_emitted_acc = 0.0_ENERGY_ACC_WP
    energy_absorbed_acc = 0.0_ENERGY_ACC_WP
    energy_escaped_acc = 0.0_ENERGY_ACC_WP
    ledger_emitted = 0.0_ENERGY_ACC_WP
    ledger_escaped = 0.0_ENERGY_ACC_WP
    ledger_absorbed = 0.0_ENERGY_ACC_WP
    if (allocated(ledger_leaf)) ledger_leaf = 0.0_wp
    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_reset
         call species_rate_path_reset()
       end block
    end if
    if (par%hei_metastable .or. par%hei_metastable_coupling) then
       block
         use hei_metastable_mod, only : hei_metastable_path_reset
         call hei_metastable_path_reset()
       end block
    end if
  end subroutine ion_score_reset

  subroutine score_ion_path(photon, il, path_lum)
    use jtally_mod,       only : jt_ion, jt_ion_on
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
    !--- secondary ionization, and their energies lie outside
    !--- the ionizing-band bin grid jt_ion is defined on.  Their only
    !--- interaction is grain absorption, which is scored into the infrared
    !--- band's own estimator; every ionizing-band estimator is left untouched.
    if (transported_band == BAND_DUST_IR) then
       call dust_ir_absorb_path(il, path_lum*photon%ionphys%dust_abs)
       return
    end if
    inu = photon%inu
    if (jt_ion_on) jt_ion(inu,il) = jt_ion(inu,il) + path_lum
    !--- EUV grain-heating integral: accumulate the packet's dust absorption
    !--- cross section, evaluated at its exact energy.
    if (allocated(dust_abs_path)) &
       dust_abs_path(il) = dust_abs_path(il) + path_lum*photon%ionphys%dust_abs
    !--- photoelectron energy above E_SECONDARY, by absorber and energy node,
    !--- at the packet's own energy.
    if (allocated(sec_ion_path)) then
       block
         real(kind=wp) :: en, sig(3), eth(3), w(N_EDEP)
         integer :: ia
         en  = photon%ionphys%energy_eV
         sig = [photon%ionphys%sigma_HI, photon%ionphys%sigma_HeI, &
                photon%ionphys%sigma_HeII]
         eth = [eth_HI, eth_HeI, eth_HeII]
         do ia = 1, 3
            if (en - eth(ia) <= E_SECONDARY) cycle
            call photoelectron_energy_weights(en - eth(ia), w)
            sec_ion_path((ia-1)*N_EDEP+1:ia*N_EDEP, il) = &
               sec_ion_path((ia-1)*N_EDEP+1:ia*N_EDEP, il) &
               + path_lum*sig(ia)*(1.0_wp - eth(ia)/en)*w
         end do
       end block
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
              photon%ionphys%sigma_metal, photon%ionphys%heat_metal, &
              photon%ionphys%n_metal, il, path_lum)
       end block
    end if
    if (par%hei_metastable .or. par%hei_metastable_coupling) then
       block
         use hei_metastable_mod, only : hei_metastable_path_score
         call hei_metastable_path_score(il, photon%ionphys%energy_eV, path_lum)
       end block
    end if
    if (par%absorption_ledger) call ledger_score(photon, il, path_lum)
  end subroutine score_ion_path

  !--- photons absorbed on one scored segment, by absorber, with the
  !--- opacity terms of ion_packet_opacity (dust: absorption only)
  subroutine ledger_score(photon, il, path_lum)
    use octree_mod,    only : amr_grid
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII
    type(photon_type), intent(in) :: photon
    integer,           intent(in) :: il
    real(kind=wp),     intent(in) :: path_lum
    real(kind=wp) :: kap(LEDGER_NABS), nH
    integer :: is
    nH = gas_nH(il)
    kap = 0.0_wp
    kap(1) = nH*gas_xHI(il)*photon%ionphys%sigma_HI*par%distance2cm
    kap(2) = nH*par%abund_He*gas_xHeI(il)*photon%ionphys%sigma_HeI*par%distance2cm
    kap(3) = nH*par%abund_He*gas_xHeII(il)*photon%ionphys%sigma_HeII*par%distance2cm
    if (par%use_metals .and. par%ion_metal_abs) then
       block
         use species_mod, only : species_cached_opacity
         kap(4) = nH*species_cached_opacity(photon%ionphys%sigma_metal, &
                  photon%ionphys%n_metal, il)*par%distance2cm
       end block
    end if
    if (par%ion_add_dust) kap(5) = amr_grid%rhokap(il)*photon%ionphys%dust_abs
    is = min(max(photon%source_channel, 0), LEDGER_NSRC-1)
    kap = path_lum*kap/(photon%ionphys%energy_eV*ev2erg)
    ledger_absorbed(:,is) = ledger_absorbed(:,is) + real(kap, ENERGY_ACC_WP)
    ledger_leaf(:,is,il) = ledger_leaf(:,is,il) + kap(1:LEDGER_NLEAF_ABS)
  end subroutine ledger_score

  !--- photons/s absorbed by absorber ia (1 H0, 2 He0, 3 He+) in every leaf
  !--- from source is, after ion_score_reduce
  subroutine ion_score_ledger_leaf(ia, is, out)
    integer,       intent(in)  :: ia, is
    real(kind=wp), intent(out) :: out(:)
    out = ledger_leaf(ia, is, 1:size(out))
  end subroutine ion_score_ledger_leaf

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
    if (par%absorption_ledger) then
       block
         integer :: is
         real(kind=wp) :: nphot
         is = min(max(photon%source_channel, 0), LEDGER_NSRC-1)
         nphot = luminosity/(photon%energy_eV*ev2erg)
         ledger_emitted(is) = ledger_emitted(is) + real(nphot, ENERGY_ACC_WP)
         ledger_escaped(is) = ledger_escaped(is) + real(nphot*transmission, ENERGY_ACC_WP)
       end block
    end if
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
    if (par%absorption_ledger) call ledger_reduce()
    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_reduce
         call species_rate_path_reduce()
       end block
    end if
    if (par%hei_metastable .or. par%hei_metastable_coupling) then
       block
         use hei_metastable_mod, only : hei_metastable_path_reduce
         call hei_metastable_path_reduce()
       end block
    end if
  end subroutine ion_score_reduce

  !--- the ledger sums over ranks, in double precision (the per-rank sums
  !--- were accumulated in the extended kind)
  subroutine ledger_reduce()
    use mpi
    real(kind=wp) :: buf((2 + LEDGER_NABS)*LEDGER_NSRC)
    integer :: ierr, n
    n = LEDGER_NSRC
    buf(1:n)         = real(ledger_emitted, wp)
    buf(n+1:2*n)     = real(ledger_escaped, wp)
    buf(2*n+1:)      = reshape(real(ledger_absorbed, wp), [LEDGER_NABS*n])
    call MPI_ALLREDUCE(MPI_IN_PLACE, buf, size(buf), MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    ledger_emitted  = real(buf(1:n), ENERGY_ACC_WP)
    ledger_escaped  = real(buf(n+1:2*n), ENERGY_ACC_WP)
    ledger_absorbed = real(reshape(buf(2*n+1:), [LEDGER_NABS, n]), ENERGY_ACC_WP)
    block
      integer :: i0, nl, nchunk
      nchunk = max(1, 100000000/(LEDGER_NLEAF_ABS*LEDGER_NSRC))
      i0 = 1
      do while (i0 <= score_nleaf)
         nl = min(nchunk, score_nleaf-i0+1)
         call MPI_ALLREDUCE(MPI_IN_PLACE, ledger_leaf(:,:,i0:i0+nl-1), &
                            LEDGER_NLEAF_ABS*LEDGER_NSRC*nl, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, MPI_COMM_WORLD, ierr)
         i0 = i0 + nl
      end do
    end block
  end subroutine ledger_reduce

  !--- one row for each source that emitted: emitted, absorbed by each
  !--- absorber, escaped [photons/s], and the closure
  !--- (absorbed + escaped)/emitted - 1, which the path estimators meet in
  !--- expectation (exactly without dust scattering, to Monte Carlo noise
  !--- with it)
  subroutine ion_score_report_ledger()
    character(len=8), parameter :: SRC_NAME(0:LEDGER_NSRC-1) = &
       ['star    ', 'dif1 HI ', 'dif2 HeI', 'dif3HeII', 'dif4 HeI', &
        'dif5HeII', 'dif6 met', 'dif7 ff ', 'dif8 val', 'star FUV']
    integer :: is
    real(kind=wp) :: em, ab(LEDGER_NABS), es
    if (.not. par%absorption_ledger .or. mpar%p_rank /= 0) return
    write(*,'(a,5(3x,a6),3x,a,3x,a)') ' LEDGER: source   emitted      ', &
       LEDGER_ABS_NAME, 'escaped    ', '(abs+esc)/em-1'
    do is = 0, LEDGER_NSRC-1
       em = real(ledger_emitted(is), wp)
       if (em <= 0.0_wp) cycle
       ab = real(ledger_absorbed(:,is), wp)
       es = real(ledger_escaped(is), wp)
       write(*,'(a,a8,7es12.4e3,es11.2)') ' LEDGER: ', SRC_NAME(is), em, ab, es, &
          (sum(ab) + es)/em - 1.0_wp
    end do
  end subroutine ion_score_report_ledger

  subroutine ion_score_apply_hhe(gHI, gHeI, gHeII, hHI, hHeI, hHeII)
    use octree_mod, only : leaf_half
    implicit none
    !--- inout: each rank writes its share of a node-shared array
    real(kind=wp), intent(inout) :: gHI(:), gHeI(:), gHeII(:)
    real(kind=wp), intent(inout) :: hHI(:), hHeI(:), hHeII(:)
    real(kind=wp) :: fac
    integer :: il

    if (.not. allocated(hhe_rate_path)) &
       error stop 'H/He path estimators are not allocated'
    !--- node-shared output: the node's ranks take disjoint leaves
    do il = shared_host_rank() + 1, score_nleaf, shared_host_size()
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


  subroutine ion_score_apply_sec(edep)
    !--- Photoelectron energy deposition rate above E_SECONDARY [erg/s per
    !--- absorber particle], edep(k, a, leaf) for energy node k and absorber
    !--- a = 1 (H I), 2 (He I), 3 (He II): sec_ion_path with the factor that
    !--- turns Sum(path_lum coeff) into Int 4 pi J_nu coeff dnu, as for the
    !--- photoheating.  Summed over k it is the part of heat_* whose
    !--- photoelectrons exceed E_SECONDARY.
    use octree_mod, only : leaf_half
    implicit none
    real(kind=wp), intent(inout) :: edep(:,:,:)   ! node-shared, this rank's share
    real(kind=wp) :: fac
    integer :: il

    if (.not. allocated(sec_ion_path)) &
       error stop 'secondary-ionization path estimators are not allocated'
    !--- node-shared output: the node's ranks take disjoint leaves
    do il = shared_host_rank() + 1, score_nleaf, shared_host_size()
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       edep(:,:,il) = reshape(sec_ion_path(:,il)*fac, [N_EDEP, 3])
    end do
  end subroutine ion_score_apply_sec


  subroutine ion_score_report_metal_rates()
    if (.not. par%use_metals) return
    if (mpar%p_rank == 0) &
       write(*,'(a)') ' ION: continuous metal direct rates: ACTIVE'
  end subroutine ion_score_report_metal_rates

end module ion_score_mod
