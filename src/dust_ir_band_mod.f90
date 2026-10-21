module dust_ir_band_mod
!---------------------------------------------------------------------------
! MoCHII: tallies of the dust thermal re-emission (infrared) band.
!
! Grains reradiate the power they absorb at h nu ~ 1e-3 - 1 eV, far below
! every gas photoionization threshold (the lowest in the registry is
! Ca I at 6.11 eV, H I at 13.598 eV).  A packet of that band therefore has
! sigma_HI = sigma_HeI = sigma_HeII = 0 and no metal cross section at all,
! and its ONLY interaction is with grains.  The transport code is shared
! with the ionizing band (transport_ion_packet evaluates every cross section
! at photon%energy_eV, so the vanishing gas terms drop out by themselves),
! but the ESTIMATORS are not: an infrared segment carries no photoionization,
! no photoheating, no Habing G0 and no ionizing-band mean intensity.
! transported_band tells the common path scorer which set of estimators a
! segment belongs to, and this module owns the infrared set:
!
!   ir_abs_path(leaf) = Sum(path_lum * sigma_abs,dust(E))
!         -> dust_ir_heating: the grain heating rate the reabsorbed infrared
!            adds to the primary (EUV/FUV) one, closing the dust self-heating
!            loop;
!   the direct split of every launched packet by its edge transmission
!            exp(-tau_edge), emitted = extinguished + transmitted, an
!            identity checked to rounding exactly as in the ionizing band;
!   the physical totals: escaped = the direct transmission plus the
!            scattered flights that leave the grid, lost = the weight of
!            walks that ended without leaving the nebula (normally zero),
!            absorbed = emitted - escaped - lost;
!   the escaping spectrum on the emission wavelength grid, direct and
!            scattered escapes together, which is the observable the
!            transported infrared produces.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  !--- which band the packets now in flight belong to.  BAND_IONIZING is the
  !--- state everywhere except inside a dust-emission transport pass.
  integer, parameter, public :: BAND_IONIZING = 1
  integer, parameter, public :: BAND_DUST_IR  = 2
  integer,            public :: transported_band = BAND_IONIZING

  public :: dust_ir_band_setup, dust_ir_band_reset, dust_ir_band_reduce
  public :: dust_ir_absorb_path, dust_ir_escape_packet
  public :: dust_ir_scattered_escape, dust_ir_lost_weight
  public :: dust_ir_heating, dust_ir_direct_split, dust_ir_energy_budget
  public :: dust_ir_escaped_spectrum

  integer, parameter :: ENERGY_ACC_WP = selected_real_kind(18,300)

  real(kind=wp), allocatable :: ir_abs_path(:)      ! (nleaf)
  integer :: ir_nleaf = 0

  !--- energy totals [erg/s], rank local in the extended kind, then reduced;
  !--- the same definitions as the ionizing band's (ion_score_mod): the
  !--- direct split emitted = extinct + transmitted, the scattered escape,
  !--- and the weight lost from the edge walk (out of the transmitted part)
  !--- and from the scattered flights (out of the extinguished part).
  real(kind=ENERGY_ACC_WP) :: ir_emitted_acc  = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ir_extinct_acc  = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ir_transmitted_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ir_scattered_escaped_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ir_lost_direct_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ir_lost_scattered_acc = 0.0_ENERGY_ACC_WP
  real(kind=wp) :: ir_emitted = 0.0_wp
  real(kind=wp) :: ir_extinct = 0.0_wp
  real(kind=wp) :: ir_transmitted = 0.0_wp
  real(kind=wp) :: ir_scattered_escaped = 0.0_wp
  real(kind=wp) :: ir_lost_direct = 0.0_wp
  real(kind=wp) :: ir_lost_scattered = 0.0_wp
  integer(kind=int64) :: ir_lost_flights = 0_int64

  !--- escaping spectrum: luminosity per node of the emission wavelength grid,
  !--- with geometric-midpoint bin edges (lam_edge) so the histogram divides
  !--- back to L_lambda on the same nodes the emission spectrum uses.
  real(kind=wp), allocatable :: lam_esc(:), lam_edge(:), Lesc_bin(:)
  integer :: ir_nlam = 0

contains

  !=========================================================================
  subroutine dust_ir_band_setup(nleaf, lam)
    implicit none
    integer,       intent(in) :: nleaf
    real(kind=wp), intent(in) :: lam(:)      ! emission grid [um], ascending
    integer :: n, k

    n = size(lam)
    if (allocated(ir_abs_path)) deallocate(ir_abs_path)
    allocate(ir_abs_path(nleaf))
    ir_nleaf = nleaf

    if (allocated(lam_esc)) deallocate(lam_esc, lam_edge, Lesc_bin)
    allocate(lam_esc(n), lam_edge(n+1), Lesc_bin(n))
    ir_nlam = n
    lam_esc = lam
    do k = 2, n
       lam_edge(k) = sqrt(lam(k-1)*lam(k))
    end do
    lam_edge(1)   = lam(1)*sqrt(lam(1)/lam(2))
    lam_edge(n+1) = lam(n)*sqrt(lam(n)/lam(n-1))
    call dust_ir_band_reset()
  end subroutine dust_ir_band_setup

  !=========================================================================
  subroutine dust_ir_band_reset()
    implicit none
    if (allocated(ir_abs_path)) ir_abs_path = 0.0_wp
    if (allocated(Lesc_bin))    Lesc_bin    = 0.0_wp
    ir_emitted_acc = 0.0_ENERGY_ACC_WP
    ir_extinct_acc = 0.0_ENERGY_ACC_WP
    ir_transmitted_acc = 0.0_ENERGY_ACC_WP
    ir_scattered_escaped_acc = 0.0_ENERGY_ACC_WP
    ir_lost_direct_acc = 0.0_ENERGY_ACC_WP
    ir_lost_scattered_acc = 0.0_ENERGY_ACC_WP
    ir_emitted = 0.0_wp;  ir_extinct = 0.0_wp;  ir_transmitted = 0.0_wp
    ir_scattered_escaped = 0.0_wp
    ir_lost_direct = 0.0_wp;  ir_lost_scattered = 0.0_wp
    ir_lost_flights = 0_int64
  end subroutine dust_ir_band_reset

  !=========================================================================
  ! One path segment's grain absorption, path_lum * sigma_abs,dust(E)
  ! [erg/s * code length].  Called from the common path scorer while
  ! transported_band = BAND_DUST_IR.
  !=========================================================================
  subroutine dust_ir_absorb_path(il, abs_lum)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: abs_lum
    ir_abs_path(il) = ir_abs_path(il) + abs_lum
  end subroutine dust_ir_absorb_path

  !=========================================================================
  ! Direct split of one launched infrared packet by its edge transmission
  ! exp(-tau_edge): emitted = extinguished + transmitted, per packet by
  ! construction.  The transmitted part escapes and enters the spectrum
  ! histogram at the packet's own wavelength, unless the edge walk was lost
  ! (edge_lost: MAX_WALK_STEPS or a failed xy-periodic re-entry), when it is
  ! booked as lost instead; the path estimates that walk already scored stay.
  !
  ! With par%ion_dust_scatter the edge optical depth is the EXTINCTION one,
  ! so "extinguished" then also holds the scattered-out energy; what of it
  ! leaves the grid comes back through dust_ir_scattered_escape, so the
  ! escaped total and the histogram hold both escapes.
  !=========================================================================
  subroutine dust_ir_escape_packet(photon, tau_edge, edge_lost)
    implicit none
    type(photon_type), intent(in) :: photon
    real(kind=wp),     intent(in) :: tau_edge
    logical,           intent(in) :: edge_lost
    real(kind=wp) :: luminosity, transmission, lesc

    luminosity   = photon%Lpacket*photon%wgt
    transmission = exp(-max(tau_edge, 0.0_wp))
    ir_emitted_acc = ir_emitted_acc + real(luminosity, ENERGY_ACC_WP)
    ir_extinct_acc = ir_extinct_acc &
       + real(luminosity*(1.0_wp - transmission), ENERGY_ACC_WP)
    lesc = luminosity*transmission
    ir_transmitted_acc = ir_transmitted_acc + real(lesc, ENERGY_ACC_WP)
    if (edge_lost) then
       ir_lost_direct_acc = ir_lost_direct_acc + real(lesc, ENERGY_ACC_WP)
       ir_lost_flights = ir_lost_flights + 1_int64
       return
    end if
    call escaping_spectrum_add(photon%lambda, lesc)
  end subroutine dust_ir_escape_packet

  !=========================================================================
  ! A dust-scattered infrared flight left the grid through an open face
  ! (called through ion_score_mod::ion_score_scattered_escape): its
  ! remaining weight escapes, into the escaped total and the spectrum.
  !=========================================================================
  subroutine dust_ir_scattered_escape(photon)
    implicit none
    type(photon_type), intent(in) :: photon
    real(kind=wp) :: lesc
    lesc = photon%Lpacket*photon%wgt
    ir_scattered_escaped_acc = ir_scattered_escaped_acc &
       + real(lesc, ENERGY_ACC_WP)
    call escaping_spectrum_add(photon%lambda, lesc)
  end subroutine dust_ir_scattered_escape

  !=========================================================================
  ! The remaining weight [erg/s] of a scattered infrared flight that left
  ! the calculation without leaving the nebula (ion_score_lost_weight).
  !=========================================================================
  subroutine dust_ir_lost_weight(luminosity)
    implicit none
    real(kind=wp), intent(in) :: luminosity
    ir_lost_scattered_acc = ir_lost_scattered_acc &
       + real(luminosity, ENERGY_ACC_WP)
    ir_lost_flights = ir_lost_flights + 1_int64
  end subroutine dust_ir_lost_weight

  !=========================================================================
  ! Escaping luminosity lesc [erg/s] at wavelength lam [um] into the
  ! histogram on the emission grid.  Every packet wavelength lies inside
  ! [lam(1), lam(n)], within the outer bin edges, so the histogram sums to
  ! the escaped total.
  !=========================================================================
  subroutine escaping_spectrum_add(lam, lesc)
    implicit none
    real(kind=wp), intent(in) :: lam, lesc
    integer :: lo, hi, mid
    if (ir_nlam <= 0 .or. lesc <= 0.0_wp) return
    if (lam <= lam_edge(1) .or. lam >= lam_edge(ir_nlam+1)) return
    lo = 1;  hi = ir_nlam
    do while (lo < hi)                       ! first bin with lam < lam_edge(k+1)
       mid = (lo + hi)/2
       if (lam < lam_edge(mid+1)) then
          hi = mid
       else
          lo = mid + 1
       end if
    end do
    Lesc_bin(lo) = Lesc_bin(lo) + lesc
  end subroutine escaping_spectrum_add

  !=========================================================================
  subroutine dust_ir_band_reduce()
    use mpi
    implicit none
    integer :: ierr, nchunk, i0, n
    real(kind=wp) :: totals(6)

    nchunk = 100000000
    i0 = 1
    do while (i0 <= ir_nleaf)
       n = min(nchunk, ir_nleaf-i0+1)
       call MPI_ALLREDUCE(MPI_IN_PLACE, ir_abs_path(i0:i0+n-1), n, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       i0 = i0 + n
    end do
    if (ir_nlam > 0) &
       call MPI_ALLREDUCE(MPI_IN_PLACE, Lesc_bin, ir_nlam, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    totals = real([ir_emitted_acc, ir_extinct_acc, ir_transmitted_acc, &
                   ir_scattered_escaped_acc, ir_lost_direct_acc, &
                   ir_lost_scattered_acc], wp)
    call MPI_ALLREDUCE(MPI_IN_PLACE, totals, 6, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    ir_emitted = totals(1);  ir_extinct = totals(2);  ir_transmitted = totals(3)
    ir_scattered_escaped = totals(4)
    ir_lost_direct = totals(5);  ir_lost_scattered = totals(6)
    call MPI_ALLREDUCE(MPI_IN_PLACE, ir_lost_flights, 1, MPI_INTEGER8, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine dust_ir_band_reduce

  !=========================================================================
  ! Grain heating [erg/s/cm^3] by the reabsorbed infrared.  Same transport
  ! form as the EUV/FUV grain heating in ion_score_mod: ir_abs_path times the
  ! factor below is 4 pi J_abs, and rhokap/distance2cm is the grain
  ! extinction per cm at the reference wavelength that sigma_abs,dust scales.
  !=========================================================================
  subroutine dust_ir_heating(heat_ir)
    use octree_mod, only : leaf_half, amr_grid
    implicit none
    real(kind=wp), intent(out) :: heat_ir(:)
    real(kind=wp) :: fac
    integer :: il
    do il = 1, ir_nleaf
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       heat_ir(il) = ir_abs_path(il)*fac*amr_grid%rhokap(il)/par%distance2cm
    end do
  end subroutine dust_ir_heating

  !=========================================================================
  ! The direct split of the reduced pass [erg/s]: emitted = extinguished +
  ! transmitted to rounding (dust_emis_transport_mod aborts otherwise).
  subroutine dust_ir_direct_split(emitted, extinguished, transmitted)
    implicit none
    real(kind=wp), intent(out) :: emitted, extinguished, transmitted
    emitted = ir_emitted;  extinguished = ir_extinct
    transmitted = ir_transmitted
  end subroutine dust_ir_direct_split

  !=========================================================================
  ! Physical totals of the reduced pass [erg/s]: escaped (direct and
  ! scattered, the integral of the escaping spectrum), lost (normally zero;
  ! nlost walks), absorbed = emitted - escaped - lost, formed from the
  ! extinguished part so that a run without scattering and losses returns
  ! the direct split itself.
  !=========================================================================
  subroutine dust_ir_energy_budget(emitted, absorbed, escaped, lost, nlost)
    implicit none
    real(kind=wp), intent(out) :: emitted, absorbed, escaped, lost
    integer(kind=int64), intent(out) :: nlost
    emitted  = ir_emitted
    escaped  = ir_transmitted - ir_lost_direct + ir_scattered_escaped
    lost     = ir_lost_direct + ir_lost_scattered
    absorbed = ir_extinct - ir_scattered_escaped - ir_lost_scattered
    nlost    = ir_lost_flights
  end subroutine dust_ir_energy_budget

  !=========================================================================
  ! Escaping spectrum L_lambda [erg/s/um] on the emission grid.
  !=========================================================================
  subroutine dust_ir_escaped_spectrum(lam, Llam, dlam)
    implicit none
    real(kind=wp), allocatable, intent(out) :: lam(:), Llam(:), dlam(:)
    integer :: k
    allocate(lam(ir_nlam), Llam(ir_nlam), dlam(ir_nlam))
    lam = lam_esc
    do k = 1, ir_nlam
       dlam(k) = lam_edge(k+1) - lam_edge(k)
       Llam(k) = Lesc_bin(k)/dlam(k)
    end do
  end subroutine dust_ir_escaped_spectrum

end module dust_ir_band_mod
