module dust_ir_band_mod
!---------------------------------------------------------------------------
! MoCHII: tallies of the dust thermal re-emission (infrared) band.
!
! Grains reradiate the power they absorb at h nu ~ 1e-3 - 1 eV, far below
! every gas photoionization threshold (the lowest in the registry is
! Ca I at 6.11 eV, H I at 13.598 eV).  A packet of that band therefore has
! sigma_HI = sigma_HeI = sigma_HeII = 0 and no metal cross section at all,
! and its ONLY interaction is with grains.  The transport machinery is shared
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
!   the emitted / extinguished / escaping luminosities, split per packet by
!            the edge transmission exp(-tau_edge) exactly as the ionizing
!            band's own budget;
!   the escaping spectrum on the emission wavelength grid, which is the
!            observable the transported infrared produces.
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
  public :: dust_ir_heating, dust_ir_energy_budget
  public :: dust_ir_escaped_spectrum

  integer, parameter :: ENERGY_ACC_WP = selected_real_kind(18,300)

  real(kind=wp), allocatable :: ir_abs_path(:)      ! (nleaf)
  integer :: ir_nleaf = 0

  real(kind=ENERGY_ACC_WP) :: ir_emitted_acc  = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ir_extinct_acc  = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ir_escaped_acc  = 0.0_ENERGY_ACC_WP
  real(kind=wp) :: ir_emitted = 0.0_wp
  real(kind=wp) :: ir_extinct = 0.0_wp
  real(kind=wp) :: ir_escaped = 0.0_wp

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
    ir_escaped_acc = 0.0_ENERGY_ACC_WP
    ir_emitted = 0.0_wp;  ir_extinct = 0.0_wp;  ir_escaped = 0.0_wp
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
  ! Energy ledger of one launched infrared packet, split by its edge
  ! transmission exp(-tau_edge): emitted = extinguished + escaping, per
  ! packet by construction.  The escaping part also enters the spectrum
  ! histogram at the packet's own wavelength.
  !
  ! With par%ion_dust_scatter the edge optical depth is the EXTINCTION one,
  ! so the "extinguished" column then also holds the scattered-out energy and
  ! the histogram is the unscattered escape.  The grain albedo in the thermal
  ! infrared is ~1e-3 or less (a 100 um photon has albedo < 1e-4 in the
  ! Draine 2003 mixture), so the two definitions differ by that much; the
  ! grain HEATING never uses this ledger, it comes from the path estimator
  ! above, which is exact with or without scattering.
  !=========================================================================
  subroutine dust_ir_escape_packet(photon, tau_edge)
    implicit none
    type(photon_type), intent(in) :: photon
    real(kind=wp),     intent(in) :: tau_edge
    real(kind=wp) :: luminosity, transmission, lesc
    integer :: lo, hi, mid

    luminosity   = photon%Lpacket*photon%wgt
    transmission = exp(-max(tau_edge, 0.0_wp))
    ir_emitted_acc = ir_emitted_acc + real(luminosity, ENERGY_ACC_WP)
    ir_extinct_acc = ir_extinct_acc &
       + real(luminosity*(1.0_wp - transmission), ENERGY_ACC_WP)
    lesc = luminosity*transmission
    ir_escaped_acc = ir_escaped_acc + real(lesc, ENERGY_ACC_WP)

    if (ir_nlam <= 0 .or. lesc <= 0.0_wp) return
    if (photon%lambda <= lam_edge(1) .or. photon%lambda >= lam_edge(ir_nlam+1)) return
    lo = 1;  hi = ir_nlam
    do while (lo < hi)                       ! first bin with lam < lam_edge(k+1)
       mid = (lo + hi)/2
       if (photon%lambda < lam_edge(mid+1)) then
          hi = mid
       else
          lo = mid + 1
       end if
    end do
    Lesc_bin(lo) = Lesc_bin(lo) + lesc
  end subroutine dust_ir_escape_packet

  !=========================================================================
  subroutine dust_ir_band_reduce()
    use mpi
    implicit none
    integer :: ierr, nchunk, i0, n
    real(kind=wp) :: totals(3)

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
    totals = real([ir_emitted_acc, ir_extinct_acc, ir_escaped_acc], wp)
    call MPI_ALLREDUCE(MPI_IN_PLACE, totals, 3, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    ir_emitted = totals(1);  ir_extinct = totals(2);  ir_escaped = totals(3)
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
  subroutine dust_ir_energy_budget(emitted, extinguished, escaped)
    implicit none
    real(kind=wp), intent(out) :: emitted, extinguished, escaped
    emitted = ir_emitted;  extinguished = ir_extinct;  escaped = ir_escaped
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
