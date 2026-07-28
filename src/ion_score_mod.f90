module ion_score_mod
!---------------------------------------------------------------------------
! Common ionizing-packet path scorer.
!
! Every direct analytic segment and every sampled scattered segment supplies
! the same path quantity
!
!   path_lum = Lpacket * weight * effective_path_length
!
! [erg/s * code length].  During the grouped-energy migration this routine
! continues to fill the authoritative jt_ion(inu,leaf) tally and, in parallel,
! accumulates direct H/He ionization/heating estimators at ion_e(inu).
! The shadow arrays are validation-only: gas_rates_mod still derives all
! solver rates from jt_ion.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: ion_score_setup, ion_score_resize, ion_score_reset
  public :: score_ion_path, ion_score_reduce, ion_score_compare_hhe
  public :: ion_score_compare_metals

  integer, parameter :: I_GHI = 1, I_GHEI = 2, I_GHEII = 3
  integer, parameter :: I_HHI = 4, I_HHEI = 5, I_HHEII = 6
  integer, parameter :: NHHE = 6
  real(kind=wp), parameter :: SHADOW_RTOL = 5.0e-12_wp

  real(kind=wp), allocatable :: hhe_coeff(:,:)   ! (NHHE,nnu_band)
  real(kind=wp), allocatable :: hhe_shadow(:,:)  ! (NHHE,nleaf), rank local then reduced
  integer :: score_nleaf = 0
  integer(kind=int64) :: direct_segments = 0_int64
  integer(kind=int64) :: scattered_segments = 0_int64

contains

  subroutine ion_score_setup(nleaf)
    use ion_band_mod, only : ion_e, nnu_band
    use photo_xsec,   only : sigma_HI, sigma_HeI, sigma_HeII
    implicit none
    integer, intent(in) :: nleaf
    real(kind=wp) :: energy, sHI, sHeI, sHeII
    integer :: inu

    if (allocated(hhe_coeff)) deallocate(hhe_coeff)
    if (allocated(hhe_shadow)) deallocate(hhe_shadow)
    score_nleaf = nleaf
    if (.not. par%ion_shadow_rates) return
    allocate(hhe_coeff(NHHE,nnu_band), hhe_shadow(NHHE,nleaf))
    hhe_shadow = 0.0_wp

    do inu = 1, nnu_band
       energy = ion_e(inu)
       sHI    = sigma_HI(energy)
       sHeI   = sigma_HeI(energy)
       sHeII  = sigma_HeII(energy)
       hhe_coeff(I_GHI,inu)   = sHI/(energy*ev2erg)
       hhe_coeff(I_GHEI,inu)  = sHeI/(energy*ev2erg)
       hhe_coeff(I_GHEII,inu) = sHeII/(energy*ev2erg)
       hhe_coeff(I_HHI,inu)   = sHI*(1.0_wp - eth_HI/energy)
       hhe_coeff(I_HHEI,inu)  = sHeI*(1.0_wp - eth_HeI/energy)
       hhe_coeff(I_HHEII,inu) = sHeII*(1.0_wp - eth_HeII/energy)
    end do
    if (par%use_metals) then
       block
         use species_mod, only : species_shadow_setup
         call species_shadow_setup(nleaf)
       end block
    end if
  end subroutine ion_score_setup

  subroutine ion_score_resize(nleaf)
    integer, intent(in) :: nleaf
    call ion_score_setup(nleaf)
  end subroutine ion_score_resize

  subroutine ion_score_reset()
    use ion_packet_mod, only : packet_opacity_shadow_reset
    call packet_opacity_shadow_reset()
    if (.not. par%ion_shadow_rates) return
    if (allocated(hhe_shadow)) hhe_shadow = 0.0_wp
    direct_segments = 0_int64
    scattered_segments = 0_int64
    if (par%use_metals) then
       block
         use species_mod, only : species_shadow_reset
         call species_shadow_reset()
       end block
    end if
  end subroutine ion_score_reset

  subroutine score_ion_path(photon, il, path_lum)
    use jtally_mod, only : jt_ion
    implicit none
    type(photon_type), intent(in) :: photon
    integer,           intent(in) :: il
    real(kind=wp),     intent(in) :: path_lum
    integer :: inu

    if (path_lum == 0.0_wp) return
    inu = photon%inu
    jt_ion(inu,il) = jt_ion(inu,il) + path_lum
    if (.not. par%ion_shadow_rates) return
    if (photon%nscatt > 0) then
       scattered_segments = scattered_segments + 1_int64
    else
       direct_segments = direct_segments + 1_int64
    end if
    hhe_shadow(:,il) = hhe_shadow(:,il) + path_lum*hhe_coeff(:,inu)
    if (par%use_metals) then
       block
         use species_mod, only : species_shadow_score
         call species_shadow_score(inu, il, path_lum)
       end block
    end if
  end subroutine score_ion_path

  subroutine ion_score_reduce()
    use mpi
    use ion_packet_mod, only : packet_opacity_shadow_reduce
    implicit none
    integer :: ierr, nchunk, i0, n
    integer(kind=int64) :: counts(2)

    call packet_opacity_shadow_reduce()
    if (.not. par%ion_shadow_rates) return
    nchunk = max(1, 100000000/NHHE)
    i0 = 1
    do while (i0 <= score_nleaf)
       n = min(nchunk, score_nleaf-i0+1)
       call MPI_ALLREDUCE(MPI_IN_PLACE, hhe_shadow(:,i0:i0+n-1), NHHE*n, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       i0 = i0 + n
    end do
    counts = [direct_segments, scattered_segments]
    call MPI_ALLREDUCE(MPI_IN_PLACE, counts, 2, MPI_INTEGER8, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    direct_segments = counts(1)
    scattered_segments = counts(2)
    if (par%use_metals) then
       block
         use species_mod, only : species_shadow_reduce
         call species_shadow_reduce()
       end block
    end if
  end subroutine ion_score_reduce

  subroutine ion_score_compare_hhe(gHI, gHeI, gHeII, hHI, hHeI, hHeII)
    use octree_mod, only : leaf_half
    use mpi
    use ion_packet_mod, only : packet_opacity_shadow_report
    implicit none
    real(kind=wp), intent(in) :: gHI(:), gHeI(:), gHeII(:)
    real(kind=wp), intent(in) :: hHI(:), hHeI(:), hHeII(:)
    real(kind=wp) :: max_abs(NHHE), scale(NHHE), rel(NHHE)
    real(kind=wp) :: reference(NHHE), shadow, fac
    integer :: il, i, ierr
    logical :: passed

    if (.not. par%ion_shadow_rates) return
    call packet_opacity_shadow_report(SHADOW_RTOL)
    max_abs = 0.0_wp
    scale = 0.0_wp
    do il = 1, score_nleaf
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       reference = [gHI(il), gHeI(il), gHeII(il), hHI(il), hHeI(il), hHeII(il)]
       do i = 1, NHHE
          shadow = hhe_shadow(i,il)*fac
          max_abs(i) = max(max_abs(i), abs(shadow-reference(i)))
          scale(i) = max(scale(i), abs(reference(i)))
       end do
    end do
    do i = 1, NHHE
       if (scale(i) > 0.0_wp) then
          rel(i) = max_abs(i)/scale(i)
       else
          rel(i) = max_abs(i)
       end if
    end do
    passed = maxval(rel) <= SHADOW_RTOL
    if (mpar%p_rank == 0) then
       write(*,'(a,2(i0,1x))') ' ION: scored path segments direct/scattered = ', &
                               direct_segments, scattered_segments
       write(*,'(a,6es11.3)') ' ION: H/He shadow relative errors: ', rel
       write(*,'(a,a)') ' ION: grouped H/He direct-rate shadow: ', &
                        merge('PASS', 'FAIL', passed)
    end if
    if (.not. passed) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end subroutine ion_score_compare_hhe

  subroutine ion_score_compare_metals()
    if (.not. par%ion_shadow_rates .or. .not. par%use_metals) return
    block
      use species_mod, only : species_shadow_compare
      call species_shadow_compare(SHADOW_RTOL)
    end block
  end subroutine ion_score_compare_metals

end module ion_score_mod
