module ion_band_mod
!---------------------------------------------------------------------------
! MoCHII: ionizing frequency band + ionizing-photon generation.
!
! The ionizing band covers photon energies [par%eion_min, par%eion_max] eV
! with par%nnu_ion log-spaced bins (par%nnu_ion is the immutable ionizing
! count; 10-20 bins reach percent-level rate integrals).  With par%add_fuv
! extra FUV bins are prepended and nnu_band is the total consumed
! downstream.  The SED grid covers dust wavelengths and is untouched; the
! ionizing band carries its own grids, source sampling, opacity, and J tally.
!
! Source spectrum in the band is resolved by component.  A single point
! source (and the global default) uses par%ion_spectrum (2-column file:
! E [eV], L_E) or, when empty, a Planck function B_nu(par%tstar).
! With several internal point sources, par%src_spectrum_file (one multi-column
! file: E then one L_E per source) column i or par%src_tstar(i) gives source i
! its own spectrum; the external field uses par%ext_spectrum / par%ext_tstar
! (see comp_shape / point_shape).
!
! par%spectrum_type fixes the column units of every file slot: 'shape'
! (arbitrary, renormalized to the scale = legacy) or a PHYSICAL type ('per_ev',
! 'per_hz', 'per_ang', 'per_um').  A physical file is ABSOLUTE: bin luminosities
! come from the file integral directly (rescaled only when the scale is set);
! an unset scale DERIVES the luminosity from the file.  par%ext_spectrum may
! also name an analytic ISRF preset ('draine'/'habing'/'mathis', FUV-only,
! absolute, needs add_fuv; see preset_je).
!
! Bin luminosities are integrated on a 32-point sub-grid per bin.  'shape'
! spectra are normalized so the ionizing segment carries the scale
! (par%luminosity / src_lum(i); external pi*J*A).  Packets sample bins from the
! luminosity CDF and carry Lpacket = ion_Ltot / nphotons.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: ion_setup, gen_ion_photon, gen_ion_photon_qmc, ion_ext_preset_id
  public :: ion_e, ion_de, ion_nu, ion_dnu, ion_lum, ion_Ltot
  public :: ion_eedge, ion_bin_of, ion_qmc_ndim
  public :: nnu_band, nfuv_band

  real(kind=wp), allocatable :: ion_e(:)     ! bin center [eV]
  real(kind=wp), allocatable :: ion_de(:)    ! bin width  [eV]
  real(kind=wp), allocatable :: ion_nu(:)    ! bin center [Hz]
  real(kind=wp), allocatable :: ion_dnu(:)   ! bin width  [Hz]
  real(kind=wp), allocatable :: ion_lum(:)   ! total bin luminosity (sum over components) [erg/s]
  real(kind=wp), allocatable :: ion_cdf(:)   ! sampling CDF over bins (single-component fast path)
  !--- band bin EDGES [eV], size nnu_band+1, ascending.  The single source of
  !--- truth for mapping a photon energy to a bin: filled at ion_setup for BOTH
  !--- the aligned and the legacy-log grids (ion_bin_of does a binary search of
  !--- it).  Ascending because energies increase with bin index.
  real(kind=wp), allocatable, protected :: ion_eedge(:)

  !--- Source-component model.  The band is fed by par%nsource internal point
  !--- sources plus, independently, one external field (ON when
  !--- par%ext_intensity>0).  With more than one component the packets are split
  !--- among the components in proportion to their band luminosity (src_cdf),
  !--- and the frequency bin is drawn from that component's own spectrum
  !--- (ion_cdf_src).  With exactly one component the single fast path is taken
  !--- (ion_cdf over the global spectrum), preserving the legacy RNG stream.
  integer :: ncomp = 1                          ! number of source components
  logical :: multi_src = .false.                ! .true. when ncomp >= 2
  integer,       allocatable :: comp_kind(:)    ! 0 = external, i>=1 = point source i
  real(kind=wp), allocatable :: src_cdf(:)      ! component-selection CDF (band totals)
  real(kind=wp), allocatable :: ion_cdf_src(:,:)! (nnu_band, ncomp) each component's bin CDF
  !--- total band luminosity carried by the packets: par%luminosity
  !--- ([eion_min, eion_max]) plus, with par%add_fuv, the FUV part of the
  !--- same source spectrum.
  real(kind=wp), protected :: ion_Ltot = 0.0_wp

  !--- total number of band bins consumed downstream: nnu_band = the
  !--- ionizing bins (par%nnu_ion) plus, with par%add_fuv, nfuv_band FUV
  !--- bins prepended at the low-energy end.  par%nnu_ion stays the input
  !--- ionizing bin count (immutable); every band array and loop uses nnu_band.
  integer, protected :: nnu_band  = 0   ! total band bins (ionizing + FUV)
  integer, protected :: nfuv_band = 0   ! FUV bins at the low-energy end

  !--- Habing (1968) FUV energy density in 6-13.6 eV [erg/cm^3]; G0 = 1 unit.
  real(kind=wp), parameter :: HABING_U_FUV = 5.29e-14_wp

contains

  !=========================================================================
  subroutine ion_setup(metal_eth, nmetal)
    use mpi
    use octree_mod, only : amr_grid
    implicit none
    !--- metal_eth(1:nmetal): active metal photoionization thresholds [eV],
    !--- gathered by the caller from the species registry.  Only consulted
    !--- when par%ion_align_edges is on; absent for the H/He-only path.
    real(kind=wp), intent(in), optional :: metal_eth(:)
    integer,       intent(in), optional :: nmetal
    integer, parameter :: NSUB = 32
    real(kind=wp), allocatable :: eedge(:), meth(:)
    real(kind=wp) :: lo, hi
    real(kind=wp) :: Lnorm, Jband, uFUV
    integer :: nnu, nfuv, i, ierr, nmet, presid
    logical :: is_abs, derived

    !--- FUV extension: nnu_fuv extra log bins on [efuv_min, eion_min]
    !--- BELOW the par%nnu_ion ionizing bins.  par%nnu_ion stays the
    !--- immutable ionizing count; nnu_band = par%nnu_ion + nfuv is the
    !--- total bin count consumed by every band array.
    nfuv = 0
    if (par%add_fuv) then
       if (par%efuv_min >= par%eion_min) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: add_fuv requires par%efuv_min < par%eion_min.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       nfuv = par%nnu_fuv
    end if
    nnu = par%nnu_ion + nfuv     ! total band bins (par%nnu_ion is the immutable ionizing count)
    nnu_band  = nnu
    nfuv_band = nfuv
    allocate(eedge(nnu+1), ion_e(nnu), ion_de(nnu), ion_nu(nnu), &
             ion_dnu(nnu), ion_lum(nnu), ion_cdf(nnu))

    lo = log(par%eion_min);  hi = log(par%eion_max)
    if (nfuv > 0) then
       do i = 1, nfuv                   ! FUV segment [efuv_min, eion_min)
          eedge(i) = exp(log(par%efuv_min) + (lo - log(par%efuv_min)) &
                     *real(i-1,wp)/real(nfuv,wp))
       end do
    end if
    if (par%ion_align_edges) then
       !--- ionizing segment [eion_min, eion_max]: bin edges pinned to the
       !--- ionization thresholds so no bin straddles one.
       nmet = 0
       if (present(nmetal))   nmet = nmetal
       if (nmet > 0 .and. present(metal_eth)) then
          allocate(meth(nmet));  meth = metal_eth(1:nmet)
       else
          allocate(meth(1));  meth = 0.0_wp;  nmet = 0
       end if
       call build_aligned_edges(eedge, nfuv, nnu, meth, nmet)
       deallocate(meth)
    else
       do i = nfuv+1, nnu+1             ! ionizing segment [eion_min, eion_max]
          eedge(i) = exp(lo + (hi - lo)*real(i-1-nfuv,wp)/real(nnu-nfuv,wp))
       end do
    end if
    do i = 1, nnu
       ion_e(i)   = sqrt(eedge(i)*eedge(i+1))          ! geometric bin center
       ion_de(i)  = eedge(i+1) - eedge(i)
       ion_nu(i)  = ion_e(i)  * ev2erg / h_planck_cgs
       ion_dnu(i) = ion_de(i) * ev2erg / h_planck_cgs
    end do

    !--- band bin luminosities and sampling CDFs.  Determine the source
    !--- components: par%nsource internal point sources plus one external field
    !--- (ON when par%ext_intensity>0, an ISRF preset is named, or a physical-
    !--- type par%ext_spectrum file is given).  With exactly one component the
    !--- single fast path runs (legacy RNG stream); otherwise the multi path
    !--- splits the packets among the components (each with its own spectrum).
    ncomp = par%nsource
    if (external_on()) ncomp = ncomp + 1
    multi_src = ncomp >= 2

    Lnorm = 0.0_wp;  Jband = 0.0_wp;  uFUV = 0.0_wp;  presid = 0;  derived = .false.
    if (.not. multi_src) then
       if (par%nsource == 0) then
          !===================== external-only fast path =====================
          !--- absolute (preset / physical file) or legacy (shape / Planck)
          !--- external field; build_ext_bins fills ion_lum with bin
          !--- luminosities [erg/s] directly and returns the diagnostics.
          call build_ext_bins(ion_lum, eedge, nnu, nfuv, NSUB, Lnorm, Jband, uFUV, presid)
       else
          !===================== single point source fast path ================
          !--- src_spectrum_file column 1 (if set) > global ion_spectrum >
          !--- global tstar Planck.  The legacy scalars rule here, so
          !--- src_tstar(1) is NOT consulted (use par%tstar with one source).
          call point_shape(ion_lum, eedge, nnu, NSUB, 1, .false., is_abs)
          call finalize_source(ion_lum, nnu, nfuv, is_abs, par%luminosity, Lnorm)
          derived = is_abs .and. (par%luminosity <= 0.0_wp)
       end if
       ion_Ltot = sum(ion_lum)
       ion_cdf(1) = ion_lum(1)
       do i = 2, nnu
          ion_cdf(i) = ion_cdf(i-1) + ion_lum(i)
       end do
       ion_cdf = ion_cdf / ion_cdf(nnu)

       if (mpar%p_rank == 0) call log_fast_path(nnu, nfuv, Lnorm, Jband, uFUV, presid, derived)
    else
       !======================== multi-component path ========================
       call setup_multi_components(eedge, nnu, nfuv, NSUB)
    end if
    !--- keep the full band edge array (both grid types) as the one source of
    !--- truth for the energy -> bin mapping (ion_bin_of).
    if (allocated(ion_eedge)) deallocate(ion_eedge)
    allocate(ion_eedge(nnu+1))
    ion_eedge = eedge
    deallocate(eedge)
  end subroutine ion_setup

  !=========================================================================
  ! Bin index for a photon energy eph [eV]: the bin whose edges bracket eph
  ! (ion_eedge(i) <= eph < ion_eedge(i+1)), so an energy exactly on a threshold
  ! edge falls in the bin ABOVE that edge.  Binary search of ion_eedge, then
  ! clamped to the ionizing bins [nfuv_band+1, nnu_band].  Used for the diffuse
  ! recombination photons (all ionizing) on BOTH the aligned and the log grid.
  !=========================================================================
  integer function ion_bin_of(eph)
    implicit none
    real(kind=wp), intent(in) :: eph
    integer :: lo, hi, mid
    lo = 1;  hi = nnu_band
    do while (lo < hi)
       mid = (lo + hi + 1)/2
       if (ion_eedge(mid) <= eph) then
          lo = mid
       else
          hi = mid - 1
       end if
    end do
    ion_bin_of = min(max(lo, nfuv_band+1), nnu_band)
  end function ion_bin_of

  !=========================================================================
  ! Number of quasi-random launch dimensions consumed by the current source
  ! configuration: 3 for a single internal point source (the stage-1 layout
  ! u = (frequency, mu, phi)); 7 for the fixed superset layout used by every
  ! other configuration (single external, multiple points, mixed).
  !=========================================================================
  integer function ion_qmc_ndim()
    implicit none
    if ((.not. multi_src) .and. trim(par%source_geometry) == 'point') then
       ion_qmc_ndim = 3
    else
       ion_qmc_ndim = 7
    end if
  end function ion_qmc_ndim

  !=========================================================================
  ! Multi-component band setup.  Builds one bin-luminosity array lum_c per
  ! component (each from its own spectrum), accumulates ion_lum = sum_c lum_c,
  ! and forms the component-selection CDF (src_cdf, over band totals) plus each
  ! component's bin CDF (ion_cdf_src).  Internal point sources are normalized to
  ! src_lum(i) (physical file with src_lum unset = derived from the file); the
  ! external field is built by build_ext_bins (absolute preset/physical file or
  ! legacy shape).  The single global ion_cdf is kept meaningful too.
  !=========================================================================
  subroutine setup_multi_components(eedge, nnu, nfuv, nsub)
    implicit none
    real(kind=wp), intent(in) :: eedge(:)
    integer,       intent(in) :: nnu, nfuv, nsub
    real(kind=wp), allocatable :: lum_c(:), comp_L(:), src_used(:)
    logical,       allocatable :: src_derived(:)
    real(kind=wp) :: Lenter, Jband, uFUV
    logical :: ext_on, is_abs
    integer :: is, ic, presid

    ext_on = external_on()
    if (allocated(comp_kind))   deallocate(comp_kind)
    if (allocated(src_cdf))     deallocate(src_cdf)
    if (allocated(ion_cdf_src)) deallocate(ion_cdf_src)
    allocate(comp_kind(ncomp), comp_L(ncomp), src_cdf(ncomp), &
             ion_cdf_src(nnu, ncomp), lum_c(nnu))
    allocate(src_used(max(par%nsource,1)), src_derived(max(par%nsource,1)))
    src_used = 0.0_wp;  src_derived = .false.
    ion_lum(:) = 0.0_wp
    Lenter = 0.0_wp;  Jband = 0.0_wp;  uFUV = 0.0_wp;  presid = 0

    ic = 0
    !--- internal point sources: each with its own resolved spectrum.
    do is = 1, par%nsource
       ic = ic + 1
       comp_kind(ic) = is
       call point_shape(lum_c, eedge, nnu, nsub, is, .true., is_abs)
       call finalize_source(lum_c, nnu, nfuv, is_abs, par%src_lum(is), src_used(is))
       src_derived(is) = is_abs .and. (par%src_lum(is) <= 0.0_wp)
       comp_L(ic) = sum(lum_c)                     ! band total (incl. FUV)
       ion_lum = ion_lum + lum_c
       call fill_cdf(ion_cdf_src(:,ic), lum_c, nnu)
    end do

    !--- external isotropic field (absolute preset/physical file or legacy).
    if (ext_on) then
       ic = ic + 1
       comp_kind(ic) = 0
       call build_ext_bins(lum_c, eedge, nnu, nfuv, nsub, Lenter, Jband, uFUV, presid)
       comp_L(ic) = sum(lum_c)
       ion_lum = ion_lum + lum_c
       call fill_cdf(ion_cdf_src(:,ic), lum_c, nnu)
    end if

    ion_Ltot = sum(ion_lum)
    !--- global bin CDF (kept meaningful for anything reading ion_lum).
    call fill_cdf(ion_cdf, ion_lum, nnu)
    !--- component-selection CDF over the band totals.
    src_cdf(1) = comp_L(1)
    do ic = 2, ncomp
       src_cdf(ic) = src_cdf(ic-1) + comp_L(ic)
    end do
    src_cdf = src_cdf / src_cdf(ncomp)

    if (mpar%p_rank == 0) then
       call band_log(nnu, nfuv)
       write(*,'(a,i0,a)') ' ION: ', ncomp, ' source components:'
       do is = 1, par%nsource
          if (src_derived(is)) then
             write(*,'(a,i0,a,3(1x,f9.4),a,es12.4)') &
                '   point ', is, ': (x,y,z) =', par%src_x(is), par%src_y(is), &
                par%src_z(is), '  L_ion (derived) =', src_used(is)
          else
             write(*,'(a,i0,a,3(1x,f9.4),a,es12.4)') &
                '   point ', is, ': (x,y,z) =', par%src_x(is), par%src_y(is), &
                par%src_z(is), '  L_ion =', src_used(is)
          end if
       end do
       if (ext_on) then
          if (presid > 0) then
             write(*,'(4a,es12.4)') '   external (', trim(par%ext_geometry), &
                '): ISRF preset ', trim(par%ext_spectrum), '  L_enter =', Lenter
             write(*,'(a,es12.4,a,f9.3)') '     u_FUV [erg/cm^3] =', uFUV, &
                ',  G0 (Habing) =', uFUV/HABING_U_FUV
          else
             write(*,'(3a,es12.4,a,es12.4)') '   external (', trim(par%ext_geometry), &
                '): J =', Jband, '  L_enter =', Lenter
          end if
       end if
       write(*,'(a,es12.4)') ' ION: band total (all components) [erg/s] = ', ion_Ltot
    end if
    deallocate(lum_c, comp_L, src_used, src_derived)
  end subroutine setup_multi_components

  !=========================================================================
  ! Cumulative sampling CDF over bins from a bin-luminosity array (normalized
  ! to end at 1).  Shared by the fast path and the multi-component path.
  !=========================================================================
  subroutine fill_cdf(cdf, lum, nnu)
    implicit none
    real(kind=wp), intent(out) :: cdf(:)
    real(kind=wp), intent(in)  :: lum(:)
    integer,       intent(in)  :: nnu
    integer :: i
    cdf(1) = lum(1)
    do i = 2, nnu
       cdf(i) = cdf(i-1) + lum(i)
    end do
    cdf(1:nnu) = cdf(1:nnu) / cdf(nnu)
  end subroutine fill_cdf

  !=========================================================================
  ! rank-0 band-grid log line (shared by both setup paths).
  !=========================================================================
  subroutine band_log(nnu, nfuv)
    implicit none
    integer, intent(in) :: nnu, nfuv
    if (nfuv > 0) then
       write(*,'(a,i4,a,i4,a,f8.3,a,f8.3,a,f8.3,a)') ' ION: band: ', &
          nnu-nfuv, ' ionizing +', nfuv, ' FUV bins, ', par%efuv_min, &
          ' /', par%eion_min, ' -', par%eion_max, ' eV'
    else
       write(*,'(a,i4,a,f8.3,a,f8.3,a)') ' ION: band: ', nnu, ' bins, ', &
          par%eion_min, ' - ', par%eion_max, ' eV'
    end if
  end subroutine band_log

  !=========================================================================
  ! Resolve a component spectrum (file > component tstar > global file >
  ! global tstar) and fill its band bin array on the NSUB sub-grid.  is_abs
  ! is .true. when the source is a PHYSICAL-type file (arr then holds absolute
  ! bin integrals in the file's density units); .false. for a 'shape' file or
  ! a Planck function (arr is an arbitrary-unit shape, renormalized later).
  !=========================================================================
  subroutine comp_shape(arr, eedge, nnu, nsub, specfile, tstar, is_abs)
    use mpi
    implicit none
    real(kind=wp),    intent(out) :: arr(:)
    real(kind=wp),    intent(in)  :: eedge(:)
    integer,          intent(in)  :: nnu, nsub
    character(len=*), intent(in)  :: specfile
    real(kind=wp),    intent(in)  :: tstar
    logical,          intent(out) :: is_abs
    character(len=256) :: eff_file
    real(kind=wp), allocatable :: etab(:), ftab(:)
    real(kind=wp) :: eff_T, e1, e2, es, fsum
    integer :: i, k, ierr, ntab
    is_abs = .false.
    if (len_trim(specfile) > 0) then
       eff_file = specfile;         eff_T = -1.0_wp
    else if (tstar > 0.0_wp) then
       eff_file = '';               eff_T = tstar
    else if (len_trim(par%ion_spectrum) > 0) then
       eff_file = par%ion_spectrum; eff_T = -1.0_wp
    else
       eff_file = '';               eff_T = par%tstar
    end if
    if (len_trim(eff_file) > 0) then
       if (spec_is_physical()) then
          !--- physical-type file: convert to (E [eV] ascending, X_E per eV)
          !--- and bin -> arr = absolute bin integrals (density units).
          call read_phys_table(eff_file, 1, 1, etab, ftab, ntab)
          call bin_interp(etab, ftab, ntab, eedge, nnu, nsub, arr)
          deallocate(etab, ftab)
          is_abs = .true.
       else
          call bin_lum_from_file(eff_file, eedge, nnu, nsub, arr)
       end if
    else
       if (eff_T <= 0.0_wp) then     ! no file and no positive temperature anywhere
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: use_ion_band needs a spectrum: set par%tstar > 0, '// &
             'par%ion_spectrum, or a component spectrum (src_spectrum_file / '// &
             'src_tstar / ext_spectrum / ext_tstar).'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       do i = 1, nnu
          e1 = eedge(i);  e2 = eedge(i+1)
          fsum = 0.0_wp
          do k = 1, nsub
             es = e1 + (e2 - e1)*(real(k,wp) - 0.5_wp)/real(nsub,wp)
             fsum = fsum + planck_nu(es, eff_T)
          end do
          arr(i) = fsum * (e2 - e1)/real(nsub,wp)
       end do
    end if
  end subroutine comp_shape

  !=========================================================================
  ! Resolve an internal point source spectrum and fill its band bin array.
  ! Priority: par%src_spectrum_file column is (if set) > src_tstar(is) Planck
  ! (only when use_src_tstar) > global par%ion_spectrum file > global par%tstar
  ! Planck.  The single-component fast path passes use_src_tstar = .false. (a
  ! lone source uses par%tstar, not src_tstar - legacy scalars rule); the multi-
  ! component path passes .true.  is_abs marks a physical-type (absolute) file.
  !=========================================================================
  subroutine point_shape(arr, eedge, nnu, nsub, is, use_src_tstar, is_abs)
    implicit none
    real(kind=wp), intent(out) :: arr(:)
    real(kind=wp), intent(in)  :: eedge(:)
    integer,       intent(in)  :: nnu, nsub, is
    logical,       intent(in)  :: use_src_tstar
    logical,       intent(out) :: is_abs
    real(kind=wp), allocatable :: etab(:), ftab(:)
    real(kind=wp) :: tst
    integer :: ntab
    is_abs = .false.
    if (len_trim(par%src_spectrum_file) > 0) then
       if (spec_is_physical()) then
          call read_phys_table(par%src_spectrum_file, is, par%nsource, etab, ftab, ntab)
          call bin_interp(etab, ftab, ntab, eedge, nnu, nsub, arr)
          deallocate(etab, ftab)
          is_abs = .true.
       else
          call bin_lum_from_multicol(par%src_spectrum_file, is, par%nsource, &
                                     eedge, nnu, nsub, arr)
       end if
    else
       if (use_src_tstar) then
          tst = par%src_tstar(is)
       else
          tst = -1.0_wp                 ! global fallback only (fast path)
       end if
       call comp_shape(arr, eedge, nnu, nsub, '', tst, is_abs)
    end if
  end subroutine point_shape

  !=========================================================================
  ! .true. when par%spectrum_type is a physical (absolute) type.
  !=========================================================================
  logical function spec_is_physical()
    select case (trim(par%spectrum_type))
    case ('per_ev', 'per_hz', 'per_ang', 'per_um')
       spec_is_physical = .true.
    case default
       spec_is_physical = .false.
    end select
  end function spec_is_physical

  !=========================================================================
  ! The external field is ON when par%ext_intensity>0, an ISRF preset is named
  ! in par%ext_spectrum, or par%ext_spectrum is a physical-type file (absolute
  ! J), so a preset / absolute field needs no ext_intensity.
  !=========================================================================
  logical function external_on()
    external_on = (par%ext_intensity > 0.0_wp) .or. &
                  (ion_ext_preset_id(par%ext_spectrum) > 0) .or. &
                  (spec_is_physical() .and. len_trim(par%ext_spectrum) > 0)
  end function external_on

  !=========================================================================
  ! ISRF preset id for a par%ext_spectrum value (case-insensitive on the
  ! trimmed string): 1 = draine, 2 = habing, 3 = mathis; 0 = a file path.
  ! (Public so setup can validate the composable-external conditions.)
  !=========================================================================
  integer function ion_ext_preset_id(name)
    implicit none
    character(len=*), intent(in) :: name
    character(len=32) :: s
    integer :: i, c
    s = adjustl(name)
    do i = 1, len_trim(s)
       c = iachar(s(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) s(i:i) = achar(c + 32)
    end do
    select case (trim(s))
    case ('draine');  ion_ext_preset_id = 1
    case ('habing');  ion_ext_preset_id = 2
    case ('mathis');  ion_ext_preset_id = 3
    case default;     ion_ext_preset_id = 0
    end select
  end function ion_ext_preset_id

  !=========================================================================
  ! Normalize a source bin array to its scale.  For a 'shape'/Planck source
  ! (is_abs=.false.) the ionizing segment is rescaled to the scale (legacy).
  ! For a physical-type source (is_abs=.true.) the bins are absolute: an unset
  ! scale (<=0) is DERIVED (bins kept, Lused = ionizing-segment integral), a
  ! set scale rescales the ionizing segment to it.  Lused = the resulting
  ! ionizing-band luminosity.
  !=========================================================================
  subroutine finalize_source(arr, nnu, nfuv, is_abs, scale, Lused)
    use mpi
    implicit none
    real(kind=wp), intent(inout) :: arr(:)
    integer,       intent(in)    :: nnu, nfuv
    logical,       intent(in)    :: is_abs
    real(kind=wp), intent(in)    :: scale
    real(kind=wp), intent(out)   :: Lused
    real(kind=wp) :: lion, ltot
    integer :: ierr
    lion = sum(arr(nfuv+1:nnu))
    if (is_abs) then
       if (scale > 0.0_wp) then          ! rescale the ionizing segment to scale
          if (lion <= 0.0_wp) then
             if (mpar%p_rank == 0) write(*,'(a)') &
                'ERROR: cannot rescale a source with no ionizing-band luminosity '// &
                'to an ionizing luminosity (the spectrum is FUV-only).'
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          end if
          if (mpar%p_rank == 0) write(*,'(a,es12.4,a,es12.4)') &
             ' ION: NOTE: physical source rescaled from file integral', lion, &
             ' to', scale
          arr = arr / lion * scale
          Lused = scale
       else                              ! derive: keep the absolute file bins
          if (lion > 0.0_wp) then
             Lused = lion
          else
             ltot = sum(arr)             ! FUV-only source: keep the FUV bins
             if (ltot <= 0.0_wp) then
                if (mpar%p_rank == 0) write(*,'(a)') &
                   'ERROR: the source spectrum has no luminosity inside the band.'
                call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
             end if
             Lused = 0.0_wp
          end if
       end if
    else                                 ! legacy shape: renormalize to scale
       if (lion <= 0.0_wp) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: the source spectrum has no luminosity inside the ionizing band.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       arr = arr / lion * scale
       Lused = scale
    end if
  end subroutine finalize_source

  !=========================================================================
  ! Build the external field's bin LUMINOSITIES [erg/s] into arr.  Three
  ! spectrum kinds: an ISRF preset (analytic J_E), a physical-type file
  ! (absolute J_E), or a legacy 'shape'/Planck spectrum.  For the absolute
  ! kinds each bin gets L_b = pi*A_surface*J_b (the interior mean intensity of
  ! an isotropic field entering an area A is J_b); if par%ext_intensity is also
  ! set the field is rescaled so the band-integrated interior J equals it.  For
  ! the legacy kind the ionizing shape is renormalized and scaled by
  ! pi*ext_intensity*A (the original behavior).  Returns Lenter (ionizing-
  ! segment power), Jband (band-integrated interior J), uFUV (FUV energy
  ! density (4pi/c) int_FUV J dE), and presid (>0 = preset).
  !=========================================================================
  subroutine build_ext_bins(arr, eedge, nnu, nfuv, nsub, Lenter, Jband, uFUV, presid)
    use mpi
    use octree_mod, only : amr_grid
    implicit none
    real(kind=wp), intent(out) :: arr(:)
    real(kind=wp), intent(in)  :: eedge(:)
    integer,       intent(in)  :: nnu, nfuv, nsub
    real(kind=wp), intent(out) :: Lenter, Jband, uFUV
    integer,       intent(out) :: presid
    real(kind=wp), allocatable :: etab(:), ftab(:)
    real(kind=wp) :: asurf, xr, yr, zr, Jtot, dummyT
    logical :: is_abs, dum
    integer :: ntab, ierr

    !--- illuminated surface area [cm^2].
    if (trim(par%ext_geometry) == 'sph') then
       asurf = fourpi * (par%rmax * par%distance2cm)**2
       if (mpar%p_rank == 0 .and. par%rmax > 0.5_wp* &
           min(amr_grid%xrange, amr_grid%yrange, amr_grid%zrange)) &
          write(*,'(a)') ' ION: WARNING: ext_geometry=sph rmax exceeds the box '// &
          'half-extent; part of the sphere falls outside the box (those '// &
          'entering packets are dropped).'
    else
       xr = amr_grid%xrange * par%distance2cm
       yr = amr_grid%yrange * par%distance2cm
       zr = amr_grid%zrange * par%distance2cm
       asurf = 2.0_wp*(xr*yr + yr*zr + zr*xr)
    end if

    presid = ion_ext_preset_id(par%ext_spectrum)
    if (presid > 0) then
       call bin_preset(presid, par%ext_scale, eedge, nnu, nsub, arr)   ! J bins
       is_abs = .true.
    else if (spec_is_physical() .and. len_trim(par%ext_spectrum) > 0) then
       call read_phys_table(par%ext_spectrum, 1, 1, etab, ftab, ntab)
       call bin_interp(etab, ftab, ntab, eedge, nnu, nsub, arr)        ! J bins
       deallocate(etab, ftab)
       is_abs = .true.
    else
       dummyT = par%ext_tstar
       call comp_shape(arr, eedge, nnu, nsub, par%ext_spectrum, dummyT, dum)
       is_abs = .false.
    end if

    if (is_abs) then
       Jtot = sum(arr)                            ! band-integrated interior J
       if (Jtot <= 0.0_wp) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: the external spectrum has no luminosity inside the band.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       if (par%ext_intensity > 0.0_wp) arr = arr * (par%ext_intensity / Jtot)  ! rescale to J
       arr = pi * asurf * arr                     ! J bins -> bin luminosities
    else
       call finalize_source(arr, nnu, nfuv, .false., 1.0_wp, dummyT)   ! ionizing shape -> 1
       arr = arr * (pi * par%ext_intensity * asurf)
    end if

    !--- diagnostics from the final bin luminosities (interior J_b = L_b/(pi A)).
    Jband  = sum(arr) / (pi * asurf)
    uFUV   = 0.0_wp
    if (nfuv > 0) uFUV = fourpi / clight_cgs * sum(arr(1:nfuv)) / (pi * asurf)
    Lenter = sum(arr(nfuv+1:nnu))
  end subroutine build_ext_bins

  !=========================================================================
  ! Threshold-aligned ionizing bin edges.  Fills eedge(nfuv+1 .. nnu+1)
  ! (the ionizing segment [eion_min, eion_max]) so that every ionization
  ! threshold strictly inside the band falls ON a bin edge.  The ionizing
  ! segment carries nion = nnu - nfuv bins, split among the sub-segments
  ! between consecutive thresholds; each sub-segment is log-uniform.
  !=========================================================================
  subroutine build_aligned_edges(eedge, nfuv, nnu, metal_eth, nmetal)
    use mpi
    implicit none
    real(kind=wp), intent(inout) :: eedge(:)
    integer,       intent(in)    :: nfuv, nnu, nmetal
    real(kind=wp), intent(in)    :: metal_eth(:)
    real(kind=wp) :: cand(256), thr(256), points(256), w(64)
    real(kind=wp) :: tmp, lo, hi, wsum, frac
    integer :: ncand, nthr, nseg, nion, i, j, ipt, ib, r
    integer :: nb(64), off, best, ierr, ntot, itmp, ndrop
    integer :: cpri(256), tpri(256)
    logical :: dup
    real(kind=wp), parameter :: MERGE_TOL = 0.01_wp   ! |ln(T2/T1)| merge

    lo = par%eion_min;  hi = par%eion_max
    nion = nnu - nfuv

    !--- 1) collect candidate thresholds STRICTLY inside (eion_min, eion_max):
    !---    H I / He I / He II plus the active metal photoionization edges
    !---    (gathered by the caller from the species registry).  H/He carry
    !---    priority 1 so that a near-coincident metal edge does not displace
    !---    them in the merge (e.g. C II 24.383 vs He I 24.587 within 1%).
    ncand = 0
    call add_cand(cand, cpri, ncand, eth_HI,   1, lo, hi)
    call add_cand(cand, cpri, ncand, eth_HeI,  1, lo, hi)
    call add_cand(cand, cpri, ncand, eth_HeII, 1, lo, hi)
    do i = 1, nmetal
       if (metal_eth(i) > 0.0_wp) &
          call add_cand(cand, cpri, ncand, metal_eth(i), 0, lo, hi)
    end do

    !--- 2) sort ascending (insertion sort; ncand is small; carry priority)
    do i = 2, ncand
       tmp = cand(i);  itmp = cpri(i);  j = i - 1
       do while (j >= 1)
          if (cand(j) <= tmp) exit
          cand(j+1) = cand(j);  cpri(j+1) = cpri(j);  j = j - 1
       end do
       cand(j+1) = tmp;  cpri(j+1) = itmp
    end do

    !--- 2b) merge near-coincident thresholds.  Within a merged cluster keep
    !---    the higher-priority member (H/He over metal); ties keep the first.
    nthr = 0
    do i = 1, ncand
       dup = .false.
       if (nthr >= 1) then
          if (abs(log(cand(i)/thr(nthr))) < MERGE_TOL) dup = .true.
       end if
       if (dup) then
          if (cpri(i) > tpri(nthr)) then     ! promote H/He onto the edge
             thr(nthr)  = cand(i)
             tpri(nthr) = cpri(i)
          end if
       else
          nthr = nthr + 1
          thr(nthr)  = cand(i)
          tpri(nthr) = cpri(i)
       end if
    end do

    !--- 2c) graceful degradation: if the thresholds outnumber the ionizing
    !---    bins (nseg = nthr+1 > nion), drop the lowest-priority (metal)
    !---    thresholds - keeping H/He always - until they fit, removing the
    !---    most redundant metal first (smallest log-gap to a neighbor; the
    !---    band edges lo/hi are the outer neighbors).  So the aligned band is
    !---    safe at any nnu_ion: it resolves the big He edges plus as many
    !---    metal edges as the bins allow.
    ndrop = 0
    do while (nthr + 1 > nion)
       best = 0;  frac = 1.0e30_wp
       do i = 1, nthr
          if (tpri(i) /= 0) cycle            ! never drop H/He (priority 1)
          if (i == 1) then;  tmp = log(thr(i)/lo)
          else;              tmp = log(thr(i)/thr(i-1));  end if
          if (i == nthr) then;  tmp = min(tmp, log(hi/thr(i)))
          else;                 tmp = min(tmp, log(thr(i+1)/thr(i)));  end if
          if (tmp < frac) then;  frac = tmp;  best = i;  end if
       end do
       if (best == 0) exit                   ! only H/He remain - cannot reduce
       do i = best, nthr-1
          thr(i) = thr(i+1);  tpri(i) = tpri(i+1)
       end do
       nthr = nthr - 1;  ndrop = ndrop + 1
    end do
    if (ndrop > 0 .and. mpar%p_rank == 0) write(*,'(a,i0,a,i0,a)') &
       ' ION: align_edges: nnu_ion=', nion, ' too small for all edges; dropped ', &
       ndrop, ' metal threshold(s) (H/He kept; raise nnu_ion to align all).'

    !--- 3) split points: eion_min, interior thresholds, eion_max
    points(1) = lo
    do i = 1, nthr
       points(i+1) = thr(i)
    end do
    points(nthr+2) = hi
    nseg = nthr + 1

    !--- 5) need at least one bin per sub-segment.  After degradation this
    !---    only trips for an absurd nnu_ion below the H/He edge count.
    if (nion < nseg) then
       if (mpar%p_rank == 0) write(*,'(a,i0,a,i0,a,i0,a)') &
          'ERROR: ion_align_edges: nnu_ion=', nion, ' cannot fit ', nthr, &
          ' H/He thresholds; raise nnu_ion to >= ', nseg, '.'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if

    !--- 4) distribute nion bins over the nseg sub-segments by log width,
    !---    largest-remainder (Hamilton) with a floor of 1 per sub-segment.
    wsum = 0.0_wp
    do j = 1, nseg
       w(j) = log(points(j+1)) - log(points(j))
       wsum = wsum + w(j)
    end do
    ntot = 0
    do j = 1, nseg
       nb(j) = max(1, int(real(nion,wp)*w(j)/wsum))   ! floor(ideal), >= 1
       ntot  = ntot + nb(j)
    end do
    r = nion - ntot
    do while (r > 0)                 ! hand out by largest fractional remainder
       best = 1;  frac = -1.0_wp
       do j = 1, nseg
          tmp = real(nion,wp)*w(j)/wsum - real(nb(j),wp)
          if (tmp > frac) then
             frac = tmp;  best = j
          end if
       end do
       nb(best) = nb(best) + 1;  r = r - 1
    end do
    do while (r < 0)                 ! too many (floor-of-1 excess): remove
       best = 0;  frac = -1.0e30_wp
       do j = 1, nseg
          if (nb(j) <= 1) cycle
          tmp = real(nb(j),wp) - real(nion,wp)*w(j)/wsum   ! excess over ideal
          if (tmp > frac) then
             frac = tmp;  best = j
          end if
       end do
       if (best == 0) exit           ! cannot reduce further (all at 1)
       nb(best) = nb(best) - 1;  r = r + 1
    end do

    !--- 6) fill eedge inside each sub-segment (log-uniform); sub-segment
    !---    boundaries are the split points = thresholds.
    off = nfuv + 1                   ! eedge index of points(1) = eion_min
    eedge(off) = points(1)
    ipt = off
    do j = 1, nseg
       do ib = 1, nb(j)
          ipt = ipt + 1
          eedge(ipt) = exp(log(points(j)) + (log(points(j+1)) - log(points(j))) &
                       *real(ib,wp)/real(nb(j),wp))
       end do
    end do

    !--- 7) rank-0 log
    if (mpar%p_rank == 0) then
       write(*,'(a,i0,a,i0,a)') ' ION: align_edges on: ', nthr, &
          ' interior thresholds, ', nseg, ' sub-segments (bins each:'
       write(*,'(a,20(1x,i0))') '      ', (nb(j), j=1,nseg)
    end if
  end subroutine build_aligned_edges

  !=========================================================================
  ! Append E (with priority pri) to the candidate list if strictly inside
  ! (lo, hi).  Priority resolves near-coincident merges (H/He = 1 > metal 0).
  !=========================================================================
  subroutine add_cand(cand, cpri, ncand, E, pri, lo, hi)
    implicit none
    real(kind=wp), intent(inout) :: cand(:)
    integer,       intent(inout) :: cpri(:)
    integer,       intent(inout) :: ncand
    real(kind=wp), intent(in)    :: E, lo, hi
    integer,       intent(in)    :: pri
    if (E > lo .and. E < hi) then
       ncand = ncand + 1
       cand(ncand) = E
       cpri(ncand) = pri
    end if
  end subroutine add_cand

  !=========================================================================
  ! Planck B_nu at photon energy E [eV] and temperature T [K], arbitrary
  ! normalization per unit E (the constant prefactor cancels in ion_lum).
  !=========================================================================
  real(kind=wp) function planck_nu(E, T) result(B)
    real(kind=wp), intent(in) :: E, T
    real(kind=wp) :: x
    x = E * ev2erg / (kboltz_cgs * T)
    if (x < 700.0_wp) then
       B = E**3 / (exp(x) - 1.0_wp)
    else
       B = E**3 * exp(-x)          ! Wien tail, avoids overflow
    end if
  end function planck_nu

  !=========================================================================
  ! 2-column spectrum file (E [eV], L_E [arb per eV]; '#' comments), linear
  ! interpolation onto the sub-grid; zero outside the tabulated range.
  !=========================================================================
  subroutine bin_lum_from_file(fname, eedge, nnu, nsub, arr)
    use mpi
    implicit none
    character(len=*), intent(in) :: fname
    real(kind=wp), intent(in) :: eedge(:)
    integer,       intent(in) :: nnu, nsub
    real(kind=wp), intent(out):: arr(:)
    real(kind=wp), allocatable :: etab(:), ftab(:)
    character(len=256) :: line
    real(kind=wp) :: e1, e2, es, fsum, ei, fi
    integer :: unit, ios, ntab, i, k, j, ierr

    ntab = 0
    open(newunit=unit, file=trim(fname), status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(2a)') 'ERROR: cannot open ', trim(fname)
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       ntab = ntab + 1
    end do
    allocate(etab(ntab), ftab(ntab))
    rewind(unit)
    i = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       i = i + 1
       read(line,*) etab(i), ftab(i)
    end do
    close(unit)

    do i = 1, nnu
       e1 = eedge(i);  e2 = eedge(i+1)
       fsum = 0.0_wp
       do k = 1, nsub
          es = e1 + (e2 - e1)*(real(k,wp) - 0.5_wp)/real(nsub,wp)
          fi = 0.0_wp
          if (es >= etab(1) .and. es <= etab(ntab)) then
             do j = 1, ntab-1
                if (es <= etab(j+1)) then
                   ei = (es - etab(j)) / (etab(j+1) - etab(j))
                   fi = ftab(j)*(1.0_wp - ei) + ftab(j+1)*ei
                   exit
                end if
             end do
          end if
          fsum = fsum + fi
       end do
       arr(i) = fsum * (e2 - e1)/real(nsub,wp)
    end do
    deallocate(etab, ftab)
  end subroutine bin_lum_from_file

  !=========================================================================
  ! Multi-column spectrum file for the internal point sources: column 1 =
  ! E [eV] (ascending), columns 2.. = L_E of source 1.. ('#' comments).  Bins
  ! the icol-th source (file column 1+icol) exactly like bin_lum_from_file
  ! (linear interpolation onto the NSUB sub-grid, zero outside the tabulated
  ! range).  Aborts when the file has fewer than nsrc source columns; a single
  ! rank-0 note (issued for icol=1) records any extra columns beyond nsrc.
  !=========================================================================
  subroutine bin_lum_from_multicol(fname, icol, nsrc, eedge, nnu, nsub, arr)
    use mpi
    implicit none
    character(len=*), intent(in) :: fname
    integer,       intent(in) :: icol, nsrc, nnu, nsub
    real(kind=wp), intent(in) :: eedge(:)
    real(kind=wp), intent(out):: arr(:)
    real(kind=wp), allocatable :: etab(:), ftab(:), row(:)
    character(len=1024) :: line
    real(kind=wp) :: e1, e2, es, fsum, ei, fi, extra
    integer :: unit, ios, ios2, ntab, i, k, j, ierr
    logical :: has_extra

    ntab = 0
    open(newunit=unit, file=trim(fname), status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(2a)') 'ERROR: cannot open ', trim(fname)
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       ntab = ntab + 1
    end do
    allocate(etab(ntab), ftab(ntab), row(nsrc+1))
    rewind(unit)
    has_extra = .false.
    i = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       i = i + 1
       read(line,*,iostat=ios2) row(1:nsrc+1)          ! E + one L_E per source
       if (ios2 /= 0) then
          if (mpar%p_rank == 0) write(*,'(3a,i0,a)') &
             'ERROR: ', trim(fname), ' needs 1 + nsource(=', nsrc, &
             ') columns (E [eV] and one L_E for each source); it has too few.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       if (i == 1) then                                ! detect extra columns once
          read(line,*,iostat=ios2) row(1:nsrc+1), extra
          if (ios2 == 0) has_extra = .true.
       end if
       etab(i) = row(1)
       ftab(i) = row(1+icol)
    end do
    close(unit)
    if (has_extra .and. icol == 1 .and. mpar%p_rank == 0) &
       write(*,'(3a,i0,a)') ' ION: NOTE: ', trim(fname), &
       ' has more than nsource(=', nsrc, ') data columns; using the first nsource.'

    do i = 1, nnu
       e1 = eedge(i);  e2 = eedge(i+1)
       fsum = 0.0_wp
       do k = 1, nsub
          es = e1 + (e2 - e1)*(real(k,wp) - 0.5_wp)/real(nsub,wp)
          fi = 0.0_wp
          if (es >= etab(1) .and. es <= etab(ntab)) then
             do j = 1, ntab-1
                if (es <= etab(j+1)) then
                   ei = (es - etab(j)) / (etab(j+1) - etab(j))
                   fi = ftab(j)*(1.0_wp - ei) + ftab(j+1)*ei
                   exit
                end if
             end do
          end if
          fsum = fsum + fi
       end do
       arr(i) = fsum * (e2 - e1)/real(nsub,wp)
    end do
    deallocate(etab, ftab, row)
  end subroutine bin_lum_from_multicol

  !=========================================================================
  ! Integrate a tabulated (E [eV] ascending, X_E per eV) spectrum onto the band
  ! bins on the NSUB midpoint sub-grid (linear interpolation, zero outside the
  ! table).  arr(i) = int_bin X_E dE.  Same quadrature as bin_lum_from_file, so
  ! a physical file bins exactly like a shape file (only the table units differ).
  !=========================================================================
  subroutine bin_interp(etab, ftab, ntab, eedge, nnu, nsub, arr)
    implicit none
    real(kind=wp), intent(in)  :: etab(:), ftab(:), eedge(:)
    integer,       intent(in)  :: ntab, nnu, nsub
    real(kind=wp), intent(out) :: arr(:)
    real(kind=wp) :: e1, e2, es, fsum, ei, fi
    integer :: i, k, j
    do i = 1, nnu
       e1 = eedge(i);  e2 = eedge(i+1)
       fsum = 0.0_wp
       do k = 1, nsub
          es = e1 + (e2 - e1)*(real(k,wp) - 0.5_wp)/real(nsub,wp)
          fi = 0.0_wp
          if (es >= etab(1) .and. es <= etab(ntab)) then
             do j = 1, ntab-1
                if (es <= etab(j+1)) then
                   ei = (es - etab(j)) / (etab(j+1) - etab(j))
                   fi = ftab(j)*(1.0_wp - ei) + ftab(j+1)*ei
                   exit
                end if
             end do
          end if
          fsum = fsum + fi
       end do
       arr(i) = fsum * (e2 - e1)/real(nsub,wp)
    end do
  end subroutine bin_interp

  !=========================================================================
  ! Read a physical-type spectrum file and return its icol-th value column as
  ! (E [eV] ascending, X_E per eV), converted from par%spectrum_type:
  !   'per_ev'  col1 E [eV],      X_E = value
  !   'per_hz'  col1 nu [Hz],     E = h*nu/e, X_E = X_nu*e/h
  !   'per_ang' col1 lambda [A],  E = hc/lambda, X_E = X_lam*hc/E^2   (hc [eV A])
  !   'per_um'  col1 lambda [um], E = hc/lambda, X_E = X_lam*hc/E^2   (hc [eV um])
  ! ncol = number of value columns (1 for a 2-column file; nsource for the
  ! multi-column source file).  Aborts when a row has fewer than ncol+1 columns;
  ! a rank-0 note (once, icol=1) records extra columns.  Wavelength tables are
  ! E-descending after conversion, so the result is sorted ascending in E.
  !=========================================================================
  subroutine read_phys_table(fname, icol, ncol, etab, ftab, ntab)
    use mpi
    implicit none
    character(len=*), intent(in) :: fname
    integer,          intent(in) :: icol, ncol
    real(kind=wp), allocatable, intent(out) :: etab(:), ftab(:)
    integer,          intent(out) :: ntab
    real(kind=wp), allocatable :: row(:)
    character(len=2048) :: line
    real(kind=wp) :: c1, cv, ee, xe, extra, tmp, hc
    integer :: unit, ios, ios2, i, j, ierr
    logical :: has_extra

    ntab = 0
    open(newunit=unit, file=trim(fname), status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(2a)') 'ERROR: cannot open ', trim(fname)
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       ntab = ntab + 1
    end do
    allocate(etab(ntab), ftab(ntab), row(ncol+1))
    rewind(unit)
    has_extra = .false.
    i = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       i = i + 1
       read(line,*,iostat=ios2) row(1:ncol+1)
       if (ios2 /= 0) then
          if (mpar%p_rank == 0) write(*,'(3a,i0,a)') &
             'ERROR: ', trim(fname), ' needs 1 + ', ncol, &
             ' columns for spectrum_type /= shape; it has too few.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       if (i == 1) then
          read(line,*,iostat=ios2) row(1:ncol+1), extra
          if (ios2 == 0) has_extra = .true.
       end if
       c1 = row(1);  cv = row(1+icol)
       select case (trim(par%spectrum_type))
       case ('per_ev')
          ee = c1
          xe = cv
       case ('per_hz')
          ee = h_planck_cgs * c1 / ev2erg
          xe = cv * ev2erg / h_planck_cgs
       case ('per_ang')
          hc = hc_evAng                     ! eV * Angstrom
          ee = hc / c1
          xe = cv * hc / (ee*ee)
       case ('per_um')
          hc = hc_evAng * 1.0e-4_wp         ! eV * micron
          ee = hc / c1
          xe = cv * hc / (ee*ee)
       case default
          ee = c1;  xe = cv                 ! unreachable (physical types only)
       end select
       etab(i) = ee;  ftab(i) = xe
    end do
    close(unit)
    if (has_extra .and. icol == 1 .and. mpar%p_rank == 0) &
       write(*,'(3a,i0,a)') ' ION: NOTE: ', trim(fname), &
       ' has more than ', ncol, ' value column(s); using the first as needed.'

    !--- sort ascending in E (wavelength tables arrive descending).
    do i = 2, ntab
       do j = i, 2, -1
          if (etab(j-1) <= etab(j)) exit
          tmp = etab(j-1);  etab(j-1) = etab(j);  etab(j) = tmp
          tmp = ftab(j-1);  ftab(j-1) = ftab(j);  ftab(j) = tmp
       end do
    end do
    deallocate(row)
  end subroutine read_phys_table

  !=========================================================================
  ! Integrate an analytic ISRF preset J_E [erg/s/cm^2/sr/eV] onto the band bins
  ! on the NSUB midpoint sub-grid, times the dimensionless par%ext_scale.
  ! arr(i) = ext_scale * int_bin J_E dE  [erg/s/cm^2/sr].
  !=========================================================================
  subroutine bin_preset(id, scale, eedge, nnu, nsub, arr)
    implicit none
    integer,       intent(in)  :: id, nnu, nsub
    real(kind=wp), intent(in)  :: scale, eedge(:)
    real(kind=wp), intent(out) :: arr(:)
    real(kind=wp) :: e1, e2, es, fsum
    integer :: i, k
    do i = 1, nnu
       e1 = eedge(i);  e2 = eedge(i+1)
       fsum = 0.0_wp
       do k = 1, nsub
          es = e1 + (e2 - e1)*(real(k,wp) - 0.5_wp)/real(nsub,wp)
          fsum = fsum + preset_je(es, id)
       end do
       arr(i) = fsum * (e2 - e1)/real(nsub,wp) * scale
    end do
  end subroutine bin_preset

  !=========================================================================
  ! ISRF preset mean-intensity density J_E [erg/s/cm^2/sr/eV] at photon energy
  ! E [eV] (before par%ext_scale).  id: 1 = Draine (1978), 2 = Habing (Draine
  ! shape / 1.71), 3 = Mathis, Mezger & Panagia (1983).  Zero outside the FUV.
  !=========================================================================
  real(kind=wp) function preset_je(E, id) result(Je)
    implicit none
    real(kind=wp), intent(in) :: E
    integer,       intent(in) :: id
    select case (id)
    case (1);  Je = draine_je(E)
    case (2);  Je = draine_je(E) / 1.71_wp
    case (3);  Je = mathis_je(E)
    case default;  Je = 0.0_wp
    end select
  end function preset_je

  !=========================================================================
  ! Draine (1978) FUV interstellar field, J_E [erg/s/cm^2/sr/eV] on 5-13.6 eV
  ! (zero outside).  The polynomial is the energy-weighted photon fit
  ! F_ph(E) = E*N(E) = 1.658e6 E^2 - 2.152e5 E^3 + 6.919e3 E^4, so the energy
  ! mean intensity is J_E = ev2erg*F_ph.  This reproduces the canonical FUV
  ! energy density u = (4pi/c) int_{6}^{13.6} J_E dE = 8.94e-14 erg/cm^3
  ! (G0 = 1.69 in Habing units).
  !=========================================================================
  real(kind=wp) function draine_je(E) result(Je)
    implicit none
    real(kind=wp), intent(in) :: E
    real(kind=wp) :: Fph
    if (E < 5.0_wp .or. E > 13.6_wp) then
       Je = 0.0_wp
    else
       Fph = 1.658e6_wp*E**2 - 2.152e5_wp*E**3 + 6.919e3_wp*E**4
       Je  = ev2erg * Fph
    end if
  end function draine_je

  !=========================================================================
  ! Mathis, Mezger & Panagia (1983) ISRF, J_E [erg/s/cm^2/sr/eV].  The
  ! mean-intensity shape J_lambda is evaluated in SI [W m^-2 m^-1 sr^-1] with
  ! the coefficients as in SEDust radfield.f90 (mathis_jl_si), then converted:
  ! J_E = J_lambda(SI) * 10 * (hc [cm eV]) / E^2, where the 10 folds
  ! W m^-2 m^-1 -> erg s^-1 cm^-2 cm^-1 and hc/E^2 = |dlambda/dE|.
  !=========================================================================
  real(kind=wp) function mathis_je(E) result(Je)
    implicit none
    real(kind=wp), intent(in) :: E
    real(kind=wp) :: lam_um
    lam_um = hc_evAng * 1.0e-4_wp / E                  ! 1.239842 / E  [um]
    Je = mathis_jl_si(lam_um) * (hc_evAng * 1.0e-7_wp) / (E*E)
  end function mathis_je

  !=========================================================================
  ! Mathis, Mezger & Panagia (1983) mean intensity J_lambda in SI
  ! [W m^-2 m^-1 sr^-1] as a function of lambda [um]; coefficients as in
  ! SEDust radfield.f90 (UV piecewise power laws + three diluted blackbodies,
  ! Draine-corrected 4000 K dilution).  CMB is omitted (negligible in the band).
  !=========================================================================
  real(kind=wp) function mathis_jl_si(lam_um) result(Jl)
    implicit none
    real(kind=wp), intent(in) :: lam_um
    if (lam_um < 0.0912_wp) then
       Jl = 0.0_wp
    else if (lam_um < 0.110_wp) then
       Jl = 3069.0_wp * lam_um**3.4172_wp
    else if (lam_um < 0.134_wp) then
       Jl = 1.627_wp
    else if (lam_um < 0.250_wp) then
       Jl = 0.0566_wp * lam_um**(-1.6678_wp)
    else
       Jl = 1.0e-14_wp  * bbody_si(7500.0_wp, lam_um) &
          + 1.65e-13_wp * bbody_si(4000.0_wp, lam_um) &
          + 4.0e-13_wp  * bbody_si(3000.0_wp, lam_um)
    end if
  end function mathis_jl_si

  !=========================================================================
  ! Planck spectral radiance B_lambda(T) in SI [W m^-2 m^-1 sr^-1], lambda [um];
  ! constants as in SEDust radfield.f90 (bbody).
  !=========================================================================
  real(kind=wp) function bbody_si(T, lam_um) result(B)
    implicit none
    real(kind=wp), intent(in) :: T, lam_um
    real(kind=wp), parameter :: c_si  = 2.99792458e8_wp
    real(kind=wp), parameter :: h_si  = 6.62606957e-34_wp
    real(kind=wp), parameter :: kB_si = 1.3806488e-23_wp
    real(kind=wp), parameter :: hc2   = 2.0_wp*h_si*c_si**2
    real(kind=wp), parameter :: hckB  = h_si*c_si/kB_si
    real(kind=wp) :: lam_m, x
    if (T <= 0.0_wp .or. lam_um <= 0.0_wp) then
       B = 0.0_wp;  return
    end if
    lam_m = lam_um * 1.0e-6_wp
    x = hckB / (T*lam_m)
    if (x >= 700.0_wp) then
       B = 0.0_wp
    else
       B = hc2 / lam_m**5 / (exp(x) - 1.0_wp)
    end if
  end function bbody_si

  !=========================================================================
  ! rank-0 log for the single-component fast path (point source or external-
  ! only field), including the derived-luminosity and preset diagnostics.
  !=========================================================================
  subroutine log_fast_path(nnu, nfuv, Lion, Jband, uFUV, presid, derived)
    implicit none
    integer,       intent(in) :: nnu, nfuv, presid
    real(kind=wp), intent(in) :: Lion, Jband, uFUV
    logical,       intent(in) :: derived
    call band_log(nnu, nfuv)
    if (par%nsource == 0) then
       !--- external-only field.
       if (presid > 0) then
          write(*,'(3a,es12.4)') ' ION: external ISRF preset = ', &
             trim(par%ext_spectrum), ',  ext_scale =', par%ext_scale
       else if (len_trim(par%ext_spectrum) > 0) then
          write(*,'(4a)') ' ION: external spectrum file = ', trim(par%ext_spectrum), &
             ',  type = ', trim(par%spectrum_type)
       else if (par%ext_tstar > 0.0_wp) then
          write(*,'(a,es12.4)') ' ION: external Planck, ext_tstar [K] = ', par%ext_tstar
       else if (len_trim(par%ion_spectrum) > 0) then
          write(*,'(4a)') ' ION: spectrum file = ', trim(par%ion_spectrum), &
             ',  type = ', trim(par%spectrum_type)
       else
          write(*,'(a,es12.4)') ' ION: Planck spectrum, tstar [K] = ', par%tstar
       end if
       write(*,'(3a)') ' ION: external field, geometry = ', &
          trim(par%ext_geometry), ' (isotropic)'
       if (par%ext_intensity > 0.0_wp) &
          write(*,'(a,es12.4)') ' ION: ext_intensity J [erg/s/cm^2/sr] = ', par%ext_intensity
       write(*,'(a,es12.4)') ' ION: entering ionizing power [erg/s]   = ', Lion
       write(*,'(a,es12.4)') ' ION: band-integrated J [erg/s/cm^2/sr] = ', Jband
       if (presid > 0 .or. uFUV > 0.0_wp) &
          write(*,'(a,es12.4,a,f9.3)') ' ION: u_FUV [erg/cm^3] =', uFUV, &
             ',  G0 (Habing) =', uFUV/HABING_U_FUV
    else
       !--- single point source.
       if (len_trim(par%src_spectrum_file) > 0) then
          write(*,'(4a)') ' ION: spectrum = column 1 of ', trim(par%src_spectrum_file), &
             ',  type = ', trim(par%spectrum_type)
       else if (len_trim(par%ion_spectrum) > 0) then
          write(*,'(4a)') ' ION: spectrum file = ', trim(par%ion_spectrum), &
             ',  type = ', trim(par%spectrum_type)
       else
          write(*,'(a,es12.4)') ' ION: Planck spectrum, tstar [K] = ', par%tstar
       end if
       if (derived) then
          write(*,'(a,es12.4)') ' ION: derived ionizing luminosity [erg/s] = ', Lion
       else
          write(*,'(a,es12.4)') ' ION: ionizing luminosity [erg/s] = ', Lion
       end if
    end if
    if (nfuv > 0 .and. Lion > 0.0_wp) write(*,'(a,es12.4,a,es10.3)') &
       ' ION: band total with FUV        = ', ion_Ltot, &
       ',  L_FUV/L_ion = ', (ion_Ltot - Lion)/Lion
  end subroutine log_fast_path

  !=========================================================================
  ! Ionizing source packet.  With a single source component the legacy fast
  ! path is taken (position/direction from par%source_geometry, bin from
  ! ion_cdf).  With several components a component is drawn in proportion to
  ! its band luminosity (src_cdf), the packet is emitted from it, and its
  ! frequency bin is drawn from that component's spectrum (ion_cdf_src):
  !   comp_kind = i>=1 - internal point source i (par%src_x/y/z(i)), isotropic;
  !   comp_kind = 0    - external field (par%ext_geometry 'rec'/'sph').
  ! Every packet carries the same Lpacket = ion_Ltot / nphotons, so component c
  ! collects a band luminosity equal to its share of ion_Ltot.
  !=========================================================================
  subroutine gen_ion_photon(photon)
    use random,     only : rand_number
    use octree_mod, only : amr_find_leaf
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp) :: u
    integer :: i, ic, ik

    if (.not. multi_src) then
       !--- single-component fast path (legacy RNG stream): emit by geometry.
       select case (trim(par%source_geometry))
       case ('point')
          call emit_point(photon, par%xs_point, par%ys_point, par%zs_point)
       case ('external_rec')
          call emit_external_rec(photon)
       case ('external_sph')
          call emit_external_sph(photon)
       end select
       u = rand_number()
       photon%inu = nnu_band
       do i = 1, nnu_band
          if (u <= ion_cdf(i)) then
             photon%inu = i
             exit
          end if
       end do
    else
       !--- multi-component: pick a component, emit from it, draw its bin.
       u = rand_number()
       ic = ncomp
       do i = 1, ncomp
          if (u <= src_cdf(i)) then
             ic = i
             exit
          end if
       end do
       ik = comp_kind(ic)
       if (ik >= 1) then
          call emit_point(photon, par%src_x(ik), par%src_y(ik), par%src_z(ik))
       else if (trim(par%ext_geometry) == 'sph') then
          call emit_external_sph(photon)
       else
          call emit_external_rec(photon)
       end if
       u = rand_number()
       photon%inu = nnu_band
       do i = 1, nnu_band
          if (u <= ion_cdf_src(i, ic)) then
             photon%inu = i
             exit
          end if
       end do
    end if

    photon%wgt     = 1.0_wp
    photon%Lpacket = ion_Ltot / real(par%nphotons, wp)
    photon%nscatt  = 0
    photon%inside  = .true.
    photon%icell_amr = amr_find_leaf(photon%x, photon%y, photon%z)
  end subroutine gen_ion_photon

  !=========================================================================
  ! Quasi-random ionizing launch.  The launch uniforms carry FIXED dimension
  ! semantics so every packet consumes the same coordinate meaning
  ! (docs/QUASI_RANDOM_LAUNCH.md).  Two layouts:
  !
  !  * SINGLE internal point source (stage-1, 3 dims) - preserved bit-for-bit:
  !      u(1) -> frequency bin (ion_cdf inverse);  u(2) -> mu = 2u-1;
  !      u(3) -> phi = 2 pi u.
  !  * SUPERSET (7 dims) - every other configuration (single external,
  !      multiple points, mixed):
  !      u(1) source component (src_cdf; unused with one component);
  !      u(2) frequency bin (the component's CDF);
  !      u(3) polar / incidence angle;   u(4) azimuth;
  !      u(5) rectangular entry face / sphere entry-point cos;
  !      u(6) first surface coordinate / sphere entry-point azimuth;
  !      u(7) second surface coordinate (rectangular only).
  !
  ! The emission expressions (isotropic point, Lambert external rec/sph, entry-
  ! surface pick, nudges, from_external/snx..snz flags) mirror the legacy
  ! gen_ion_photon path EXACTLY; only the uniforms are supplied instead of drawn.
  ! No Mersenne Twister draw is consumed, so the launch set is indexed by the
  ! global photon number and is independent of the MPI task count.
  !=========================================================================
  subroutine gen_ion_photon_qmc(photon, u)
    use octree_mod, only : amr_find_leaf
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp),     intent(in)    :: u(:)
    integer :: i, ic, ik

    if ((.not. multi_src) .and. trim(par%source_geometry) == 'point') then
       !--- single point source: stage-1 layout (u1 = frequency).
       call emit_point_qmc(photon, par%xs_point, par%ys_point, par%zs_point, &
                           u(2), u(3))
       photon%inu = bin_from_cdf(ion_cdf, u(1))
    else if (.not. multi_src) then
       !--- single external component (nsource=0): superset angle/surface dims,
       !--- frequency from the global CDF (u2).  d1 (component) is unused.
       if (trim(par%source_geometry) == 'external_sph') then
          call emit_external_sph_qmc(photon, u(3), u(4), u(5), u(6))
       else
          call emit_external_rec_qmc(photon, u(3), u(4), u(5), u(6), u(7))
       end if
       photon%inu = bin_from_cdf(ion_cdf, u(2))
    else
       !--- multi-component: pick a component (u1), emit from it, draw its bin.
       ic = ncomp
       do i = 1, ncomp
          if (u(1) <= src_cdf(i)) then
             ic = i;  exit
          end if
       end do
       ik = comp_kind(ic)
       if (ik >= 1) then
          call emit_point_qmc(photon, par%src_x(ik), par%src_y(ik), &
                              par%src_z(ik), u(3), u(4))
       else if (trim(par%ext_geometry) == 'sph') then
          call emit_external_sph_qmc(photon, u(3), u(4), u(5), u(6))
       else
          call emit_external_rec_qmc(photon, u(3), u(4), u(5), u(6), u(7))
       end if
       photon%inu = bin_from_cdf(ion_cdf_src(:,ic), u(2))
    end if

    photon%wgt     = 1.0_wp
    photon%Lpacket = ion_Ltot / real(par%nphotons, wp)
    photon%nscatt  = 0
    photon%inside  = .true.
    photon%icell_amr = amr_find_leaf(photon%x, photon%y, photon%z)
  end subroutine gen_ion_photon_qmc

  !=========================================================================
  ! Frequency bin from a monotone CDF and a launch uniform uf: the smallest bin
  ! i with uf <= cdf(i).  Same linear inverse as gen_ion_photon.
  !=========================================================================
  integer function bin_from_cdf(cdf, uf) result(inu)
    implicit none
    real(kind=wp), intent(in) :: cdf(:)
    real(kind=wp), intent(in) :: uf
    integer :: i
    inu = nnu_band
    do i = 1, nnu_band
       if (uf <= cdf(i)) then
          inu = i;  exit
       end if
    end do
  end function bin_from_cdf

  !=========================================================================
  ! Emit an internal point source at (xs,ys,zs) with an isotropic 4pi
  ! direction.
  !=========================================================================
  subroutine emit_point(photon, xs, ys, zs)
    use random, only : rand_number
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp),     intent(in)    :: xs, ys, zs
    real(kind=wp) :: cost, sint, phi
    photon%from_external = .false.
    photon%x = xs
    photon%y = ys
    photon%z = zs
    cost = 2.0_wp*rand_number() - 1.0_wp
    sint = sqrt(1.0_wp - cost*cost)
    phi  = twopi*rand_number()
    photon%kx = sint*cos(phi)
    photon%ky = sint*sin(phi)
    photon%kz = cost
  end subroutine emit_point

  !=========================================================================
  ! Emit an isotropic external field entering the box faces.  Pick a face
  ! weighted by area, a uniform point on it, then a cosine-weighted (Lambert)
  ! direction into the inward hemisphere -> isotropic interior field.
  !=========================================================================
  subroutine emit_external_rec(photon)
    use random,     only : rand_number
    use octree_mod, only : amr_grid
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp) :: mu, sinm, phi, cphi, sphi, nudge, span
    real(kind=wp) :: farea(6), fcum(6), fpick
    integer :: i, iface
    span  = min(amr_grid%xrange, amr_grid%yrange, amr_grid%zrange)
    nudge = 1.0e-6_wp * span
    farea(1) = amr_grid%yrange*amr_grid%zrange     ! +x face
    farea(2) = farea(1)                            ! -x face
    farea(3) = amr_grid%zrange*amr_grid%xrange     ! +y face
    farea(4) = farea(3)                            ! -y face
    farea(5) = amr_grid%xrange*amr_grid%yrange     ! +z face
    farea(6) = farea(5)                            ! -z face
    fcum(1) = farea(1)
    do i = 2, 6
       fcum(i) = fcum(i-1) + farea(i)
    end do
    fpick = rand_number()*fcum(6)
    iface = 6
    do i = 1, 6
       if (fpick <= fcum(i)) then
          iface = i;  exit
       end if
    end do
    mu   = sqrt(rand_number())              ! cos(angle from inward normal)
    sinm = sqrt(1.0_wp - mu*mu)
    phi  = twopi*rand_number()
    cphi = cos(phi);  sphi = sin(phi)
    photon%from_external = .true.
    select case (iface)
    case (1)   ! +x face, inward normal (-1,0,0)
       photon%x = amr_grid%xmax - nudge
       photon%y = amr_grid%ymin + amr_grid%yrange*rand_number()
       photon%z = amr_grid%zmin + amr_grid%zrange*rand_number()
       photon%kx = -mu;  photon%ky = sinm*cphi;  photon%kz = sinm*sphi
       photon%snx = -1.0_wp;  photon%sny = 0.0_wp;  photon%snz = 0.0_wp
    case (2)   ! -x face, inward normal (+1,0,0)
       photon%x = amr_grid%xmin + nudge
       photon%y = amr_grid%ymin + amr_grid%yrange*rand_number()
       photon%z = amr_grid%zmin + amr_grid%zrange*rand_number()
       photon%kx =  mu;  photon%ky = sinm*cphi;  photon%kz = sinm*sphi
       photon%snx =  1.0_wp;  photon%sny = 0.0_wp;  photon%snz = 0.0_wp
    case (3)   ! +y face, inward normal (0,-1,0)
       photon%x = amr_grid%xmin + amr_grid%xrange*rand_number()
       photon%y = amr_grid%ymax - nudge
       photon%z = amr_grid%zmin + amr_grid%zrange*rand_number()
       photon%kx = sinm*cphi;  photon%ky = -mu;  photon%kz = sinm*sphi
       photon%snx = 0.0_wp;  photon%sny = -1.0_wp;  photon%snz = 0.0_wp
    case (4)   ! -y face, inward normal (0,+1,0)
       photon%x = amr_grid%xmin + amr_grid%xrange*rand_number()
       photon%y = amr_grid%ymin + nudge
       photon%z = amr_grid%zmin + amr_grid%zrange*rand_number()
       photon%kx = sinm*cphi;  photon%ky =  mu;  photon%kz = sinm*sphi
       photon%snx = 0.0_wp;  photon%sny =  1.0_wp;  photon%snz = 0.0_wp
    case (5)   ! +z face, inward normal (0,0,-1)
       photon%x = amr_grid%xmin + amr_grid%xrange*rand_number()
       photon%y = amr_grid%ymin + amr_grid%yrange*rand_number()
       photon%z = amr_grid%zmax - nudge
       photon%kx = sinm*cphi;  photon%ky = sinm*sphi;  photon%kz = -mu
       photon%snx = 0.0_wp;  photon%sny = 0.0_wp;  photon%snz = -1.0_wp
    case (6)   ! -z face, inward normal (0,0,+1)
       photon%x = amr_grid%xmin + amr_grid%xrange*rand_number()
       photon%y = amr_grid%ymin + amr_grid%yrange*rand_number()
       photon%z = amr_grid%zmin + nudge
       photon%kx = sinm*cphi;  photon%ky = sinm*sphi;  photon%kz =  mu
       photon%snx = 0.0_wp;  photon%sny = 0.0_wp;  photon%snz =  1.0_wp
    end select
  end subroutine emit_external_rec

  !=========================================================================
  ! Emit an isotropic external field entering a bounding sphere of radius
  ! par%rmax centered on the box.  Uniform entry point on the sphere; cosine-
  ! weighted inward direction about the inward radial -rhat.
  !=========================================================================
  subroutine emit_external_sph(photon)
    use random,     only : rand_number
    use octree_mod, only : amr_grid
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp) :: cost, sint, phi, mu, sinm, cphi, sphi, nudge
    real(kind=wp) :: xc, yc, zc, rhx, rhy, rhz
    real(kind=wp) :: tx, ty, tz, bx, by, bz
    xc = 0.5_wp*(amr_grid%xmin + amr_grid%xmax)
    yc = 0.5_wp*(amr_grid%ymin + amr_grid%ymax)
    zc = 0.5_wp*(amr_grid%zmin + amr_grid%zmax)
    cost = 2.0_wp*rand_number() - 1.0_wp    ! uniform point on the unit sphere
    sint = sqrt(1.0_wp - cost*cost)
    phi  = twopi*rand_number()
    rhx = sint*cos(phi);  rhy = sint*sin(phi);  rhz = cost   ! outward radial
    photon%x = xc + par%rmax*rhx
    photon%y = yc + par%rmax*rhy
    photon%z = zc + par%rmax*rhz
    !--- orthonormal frame (t,b) perpendicular to the inward normal -rhat.
    call ortho_frame(-rhx, -rhy, -rhz, tx, ty, tz, bx, by, bz)
    mu   = sqrt(rand_number())
    sinm = sqrt(1.0_wp - mu*mu)
    phi  = twopi*rand_number()
    cphi = cos(phi);  sphi = sin(phi)
    photon%kx = -mu*rhx + sinm*(cphi*tx + sphi*bx)
    photon%ky = -mu*rhy + sinm*(cphi*ty + sphi*by)
    photon%kz = -mu*rhz + sinm*(cphi*tz + sphi*bz)
    photon%from_external = .true.
    photon%snx = -rhx;  photon%sny = -rhy;  photon%snz = -rhz   ! inward normal
    nudge = 1.0e-6_wp * min(amr_grid%xrange, amr_grid%yrange, amr_grid%zrange)
    photon%x = photon%x - nudge*rhx      ! nudge strictly inside the sphere
    photon%y = photon%y - nudge*rhy
    photon%z = photon%z - nudge*rhz
  end subroutine emit_external_sph

  !=========================================================================
  ! Quasi-random emit helpers.  Each mirrors its legacy counterpart EXACTLY
  ! (same expressions, frames, nudges, from_external/snx..snz flags) with the
  ! Mersenne Twister draws replaced by the supplied launch uniforms.  The
  ! legacy routines are left untouched.
  !=========================================================================

  !--- isotropic internal point source at (xs,ys,zs): mu = 2 u_mu - 1,
  !--- phi = 2 pi u_phi (matches emit_point).
  subroutine emit_point_qmc(photon, xs, ys, zs, u_mu, u_phi)
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp),     intent(in)    :: xs, ys, zs, u_mu, u_phi
    real(kind=wp) :: cost, sint, phi
    photon%from_external = .false.
    photon%x = xs
    photon%y = ys
    photon%z = zs
    cost = 2.0_wp*u_mu - 1.0_wp
    sint = sqrt(1.0_wp - cost*cost)
    phi  = twopi*u_phi
    photon%kx = sint*cos(phi)
    photon%ky = sint*sin(phi)
    photon%kz = cost
  end subroutine emit_point_qmc

  !--- isotropic external field entering the box faces (matches
  !--- emit_external_rec): u_face -> area-weighted entry face, u_mu -> Lambert
  !--- incidence mu = sqrt(u_mu), u_phi -> incidence azimuth, u_c1/u_c2 -> the
  !--- two in-face coordinates (in the same positions rand_number() appeared).
  subroutine emit_external_rec_qmc(photon, u_mu, u_phi, u_face, u_c1, u_c2)
    use octree_mod, only : amr_grid
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp),     intent(in)    :: u_mu, u_phi, u_face, u_c1, u_c2
    real(kind=wp) :: mu, sinm, phi, cphi, sphi, nudge, span
    real(kind=wp) :: farea(6), fcum(6), fpick
    integer :: i, iface
    span  = min(amr_grid%xrange, amr_grid%yrange, amr_grid%zrange)
    nudge = 1.0e-6_wp * span
    farea(1) = amr_grid%yrange*amr_grid%zrange     ! +x face
    farea(2) = farea(1)                            ! -x face
    farea(3) = amr_grid%zrange*amr_grid%xrange     ! +y face
    farea(4) = farea(3)                            ! -y face
    farea(5) = amr_grid%xrange*amr_grid%yrange     ! +z face
    farea(6) = farea(5)                            ! -z face
    fcum(1) = farea(1)
    do i = 2, 6
       fcum(i) = fcum(i-1) + farea(i)
    end do
    fpick = u_face*fcum(6)
    iface = 6
    do i = 1, 6
       if (fpick <= fcum(i)) then
          iface = i;  exit
       end if
    end do
    mu   = sqrt(u_mu)                       ! cos(angle from inward normal)
    sinm = sqrt(1.0_wp - mu*mu)
    phi  = twopi*u_phi
    cphi = cos(phi);  sphi = sin(phi)
    photon%from_external = .true.
    select case (iface)
    case (1)   ! +x face, inward normal (-1,0,0)
       photon%x = amr_grid%xmax - nudge
       photon%y = amr_grid%ymin + amr_grid%yrange*u_c1
       photon%z = amr_grid%zmin + amr_grid%zrange*u_c2
       photon%kx = -mu;  photon%ky = sinm*cphi;  photon%kz = sinm*sphi
       photon%snx = -1.0_wp;  photon%sny = 0.0_wp;  photon%snz = 0.0_wp
    case (2)   ! -x face, inward normal (+1,0,0)
       photon%x = amr_grid%xmin + nudge
       photon%y = amr_grid%ymin + amr_grid%yrange*u_c1
       photon%z = amr_grid%zmin + amr_grid%zrange*u_c2
       photon%kx =  mu;  photon%ky = sinm*cphi;  photon%kz = sinm*sphi
       photon%snx =  1.0_wp;  photon%sny = 0.0_wp;  photon%snz = 0.0_wp
    case (3)   ! +y face, inward normal (0,-1,0)
       photon%x = amr_grid%xmin + amr_grid%xrange*u_c1
       photon%y = amr_grid%ymax - nudge
       photon%z = amr_grid%zmin + amr_grid%zrange*u_c2
       photon%kx = sinm*cphi;  photon%ky = -mu;  photon%kz = sinm*sphi
       photon%snx = 0.0_wp;  photon%sny = -1.0_wp;  photon%snz = 0.0_wp
    case (4)   ! -y face, inward normal (0,+1,0)
       photon%x = amr_grid%xmin + amr_grid%xrange*u_c1
       photon%y = amr_grid%ymin + nudge
       photon%z = amr_grid%zmin + amr_grid%zrange*u_c2
       photon%kx = sinm*cphi;  photon%ky =  mu;  photon%kz = sinm*sphi
       photon%snx = 0.0_wp;  photon%sny =  1.0_wp;  photon%snz = 0.0_wp
    case (5)   ! +z face, inward normal (0,0,-1)
       photon%x = amr_grid%xmin + amr_grid%xrange*u_c1
       photon%y = amr_grid%ymin + amr_grid%yrange*u_c2
       photon%z = amr_grid%zmax - nudge
       photon%kx = sinm*cphi;  photon%ky = sinm*sphi;  photon%kz = -mu
       photon%snx = 0.0_wp;  photon%sny = 0.0_wp;  photon%snz = -1.0_wp
    case (6)   ! -z face, inward normal (0,0,+1)
       photon%x = amr_grid%xmin + amr_grid%xrange*u_c1
       photon%y = amr_grid%ymin + amr_grid%yrange*u_c2
       photon%z = amr_grid%zmin + nudge
       photon%kx = sinm*cphi;  photon%ky = sinm*sphi;  photon%kz =  mu
       photon%snx = 0.0_wp;  photon%sny = 0.0_wp;  photon%snz =  1.0_wp
    end select
  end subroutine emit_external_rec_qmc

  !--- isotropic external field entering a bounding sphere (matches
  !--- emit_external_sph): u_cost/u_saz -> uniform entry point on the sphere,
  !--- u_mu -> Lambert incidence mu = sqrt(u_mu), u_phi -> incidence azimuth.
  subroutine emit_external_sph_qmc(photon, u_mu, u_phi, u_cost, u_saz)
    use octree_mod, only : amr_grid
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp),     intent(in)    :: u_mu, u_phi, u_cost, u_saz
    real(kind=wp) :: cost, sint, phi, mu, sinm, cphi, sphi, nudge
    real(kind=wp) :: xc, yc, zc, rhx, rhy, rhz
    real(kind=wp) :: tx, ty, tz, bx, by, bz
    xc = 0.5_wp*(amr_grid%xmin + amr_grid%xmax)
    yc = 0.5_wp*(amr_grid%ymin + amr_grid%ymax)
    zc = 0.5_wp*(amr_grid%zmin + amr_grid%zmax)
    cost = 2.0_wp*u_cost - 1.0_wp           ! uniform point on the unit sphere
    sint = sqrt(1.0_wp - cost*cost)
    phi  = twopi*u_saz
    rhx = sint*cos(phi);  rhy = sint*sin(phi);  rhz = cost   ! outward radial
    photon%x = xc + par%rmax*rhx
    photon%y = yc + par%rmax*rhy
    photon%z = zc + par%rmax*rhz
    !--- orthonormal frame (t,b) perpendicular to the inward normal -rhat.
    call ortho_frame(-rhx, -rhy, -rhz, tx, ty, tz, bx, by, bz)
    mu   = sqrt(u_mu)
    sinm = sqrt(1.0_wp - mu*mu)
    phi  = twopi*u_phi
    cphi = cos(phi);  sphi = sin(phi)
    photon%kx = -mu*rhx + sinm*(cphi*tx + sphi*bx)
    photon%ky = -mu*rhy + sinm*(cphi*ty + sphi*by)
    photon%kz = -mu*rhz + sinm*(cphi*tz + sphi*bz)
    photon%from_external = .true.
    photon%snx = -rhx;  photon%sny = -rhy;  photon%snz = -rhz   ! inward normal
    nudge = 1.0e-6_wp * min(amr_grid%xrange, amr_grid%yrange, amr_grid%zrange)
    photon%x = photon%x - nudge*rhx      ! nudge strictly inside the sphere
    photon%y = photon%y - nudge*rhy
    photon%z = photon%z - nudge*rhz
  end subroutine emit_external_sph_qmc

  !=========================================================================
  ! Orthonormal frame: given a unit vector n = (nx,ny,nz), return two unit
  ! vectors t, b that complete a right-handed orthonormal basis (t, b, n).
  ! The seed axis is the coordinate axis most orthogonal to n (stable).
  !=========================================================================
  subroutine ortho_frame(nx, ny, nz, tx, ty, tz, bx, by, bz)
    implicit none
    real(kind=wp), intent(in)  :: nx, ny, nz
    real(kind=wp), intent(out) :: tx, ty, tz, bx, by, bz
    real(kind=wp) :: inv
    if (abs(nx) <= abs(ny) .and. abs(nx) <= abs(nz)) then
       tx = 0.0_wp;  ty = nz;  tz = -ny          ! t = n x x_hat
    else if (abs(ny) <= abs(nz)) then
       tx = -nz;  ty = 0.0_wp;  tz = nx          ! t = n x y_hat
    else
       tx = ny;  ty = -nx;  tz = 0.0_wp          ! t = n x z_hat
    end if
    inv = 1.0_wp / sqrt(tx*tx + ty*ty + tz*tz)
    tx = tx*inv;  ty = ty*inv;  tz = tz*inv
    bx = ny*tz - nz*ty                           ! b = n x t (unit)
    by = nz*tx - nx*tz
    bz = nx*ty - ny*tx
  end subroutine ortho_frame

end module ion_band_mod
