! MoCHII: copied from MoCafe_v2.00/src/grid_mod_amr.f90 (2026-07-11)
module grid_mod_amr
!---------------------------------------------------------------------------
! AMR grid setup for MoCafe (dust-only, par%grid_type = 'amr').
!
! Reads a generic AMR file (FITS/HDF5/text) of leaf cells, builds the octree
! (octree_mod), computes the grey dust opacity of each leaf according to par%dust_model,
! and normalizes it to the system target par%taumax (radial pole) or
! par%tauhomo (volume average).  A small vestigial Cartesian box grid is then
! built (grid_create) only so the observer/output code runs unchanged --
! the transport uses amr_grid, not grid%rhokap.  Dust-only slim of LaRT's
! grid_mod_amr.f90 (no frequency grid, emissivity, velocity, or Voigt); see
! AMR_CLUMPS_PLAN.md Part A.
!
! The box is recentred on the origin (MoCafe convention: box centered at 0).
!---------------------------------------------------------------------------
  use grid_mod
  use octree_mod
  use read_generic_amr_mod
  use physics_amr_mod
  use gas_state_mod, only : gas_state_setup   ! MoCHII
  implicit none

  public :: grid_create_amr, grid_destroy_amr

contains

  !=========================================================================
  subroutine grid_create_amr(grid)
  use define
  use mpi
  implicit none
  type(grid_type), intent(inout) :: grid

  real(wp), allocatable :: xleaf(:), yleaf(:), zleaf(:)
  integer,  allocatable :: lev(:)
  real(wp), allocatable :: nH(:)
  real(wp), allocatable :: Zarr(:), xHIarr(:), ndustarr(:)
  logical  :: have_Z, have_xHI, have_ndust
  integer  :: nleaf, il, ierr
  real(wp) :: boxlen, ox, oy, oz, cxb, cyb, czb, half
  real(wp) :: Zuse, rho, taupole, tauhomo_real, vol, volsum, kapsum, opac_norm
  character(len=128) :: saved_distance_unit
  real(wp) :: saved_distance2cm

  if (trim(par%amr_type) == 'ramses') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: amr_type = ''ramses'' is not read directly; convert with '// &
        'python/AMR_grid/convert_ramses_to_generic.py first.'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if
  if (len_trim(par%amr_file) == 0) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: par%amr_file is empty (grid_type=amr).'
     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end if

  !--- read leaf data (every rank; identical read-only data).
  call generic_amr_read(trim(par%amr_file), xleaf, yleaf, zleaf, lev, &
       nH, nleaf, boxlen, &
       metallicity=Zarr, xHI=xHIarr, ndust=ndustarr, &
       origin_x=ox, origin_y=oy, origin_z=oz)
  have_Z     = allocated(Zarr)
  have_xHI   = allocated(xHIarr)
  have_ndust = allocated(ndustarr)

  !--- recenter the box on the origin: shift leaf coords by the box center.
  cxb = ox + 0.5_wp*boxlen;  cyb = oy + 0.5_wp*boxlen;  czb = oz + 0.5_wp*boxlen
  xleaf = xleaf - cxb;  yleaf = yleaf - cyb;  zleaf = zleaf - czb
  half  = 0.5_wp * boxlen
  par%rmax = half;  par%xmax = half;  par%ymax = half;  par%zmax = half
  par%xyz_symmetry = .false.;  par%xy_periodic = .false.;  par%z_symmetry = .false.

  if (mpar%p_rank == 0) then
     write(*,'(a,i12)')   ' AMR: nleaf      = ', nleaf
     write(*,'(a,f12.5)') ' AMR: boxlen     = ', boxlen
     write(*,'(a,a)')     ' AMR: dust_model = ', trim(par%dust_model)
  end if

  !--- build the octree + neighbor table (shared memory), or — for
  !--- grid_type='uniform' — the raster-ordered Cartesian grid with DDA
  !--- traversal (no tree, no neighbor table; PLAN section 4).
  if (trim(par%grid_type) == 'uniform') then
     block
       integer,  allocatable :: perm(:)
       real(wp), allocatable :: tmp(:)
       real(wp) :: dc
       integer :: nxu, ix, iy, iz, idx
       if (any(lev /= lev(1))) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: grid_type=''uniform'' needs a single-level leaf list.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       nxu = 2**lev(1)
       if (nleaf /= nxu**3) then
          if (mpar%p_rank == 0) write(*,'(a,i0,a,i0)') &
             'ERROR: uniform grid expects nleaf = (2^level)^3 = ', &
             nxu**3, ', got ', nleaf
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       !--- permutation file row -> raster slot, then reorder every
       !--- leaf array so that leaf index = raster index.
       dc = boxlen/real(nxu, wp)
       allocate(perm(nleaf), tmp(nleaf))
       perm = 0
       do il = 1, nleaf
          ix = min(max(int((xleaf(il) + half)/dc), 0), nxu-1)
          iy = min(max(int((yleaf(il) + half)/dc), 0), nxu-1)
          iz = min(max(int((zleaf(il) + half)/dc), 0), nxu-1)
          idx = 1 + ix + nxu*(iy + nxu*iz)
          if (perm(idx) /= 0) then
             if (mpar%p_rank == 0) write(*,'(a)') &
                'ERROR: uniform grid: two leaves map to one raster cell.'
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          end if
          perm(idx) = il
       end do
       tmp = xleaf;  xleaf = tmp(perm)
       tmp = yleaf;  yleaf = tmp(perm)
       tmp = zleaf;  zleaf = tmp(perm)
       tmp = nH;     nH    = tmp(perm)
       if (have_Z)     then;  tmp = Zarr;     Zarr     = tmp(perm);  end if
       if (have_xHI)   then;  tmp = xHIarr;   xHIarr   = tmp(perm);  end if
       if (have_ndust) then;  tmp = ndustarr; ndustarr = tmp(perm);  end if
       deallocate(perm, tmp)
       call amr_build_uniform(nxu, nxu, nxu, -half, half, -half, half, &
                              -half, half)
     end block
  else
     call amr_build_tree(xleaf, yleaf, zleaf, lev, nleaf, &
                         -half, half, -half, half, -half, half)
     call amr_build_neighbors
  end if
  call amr_alloc_phys()

  !--- MoCHII: persist the gas leaf state (nH, ion fractions) BEFORE the
  !--- transient read arrays are deallocated below.  MoCafe deallocates
  !--- nH/xHI after the dust-density step; the gas physics needs them.
  if (par%use_ion_band) then
     if (have_xHI) then
        call gas_state_setup(nH, nleaf, xHI=xHIarr)
     else
        call gas_state_setup(nH, nleaf)
     end if
  end if

  !--- grey dust opacity of each leaf (h_rank=0 fills the shared array).
  if (mpar%h_rank == 0) then
     do il = 1, nleaf
        select case (trim(par%dust_model))
        case ('none')   ! MoCHII: gas-only run (G0 gate); no dust opacity
           rho = 0.0_wp
        case ('laursen09_live')
           !--- MoCHII PLAN section 7 item 3: dust tied to the COMPUTED
           !--- ionization state; initial value from the initial state
           !--- (gas_state_setup ran above), refreshed each iteration in
           !--- gas_opacity_fill (with the PAH survival split there).
           !--- Uses par%Z_global (no live Z column).
           block
             use gas_state_mod, only : gas_xHI
             real(wp) :: xh, fac
             xh  = gas_xHI(il)
             fac = (1.0_wp - par%f_pah)*(xh + par%f_ion_dust*(1.0_wp - xh)) &
                   + par%f_pah*(xh + par%f_ion_pah*(1.0_wp - xh))
             rho = (par%Z_global/max(par%Z_ref,1.0e-30_wp))*nH(il)*fac &
                   * par%cext_dust * par%DGR * par%distance2cm
           end block
        case ('from_file')
           if (.not. have_ndust) then
              if (mpar%p_rank == 0) write(*,'(a)') &
                 'ERROR: dust_model=''from_file'' requires an ndust column.'
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           end if
           rho = ndustarr(il) * par%cext_dust * par%distance2cm
        case ('laursen09')
           if (.not. have_xHI) then
              if (mpar%p_rank == 0) write(*,'(a)') &
                 'ERROR: dust_model=''laursen09'' requires an xHI column '// &
                 '(T->xHI CIE is deferred with temperature).'
              call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           end if
           if (have_Z) then
              Zuse = Zarr(il)
           else
              Zuse = par%Z_global
           end if
           rho = laursen09_ndust(nH(il), xHIarr(il), Zuse, par%Z_ref, par%f_ion_dust) &
                 * par%cext_dust * par%distance2cm
        case default   ! 'global_dgr'
           rho = nH(il) * par%cext_dust * par%DGR * par%distance2cm
        end select
        amr_grid%rhokap(il) = rho
     end do
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- normalize to the system target (rank 0 computes the factor; broadcast).
  opac_norm = 1.0_wp
  if (mpar%p_rank == 0) then
     if (par%taumax > 0.0_wp) then
        taupole = amr_pole_tau()
        if (taupole > 0.0_wp) opac_norm = par%taumax / taupole
     else if (par%tauhomo > 0.0_wp) then
        volsum = 0.0_wp;  kapsum = 0.0_wp
        do il = 1, nleaf
           vol = (2.0_wp * amr_grid%ch(amr_grid%icell_of_leaf(il)))**3
           if (amr_grid%rhokap(il) > 0.0_wp) then
              volsum = volsum + vol
              kapsum = kapsum + amr_grid%rhokap(il) * vol
           end if
        end do
        if (kapsum > 0.0_wp .and. volsum > 0.0_wp) &
           opac_norm = par%tauhomo / (kapsum / volsum * half)
     end if
  end if
  call MPI_BCAST(opac_norm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (opac_norm /= 1.0_wp .and. mpar%h_rank == 0) then
     do il = 1, nleaf
        amr_grid%rhokap(il) = amr_grid%rhokap(il) * opac_norm
     end do
  end if
  call MPI_BARRIER(mpar%hostcomm, ierr)

  !--- realized system scalars (after normalization), reported to par.
  saved_distance_unit = par%distance_unit
  saved_distance2cm   = par%distance2cm
  if (par%nx < 2) par%nx = 2
  if (par%ny < 2) par%ny = 2
  if (par%nz < 2) par%nz = 2
  call grid_create(grid)                 ! vestigial box grid for observer/output
  par%distance_unit = saved_distance_unit
  par%distance2cm   = saved_distance2cm
  if (mpar%h_rank == 0) grid%rhokap(:,:,:) = 0.0_wp
  call MPI_BARRIER(mpar%hostcomm, ierr)

  if (mpar%p_rank == 0) then
     taupole      = amr_pole_tau()
     volsum = 0.0_wp;  kapsum = 0.0_wp
     do il = 1, nleaf
        vol = (2.0_wp * amr_grid%ch(amr_grid%icell_of_leaf(il)))**3
        if (amr_grid%rhokap(il) > 0.0_wp) then
           volsum = volsum + vol
           kapsum = kapsum + amr_grid%rhokap(il) * vol
        end if
     end do
     tauhomo_real = 0.0_wp
     if (volsum > 0.0_wp) tauhomo_real = kapsum / volsum * half
     par%taumax  = taupole
     par%tauhomo = tauhomo_real
     write(*,'(a,es14.5)') ' AMR derived: taumax (pole)   = ', par%taumax
     write(*,'(a,es14.5)') ' AMR derived: tauhomo (volavg)= ', par%tauhomo
     write(*,'(a,3i5)')    ' AMR box grid: nx,ny,nz = ', grid%nx, grid%ny, grid%nz
  end if
  call MPI_BCAST(par%taumax,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%tauhomo, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if (allocated(xleaf)) deallocate(xleaf, yleaf, zleaf, lev, nH)
  if (allocated(Zarr))     deallocate(Zarr)
  if (allocated(xHIarr))   deallocate(xHIarr)
  if (allocated(ndustarr)) deallocate(ndustarr)
  end subroutine grid_create_amr

  !=========================================================================
  ! Radial dust optical depth from the box center to the +z edge (the AMR
  ! analog of the Cartesian "taupole").  A small transverse offset keeps the
  ! ray off the axis-aligned cell faces.  Called on rank 0 only.
  !=========================================================================
  real(wp) function amr_pole_tau() result(tau)
  implicit none
  integer  :: il, il_new, icell, iface
  real(wp) :: x, y, z, kx, ky, kz, t_exit, off
  tau = 0.0_wp
  off = amr_grid%L_box / real(2**(amr_grid%levelmax+3), wp)   ! << smallest leaf
  x = off;  y = off;  z = 0.0_wp
  kx = 0.0_wp;  ky = 0.0_wp;  kz = 1.0_wp
  il = amr_find_leaf(x, y, z)
  do while (il > 0)
     icell = amr_grid%icell_of_leaf(il)
     call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
     tau = tau + amr_grid%rhokap(il) * t_exit
     x = x + t_exit*kx;  y = y + t_exit*ky;  z = z + t_exit*kz
     if (z >= amr_grid%zmax) exit
     il_new = amr_next_leaf(icell, iface, x, y, z)
     if (il_new <= 0) exit
     il = il_new
  end do
  end function amr_pole_tau

  !=========================================================================
  subroutine grid_destroy_amr(grid)
  use define
  implicit none
  type(grid_type), intent(inout) :: grid
  ! grid_destroy frees ALL shared-memory windows (box grid + octree); just
  ! nullify the AMR pointers afterwards.
  call grid_destroy(grid)
  nullify(amr_grid%parent, amr_grid%children, amr_grid%level, amr_grid%ileaf, &
          amr_grid%icell_of_leaf, amr_grid%cx, amr_grid%cy, amr_grid%cz, &
          amr_grid%ch, amr_grid%neighbor, amr_grid%rhokap)
  amr_grid%ncells = 0;  amr_grid%nleaf = 0
  end subroutine grid_destroy_amr

end module grid_mod_amr
