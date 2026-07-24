! MoCHII: copied from MoCafe_v2.00/src/grid_mod_amr.f90 (2026-07-11)
module grid_mod_amr
!---------------------------------------------------------------------------
! AMR grid setup (dust-only, par%grid_type = 'amr').
!
! Reads a generic AMR file (FITS/HDF5/text) of leaf cells, builds the octree
! (octree_mod), computes the grey dust opacity of each leaf according to par%dust_model,
! and normalizes it to the system target par%taumax (radial pole) or
! par%tauhomo (volume average).  A small vestigial Cartesian box grid is then
! built (grid_create) only so the observer/output code runs unchanged --
! the transport uses amr_grid, not grid%rhokap.
!
! The box is recentred on the origin (box centered at 0).
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
  use read_mod, only : get_dimension, read_3D
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

  if (len_trim(par%amr_file) > 0) then
     !--- read leaf data from a file (every rank; identical read-only data).
     call generic_amr_read(trim(par%amr_file), xleaf, yleaf, zleaf, lev, &
          nH, nleaf, boxlen, &
          metallicity=Zarr, xHI=xHIarr, ndust=ndustarr, &
          origin_x=ox, origin_y=oy, origin_z=oz)
     !--- recenter the box on the origin: shift leaf coords by the box center.
     cxb = ox + 0.5_wp*boxlen;  cyb = oy + 0.5_wp*boxlen;  czb = oz + 0.5_wp*boxlen
     xleaf = xleaf - cxb;  yleaf = yleaf - cyb;  zleaf = zleaf - czb
     half  = 0.5_wp * boxlen
     par%xmax = half;  par%ymax = half;  par%zmax = half
  else
     !--- MoCHII: build a single-level Cartesian ('car') grid.  With
     !--- par%density_file set, nx/ny/nz and the per-leaf nH come from a 3D
     !--- density cube (FITS/HDF5, NAXIS1/2/3 = nx/ny/nz, values = nH [cm^-3]);
     !--- otherwise nx/ny/nz are from the namelist and nH is uniform (set by
     !--- the density model below).  The box is par%xmax/ymax/zmax; the car
     !--- traversal uses one cell size, so the cells must be cubic.
     block
       real(wp) :: dcx, dcy, dcz
       integer  :: ix, iy, iz, dstat
       real(wp), allocatable :: cube(:,:,:)
       if (len_trim(par%density_file) > 0) then
          dstat = 0
          call get_dimension(trim(par%density_file), par%nx, par%ny, par%nz, &
                             dstat, reduce_factor=par%reduce_factor)
          if (dstat /= 0) then
             if (mpar%p_rank == 0) write(*,'(3a)') &
                'ERROR: cannot read cube dimensions from par%density_file (', &
                trim(par%density_file), ').'
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          end if
       end if
       dcx = 2.0_wp*par%xmax/real(par%nx, wp)
       dcy = 2.0_wp*par%ymax/real(par%ny, wp)
       dcz = 2.0_wp*par%zmax/real(par%nz, wp)
       if (abs(dcy-dcx) > 1.0e-9_wp*dcx .or. abs(dcz-dcx) > 1.0e-9_wp*dcx) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: car grid needs cubic cells: 2*xmax/nx = 2*ymax/ny = 2*zmax/nz.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       nleaf  = par%nx*par%ny*par%nz
       boxlen = 2.0_wp*par%xmax
       half   = par%xmax
       allocate(xleaf(nleaf), yleaf(nleaf), zleaf(nleaf), lev(nleaf), nH(nleaf))
       lev = 0;  nH = 1.0_wp
       il = 0
       do iz = 0, par%nz-1
          do iy = 0, par%ny-1
             do ix = 0, par%nx-1
                il = il + 1     ! = 1 + ix + nx*(iy + ny*iz): raster order
                xleaf(il) = -par%xmax + (real(ix,wp)+0.5_wp)*dcx
                yleaf(il) = -par%ymax + (real(iy,wp)+0.5_wp)*dcy
                zleaf(il) = -par%zmax + (real(iz,wp)+0.5_wp)*dcz
             end do
          end do
       end do
       if (len_trim(par%density_file) > 0) then
          allocate(cube(par%nx, par%ny, par%nz))
          call read_3D(trim(par%density_file), cube, &
                       reduce_factor=par%reduce_factor, centering=par%centering)
          il = 0
          do iz = 0, par%nz-1
             do iy = 0, par%ny-1
                do ix = 0, par%nx-1
                   il = il + 1                        ! same raster order (x fastest)
                   nH(il) = cube(ix+1, iy+1, iz+1)
                end do
             end do
          end do
          deallocate(cube)
          if (mpar%p_rank == 0) write(*,'(3a,3(i0,a),es10.3,a,es10.3)') &
             ' GRID: nH from density cube ', trim(par%density_file), &
             ' [', par%nx, 'x', par%ny, 'x', par%nz, '], nH min/max = ', &
             minval(nH), ' /', maxval(nH)
       end if
     end block
  end if
  have_Z     = allocated(Zarr)
  have_xHI   = allocated(xHIarr)
  have_ndust = allocated(ndustarr)
  !--- honor the namelist par%xy_periodic (the slab boundary condition); the
  !--- reflecting symmetries are not used by MoCHII grids and stay off.
  par%xyz_symmetry = .false.;  par%z_symmetry = .false.

  !--- MoCHII gas-density model (both 'car' and 'amr'): a uniform density
  !--- (par%nH_const) and/or a spherical sphere/shell cut (par%rmax/rmin) let
  !--- a sharp geometry be resolved by the octree while the density stays
  !--- uniform.  Applied on the leaf centers before the grid is built.
  if (par%nH_const >= 0.0_wp .and. len_trim(par%density_file) == 0) nH = par%nH_const
  if (par%rmax > 0.0_wp .or. par%rmin > 0.0_wp) then
     do il = 1, nleaf
        rho = sqrt(xleaf(il)**2 + yleaf(il)**2 + zleaf(il)**2)
        if (par%rmax > 0.0_wp .and. rho > par%rmax) nH(il) = 0.0_wp
        if (par%rmin > 0.0_wp .and. rho < par%rmin) nH(il) = 0.0_wp
     end do
  end if

  if (mpar%p_rank == 0) then
     write(*,'(a,a)')     ' GRID: type       = ', trim(par%grid_type)
     write(*,'(a,i12)')   ' GRID: nleaf      = ', nleaf
     write(*,'(a,f12.5)') ' GRID: boxlen     = ', boxlen
     write(*,'(a,a)')     ' GRID: dust_model = ', trim(par%dust_model)
     if (par%nH_const >= 0.0_wp) write(*,'(a,es12.4)') ' GRID: nH_const   = ', par%nH_const
     if (par%rmax > 0.0_wp)      write(*,'(a,f12.5)')  ' GRID: rmax       = ', par%rmax
     if (par%rmin > 0.0_wp)      write(*,'(a,f12.5)')  ' GRID: rmin       = ', par%rmin
  end if

  !--- build the octree + neighbor table (shared memory), or — for
  !--- grid_type='car' — the raster-ordered Cartesian grid with DDA
  !--- traversal (no tree, no neighbor table).
  if (trim(par%grid_type) == 'car') then
     if (len_trim(par%amr_file) > 0) then
        !--- file 'car' grid: a single-level octree-style leaf list, permuted
        !--- from file row to raster slot so that leaf index = raster index.
        block
          integer,  allocatable :: perm(:)
          real(wp), allocatable :: tmp(:)
          real(wp) :: dc
          integer :: nxu, ix, iy, iz, idx
          if (any(lev /= lev(1))) then
             if (mpar%p_rank == 0) write(*,'(a)') &
                'ERROR: grid_type=''car'' from a file needs a single-level leaf list.'
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          end if
          nxu = 2**lev(1)
          if (nleaf /= nxu**3) then
             if (mpar%p_rank == 0) write(*,'(a,i0,a,i0)') &
                'ERROR: car grid from a file expects nleaf = (2^level)^3 = ', &
                nxu**3, ', got ', nleaf
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          end if
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
                   'ERROR: car grid: two leaves map to one raster cell.'
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
          call amr_build_car(nxu, nxu, nxu, -half, half, -half, half, &
                             -half, half)
        end block
     else
        !--- namelist 'car' grid: already in raster order.
        call amr_build_car(par%nx, par%ny, par%nz, &
                           -par%xmax, par%xmax, -par%ymax, par%ymax, &
                           -par%zmax, par%zmax)
     end if
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
        case ('none')   ! MoCHII: gas-only run; no dust opacity
           rho = 0.0_wp
        case ('laursen09_live')
           !--- MoCHII: dust tied to the COMPUTED
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
           vol = (2.0_wp * leaf_half(il))**3
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
        vol = (2.0_wp * leaf_half(il))**3
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
     icell = leaf_cell(il)
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
