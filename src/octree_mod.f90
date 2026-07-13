! MoCHII: copied from MoCafe_v2.00/src/octree_mod.f90 (2026-07-11)
module octree_mod
!---------------------------------------------------------------------------
! AMR octree for MoCafe (dust-only), ported from LaRT_v2.00/octree_mod.f90.
!
! Tree storage is a flat array of cells (internal + leaf), 1-indexed.
! Root = cell 1.  Octant ordering: child index = 1 + ix + 2*iy + 4*iz with
! ix=1 if x >= cell_center_x else 0 (similarly iy, iz).  Tree-structure and
! the dust opacity of each leaf are MPI-3 shared memory (one copy per node).
!
! This is the dust-only slim of the LaRT octree: the Lyman-alpha physical-cell
! data (HI line opacity, Doppler width, Voigt parameter, velocity field) and
! the frequency/spectral handling are stripped.  Each leaf carries a single
! grey dust opacity rhokap(il) (extinction per unit code length).  The octree
! geometry (build, neighbor table, leaf finder, ray traversal incl. the
! pole-gap fix) is preserved verbatim.  See AMR_CLUMPS_PLAN.md Part A.
!---------------------------------------------------------------------------
  use define
  use memory_mod
  implicit none
  public

  type amr_grid_type
    ! ------ domain geometry ------
    real(wp) :: xmin = 0.0_wp, xmax = 1.0_wp, xrange = 1.0_wp
    real(wp) :: ymin = 0.0_wp, ymax = 1.0_wp, yrange = 1.0_wp
    real(wp) :: zmin = 0.0_wp, zmax = 1.0_wp, zrange = 1.0_wp
    real(wp) :: L_box = 1.0_wp

    ! ------ octree structure ------
    integer  :: ncells   = 0
    integer  :: nleaf    = 0
    integer  :: levelmax = 0

    ! Arrays of size ncells (tree topology) -- shared memory
    integer,  pointer :: parent(:)     => null()  ! parent index; 0 for root
    integer,  pointer :: children(:,:) => null()  ! (8, ncells); 0 = no child
    integer,  pointer :: level(:)      => null()  ! cell AMR level (0 = root)
    real(wp), pointer :: cx(:)         => null()  ! cell center x
    real(wp), pointer :: cy(:)         => null()  ! cell center y
    real(wp), pointer :: cz(:)         => null()  ! cell center z
    real(wp), pointer :: ch(:)         => null()  ! cell half-width

    ! Leaf-index mapping -- shared memory
    integer,  pointer :: ileaf(:)         => null()  ! ncells: >0 leaf index; 0 internal
    integer,  pointer :: icell_of_leaf(:) => null()  ! nleaf -> cell index

    ! ------ face-neighbor table (6 x ncells) -- shared memory ------
    !   iface 1=+x 2=-x 3=+y 4=-y 5=+z 6=-z; 0 => outside the box.
    integer, pointer :: neighbor(:,:) => null()

    ! ------ dust opacity of each leaf (size nleaf) -- shared memory ------
    real(wp), pointer :: rhokap(:) => null()  ! grey dust opacity per code length

    ! ------ uniform-grid mode (grid_type='uniform') ------
    ! A single-level Cartesian grid in raster order (leaf = cell =
    ! 1 + ix + nx*(iy + ny*iz), 0-based ix/iy/iz): no tree topology and
    ! no neighbor table are built; the leaf finder and the face-neighbor
    ! step become integer arithmetic (the DDA traversal — amr_cell_exit
    ! already computes the exact plane-crossing step).
    logical  :: uniform = .false.
    integer  :: nx = 0, ny = 0, nz = 0
    real(wp) :: dcell = 0.0_wp
  end type amr_grid_type

  ! Global AMR grid instance.
  type(amr_grid_type) :: amr_grid

contains

  !=========================================================================
  integer function amr_find_leaf(x, y, z) result(ileaf)
    real(wp), intent(in) :: x, y, z
    integer :: icell, ioct
    ileaf = 0
    if (x < amr_grid%xmin .or. x > amr_grid%xmax .or. &
        y < amr_grid%ymin .or. y > amr_grid%ymax .or. &
        z < amr_grid%zmin .or. z > amr_grid%zmax) return
    if (amr_grid%uniform) then
      ileaf = uni_index(x, y, z)
      return
    end if
    icell = 1
    do
      if (amr_grid%ileaf(icell) > 0) then
        ileaf = amr_grid%ileaf(icell)
        return
      end if
      ioct = 1
      if (x >= amr_grid%cx(icell)) ioct = ioct + 1
      if (y >= amr_grid%cy(icell)) ioct = ioct + 2
      if (z >= amr_grid%cz(icell)) ioct = ioct + 4
      icell = amr_grid%children(ioct, icell)
      if (icell == 0) return
    end do
  end function amr_find_leaf

  !=========================================================================
  !=========================================================================
  ! Geometry accessors: analytic in uniform mode (the cx/cy/cz/ch and
  ! leaf-map arrays are NOT allocated there — the "true" uniform
  ! backend), plain array reads on the octree.  The analytic center
  ! expressions are the SAME ones amr_build_uniform previously used to
  ! fill the arrays, so values are bit-identical.
  !=========================================================================
  pure subroutine uni_ixyz(icell, ix, iy, iz)
    integer, intent(in)  :: icell
    integer, intent(out) :: ix, iy, iz
    integer :: rem
    iz  = (icell-1)/(amr_grid%nx*amr_grid%ny)
    rem = (icell-1) - iz*amr_grid%nx*amr_grid%ny
    iy  = rem/amr_grid%nx
    ix  = rem - iy*amr_grid%nx
  end subroutine uni_ixyz

  pure real(wp) function cell_cx(icell)
    integer, intent(in) :: icell
    integer :: ix, iy, iz
    if (amr_grid%uniform) then
       call uni_ixyz(icell, ix, iy, iz)
       cell_cx = amr_grid%xmin + (real(ix,wp) + 0.5_wp)*amr_grid%dcell
    else
       cell_cx = amr_grid%cx(icell)
    end if
  end function cell_cx

  pure real(wp) function cell_cy(icell)
    integer, intent(in) :: icell
    integer :: ix, iy, iz
    if (amr_grid%uniform) then
       call uni_ixyz(icell, ix, iy, iz)
       cell_cy = amr_grid%ymin + (real(iy,wp) + 0.5_wp)*amr_grid%dcell
    else
       cell_cy = amr_grid%cy(icell)
    end if
  end function cell_cy

  pure real(wp) function cell_cz(icell)
    integer, intent(in) :: icell
    integer :: ix, iy, iz
    if (amr_grid%uniform) then
       call uni_ixyz(icell, ix, iy, iz)
       cell_cz = amr_grid%zmin + (real(iz,wp) + 0.5_wp)*amr_grid%dcell
    else
       cell_cz = amr_grid%cz(icell)
    end if
  end function cell_cz

  pure real(wp) function cell_ch(icell)
    integer, intent(in) :: icell
    if (amr_grid%uniform) then
       cell_ch = 0.5_wp*amr_grid%dcell
    else
       cell_ch = amr_grid%ch(icell)
    end if
  end function cell_ch

  pure integer function leaf_cell(il)      ! leaf index -> cell index
    integer, intent(in) :: il
    if (amr_grid%uniform) then
       leaf_cell = il
    else
       leaf_cell = amr_grid%icell_of_leaf(il)
    end if
  end function leaf_cell

  !--- leaf-level conveniences (the idioms the gas modules use)
  pure real(wp) function leaf_cx(il)
    integer, intent(in) :: il
    leaf_cx = cell_cx(leaf_cell(il))
  end function leaf_cx
  pure real(wp) function leaf_cy(il)
    integer, intent(in) :: il
    leaf_cy = cell_cy(leaf_cell(il))
  end function leaf_cy
  pure real(wp) function leaf_cz(il)
    integer, intent(in) :: il
    leaf_cz = cell_cz(leaf_cell(il))
  end function leaf_cz
  pure real(wp) function leaf_half(il)     ! cell half-width of a leaf
    integer, intent(in) :: il
    leaf_half = cell_ch(leaf_cell(il))
  end function leaf_half

  !=========================================================================
  ! Uniform-grid raster index of (x,y,z); 1-based, 0 outside.
  !=========================================================================
  integer function uni_index(x, y, z) result(il)
    real(wp), intent(in) :: x, y, z
    integer :: ix, iy, iz
    ix = int((x - amr_grid%xmin)/amr_grid%dcell)
    iy = int((y - amr_grid%ymin)/amr_grid%dcell)
    iz = int((z - amr_grid%zmin)/amr_grid%dcell)
    ix = min(max(ix, 0), amr_grid%nx-1)
    iy = min(max(iy, 0), amr_grid%ny-1)
    iz = min(max(iz, 0), amr_grid%nz-1)
    il = 1 + ix + amr_grid%nx*(iy + amr_grid%ny*iz)
  end function uni_index

  !=========================================================================
  integer function amr_find_cell(x, y, z) result(icell_out)
    real(wp), intent(in) :: x, y, z
    integer :: icell, ioct
    icell_out = 0
    if (x < amr_grid%xmin .or. x > amr_grid%xmax .or. &
        y < amr_grid%ymin .or. y > amr_grid%ymax .or. &
        z < amr_grid%zmin .or. z > amr_grid%zmax) return
    if (amr_grid%uniform) then
      icell_out = uni_index(x, y, z)
      return
    end if
    icell = 1
    do
      if (amr_grid%ileaf(icell) > 0) then
        icell_out = icell
        return
      end if
      ioct = 1
      if (x >= amr_grid%cx(icell)) ioct = ioct + 1
      if (y >= amr_grid%cy(icell)) ioct = ioct + 2
      if (z >= amr_grid%cz(icell)) ioct = ioct + 4
      icell = amr_grid%children(ioct, icell)
      if (icell == 0) then
        icell_out = 0
        return
      end if
    end do
  end function amr_find_cell

  !=========================================================================
  ! Deepest cell (leaf or internal) containing (x,y,z); returns the internal
  ! gap-cell (not 0) when descent reaches an octant with no child.
  !=========================================================================
  integer function amr_find_enclosing_cell(x, y, z) result(icell_out)
    real(wp), intent(in) :: x, y, z
    integer :: icell, ioct, child
    icell_out = 0
    if (x < amr_grid%xmin .or. x > amr_grid%xmax .or. &
        y < amr_grid%ymin .or. y > amr_grid%ymax .or. &
        z < amr_grid%zmin .or. z > amr_grid%zmax) return
    if (amr_grid%uniform) then
      icell_out = uni_index(x, y, z)
      return
    end if
    icell = 1
    do
      if (amr_grid%ileaf(icell) > 0) then
        icell_out = icell
        return
      end if
      ioct = 1
      if (x >= amr_grid%cx(icell)) ioct = ioct + 1
      if (y >= amr_grid%cy(icell)) ioct = ioct + 2
      if (z >= amr_grid%cz(icell)) ioct = ioct + 4
      child = amr_grid%children(ioct, icell)
      if (child == 0) then
        icell_out = icell
        return
      end if
      icell = child
    end do
  end function amr_find_enclosing_cell

  !=========================================================================
  ! Step through a virtual sub-cell (half the size of icell) centered in the
  ! sub-octant containing (x,y,z).  Pole-gap traversal helper.
  !=========================================================================
  subroutine amr_gap_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
    integer,  intent(in)  :: icell
    real(wp), intent(in)  :: x, y, z, kx, ky, kz
    real(wp), intent(out) :: t_exit
    integer,  intent(out) :: iface
    real(wp) :: cx, cy, cz, h_sub, cx_sub, cy_sub, cz_sub, t(6)
    integer  :: ix, iy, iz
    cx = amr_grid%cx(icell);  cy = amr_grid%cy(icell);  cz = amr_grid%cz(icell)
    h_sub = amr_grid%ch(icell) * 0.5_wp
    ix = 0;  if (x >= cx) ix = 1
    iy = 0;  if (y >= cy) iy = 1
    iz = 0;  if (z >= cz) iz = 1
    cx_sub = cx + real(2*ix - 1, wp) * h_sub
    cy_sub = cy + real(2*iy - 1, wp) * h_sub
    cz_sub = cz + real(2*iz - 1, wp) * h_sub
    if (kx > 0.0_wp) then
      t(1) = (cx_sub + h_sub - x) / kx;  t(2) = hugest
    else if (kx < 0.0_wp) then
      t(1) = hugest;  t(2) = (cx_sub - h_sub - x) / kx
    else
      t(1) = hugest;  t(2) = hugest
    end if
    if (ky > 0.0_wp) then
      t(3) = (cy_sub + h_sub - y) / ky;  t(4) = hugest
    else if (ky < 0.0_wp) then
      t(3) = hugest;  t(4) = (cy_sub - h_sub - y) / ky
    else
      t(3) = hugest;  t(4) = hugest
    end if
    if (kz > 0.0_wp) then
      t(5) = (cz_sub + h_sub - z) / kz;  t(6) = hugest
    else if (kz < 0.0_wp) then
      t(5) = hugest;  t(6) = (cz_sub - h_sub - z) / kz
    else
      t(5) = hugest;  t(6) = hugest
    end if
    iface  = minloc(t, dim=1)
    t_exit = t(iface)
  end subroutine amr_gap_exit

  !=========================================================================
  ! Distance from (x,y,z) along (kx,ky,kz) to the nearest face of cell icell.
  ! iface: 1=+x, 2=-x, 3=+y, 4=-y, 5=+z, 6=-z.
  !=========================================================================
  subroutine amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
    integer,  intent(in)  :: icell
    real(wp), intent(in)  :: x, y, z, kx, ky, kz
    real(wp), intent(out) :: t_exit
    integer,  intent(out) :: iface
    real(wp) :: h, cx, cy, cz, t(6)
    cx = cell_cx(icell);  cy = cell_cy(icell);  cz = cell_cz(icell)
    h  = cell_ch(icell)
    if (kx > 0.0_wp) then
      t(1) = (cx + h - x) / kx;  t(2) = hugest
    else if (kx < 0.0_wp) then
      t(1) = hugest;             t(2) = (cx - h - x) / kx
    else
      t(1) = hugest;             t(2) = hugest
    end if
    if (ky > 0.0_wp) then
      t(3) = (cy + h - y) / ky;  t(4) = hugest
    else if (ky < 0.0_wp) then
      t(3) = hugest;             t(4) = (cy - h - y) / ky
    else
      t(3) = hugest;             t(4) = hugest
    end if
    if (kz > 0.0_wp) then
      t(5) = (cz + h - z) / kz;  t(6) = hugest
    else if (kz < 0.0_wp) then
      t(5) = hugest;             t(6) = (cz - h - z) / kz
    else
      t(5) = hugest;             t(6) = hugest
    end if
    iface  = minloc(t, dim=1)
    t_exit = t(iface)
  end subroutine amr_cell_exit

  !=========================================================================
  ! Build an octree from flat arrays of leaf positions and AMR levels.
  ! Built locally on every rank (identical broadcast data), then transferred
  ! to MPI-3 shared memory (one physical copy per node; h_rank==0 writes).
  !=========================================================================
  subroutine amr_build_tree(xleaf, yleaf, zleaf, leaf_level, nleaf, &
                             box_xmin, box_xmax, box_ymin, box_ymax, box_zmin, box_zmax)
    use mpi
    integer,  intent(in) :: nleaf
    real(wp), intent(in) :: xleaf(nleaf), yleaf(nleaf), zleaf(nleaf)
    integer,  intent(in) :: leaf_level(nleaf)
    real(wp), intent(in) :: box_xmin, box_xmax, box_ymin, box_ymax, box_zmin, box_zmax

    integer,  allocatable :: t_parent(:), t_children(:,:), t_level(:), t_ileaf(:)
    real(wp), allocatable :: t_cx(:), t_cy(:), t_cz(:), t_ch(:)
    integer,  allocatable :: t_icell_of_leaf(:)
    integer  :: ncells_max, ncells
    integer  :: il, icell, lev, ioct, nc, ix, iy, iz, ierr

    ncells_max = max(nleaf * 2, 16)
    allocate(t_parent(ncells_max),     source=0)
    allocate(t_children(8,ncells_max), source=0)
    allocate(t_level(ncells_max),      source=0)
    allocate(t_ileaf(ncells_max),      source=0)
    allocate(t_cx(ncells_max), source=0.0_wp)
    allocate(t_cy(ncells_max), source=0.0_wp)
    allocate(t_cz(ncells_max), source=0.0_wp)
    allocate(t_ch(ncells_max), source=0.0_wp)
    allocate(t_icell_of_leaf(nleaf), source=0)

    amr_grid%xmin = box_xmin;  amr_grid%xmax = box_xmax;  amr_grid%xrange = box_xmax - box_xmin
    amr_grid%ymin = box_ymin;  amr_grid%ymax = box_ymax;  amr_grid%yrange = box_ymax - box_ymin
    amr_grid%zmin = box_zmin;  amr_grid%zmax = box_zmax;  amr_grid%zrange = box_zmax - box_zmin
    amr_grid%L_box = box_xmax - box_xmin

    ncells          = 1
    t_parent(1)     = 0
    t_children(:,1) = 0
    t_level(1)      = 0
    t_ileaf(1)      = 0
    t_cx(1) = (box_xmin + box_xmax) * 0.5_wp
    t_cy(1) = (box_ymin + box_ymax) * 0.5_wp
    t_cz(1) = (box_zmin + box_zmax) * 0.5_wp
    t_ch(1) = amr_grid%L_box * 0.5_wp

    do il = 1, nleaf
      icell = 1
      do lev = 0, leaf_level(il) - 1
        ioct = 1
        if (xleaf(il) >= t_cx(icell)) ioct = ioct + 1
        if (yleaf(il) >= t_cy(icell)) ioct = ioct + 2
        if (zleaf(il) >= t_cz(icell)) ioct = ioct + 4
        if (t_children(ioct, icell) == 0) then
          ncells = ncells + 1
          if (ncells > ncells_max) then
            call grow_local_arrays(ncells_max*2, ncells_max, &
                t_parent, t_children, t_level, t_ileaf, t_cx, t_cy, t_cz, t_ch)
            ncells_max = size(t_parent)
          end if
          nc = ncells
          t_parent(nc)      = icell
          t_children(:, nc) = 0
          t_level(nc)       = lev + 1
          t_ileaf(nc)       = 0
          t_ch(nc)          = t_ch(icell) * 0.5_wp
          ix = mod(ioct - 1, 2)
          iy = mod((ioct - 1) / 2, 2)
          iz = (ioct - 1) / 4
          t_cx(nc) = t_cx(icell) + real(2*ix - 1, wp) * t_ch(nc)
          t_cy(nc) = t_cy(icell) + real(2*iy - 1, wp) * t_ch(nc)
          t_cz(nc) = t_cz(icell) + real(2*iz - 1, wp) * t_ch(nc)
          t_children(ioct, icell) = nc
        end if
        icell = t_children(ioct, icell)
      end do
      t_ileaf(icell)      = il
      t_icell_of_leaf(il) = icell
    end do

    amr_grid%ncells   = ncells
    amr_grid%nleaf    = nleaf
    amr_grid%levelmax = maxval(leaf_level)

    call create_shared_mem(amr_grid%parent,        [ncells])
    call create_shared_mem(amr_grid%children,      [8, ncells])
    call create_shared_mem(amr_grid%level,         [ncells])
    call create_shared_mem(amr_grid%ileaf,         [ncells])
    call create_shared_mem(amr_grid%cx,            [ncells])
    call create_shared_mem(amr_grid%cy,            [ncells])
    call create_shared_mem(amr_grid%cz,            [ncells])
    call create_shared_mem(amr_grid%ch,            [ncells])
    call create_shared_mem(amr_grid%icell_of_leaf, [nleaf])

    if (mpar%h_rank == 0) then
      amr_grid%parent(1:ncells)       = t_parent(1:ncells)
      amr_grid%children(:, 1:ncells)  = t_children(:, 1:ncells)
      amr_grid%level(1:ncells)        = t_level(1:ncells)
      amr_grid%ileaf(1:ncells)        = t_ileaf(1:ncells)
      amr_grid%cx(1:ncells)           = t_cx(1:ncells)
      amr_grid%cy(1:ncells)           = t_cy(1:ncells)
      amr_grid%cz(1:ncells)           = t_cz(1:ncells)
      amr_grid%ch(1:ncells)           = t_ch(1:ncells)
      amr_grid%icell_of_leaf(1:nleaf) = t_icell_of_leaf(1:nleaf)
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)

    deallocate(t_parent, t_children, t_level, t_ileaf, &
               t_cx, t_cy, t_cz, t_ch, t_icell_of_leaf)
  end subroutine amr_build_tree

  !=========================================================================
  ! Allocate the dust opacity of each leaf as shared memory.  h_rank==0 fills it
  ! and calls MPI_BARRIER(mpar%hostcomm) after this returns.
  !=========================================================================
  subroutine amr_alloc_phys()
    call create_shared_mem(amr_grid%rhokap, [amr_grid%nleaf])
  end subroutine amr_alloc_phys

  !=========================================================================
  ! Precompute the face-neighbor table neighbor(6, ncells) as shared memory.
  ! Must be called after amr_build_tree.
  !=========================================================================
  subroutine amr_build_neighbors
    use mpi
    integer  :: icell, ierr, iface, ineigh
    real(wp) :: cx, cy, cz, h, hp

    call create_shared_mem(amr_grid%neighbor, [6, amr_grid%ncells])

    if (mpar%h_rank == 0) then
      amr_grid%neighbor = 0
      do icell = 1, amr_grid%ncells
        cx = amr_grid%cx(icell);  cy = amr_grid%cy(icell);  cz = amr_grid%cz(icell)
        h  = amr_grid%ch(icell)
        hp = 2.0_wp * h   ! query the neighbor center (1 cell width past the face)
        if (cx + hp <= amr_grid%xmax) &
          amr_grid%neighbor(1, icell) = amr_find_cell_at_level(cx+hp, cy,    cz,    amr_grid%level(icell))
        if (cx - hp >= amr_grid%xmin) &
          amr_grid%neighbor(2, icell) = amr_find_cell_at_level(cx-hp, cy,    cz,    amr_grid%level(icell))
        if (cy + hp <= amr_grid%ymax) &
          amr_grid%neighbor(3, icell) = amr_find_cell_at_level(cx,    cy+hp, cz,    amr_grid%level(icell))
        if (cy - hp >= amr_grid%ymin) &
          amr_grid%neighbor(4, icell) = amr_find_cell_at_level(cx,    cy-hp, cz,    amr_grid%level(icell))
        if (cz + hp <= amr_grid%zmax) &
          amr_grid%neighbor(5, icell) = amr_find_cell_at_level(cx,    cy,    cz+hp, amr_grid%level(icell))
        if (cz - hp >= amr_grid%zmin) &
          amr_grid%neighbor(6, icell) = amr_find_cell_at_level(cx,    cy,    cz-hp, amr_grid%level(icell))
        do iface = 1, 6
          ineigh = amr_grid%neighbor(iface, icell)
          if (ineigh > 0 .and. ineigh /= icell .and. is_ancestor(ineigh, icell)) &
            amr_grid%neighbor(iface, icell) = 0
        end do
      end do
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine amr_build_neighbors

  pure logical function is_ancestor(anc, desc)
    integer, intent(in) :: anc, desc
    integer :: c
    is_ancestor = .false.
    c = desc
    do while (c > 0)
      c = amr_grid%parent(c)
      if (c == anc) then
        is_ancestor = .true.
        return
      end if
    end do
  end function is_ancestor

  !=========================================================================
  ! Leaf entered after crossing face iface of cell icell.  (x,y,z) is the
  ! photon position AT the face; the iface-normal sub-octant bit is set
  ! topologically (round-off safe).  Returns 0 if the photon left the box.
  !=========================================================================
  integer function amr_next_leaf(icell, iface, x, y, z) result(il_new)
    integer,  intent(in) :: icell, iface
    real(wp), intent(in) :: x, y, z
    integer :: ineigh, child, ioct
    il_new = 0
    if (amr_grid%uniform) then
      il_new = uni_next(icell, iface)
      return
    end if
    ineigh = amr_grid%neighbor(iface, icell)
    if (ineigh == 0) return
    do while (amr_grid%ileaf(ineigh) == 0)
      ioct = 1
      select case (iface)
      case (1)
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (2)
        ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (3)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (4)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (5)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
      case (6)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        ioct = ioct + 4
      end select
      child = amr_grid%children(ioct, ineigh)
      if (child == 0) exit
      ineigh = child
    end do
    il_new = amr_grid%ileaf(ineigh)
  end function amr_next_leaf

  !=========================================================================
  ! Like amr_next_leaf but distinguishes "outside box" from "in a gap".
  !   il_new > 0                       : valid leaf
  !   il_new = 0, icell_gap > 0        : in a gap (h_gap = gap-cell half-width)
  !   il_new = 0, icell_gap = 0        : outside the box
  !=========================================================================
  subroutine amr_next_leaf_or_gap(icell, iface, x, y, z, il_new, icell_gap, h_gap)
    integer,  intent(in)  :: icell, iface
    real(wp), intent(in)  :: x, y, z
    integer,  intent(out) :: il_new, icell_gap
    real(wp), intent(out) :: h_gap
    integer :: ineigh, child, ioct
    il_new = 0;  icell_gap = 0;  h_gap = 0.0_wp
    if (amr_grid%uniform) then
      il_new = uni_next(icell, iface)
      return
    end if
    ineigh = amr_grid%neighbor(iface, icell)
    if (ineigh == 0) return
    do while (amr_grid%ileaf(ineigh) == 0)
      ioct = 1
      select case (iface)
      case (1)
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (2)
        ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (3)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (4)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        ioct = ioct + 2
        if (z >= amr_grid%cz(ineigh)) ioct = ioct + 4
      case (5)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
      case (6)
        if (x >= amr_grid%cx(ineigh)) ioct = ioct + 1
        if (y >= amr_grid%cy(ineigh)) ioct = ioct + 2
        ioct = ioct + 4
      end select
      child = amr_grid%children(ioct, ineigh)
      if (child == 0) then
        icell_gap = ineigh
        h_gap     = amr_grid%ch(ineigh)
        return
      end if
      ineigh = child
    end do
    il_new = amr_grid%ileaf(ineigh)
  end subroutine amr_next_leaf_or_gap

  !=========================================================================
  ! Uniform-grid face-neighbor step: pure integer arithmetic (the DDA
  ! traversal step; no neighbor table).  Returns 0 outside the box.
  !=========================================================================
  integer function uni_next(icell, iface) result(il_new)
    integer, intent(in) :: icell, iface
    integer :: ix, iy, iz, rem
    iz  = (icell-1)/(amr_grid%nx*amr_grid%ny)
    rem = (icell-1) - iz*amr_grid%nx*amr_grid%ny
    iy  = rem/amr_grid%nx
    ix  = rem - iy*amr_grid%nx
    select case (iface)
    case (1);  ix = ix + 1
    case (2);  ix = ix - 1
    case (3);  iy = iy + 1
    case (4);  iy = iy - 1
    case (5);  iz = iz + 1
    case (6);  iz = iz - 1
    end select
    if (ix < 0 .or. ix >= amr_grid%nx .or. iy < 0 .or. iy >= amr_grid%ny &
        .or. iz < 0 .or. iz >= amr_grid%nz) then
       il_new = 0
    else
       il_new = 1 + ix + amr_grid%nx*(iy + amr_grid%ny*iz)
    end if
  end function uni_next

  !=========================================================================
  ! Build the uniform-grid geometry (grid_type='uniform'): cx/cy/cz/ch and
  ! identity leaf maps in raster order; no tree topology, no neighbor
  ! table.  The caller must pass the state arrays already raster-ordered.
  !=========================================================================
  subroutine amr_build_uniform(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax)
    use memory_mod, only : create_shared_mem
    use mpi
    implicit none
    integer,  intent(in) :: nx, ny, nz
    real(wp), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
    integer :: il, ix, iy, iz, ierr

    amr_grid%uniform = .true.
    amr_grid%nx = nx;  amr_grid%ny = ny;  amr_grid%nz = nz
    amr_grid%xmin = xmin;  amr_grid%xmax = xmax;  amr_grid%xrange = xmax - xmin
    amr_grid%ymin = ymin;  amr_grid%ymax = ymax;  amr_grid%yrange = ymax - ymin
    amr_grid%zmin = zmin;  amr_grid%zmax = zmax;  amr_grid%zrange = zmax - zmin
    amr_grid%L_box = xmax - xmin
    amr_grid%dcell = (xmax - xmin)/real(nx, wp)
    amr_grid%ncells = nx*ny*nz
    amr_grid%nleaf  = nx*ny*nz
    amr_grid%levelmax = nint(log(real(nx,wp))/log(2.0_wp))

    !--- the true uniform backend: NO geometry arrays — centers, half
    !--- widths, and the leaf<->cell maps come analytically from the
    !--- raster index (cell_cx/cy/cz/ch, leaf_cell accessors).
    call MPI_BARRIER(mpar%hostcomm, ierr)
    if (mpar%p_rank == 0) write(*,'(a,3i5,a)') &
       ' UNI: uniform grid ', nx, ny, nz, &
       ' (raster order, analytic geometry — no tree arrays)'
  end subroutine amr_build_uniform

  ! ---- private helpers ----

  integer function amr_find_cell_at_level(x, y, z, target_level) result(icell)
    real(wp), intent(in) :: x, y, z
    integer,  intent(in) :: target_level
    integer :: ioct, child
    icell = 0
    if (x < amr_grid%xmin .or. x > amr_grid%xmax .or. &
        y < amr_grid%ymin .or. y > amr_grid%ymax .or. &
        z < amr_grid%zmin .or. z > amr_grid%zmax) return
    icell = 1
    do
      if (amr_grid%level(icell) >= target_level) return
      if (amr_grid%ileaf(icell) > 0)             return
      ioct = 1
      if (x >= amr_grid%cx(icell)) ioct = ioct + 1
      if (y >= amr_grid%cy(icell)) ioct = ioct + 2
      if (z >= amr_grid%cz(icell)) ioct = ioct + 4
      child = amr_grid%children(ioct, icell)
      if (child == 0) return
      icell = child
    end do
  end function amr_find_cell_at_level

  subroutine grow_local_arrays(new_max, old_max, &
      t_par, t_chi, t_lev, t_ile, t_cx, t_cy, t_cz, t_ch)
    integer,  intent(in) :: new_max, old_max
    integer,  allocatable, intent(inout) :: t_par(:), t_chi(:,:), t_lev(:), t_ile(:)
    real(wp), allocatable, intent(inout) :: t_cx(:), t_cy(:), t_cz(:), t_ch(:)
    integer,  allocatable :: tmp_i1(:), tmp_i2(:,:)
    real(wp), allocatable :: tmp_r(:)

    allocate(tmp_i1(new_max));  tmp_i1(1:old_max) = t_par;  tmp_i1(old_max+1:) = 0
    call move_alloc(tmp_i1, t_par)
    allocate(tmp_i2(8, new_max));  tmp_i2(:, 1:old_max) = t_chi;  tmp_i2(:, old_max+1:) = 0
    call move_alloc(tmp_i2, t_chi)
    allocate(tmp_i1(new_max));  tmp_i1(1:old_max) = t_lev;  tmp_i1(old_max+1:) = 0
    call move_alloc(tmp_i1, t_lev)
    allocate(tmp_i1(new_max));  tmp_i1(1:old_max) = t_ile;  tmp_i1(old_max+1:) = 0
    call move_alloc(tmp_i1, t_ile)
    allocate(tmp_r(new_max));  tmp_r(1:old_max) = t_cx;  tmp_r(old_max+1:) = 0.0_wp
    call move_alloc(tmp_r, t_cx)
    allocate(tmp_r(new_max));  tmp_r(1:old_max) = t_cy;  tmp_r(old_max+1:) = 0.0_wp
    call move_alloc(tmp_r, t_cy)
    allocate(tmp_r(new_max));  tmp_r(1:old_max) = t_cz;  tmp_r(old_max+1:) = 0.0_wp
    call move_alloc(tmp_r, t_cz)
    allocate(tmp_r(new_max));  tmp_r(1:old_max) = t_ch;  tmp_r(old_max+1:) = 0.0_wp
    call move_alloc(tmp_r, t_ch)
  end subroutine grow_local_arrays

end module octree_mod
