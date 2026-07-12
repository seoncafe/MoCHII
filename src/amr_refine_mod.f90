module amr_refine_mod
!---------------------------------------------------------------------------
! MoCHII: solution-driven octree re-refinement (Stage G4).
!
! The I-front is the one place needing resolution (docs/PLAN.md section 4).
! amr_refine_front rebuilds the octree from the CURRENT solution:
!   - an old leaf is flagged as front when eps < x_HI < 1-eps
!     (par%refine_eps) or a face neighbor differs by more than
!     par%refine_dx in x_HI (gas cells only);
!   - the new tree is level par%refine_lmax wherever a cell overlaps a
!     front-flagged old leaf, par%refine_lbase elsewhere;
!   - the gas state (and rhokap) map over by position: each new leaf
!     samples the old leaf containing its center (amr_build_tree is cheap
!     to re-run; the state maps by position).
! All ranks build identical trees from the shared state; the shared-memory
! windows of the OLD tree/arrays are not recycled (one refinement event
! leaks the old grid until MPI_FINALIZE — acceptable; window recycling is
! a later cleanup).  Call between iterations, after the equilibrium
! update; the caller continues with the refilled opacity and a fresh tally.
!---------------------------------------------------------------------------
  use define
  use octree_mod
  use gas_state_mod
  implicit none
  private

  public :: amr_refine_front

  !--- old-tree snapshot used during construction
  real(kind=wp), allocatable :: old_cx(:), old_cy(:), old_cz(:), old_ch(:)
  integer,       allocatable :: old_children(:,:), old_ileaf(:)
  logical,       allocatable :: front(:)

  !--- new leaf list (grown dynamically)
  integer :: nnew = 0, nmax = 0
  real(kind=wp), allocatable :: nx_(:), ny_(:), nz_(:)
  integer,       allocatable :: nlev_(:)

contains

  !=========================================================================
  subroutine amr_refine_front()
    use mpi
    use gas_opacity_mod, only : gas_opacity_setup, kap_ion
    use jtally_mod,      only : jtally_ion_resize, jt_ion
    use species_mod,     only : species_resize
    use gas_state_mod,   only : gas_nH_p => gas_nH, gas_xHI_p => gas_xHI, &
                                gas_xHeI_p => gas_xHeI, &
                                gas_xHeII_p => gas_xHeII, &
                                gas_ne_p => gas_ne, gas_Te_p => gas_Te
    use memory_mod,      only : destroy_shared_mem
    implicit none
    real(kind=wp), allocatable :: s_nH(:), s_x1(:), s_x2(:), s_x3(:), &
                                  s_ne(:), s_te(:), s_rk(:)
    real(kind=wp) :: half, xr, dx
    integer :: il, jl, icell, jcell, iface, nold, nfront, ierr

    nold = amr_grid%nleaf
    half = 0.5_wp*amr_grid%L_box

    !--- 1. front flags on the old leaves
    allocate(front(nold))
    front = .false.
    do il = 1, nold
       if (gas_nH(il) <= 0.0_wp) cycle
       xr = gas_xHI(il)
       if (xr > par%refine_eps .and. xr < 1.0_wp - par%refine_eps) then
          front(il) = .true.
          cycle
       end if
       icell = amr_grid%icell_of_leaf(il)
       do iface = 1, 6
          jcell = amr_grid%neighbor(iface, icell)
          if (jcell <= 0) cycle
          jl = amr_grid%ileaf(jcell)
          if (jl <= 0) cycle
          if (gas_nH(jl) <= 0.0_wp) cycle
          dx = abs(gas_xHI(jl) - xr)
          if (dx > par%refine_dx) then
             front(il) = .true.
             exit
          end if
       end do
    end do
    nfront = count(front)

    !--- 2. snapshot the old tree topology (local copies; the shared arrays
    !---    stay alive but we avoid touching them mid-rebuild)
    allocate(old_cx(amr_grid%ncells), old_cy(amr_grid%ncells), &
             old_cz(amr_grid%ncells), old_ch(amr_grid%ncells), &
             old_children(8, amr_grid%ncells), old_ileaf(amr_grid%ncells))
    old_cx = amr_grid%cx(1:amr_grid%ncells)
    old_cy = amr_grid%cy(1:amr_grid%ncells)
    old_cz = amr_grid%cz(1:amr_grid%ncells)
    old_ch = amr_grid%ch(1:amr_grid%ncells)
    old_children = amr_grid%children(:, 1:amr_grid%ncells)
    old_ileaf    = amr_grid%ileaf(1:amr_grid%ncells)

    !--- 3. build the new leaf list
    nnew = 0;  nmax = max(2*nold, 1024)
    allocate(nx_(nmax), ny_(nmax), nz_(nmax), nlev_(nmax))
    call rec_build(0.0_wp, 0.0_wp, 0.0_wp, half, 0)

    !--- 4. sample the state at the new leaf centers from the old tree
    allocate(s_nH(nnew), s_x1(nnew), s_x2(nnew), s_x3(nnew), &
             s_ne(nnew), s_te(nnew), s_rk(nnew))
    do il = 1, nnew
       jl = amr_find_leaf(nx_(il), ny_(il), nz_(il))
       s_nH(il) = gas_nH(jl);   s_x1(il) = gas_xHI(jl)
       s_x2(il) = gas_xHeI(jl); s_x3(il) = gas_xHeII(jl)
       s_ne(il) = gas_ne(jl);   s_te(il) = gas_Te(jl)
       s_rk(il) = amr_grid%rhokap(jl)
    end do

    !--- 4b. recycle the OLD grid's shared-memory windows: nothing below
    !--- reads the old tree or state (everything needed is in the local
    !--- old_*/s_* copies), so free them BEFORE the rebuild — collective,
    !--- every rank runs this routine.  Without the recycling every
    !--- re-refinement leaks a full grid's shared memory until finalize.
    !--- jt_ion is a per-rank array, not a window — jtally_ion_resize
    !--- deallocates it properly.
    call destroy_shared_mem(amr_grid%parent)
    call destroy_shared_mem(amr_grid%children)
    call destroy_shared_mem(amr_grid%level)
    call destroy_shared_mem(amr_grid%cx)
    call destroy_shared_mem(amr_grid%cy)
    call destroy_shared_mem(amr_grid%cz)
    call destroy_shared_mem(amr_grid%ch)
    call destroy_shared_mem(amr_grid%ileaf)
    call destroy_shared_mem(amr_grid%icell_of_leaf)
    call destroy_shared_mem(amr_grid%neighbor)
    call destroy_shared_mem(amr_grid%rhokap)
    call destroy_shared_mem(gas_nH_p)
    call destroy_shared_mem(gas_xHI_p)
    call destroy_shared_mem(gas_xHeI_p)
    call destroy_shared_mem(gas_xHeII_p)
    call destroy_shared_mem(gas_ne_p)
    call destroy_shared_mem(gas_Te_p)
    call destroy_shared_mem(kap_ion)
    if (mpar%p_rank == 0) write(*,'(a)') &
       ' AMR: recycled the shared-memory windows of the old grid'

    !--- 5. rebuild the octree + neighbors + rhokap
    call amr_build_tree(nx_(1:nnew), ny_(1:nnew), nz_(1:nnew), &
                        nlev_(1:nnew), nnew, -half, half, -half, half, &
                        -half, half)
    call amr_build_neighbors
    call amr_alloc_phys()
    if (mpar%h_rank == 0) amr_grid%rhokap(1:nnew) = s_rk
    call MPI_BARRIER(mpar%hostcomm, ierr)

    !--- 6. recreate the gas state and dependent arrays
    call gas_state_recreate(nnew, s_nH, s_x1, s_x2, s_x3, s_ne, s_te)
    !--- species BEFORE opacity: with par%ion_metal_abs the opacity fill
    !--- reads the Gamma blocks at the NEW leaf count.
    if (par%use_metals) call species_resize(nnew)
    call gas_opacity_setup()
    call jtally_ion_resize(nnew)

    if (mpar%p_rank == 0) write(*,'(a,i0,a,i0,a,i0,a,i0,a)') &
       ' AMR: re-refined on the I-front: ', nold, ' -> ', nnew, &
       ' leaves (', nfront, ' front cells, lmax = ', par%refine_lmax, ')'

    deallocate(front, old_cx, old_cy, old_cz, old_ch, old_children, &
               old_ileaf, nx_, ny_, nz_, nlev_, s_nH, s_x1, s_x2, s_x3, &
               s_ne, s_te, s_rk)

  end subroutine amr_refine_front

  !=========================================================================
  recursive subroutine rec_build(x, y, z, h, level)
    implicit none
    real(kind=wp), intent(in) :: x, y, z, h
    integer,       intent(in) :: level
    real(kind=wp) :: hc
    integer :: ix, iy, iz

    if ( (level < par%refine_lmax .and. cell_has_front(x, y, z, h)) .or. &
         level < par%refine_lbase ) then
       hc = 0.5_wp*h
       do iz = 0, 1
          do iy = 0, 1
             do ix = 0, 1
                call rec_build(x + real(2*ix-1,wp)*hc, &
                               y + real(2*iy-1,wp)*hc, &
                               z + real(2*iz-1,wp)*hc, hc, level+1)
             end do
          end do
       end do
    else
       nnew = nnew + 1
       if (nnew > nmax) call grow_new()
       nx_(nnew) = x;  ny_(nnew) = y;  nz_(nnew) = z;  nlev_(nnew) = level
    end if
  end subroutine rec_build

  !=========================================================================
  ! Does the cell (center, half width h) overlap any front-flagged old
  ! leaf?  Old-tree descent with pruning (explicit stack).
  !=========================================================================
  logical function cell_has_front(x, y, z, h) result(hit)
    implicit none
    real(kind=wp), intent(in) :: x, y, z, h
    integer :: stack(512), nstk, icell, ioct, jc

    hit = .false.
    nstk = 1;  stack(1) = 1
    do while (nstk > 0)
       icell = stack(nstk);  nstk = nstk - 1
       if (abs(old_cx(icell) - x) > old_ch(icell) + h) cycle
       if (abs(old_cy(icell) - y) > old_ch(icell) + h) cycle
       if (abs(old_cz(icell) - z) > old_ch(icell) + h) cycle
       if (old_ileaf(icell) > 0) then
          if (front(old_ileaf(icell))) then
             hit = .true.
             return
          end if
          cycle
       end if
       do ioct = 1, 8
          jc = old_children(ioct, icell)
          if (jc > 0) then
             nstk = nstk + 1
             stack(nstk) = jc
          end if
       end do
    end do
  end function cell_has_front

  !=========================================================================
  subroutine grow_new()
    implicit none
    real(kind=wp), allocatable :: tx(:), ty(:), tz(:)
    integer,       allocatable :: tl(:)
    call move_alloc(nx_, tx);  call move_alloc(ny_, ty)
    call move_alloc(nz_, tz);  call move_alloc(nlev_, tl)
    nmax = 2*nmax
    allocate(nx_(nmax), ny_(nmax), nz_(nmax), nlev_(nmax))
    nx_(1:size(tx)) = tx;  ny_(1:size(ty)) = ty
    nz_(1:size(tz)) = tz;  nlev_(1:size(tl)) = tl
    deallocate(tx, ty, tz, tl)
  end subroutine grow_new

end module amr_refine_mod
