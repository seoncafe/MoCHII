module amr_refine_mod
!---------------------------------------------------------------------------
! MoCHII: solution-driven octree re-refinement.
!
! The I-front is the one place needing resolution.
! amr_refine_front rebuilds the octree from the CURRENT solution:
!   - an old leaf is flagged as front when eps < x_HI < 1-eps
!     (par%refine_eps) or a face neighbor differs by more than
!     par%refine_dx in x_HI (gas cells only);
!   - the new tree is level par%refine_lmax wherever a cell overlaps a
!     front-flagged old leaf, par%refine_lbase elsewhere;
!   - the gas state (and rhokap) are remapped by OVERLAP VOLUME, so the
!     box integrals of H nuclei, each ion stage, electrons, thermal energy
!     and the dust extinction survive the rebuild.  Where a leaf is only
!     subdivided the sum collapses to the parent's values; the weighting
!     matters when heterogeneous leaves are coarsened back together, which
!     is what refine_lbase does away from the front.  Every event prints
!     the before/after integrals and aborts if one of them moved.
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
    use gas_opacity_mod, only : gas_opacity_setup
    use jtally_mod,      only : jtally_ion_resize
    use ion_score_mod,   only : ion_score_resize
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

    !--- 4. remap the state onto the new leaves by volume overlap.
    !--- Sampling the one old leaf containing the new leaf center is exact
    !--- while a leaf is only subdivided, but it is not conservative when
    !--- heterogeneous old leaves are coarsened into a larger new one: the
    !--- center cell's state would stand in for the whole region, and H
    !--- nuclei, the ion stage counts, electrons, thermal energy and the
    !--- dust extinction integral would all jump.  That matters exactly at
    !--- the discontinuities re-refinement exists to resolve.  remap_leaf
    !--- collapses to the old behavior under pure refinement.
    allocate(s_nH(nnew), s_x1(nnew), s_x2(nnew), s_x3(nnew), &
             s_ne(nnew), s_te(nnew), s_rk(nnew))
    do il = 1, nnew
       call remap_leaf(nx_(il), ny_(il), nz_(il), half/2.0_wp**nlev_(il), &
                       s_nH(il), s_x1(il), s_x2(il), s_x3(il), &
                       s_ne(il), s_te(il), s_rk(il))
    end do

    !--- 4a. conservation report.  These are integrals over the whole box
    !--- and are preserved exactly by the overlap weights, so anything above
    !--- roundoff is a coding error rather than a modeling choice.
    call conservation_report(nold, half, s_nH, s_x1, s_x2, s_x3, &
                             s_ne, s_te, s_rk)

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
    !--- species BEFORE opacity: the stage-fraction refresh inside
    !--- gas_opacity_fill reads the Gamma blocks at the NEW leaf count.
    if (par%use_metals) call species_resize(nnew)
    call gas_opacity_setup()
    call jtally_ion_resize(nnew)
    call ion_score_resize(nnew)

    if (mpar%p_rank == 0) write(*,'(a,i0,a,i0,a,i0,a,i0,a)') &
       ' AMR: re-refined on the I-front: ', nold, ' -> ', nnew, &
       ' leaves (', nfront, ' front cells, lmax = ', par%refine_lmax, ')'

    deallocate(front, old_cx, old_cy, old_cz, old_ch, old_children, &
               old_ileaf, nx_, ny_, nz_, nlev_, s_nH, s_x1, s_x2, s_x3, &
               s_ne, s_te, s_rk)

  end subroutine amr_refine_front

  !=========================================================================
  ! Overlap of two intervals given as (center, half width).
  !=========================================================================
  pure real(kind=wp) function overlap_1d(ca, ha, cb, hb) result(ov)
    real(kind=wp), intent(in) :: ca, ha, cb, hb
    ov = min(ca + ha, cb + hb) - max(ca - ha, cb - hb)
    if (ov < 0.0_wp) ov = 0.0_wp
  end function overlap_1d

  !=========================================================================
  ! State of one new leaf (center x,y,z, half width h) from the old tree,
  ! weighted by the overlap volume.  Both trees are octrees on the same root
  ! box, so a new leaf and an old leaf are either disjoint or one contains
  ! the other; the overlap is the product of the three 1-D overlaps in both
  ! cases.  When the new leaf lies inside a single old leaf the sum collapses
  ! to that leaf's values, so pure refinement reproduces center sampling.
  !
  ! The weight of each field is set by the quantity that must be conserved:
  !   H nuclei          sum(nH V)        -> nH by volume
  !   ion stage counts  sum(nH x_i V)    -> x_i by nH V
  !   electrons         sum(ne V)        -> ne by volume
  !   dust extinction   sum(rhokap V)    -> rhokap by volume
  !   thermal energy    sum(n_tot Te V)  -> Te by n_tot V, with
  !                                         n_tot = nH(1 + He_abund) + ne
  ! Where a region holds no material the fractions and Te carry no thermal
  ! content, and the plain volume average is used so the value stays inside
  ! the range of the old leaves instead of being invented.
  !=========================================================================
  subroutine remap_leaf(x, y, z, h, nH, x1, x2, x3, ne, te, rk)
    use mpi
    implicit none
    real(kind=wp), intent(in)  :: x, y, z, h
    real(kind=wp), intent(out) :: nH, x1, x2, x3, ne, te, rk
    integer, parameter :: MAXSTK = 512
    integer :: stack(MAXSTK), nstk, icell, ioct, jc, jl, ierr
    real(kind=wp) :: ov, hw, wm, wt, vol
    real(kind=wp) :: sv, snh, sne, srk                  ! volume weighted
    real(kind=wp) :: sm, smx1, smx2, smx3               ! nH V weighted
    real(kind=wp) :: st, ste                            ! n_tot V weighted
    real(kind=wp) :: svx1, svx2, svx3, svte             ! volume fallback

    sv = 0.0_wp;  snh = 0.0_wp;  sne = 0.0_wp;  srk = 0.0_wp
    sm = 0.0_wp;  smx1 = 0.0_wp; smx2 = 0.0_wp; smx3 = 0.0_wp
    st = 0.0_wp;  ste = 0.0_wp
    svx1 = 0.0_wp; svx2 = 0.0_wp; svx3 = 0.0_wp; svte = 0.0_wp

    nstk = 1;  stack(1) = 1
    do while (nstk > 0)
       icell = stack(nstk);  nstk = nstk - 1
       hw = old_ch(icell)
       if (abs(old_cx(icell) - x) >= hw + h) cycle
       if (abs(old_cy(icell) - y) >= hw + h) cycle
       if (abs(old_cz(icell) - z) >= hw + h) cycle
       if (old_ileaf(icell) > 0) then
          jl = old_ileaf(icell)
          ov = overlap_1d(old_cx(icell), hw, x, h) &
             * overlap_1d(old_cy(icell), hw, y, h) &
             * overlap_1d(old_cz(icell), hw, z, h)
          if (ov <= 0.0_wp) cycle
          sv   = sv  + ov
          snh  = snh + gas_nH(jl)*ov
          sne  = sne + gas_ne(jl)*ov
          srk  = srk + amr_grid%rhokap(jl)*ov
          wm   = gas_nH(jl)*ov
          sm   = sm   + wm
          smx1 = smx1 + wm*gas_xHI(jl)
          smx2 = smx2 + wm*gas_xHeI(jl)
          smx3 = smx3 + wm*gas_xHeII(jl)
          wt   = (gas_nH(jl)*(1.0_wp + par%He_abund) + gas_ne(jl))*ov
          st   = st   + wt
          ste  = ste  + wt*gas_Te(jl)
          svx1 = svx1 + ov*gas_xHI(jl)
          svx2 = svx2 + ov*gas_xHeI(jl)
          svx3 = svx3 + ov*gas_xHeII(jl)
          svte = svte + ov*gas_Te(jl)
          cycle
       end if
       do ioct = 1, 8
          jc = old_children(ioct, icell)
          if (jc > 0) then
             nstk = nstk + 1
             if (nstk > MAXSTK) then
                if (mpar%p_rank == 0) write(*,'(a)') &
                   'ERROR: amr remap descent stack overflow'
                call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
             end if
             stack(nstk) = jc
          end if
       end do
    end do

    vol = (2.0_wp*h)**3
    nH = snh/vol
    ne = sne/vol
    rk = srk/vol
    if (sm > 0.0_wp) then
       x1 = smx1/sm;  x2 = smx2/sm;  x3 = smx3/sm
    else if (sv > 0.0_wp) then
       x1 = svx1/sv;  x2 = svx2/sv;  x3 = svx3/sv
    else
       x1 = 1.0_wp;  x2 = 1.0_wp;  x3 = 0.0_wp
    end if
    if (st > 0.0_wp) then
       te = ste/st
    else if (sv > 0.0_wp) then
       te = svte/sv
    else
       te = par%te_min
    end if
  end subroutine remap_leaf

  !=========================================================================
  ! Box integrals before and after the rebuild.  The overlap weights make
  ! each of these exact, so a deviation above roundoff means the remap is
  ! wrong, not that the model changed.
  !=========================================================================
  subroutine conservation_report(nold, half, s_nH, s_x1, s_x2, s_x3, &
                                 s_ne, s_te, s_rk)
    use mpi
    use octree_mod, only : leaf_half
    implicit none
    integer,       intent(in) :: nold
    real(kind=wp), intent(in) :: half
    real(kind=wp), intent(in) :: s_nH(:), s_x1(:), s_x2(:), s_x3(:), &
                                 s_ne(:), s_te(:), s_rk(:)
    integer, parameter :: NQ = 7
    character(len=12), parameter :: qname(NQ) = [character(len=12) :: &
       'H nuclei', 'H I', 'He I', 'He II', 'electrons', 'extinction', &
       'thermal']
    real(kind=wp) :: a(NQ), b(NQ), vol, dev, worst
    integer :: il, ierr

    a = 0.0_wp
    do il = 1, nold
       vol  = (2.0_wp*leaf_half(il))**3
       a(1) = a(1) + gas_nH(il)*vol
       a(2) = a(2) + gas_nH(il)*gas_xHI(il)*vol
       a(3) = a(3) + gas_nH(il)*par%He_abund*gas_xHeI(il)*vol
       a(4) = a(4) + gas_nH(il)*par%He_abund*gas_xHeII(il)*vol
       a(5) = a(5) + gas_ne(il)*vol
       a(6) = a(6) + amr_grid%rhokap(il)*vol
       a(7) = a(7) + (gas_nH(il)*(1.0_wp + par%He_abund) + gas_ne(il)) &
                     *gas_Te(il)*vol
    end do

    b = 0.0_wp
    do il = 1, nnew
       vol  = (2.0_wp*half/2.0_wp**nlev_(il))**3
       b(1) = b(1) + s_nH(il)*vol
       b(2) = b(2) + s_nH(il)*s_x1(il)*vol
       b(3) = b(3) + s_nH(il)*par%He_abund*s_x2(il)*vol
       b(4) = b(4) + s_nH(il)*par%He_abund*s_x3(il)*vol
       b(5) = b(5) + s_ne(il)*vol
       b(6) = b(6) + s_rk(il)*vol
       b(7) = b(7) + (s_nH(il)*(1.0_wp + par%He_abund) + s_ne(il)) &
                     *s_te(il)*vol
    end do

    worst = 0.0_wp
    do il = 1, NQ
       dev = 0.0_wp
       if (abs(a(il)) > 0.0_wp) dev = abs(b(il) - a(il))/abs(a(il))
       worst = max(worst, dev)
       if (mpar%p_rank == 0) write(*,'(3a,es14.6,a,es14.6,a,es9.2)') &
          ' AMR: conserved ', qname(il), ': ', a(il), ' -> ', b(il), &
          '   rel ', dev
    end do
    if (worst > 1.0e-8_wp) then
       if (mpar%p_rank == 0) write(*,'(a,es10.3,a)') &
          'ERROR: AMR re-refinement lost ', worst, &
          ' of a conserved box integral; the remap is wrong.'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
  end subroutine conservation_report

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
