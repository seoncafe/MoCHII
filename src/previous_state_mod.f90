module previous_state_mod
!---------------------------------------------------------------------------
! MoCHII: start the gas iteration from the state a previous run wrote
! (par%restart_file = that run's <base>_rates file, HDF5 or FITS).
!
! What a transport pass depends on, beyond the grid and the sources, is
!   x_HI, x_HII, x_HeI, x_HeII
!                         the H and He opacity, and the start of the next
!                         solve; x_HII is read where the file carries it
!                         (written since 2026-09-26), else taken as
!                         1 - x_HI with a printed note -- beyond the
!                         ionization front that has lost x_HII's digits, so
!                         the first transport pass (metal stages, n_HII)
!                         sees that state until the first cell solve
!                         recovers x_HII; likewise x_HeIII (written since
!                         2026-09-26), else max(0, 1 - x_HeI - x_HeII);
!   T_e, n_e              the diffuse-field emission (recombination spectra,
!                         alpha(T_e), n_e n_ion) of the first diffuse pass;
!   Gamma_<el>_stages     the metal photoionization rates: the stage-fraction
!                         cache behind the metal opacity is built from them
!                         (species_stage_cache_refresh), and a fresh run
!                         starts them at zero, i.e. with metals in their low
!                         stages -- where x_HI ~ 1e-3 that metal absorption is
!                         as large as hydrogen's.
!   Phi_HeI_23S           the He I 2^3S photoionization rate of the last
!                         pass (with par%hei_metastable_coupling): the 2^3S
!                         population that splits the He I diffuse channel of
!                         the first pass is solved with it.  Without it (a
!                         file written before 2026-10-12, or a run without
!                         the coupling) the first pass takes Phi_3 = 0, as a
!                         fresh run does, and a restart is not the
!                         continuation of the run that wrote the file.
!                         It is held here and set after ion_score_setup
!                         (restore_hei_phi3), which clears the leaf rates.
! Everything else is rebuilt every pass from these: the H/He rates and heat,
! the diffuse emission, the metal stage fractions, the PAH survival factor.
! All of them are written to the rates file already.
!
! Leaves are matched by position: each new leaf takes the values of the
! previous leaf that contains its center, so the new grid may differ from the
! old one (refined, coarsened, extended).  A new leaf with no previous
! counterpart -- outside the old domain, or gas where the old leaf was
! vacuum (n_e = 0) -- keeps the par%*_init values.  The lookup is a hash on
! the integer position of the containing cell, one table per distinct leaf
! size, so a deep octree costs no dense lattice.
!
! The convergence history (fixed_point_distance_mod) starts empty: the
! stopping test needs NWIN = 8 iterations of the new run before it can pass.
!---------------------------------------------------------------------------
  use define
  use mpi
  implicit none
  private

  public :: load_previous_state, restore_hei_phi3
  !--- Phi_3 on the new leaves, held until restore_hei_phi3
  real(kind=wp), allocatable :: restart_phi3(:)

  integer, parameter :: MAX_SIZES = 32
  integer(kind=int64), parameter :: EMPTY = -1_int64

contains

  subroutine load_previous_state(fname)
    use iofile_mod
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHII, gas_xHeI, gas_xHeII, &
                              gas_xHeIII, gas_ne, gas_Te, gas_nleaf, &
                              gas_state_sync
    use species_mod,   only : n_elements, elem_name, elem_nstage, &
                              species_set_photoionization_rates
    use utility,       only : fatal_error
    character(len=*), intent(in) :: fname

    type(io_file_type) :: file
    character(len=64)  :: extname
    real(kind=wp), allocatable :: oxyz(:,:), osize(:), ox1(:), ox1p(:), ox2(:), &
                                  ox3(:), ox4(:), one(:), ote(:), g(:,:), &
                                  ophi3(:)
    real(kind=wp), allocatable :: nx1(:), nx1p(:), nx2(:), nx3(:), nx4(:), &
                                  nne(:), nte(:)
    type :: rate_block
       real(kind=wp), allocatable :: g(:,:)
       character(len=8) :: name = ''
    end type rate_block
    type(rate_block), allocatable :: orates(:)
    integer, allocatable :: match(:)
    integer :: status, n1, n2, nold, il, j, ie, k, nr, nmatched, ngas_new, &
               nunmatched_gas, ierr, nst
    logical :: have(8), have_xhii, have_xheiii, have_phi3
    logical :: restored_rates(n_elements)

    nold = 0;  nr = 0;  have = .false.;  have_xhii = .false.
    have_xheiii = .false.;  have_phi3 = .false.
    ! The old run can contain more elements than this run. Read its full
    ! supported roster so earlier, now-disabled elements do not hide the
    ! rate block of an active element later in the file.
    allocate(orates(MAX_METAL_ELEMENTS))
    restored_rates = .false.

    !--- world rank 0 reads the previous state
    if (mpar%p_rank == 0) then
       status = 0
       call io_open_old(file, trim(fname), status)
       if (status /= 0) call fatal_error('restart_file: cannot open '//trim(fname))
       do
          extname = ''
          call io_get_keyword(file, 'EXTNAME', extname, status)
          if (status /= 0) exit
          n1 = 0;  n2 = 1
          call io_get_keyword(file, 'NAXIS1', n1, status)
          if (index(extname, '_stages') > 0 .or. trim(extname) == 'LeafXYZ') &
             call io_get_keyword(file, 'NAXIS2', n2, status)
          select case (trim(extname))
          case ('LeafXYZ')
             allocate(oxyz(n1, n2));  call io_read_image(file, oxyz, status)
             nold = n1;  have(1) = .true.
          case ('LeafSize')
             allocate(osize(n1));  call io_read_image(file, osize, status);  have(2) = .true.
          case ('x_HI')
             allocate(ox1(n1));  call io_read_image(file, ox1, status);  have(3) = .true.
          case ('x_HII')
             allocate(ox1p(n1));  call io_read_image(file, ox1p, status)
             have_xhii = .true.
          case ('x_HeI')
             allocate(ox2(n1));  call io_read_image(file, ox2, status);  have(4) = .true.
          case ('x_HeII')
             allocate(ox3(n1));  call io_read_image(file, ox3, status);  have(5) = .true.
          case ('x_HeIII')
             allocate(ox4(n1));  call io_read_image(file, ox4, status)
             have_xheiii = .true.
          case ('n_e')
             allocate(one(n1));  call io_read_image(file, one, status);  have(6) = .true.
          case ('T_e')
             allocate(ote(n1));  call io_read_image(file, ote, status);  have(7) = .true.
          case ('Phi_HeI_23S')
             allocate(ophi3(n1));  call io_read_image(file, ophi3, status)
             have_phi3 = .true.
          case default
             if (extname(1:6) == 'Gamma_' .and. index(extname, '_stages') > 0 &
                 .and. nr < size(orates)) then
                nr = nr + 1
                orates(nr)%name = extname(7:index(extname, '_stages')-1)
                allocate(orates(nr)%g(n1, n2))
                call io_read_image(file, orates(nr)%g, status)
             end if
          end select
          if (status /= 0) call fatal_error('restart_file: cannot read '// &
                                            trim(extname)//' in '//trim(fname))
          call io_move_to_next_section(file, status)
          if (status /= 0) exit
       end do
       status = 0
       call io_close(file, status)
       if (.not. all(have(1:7))) call fatal_error('restart_file: '//trim(fname)// &
          ' lacks one of LeafXYZ, LeafSize, x_HI, x_HeI, x_HeII, n_e, T_e')
       if (.not. have_xhii) then
          !--- a file written before x_HII was stored (2026-09-26)
          allocate(ox1p(size(ox1)));  ox1p = 1.0_wp - ox1
          write(*,'(a)') ' RESTART: NOTE: '//trim(fname)//' carries no x_HII'// &
             ' (written before 2026-09-26); x_HII = 1 - x_HI, which beyond'// &
             ' the ionization front has lost its digits: the first transport'// &
             ' pass sees that state, the first cell solve recovers x_HII'
       end if
       if (.not. have_xheiii) then
          !--- a file written before x_HeIII was stored (2026-09-26)
          allocate(ox4(size(ox2)));  ox4 = max(0.0_wp, 1.0_wp - (ox2 + ox3))
          write(*,'(a)') ' RESTART: NOTE: '//trim(fname)//' carries no x_HeIII'// &
             ' (written before 2026-09-26); x_HeIII = max(0, 1 - x_HeI - x_HeII),'// &
             ' which where He2+ is scarce is rounding: the first transport'// &
             ' pass sees that state, the first cell solve recovers x_HeIII'
       end if
    end if

    !--- rank 0 matches every new leaf to the previous leaf containing its
    !--- center and assembles the new-grid arrays; they are then broadcast.
    allocate(nx1(gas_nleaf), nx1p(gas_nleaf), nx2(gas_nleaf), nx3(gas_nleaf), &
             nx4(gas_nleaf), nne(gas_nleaf), nte(gas_nleaf), match(gas_nleaf))
    if (mpar%p_rank == 0) then
       call match_leaves(oxyz, osize, match)
       nmatched = 0;  ngas_new = 0;  nunmatched_gas = 0
       do il = 1, gas_nleaf
          j = match(il)
          !--- a previous vacuum leaf carries no state for new gas
          if (j > 0) then
             if (one(j) <= 0.0_wp .and. gas_nH(il) > 0.0_wp) j = 0
          end if
          match(il) = j
          if (gas_nH(il) > 0.0_wp) ngas_new = ngas_new + 1
          if (j > 0) then
             nmatched = nmatched + 1
             nx1(il) = ox1(j);  nx1p(il) = ox1p(j)
             nx2(il) = ox2(j);  nx3(il) = ox3(j);  nx4(il) = ox4(j)
             nne(il) = one(j);  nte(il) = ote(j)
          else
             if (gas_nH(il) > 0.0_wp) nunmatched_gas = nunmatched_gas + 1
             nx1(il) = gas_xHI(il);  nx1p(il) = gas_xHII(il)
             nx2(il) = gas_xHeI(il);  nx3(il) = gas_xHeII(il)
             nx4(il) = gas_xHeIII(il)
             nne(il) = gas_ne(il);   nte(il) = gas_Te(il)
          end if
       end do
       write(*,'(a,a)')  ' RESTART: gas state from ', trim(fname)
       write(*,'(a,i10,a,i10,a,i10,a)') ' RESTART: ', nold, ' previous leaves; ', &
          nmatched, ' new leaves matched; ', nunmatched_gas, &
          ' gas leaves unmatched (kept at the *_init values)'
    end if
    call MPI_BCAST(match, gas_nleaf, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nx1, gas_nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nx1p, gas_nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nx2, gas_nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nx3, gas_nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nx4, gas_nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nne, gas_nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nte, gas_nleaf, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (mpar%h_rank == 0) then
       gas_xHI(1:gas_nleaf)   = nx1;  gas_xHII(1:gas_nleaf) = nx1p
       gas_xHeI(1:gas_nleaf)  = nx2
       gas_xHeII(1:gas_nleaf) = nx3;  gas_xHeIII(1:gas_nleaf) = nx4
       gas_ne(1:gas_nleaf)    = nne
       gas_Te(1:gas_nleaf)    = nte
    end if
    call gas_state_sync()

    !--- the He I 2^3S photoionization rate (0 on unmatched leaves)
    call MPI_BCAST(have_phi3, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (par%hei_metastable_coupling) then
       allocate(restart_phi3(gas_nleaf));  restart_phi3 = 0.0_wp
       if (have_phi3) then
          if (mpar%p_rank == 0) then
             do il = 1, gas_nleaf
                if (match(il) > 0) restart_phi3(il) = ophi3(match(il))
             end do
          end if
          call MPI_BCAST(restart_phi3, gas_nleaf, MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_WORLD, ierr)
          if (mpar%p_rank == 0) write(*,'(a)') &
             ' RESTART: He I 2^3S photoionization rate Phi_3'
       else if (mpar%p_rank == 0) then
          write(*,'(a)') ' RESTART: NOTE: '//trim(fname)//' carries no'// &
             ' Phi_HeI_23S; the first pass takes Phi_3 = 0, as a fresh run'// &
             ' does, so the 2^3S split of the He I diffuse channel differs'// &
             ' from the run that wrote the file for that pass'
       end if
    end if

    !--- metal photoionization rates, element by element (by name)
    call MPI_BCAST(nr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    do k = 1, nr
       call MPI_BCAST(orates(k)%name, 8, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
       n1 = 0;  n2 = 0
       if (mpar%p_rank == 0) then
          n1 = size(orates(k)%g, 1);  n2 = size(orates(k)%g, 2)
       end if
       call MPI_BCAST(n1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       ie = 0
       do j = 1, n_elements
          if (trim(elem_name(j)) == trim(orates(k)%name)) ie = j
       end do
       if (ie == 0) cycle                     ! element not in this run
       nst = min(n1, elem_nstage(ie) - 1)
       restored_rates(ie) = n1 >= elem_nstage(ie)-1
       allocate(g(elem_nstage(ie) - 1, gas_nleaf));  g = 0.0_wp
       if (mpar%p_rank == 0) then
          do il = 1, gas_nleaf
             if (match(il) > 0) g(1:nst, il) = orates(k)%g(1:nst, match(il))
          end do
       end if
       call MPI_BCAST(g, size(g), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       call species_set_photoionization_rates(ie, g)
       if (mpar%p_rank == 0) write(*,'(a,a,a,i0,a)') ' RESTART: ', &
          trim(orates(k)%name), ' photoionization rates, ', nst, ' transitions'
       deallocate(g)
    end do
    if (any(.not. restored_rates) .and. mpar%p_rank == 0) then
       write(*,'(a)') ' RESTART: incomplete metal transition rates; added rates start at zero.'
       write(*,'(a)') ' RESTART: fresh transport and gas reconvergence are required before output.'
       if (par%gas_niter < 1) &
          call fatal_error('restart with missing metal rates requires gas_niter >= 1')
       ! The ordinary main-loop convergence check now refuses an output
       ! whose warm start has not reconverged after the stage change.
       par%require_convergence = .true.
    end if
  end subroutine load_previous_state

  !--- set the restored Phi_3 on the leaves once ion_score_setup has
  !--- allocated (and cleared) the transport estimators; its heating is
  !--- read only after the next transport, which replaces both.
  subroutine restore_hei_phi3()
    use hei_metastable_mod, only : hei_phi3_assign
    real(kind=wp), allocatable :: zero(:)
    if (.not. allocated(restart_phi3)) return
    allocate(zero(size(restart_phi3)));  zero = 0.0_wp
    call hei_phi3_assign(restart_phi3, zero)
    deallocate(restart_phi3, zero)
  end subroutine restore_hei_phi3

  !=========================================================================
  ! match(il) = the previous leaf whose cube contains the center of new leaf
  ! il, or 0.  Previous leaves are grouped by size; within a size class a
  ! leaf is keyed by the integer position floor((x - xmin)/s) of its cube,
  ! and a new center is looked up class by class.
  !=========================================================================
  subroutine match_leaves(oxyz, osize, match)
    use octree_mod,    only : leaf_cx, leaf_cy, leaf_cz
    use gas_state_mod, only : gas_nleaf
    use utility,       only : fatal_error
    real(kind=wp), intent(in)  :: oxyz(:,:), osize(:)
    integer,       intent(out) :: match(:)
    real(kind=wp) :: sizes(MAX_SIZES), org(3), p(3)
    integer :: nsize, nold, j, c, il, ijk(3), slot, nslot
    integer(kind=int64), allocatable :: keys(:)
    integer, allocatable :: vals(:)
    integer(kind=int64) :: key

    nold = size(osize)
    nsize = 0
    do j = 1, nold
       if (any(abs(sizes(1:nsize) - osize(j)) <= 1.0e-6_wp*osize(j))) cycle
       if (nsize == MAX_SIZES) call fatal_error('restart_file: more than 32 distinct leaf sizes')
       nsize = nsize + 1;  sizes(nsize) = osize(j)
    end do
    org = [minval(oxyz(:,1) - 0.5_wp*osize), minval(oxyz(:,2) - 0.5_wp*osize), &
           minval(oxyz(:,3) - 0.5_wp*osize)]

    !--- open-addressing hash table, load factor <= 1/2
    nslot = 1
    do while (nslot < 2*nold)
       nslot = 2*nslot
    end do
    allocate(keys(0:nslot-1), vals(0:nslot-1));  keys = EMPTY;  vals = 0
    do j = 1, nold
       c = size_class(osize(j))
       ijk = floor((oxyz(j,:) - org)/osize(j))
       key = pack_key(c, ijk)
       slot = int(modulo(hash(key), int(nslot, int64)))
       do while (keys(slot) /= EMPTY)
          slot = modulo(slot + 1, nslot)
       end do
       keys(slot) = key;  vals(slot) = j
    end do

    match = 0
    do il = 1, gas_nleaf
       p = [leaf_cx(il), leaf_cy(il), leaf_cz(il)]
       do c = 1, nsize
          if (any(p < org)) exit
          ijk = floor((p - org)/sizes(c))
          if (any(ijk >= 2**19)) cycle
          key = pack_key(c, ijk)
          slot = int(modulo(hash(key), int(nslot, int64)))
          do while (keys(slot) /= EMPTY)
             if (keys(slot) == key) then
                match(il) = vals(slot);  exit
             end if
             slot = modulo(slot + 1, nslot)
          end do
          if (match(il) > 0) exit
       end do
    end do
    deallocate(keys, vals)

  contains

    integer function size_class(s)
      real(kind=wp), intent(in) :: s
      integer :: m
      size_class = 0
      do m = 1, nsize
         if (abs(sizes(m) - s) <= 1.0e-6_wp*s) size_class = m
      end do
    end function size_class

  end subroutine match_leaves

  !--- class (5 bits) and three 19-bit cell indices in one 62-bit key
  pure integer(kind=int64) function pack_key(c, ijk) result(key)
    integer, intent(in) :: c, ijk(3)
    key = ishft(int(c, int64), 57) + ishft(int(ijk(1), int64), 38) &
        + ishft(int(ijk(2), int64), 19) + int(ijk(3), int64)
  end function pack_key

  pure integer(kind=int64) function hash(key) result(h)
    integer(kind=int64), intent(in) :: key
    !--- multiplicative mixing (Knuth); the modulo by the table size follows
    h = ieor(key, ishft(key, -29))*2654435761_int64
    h = ieor(h, ishft(h, -32))
    if (h < 0_int64) h = -(h + 1_int64)
  end function hash

end module previous_state_mod
