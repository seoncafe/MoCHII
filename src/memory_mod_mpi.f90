! Origin: MoCafe_v2.00/src/memory_mod_mpi.f90, copied on 2026-07-11.
! MoCHII retains its host-communicator decomposition. Shared allocation,
! ownership, publication, release, and named sums were rewritten for MoCHII.
module memory_mod
  use utility, only : abort_run
  use iso_c_binding, only : c_ptr, c_null_ptr, c_f_pointer, c_loc, c_associated, c_int8_t
  use iso_fortran_env, only : real64, int32, int64
  use mpi
  use mpi_f08, only : MPI_F_SYNC_REG
  implicit none
  private

  ! MPI communicators and windows retain MPI_ERRORS_ARE_FATAL.
  ! Shared storage has one explicit owner and one canonical pointer.
  ! Metadata copies are publication views, never independent allocations.
  ! Slot generations reject stale views after release or slot reuse.
  type, public :: shared_owner
    private
    integer :: window = MPI_WIN_NULL, slot = 0, rank = 0, value_kind = 0
    integer(int64) :: generation = 0_int64, shape(3) = 0_int64
    integer(MPI_ADDRESS_KIND) :: bytes = 0_MPI_ADDRESS_KIND
    type(c_ptr) :: base = c_null_ptr
  end type shared_owner
  integer(int64), allocatable :: live_generation(:)
  integer(int64) :: next_generation = 0_int64
  logical :: mpi_parameters_defined = .false.
  integer :: nprocs, p_rank, h_rank, hostcomm, same_hrank_comm, num_hosts
  integer, parameter :: SUM_ALL_RANKS = 1, SUM_SHARED_RESULT = 2, SUM_ROOT = 3
  integer(int64), parameter :: SUM_CHUNK = 100000000_int64
  public :: allocate_shared, free_shared, publish_shared, shared_is_live, shared_shape, shared_extent
  public :: shared_memory_finish
  public :: shared_host_rank, shared_host_size, shared_host_comm
  public :: sum_to_all_ranks, sum_to_root, sum_into_shared, sum_shared_hosts

  ! Allocation and release are collective over the host, in identical order.
  ! Never overwrite a live owner or an associated canonical pointer.
  ! Allocation publishes leader-zeroed storage before returning.
  ! Before reusing storage, finish ALL previous readers: call publish_shared
  ! before the first store unless an existing collective already orders them.
  ! After disjoint-writer stores, publish the entire physical update together.
  ! Publication orders readers; it does not make concurrent additions atomic.
  ! Release only through the canonical owner/pointer pair after all readers
  ! finish. It invalidates the token, closes the epoch and nullifies that
  ! pointer. Other pointer aliases cannot be nullified automatically; their
  ! module must drop them before release. There is no global implicit free.
  interface allocate_shared
    module procedure allocate_shared_1d_real64
    module procedure allocate_shared_2d_real64
    module procedure allocate_shared_3d_real64
    module procedure allocate_shared_1d_int32
    module procedure allocate_shared_2d_int32
  end interface
  interface free_shared
    module procedure free_shared_1d_real64
    module procedure free_shared_2d_real64
    module procedure free_shared_3d_real64
    module procedure free_shared_1d_int32
    module procedure free_shared_2d_int32
  end interface
  interface sum_to_all_ranks
    module procedure sum_to_all_ranks_1D_real64
    module procedure sum_to_all_ranks_2D_real64
  end interface
  interface sum_to_root
    module procedure sum_to_root_1D_real64
    module procedure sum_to_root_2D_real64
    module procedure sum_to_root_3D_real64
    module procedure sum_to_root_4D_real64
  end interface
  interface sum_into_shared
    module procedure sum_into_shared_1d_real64, sum_into_shared_2d_real64, sum_into_shared_3d_real64
  end interface
  ! sum_into_shared reads private scores without changing them. Its result
  ! is the whole allocation identified by an owner, never an array section.
  ! It finishes old readers before overwriting and publishes before return.
  ! Named sums accept private real64 storage, never node-shared arrays.
  ! CONTIGUOUS dummies allow packing/copy back of private strided sections.
  ! Only the named destinations hold totals; reset other contributions.
contains
  subroutine init_mpi_parameters()
    integer :: ierr
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, p_rank, ierr)
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, hostcomm, ierr)
    call MPI_COMM_RANK(hostcomm, h_rank, ierr)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, h_rank, p_rank, same_hrank_comm, ierr)
    call MPI_COMM_SIZE(same_hrank_comm, num_hosts, ierr)
    mpi_parameters_defined = .true.
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine init_mpi_parameters
  !---------------------------
  !--- MoCHII (2026-10-12): the rank, size and communicator of the node that
  !--- shares the windows allocate_shared creates, so a caller filling a
  !--- shared array divides its elements over exactly the ranks that share it
  integer function shared_host_rank() result(r)
    if (.not. mpi_parameters_defined) call init_mpi_parameters()
    r = h_rank
  end function shared_host_rank
  integer function shared_host_size() result(n)
    integer :: ierr
    if (.not. mpi_parameters_defined) call init_mpi_parameters()
    call MPI_COMM_SIZE(hostcomm, n, ierr)
  end function shared_host_size
  integer function shared_host_comm() result(c)
    if (.not. mpi_parameters_defined) call init_mpi_parameters()
    c = hostcomm
  end function shared_host_comm

  logical function shared_is_live(owner) result(live)
    type(shared_owner), intent(in) :: owner
    live = .false.
    if (owner%slot == 0) return
    if (.not. allocated(live_generation)) return
    if (owner%slot < 1 .or. owner%slot > size(live_generation)) return
    live = owner%generation /= 0_int64 .and. live_generation(owner%slot) == owner%generation
  end function shared_is_live

  subroutine require_shared_owner(owner)
    type(shared_owner), intent(in) :: owner
    if (.not. shared_is_live(owner)) call abort_run('inactive or stale shared owner')
  end subroutine require_shared_owner

  function shared_shape(owner) result(shape)
    type(shared_owner), intent(in) :: owner
    integer(int64), allocatable :: shape(:)
    call require_shared_owner(owner)
    shape = owner%shape(1:owner%rank)
  end function shared_shape

  function shared_extent(owner) result(bytes)
    type(shared_owner), intent(in) :: owner
    integer(MPI_ADDRESS_KIND) :: bytes
    call require_shared_owner(owner)
    bytes = owner%bytes
  end function shared_extent

  subroutine shared_allocation(owner, shape, value_kind, item_bytes)
    type(shared_owner), intent(inout) :: owner
    integer, intent(in) :: shape(:), value_kind, item_bytes
    integer(int64), allocatable :: larger(:)
    integer(int64) :: elements, count_info(6), minimum(6), maximum(6)
    integer(MPI_ADDRESS_KIND) :: local_bytes, queried_bytes
    integer :: k, slot, disp_unit, ierr
    if (.not. mpi_parameters_defined) call init_mpi_parameters()
    if (owner%slot /= 0 .or. owner%window /= MPI_WIN_NULL) &
       call abort_run('allocate_shared: owner already live or stale; release it first')
    if (size(shape) < 1 .or. size(shape) > 3 .or. any(shape < 0)) &
       call abort_run('allocate_shared: invalid shape')
    elements = 1_int64
    do k = 1, size(shape)
       if (shape(k) > 0) then
          if (elements > huge(elements) / int(shape(k), int64)) &
             call abort_run('allocate_shared: element count overflow')
       end if
       elements = elements * int(shape(k), int64)
    end do
    if (elements > int(huge(local_bytes), int64) / int(item_bytes, int64)) &
       call abort_run('allocate_shared: byte extent overflow')
    owner%shape = 0_int64
    owner%shape(1:size(shape)) = int(shape, int64)
    owner%bytes = int(elements, MPI_ADDRESS_KIND) * int(item_bytes, MPI_ADDRESS_KIND)
    count_info = [int(size(shape), int64), int(value_kind, int64), owner%shape, int(owner%bytes, int64)]
    call MPI_ALLREDUCE(count_info, minimum, 6, MPI_INTEGER8, MPI_MIN, hostcomm, ierr)
    call MPI_ALLREDUCE(count_info, maximum, 6, MPI_INTEGER8, MPI_MAX, hostcomm, ierr)
    if (any(minimum /= maximum)) call abort_run('allocate_shared: host ranks disagree on shape or type')
    if (.not. allocated(live_generation)) then
       allocate(live_generation(64))
       live_generation = 0_int64
    end if
    slot = 0
    do k = 1, size(live_generation)
       if (live_generation(k) == 0_int64) then
          slot = k
          exit
       end if
    end do
    if (slot == 0) then
       k = size(live_generation)
       allocate(larger(2 * k))
       larger = 0_int64
       larger(1:k) = live_generation
       call move_alloc(larger, live_generation)
       slot = k + 1
    end if
    if (next_generation == huge(next_generation)) call abort_run('shared owner generation overflow')
    next_generation = next_generation + 1_int64
    owner%slot = slot
    owner%generation = next_generation
    live_generation(slot) = next_generation
    owner%rank = size(shape)
    owner%value_kind = value_kind
    local_bytes = 0_MPI_ADDRESS_KIND
    if (h_rank == 0) local_bytes = owner%bytes
    call MPI_WIN_ALLOCATE_SHARED(local_bytes, item_bytes, MPI_INFO_NULL, hostcomm, &
       owner%base, owner%window, ierr)
    if (h_rank /= 0) then
       call MPI_WIN_SHARED_QUERY(owner%window, 0, queried_bytes, disp_unit, owner%base, ierr)
       ! Intel MPI reports displacement unit zero for an empty segment on
       ! nonleaders. There are no elements whose displacement unit is used.
       if (queried_bytes /= owner%bytes .or. (owner%bytes > 0 .and. disp_unit /= item_bytes)) &
          call abort_run('allocate_shared: unexpected shared segment extent')
    end if
    call MPI_WIN_FENCE(0, owner%window, ierr)
  end subroutine shared_allocation

  subroutine shared_register_sync(owner)
    type(shared_owner), intent(in) :: owner
    integer(c_int8_t), pointer, contiguous, asynchronous :: memory(:)
    if (owner%bytes == 0_MPI_ADDRESS_KIND) return
    call c_f_pointer(owner%base, memory, [owner%bytes])
    ! A contiguous nonpointer section avoids pointer association on the
    ! Fortran buffer argument. Synchronize the full allocation, without packing.
    call MPI_F_SYNC_REG(memory(:))
  end subroutine shared_register_sync

  subroutine publish_shared(owners)
    type(shared_owner), intent(in) :: owners(:)
    integer :: k, ierr
    if (.not. mpi_parameters_defined) call init_mpi_parameters()
    do k = 1, size(owners)
       if (owners(k)%slot == 0) cycle
       call require_shared_owner(owners(k))
       call shared_register_sync(owners(k))
       call MPI_WIN_FENCE(0, owners(k)%window, ierr)
    end do
    call MPI_BARRIER(hostcomm, ierr)
    do k = 1, size(owners)
       if (owners(k)%slot /= 0) call shared_register_sync(owners(k))
    end do
  end subroutine publish_shared

  subroutine shared_release(owner)
    type(shared_owner), intent(inout) :: owner
    integer :: ierr
    call require_shared_owner(owner)
    call shared_register_sync(owner)
    call MPI_WIN_FENCE(MPI_MODE_NOSUCCEED, owner%window, ierr)
    call MPI_WIN_FREE(owner%window, ierr)
    live_generation(owner%slot) = 0_int64
    owner = shared_owner()
  end subroutine shared_release

  subroutine shared_pointer_check(owner, rank, value_kind, shape, base)
    type(shared_owner), intent(in) :: owner
    integer, intent(in) :: rank, value_kind, shape(:)
    type(c_ptr), intent(in) :: base
    call require_shared_owner(owner)
    if (owner%rank /= rank .or. owner%value_kind /= value_kind) &
       call abort_run('free_shared: pointer type or rank does not match owner')
    if (any(owner%shape(1:rank) /= int(shape, int64))) &
       call abort_run('free_shared: pointer shape does not match owner')
    if (owner%bytes > 0 .and. .not. c_associated(owner%base, base)) &
       call abort_run('free_shared: pointer does not belong to owner')
  end subroutine shared_pointer_check

  subroutine allocate_shared_1d_real64(owner, array, shape)
    type(shared_owner), intent(inout) :: owner
    real(real64), pointer, intent(inout) :: array(:)
    integer, intent(in) :: shape(1)
    if (associated(array)) call abort_run('allocate_shared: canonical pointer already associated')
    call shared_allocation(owner, shape, 1, storage_size(0.0_real64) / 8)
    call c_f_pointer(owner%base, array, shape)
    if (h_rank == 0 .and. owner%bytes > 0) array = 0.0_real64
    call publish_shared([owner])
  end subroutine allocate_shared_1d_real64

  subroutine free_shared_1d_real64(owner, array)
    type(shared_owner), intent(inout) :: owner
    real(real64), pointer, intent(inout) :: array(:)
    type(c_ptr) :: base
    if (owner%slot == 0) then
       if (associated(array)) call abort_run('free_shared: associated pointer without owner')
       return
    end if
    call require_shared_owner(owner)
    if (owner%rank /= 1 .or. owner%value_kind /= 1) &
       call abort_run('free_shared: pointer type or rank does not match owner')
    base = c_null_ptr
    if (associated(array)) then
       if (.not. is_contiguous(array)) call abort_run('free_shared: canonical pointer is not contiguous')
       if (owner%bytes > 0) base = c_loc(array(1))
       call shared_pointer_check(owner, 1, 1, shape(array), base)
    else if (owner%bytes > 0) then
       call abort_run('free_shared: missing canonical pointer')
    end if
    call shared_release(owner)
    nullify(array)
  end subroutine free_shared_1d_real64

  subroutine allocate_shared_2d_real64(owner, array, shape)
    type(shared_owner), intent(inout) :: owner
    real(real64), pointer, intent(inout) :: array(:, :)
    integer, intent(in) :: shape(2)
    if (associated(array)) call abort_run('allocate_shared: canonical pointer already associated')
    call shared_allocation(owner, shape, 1, storage_size(0.0_real64) / 8)
    call c_f_pointer(owner%base, array, shape)
    if (h_rank == 0 .and. owner%bytes > 0) array = 0.0_real64
    call publish_shared([owner])
  end subroutine allocate_shared_2d_real64

  subroutine free_shared_2d_real64(owner, array)
    type(shared_owner), intent(inout) :: owner
    real(real64), pointer, intent(inout) :: array(:, :)
    type(c_ptr) :: base
    if (owner%slot == 0) then
       if (associated(array)) call abort_run('free_shared: associated pointer without owner')
       return
    end if
    call require_shared_owner(owner)
    if (owner%rank /= 2 .or. owner%value_kind /= 1) &
       call abort_run('free_shared: pointer type or rank does not match owner')
    base = c_null_ptr
    if (associated(array)) then
       if (.not. is_contiguous(array)) call abort_run('free_shared: canonical pointer is not contiguous')
       if (owner%bytes > 0) base = c_loc(array(1, 1))
       call shared_pointer_check(owner, 2, 1, shape(array), base)
    else if (owner%bytes > 0) then
       call abort_run('free_shared: missing canonical pointer')
    end if
    call shared_release(owner)
    nullify(array)
  end subroutine free_shared_2d_real64

  subroutine allocate_shared_3d_real64(owner, array, shape)
    type(shared_owner), intent(inout) :: owner
    real(real64), pointer, intent(inout) :: array(:, :, :)
    integer, intent(in) :: shape(3)
    if (associated(array)) call abort_run('allocate_shared: canonical pointer already associated')
    call shared_allocation(owner, shape, 1, storage_size(0.0_real64) / 8)
    call c_f_pointer(owner%base, array, shape)
    if (h_rank == 0 .and. owner%bytes > 0) array = 0.0_real64
    call publish_shared([owner])
  end subroutine allocate_shared_3d_real64

  subroutine free_shared_3d_real64(owner, array)
    type(shared_owner), intent(inout) :: owner
    real(real64), pointer, intent(inout) :: array(:, :, :)
    type(c_ptr) :: base
    if (owner%slot == 0) then
       if (associated(array)) call abort_run('free_shared: associated pointer without owner')
       return
    end if
    call require_shared_owner(owner)
    if (owner%rank /= 3 .or. owner%value_kind /= 1) &
       call abort_run('free_shared: pointer type or rank does not match owner')
    base = c_null_ptr
    if (associated(array)) then
       if (.not. is_contiguous(array)) call abort_run('free_shared: canonical pointer is not contiguous')
       if (owner%bytes > 0) base = c_loc(array(1, 1, 1))
       call shared_pointer_check(owner, 3, 1, shape(array), base)
    else if (owner%bytes > 0) then
       call abort_run('free_shared: missing canonical pointer')
    end if
    call shared_release(owner)
    nullify(array)
  end subroutine free_shared_3d_real64

  subroutine allocate_shared_1d_int32(owner, array, shape)
    type(shared_owner), intent(inout) :: owner
    integer(int32), pointer, intent(inout) :: array(:)
    integer, intent(in) :: shape(1)
    if (associated(array)) call abort_run('allocate_shared: canonical pointer already associated')
    call shared_allocation(owner, shape, 2, storage_size(0_int32) / 8)
    call c_f_pointer(owner%base, array, shape)
    if (h_rank == 0 .and. owner%bytes > 0) array = 0_int32
    call publish_shared([owner])
  end subroutine allocate_shared_1d_int32

  subroutine free_shared_1d_int32(owner, array)
    type(shared_owner), intent(inout) :: owner
    integer(int32), pointer, intent(inout) :: array(:)
    type(c_ptr) :: base
    if (owner%slot == 0) then
       if (associated(array)) call abort_run('free_shared: associated pointer without owner')
       return
    end if
    call require_shared_owner(owner)
    if (owner%rank /= 1 .or. owner%value_kind /= 2) &
       call abort_run('free_shared: pointer type or rank does not match owner')
    base = c_null_ptr
    if (associated(array)) then
       if (.not. is_contiguous(array)) call abort_run('free_shared: canonical pointer is not contiguous')
       if (owner%bytes > 0) base = c_loc(array(1))
       call shared_pointer_check(owner, 1, 2, shape(array), base)
    else if (owner%bytes > 0) then
       call abort_run('free_shared: missing canonical pointer')
    end if
    call shared_release(owner)
    nullify(array)
  end subroutine free_shared_1d_int32

  subroutine allocate_shared_2d_int32(owner, array, shape)
    type(shared_owner), intent(inout) :: owner
    integer(int32), pointer, intent(inout) :: array(:, :)
    integer, intent(in) :: shape(2)
    if (associated(array)) call abort_run('allocate_shared: canonical pointer already associated')
    call shared_allocation(owner, shape, 2, storage_size(0_int32) / 8)
    call c_f_pointer(owner%base, array, shape)
    if (h_rank == 0 .and. owner%bytes > 0) array = 0_int32
    call publish_shared([owner])
  end subroutine allocate_shared_2d_int32

  subroutine free_shared_2d_int32(owner, array)
    type(shared_owner), intent(inout) :: owner
    integer(int32), pointer, intent(inout) :: array(:, :)
    type(c_ptr) :: base
    if (owner%slot == 0) then
       if (associated(array)) call abort_run('free_shared: associated pointer without owner')
       return
    end if
    call require_shared_owner(owner)
    if (owner%rank /= 2 .or. owner%value_kind /= 2) &
       call abort_run('free_shared: pointer type or rank does not match owner')
    base = c_null_ptr
    if (associated(array)) then
       if (.not. is_contiguous(array)) call abort_run('free_shared: canonical pointer is not contiguous')
       if (owner%bytes > 0) base = c_loc(array(1, 1))
       call shared_pointer_check(owner, 2, 2, shape(array), base)
    else if (owner%bytes > 0) then
       call abort_run('free_shared: missing canonical pointer')
    end if
    call shared_release(owner)
    nullify(array)
  end subroutine free_shared_2d_int32

  ! Called after every module has released its canonical owner/pointer pairs.
  ! This is a leak check and communicator release, never an implicit array free.
  subroutine shared_memory_finish()
    integer :: ierr
    if (allocated(live_generation)) then
       if (any(live_generation /= 0_int64)) call abort_run('shared_memory_finish: live owners remain')
       deallocate(live_generation)
    end if
    if (mpi_parameters_defined) then
       call MPI_COMM_FREE(same_hrank_comm, ierr)
       call MPI_COMM_FREE(hostcomm, ierr)
       mpi_parameters_defined = .false.
    end if
  end subroutine shared_memory_finish

  subroutine shared_sum_shape(owner, actual_shape)
    type(shared_owner), intent(in) :: owner
    integer, intent(in) :: actual_shape(:)
    call require_shared_owner(owner)
    if (owner%value_kind /= 1 .or. owner%rank /= size(actual_shape)) &
       call abort_run('sum_into_shared: result kind or rank mismatch')
    if (any(owner%shape(1:owner%rank) /= int(actual_shape, int64))) &
       call abort_run('sum_into_shared: result shape mismatch')
  end subroutine shared_sum_shape

  subroutine shared_sum_real64(scores, n, owner)
    integer(int64), intent(in) :: n
    real(real64), target, intent(in) :: scores(n)
    type(shared_owner), intent(in) :: owner
    real(real64), pointer, contiguous :: result(:)
    ! This entry also orders the last readers of the previous result.
    call publish_shared([owner])
    if (n > 0_int64) then
       if (c_associated(c_loc(scores(1)), owner%base)) &
          call abort_run('sum_into_shared: scores must be private, not the shared result')
       call c_f_pointer(owner%base, result, [n])
       call sum_real64(result, n, SUM_SHARED_RESULT, scores)
    end if
    call publish_shared([owner])
  end subroutine shared_sum_real64

  subroutine sum_into_shared_1d_real64(scores, owner)
    real(real64), contiguous, target, intent(in) :: scores(:)
    type(shared_owner), intent(in) :: owner
    call shared_sum_shape(owner, shape(scores))
    call shared_sum_real64(scores, size(scores, kind = int64), owner)
  end subroutine sum_into_shared_1d_real64

  subroutine sum_into_shared_2d_real64(scores, owner)
    real(real64), contiguous, target, intent(in) :: scores(:, :)
    type(shared_owner), intent(in) :: owner
    call shared_sum_shape(owner, shape(scores))
    call shared_sum_real64(scores, size(scores, kind = int64), owner)
  end subroutine sum_into_shared_2d_real64

  subroutine sum_into_shared_3d_real64(scores, owner)
    real(real64), contiguous, target, intent(in) :: scores(:, :, :)
    type(shared_owner), intent(in) :: owner
    call shared_sum_shape(owner, shape(scores))
    call shared_sum_real64(scores, size(scores, kind = int64), owner)
  end subroutine sum_into_shared_3d_real64

  ! Combine disjoint writes in whole shared allocations across hosts.
  ! Collective over all ranks. Foreign-host contributions must be zero.
  ! No within-host reduction: all host ranks see the same allocation.
  subroutine sum_shared_hosts(owners)
    type(shared_owner), intent(in) :: owners(:)
    real(real64), pointer, contiguous :: values(:)
    integer(int64) :: n
    integer :: k
    call publish_shared(owners)
    do k = 1, size(owners)
       if (owners(k)%slot == 0) cycle
       call require_shared_owner(owners(k))
       if (owners(k)%value_kind /= 1) call abort_run('sum_shared_hosts: result must be real64')
       n = product(owners(k)%shape(1:owners(k)%rank))
       if (n == 0_int64) cycle
       call c_f_pointer(owners(k)%base, values, [n])
       call sum_real64(values, n, SUM_SHARED_RESULT, host_shared = .true.)
    end do
    call publish_shared(owners)
  end subroutine sum_shared_hosts

  ! Sum private contributions within hosts, then between host leaders, or
  ! combine an already shared host contribution without a local reduction.
  ! Only sum_to_all_ranks broadcasts the result back within each host.
  ! MPI ignores the receive argument on nonroot ranks of MPI_REDUCE.
  subroutine sum_real64(buf, n, destination, contribution, host_shared)
    integer(kind = int64), intent(in) :: n
    integer, intent(in) :: destination
    real(kind = real64), intent(inout) :: buf(n)
    real(kind = real64), intent(in), optional :: contribution(n)
    logical, intent(in), optional :: host_shared
    logical :: shared_host_values
    real(kind = real64) :: unused(1)
    integer(kind = int64) :: i0
    integer :: m, ierr

    if (n < 0_int64) call abort_run('negative element count in sum_real64')
    if (.not. mpi_parameters_defined) call init_mpi_parameters()
    shared_host_values = .false.
    if (present(host_shared)) shared_host_values = host_shared
    i0 = 1_int64
    do while (i0 <= n)
       m = int(min(SUM_CHUNK, n - i0 + 1_int64))
       if (h_rank == 0) then
          if (.not. shared_host_values) then
             if (present(contribution)) then
                call MPI_REDUCE(contribution(i0), buf(i0), m, MPI_REAL8, MPI_SUM, 0, hostcomm, ierr)
             else
                call MPI_REDUCE(MPI_IN_PLACE, buf(i0), m, MPI_REAL8, MPI_SUM, 0, hostcomm, ierr)
             end if
          end if
          if (num_hosts > 1) then
             select case (destination)
             case (SUM_ROOT)
                if (p_rank == 0) then
                   call MPI_REDUCE(MPI_IN_PLACE, buf(i0), m, MPI_REAL8, MPI_SUM, 0, same_hrank_comm, ierr)
                else
                   call MPI_REDUCE(buf(i0), unused, m, MPI_REAL8, MPI_SUM, 0, same_hrank_comm, ierr)
                end if
             case (SUM_ALL_RANKS, SUM_SHARED_RESULT)
                call MPI_ALLREDUCE(MPI_IN_PLACE, buf(i0), m, MPI_REAL8, MPI_SUM, same_hrank_comm, ierr)
             end select
          end if
       else if (.not. shared_host_values) then
          if (present(contribution)) then
             call MPI_REDUCE(contribution(i0), unused, m, MPI_REAL8, MPI_SUM, 0, hostcomm, ierr)
          else
             call MPI_REDUCE(buf(i0), unused, m, MPI_REAL8, MPI_SUM, 0, hostcomm, ierr)
          end if
       end if
       if (destination == SUM_ALL_RANKS) call MPI_BCAST(buf(i0), m, MPI_REAL8, 0, hostcomm, ierr)
       i0 = i0 + int(m, int64)
    end do
  end subroutine sum_real64

  subroutine sum_to_all_ranks_1D_real64(array)
    real(kind = real64), contiguous, intent(inout) :: array(:)
    call sum_real64(array, size(array, kind = int64), SUM_ALL_RANKS)
  end subroutine sum_to_all_ranks_1D_real64

  subroutine sum_to_all_ranks_2D_real64(array)
    real(kind = real64), contiguous, intent(inout) :: array(:, :)
    call sum_real64(array, size(array, kind = int64), SUM_ALL_RANKS)
  end subroutine sum_to_all_ranks_2D_real64


  subroutine sum_to_root_1D_real64(array)
    real(kind = real64), contiguous, intent(inout) :: array(:)
    call sum_real64(array, size(array, kind = int64), SUM_ROOT)
  end subroutine sum_to_root_1D_real64

  subroutine sum_to_root_2D_real64(array)
    real(kind = real64), contiguous, intent(inout) :: array(:, :)
    call sum_real64(array, size(array, kind = int64), SUM_ROOT)
  end subroutine sum_to_root_2D_real64

  subroutine sum_to_root_3D_real64(array)
    real(kind = real64), contiguous, intent(inout) :: array(:, :, :)
    call sum_real64(array, size(array, kind = int64), SUM_ROOT)
  end subroutine sum_to_root_3D_real64

  subroutine sum_to_root_4D_real64(array)
    real(kind = real64), contiguous, intent(inout) :: array(:, :, :, :)
    call sum_real64(array, size(array, kind = int64), SUM_ROOT)
  end subroutine sum_to_root_4D_real64

end module memory_mod
