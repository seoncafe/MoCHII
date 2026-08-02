module photon_schedule_mod
!---------------------------------------------------------------------------
! MoCHII: how the global photon indices of one transport pass are split over
! the MPI ranks.  This is infrastructure, not physics -- the loop body is the
! same either way, and only *which* rank draws *which* global index changes:
!
!     call photon_schedule_begin(ntotal)
!     do while (photon_schedule_next(ip))
!        ... transport packet ip ...
!     end do
!     call photon_schedule_end()
!
! Two schedules sit behind that interface, selected by par%use_master_slave.
!
!   static (.false.)    Every rank walks its own arithmetic sequence
!                       ip = p_rank+1, p_rank+1+nproc, ... -- the same indices
!                       in the same order as the
!                          do ip = mpar%p_rank+1, ntotal, mpar%nproc
!                       it replaces, so results are unchanged bit for bit.
!
!   dynamic (.true.)    Rank 0 is the master: it serves batches of
!                       par%num_send_at_once consecutive indices to whichever
!                       rank asks next, so a rank running slow takes a smaller
!                       share instead of holding up the reduction that closes
!                       the pass.  The master transports nothing -- its first
!                       call to photon_schedule_next runs the whole dispatch
!                       loop and then reports .false.
!
! Why: measured on hii40_bench at 70 ranks, every rank carries the same photon
! count by construction and the same path-segment count to 0.17%, yet the ranks
! take 17.1-35.8 s to do it (hyperthread contention on a 36-core node).  The
! static schedule pays the maximum of those times; served on demand the same
! pass costs close to their harmonic mean (docs/MPI_LOAD_BALANCE_2026-08-02.md).
! Measured here on the same case, the transport pass drops from a median 52.1 s
! to 43.4 s at 70 ranks, and its run-to-run spread from 28% to 1.5%.
!
! Reproducibility: with par%launch_sequence = 'sobol' a packet's launch is
! fixed by its global index, so the two schedules transport the *same* set of
! packets and differ only in the order the contributions are summed.  With the
! default 'random' each rank draws from its own stream, so the schedule also
! selects the realization: statistically equivalent, not numerically identical.
!
! Deadlock boundaries handled here: nproc = 1 (no worker exists, so a dynamic
! request degrades to the static schedule and rank 0 transports everything),
! fewer batches than workers, and ntotal = 0 (a pass with no diffuse packets),
! where every worker is released immediately.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: photon_schedule_begin, photon_schedule_next, photon_schedule_end

  !--- message tags of the master-worker protocol.  These are the only
  !--- point-to-point messages in MoCHII; everything else is collective.
  integer, parameter :: TAG_WORK = 1     ! payload = last index of a batch
  integer, parameter :: TAG_STOP = 0     ! no more work for this rank
  integer, parameter :: MASTER   = 0

  !--- schedule state of the pass in progress; set by photon_schedule_begin
  !--- and advanced by photon_schedule_next.
  logical             :: sched_dynamic = .false.  ! serve batches on demand
  logical             :: sched_master  = .false.  ! this rank dispatches only
  logical             :: sched_holding = .false.  ! a batch is unacknowledged
  logical             :: sched_done    = .false.  ! nothing left for this rank
  integer(kind=int64) :: sched_ntotal  = 0        ! packets in the pass
  integer(kind=int64) :: sched_ip      = 0        ! next index to hand back
  integer(kind=int64) :: sched_hi      = -1       ! last index of the batch
  integer(kind=int64) :: sched_report  = 0        ! progress stride, 0 = quiet
  !--- batch length in use, taken once from par%num_send_at_once so that the
  !--- master and the workers cannot read a different value for it.
  integer(kind=int64) :: sched_nsend   = 1

contains

  !=========================================================================
  ! Open a pass of ntotal packets.  report_stride > 0 asks the dynamic master
  ! to print the same progress line the static schedule prints from the loop
  ! body, at multiples of that many dispatched packets; the static schedule
  ! ignores it (rank 0 still reports its own count there).
  !=========================================================================
  subroutine photon_schedule_begin(ntotal, report_stride)
    implicit none
    integer(kind=int64), intent(in)           :: ntotal
    integer(kind=int64), intent(in), optional :: report_stride

    sched_ntotal = max(ntotal, 0_int64)
    sched_report = 0_int64
    if (present(report_stride)) sched_report = max(report_stride, 0_int64)
    !--- a served schedule needs at least one rank besides the master.
    sched_dynamic = par%use_master_slave .and. mpar%nproc > 1
    sched_master  = sched_dynamic .and. mpar%p_rank == MASTER
    sched_nsend   = max(int(par%num_send_at_once, int64), 1_int64)
    sched_done    = .false.
    sched_holding = .false.
    sched_hi      = -1_int64
    if (sched_dynamic) then
       sched_ip = 1_int64                              ! set from the batch
    else
       sched_ip = int(mpar%p_rank, int64) + 1_int64    ! cyclic sequence start
    end if
  end subroutine photon_schedule_begin

  !=========================================================================
  ! Next global packet index for this rank, .false. when the rank is finished.
  !=========================================================================
  function photon_schedule_next(ip) result(more)
    use mpi
    implicit none
    integer(kind=int64), intent(out) :: ip
    logical :: more
    integer :: ierr, ans
    integer :: status(MPI_STATUS_SIZE)
    integer(kind=int64) :: batch_end

    ip   = 0_int64
    more = .false.
    if (sched_done) return

    !--- static: the arithmetic sequence, identical to the loop it replaces.
    if (.not. sched_dynamic) then
       if (sched_ip > sched_ntotal) then
          sched_done = .true.
          return
       end if
       ip       = sched_ip
       sched_ip = sched_ip + int(mpar%nproc, int64)
       more     = .true.
       return
    end if

    !--- dynamic master: hand out every batch here, then transport nothing.
    if (sched_master) then
       call dispatch_batches()
       sched_done = .true.
       return
    end if

    !--- dynamic worker: walk the held batch, then acknowledge and ask again.
    do
       if (sched_holding .and. sched_ip <= sched_hi) then
          ip       = sched_ip
          sched_ip = sched_ip + 1_int64
          more     = .true.
          return
       end if
       if (sched_holding) then
          ans = 1
          call MPI_SEND(ans, 1, MPI_INTEGER, MASTER, TAG_WORK, &
                        MPI_COMM_WORLD, ierr)
          sched_holding = .false.
       end if
       call MPI_RECV(batch_end, 1, MPI_INTEGER8, MASTER, MPI_ANY_TAG, &
                     MPI_COMM_WORLD, status, ierr)
       if (status(MPI_TAG) == TAG_STOP) then
          sched_done = .true.
          return
       end if
       !--- the master sends the running total, so the batch is the
       !--- num_send_at_once indices ending there, clipped at ntotal (the
       !--- last batch is partial, and a batch can be empty when there are
       !--- fewer indices than the initial hand-out assumed).
       sched_ip      = batch_end - sched_nsend + 1_int64
       sched_hi      = min(batch_end, sched_ntotal)
       sched_holding = .true.
    end do
  end function photon_schedule_next

  !=========================================================================
  ! Close the pass.  The dynamic protocol is already matched send for receive
  ! when the master leaves its dispatch loop; the barrier only keeps the two
  ! dispatch phases of one pass (stellar, then diffuse) from overlapping.
  !=========================================================================
  subroutine photon_schedule_end()
    use mpi
    implicit none
    integer :: ierr

    if (sched_dynamic) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    sched_done    = .true.
    sched_holding = .false.
  end subroutine photon_schedule_end

  !=========================================================================
  ! The master's dispatch loop: one initial batch to each worker that has one
  ! waiting, immediate release for the rest, then one reply per acknowledgment
  ! until every batch has been worked and every worker released.
  !=========================================================================
  subroutine dispatch_batches()
    use mpi
    use utility, only : time_stamp
    implicit none
    integer :: ierr, irank, worker, ans, nworkers, n_init
    integer :: status(MPI_STATUS_SIZE)
    integer(kind=int64) :: nsend, nbatch, numsent, ndone, next_report, nout
    real(kind=wp) :: dtime

    nsend    = sched_nsend
    nworkers = mpar%nproc - 1
    nbatch   = (sched_ntotal + nsend - 1_int64)/nsend
    n_init   = int(min(int(nworkers, int64), nbatch))

    numsent = 0_int64
    do irank = 1, n_init
       numsent = numsent + nsend
       call MPI_SEND(numsent, 1, MPI_INTEGER8, irank, TAG_WORK, &
                     MPI_COMM_WORLD, ierr)
    end do
    !--- workers with no initial batch (fewer batches than workers, ntotal = 0
    !--- included) are released at once; otherwise they block forever in a
    !--- receive the master would never answer.
    do irank = n_init + 1, nworkers
       call MPI_SEND(numsent, 1, MPI_INTEGER8, irank, TAG_STOP, &
                     MPI_COMM_WORLD, ierr)
    end do

    ndone       = 0_int64
    next_report = sched_report
    nout        = int(n_init, int64)
    do while (nout > 0_int64)
       call MPI_RECV(ans, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                     MPI_COMM_WORLD, status, ierr)
       worker = status(MPI_SOURCE)
       nout   = nout - 1_int64
       !--- cap the printed counter at ntotal (the last batch is partial)
       ndone  = min(ndone + nsend, sched_ntotal)
       if (sched_report > 0_int64 .and. ndone >= next_report) then
          call time_stamp(dtime)
          write(6,'(es14.3,a,f8.3,a)') real(ndone,wp), ' photons  @ ', &
             dtime/60.0_wp, ' mins'
          do while (next_report <= ndone)
             next_report = next_report + sched_report
          end do
       end if
       if (numsent < sched_ntotal) then
          numsent = numsent + nsend
          call MPI_SEND(numsent, 1, MPI_INTEGER8, worker, TAG_WORK, &
                        MPI_COMM_WORLD, ierr)
          nout = nout + 1_int64
       else
          call MPI_SEND(numsent, 1, MPI_INTEGER8, worker, TAG_STOP, &
                        MPI_COMM_WORLD, ierr)
       end if
    end do
  end subroutine dispatch_batches

end module photon_schedule_mod
