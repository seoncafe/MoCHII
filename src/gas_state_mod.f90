module gas_state_mod
!---------------------------------------------------------------------------
! MoCHII: persistent gas leaf state.
!
! nH/xHI are read from the AMR file and used for the dust-density step;
! MoCHII persists them as shared-memory leaf arrays and adds the He
! ionization fractions and T_e/n_e.
!
! Hydrogen's state is carried as BOTH fractions, x_HI and x_HII, each to its
! own relative precision; x_HI + x_HII = 1 to rounding.  Beyond the
! ionization front x_HII falls to 1e-10 .. 1e-17 (and below), where
! 1 - x_HI keeps none of its digits or rounds to zero, so every consumer of
! n_HII reads gas_xHII instead of forming it from x_HI.  The cell solves
! produce both (ion_balance_mod::hydrogen_ionization_balance).
!
! Helium likewise carries all three fractions, x_HeI, x_HeII, x_HeIII, each
! to its own relative precision (x_HeI + x_HeII + x_HeIII = 1 to rounding):
! the cell solve forms them from the stage ratios, and in radiation-starved
! gas 1 - x_HeI - x_HeII rounds to one unit of rounding (1.1e-16) where
! x_HeIII is << 1e-30, which the He2+ + H exchange term then carries.
!
! The state is INITIALIZED (from the file columns where present,
! else from par%xHI_init etc.); when gas_niter = 0 it is never updated —
! no equilibrium solve, no opacity feedback.  gas_state_setup is called
! from grid_create_amr before the transient read arrays are deallocated.
!
! Where the cell solves run.  The state is one node-shared copy per node;
! the ionization/thermal solves of the leaves are divided over ALL ranks of
! MPI_COMM_WORLD, rank p taking the leaves gas_leaf_owner_first() +
! k gas_leaf_owner_stride() (cyclic, as the SEDust solve and the dust
! emission launch divide theirs), so adding nodes divides the work.  Each
! rank fills its own leaves of private arrays, zero elsewhere, and
! gas_state_publish sums them over MPI_COMM_WORLD: every leaf has exactly one
! nonzero contribution, so the sum is exact and the same on every rank
! whatever the reduction order, and each node's rank 0 then writes it into
! that node's shared copy.  Until 2026-09-26 the thermal solve was divided
! over the ranks of one node only, every node repeating all leaves, and the
! fixed-Te solve ran on each node's rank 0 alone.
!---------------------------------------------------------------------------
  use define
  use memory_mod, only : create_shared_mem, get_window
  use mpi
  implicit none
  private

  public :: gas_state_setup, gas_state_recreate, gas_state_sync
  public :: gas_state_publish, gas_leaf_owner_first, gas_leaf_owner_stride
  public :: gas_nH, gas_xHI, gas_xHII, gas_xHeI, gas_xHeII, gas_xHeIII
  public :: gas_ne, gas_Te
  public :: gas_nleaf

  real(kind=wp), pointer :: gas_nH(:)    => null()  ! H number density [cm^-3]
  real(kind=wp), pointer :: gas_xHI(:)   => null()  ! n_HI / n_H
  real(kind=wp), pointer :: gas_xHII(:)  => null()  ! n_HII / n_H
  real(kind=wp), pointer :: gas_xHeI(:)  => null()  ! n_HeI / n_He
  real(kind=wp), pointer :: gas_xHeII(:) => null()  ! n_HeII / n_He
  real(kind=wp), pointer :: gas_xHeIII(:) => null() ! n_HeIII / n_He
  real(kind=wp), pointer :: gas_ne(:)    => null()  ! electron density [cm^-3]
  real(kind=wp), pointer :: gas_Te(:)    => null()  ! electron temperature [K]
  integer :: gas_nleaf = 0
  integer :: win_xHI = MPI_WIN_NULL, win_xHII = MPI_WIN_NULL
  integer :: win_xHeI = MPI_WIN_NULL
  integer :: win_xHeII = MPI_WIN_NULL, win_xHeIII = MPI_WIN_NULL
  integer :: win_ne = MPI_WIN_NULL
  integer :: win_Te = MPI_WIN_NULL

contains

  !=========================================================================
  ! nH is mandatory (from the AMR file); xHI is optional (file column) —
  ! when absent every leaf starts at par%xHI_init.  x_HII starts at 1 - x_HI
  ! (the only information the file or the namelist gives).  He fractions
  ! start at par%xHeI_init / par%xHeII_init (no He columns in the generic
  ! format), x_HeIII at max(0, 1 - xHeI_init - xHeII_init).
  !=========================================================================
  subroutine gas_state_setup(nH, nleaf, xHI)
    implicit none
    real(kind=wp), intent(in)           :: nH(:)
    integer,       intent(in)           :: nleaf
    real(kind=wp), intent(in), optional :: xHI(:)
    integer :: il

    gas_nleaf = nleaf
    call create_shared_mem(gas_nH,    [nleaf])
    call create_shared_mem(gas_xHI,   [nleaf])
    call create_shared_mem(gas_xHII,  [nleaf])
    call create_shared_mem(gas_xHeI,  [nleaf])
    call create_shared_mem(gas_xHeII, [nleaf])
    call create_shared_mem(gas_xHeIII, [nleaf])
    call create_shared_mem(gas_ne,    [nleaf])
    call create_shared_mem(gas_Te,    [nleaf])
    call cache_state_windows()

    if (mpar%h_rank == 0) then
       do il = 1, nleaf
          gas_nH(il) = nH(il)
       end do
       if (present(xHI)) then
          do il = 1, nleaf
             gas_xHI(il) = xHI(il)
          end do
       else
          gas_xHI(:) = par%xHI_init
       end if
       do il = 1, nleaf
          gas_xHII(il) = 1.0_wp - gas_xHI(il)
       end do
       gas_xHeI(:)  = par%xHeI_init
       gas_xHeII(:) = par%xHeII_init
       gas_xHeIII(:) = max(0.0_wp, 1.0_wp - par%xHeI_init - par%xHeII_init)
       do il = 1, nleaf
          gas_ne(il) = gas_nH(il) * ( gas_xHII(il) &
                       + par%abund_He*(gas_xHeII(il) + 2.0_wp*gas_xHeIII(il)) )
       end do
       gas_Te(:) = par%te_fixed
    end if
    call gas_state_sync()

    if (mpar%p_rank == 0) then
       write(*,'(a,i12)')    ' GAS: leaf state allocated, nleaf = ', nleaf
       if (present(xHI)) then
          write(*,'(a)')     ' GAS: x_HI initialized from the AMR file column'
       else
          write(*,'(a,f8.4)')' GAS: x_HI initialized to ', par%xHI_init
       end if
       write(*,'(a,2f8.4)')  ' GAS: x_HeI, x_HeII initialized to ', &
                             par%xHeI_init, par%xHeII_init
       write(*,'(a,f8.4)')   ' GAS: He abundance n_He/n_H = ', par%abund_He
    end if
  end subroutine gas_state_setup

  !=========================================================================
  ! Recreate the state on a re-refined tree: fresh shared windows
  ! (the old ones are leaked until finalize), filled from position-mapped
  ! local arrays.  Collective; h_rank 0 writes.
  !=========================================================================
  subroutine gas_state_recreate(nleaf, nHv, x1, x1p, x2, x3, x4, nev, tev)
    implicit none
    integer,       intent(in) :: nleaf
    real(kind=wp), intent(in) :: nHv(nleaf), x1(nleaf), x1p(nleaf), &
                                 x2(nleaf), x3(nleaf), x4(nleaf), &
                                 nev(nleaf), tev(nleaf)
    gas_nleaf = nleaf
    call create_shared_mem(gas_nH,    [nleaf])
    call create_shared_mem(gas_xHI,   [nleaf])
    call create_shared_mem(gas_xHII,  [nleaf])
    call create_shared_mem(gas_xHeI,  [nleaf])
    call create_shared_mem(gas_xHeII, [nleaf])
    call create_shared_mem(gas_xHeIII, [nleaf])
    call create_shared_mem(gas_ne,    [nleaf])
    call create_shared_mem(gas_Te,    [nleaf])
    call cache_state_windows()
    if (mpar%h_rank == 0) then
       gas_nH(1:nleaf)    = nHv
       gas_xHI(1:nleaf)   = x1
       gas_xHII(1:nleaf)  = x1p
       gas_xHeI(1:nleaf)  = x2
       gas_xHeII(1:nleaf) = x3
       gas_xHeIII(1:nleaf) = x4
       gas_ne(1:nleaf)    = nev
       gas_Te(1:nleaf)    = tev
    end if
    call gas_state_sync()
  end subroutine gas_state_recreate

  !--- the first leaf this rank solves, and the step to its next one
  integer function gas_leaf_owner_first() result(i)
    i = mpar%p_rank + 1
  end function gas_leaf_owner_first

  integer function gas_leaf_owner_stride() result(i)
    i = mpar%nproc
  end function gas_leaf_owner_stride

  !=========================================================================
  ! Sum the ranks' solved leaves over MPI_COMM_WORLD (each array holds this
  ! rank's own leaves, zero elsewhere) and write the result into every
  ! node's shared state.  te_new is absent for the fixed-T_e solve, which
  ! leaves T_e as it is.  The arrays come back holding the full state.
  !=========================================================================
  subroutine gas_state_publish(xHI_new, xHII_new, xHeI_new, xHeII_new, &
                               xHeIII_new, ne_new, te_new)
    real(kind=wp), intent(inout)           :: xHI_new(:), xHII_new(:), &
                                              xHeI_new(:), xHeII_new(:), &
                                              xHeIII_new(:), ne_new(:)
    real(kind=wp), intent(inout), optional :: te_new(:)
    call sum_over_world(xHI_new)
    call sum_over_world(xHII_new)
    call sum_over_world(xHeI_new)
    call sum_over_world(xHeII_new)
    call sum_over_world(xHeIII_new)
    call sum_over_world(ne_new)
    if (present(te_new)) call sum_over_world(te_new)
    if (mpar%h_rank == 0) then
       gas_xHI(1:gas_nleaf)   = xHI_new(1:gas_nleaf)
       gas_xHII(1:gas_nleaf)  = xHII_new(1:gas_nleaf)
       gas_xHeI(1:gas_nleaf)  = xHeI_new(1:gas_nleaf)
       gas_xHeII(1:gas_nleaf) = xHeII_new(1:gas_nleaf)
       gas_xHeIII(1:gas_nleaf) = xHeIII_new(1:gas_nleaf)
       gas_ne(1:gas_nleaf)    = ne_new(1:gas_nleaf)
       if (present(te_new)) gas_Te(1:gas_nleaf) = te_new(1:gas_nleaf)
    end if
    call gas_state_sync()
  end subroutine gas_state_publish

  !--- in-place MPI_SUM over MPI_COMM_WORLD, in chunks of at most 1e8 values
  subroutine sum_over_world(a)
    real(kind=wp), intent(inout) :: a(:)
    integer, parameter :: NCHUNK = 100000000
    integer :: i0, n, ierr
    i0 = 1
    do while (i0 <= size(a))
       n = min(NCHUNK, size(a) - i0 + 1)
       call MPI_ALLREDUCE(MPI_IN_PLACE, a(i0:i0+n-1), n, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
       i0 = i0 + n
    end do
  end subroutine sum_over_world

  !=========================================================================
  ! Publish direct loads/stores to the node-shared gas-state windows.
  !
  ! The inherited shared-memory allocator keeps a fence epoch open for each
  ! window.  Advancing that epoch collectively publishes the direct stores;
  ! the host barrier then keeps all state windows at the same update boundary.
  !=========================================================================
  subroutine gas_state_sync()
    integer :: ierr

    if (win_xHI   /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, win_xHI,   ierr)
    if (win_xHII  /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, win_xHII,  ierr)
    if (win_xHeI  /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, win_xHeI,  ierr)
    if (win_xHeII /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, win_xHeII, ierr)
    if (win_xHeIII /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, win_xHeIII, ierr)
    if (win_ne    /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, win_ne,    ierr)
    if (win_Te    /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, win_Te,    ierr)
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine gas_state_sync

  subroutine cache_state_windows()
    call get_window(gas_xHI,   win_xHI)
    call get_window(gas_xHII,  win_xHII)
    call get_window(gas_xHeI,  win_xHeI)
    call get_window(gas_xHeII, win_xHeII)
    call get_window(gas_xHeIII, win_xHeIII)
    call get_window(gas_ne,    win_ne)
    call get_window(gas_Te,    win_Te)
  end subroutine cache_state_windows

end module gas_state_mod
