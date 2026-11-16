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
! rank writes its leaves into separate node-shared new-state buffers.
! Leaves owned by another host remain zero. After all old-state readers
! and change measures finish, gas_state_publish combines the host buffers
! among their leaders, copies into the unchanged canonical current-state
! allocations, and publishes the entire update. Packet scoring stays private.
!---------------------------------------------------------------------------
  use define
  use utility, only : abort_run
  use memory_mod, only : shared_owner, allocate_shared, free_shared, publish_shared, shared_is_live, &
                        sum_shared_hosts
  implicit none
  private

  public :: gas_state_release, gas_state_new_reset
  public :: xHI_new, xHII_new, xHeI_new, xHeII_new, xHeIII_new, ne_new, te_new

  public :: gas_state_setup, gas_state_recreate, gas_state_sync
  public :: gas_state_publish, gas_leaf_owner_first, gas_leaf_owner_stride
  public :: gas_nH, gas_xHI, gas_xHII, gas_xHeI, gas_xHeII, gas_xHeIII
  public :: gas_ne, gas_Te
  public :: gas_nleaf
  public :: helium_starts_as_hydrogen

  real(kind=wp), pointer :: gas_nH(:)    => null()  ! H number density [cm^-3]
  real(kind=wp), pointer :: gas_xHI(:)   => null()  ! n_HI / n_H
  real(kind=wp), pointer :: gas_xHII(:)  => null()  ! n_HII / n_H
  real(kind=wp), pointer :: gas_xHeI(:)  => null()  ! n_HeI / n_He
  real(kind=wp), pointer :: gas_xHeII(:) => null()  ! n_HeII / n_He
  real(kind=wp), pointer :: gas_xHeIII(:) => null() ! n_HeIII / n_He
  real(kind=wp), pointer :: gas_ne(:)    => null()  ! electron density [cm^-3]
  real(kind=wp), pointer :: gas_Te(:)    => null()  ! electron temperature [K]
  integer :: gas_nleaf = 0
  !--- set by setup.f90 (read_parameters_with_implicit_defaults): conv_crit =
  !--- 'cell_fraction' with xHeI_init and xHeII_init both left out of the
  !--- input, so gas_state_setup starts helium as hydrogen in each leaf
  logical :: helium_starts_as_hydrogen = .false.
  type(shared_owner) :: state_owner(8)
  ! Separate canonical allocations preserve the old state during every solve.
  real(kind=wp), pointer, contiguous :: xHI_new(:) => null(), xHII_new(:) => null()
  real(kind=wp), pointer, contiguous :: xHeI_new(:) => null(), xHeII_new(:) => null()
  real(kind=wp), pointer, contiguous :: xHeIII_new(:) => null(), ne_new(:) => null()
  real(kind=wp), pointer, contiguous :: te_new(:) => null()
  type(shared_owner) :: new_state_owner(7)
  logical :: new_state_ready = .false.

contains

  !=========================================================================
  ! nH is mandatory (from the AMR file); xHI is optional (file column) —
  ! when absent every leaf starts at par%xHI_init.  x_HII starts at 1 - x_HI
  ! (the only information the file or the namelist gives).  He fractions
  ! start at par%xHeI_init / par%xHeII_init (no He columns in the generic
  ! format), x_HeIII at max(0, 1 - xHeI_init - xHeII_init); with
  ! helium_starts_as_hydrogen (conv_crit = 'cell_fraction', both He items
  ! left out of the input) at x_HeI = x_HI, x_HeII = x_HII, x_HeIII = 0 in
  ! each leaf, as MOCASSIN starts.  T_e starts at par%te_start with
  ! solve_te and is held at par%te_fixed without it.
  !=========================================================================
  subroutine gas_state_setup(nH, nleaf, xHI)
    implicit none
    real(kind=wp), intent(in)           :: nH(:)
    integer,       intent(in)           :: nleaf
    real(kind=wp), intent(in), optional :: xHI(:)
    integer :: il

    gas_nleaf = nleaf
    call allocate_shared(state_owner(1), gas_nH, [nleaf])
    call allocate_shared(state_owner(2), gas_xHI, [nleaf])
    call allocate_shared(state_owner(3), gas_xHII, [nleaf])
    call allocate_shared(state_owner(4), gas_xHeI, [nleaf])
    call allocate_shared(state_owner(5), gas_xHeII, [nleaf])
    call allocate_shared(state_owner(6), gas_xHeIII, [nleaf])
    call allocate_shared(state_owner(7), gas_ne, [nleaf])
    call allocate_shared(state_owner(8), gas_Te, [nleaf])

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
       if (helium_starts_as_hydrogen) then
          do il = 1, nleaf
             gas_xHeI(il)   = gas_xHI(il)
             gas_xHeII(il)  = gas_xHII(il)
             gas_xHeIII(il) = 0.0_wp
          end do
       else
          gas_xHeI(:)  = par%xHeI_init
          gas_xHeII(:) = par%xHeII_init
          gas_xHeIII(:) = max(0.0_wp, 1.0_wp - par%xHeI_init - par%xHeII_init)
       end if
       do il = 1, nleaf
          gas_ne(il) = gas_nH(il) * ( gas_xHII(il) &
                       + par%abund_He*(gas_xHeII(il) + 2.0_wp*gas_xHeIII(il)) )
       end do
       if (par%solve_te) then
          gas_Te(:) = par%te_start
       else
          gas_Te(:) = par%te_fixed
       end if
    end if
    call gas_state_sync()

    if (mpar%p_rank == 0) then
       write(*,'(a,i12)')    ' GAS: leaf state allocated, nleaf = ', nleaf
       if (present(xHI)) then
          write(*,'(a)')     ' GAS: x_HI initialized from the AMR file column'
       else
          write(*,'(a,es13.6)') ' GAS: x_HI initialized to ', par%xHI_init
       end if
       if (helium_starts_as_hydrogen .and. present(xHI)) then
          write(*,'(a)')     ' GAS: x_HeI = x_HI, x_HeII = 1 - x_HI, '// &
                             'x_HeIII = 0 in each leaf (helium starts as hydrogen)'
       else if (helium_starts_as_hydrogen) then
          write(*,'(a,2es13.6,a)') ' GAS: x_HeI, x_HeII initialized to ', &
                             par%xHI_init, 1.0_wp - par%xHI_init, &
                             ' (helium starts as hydrogen)'
       else
          write(*,'(a,2es13.6)') ' GAS: x_HeI, x_HeII initialized to ', &
                             par%xHeI_init, par%xHeII_init
       end if
       if (par%solve_te) then
          write(*,'(a,f10.1,a)') ' GAS: T_e initialized to ', par%te_start, &
             ' K (te_start, the start of the thermal solve)'
       else
          write(*,'(a,f10.1,a)') ' GAS: T_e held at ', par%te_fixed, &
             ' K in every leaf (te_fixed; solve_te = F)'
       end if
       write(*,'(a,f8.4)')   ' GAS: He abundance n_He/n_H = ', par%abund_He
    end if
  end subroutine gas_state_setup

  !=========================================================================
  ! Recreate the state on a re-refined tree: fresh shared windows
  ! (the caller freed the old ones, amr_refine_mod step 4b), filled from position-mapped
  ! local arrays.  Collective; h_rank 0 writes.
  !=========================================================================
  subroutine gas_state_recreate(nleaf, nHv, x1, x1p, x2, x3, x4, nev, tev)
    implicit none
    integer,       intent(in) :: nleaf
    real(kind=wp), intent(in) :: nHv(nleaf), x1(nleaf), x1p(nleaf), &
                                 x2(nleaf), x3(nleaf), x4(nleaf), &
                                 nev(nleaf), tev(nleaf)
    gas_nleaf = nleaf
    call allocate_shared(state_owner(1), gas_nH, [nleaf])
    call allocate_shared(state_owner(2), gas_xHI, [nleaf])
    call allocate_shared(state_owner(3), gas_xHII, [nleaf])
    call allocate_shared(state_owner(4), gas_xHeI, [nleaf])
    call allocate_shared(state_owner(5), gas_xHeII, [nleaf])
    call allocate_shared(state_owner(6), gas_xHeIII, [nleaf])
    call allocate_shared(state_owner(7), gas_ne, [nleaf])
    call allocate_shared(state_owner(8), gas_Te, [nleaf])
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
  ! Finish previous buffer readers, then clear the whole host contribution.
  ! World-rank ownership is unchanged: each writer stores only its own leaves.
  ! Foreign-host leaves stay zero until sum_shared_hosts exchanges the sums.
  ! Current state is read-only throughout this reset and the following solve.
  !=========================================================================
  subroutine gas_state_new_reset(solve_temperature)
    logical, intent(in) :: solve_temperature
    if (.not. shared_is_live(new_state_owner(1))) then
       call allocate_shared(new_state_owner(1), xHI_new, [gas_nleaf])
       call allocate_shared(new_state_owner(2), xHII_new, [gas_nleaf])
       call allocate_shared(new_state_owner(3), xHeI_new, [gas_nleaf])
       call allocate_shared(new_state_owner(4), xHeII_new, [gas_nleaf])
       call allocate_shared(new_state_owner(5), xHeIII_new, [gas_nleaf])
       call allocate_shared(new_state_owner(6), ne_new, [gas_nleaf])
    end if
    if (solve_temperature) then
       if (.not. shared_is_live(new_state_owner(7))) &
          call allocate_shared(new_state_owner(7), te_new, [gas_nleaf])
    else
       if (shared_is_live(new_state_owner(7))) call free_shared(new_state_owner(7), te_new)
    end if
    call publish_shared(new_state_owner)
    if (mpar%h_rank == 0) then
       xHI_new = 0.0_wp
       xHII_new = 0.0_wp
       xHeI_new = 0.0_wp
       xHeII_new = 0.0_wp
       xHeIII_new = 0.0_wp
       ne_new = 0.0_wp
       if (associated(te_new)) te_new = 0.0_wp
    end if
    call publish_shared(new_state_owner)
    new_state_ready = .true.
  end subroutine gas_state_new_reset

  !=========================================================================
  ! Called only after solves and all old-state change measures have finished.
  ! Publish disjoint writers, combine the zero-padded host contributions,
  ! then copy into stable current-state allocations and publish that update.
  ! No temperature buffer in a fixed-temperature solve: T_e stays untouched.
  !=========================================================================
  subroutine gas_state_publish()
    integer :: il
    if (.not. new_state_ready) call abort_run('gas_state_publish: new state was not reset')
    call sum_shared_hosts(new_state_owner)
    call gas_state_sync()
    if (mpar%h_rank == 0) then
       ! Scalar stores need no array overlap temporary.
       do il = 1, gas_nleaf
          gas_xHI(il) = xHI_new(il)
          gas_xHII(il) = xHII_new(il)
          gas_xHeI(il) = xHeI_new(il)
          gas_xHeII(il) = xHeII_new(il)
          gas_xHeIII(il) = xHeIII_new(il)
          gas_ne(il) = ne_new(il)
          if (associated(te_new)) gas_Te(il) = te_new(il)
       end do
    end if
    call gas_state_sync()
    new_state_ready = .false.
  end subroutine gas_state_publish

  subroutine gas_state_new_release()
    if (shared_is_live(new_state_owner(1))) call free_shared(new_state_owner(1), xHI_new)
    if (shared_is_live(new_state_owner(2))) call free_shared(new_state_owner(2), xHII_new)
    if (shared_is_live(new_state_owner(3))) call free_shared(new_state_owner(3), xHeI_new)
    if (shared_is_live(new_state_owner(4))) call free_shared(new_state_owner(4), xHeII_new)
    if (shared_is_live(new_state_owner(5))) call free_shared(new_state_owner(5), xHeIII_new)
    if (shared_is_live(new_state_owner(6))) call free_shared(new_state_owner(6), ne_new)
    if (shared_is_live(new_state_owner(7))) call free_shared(new_state_owner(7), te_new)
    new_state_ready = .false.
  end subroutine gas_state_new_release

  !=========================================================================
  ! Publish one complete gas-state update and order its readers.
  !=========================================================================
  subroutine gas_state_sync()
    call publish_shared(state_owner)
  end subroutine gas_state_sync

  subroutine gas_state_release()
    ! Complete old-state reads before releasing the canonical pointers.
    call gas_state_sync()
    call gas_state_new_release()
    if (shared_is_live(state_owner(1))) call free_shared(state_owner(1), gas_nH)
    if (shared_is_live(state_owner(2))) call free_shared(state_owner(2), gas_xHI)
    if (shared_is_live(state_owner(3))) call free_shared(state_owner(3), gas_xHII)
    if (shared_is_live(state_owner(4))) call free_shared(state_owner(4), gas_xHeI)
    if (shared_is_live(state_owner(5))) call free_shared(state_owner(5), gas_xHeII)
    if (shared_is_live(state_owner(6))) call free_shared(state_owner(6), gas_xHeIII)
    if (shared_is_live(state_owner(7))) call free_shared(state_owner(7), gas_ne)
    if (shared_is_live(state_owner(8))) call free_shared(state_owner(8), gas_Te)
    gas_nleaf = 0
  end subroutine gas_state_release

end module gas_state_mod
