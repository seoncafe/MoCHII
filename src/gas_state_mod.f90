module gas_state_mod
!---------------------------------------------------------------------------
! MoCHII: persistent gas leaf state.
!
! nH/xHI are read from the AMR file and used for the dust-density step;
! MoCHII persists them as shared-memory leaf arrays and adds the He
! ionization fractions and T_e/n_e.
!
! The state is INITIALIZED (from the file columns where present,
! else from par%xHI_init etc.); when gas_niter = 0 it is never updated —
! no equilibrium solve, no opacity feedback.  gas_state_setup is called
! from grid_create_amr before the transient read arrays are deallocated.
!---------------------------------------------------------------------------
  use define
  use memory_mod, only : create_shared_mem
  implicit none
  private

  public :: gas_state_setup, gas_state_recreate
  public :: gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, gas_Te, gas_nleaf

  real(kind=wp), pointer :: gas_nH(:)    => null()  ! H number density [cm^-3]
  real(kind=wp), pointer :: gas_xHI(:)   => null()  ! n_HI / n_H
  real(kind=wp), pointer :: gas_xHeI(:)  => null()  ! n_HeI / n_He
  real(kind=wp), pointer :: gas_xHeII(:) => null()  ! n_HeII / n_He
  real(kind=wp), pointer :: gas_ne(:)    => null()  ! electron density [cm^-3]
  real(kind=wp), pointer :: gas_Te(:)    => null()  ! electron temperature [K]
  integer :: gas_nleaf = 0

contains

  !=========================================================================
  ! nH is mandatory (from the AMR file); xHI is optional (file column) —
  ! when absent every leaf starts at par%xHI_init.  He fractions start at
  ! par%xHeI_init / par%xHeII_init (no He columns in the generic format).
  !=========================================================================
  subroutine gas_state_setup(nH, nleaf, xHI)
    use mpi
    implicit none
    real(kind=wp), intent(in)           :: nH(:)
    integer,       intent(in)           :: nleaf
    real(kind=wp), intent(in), optional :: xHI(:)
    integer :: il, ierr

    gas_nleaf = nleaf
    call create_shared_mem(gas_nH,    [nleaf])
    call create_shared_mem(gas_xHI,   [nleaf])
    call create_shared_mem(gas_xHeI,  [nleaf])
    call create_shared_mem(gas_xHeII, [nleaf])
    call create_shared_mem(gas_ne,    [nleaf])
    call create_shared_mem(gas_Te,    [nleaf])

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
       gas_xHeI(:)  = par%xHeI_init
       gas_xHeII(:) = par%xHeII_init
       do il = 1, nleaf
          gas_ne(il) = gas_nH(il) * ( (1.0_wp - gas_xHI(il)) &
                       + par%He_abund*(gas_xHeII(il) &
                       + 2.0_wp*max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il))) )
       end do
       gas_Te(:) = par%te_fixed
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)

    if (mpar%p_rank == 0) then
       write(*,'(a,i12)')    ' GAS: leaf state allocated, nleaf = ', nleaf
       if (present(xHI)) then
          write(*,'(a)')     ' GAS: x_HI initialized from the AMR file column'
       else
          write(*,'(a,f8.4)')' GAS: x_HI initialized to ', par%xHI_init
       end if
       write(*,'(a,2f8.4)')  ' GAS: x_HeI, x_HeII initialized to ', &
                             par%xHeI_init, par%xHeII_init
       write(*,'(a,f8.4)')   ' GAS: He abundance n_He/n_H = ', par%He_abund
    end if
  end subroutine gas_state_setup

  !=========================================================================
  ! Recreate the state on a re-refined tree: fresh shared windows
  ! (the old ones are leaked until finalize), filled from position-mapped
  ! local arrays.  Collective; h_rank 0 writes.
  !=========================================================================
  subroutine gas_state_recreate(nleaf, nHv, x1, x2, x3, nev, tev)
    use mpi
    implicit none
    integer,       intent(in) :: nleaf
    real(kind=wp), intent(in) :: nHv(nleaf), x1(nleaf), x2(nleaf), &
                                 x3(nleaf), nev(nleaf), tev(nleaf)
    integer :: ierr

    gas_nleaf = nleaf
    call create_shared_mem(gas_nH,    [nleaf])
    call create_shared_mem(gas_xHI,   [nleaf])
    call create_shared_mem(gas_xHeI,  [nleaf])
    call create_shared_mem(gas_xHeII, [nleaf])
    call create_shared_mem(gas_ne,    [nleaf])
    call create_shared_mem(gas_Te,    [nleaf])
    if (mpar%h_rank == 0) then
       gas_nH(1:nleaf)    = nHv
       gas_xHI(1:nleaf)   = x1
       gas_xHeI(1:nleaf)  = x2
       gas_xHeII(1:nleaf) = x3
       gas_ne(1:nleaf)    = nev
       gas_Te(1:nleaf)    = tev
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine gas_state_recreate

end module gas_state_mod
