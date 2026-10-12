program tabulate_freefree_gaunt
!---------------------------------------------------------------------------
! Writes free_free_mod::free_free_gaunt on a (log T, E) grid for Z = 1..6,
! the Gaunt factor of the code's free-free emission, for the free-free
! absorption measurement (measure_freefree_absorption.py).
!---------------------------------------------------------------------------
  use define
  use free_free_mod, only : free_free_gaunt
  use mpi
  implicit none
  integer, parameter :: NT = 161, NE = 241
  real(kind=wp) :: lt, T, E, u
  integer :: it, ie, iz, ierr, unit
  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%data_dir = 'data'
  open(newunit=unit, file='tests/freefree_absorption/freefree_gaunt_table.txt', status='replace')
  write(unit,'(a)') '# log10T  E[eV]  g_ff(Z=1..6)  (free_free_mod::free_free_gaunt)'
  do it = 1, NT
     lt = 3.0_wp + 2.0_wp*real(it-1,wp)/real(NT-1,wp)      ! 1e3 - 1e5 K
     T = 10.0_wp**lt
     do ie = 1, NE
        E = 13.0_wp*(160.0_wp/13.0_wp)**(real(ie-1,wp)/real(NE-1,wp))
        u = E*ev2erg/(kboltz_cgs*T)
        write(unit,'(f8.4,f10.4,6es14.6)') lt, E, &
           (free_free_gaunt(real(iz,wp), T, u), iz = 1, 6)
     end do
  end do
  close(unit)
  call MPI_FINALIZE(ierr)
end program tabulate_freefree_gaunt
