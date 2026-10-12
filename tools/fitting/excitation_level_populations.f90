! Build outside the tree from copies of src/define.f90 and src/nlevel_mod.f90:
!   ifx -O2 -cpp -c define.f90 nlevel_mod.f90
!   ifx -O2 -cpp -o excitation_level_populations.x excitation_level_populations.f90 \
!       define.o nlevel_mod.o -qmkl=sequential
! Run: excitation_level_populations.x nlevel_<el>_<stage>.txt T_K n_e_cm3
! Writes to standard output the line cooling per (n_e n_ion) [erg cm^3 s^-1]
! and the level populations (normalized to the ion) of the loaded atom.
program excitation_level_populations
! Statistical-equilibrium populations of one atom at one (T, n_e).
  use define, only: wp
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_cooling_prepT, &
       nlevel_cooling_solve, nlevel_populations
  implicit none
  type(nlevel_atom_type) :: atom
  real(kind=wp), allocatable :: rc(:,:), rr(:,:), pop(:)
  real(kind=wp) :: t, ne
  character(len=1024) :: atom_file, argument
  integer :: i
  logical :: ok
  call get_command_argument(1, atom_file)
  call get_command_argument(2, argument); read(argument,*) t
  call get_command_argument(3, argument); read(argument,*) ne
  call nlevel_load(trim(atom_file), atom, ok)
  if (.not. ok) error stop 'Atom did not load'
  allocate(rc(atom%nlev,atom%nlev), rr(atom%nlev,atom%nlev), pop(atom%nlev))
  call nlevel_cooling_prepT(atom, t, rc, rr)
  call nlevel_populations(atom, rc, rr, ne, pop)
  write(*,'(a,3i9)') '# nlev nrad nups', atom%nlev, atom%nrad, atom%nups
  write(*,'(a,es24.16e3)') '# cooling ', nlevel_cooling_solve(atom, rc, rr, ne)
  do i = 1, atom%nlev
     write(*,'(i5,1x,es24.16e3)') i, pop(i)
  end do
end program excitation_level_populations
