! Build outside the tree from copies of src/define.f90 and src/nlevel_mod.f90:
!   ifx -O2 -cpp -c define.f90 nlevel_mod.f90
!   ifx -O2 -cpp -o excitation_atom_cooling.x excitation_atom_cooling.f90 define.o nlevel_mod.o -qmkl
! Run: excitation_atom_cooling.x nlevel_<el>_<stage>.txt OUT; writes OUT (cooling at 81
! temperatures, 1e3-1e5 K) and OUT.collision_balance (radiation-off populations).
! Used by validate_chianti_excitation_model.py (--fortran-cooling OUT).
program excitation_atom_cooling
! Load one atom with nlevel_load and write the line cooling per (n_e n_ion)
! at 81 temperatures over 1e3-1e5 K for n_e = 1e-6, 1e2, 1e4, 1e6 cm^-3.
  use define, only: wp
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_cooling_prepT, &
       nlevel_cooling_solve, nlevel_populations
  implicit none
  type(nlevel_atom_type) :: atom
  real(kind=wp), allocatable :: rc(:,:), rr(:,:), pop(:), zero(:,:)
  real(kind=wp), parameter :: balance_t(3)=[3e3_wp,1e4_wp,1e5_wp]
  real(kind=wp), parameter :: density(4)=[1e-6_wp,1e2_wp,1e4_wp,1e6_wp]
  real(kind=wp) :: t
  character(len=1024) :: atom_file, output
  integer :: i, j, n, unit
  logical :: ok
  call get_command_argument(1,atom_file)
  call get_command_argument(2,output)
  call nlevel_load(trim(atom_file),atom,ok)
  if (.not.ok) error stop 'Atom did not load'
  allocate(rc(atom%nlev,atom%nlev),rr(atom%nlev,atom%nlev))
  n=81
  open(newunit=unit,file=trim(output),status='replace')
  write(unit,'(a,3i9)') '# nlev nrad nups', atom%nlev, atom%nrad, atom%nups
  do i=1,n
     t=10.0_wp**(3.0_wp+2.0_wp*real(i-1,wp)/real(n-1,wp))
     call nlevel_cooling_prepT(atom,t,rc,rr)
     write(unit,'(*(es26.16e3,1x))') t,(nlevel_cooling_solve(atom,rc,rr,density(j)),j=1,4)
  end do
  close(unit)
  ! Collision-only balance (radiation off): must be the Boltzmann distribution.
  allocate(pop(atom%nlev),zero(atom%nlev,atom%nlev))
  zero=0.0_wp
  open(newunit=unit,file=trim(output)//'.collision_balance',status='replace')
  do i=1,3
     call nlevel_cooling_prepT(atom,balance_t(i),rc,rr)
     call nlevel_populations(atom,rc,zero,1.0_wp,pop)
     write(unit,'(*(es26.16e3,1x))') balance_t(i),pop
  end do
  close(unit)
end program excitation_atom_cooling
