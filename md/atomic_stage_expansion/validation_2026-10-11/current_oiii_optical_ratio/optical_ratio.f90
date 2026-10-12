program optical_ratio
  use mpi
  use define, only: wp
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_emissivities, nlevel_line_ident
  implicit none
  type(nlevel_atom_type) :: atom
  real(kind=wp), allocatable :: emissivity(:)
  real(kind=wp), parameter :: temperatures(3)=[1000._wp,10000._wp,100000._wp]
  real(kind=wp), parameter :: densities(4)=[1e-6_wp,100._wp,1e6_wp,1e10_wp]
  character(len=1024) :: filename, output
  integer :: ierr, i, j, k, short_line, long_line, unit
  logical :: ok
  call MPI_INIT(ierr)
  call get_command_argument(1,filename)
  call get_command_argument(2,output)
  call nlevel_load(trim(filename),atom,ok)
  if (.not.ok) error stop 'O III atom did not load'
  short_line=0;long_line=0
  do k=1,atom%nrad
    if (atom%rl(k)==2 .and. atom%ru(k)==4) short_line=k
    if (atom%rl(k)==3 .and. atom%ru(k)==4) long_line=k
  end do
  if (short_line==0 .or. long_line==0) error stop 'Optical O III transitions missing'
  allocate(emissivity(atom%nrad))
  open(newunit=unit,file=trim(output),status='new')
  write(unit,'(a,i0)') '# levels ',atom%nlev
  do i=1,3
    do j=1,4
      call nlevel_emissivities(atom,temperatures(i),densities(j),emissivity)
      if (min(emissivity(short_line),emissivity(long_line))<=0._wp) error stop 'Nonpositive line energy'
      write(unit,'(6(es26.16e3,1x))') temperatures(i),densities(j), &
        nlevel_line_ident(atom,short_line),nlevel_line_ident(atom,long_line), &
        emissivity(short_line),emissivity(long_line)
    end do
  end do
  close(unit)
  call MPI_FINALIZE(ierr)
end program optical_ratio
