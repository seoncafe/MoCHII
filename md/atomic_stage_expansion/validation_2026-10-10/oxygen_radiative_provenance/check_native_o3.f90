program check_native_o3
  use mpi
  use define, only: wp
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_emissivities, nlevel_cooling_at
  implicit none
  type(nlevel_atom_type) :: atom
  real(kind=wp), allocatable :: emis(:)
  real(kind=wp), parameter :: density(4)=[1e-6_wp,100.0_wp,1e6_wp,1e10_wp]
  real(kind=wp) :: t, lam, ratio
  character(len=1024) :: fname, output
  integer :: ierr, i, j, k, line4959, line5007, unit
  logical :: ok
  call MPI_INIT(ierr)
  call get_command_argument(1,fname)
  call get_command_argument(2,output)
  call nlevel_load(trim(fname),atom,ok)
  if (.not.ok) error stop 'O III atom did not load'
  allocate(emis(atom%nrad))
  line4959=0; line5007=0
  do k=1,atom%nrad
     if(atom%rl(k)==2.and.atom%ru(k)==4) line4959=k
     if(atom%rl(k)==3.and.atom%ru(k)==4) line5007=k
  end do
  if(line4959==0.or.line5007==0) error stop 'Missing optical transitions'
  open(newunit=unit,file=trim(output),status='new')
  do i=0,800
     t=10.0_wp**(3.0_wp+real(i,wp)/400.0_wp)
     do j=1,4
        call nlevel_emissivities(atom,t,density(j),emis)
        lam=nlevel_cooling_at(atom,t,density(j))
        if(abs(sum(emis)/density(j)-lam)>1e-12_wp*abs(lam)) &
             error stop 'Line emission and cooling differ'
        ratio=emis(line5007)/emis(line4959)
        write(unit,'(*(es26.16e3,1x))') t,density(j),lam,ratio
     end do
  end do
  close(unit)
  call MPI_FINALIZE(ierr)
end program check_native_o3
