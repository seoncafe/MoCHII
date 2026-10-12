program check_charge_transfer_ranges
  use define
  use species_mod
  use mpi
  implicit none
  integer :: ierr, ie, it, ios
  real(kind=wp) :: temperature, ionizing, recombining
  character(len=PATH_LEN) :: directory, input, output
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD; mpar%h_rank=0
  call get_command_argument(1,directory)
  call get_command_argument(2,input)
  call get_command_argument(3,output)
  par%atomic_dir=trim(directory); par%cooling_model='low_density'
  par%abund_C=1e-4_wp; par%abund_N=1e-4_wp; par%abund_O=1e-4_wp
  par%abund_Ne=1e-4_wp; par%abund_S=1e-4_wp; par%abund_Ar=1e-4_wp
  par%abund_Mg=1e-4_wp; par%abund_Fe=1e-4_wp; par%abund_Si=1e-4_wp
  par%abund_Cl=1e-4_wp; par%abund_Ca=1e-4_wp
  call species_setup(1)
  open(unit=71,file=trim(input),status='old')
  open(unit=72,file=trim(output),status='replace')
  do
     read(71,*,iostat=ios) ie,it,temperature
     if (ios /= 0) exit
     call species_cx_coefficients(ie,it,temperature,ionizing,recombining)
     write(72,'(2i5,3es27.17e3)') ie,it,temperature,ionizing,recombining
  end do
  close(71); close(72)
  call MPI_FINALIZE(ierr)
end program check_charge_transfer_ranges
