program check_fe_capture_cooling
! Isolated ground-level RR+DR check, with CI and charge transfer disabled.
! This fixture measures a coefficient; it is not an Fe plasma model.
  use define
  use species_mod
  use mpi
  implicit none
  integer :: ierr, j, q
  real(kind=wp) :: t, alpha, nem, direct_ne, cxi, cxr
  real(kind=wp) :: frac(MAX_ION_STAGES), gamma(MAX_ION_STAGES-1,1)
  character(len=PATH_LEN) :: directory, output
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD; mpar%h_rank=0
  call get_command_argument(1,directory)
  call get_command_argument(2,output)
  par%atomic_dir=trim(directory); par%cooling_model='low_density'
  par%abund_C=0; par%abund_N=0; par%abund_O=0; par%abund_Ne=0
  par%abund_S=0; par%abund_Ar=0; par%abund_Mg=0; par%abund_Si=0
  par%abund_Cl=0; par%abund_Ca=0; par%abund_Fe=1e-4_wp
  call species_setup(1)
  if(n_elements/=1.or.elem_nstage(1)/=7) error stop 'Fe rate fixture roster'
  gamma=1e-8_wp
  call species_set_photoionization_rates(1,gamma)
  open(unit=71,file=trim(output),status='new')
  do j=0,300
     t=1000.0_wp*100.0_wp**(real(j,wp)/300.0_wp)
     call species_fractions(1,1,t,10.0_wp,0.0_wp,0.0_wp,frac)
     alpha=metal_ionization_cooling(1,t,100.0_wp,10.0_wp,0.0_wp,0.0_wp)
     if(alpha<=0) error stop 'positive captured-electron energy loss'
     write(71,'(9es27.17e3)') t,alpha,frac
  end do
  close(71)
  print *, 'PASS: 301 captured-electron energy loss values'
  call MPI_FINALIZE(ierr)
end program check_fe_capture_cooling
