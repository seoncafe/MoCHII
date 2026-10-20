program check_neutral_metal_initialization
  use define
  use mpi
  use gas_state_mod
  use species_mod, only: species_setup, species_stage_cache_refresh, &
       species_set_photoionization_rates, species_fractions, species_stage_fraction, &
       species_ne_prepare, species_metal_sums
  implicit none
  character(len=PATH_LEN) :: directory
  real(kind=wp) :: gamma(6,3), fraction(MAX_ION_STAGES), expected(7)
  real(kind=wp) :: electrons, cx_in, cx_out
  integer :: ierr, i
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpar%nproc,ierr)
  mpar%hostcomm = MPI_COMM_WORLD
  mpar%h_rank = mpar%p_rank
  call get_command_argument(1,directory)
  par%atomic_dir = trim(directory)
  par%cooling_model = 'low_density'
  par%abund_C = 1.0e-4_wp
  par%metal_collisional_ionization = .false.
  par%xHI_init = 1.0_wp
  par%xHeI_init = 1.0_wp
  par%xHeII_init = 0.0_wp
  par%te_fixed = 1.0e4_wp
  call gas_state_setup([100.0_wp,10.0_wp,0.0_wp],3)
  call species_setup(3)
  call species_stage_cache_refresh()
  do i=1,3
     if (species_stage_fraction(1,1,i) /= 1.0_wp) error stop 'Neutral initial opacity missing'
     if (sum([(species_stage_fraction(1,ierr,i),ierr=2,7)]) /= 0.0_wp) error stop 'Charged initial metals'
  end do
  call species_fractions(1,1,1.0e4_wp,0.0_wp,100.0_wp,0.0_wp,fraction)
  if (fraction(1) /= 1.0_wp .or. sum(fraction(2:)) /= 0.0_wp) error stop 'Direct neutral limit'
  call species_ne_prepare(1,1.0e4_wp)
  call species_metal_sums(100.0_wp,0.0_wp,100.0_wp,0.0_wp,electrons,cx_in,cx_out)
  if (electrons /= 0.0_wp .or. cx_in /= 0.0_wp .or. cx_out /= 0.0_wp) error stop 'Neutral cached balance'

  ! The synthetic atom has alpha=2.5e-12 at 1e4 K and no charge exchange.
  ! At ne=100 and Gamma=2e-10 the analytic stage ratio is 0.8.
  gamma = 2.0e-10_wp
  call species_set_photoionization_rates(1,gamma)
  if (mpar%h_rank == 0) gas_ne(1:2) = 100.0_wp
  call gas_state_sync()
  call species_stage_cache_refresh()
  do i=1,7
     expected(i) = 0.8_wp**(i-1)
  end do
  expected = expected/sum(expected)
  if (maxval(abs([(species_stage_fraction(1,i,1),i=1,7)]-expected)) > 1e-13_wp) &
       error stop 'Illuminated opacity'
  if (species_stage_fraction(1,1,3) /= 1.0_wp) error stop 'Vacuum convention'
  call species_ne_prepare(1,1.0e4_wp)
  call species_metal_sums(100.0_wp,100.0_wp,100.0_wp,0.0_wp,electrons,cx_in,cx_out)
  if (abs(electrons-.01_wp*sum(expected*[(real(i-1,wp),i=1,7)])) > 1e-14_wp) &
       error stop 'Illuminated cached electron balance'
  if (mpar%p_rank == 0) write(*,'(a)') 'PASS: neutral initialization, illuminated equilibrium, and vacuum'
  call MPI_FINALIZE(ierr)
end program check_neutral_metal_initialization
