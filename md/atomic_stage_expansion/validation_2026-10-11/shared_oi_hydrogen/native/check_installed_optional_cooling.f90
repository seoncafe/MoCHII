program check_installed_optional_cooling
  use define
  use species_mod
  use nlevel_mod
  use mpi
  implicit none
  type(nlevel_atom_type) :: atom
  character(len=PATH_LEN) :: directory, option_file
  character(len=8) :: element
  real(kind=wp), parameter :: temperatures(3)=[1000.0_wp,3000.0_wp,10000.0_wp]
  real(kind=wp), parameter :: densities(3)=[1e-4_wp,1.0_wp,1e4_wp]
  real(kind=wp) :: frac(7), gamma(2,1), cooling, increment, expected, T, ne
  real(kind=wp), allocatable :: emission(:), electrons(:)
  integer :: ierr, i, j
  logical :: ok
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD
  call get_command_argument(1,directory)
  call get_command_argument(2,element)
  call get_command_argument(3,option_file)
  if (par%additional_line_cooling .or. par%neutral_hydrogen_excitation) stop 1
  block
    use setup_mod, only: read_parameters_with_implicit_defaults
    call read_parameters_with_implicit_defaults(trim(option_file))
  end block
  par%abund_C=0; par%abund_N=0; par%abund_O=0; par%abund_Ne=0
  par%abund_S=0; par%abund_Ar=0; par%abund_Mg=0; par%abund_Fe=0
  par%abund_Si=0; par%abund_Cl=0; par%abund_Ca=0
  if (element=='ne') par%abund_Ne=1e-4_wp
  if (element=='o') par%abund_O=1e-4_wp
  par%atomic_dir=trim(directory)
  par%cooling_model='low_density'
  call species_setup(1)
  gamma=0
  call species_set_photoionization_rates(1,gamma)
  call nlevel_load(trim(directory)//'/nlevel_'//trim(element)//'_1.txt',atom,ok)
  if (.not. ok) stop 2
  allocate(emission(atom%nrad),electrons(atom%nrad))
  do i=1,3
    T=temperatures(i)
    if (element=='ne') T=100*T
    do j=1,3
      ne=densities(j)
      call species_fractions(1,1,T,ne,1000.0_wp,0.0_wp,frac)
      if (abs(frac(1)-1)>1e-12_wp) stop 3
      cooling=metal_cooling(1,T,1000.0_wp,ne,1000.0_wp,0.0_wp)
      increment=metal_cooling_H(1,T,1000.0_wp,ne,1000.0_wp,0.0_wp)
      call nlevel_emissivities(atom,T,ne,electrons)
      emission=electrons
      if (par%neutral_hydrogen_excitation) call nlevel_emissivities(atom,T,ne,emission,1000.0_wp)
      expected=0.1_wp*sum(emission-electrons)
      if (abs(increment-expected)>1e-10_wp*max(abs(expected),1e-35_wp)) stop 4
      if (element=='ne') then
        if (.not. par%additional_line_cooling .and. cooling/=0) stop 5
        if (par%additional_line_cooling .and. cooling<=0) stop 6
      end if
      if (.not. par%neutral_hydrogen_excitation .and. increment/=0) stop 7
      write(*,'(a,5(1x,es24.16))') 'COOLING',T,ne,cooling,increment,sum(emission)
    end do
  end do
  call MPI_FINALIZE(ierr)
end program check_installed_optional_cooling
