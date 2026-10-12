program check_neutral_hydrogen_excitation
! Atomic-hydrogen impact on the ground-term fine structure (C I, C II, O I;
! par%neutral_hydrogen_excitation): species_mod::metal_cooling_H must be,
! for every stage whose atom carries hydrogen rates, n_ion times the
! radiated power with both colliders in one statistical equilibrium minus
! the power with electrons alone (nlevel_emissivities), and zero with the
! switch off.  Neutral gas (n_HII = 0), no photoionization.
  use define
  use species_mod
  use nlevel_mod
  use setup_mod, only: read_parameters_with_implicit_defaults
  use mpi
  implicit none
  type(nlevel_atom_type) :: atom
  character(len=PATH_LEN) :: directory, option_file, filename
  character(len=8) :: element
  real(kind=wp), parameter :: temperatures(3)=[100.0_wp,300.0_wp,1000.0_wp]
  real(kind=wp), parameter :: densities(3)=[1e-4_wp,1.0_wp,1e4_wp]
  real(kind=wp), parameter :: nH=1000.0_wp
  real(kind=wp) :: frac(MAX_ION_STAGES), cooling, increment, expected, T, ne
  real(kind=wp), allocatable :: gamma(:,:), emission(:), electrons(:)
  integer :: ierr, i, j, stage, nhydrogen
  logical :: ok
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD
  call get_command_argument(1,directory)
  call get_command_argument(2,element)
  call get_command_argument(3,option_file)
  ! Compiled default: on.
  if (.not. par%neutral_hydrogen_excitation) stop 1
  call read_parameters_with_implicit_defaults(trim(option_file))
  par%abund_C=0; par%abund_N=0; par%abund_O=0; par%abund_Ne=0
  par%abund_S=0; par%abund_Ar=0; par%abund_Mg=0; par%abund_Fe=0
  par%abund_Si=0; par%abund_Cl=0; par%abund_Ca=0
  if (element=='c') par%abund_C=1e-4_wp
  if (element=='o') par%abund_O=1e-4_wp
  par%atomic_dir=trim(directory)
  par%cooling_model='low_density'
  call species_setup(1)
  allocate(gamma(elem_nstage(1)-1,1))
  gamma=0
  call species_set_photoionization_rates(1,gamma)
  do i=1,3
    T=temperatures(i)
    do j=1,3
      ne=densities(j)
      call species_fractions(1,1,T,ne,nH,0.0_wp,frac)
      cooling=metal_cooling(1,T,nH,ne,nH,0.0_wp)
      increment=metal_cooling_H(1,T,nH,ne,nH,0.0_wp)
      expected=0
      nhydrogen=0
      do stage=1,elem_nstage(1)
        filename=metal_excitation_data_file('nlevel',trim(element),stage)
        call nlevel_load(trim(filename),atom,ok)
        if (.not. ok) cycle
        if (.not. nlevel_has_hcoll(atom)) cycle
        nhydrogen=nhydrogen+1
        if (allocated(emission)) deallocate(emission,electrons)
        allocate(emission(atom%nrad),electrons(atom%nrad))
        call nlevel_emissivities(atom,T,ne,electrons)
        call nlevel_emissivities(atom,T,ne,emission,nH)
        if (par%neutral_hydrogen_excitation) &
           expected=expected+elem_abund(1)*nH*frac(stage)*sum(emission-electrons)
      end do
      if (nhydrogen==0) stop 2
      if (abs(increment-expected)>1e-10_wp*max(abs(expected),1e-300_wp)) stop 3
      if (.not. par%neutral_hydrogen_excitation .and. increment/=0) stop 4
      write(*,'(a,i2,4(1x,es24.16))') 'COOLING',nhydrogen,T,ne,cooling,increment
    end do
  end do
  call MPI_FINALIZE(ierr)
end program check_neutral_hydrogen_excitation
