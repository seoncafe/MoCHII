program check_atomic_cooling_without_fit
! A stage whose excitation model has a level atom but no fitted cooling file
! cools by the absolute n-level table (species_setup, nlevel_cooling_add).
! In a model directory that holds level atoms and no cooling file at all,
! species_mod::metal_cooling of a one-element gas, divided by n_e n_X, must
! equal sum_i f(i) Lambda_i over the stages with an atom, Lambda_i the
! direct statistical-equilibrium cooling of that atom (nlevel_cooling_at) at
! n_e for 'local_ne' and at the 1e-6 cm^-3 low-density reference for
! 'low_density', to the table's interpolation accuracy.
  use define
  use mpi
  use setup_mod, only: read_parameters_with_implicit_defaults
  use species_mod
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_cooling_at, metal_excitation_data_file
  implicit none
  character(len=PATH_LEN) :: directory,input,output,filename
  character(len=8) :: element
  real(kind=wp),parameter :: densities(4)=[1e-6_wp,100.0_wp,1e6_wp,1e10_wp]
  real(kind=wp),parameter :: nH=100.0_wp, abundance=1e-4_wp
  type(nlevel_atom_type) :: atom(MAX_ION_STAGES)
  logical :: has_atom(MAX_ION_STAGES)
  real(kind=wp) :: T,ne,cool,frac(MAX_ION_STAGES),coefficient,direct,line_ne
  real(kind=wp),allocatable :: gamma(:,:)
  integer :: ierr,i,j,k,unit
  logical :: exists
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD
  call get_command_argument(1,directory)
  call get_command_argument(2,element)
  call get_command_argument(3,input)
  call get_command_argument(4,output)
  call read_parameters_with_implicit_defaults(trim(input))
  par%atomic_dir=trim(directory)
  par%abund_C=0;par%abund_N=0;par%abund_O=0;par%abund_Ne=0
  par%abund_S=0;par%abund_Ar=0;par%abund_Mg=0;par%abund_Fe=0
  par%abund_Si=0;par%abund_Cl=0;par%abund_Ca=0
  select case(trim(element))
  case('c');par%abund_C=abundance
  case('n');par%abund_N=abundance
  case('o');par%abund_O=abundance
  case('ne');par%abund_Ne=abundance
  case('s');par%abund_S=abundance
  case('ar');par%abund_Ar=abundance
  case('mg');par%abund_Mg=abundance
  case('fe');par%abund_Fe=abundance
  case('si');par%abund_Si=abundance
  case('cl');par%abund_Cl=abundance
  case('ca');par%abund_Ca=abundance
  case default;error stop 'unknown element'
  end select
  ! Atomic-hydrogen impact is the separate consumer metal_cooling_H.
  par%neutral_hydrogen_excitation=.false.
  call species_setup(1)
  if(n_elements/=1) stop 1
  has_atom=.false.
  do k=1,elem_nstage(1)
     filename=metal_excitation_data_file('cooling',trim(element),k)
     inquire(file=trim(filename),exist=exists)
     if(exists) error stop 'the fixture has a fitted cooling file'
     filename=metal_excitation_data_file('nlevel',trim(element),k)
     call nlevel_load(trim(filename),atom(k),has_atom(k))
  end do
  if(.not.any(has_atom)) error stop 'the fixture has no level atom'
  allocate(gamma(elem_nstage(1)-1,1))
  gamma=2e-8_wp
  call species_set_photoionization_rates(1,gamma)
  if(metal_cooling(1,10000.0_wp,nH,0.0_wp,50.0_wp,50.0_wp)/=0) stop 2
  if(metal_cooling(1,10000.0_wp,0.0_wp,1.0_wp,0.0_wp,0.0_wp)/=0) stop 3
  open(newunit=unit,file=trim(output),status='new')
  do i=0,32
    T=10.0_wp**(3.0_wp+real(i,wp)/16.0_wp)
    do j=1,4
      ne=densities(j)
      call species_fractions(1,1,T,ne,50.0_wp,50.0_wp,frac)
      cool=metal_cooling(1,T,nH,ne,50.0_wp,50.0_wp)
      if(cool<=0) stop 4
      coefficient=cool/(ne*abundance*nH)
      line_ne=ne
      if(trim(par%cooling_model)=='low_density') line_ne=densities(1)
      direct=0
      do k=1,elem_nstage(1)
         if(has_atom(k)) direct=direct+frac(k)*nlevel_cooling_at(atom(k),T,line_ne)
      end do
      write(unit,'(4(es26.16e3,1x))') T,ne,coefficient,direct
    enddo
  enddo
  close(unit)
  write(*,'(a)') 'PASS: absolute atomic cooling without a coefficient file'
  call MPI_FINALIZE(ierr)
end program check_atomic_cooling_without_fit
