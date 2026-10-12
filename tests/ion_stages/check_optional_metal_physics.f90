program check_optional_metal_physics
! Electron-impact ionization of the metals (par%metal_collisional_ionization)
! on a synthetic seven-stage carbon (finite positive rates, not carbon data):
! the namelist switch, its effect on the stage fractions, the direct and
! cached electron closure, and the ionization-energy loss
!   n_e n_C sum_it [ f(it+1) beta_it + f(it) C_it E_th,it ],
! printed for an independent evaluation (check_optional_metal_physics.py);
! and the terminal-stage refusal: bare C VII has no edge, C VI does.
  use define
  use species_mod
  use setup_mod, only: read_parameters_with_implicit_defaults
  use mpi
  implicit none
  character(len=PATH_LEN) :: directory, input_file
  character(len=8) :: mode
  real(kind=wp) :: gamma(6,1), frac(7), previous(7), cool, nem, expected
  real(kind=wp), parameter :: eedge(3)=[13.598_wp, 500.0_wp, 1000.0_wp], lum(2)=[1.0_wp, 1.0_wp]
  integer :: ierr, j, selection
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  mpar%hostcomm = MPI_COMM_WORLD
  call get_command_argument(1, directory)
  call get_command_argument(2, input_file)
  call get_command_argument(3, mode)
  if (.not. par%metal_collisional_ionization) stop 1
  if (MAX_METAL_TRANSITIONS /= 66) stop 2
  call read_parameters_with_implicit_defaults(trim(input_file))
  if (par%metal_collisional_ionization) stop 3
  par%abund_C=1e-4_wp; par%abund_N=0; par%abund_O=0; par%abund_Ne=0
  par%abund_S=0; par%abund_Ar=0; par%abund_Mg=0; par%abund_Fe=0
  par%abund_Si=0; par%abund_Cl=0; par%abund_Ca=0
  par%atomic_dir=trim(directory)
  call species_setup(1)
  ! A band to 1000 eV: accepted for bare C VII, refused for a C VI terminal.
  call species_terminal_photoionization_check(1000.0_wp, eedge, lum)
  if (trim(mode) == 'bound') error stop 'C VI terminal stage accepted'
  gamma=2e-10_wp
  call species_set_photoionization_rates(1,gamma(:,1:1))
  do selection=0,1
    par%metal_collisional_ionization=selection==1
    call species_fractions(1,1,1e5_wp,100.0_wp,0.0_wp,100.0_wp,frac)
    call species_ne_prepare(1,1e5_wp)
    nem=species_ne_cached(100.0_wp,100.0_wp,0.0_wp,100.0_wp)
    expected=0
    do j=1,7
      expected=expected+100.0_wp*elem_abund(1)*real(j-1,wp)*frac(j)
    end do
    if (abs(nem-expected)>1e-13_wp) stop 4
    if (selection==0) previous=frac
    if (selection==1 .and. maxval(abs(frac-previous))<1e-4_wp) stop 5
    cool=metal_ionization_cooling(1,1e5_wp,100.0_wp,100.0_wp,0.0_wp,100.0_wp)
    write(*,'(a,1x,i1,9(1x,es24.16))') 'METAL_STATE',selection,frac,nem,cool
  end do
  call MPI_FINALIZE(ierr)
end program check_optional_metal_physics
