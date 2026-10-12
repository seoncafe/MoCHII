program check_sulfur_charge_exchange_model
! S + H+ <-> S+ + H resolved by S+ product term (Zhao et al. 2005 Table III;
! element_s.txt form 9): the reverse rate and the mean kinetic energies the
! species code forms at (T, n_e), beside the 4S / 2D / 2P term populations
! of the S II atom of the selected excitation model, for the independent
! detailed-balance evaluation in check_sulfur_charge_exchange_model.py.
  use define
  use mpi
  use species_mod, only: species_setup, species_cx_coefficients, species_cx_kinetic_energy
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_cooling_prepT, nlevel_populations, &
       metal_excitation_data_file
  implicit none
  type(nlevel_atom_type) :: atom
  character(len=PATH_LEN) :: directory, output, model
  real(kind=wp), allocatable :: coll(:,:), rad(:,:), pop(:)
  real(kind=wp), parameter :: density(4)=[1e-6_wp,100.0_wp,1e6_wp,1e10_wp]
  real(kind=wp) :: T, ki, kr, di, dr, split, terms(3)
  integer :: ierr, i, j, unit
  logical :: ok
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD
  if(command_argument_count()/=3) error stop 'Expected ATOMIC_DIR MODEL OUTPUT'
  call get_command_argument(1,directory)
  call get_command_argument(2,model)
  call get_command_argument(3,output)
  par%atomic_dir=trim(directory);par%use_metals=.true.;par%cooling_model='low_density'
  par%metal_excitation_model=trim(model)
  par%abund_C=0;par%abund_N=0;par%abund_O=0;par%abund_Ne=0
  par%abund_S=1e-4_wp;par%abund_Ar=0;par%abund_Mg=0;par%abund_Fe=0
  par%abund_Si=0;par%abund_Cl=0;par%abund_Ca=0
  call species_setup(1)
  call nlevel_load(trim(metal_excitation_data_file('nlevel','s',2)),atom,ok)
  if(.not.ok) error stop 'S II atom unreadable'
  allocate(coll(atom%nlev,atom%nlev),rad(atom%nlev,atom%nlev),pop(atom%nlev))
  open(newunit=unit,file=trim(output),status='new')
  write(unit,'(a,i0)') '# nlev ',atom%nlev
  write(unit,'(a)') '# T ne ki kr de_ion_eV de_rec_eV split pop_4S pop_2D pop_2P'
  do i=0,100
    T=10.0_wp**(3.0_wp+0.02_wp*i)
    call nlevel_cooling_prepT(atom,T,coll,rad)
    do j=1,4
      call species_cx_coefficients(1,1,T,ki,kr,density(j))
      call species_cx_kinetic_energy(1,1,di,dr,split,T,density(j))
      call nlevel_populations(atom,coll,rad,density(j),pop)
      terms=[pop(1),sum(pop(2:3)),sum(pop(4:5))]
      write(unit,'(10es27.17e3)') T,density(j),ki,kr,di,dr,split,terms
    end do
  end do
  close(unit)
  call MPI_FINALIZE(ierr)
end program check_sulfur_charge_exchange_model
! S + H+ <-> S+ + H resolved by S+ product term (Zhao et al. 2005 Table III;
! element_s.txt form 9): the reverse rate and the mean kinetic energies the
! species code forms at (T, n_e), beside the 4S / 2D / 2P term populations
! of the S II atom of the selected excitation model, for the independent
! detailed-balance evaluation in check_sulfur_charge_exchange_model.py.
