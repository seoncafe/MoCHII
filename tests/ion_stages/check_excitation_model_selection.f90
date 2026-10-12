program check_excitation_model_selection
  use define
  use setup_mod, only: read_parameters_with_implicit_defaults
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, metal_excitation_data_file
  use mpi
  implicit none
  type(nlevel_atom_type) :: atom
  character(len=PATH_LEN) :: input, root, output, filename
  character(len=2),parameter :: elements(11)=[character(len=2)::'c','n','o','ne','s','ar','mg','fe','si','cl','ca']
  integer :: ierr,ie,stage,u
  logical :: ok
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  call get_command_argument(1,input)
  call get_command_argument(2,root)
  call get_command_argument(3,output)
  call read_parameters_with_implicit_defaults(trim(input))
  par%atomic_dir=trim(root)
  open(newunit=u,file=trim(output),status='new')
  write(u,'(a)') trim(par%metal_excitation_model)
  do ie=1,11
     do stage=1,7
        filename=metal_excitation_data_file('nlevel',trim(elements(ie)),stage)
        call nlevel_load(trim(filename),atom,ok)
        if (.not.ok) cycle
        write(u,'(a,1x,5(i0,1x))') trim(elements(ie)),stage,atom%nlev,atom%nrad,atom%nups,atom%nhc
     enddo
  enddo
  close(u)
  call MPI_FINALIZE(ierr)
end program check_excitation_model_selection
