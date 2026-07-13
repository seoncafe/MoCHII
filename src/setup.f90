module setup_mod
!---------------------------------------------------------------------------
! MoCHII: trimmed driver setup (from MoCafe_v2.00/src/setup.f90, G0 slim).
!
! read_input: namelist /parameters/ par + MPI communicators + validation of
! the AMR + ionizing-band options.  setup_procedure: RNG seed + the AMR
! raytrace procedure pointers.  The dust SED / scattering / peel-off /
! observer bindings return with later stages (G1+); G0 transports only
! ionizing packets through the analytic edge walk.
!---------------------------------------------------------------------------
contains
  !+++++++++++++++++++++++++++++++++++++++++++
  subroutine read_input
  use define
  use utility
  use iofile_mod, only : io_file_extension
  use mpi
  implicit none

  character(len=128) :: model_infile, arg
  integer :: unit, ierr

  namelist /parameters/ par

  if (command_argument_count() >= 1) then
     call get_command_argument(1, model_infile)
  else
     call get_command_argument(0, arg)
     write(*,*) 'Usage: ',trim(arg),' input_file.'
     stop
  endif

  open(newunit=unit,file=trim(model_infile),status='old')
  read(unit,parameters)
  close(unit)

  par%nprint = par%no_print
  if (par%nprint >= par%no_photons) par%nprint = par%no_photons/10
  if (par%no_photons < 10) par%nprint = 1
  par%nphotons   = par%no_photons
  par%nscatt_tot = 0.0_wp

  !--- MPI-related parameters (identical to MoCafe).
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpar%nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, mpar%p_rank, &
                           MPI_INFO_NULL, mpar%hostcomm, ierr)
  call MPI_COMM_RANK(mpar%hostcomm, mpar%h_rank, ierr)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, mpar%h_rank, mpar%p_rank, mpar%SAME_HRANK_COMM, ierr)
  call MPI_COMM_SIZE(mpar%SAME_HRANK_COMM, mpar%SAME_HRANK_NPROC, ierr)

  !--- MoCHII grids: 'amr' (octree) or 'car' (single-level Cartesian, raster
  !--- order; 'uniform' is accepted as a backward-compatible alias).  Either
  !--- grid reads a leaf-list file (par%amr_file); a 'car' grid may instead be
  !--- built from the namelist (par%nx/ny/nz + par%xmax/ymax/zmax).
  if (trim(par%grid_type) == 'uniform') par%grid_type = 'car'
  if (trim(par%grid_type) /= 'amr' .and. trim(par%grid_type) /= 'car') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: MoCHII requires par%grid_type = ''amr'' or ''car''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (len_trim(par%amr_file) == 0) then
     if (trim(par%grid_type) == 'amr') then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: grid_type=''amr'' requires par%amr_file (refinement structure).'
        call MPI_FINALIZE(ierr);  stop
     else if (par%nx < 2 .or. par%xmax <= 0.0_wp) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: grid_type=''car'' without par%amr_file needs par%nx>=2 and par%xmax>0.'
        call MPI_FINALIZE(ierr);  stop
     end if
  endif
  if (trim(par%grid_type) == 'car' .and. par%refine_front) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: refine_front needs the octree (par%grid_type = ''amr'').'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%grid_type) == 'car' .and. &
      trim(par%car_walk) /= 'dda' .and. trim(par%car_walk) /= 'shared') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%car_walk must be ''dda'' (default) or ''shared''.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- ionizing band (G0: the only transport mode).
  if (.not. par%use_ion_band) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: Stage G0 requires par%use_ion_band = .true.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (par%eion_max <= par%eion_min .or. par%nnu_ion < 1) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: need par%eion_max > par%eion_min and par%nnu_ion >= 1.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (len_trim(par%ion_spectrum) == 0 .and. par%tstar <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: use_ion_band requires par%tstar > 0 or par%ion_spectrum.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- dust in the ionizing band needs the EUV extinction table.
  if (par%ion_add_dust .and. len_trim(par%ion_dust_kext) == 0) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%ion_add_dust requires par%ion_dust_kext (EUV kext table).'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- G3: the explicit diffuse field requires case A (case B would
  !--- double-count the on-the-spot absorption of ground recombinations).
  if (par%diffuse_field .and. trim(par%case_ab) /= 'A') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'NOTE: par%diffuse_field forces par%case_ab = ''A''.'
     par%case_ab = 'A'
  endif

  select case(trim(par%distance_unit))
     case ('kpc')
        par%distance2cm = kpc2cm
     case ('pc')
        par%distance2cm = pc2cm
     case ('au')
        par%distance2cm = au2cm
     case ('')
        par%distance2cm = 1.0_wp
     case default
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: unknown par%distance_unit ('''//trim(par%distance_unit)// &
           ''') — use ''kpc'', ''pc'', ''au'', or '''' (cm).'
        call MPI_FINALIZE(ierr);  stop
  end select

  if (len_trim(par%out_file) == 0) then
     par%base_name = trim(get_base_input_name(model_infile))
     par%out_file  = trim(par%base_name)//trim(io_file_extension(par%file_format))
  else
     par%base_name = trim(get_base_name(trim(par%out_file)))
  endif

  if (mpar%p_rank == 0) then
     write(*,'(a)')     ''
     write(*,'(3a)')    '+++++ MoCHII (G0): ',trim(model_infile),' +++++'
     write(*,'(2a)')    ' >>> START @ ', get_date_time()
     write(*,'(a,i14)') 'Total ionizing photons    : ', par%nphotons
  endif
  end subroutine read_input

  !---------------------------------------
  subroutine setup_procedure
  use define
  use random, only : init_random_seed
  use raytrace_amr_mod
  implicit none
  !--- Initialize Random Number Generator
  call init_random_seed(par%iseed)
  !--- AMR raytrace bindings (used by later stages; G0 calls the ionizing
  !--- edge walk directly).
  raytrace_to_tau  => raytrace_to_tau_amr
  raytrace_to_edge => raytrace_to_edge_amr
  end subroutine setup_procedure
end module setup_mod
