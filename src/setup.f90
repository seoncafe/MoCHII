module setup_mod
!---------------------------------------------------------------------------
! MoCHII: driver setup.
!
! read_input: namelist /parameters/ par + MPI communicators + validation of
! the AMR + ionizing-band options.  setup_procedure: RNG seed + the AMR
! raytrace procedure pointers.
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
  character(len=256) :: exepath
  integer :: unit, ierr, islash

  namelist /parameters/ par

  if (command_argument_count() >= 1) then
     call get_command_argument(1, model_infile)
  else
     call get_command_argument(0, arg)
     write(*,*) 'Usage: ',trim(arg),' input_file.'
     stop
  endif

  par%require_convergence = .false.   ! definite baseline; namelist may override
  open(newunit=unit,file=trim(model_infile),status='old')
  read(unit,parameters)
  close(unit)

  !--- Resolve the SEDust directory relative to the executable when the user
  !--- leaves par%sed_workdir blank, so a fresh checkout runs dust emission from
  !--- any working directory without editing paths.  argv(0) is the path to
  !--- MoCHII.x; take its directory and append SEDust/sed.
  if (len_trim(par%sed_workdir) == 0) then
     call get_command_argument(0, exepath)
     islash = index(trim(exepath), '/', back=.true.)
     if (islash > 0) then
        par%sed_workdir = exepath(1:islash-1) // '/SEDust/sed'
     else
        par%sed_workdir = 'SEDust/sed'
     end if
  end if

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
  !--- par%density_file: a uniform 3D density cube (FITS/HDF5) for the 'car'
  !--- grid; nx/ny/nz come from the file, so it is an alternative to both the
  !--- namelist nx/ny/nz and the amr leaf-list.
  if (len_trim(par%density_file) > 0 .and. trim(par%grid_type) /= 'car') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%density_file (a 3D density cube) requires par%grid_type = ''car''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (len_trim(par%density_file) > 0 .and. len_trim(par%amr_file) > 0) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: set either par%amr_file (leaf list) or par%density_file (cube), not both.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (len_trim(par%amr_file) == 0) then
     if (trim(par%grid_type) == 'amr') then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: grid_type=''amr'' requires par%amr_file (refinement structure).'
        call MPI_FINALIZE(ierr);  stop
     else if (len_trim(par%density_file) > 0) then
        if (par%xmax <= 0.0_wp) then
           if (mpar%p_rank == 0) write(*,'(a)') &
              'ERROR: grid_type=''car'' with par%density_file needs par%xmax > 0.'
           call MPI_FINALIZE(ierr);  stop
        end if
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
  if (trim(par%ci_model) /= 'voronov' .and. &
      trim(par%ci_model) /= 'dere_hybrid') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%ci_model must be ''voronov'' (default) or ''dere_hybrid''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%recomb_model) /= 'badnell_mao' .and. &
      trim(par%recomb_model) /= 'hui_gnedin') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%recomb_model must be ''badnell_mao'' (default) or ''hui_gnedin''.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- general parameter sanity (fail fast rather than run on nonsense).
  if (par%no_photons < 1.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%no_photons must be >= 1.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (par%ion_relax <= 0.0_wp .or. par%ion_relax > 1.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%ion_relax must be in (0, 1].'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%case_ab) /= 'A' .and. trim(par%case_ab) /= 'B') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%case_ab must be ''A'' or ''B''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (par%solve_te .and. par%te_min >= par%te_max) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%solve_te needs par%te_min < par%te_max.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (par%ion_peel .and. (par%nxim < 1 .or. par%nyim < 1)) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%ion_peel needs par%nxim >= 1 and par%nyim >= 1.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- ionizing band (the ionizing transport mode).
  if (.not. par%use_ion_band) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: MoCHII requires par%use_ion_band = .true.'
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

  !--- the explicit diffuse field requires case A (case B would
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
     write(*,'(3a)')    '+++++ MoCHII: ',trim(model_infile),' +++++'
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
  !--- AMR raytrace bindings (the ionizing edge walk is called directly).
  raytrace_to_tau  => raytrace_to_tau_amr
  raytrace_to_edge => raytrace_to_edge_amr
  end subroutine setup_procedure
end module setup_mod
