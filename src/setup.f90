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
  use ion_band_mod, only : ion_ext_preset_id
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

  !--- Source model.  The ionizing band is fed by any number of internal point
  !--- sources (par%nsource, positions src_x/y/z, luminosities src_lum) plus,
  !--- INDEPENDENTLY, an isotropic external field (ON when par%ext_intensity>0,
  !--- entry geometry par%ext_geometry = 'rec'|'sph').  Packets split among all
  !--- components in proportion to their band luminosity, each with its own
  !--- spectrum.  par%source_geometry is a LEGACY ALIAS: 'point' (default) leaves
  !--- the composable form untouched; 'external'|'external_rec'|'external_sph' is
  !--- the external-ONLY shorthand (forces nsource=0, external on).
  block
    integer :: is, nset
    logical :: ext_on, is_phys, ext_is_preset, ext_is_physfile, src_is_absolute

    !--- par%spectrum_type sets the column units of every spectrum file slot.
    select case (trim(par%spectrum_type))
    case ('shape', 'per_ev', 'per_hz', 'per_ang', 'per_um')
       is_phys = trim(par%spectrum_type) /= 'shape'
    case default
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: par%spectrum_type must be ''shape'' (default), ''per_ev'', ''per_hz'', '// &
          '''per_ang'', or ''per_um''.'
       call MPI_FINALIZE(ierr);  stop
    end select

    !--- the ISRF presets (draine/habing/mathis) are external-only; an internal
    !--- file slot naming one is almost certainly a mistake.
    if (ion_ext_preset_id(par%ion_spectrum) > 0 .or. &
        ion_ext_preset_id(par%src_spectrum_file) > 0) then
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: ISRF preset names (draine/habing/mathis) are valid only for '// &
          'par%ext_spectrum, not for an internal source file slot.'
       call MPI_FINALIZE(ierr);  stop
    endif
    ext_is_preset   = ion_ext_preset_id(par%ext_spectrum) > 0
    ext_is_physfile = (.not. ext_is_preset) .and. len_trim(par%ext_spectrum) > 0 .and. is_phys
    !--- the primary internal source is ABSOLUTE when a physical-type spectrum
    !--- file feeds it (its own src_spectrum_file or the global ion_spectrum).
    src_is_absolute = is_phys .and. (len_trim(par%src_spectrum_file) > 0 .or. &
                                     len_trim(par%ion_spectrum) > 0)

    if (trim(par%source_geometry) == 'external' .or. &
        trim(par%source_geometry) == 'external_rec' .or. &
        trim(par%source_geometry) == 'external_sph') then
       !--- legacy external-only alias: map onto ext_geometry + nsource=0.
       if (trim(par%source_geometry) == 'external_sph') then
          par%ext_geometry = 'sph'
       else
          par%ext_geometry = 'rec'
       end if
       if (par%ext_intensity <= 0.0_wp .and. .not. (ext_is_preset .or. ext_is_physfile)) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: external source_geometry requires par%ext_intensity > 0 '// &
             '(mean intensity J), an ISRF preset, or a physical-type ext_spectrum.'
          call MPI_FINALIZE(ierr);  stop
       endif
       par%nsource = 0
       if (mpar%p_rank == 0) write(*,'(a)') &
          'NOTE: source_geometry=''external*'' is the external-only shorthand; '// &
          'the composable form is par%ext_intensity/par%ext_geometry + par%nsource.'
    else if (trim(par%source_geometry) /= 'point') then
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: par%source_geometry must be ''point'', ''external_rec'', or ''external_sph''.'
       call MPI_FINALIZE(ierr);  stop
    endif

    !--- external field ON: mean intensity, an ISRF preset, or a physical-type
    !--- ext_spectrum file (the last two are absolute and need no ext_intensity).
    ext_on = (par%ext_intensity > 0.0_wp) .or. ext_is_preset .or. ext_is_physfile

    !--- ISRF presets are analytic and FUV-only: require add_fuv, a positive
    !--- scale, and note that ext_intensity (if set) overrides ext_scale.
    if (ext_is_preset) then
       if (par%ext_scale <= 0.0_wp) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: par%ext_scale must be > 0 for an ISRF preset.'
          call MPI_FINALIZE(ierr);  stop
       endif
       if (.not. par%add_fuv) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: the ISRF presets are FUV-only; set par%add_fuv (with '// &
             'par%efuv_min < 13.6) so the FUV bins carry the field.'
          call MPI_FINALIZE(ierr);  stop
       endif
       if (par%ext_intensity > 0.0_wp .and. par%ext_scale /= 1.0_wp &
           .and. mpar%p_rank == 0) write(*,'(a)') &
          'NOTE: par%ext_intensity rescales the preset to a target band J; '// &
          'par%ext_scale is then ignored.'
    endif
    !--- external geometry validation (only when the external field is on).
    if (ext_on) then
       if (trim(par%ext_geometry) /= 'rec' .and. trim(par%ext_geometry) /= 'sph') then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: par%ext_geometry must be ''rec'' or ''sph''.'
          call MPI_FINALIZE(ierr);  stop
       endif
       if (trim(par%ext_geometry) == 'sph' .and. par%rmax <= 0.0_wp) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: ext_geometry=''sph'' requires par%rmax > 0.'
          call MPI_FINALIZE(ierr);  stop
       endif
    endif

    !--- luminosity sentinel (default -999 = unset).  A 'shape'/Planck primary
    !--- source maps it to 1.0 (exact legacy); a physical-type (absolute) source
    !--- leaves it unset so its luminosity is DERIVED from the file integral.
    if (par%luminosity < -900.0_wp) then
       if (par%nsource >= 1 .and. src_is_absolute) then
          continue                          ! leave unset -> derive in ion_band
       else
          par%luminosity = 1.0_wp
       endif
    endif

    !--- internal point sources.
    if (par%nsource < 0 .or. par%nsource > MAX_SRC) then
       if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
          'ERROR: par%nsource must be in [0, ', MAX_SRC, '].'
       call MPI_FINALIZE(ierr);  stop
    endif
    if (par%nsource == 0) then
       if (.not. ext_on) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: no source: set par%nsource >= 1 or the external field '// &
             '(par%ext_intensity > 0, an ISRF preset, or a physical ext_spectrum).'
          call MPI_FINALIZE(ierr);  stop
       endif
    else if (par%nsource == 1) then
       !--- single point source: normalize the legacy scalars into the arrays
       !--- (only read on the multi-component path; the single-component fast
       !--- path uses the legacy scalars directly, so nothing changes there).
       if (par%src_lum(1) <= 0.0_wp) then
          par%src_x(1) = par%xs_point
          par%src_y(1) = par%ys_point
          par%src_z(1) = par%zs_point
          par%src_lum(1)      = par%luminosity
          par%src_geometry(1) = 'point'
       endif
       if (trim(par%src_geometry(1)) /= 'point') then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: the ionizing band supports only src_geometry = ''point''.'
          call MPI_FINALIZE(ierr);  stop
       endif
    else
       !--- multiple point sources.  src_lum: all set (keep) or all unset.  For
       !--- a 'shape' spectrum an unset set is the equal split of par%luminosity
       !--- (MoCafe convention); for a physical-type (absolute) source file an
       !--- unset set is DERIVED per column (src_lum left as the -999 sentinel).
       !--- A partial set is an error.
       nset = 0
       do is = 1, par%nsource
          if (par%src_lum(is) > 0.0_wp) nset = nset + 1
          if (trim(par%src_geometry(is)) /= 'point') then
             if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
                'ERROR: src_geometry(', is, ') must be ''point'' in the ionizing band.'
             call MPI_FINALIZE(ierr);  stop
          endif
       end do
       if (nset == 0) then
          if (.not. src_is_absolute) then
             do is = 1, par%nsource
                par%src_lum(is) = par%luminosity / real(par%nsource, wp)
             end do
          endif
       else if (nset /= par%nsource) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: set par%src_lum for ALL sources or for NONE '// &
             '(equal split, or derived for a physical-type source file).'
          call MPI_FINALIZE(ierr);  stop
       endif
       if (nset == par%nsource) par%luminosity = sum(par%src_lum(1:par%nsource))  ! informational total
    endif

    !--- route the external-only run (no point source) through the legacy
    !--- external fast path by encoding the geometry in source_geometry.
    if (par%nsource == 0) par%source_geometry = 'external_'//trim(par%ext_geometry)

    !--- src_spectrum_file feeds the internal point sources only; it has no
    !--- effect on an external-only run (nsource=0).
    if (len_trim(par%src_spectrum_file) > 0 .and. par%nsource == 0 &
        .and. mpar%p_rank == 0) write(*,'(a)') &
       'NOTE: par%src_spectrum_file is ignored on an external-only run (nsource=0).'
  end block

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
  if (len_trim(par%ion_spectrum) == 0 .and. par%tstar <= 0.0_wp .and. &
      len_trim(par%src_spectrum_file) == 0 .and. &
      len_trim(par%ext_spectrum) == 0 .and. par%ext_tstar <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: use_ion_band needs a spectrum: par%tstar > 0, par%ion_spectrum, '// &
        'par%src_spectrum_file, par%ext_spectrum, or par%ext_tstar > 0.'
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
