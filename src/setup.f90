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

  character(len=PATH_LEN) :: model_infile, arg
  character(len=PATH_LEN) :: exepath
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
  call read_parameters_with_implicit_defaults(trim(model_infile))

  !--- Resolve the SEDust data directory relative to the executable when the
  !--- user leaves par%sed_data_dir blank, so a fresh checkout runs dust
  !--- emission from any working directory without editing paths.  argv(0) is
  !--- the path to MoCHII.x; take its directory and append SEDust/data.
  if (len_trim(par%sed_data_dir) == 0) then
     call get_command_argument(0, exepath)
     islash = index(trim(exepath), '/', back=.true.)
     if (islash > 0) then
        par%sed_data_dir = exepath(1:islash-1) // '/SEDust/data'
     else
        par%sed_data_dir = 'SEDust/data'
     end if
  end if

  !--- Same treatment for the shipped data tables.  They used to be opened
  !--- as '../../data/x' with a 'data/x' retry, which resolves against the
  !--- CURRENT WORKING DIRECTORY: the same run then succeeds or dies
  !--- depending on where it was launched from, and the failure surfaces
  !--- only after the transport has already been paid for.
  if (len_trim(par%data_dir) == 0) then
     call get_command_argument(0, exepath)
     islash = index(trim(exepath), '/', back=.true.)
     if (islash > 0) then
        par%data_dir = exepath(1:islash-1) // '/data'
     else
        par%data_dir = 'data'
     end if
  end if
  if (len_trim(par%atomic_dir) == 0) &
     par%atomic_dir = trim(par%data_dir) // '/atomic'

  par%nprint = par%no_print
  if (par%nprint >= par%no_photons) par%nprint = par%no_photons/10
  if (par%no_photons < 10) par%nprint = 1
  par%nphotons   = par%no_photons
  par%nphotons_fuv = par%no_photons_fuv
  par%nscatt_tot = 0.0_wp

  !--- MPI-related parameters (identical to MoCafe).
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpar%nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, mpar%p_rank, &
                           MPI_INFO_NULL, mpar%hostcomm, ierr)
  call MPI_COMM_RANK(mpar%hostcomm, mpar%h_rank, ierr)
  call MPI_COMM_SIZE(mpar%hostcomm, mpar%h_size, ierr)
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
     else if (par%xmax <= 0.0_wp .or. &
              (par%nx < 2 .and. .not. par%xy_periodic)) then
        !--- nx=1 is allowed only with xy_periodic (a plane-parallel column
        !--- whose single lateral cell wraps onto itself).
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: grid_type=''car'' without par%amr_file needs par%xmax>0 and '// &
           'par%nx>=2 (par%nx>=1 allowed with par%xy_periodic).'
        call MPI_FINALIZE(ierr);  stop
     else if (par%ny < 1 .or. par%nz < 1 .or. &
              par%ymax <= 0.0_wp .or. par%zmax <= 0.0_wp) then
        !--- the remaining namelist-built dimensions must be positive.
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: grid_type=''car'' without par%amr_file needs par%ny>=1, '// &
           'par%nz>=1, par%ymax>0, and par%zmax>0.'
        call MPI_FINALIZE(ierr);  stop
     end if
  endif
  if (trim(par%grid_type) == 'car' .and. par%refine_front) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: refine_front needs the octree (par%grid_type = ''amr'').'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (par%refine_front) then
     if (par%refine_lbase > par%refine_lmax) then
        if (mpar%p_rank == 0) write(*,'(a,i0,a,i0,a)') &
           'ERROR: par%refine_lbase (', par%refine_lbase, &
           ') must not exceed par%refine_lmax (', par%refine_lmax, ').'
        call MPI_FINALIZE(ierr);  stop
     endif
     if (par%refine_lbase < 1 .or. par%refine_lmax < 1) then
        if (mpar%p_rank == 0) write(*,'(a)') &
           'ERROR: par%refine_lbase and par%refine_lmax must be >= 1.'
        call MPI_FINALIZE(ierr);  stop
     endif
     if (par%refine_iter < 1 .or. par%refine_iter >= par%gas_niter) then
        if (mpar%p_rank == 0) write(*,'(a,i0,a,i0,a)') &
           'ERROR: par%refine_iter (', par%refine_iter, &
           ') must be in [1, par%gas_niter-1] so the new grid is iterated '// &
           'on; par%gas_niter = ', par%gas_niter, '.'
        call MPI_FINALIZE(ierr);  stop
     endif
  endif
  if (trim(par%grid_type) == 'car' .and. &
      trim(par%car_walk) /= 'dda' .and. trim(par%car_walk) /= 'amr') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%car_walk must be ''dda'' (default) or ''amr''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%ci_model) /= 'voronov' .and. &
      trim(par%ci_model) /= 'dere_hybrid') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%ci_model must be ''voronov'' (default) or ''dere_hybrid''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%charge_exchange) /= 'two_way' .and. &
      trim(par%charge_exchange) /= 'metal_only') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%charge_exchange must be ''two_way'' (default) or '// &
        '''metal_only'' (reproduction only).'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%recomb_model) /= 'detailed_balance' .and. &
      trim(par%recomb_model) /= 'hui_gnedin') then
     if (mpar%p_rank == 0) then
        write(*,'(a)') 'ERROR: par%recomb_model must be ''detailed_balance'' '// &
           '(default) or ''hui_gnedin''.'
        if (trim(par%recomb_model) == 'badnell_milne') write(*,'(a)') &
           '       ''badnell_milne'' was renamed ''detailed_balance'' on '// &
           '2026-09-26 (the hydrogenic ions no longer use the Badnell total).'
     end if
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%diffuse_energy_model) /= 'exponential' .and. &
      trim(par%diffuse_energy_model) /= 'threshold'   .and. &
      trim(par%diffuse_energy_model) /= 'milne') then
     if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: par%diffuse_energy_model must be ''milne'', ' // &
          '''exponential'' or ''threshold''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%cooling_model) /= 'low_density' .and. &
      trim(par%cooling_model) /= 'local_ne') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%cooling_model must be ''local_ne'' (default) or ''low_density''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%twophoton_quenching) /= 'colliders' .and. &
      trim(par%twophoton_quenching) /= 'lumped_ne') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%twophoton_quenching must be ''colliders'' (default) or ''lumped_ne''.'
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
    case ('shape_ev', 'per_ev', 'per_hz', 'per_ang', 'per_um')
       is_phys = trim(par%spectrum_type) /= 'shape_ev'
    case default
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: par%spectrum_type must be ''shape_ev'' (default), ''per_ev'', ''per_hz'', '// &
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
    else if (trim(par%source_geometry) == 'slab') then
       !--- plane-parallel slab illumination (its own single-component source).
       if (.not. par%xy_periodic) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: source_geometry=''slab'' needs par%xy_periodic (a slab is '// &
             'periodic in x and y).'
          call MPI_FINALIZE(ierr);  stop
       endif
       if (trim(par%slab_faces) /= 'top' .and. trim(par%slab_faces) /= 'bottom' &
           .and. trim(par%slab_faces) /= 'both') then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: par%slab_faces must be ''top'', ''bottom'', or ''both''.'
          call MPI_FINALIZE(ierr);  stop
       endif
       if ((trim(par%slab_top_mode) /= 'beam' .and. trim(par%slab_top_mode) /= 'isotropic') .or. &
           (trim(par%slab_bot_mode) /= 'beam' .and. trim(par%slab_bot_mode) /= 'isotropic')) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: par%slab_top_mode / slab_bot_mode must be ''beam'' or ''isotropic''.'
          call MPI_FINALIZE(ierr);  stop
       endif
       block
         logical :: ltop, lbot
         ltop = trim(par%slab_faces) == 'top'    .or. trim(par%slab_faces) == 'both'
         lbot = trim(par%slab_faces) == 'bottom' .or. trim(par%slab_faces) == 'both'
         !--- a beam incidence angle theta -> 90 deg gives mu = cos theta -> 0,
         !--- so kz -> 0 and the packet would wrap the periodic x/y faces without
         !--- ever reaching a z-face.  Require 0 <= theta < 90 on each lit beam face.
         if (ltop .and. trim(par%slab_top_mode) == 'beam' .and. &
             (par%slab_top_theta < 0.0_wp .or. par%slab_top_theta >= 90.0_wp)) then
            if (mpar%p_rank == 0) write(*,'(a)') &
               'ERROR: par%slab_top_theta must be in [0, 90) degrees for a beam.'
            call MPI_FINALIZE(ierr);  stop
         endif
         if (lbot .and. trim(par%slab_bot_mode) == 'beam' .and. &
             (par%slab_bot_theta < 0.0_wp .or. par%slab_bot_theta >= 90.0_wp)) then
            if (mpar%p_rank == 0) write(*,'(a)') &
               'ERROR: par%slab_bot_theta must be in [0, 90) degrees for a beam.'
            call MPI_FINALIZE(ierr);  stop
         endif
         !--- each lit face must carry a positive source strength on its own (a
         !--- negative one could hide behind a positive one in the sum and push
         !--- slab_ptop outside [0, 1]).
         if (ltop .and. par%slab_top_source <= 0.0_wp) then
            if (mpar%p_rank == 0) write(*,'(a)') &
               'ERROR: par%slab_top_source must be > 0 on the lit top face.'
            call MPI_FINALIZE(ierr);  stop
         endif
         if (lbot .and. par%slab_bot_source <= 0.0_wp) then
            if (mpar%p_rank == 0) write(*,'(a)') &
               'ERROR: par%slab_bot_source must be > 0 on the lit bottom face.'
            call MPI_FINALIZE(ierr);  stop
         endif
       end block
       if (par%slab_nmu < 1) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: par%slab_nmu must be >= 1.'
          call MPI_FINALIZE(ierr);  stop
       endif
       !--- slab illumination is its own single source; composing it with
       !--- internal point sources or an external field is not yet supported.
       !--- nsource defaults to 1, so only nsource > 1 (or an external field)
       !--- signals an explicit extra-source request; nsource <= 1 is silently
       !--- overridden by the slab below.
       if (par%nsource > 1 .or. par%ext_intensity > 0.0_wp) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: source_geometry=''slab'' cannot be combined with internal '// &
             'point sources (par%nsource > 1) or an external field '// &
             '(par%ext_intensity > 0); slab+other-source composition is not yet supported.'
          call MPI_FINALIZE(ierr);  stop
       endif
       par%nsource = 0
    else if (trim(par%source_geometry) /= 'point') then
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: par%source_geometry must be ''point'', ''slab'', ''external_rec'', or ''external_sph''.'
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

    !--- luminosity sentinel (default -999 = unset).  A 'shape_ev'/Planck primary
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
    if (par%nsource == 0 .and. trim(par%source_geometry) /= 'slab') then
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
       !--- a 'shape_ev' spectrum an unset set is the equal split of par%luminosity
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
    !--- (the slab is its own geometry and is left untouched.)
    if (par%nsource == 0 .and. trim(par%source_geometry) /= 'slab') &
       par%source_geometry = 'external_'//trim(par%ext_geometry)

    !--- src_spectrum_file feeds the internal point sources only; it has no
    !--- effect on an external-only run (nsource=0).
    if (len_trim(par%src_spectrum_file) > 0 .and. par%nsource == 0 &
        .and. mpar%p_rank == 0) write(*,'(a)') &
       'NOTE: par%src_spectrum_file is ignored on an external-only run (nsource=0).'

    !--- launch sequence: 'random' (Mersenne Twister) or 'sobol' (Owen-
    !--- scrambled quasi-random launch).  'sobol' now covers every source
    !--- configuration - single point, single external (rec/sph), multiple
    !--- points, and mixed point+external - through the fixed superset launch
    !--- dimension layout; the diffuse packets ride a second decorrelated Sobol
    !--- stream; the dust-scattering / SED draws keep the Mersenne Twister.
    if (trim(par%launch_sequence) /= 'random' .and. &
        trim(par%launch_sequence) /= 'sobol') then
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: par%launch_sequence must be ''sobol'' (default) or ''random''.'
       call MPI_FINALIZE(ierr);  stop
    endif
  end block

  !--- general parameter sanity (fail fast rather than run on nonsense).
  if (par%no_photons < 1.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%no_photons must be >= 1.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%fuv_sampling) /= 'separate' .and. &
      trim(par%fuv_sampling) /= 'proportional') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%fuv_sampling must be ''separate'' (default) or ''proportional''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (par%add_fuv .and. trim(par%fuv_sampling) == 'separate' .and. &
      par%no_photons_fuv < 1.0_wp .and. par%no_photons_fuv /= -1.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%no_photons_fuv must be >= 1 with par%add_fuv and '// &
        'par%fuv_sampling = ''separate''.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (par%source_cdf_tol <= 0.0_wp .or. par%source_cdf_tol >= 1.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%source_cdf_tol must be in (0, 1).'
     call MPI_FINALIZE(ierr);  stop
  endif
  !--- Energy sampling covers every source configuration (point sources,
  !--- external field, slab), so the one requirement is that SOME spectrum is
  !--- defined for the components to sample.
  if (len_trim(par%ion_spectrum) == 0 .and. &
      len_trim(par%src_spectrum_file) == 0 .and. &
      len_trim(par%ext_spectrum) == 0 .and. &
      par%tstar <= 0.0_wp .and. par%ext_tstar <= 0.0_wp .and. &
      all(par%src_tstar <= 0.0_wp)) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: the ionizing band needs a Planck or tabulated source spectrum.'
     call MPI_FINALIZE(ierr);  stop 1
  endif
  if (par%ion_relax <= 0.0_wp .or. par%ion_relax > 1.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%ion_relax must be in (0, 1].'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%lyman_lines) /= 'on_the_spot' .and. &
      trim(par%lyman_lines) /= 'escaping') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%lyman_lines must be ''on_the_spot'' (default) or ' // &
        '''escaping'' (MoCHII does not transport the Lyman lines).'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (trim(par%recombination_continuum) /= 'on_the_spot' .and. &
      trim(par%recombination_continuum) /= 'transported' .and. &
      trim(par%recombination_continuum) /= 'escaping') then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%recombination_continuum must be ''on_the_spot'', ' // &
        '''transported'' (default) or ''escaping''.'
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
  if ((par%hei_metastable .or. par%hei_metastable_coupling) .and. &
      par%abund_He > 0.0_wp) then
     !--- Phi_3 (2^3S photoionization) starts at E_HEI_23S_SIGMA_EDGE = 4.78
     !--- eV; a band floor above it misses the soft-UV part of the integral
     !--- (rank-0 note, no abort).  Phi_3 feeds the coupled balances too, so
     !--- the note follows either switch.  For a 40 kK star the unattenuated
     !--- 4.78-13.6 eV part is ~10x the part above 13.6 eV.  An input that
     !--- leaves add_fuv and efuv_min out gets exactly this floor
     !--- (read_parameters_with_implicit_defaults), so the note fires only
     !--- when the input sets one of them otherwise.
     block
       real(kind=wp) :: efloor
       efloor = par%eion_min
       if (par%add_fuv) efloor = par%efuv_min
       if (mpar%p_rank == 0 .and. efloor > E_HEI_23S_SIGMA_EDGE) &
          write(*,'(a,f7.3,a,f4.2,a,f4.2,a,f4.2,a)') &
          'NOTE: hei_metastable band floor ', efloor, ' eV > ', &
          E_HEI_23S_SIGMA_EDGE, ' eV -> Phi_3 misses [', E_HEI_23S_SIGMA_EDGE, &
          ', floor); set add_fuv with efuv_min = ', E_HEI_23S_SIGMA_EDGE, &
          ' (the default when both are left out) for the full 2^3S '// &
          'photoionization integral.'
     end block
  endif
  if (par%ion_peel .and. par%xy_periodic) then
     !--- the point-observer peel is ill-defined under xy-periodicity (the
     !--- medium tiles into infinitely many periodic images); the emergent
     !--- observable of a periodic slab is the boundary I(mu) tally instead.
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%ion_peel is not supported with par%xy_periodic (a periodic '// &
        'slab has no finite-distance observer; use the boundary I(mu) output '// &
        'written for source_geometry=''slab'').'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- ionizing band
  if (par%eion_max <= par%eion_min .or. par%nnu_ion < 1) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: need par%eion_max > par%eion_min and par%nnu_ion >= 1.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (len_trim(par%ion_spectrum) == 0 .and. par%tstar <= 0.0_wp .and. &
      len_trim(par%src_spectrum_file) == 0 .and. &
      len_trim(par%ext_spectrum) == 0 .and. par%ext_tstar <= 0.0_wp) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: the ionizing band needs a spectrum: par%tstar > 0, par%ion_spectrum, '// &
        'par%src_spectrum_file, par%ext_spectrum, or par%ext_tstar > 0.'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- dust in the ionizing band needs extinction optics: either the EUV kext
  !--- table (par%ion_dust_kext, which OVERRIDES the model) or the SEDust grain
  !--- model (par%dust_model), extended into the EUV.  gas_opacity_mod
  !--- checks that the grain model's grid really covers the band, for every
  !--- model alike (Zubko needs no extension -- its own optics reach 0.001 um).
  if (par%ion_add_dust .and. len_trim(par%ion_dust_kext) == 0 .and. &
      mpar%p_rank == 0) write(*,'(3a)') &
     ' NOTE: par%ion_add_dust with no par%ion_dust_kext -> ionizing-band dust '// &
     'optics from the ''', trim(par%dust_model), ''' SEDust grain model.'

  if (par%grain_pe .and. .not. par%add_fuv) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%grain_pe needs par%add_fuv (the FUV field is its heating source).'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (par%grain_pe .and. .not. par%ion_add_dust) then
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: par%grain_pe needs par%ion_add_dust (grain opacity carries the field).'
     call MPI_FINALIZE(ierr);  stop
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
     !--- how the photon indices are split over the ranks: .false. = static
     !--- cyclic, .true. = rank 0 serves batches on demand.
     write(*,'(a,l14)') 'Master-slave scheduling   : ', par%use_master_slave
     write(*,'(a,i14)') 'Total ionizing photons    : ', par%nphotons
     if (par%add_fuv) then
        write(*,'(2a)')    'FUV sampling              : ', trim(par%fuv_sampling)
        if (trim(par%fuv_sampling) == 'separate') then
           if (par%nphotons_fuv >= 0_int64) then
              write(*,'(a,i14)') 'FUV photons               : ', par%nphotons_fuv
           else
              write(*,'(a)') 'FUV photons               :  set by the band luminosities'
           end if
        end if
     end if
     write(*,'(a,l14)') 'Secondary ionization      : ', par%use_sec_ion
     write(*,'(a,l14)') 'He I 2^3S coupling        : ', par%hei_metastable_coupling
     if (par%ion_add_dust) &
        write(*,'(a,l14)') 'Dust emission transported : ', par%dust_emis_transport
     write(*,'(a,l14)') 'Metal free-free cooling   : ', par%metal_freefree
     write(*,'(2a)')    'Recombination continuum   : ', trim(par%recombination_continuum)
     write(*,'(2a)')    'Lyman lines               : ', trim(par%lyman_lines)
     if (trim(par%recombination_continuum) == 'on_the_spot') then
        write(*,'(a,l14)') 'He I on the spot          : ', par%hei_on_the_spot
        write(*,'(a,l14)') 'He II on the spot         : ', par%heii_on_the_spot
        write(*,'(a,l14)') 'He III cascade transported: ', par%heii_diffuse
        !--- heii_recombination_photons_mod, "Absorption limit"
        if (par%heii_on_the_spot .and. .not. par%heii_diffuse .and. &
            par%abund_He > 0.12_wp) &
           write(*,'(a,f6.3,a)') 'NOTE: abund_He = ', par%abund_He, ' > 0.12: in a'// &
              ' fully He2+ zone the He II cascade on the spot can outnumber'// &
              ' hydrogen''s recombinations, where no local balance exists'// &
              ' (the H balance residual shows those leaves); use'// &
              ' recombination_continuum = ''transported'' there.'
     end if
     write(*,'(a,l14)') 'Charge-exchange heating   : ', par%charge_exchange_heating
     write(*,'(2a)')    'Photon launch sequence    : ', trim(par%launch_sequence)
     write(*,'(a,es14.6)') 'Source CDF tolerance      : ', par%source_cdf_tol
     if (trim(par%launch_sequence) == 'sobol') &
        write(*,'(a,i14)') 'QMC scramble seed         : ', par%qmc_seed
  endif
  end subroutine read_input

  !=========================================================================
  ! Read the namelist, then set the switches whose defaults depend on other
  ! items: the metal switches the input left out from the metal abundances it
  ! gave, and the excited-cascade diffuse channels from
  ! recombination_continuum.
  !   use_metals            not in the input: .true. if any abund_<X> the
  !                         input sets is > 0, else .false. (every abundance
  !                         defaults to 0, so an input that names no
  !                         abundance stays H/He only);
  !   metal_ne, metal_heat  not in the input: follow use_metals;
  !   hei_diffuse           not in the input: on when recombination_continuum
  !                         = 'transported';
  !   heii_diffuse          not in the input: the same, when the
  !                         band reaches the 40.813 eV He II Ly-alpha
  !                         (eion_max), else .false. with a note; set in
  !                         the input with 'on_the_spot', it transports the
  !                         He III cascade alone (heii_cascade_transported);
  !   add_fuv               not in the input: .true. when the He I 2^3S
  !                         population is used (hei_metastable or
  !                         hei_metastable_coupling), else .false.: Phi_3
  !                         is a path estimator over the transported band
  !                         and starts at 4.78 eV;
  !   efuv_min              not in the input: E_HEI_23S_SIGMA_EDGE = 4.78 eV
  !                         when the 2^3S population is used, else 6.0 eV;
  !   no_photons_fuv        not in the input: no_photons/4, or no_photons when
  !                         no source has ionizing luminosity (resolved in
  !                         ion_band_mod once the band luminosities are known).
  ! A value the input does set is kept as given.  Fortran namelist input does
  ! not report which items it contained, so the file is read twice, from two
  ! states that differ only in these items (the switches .true. then
  ! .false., the abundances, efuv_min and no_photons_fuv -1 then -2): an item
  ! the input sets comes out equal in both reads, an item it leaves out does
  ! not.
  !=========================================================================
  subroutine read_parameters_with_implicit_defaults(fname)
    use define
    use mpi
    implicit none
    character(len=*), intent(in) :: fname
    integer, parameter :: NAB = 11
    character(len=2), parameter :: ABNAME(NAB) = &
       ['C ', 'N ', 'O ', 'Ne', 'S ', 'Ar', 'Mg', 'Fe', 'Si', 'Cl', 'Ca']
    type(params_type) :: par0, par_hi
    real(kind=wp) :: ab_hi(NAB), ab_lo(NAB), ab_default(NAB)
    logical :: given_um, given_ne, given_heat, given_ab(NAB), any_positive
    logical :: given_hei, given_heii, transported
    logical :: given_fuv, given_efuv, given_nfuv, uses_23s
    !--- the He II Ly-alpha energy [eV] (heii_cascade_mod::E_HEII_LYA)
    real(kind=wp), parameter :: E_HEII_LYA = 40.813_wp
    character(len=64) :: names
    integer :: unit, k, world_rank, ierr
    namelist /parameters/ par

    par0 = par
    ab_default = abundances()

    call set_implicit_items(.true., -1.0_wp)
    open(newunit=unit, file=fname, status='old');  read(unit, parameters);  close(unit)
    par_hi = par
    ab_hi  = abundances()

    par = par0
    call set_implicit_items(.false., -2.0_wp)
    open(newunit=unit, file=fname, status='old');  read(unit, parameters);  close(unit)
    ab_lo = abundances()

    given_um   = par_hi%use_metals .eqv. par%use_metals
    given_ne   = par_hi%metal_ne   .eqv. par%metal_ne
    given_heat = par_hi%metal_heat .eqv. par%metal_heat
    given_ab   = ab_hi == ab_lo
    given_hei  = par_hi%hei_diffuse  .eqv. par%hei_diffuse
    given_heii = par_hi%heii_diffuse .eqv. par%heii_diffuse
    given_fuv  = par_hi%add_fuv      .eqv. par%add_fuv
    given_efuv = par_hi%efuv_min       == par%efuv_min
    given_nfuv = par_hi%no_photons_fuv == par%no_photons_fuv

    !--- every other item is identical in the two reads: keep the first, and
    !--- restore the defaults of the abundances the input left out
    par = par_hi
    call set_abundances(merge(ab_hi, ab_default, given_ab))

    any_positive = any(given_ab .and. ab_hi > 0.0_wp)
    if (.not. given_um)   par%use_metals = any_positive
    if (.not. given_ne)   par%metal_ne   = par%use_metals
    if (.not. given_heat) par%metal_heat = par%use_metals
    !--- the excited-cascade diffuse channels follow the diffuse field
    transported = trim(par%recombination_continuum) == 'transported'
    if (.not. given_hei)  par%hei_diffuse  = transported
    if (.not. given_heii) par%heii_diffuse = transported .and. &
                                             par%eion_max >= E_HEII_LYA
    !--- the FUV extension follows the He I 2^3S population: Phi_3 is scored
    !--- from the transported band and its threshold is 4.78 eV
    !--- no helium, no 2^3S population (solve_ion_cell and thermal_mod run the
    !--- coupling only for abund_He > 0), so no band for it either
    uses_23s = (par%hei_metastable .or. par%hei_metastable_coupling) .and. &
               par%abund_He > 0.0_wp
    if (.not. given_fuv)  par%add_fuv = uses_23s
    if (.not. given_efuv) par%efuv_min = merge(E_HEI_23S_SIGMA_EDGE, &
                                               par0%efuv_min, uses_23s)
    !--- resolved once the band luminosities are known (ion_band_mod):
    !--- no_photons/4, or no_photons when no source has ionizing luminosity
    if (.not. given_nfuv) par%no_photons_fuv = -1.0_wp

    !--- the rank communicators are built after this read; ask MPI directly
    call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, ierr)
    if (world_rank == 0 .and. .not. (given_um .and. given_ne .and. given_heat)) then
       names = ''
       do k = 1, NAB
          if (given_ab(k) .and. ab_hi(k) > 0.0_wp) names = trim(names)//' '//trim(ABNAME(k))
       end do
       write(*,'(a,3(a,l1))') ' METALS: set from the input (items it leaves out):', &
          ' use_metals = ', par%use_metals, ', metal_ne = ', par%metal_ne, &
          ', metal_heat = ', par%metal_heat
       if (len_trim(names) > 0) write(*,'(2a)') &
          ' METALS: abundances the input gives > 0:', trim(names)
    end if
    if (world_rank == 0 .and. transported .and. &
        .not. (given_hei .and. given_heii)) then
       write(*,'(a,2(a,l1))') ' DIFFUSE: excited-cascade channels on with ', &
          'the transported continuum (items the input leaves out): hei_diffuse = ', &
          par%hei_diffuse, ', heii_diffuse = ', par%heii_diffuse
       if (.not. given_heii .and. .not. par%heii_diffuse) write(*,'(a,f7.3,a)') &
          ' DIFFUSE: heii_diffuse left off: the band stops at eion_max =', &
          par%eion_max, ' eV, below the 40.813 eV He II Ly-alpha, so no He III'
    end if
    if (world_rank == 0 .and. uses_23s .and. .not. given_fuv) &
       write(*,'(a)') ' FUV: add_fuv = T for the He I 2^3S population '// &
          '(the input leaves it out)'
    if (world_rank == 0 .and. uses_23s .and. par%add_fuv .and. .not. given_efuv) &
       write(*,'(a,f7.3,a)') ' FUV: efuv_min = ', par%efuv_min, &
          ' eV, where the 2^3S photoionization cross section starts'// &
          ' (threshold 4.768 eV; the input leaves it out)'

  contains

    !--- flag for the logical items, sentinel for the real ones (a value no
    !--- input gives: every real item here is > 0 when set).
    subroutine set_implicit_items(flag, sentinel)
      logical,       intent(in) :: flag
      real(kind=wp), intent(in) :: sentinel
      par%use_metals = flag;  par%metal_ne = flag;  par%metal_heat = flag
      par%hei_diffuse = flag;  par%heii_diffuse = flag
      par%add_fuv = flag
      par%efuv_min = sentinel;  par%no_photons_fuv = sentinel
      call set_abundances([(sentinel, k = 1, NAB)])
    end subroutine set_implicit_items

    subroutine set_abundances(a)
      real(kind=wp), intent(in) :: a(NAB)
      par%abund_C  = a(1);  par%abund_N  = a(2);  par%abund_O  = a(3)
      par%abund_Ne = a(4);  par%abund_S  = a(5);  par%abund_Ar = a(6)
      par%abund_Mg = a(7);  par%abund_Fe = a(8);  par%abund_Si = a(9)
      par%abund_Cl = a(10); par%abund_Ca = a(11)
    end subroutine set_abundances

    function abundances() result(a)
      real(kind=wp) :: a(NAB)
      a = [par%abund_C, par%abund_N, par%abund_O, par%abund_Ne, par%abund_S, &
           par%abund_Ar, par%abund_Mg, par%abund_Fe, par%abund_Si, &
           par%abund_Cl, par%abund_Ca]
    end function abundances

  end subroutine read_parameters_with_implicit_defaults

  !---------------------------------------
  subroutine setup_procedure
  use define
  use random, only : init_random_seed
  use qmc_mod, only : qmc_setup
  use raytrace_amr_mod
  implicit none
  !--- Initialize Random Number Generator
  call init_random_seed(par%iseed)
  !--- Quasi-random launch generator (used only when launch_sequence='sobol').
  call qmc_setup(par%qmc_seed)
  !--- AMR raytrace bindings (the ionizing edge walk is called directly).
  raytrace_to_tau  => raytrace_to_tau_amr
  raytrace_to_edge => raytrace_to_edge_amr
  end subroutine setup_procedure
end module setup_mod
