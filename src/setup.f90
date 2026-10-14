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
  integer :: ierr, islash

  namelist /parameters/ par

  if (command_argument_count() >= 1) then
     call get_command_argument(1, model_infile)
  else
     call get_command_argument(0, arg)
     call refuse_input('no input file.  Usage: '//trim(arg)//' input_file')
  endif

  par%require_convergence = .false.   ! definite baseline; namelist may override
  call read_parameters_with_implicit_defaults(trim(model_infile))
  call refuse_unused_mocafe_items()

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
     call refuse_input( &
        'MoCHII requires par%grid_type = ''amr'' or ''car''.')
  endif
  !--- par%density_file: a uniform 3D density cube (FITS/HDF5) for the 'car'
  !--- grid; nx/ny/nz come from the file, so it is an alternative to both the
  !--- namelist nx/ny/nz and the amr leaf-list.
  if (len_trim(par%density_file) > 0 .and. trim(par%grid_type) /= 'car') then
     call refuse_input( &
        'par%density_file (a 3D density cube) requires par%grid_type = ''car''.')
  endif
  if (len_trim(par%density_file) > 0 .and. len_trim(par%amr_file) > 0) then
     call refuse_input( &
        'set either par%amr_file (leaf list) or par%density_file (cube), not both.')
  endif
  if (len_trim(par%amr_file) == 0) then
     if (trim(par%grid_type) == 'amr') then
        call refuse_input( &
           'grid_type=''amr'' requires par%amr_file (refinement structure).')
     else if (len_trim(par%density_file) > 0) then
        if (par%xmax <= 0.0_wp) then
           call refuse_input( &
              'grid_type=''car'' with par%density_file needs par%xmax > 0.')
        end if
     else if (par%xmax <= 0.0_wp .or. &
              (par%nx < 2 .and. .not. par%xy_periodic)) then
        !--- nx=1 is allowed only with xy_periodic (a plane-parallel column
        !--- whose single lateral cell wraps onto itself).
        call refuse_input( &
           'grid_type=''car'' without par%amr_file needs par%xmax>0 and '// &
           'par%nx>=2 (par%nx>=1 allowed with par%xy_periodic).')
     else if (par%ny < 1 .or. par%nz < 1 .or. &
              par%ymax <= 0.0_wp .or. par%zmax <= 0.0_wp) then
        !--- the remaining namelist-built dimensions must be positive.
        call refuse_input( &
           'grid_type=''car'' without par%amr_file needs par%ny>=1, '// &
           'par%nz>=1, par%ymax>0, and par%zmax>0.')
     end if
  endif
  if (trim(par%grid_type) == 'car' .and. par%refine_front) then
     call refuse_input( &
        'refine_front needs the octree (par%grid_type = ''amr'').')
  endif
  if (par%refine_front) then
     if (par%refine_lbase > par%refine_lmax) then
        call refuse_input('par%refine_lbase ('//integer_text(par%refine_lbase)// &
           ') must not exceed par%refine_lmax ('//integer_text(par%refine_lmax)//').')
     endif
     if (par%refine_lbase < 1 .or. par%refine_lmax < 1) then
        call refuse_input( &
           'par%refine_lbase and par%refine_lmax must be >= 1.')
     endif
     if (par%refine_iter < 1 .or. par%refine_iter >= par%gas_niter) then
        call refuse_input('par%refine_iter ('//integer_text(par%refine_iter)// &
           ') must be in [1, par%gas_niter-1] so the new grid is iterated '// &
           'on; par%gas_niter = '//integer_text(par%gas_niter)//'.')
     endif
  endif
  !--- A string item with a fixed set of values is checked against that set
  !--- whether or not this run reads it, so a misspelled value is refused
  !--- rather than falling through to a default branch or going unread.
  if (trim(par%car_walk) /= 'dda' .and. trim(par%car_walk) /= 'amr') then
     call refuse_input( &
        'par%car_walk must be ''dda'' (default) or ''amr''.')
  endif
  if (par%restart_iteration_history .and. len_trim(par%restart_file) == 0) &
     call refuse_input('par%restart_iteration_history needs par%restart_file.')
  if (.not. par%conv_drift_test .and. trim(par%conv_crit) /= 'vol' &
      .and. mpar%p_rank == 0) write(*,'(a)') &
     ' NOTE: par%conv_drift_test is read only with par%conv_crit = ''vol''.'
  if (trim(par%conv_crit) /= 'cell' .and. trim(par%conv_crit) /= 'vol') then
     call refuse_input('par%conv_crit must be ''cell'' (default; the cell '// &
        'maxima) or ''vol'' (the volume-integrated drift test).')
  endif
  !--- iofile_mod::io_file_extension maps 'hdf5'/'h5' to HDF5 and anything
  !--- else to FITS, so a misspelling would silently write FITS.
  if (strlowcase(trim(par%file_format)) /= 'hdf5' .and. &
      strlowcase(trim(par%file_format)) /= 'h5'   .and. &
      strlowcase(trim(par%file_format)) /= 'fits') then
     call refuse_input('par%file_format must be ''hdf5'' (default) or ''fits''.')
  endif
  !--- grid_mod_amr's leaf dust opacity: an unknown law used to fall through
  !--- to 'global_dgr'.
  select case (trim(par%dust_density_law))
  case ('global_dgr', 'none', 'laursen09_live', 'laursen09', 'from_file')
  case default
     call refuse_input('par%dust_density_law must be ''global_dgr'' (default), '// &
        '''none'', ''laursen09_live'', ''laursen09'' or ''from_file''.')
  end select
  if (trim(par%amr_type) == 'ramses') then
     call refuse_input('par%amr_type = ''ramses'' is not read directly; convert '// &
        'the snapshot to the generic leaf file first (MoCafe v2.00 carries '// &
        'python/AMR_grid/convert_ramses_to_generic.py; this tree does not).')
  else if (trim(par%amr_type) /= 'generic') then
     call refuse_input('par%amr_type must be ''generic'' (the leaf-list file).')
  endif
  if (trim(par%ext_geometry) /= 'rec' .and. trim(par%ext_geometry) /= 'sph') then
     call refuse_input('par%ext_geometry must be ''rec'' (default) or ''sph''.')
  endif
  if (trim(par%slab_faces) /= 'top' .and. trim(par%slab_faces) /= 'bottom' &
      .and. trim(par%slab_faces) /= 'both') then
     call refuse_input( &
        'par%slab_faces must be ''top'', ''bottom'', or ''both''.')
  endif
  if ((trim(par%slab_top_mode) /= 'beam' .and. trim(par%slab_top_mode) /= 'isotropic') .or. &
      (trim(par%slab_bot_mode) /= 'beam' .and. trim(par%slab_bot_mode) /= 'isotropic')) then
     call refuse_input( &
        'par%slab_top_mode / slab_bot_mode must be ''beam'' or ''isotropic''.')
  endif
  !--- the ionizing band launches point sources only (every slot, whether or
  !--- not par%nsource reaches it: a slot the run does not use is still a
  !--- value the input asked for)
  if (any(par%src_geometry /= 'point')) then
     call refuse_input('the ionizing band supports only par%src_geometry = ''point''.')
  endif
  if (trim(par%ci_model) /= 'voronov' .and. &
      trim(par%ci_model) /= 'dere_hybrid') then
     call refuse_input( &
        'par%ci_model must be ''voronov'' (default) or ''dere_hybrid''.')
  endif
  if (trim(par%charge_exchange) /= 'two_way' .and. &
      trim(par%charge_exchange) /= 'metal_only') then
     call refuse_input( &
        'par%charge_exchange must be ''two_way'' (default) or '// &
        '''metal_only'' (reproduction only).')
  endif
  if (trim(par%recomb_model) /= 'detailed_balance' .and. &
      trim(par%recomb_model) /= 'hui_gnedin') then
     if (trim(par%recomb_model) == 'badnell_milne') call refuse_input( &
        'par%recomb_model = ''badnell_milne'' was renamed ''detailed_balance'' '// &
        'on 2026-09-26 (the hydrogenic ions no longer use the Badnell total).')
     call refuse_input('par%recomb_model must be ''detailed_balance'' '// &
        '(default) or ''hui_gnedin''.')
  endif
  if (trim(par%metal_excitation_model) == 'extended') call refuse_input( &
     'par%metal_excitation_model = ''extended'' was removed 2026-10-11: '// &
     'use ''chianti'' (full CHIANTI 11.0.2 atoms, the default).')
  if (trim(par%metal_excitation_model) /= 'chianti' .and. &
      trim(par%metal_excitation_model) /= 'mocassin_v2') then
     call refuse_input('par%metal_excitation_model must be ''chianti'' '// &
        '(default; full CHIANTI 11.0.2 atoms) or ''mocassin_v2'' (MOCASSIN '// &
        'v2 atoms, for parity runs).')
  endif
  if (trim(par%diffuse_energy_model) /= 'exponential' .and. &
      trim(par%diffuse_energy_model) /= 'threshold'   .and. &
      trim(par%diffuse_energy_model) /= 'milne') then
     call refuse_input( &
          'par%diffuse_energy_model must be ''milne'', ' // &
          '''exponential'' or ''threshold''.')
  endif
  if (trim(par%cooling_model) /= 'low_density' .and. &
      trim(par%cooling_model) /= 'local_ne') then
     call refuse_input( &
        'par%cooling_model must be ''local_ne'' (default) or ''low_density''.')
  endif
  if (trim(par%twophoton_quenching) /= 'colliders' .and. &
      trim(par%twophoton_quenching) /= 'lumped_ne') then
     call refuse_input( &
        'par%twophoton_quenching must be ''colliders'' (default) or ''lumped_ne''.')
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
       call refuse_input( &
          'par%spectrum_type must be ''shape_ev'' (default), ''per_ev'', ''per_hz'', '// &
          '''per_ang'', or ''per_um''.')
    end select

    !--- the ISRF presets (draine/habing/mathis) are external-only; an internal
    !--- file slot naming one is almost certainly a mistake.
    if (ion_ext_preset_id(par%ion_spectrum) > 0 .or. &
        ion_ext_preset_id(par%src_spectrum_file) > 0) then
       call refuse_input( &
          'ISRF preset names (draine/habing/mathis) are valid only for '// &
          'par%ext_spectrum, not for an internal source file slot.')
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
          call refuse_input( &
             'external source_geometry requires par%ext_intensity > 0 '// &
             '(mean intensity J), an ISRF preset, or a physical-type ext_spectrum.')
       endif
       par%nsource = 0
       if (mpar%p_rank == 0) write(*,'(a)') &
          'NOTE: source_geometry=''external*'' is the external-only shorthand; '// &
          'the composable form is par%ext_intensity/par%ext_geometry + par%nsource.'
    else if (trim(par%source_geometry) == 'slab') then
       !--- plane-parallel slab illumination (its own single-component source).
       if (.not. par%xy_periodic) then
          call refuse_input( &
             'source_geometry=''slab'' needs par%xy_periodic (a slab is '// &
             'periodic in x and y).')
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
            call refuse_input( &
               'par%slab_top_theta must be in [0, 90) degrees for a beam.')
         endif
         if (lbot .and. trim(par%slab_bot_mode) == 'beam' .and. &
             (par%slab_bot_theta < 0.0_wp .or. par%slab_bot_theta >= 90.0_wp)) then
            call refuse_input( &
               'par%slab_bot_theta must be in [0, 90) degrees for a beam.')
         endif
         !--- each lit face must carry a positive source strength on its own (a
         !--- negative one could hide behind a positive one in the sum and push
         !--- slab_ptop outside [0, 1]).
         if (ltop .and. par%slab_top_source <= 0.0_wp) then
            call refuse_input( &
               'par%slab_top_source must be > 0 on the lit top face.')
         endif
         if (lbot .and. par%slab_bot_source <= 0.0_wp) then
            call refuse_input( &
               'par%slab_bot_source must be > 0 on the lit bottom face.')
         endif
       end block
       if (par%slab_nmu < 1) then
          call refuse_input( &
             'par%slab_nmu must be >= 1.')
       endif
       !--- slab illumination is its own single source; composing it with
       !--- internal point sources or an external field is not yet supported.
       !--- nsource defaults to 1, so only nsource > 1 (or an external field)
       !--- signals an explicit extra-source request; nsource <= 1 is silently
       !--- overridden by the slab below.
       if (par%nsource > 1 .or. par%ext_intensity > 0.0_wp) then
          call refuse_input( &
             'source_geometry=''slab'' cannot be combined with internal '// &
             'point sources (par%nsource > 1) or an external field '// &
             '(par%ext_intensity > 0); slab+other-source composition is not yet supported.')
       endif
       par%nsource = 0
    else if (trim(par%source_geometry) /= 'point') then
       call refuse_input( &
          'par%source_geometry must be ''point'', ''slab'', ''external_rec'', or ''external_sph''.')
    endif

    !--- external field ON: mean intensity, an ISRF preset, or a physical-type
    !--- ext_spectrum file (the last two are absolute and need no ext_intensity).
    ext_on = (par%ext_intensity > 0.0_wp) .or. ext_is_preset .or. ext_is_physfile

    !--- ISRF presets are analytic and FUV-only: require add_fuv, a positive
    !--- scale, and note that ext_intensity (if set) overrides ext_scale.
    if (ext_is_preset) then
       if (par%ext_scale <= 0.0_wp) then
          call refuse_input( &
             'par%ext_scale must be > 0 for an ISRF preset.')
       endif
       if (.not. par%add_fuv) then
          call refuse_input( &
             'the ISRF presets are FUV-only; set par%add_fuv (with '// &
             'par%efuv_min < 13.6) so the FUV bins carry the field.')
       endif
       if (par%ext_intensity > 0.0_wp .and. par%ext_scale /= 1.0_wp &
           .and. mpar%p_rank == 0) write(*,'(a)') &
          'NOTE: par%ext_intensity rescales the preset to a target band J; '// &
          'par%ext_scale is then ignored.'
    endif
    !--- the spherical entry surface needs its radius (only when the external
    !--- field is on; the value set of ext_geometry is checked above).
    if (ext_on) then
       if (trim(par%ext_geometry) == 'sph' .and. par%rmax <= 0.0_wp) then
          call refuse_input( &
             'ext_geometry=''sph'' requires par%rmax > 0.')
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
       call refuse_input('par%nsource must be in [0, '//integer_text(MAX_SRC)//'].')
    endif
    if (par%nsource == 0 .and. trim(par%source_geometry) /= 'slab') then
       if (.not. ext_on) then
          call refuse_input( &
             'no source: set par%nsource >= 1 or the external field '// &
             '(par%ext_intensity > 0, an ISRF preset, or a physical ext_spectrum).')
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
    else
       !--- multiple point sources.  src_lum: all set (keep) or all unset.  For
       !--- a 'shape_ev' spectrum an unset set is the equal split of par%luminosity
       !--- (MoCafe convention); for a physical-type (absolute) source file an
       !--- unset set is DERIVED per column (src_lum left as the -999 sentinel).
       !--- A partial set is an error.
       nset = 0
       do is = 1, par%nsource
          if (par%src_lum(is) > 0.0_wp) nset = nset + 1
       end do
       if (nset == 0) then
          if (.not. src_is_absolute) then
             do is = 1, par%nsource
                par%src_lum(is) = par%luminosity / real(par%nsource, wp)
             end do
          endif
       else if (nset /= par%nsource) then
          call refuse_input( &
             'set par%src_lum for ALL sources or for NONE '// &
             '(equal split, or derived for a physical-type source file).')
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
       call refuse_input( &
          'par%launch_sequence must be ''sobol'' (default) or ''random''.')
    endif
  end block

  !--- general parameter sanity (fail fast rather than run on nonsense).
  if (par%no_photons < 1.0_wp) then
     call refuse_input( &
        'par%no_photons must be >= 1.')
  endif
  if (trim(par%fuv_sampling) /= 'separate' .and. &
      trim(par%fuv_sampling) /= 'proportional') then
     call refuse_input( &
        'par%fuv_sampling must be ''separate'' (default) or ''proportional''.')
  endif
  if (par%add_fuv .and. trim(par%fuv_sampling) == 'separate' .and. &
      par%no_photons_fuv < 1.0_wp .and. par%no_photons_fuv /= -1.0_wp) then
     call refuse_input( &
        'par%no_photons_fuv must be >= 1 with par%add_fuv and '// &
        'par%fuv_sampling = ''separate''.')
  endif
  if (par%source_cdf_tol <= 0.0_wp .or. par%source_cdf_tol >= 1.0_wp) then
     call refuse_input( &
        'par%source_cdf_tol must be in (0, 1).')
  endif
  !--- Energy sampling covers every source configuration (point sources,
  !--- external field, slab), so the one requirement is that SOME spectrum is
  !--- defined for the components to sample.
  if (len_trim(par%ion_spectrum) == 0 .and. &
      len_trim(par%src_spectrum_file) == 0 .and. &
      len_trim(par%ext_spectrum) == 0 .and. &
      par%tstar <= 0.0_wp .and. par%ext_tstar <= 0.0_wp .and. &
      all(par%src_tstar <= 0.0_wp)) then
     call refuse_input( &
        'the ionizing band needs a Planck or tabulated source spectrum.')
  endif
  if (par%ion_relax <= 0.0_wp .or. par%ion_relax > 1.0_wp) then
     call refuse_input( &
        'par%ion_relax must be in (0, 1].')
  endif
  if (trim(par%lyman_lines) /= 'on_the_spot' .and. &
      trim(par%lyman_lines) /= 'escaping') then
     call refuse_input( &
        'par%lyman_lines must be ''on_the_spot'' (default) or ' // &
        '''escaping'' (MoCHII does not transport the Lyman lines).')
  endif
  if (trim(par%recombination_continuum) /= 'on_the_spot' .and. &
      trim(par%recombination_continuum) /= 'transported' .and. &
      trim(par%recombination_continuum) /= 'escaping') then
     call refuse_input( &
        'par%recombination_continuum must be ''on_the_spot'', ' // &
        '''transported'' (default) or ''escaping''.')
  endif
  if (par%solve_te .and. par%te_min >= par%te_max) then
     call refuse_input( &
        'par%solve_te needs par%te_min < par%te_max.')
  endif
  if (par%ion_peel .and. (par%nxim < 1 .or. par%nyim < 1)) then
     call refuse_input( &
        'par%ion_peel needs par%nxim >= 1 and par%nyim >= 1.')
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
     call refuse_input( &
        'par%ion_peel is not supported with par%xy_periodic (a periodic '// &
        'slab has no finite-distance observer; use the boundary I(mu) output '// &
        'written for source_geometry=''slab'').')
  endif

  !--- ionizing band
  if (par%eion_max <= par%eion_min .or. par%nnu_ion < 1) then
     call refuse_input( &
        'need par%eion_max > par%eion_min and par%nnu_ion >= 1.')
  endif
  if (len_trim(par%ion_spectrum) == 0 .and. par%tstar <= 0.0_wp .and. &
      len_trim(par%src_spectrum_file) == 0 .and. &
      len_trim(par%ext_spectrum) == 0 .and. par%ext_tstar <= 0.0_wp) then
     call refuse_input( &
        'the ionizing band needs a spectrum: par%tstar > 0, par%ion_spectrum, '// &
        'par%src_spectrum_file, par%ext_spectrum, or par%ext_tstar > 0.')
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
     call refuse_input( &
        'par%grain_pe needs par%add_fuv (the FUV field is its heating source).')
  endif
  if (par%grain_pe .and. .not. par%ion_add_dust) then
     call refuse_input( &
        'par%grain_pe needs par%ion_add_dust (grain opacity carries the field).')
  endif
  if (par%gas_grain_energy_exchange .and. .not. par%ion_add_dust) then
     call refuse_input('par%gas_grain_energy_exchange needs par%ion_add_dust '// &
        '(the grains and their temperature come from the dust model).')
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
        call refuse_input( &
           'unknown par%distance_unit ('''//trim(par%distance_unit)// &
           ''') — use ''kpc'', ''pc'', ''au'', or '''' (cm).')
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
  !   add_fuv               not in the input: .true. when something absorbs
  !                         below 13.6 eV -- the He I 2^3S population
  !                         (hei_metastable or hei_metastable_coupling,
  !                         abund_He > 0; Phi_3 is a path estimator over the
  !                         transported band and starts at 4.78 eV), a metal
  !                         with a first ionization potential below the H I
  !                         edge, or dust -- else .false.;
  !   efuv_min              not in the input: E_HEI_23S_SIGMA_EDGE = 4.78 eV
  !                         under the same condition, else the define.f90
  !                         default;
  !   grain_pe              not in the input: ion_add_dust .and. add_fuv
  !                         (WD01 photoelectric heating and accretion
  !                         cooling, grain_photoelectric_mod);
  !   gas_grain_energy_exchange  not in the input: ion_add_dust;
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
    use metal_ionization_potential_mod, only : metal_ionization_potential
    implicit none
    character(len=*), intent(in) :: fname
    integer, parameter :: NAB = 11
    character(len=2), parameter :: ABNAME(NAB) = &
       ['C ', 'N ', 'O ', 'Ne', 'S ', 'Ar', 'Mg', 'Fe', 'Si', 'Cl', 'Ca']
    integer, parameter :: ABZ(NAB) = [6, 7, 8, 10, 16, 18, 12, 26, 14, 17, 20]
    type(params_type) :: par0, par_hi
    real(kind=wp) :: ab_hi(NAB), ab_lo(NAB), ab_default(NAB), ab_now(NAB)
    logical :: given_um, given_ne, given_heat, given_ab(NAB), any_positive
    logical :: given_hei, given_heii, transported
    logical :: given_fuv, given_efuv, given_nfuv, uses_23s
    logical :: given_pe, given_gge
    logical :: metal_below_hi_edge, fuv_absorbers
    character(len=128) :: fuv_reasons
    !--- the He II Ly-alpha energy [eV] (heii_cascade_mod::E_HEII_LYA)
    real(kind=wp), parameter :: E_HEII_LYA = 40.813_wp
    character(len=64) :: names
    integer :: k, world_rank, ierr
    namelist /parameters/ par

    par0 = par
    ab_default = abundances()

    call set_implicit_items(.true., -1.0_wp)
    call read_namelist_file()
    par_hi = par
    ab_hi  = abundances()

    par = par0
    call set_implicit_items(.false., -2.0_wp)
    call read_namelist_file()
    ab_lo = abundances()

    given_um   = par_hi%use_metals .eqv. par%use_metals
    given_ne   = par_hi%metal_ne   .eqv. par%metal_ne
    given_heat = par_hi%metal_heat .eqv. par%metal_heat
    given_ab   = ab_hi == ab_lo
    given_hei  = par_hi%hei_diffuse  .eqv. par%hei_diffuse
    given_heii = par_hi%heii_diffuse .eqv. par%heii_diffuse
    given_fuv  = par_hi%add_fuv      .eqv. par%add_fuv
    given_pe   = par_hi%grain_pe     .eqv. par%grain_pe
    given_gge  = par_hi%gas_grain_energy_exchange .eqv. &
                 par%gas_grain_energy_exchange
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
    !--- the FUV extension is on whenever something in the run absorbs below
    !--- the H I edge: the He I 2^3S population (Phi_3, threshold 4.78 eV;
    !--- no helium, no 2^3S population -- solve_ion_cell and thermal_mod run
    !--- the coupling only for abund_He > 0), a metal with its neutral
    !--- ionization potential below 13.598 eV, or dust
    uses_23s = (par%hei_metastable .or. par%hei_metastable_coupling) .and. &
               par%abund_He > 0.0_wp
    metal_below_hi_edge = .false.
    ab_now = abundances()
    if (par%use_metals) then
       do k = 1, NAB
          if (ab_now(k) > 0.0_wp .and. &
              metal_ionization_potential(ABZ(k), 1) < eth_HI) &
             metal_below_hi_edge = .true.
       end do
    end if
    fuv_absorbers = uses_23s .or. metal_below_hi_edge .or. par%ion_add_dust
    fuv_reasons = ''
    if (uses_23s) fuv_reasons = trim(fuv_reasons)//' He I 2^3S,'
    if (metal_below_hi_edge) fuv_reasons = trim(fuv_reasons)//' metals,'
    if (par%ion_add_dust) fuv_reasons = trim(fuv_reasons)//' dust,'
    if (len_trim(fuv_reasons) > 0) fuv_reasons = fuv_reasons(1:len_trim(fuv_reasons)-1)
    if (.not. given_fuv)  par%add_fuv = fuv_absorbers
    !--- the gas-grain physics is on wherever there are grains: photoelectric
    !--- emission and electron and ion capture need the FUV band as well
    if (.not. given_pe)   par%grain_pe = par%ion_add_dust .and. par%add_fuv
    if (.not. given_gge)  par%gas_grain_energy_exchange = par%ion_add_dust
    if (.not. given_efuv) par%efuv_min = merge(E_HEI_23S_SIGMA_EDGE, &
                                               par0%efuv_min, fuv_absorbers)
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
    if (world_rank == 0 .and. par%ion_add_dust .and. &
        .not. (given_pe .and. given_gge)) &
       write(*,'(2(a,l1),a)') ' GRAIN: set with the dust (items the input '// &
          'leaves out): grain_pe = ', par%grain_pe, &
          ', gas_grain_energy_exchange = ', par%gas_grain_energy_exchange, &
          merge(' (grain_pe needs add_fuv)', '                         ', &
                .not. par%add_fuv .and. .not. given_pe)
    if (world_rank == 0 .and. fuv_absorbers .and. .not. given_fuv) &
       write(*,'(3a)') ' FUV: add_fuv = T for the absorbers below 13.6 eV:', &
          trim(fuv_reasons), ' (the input leaves it out)'
    if (world_rank == 0 .and. fuv_absorbers .and. par%add_fuv .and. .not. given_efuv) &
       write(*,'(a,f7.3,a)') ' FUV: efuv_min = ', par%efuv_min, &
          ' eV, where the 2^3S photoionization cross section starts,'// &
          ' below every metal edge (the input leaves it out)'

  contains

    !--- One read of the namelist group /parameters/.  An unknown or
    !--- misspelled item, a value that does not fit its item, or a file with
    !--- no &parameters group ends the run with the compiler's message; Intel
    !--- Fortran's names the line where it stopped reading ("line N"), which
    !--- can be the line after the offending one, so both are echoed.
    subroutine read_namelist_file()
      use utility, only : refuse_input
      integer :: unit, ios, ios_open
      character(len=1024) :: msg
      open(newunit=unit, file=fname, status='old', action='read', &
           iostat=ios_open, iomsg=msg)
      if (ios_open /= 0) call refuse_input('cannot open the input file '''// &
         fname//''': '//trim(msg))
      close(unit)
      call refuse_removed_items()
      open(newunit=unit, file=fname, status='old', action='read')
      read(unit, nml=parameters, iostat=ios, iomsg=msg)
      close(unit)
      if (ios < 0) call refuse_input('no namelist group &parameters ... / '// &
         'was found in '''//fname//''' ('//trim(msg)//').')
      if (ios > 0) call refuse_input('the namelist group &parameters in '''// &
         fname//''' could not be read (an unknown or misspelled item, or a '// &
         'value that does not fit its item):'//new_line('a')//'       '// &
         trim(msg)//input_lines_near(msg))
    end subroutine read_namelist_file

    !--- Items removed or renamed are refused with the name of what replaced
    !--- them (the namelist read alone would only say "unknown item").
    subroutine refuse_removed_items()
      use utility, only : refuse_input
      integer, parameter :: NREMOVED = 8
      character(len=32), parameter :: removed(NREMOVED) = [character(len=32) :: &
         'metal_inner_shell_absorption', 'additional_line_cooling', &
         'high_charge_exchange', 'fe_levels_full', 'case_ab', &
         'diffuse_field', 'line_case', 'pe_scale']
      character(len=200), parameter :: replacement(NREMOVED) = [character(len=200) :: &
         'removed 2026-10-11: metal absorption always sums every subshell '// &
            '(phfit2); multi-stage Auger products are the planned '// &
            'par%auger_ionization', &
         'removed 2026-10-11: every stage with an atom in '// &
            'par%metal_excitation_model contributes cooling and lines', &
         'renamed 2026-10-11: use par%langevin_charge_exchange', &
         'removed 2026-10-11: the Fe II/III atoms come from '// &
            'par%metal_excitation_model', &
         'removed 2026-09-25: use par%recombination_continuum '// &
            '(''transported'', ''on_the_spot'' or ''escaping'')', &
         'removed 2026-09-25: use par%recombination_continuum '// &
            '(''transported'', ''on_the_spot'' or ''escaping'')', &
         'removed 2026-09-25: use par%lyman_lines '// &
            '(''on_the_spot'' or ''escaping'')', &
         'removed 2026-10-11 with the Bakes & Tielens (1994) fits it '// &
            'scaled: par%grain_pe (Weingartner & Draine 2001) takes the '// &
            'grains of par%dust_model at the leaf''s dust opacity']
      character(len=1024) :: line, packed
      integer :: unit, ios, i, j, at, nline, n
      open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      nline = 0
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         nline = nline + 1
         at = index(line, '!')
         if (at > 0) line(at:) = ' '
         !--- lower case, blanks removed
         packed = ''
         n = 0
         do j = 1, len_trim(line)
            if (line(j:j) == ' ' .or. line(j:j) == char(9)) cycle
            n = n + 1
            packed(n:n) = line(j:j)
            if (packed(n:n) >= 'A' .and. packed(n:n) <= 'Z') &
               packed(n:n) = achar(iachar(packed(n:n)) + 32)
         end do
         do i = 1, NREMOVED
            at = index(packed, 'par%'//trim(removed(i)))
            if (at == 0) cycle
            j = at + 4 + len_trim(removed(i))
            if (j > n) cycle
            if (packed(j:j) /= '=' .and. packed(j:j) /= '(') cycle
            close(unit)
            call refuse_input('par%'//trim(removed(i))//' (line '// &
               integer_text(nline)//' of '''//fname//''') was '// &
               trim(replacement(i))//'.')
         end do
      end do
      close(unit)
    end subroutine refuse_removed_items

    !--- The line the read message names and the one before it, numbered,
    !--- or nothing when the message names no line.
    function input_lines_near(msg) result(text)
      character(len=*), intent(in)  :: msg
      character(len=:), allocatable :: text
      character(len=1024) :: line
      integer :: at, nline, unit, ios, i
      text = ''
      at = index(msg, 'line ', back=.true.)
      if (at == 0) return
      read(msg(at+5:), *, iostat=ios) nline
      if (ios /= 0 .or. nline < 1) return
      open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do i = 1, nline
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (i >= nline - 1) text = text//new_line('a')//'       line '// &
            integer_text(i)//': '//trim(line)
      end do
      close(unit)
    end function input_lines_near

    !--- flag for the logical items, sentinel for the real ones (a value no
    !--- input gives: every real item here is > 0 when set).
    subroutine set_implicit_items(flag, sentinel)
      logical,       intent(in) :: flag
      real(kind=wp), intent(in) :: sentinel
      par%use_metals = flag;  par%metal_ne = flag;  par%metal_heat = flag
      par%hei_diffuse = flag;  par%heii_diffuse = flag
      par%add_fuv = flag
      par%grain_pe = flag;  par%gas_grain_energy_exchange = flag
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

  !=========================================================================
  ! params_type (define.f90) still declares items of the MoCafe dust-
  ! continuum code that MoCHII has no code for.  The namelist accepts them,
  ! so a MoCafe input still parses, but an item set in the input would do
  ! nothing and say nothing; this refuses every one that differs from its
  ! declared default, by name.  The list is checked by grep over src/:
  !   referenced nowhere outside define.f90 (33): the first block below;
  !   referenced only where it has no effect (19):
  !     grid_mod_car::grid_create builds a box grid that grid_create_amr
  !       zeroes and whose taumax/tauhomo it overwrites -- geometry,
  !       density_zscale, density_rscale, density_powerlaw_index,
  !       density_zprofile, spiral_m, spiral_pitch, spiral_amp, and use_sed
  !       (whose distance-unit change grid_create_amr undoes);
  !     grid_create_amr sets xyz_symmetry and z_symmetry .false. before any
  !       transport;
  !     sed_mod::setup_sed, which nothing calls, is the only reader of
  !       nlambda, lambda_min, lambda_max, lambda0, kext_file,
  !       source_spectrum, hgg and albedo.
  ! docs/MoCHII_UserGuide.tex and md/MOCAFE_PARAMETER_DIVERGENCE.md carry
  ! the same list.  Porting one of these features means deleting its line.
  !=========================================================================
  subroutine refuse_unused_mocafe_items()
    use define
    use utility, only : refuse_input
    use, intrinsic :: ieee_arithmetic, only : ieee_class, ieee_quiet_nan, &
                                              operator(==), operator(/=)
    implicit none
    type(params_type) :: d          ! every item at its declared default
    character(len=:), allocatable :: report
    integer :: nset

    report = '';  nset = 0
    !--- referenced nowhere outside define.f90
    call item(par%sightline_tau .neqv. d%sightline_tau, 'sightline_tau')
    call item(par%nr /= d%nr, 'nr')
    call item(differs(par%source_zscale, d%source_zscale), 'source_zscale')
    call item(differs(par%source_rscale, d%source_rscale), 'source_rscale')
    call item(par%out_bitpix /= d%out_bitpix, 'out_bitpix')
    call item(par%out_bitpix_force .neqv. d%out_bitpix_force, 'out_bitpix_force')
    call item(par%use_reduced_wgt .neqv. d%use_reduced_wgt, 'use_reduced_wgt')
    call item(par%use_stokes .neqv. d%use_stokes, 'use_stokes')
    call item(par%scatt_mat_file /= d%scatt_mat_file, 'scatt_mat_file')
    call item(par%save_jlam .neqv. d%save_jlam, 'save_jlam')
    call item(par%use_dustemis .neqv. d%use_dustemis, 'use_dustemis', &
         'dust emission is par%dust_sed and par%dust_emis_transport')
    call item(par%dust_emission_method /= d%dust_emission_method, &
         'dust_emission_method', 'dust emission is par%dust_sed')
    call item(par%dust_fast_table .neqv. d%dust_fast_table, 'dust_fast_table', &
         'dust emission is par%dust_sed')
    call item(par%dust_nU /= d%dust_nU, 'dust_nU', 'dust emission is par%dust_sed')
    call item(par%dust_single_teq .neqv. d%dust_single_teq, 'dust_single_teq', &
         'dust emission is par%dust_sed')
    call item(any(par%src_spectrum /= d%src_spectrum), 'src_spectrum', &
         'the source spectra are par%src_spectrum_file, par%ion_spectrum and '// &
         'par%src_tstar')
    call item(any(par%src_zscale /= d%src_zscale), 'src_zscale')
    call item(any(par%src_rscale /= d%src_rscale), 'src_rscale')
    call item(any(par%src_sersic_index /= d%src_sersic_index), 'src_sersic_index')
    call item(any(par%src_reff /= d%src_reff), 'src_reff')
    call item(any(par%src_axial_ratio /= d%src_axial_ratio), 'src_axial_ratio')
    call item(any(par%src_boxiness /= d%src_boxiness), 'src_boxiness')
    call item(par%allsky .neqv. d%allsky, 'allsky')
    call item(par%allsky_nside /= d%allsky_nside, 'allsky_nside')
    call item(par%allsky_x /= d%allsky_x, 'allsky_x')
    call item(par%allsky_y /= d%allsky_y, 'allsky_y')
    call item(par%allsky_z /= d%allsky_z, 'allsky_z')
    call item(par%use_mrw .neqv. d%use_mrw, 'use_mrw')
    call item(par%mrw_gamma /= d%mrw_gamma, 'mrw_gamma')
    call item(par%use_amr_grid .neqv. d%use_amr_grid, 'use_amr_grid', &
         'the octree is par%grid_type = ''amr''')
    call item(par%cone_opening /= d%cone_opening, 'cone_opening')
    call item(par%output_normalization /= d%output_normalization, &
         'output_normalization')
    call item(par%radiation_angular_PDF_file /= d%radiation_angular_PDF_file, &
         'radiation_angular_PDF_file')
    !--- referenced only where it has no effect (see the header)
    call item(par%geometry /= d%geometry, 'geometry', &
         'the density is par%amr_file, par%density_file or par%nH_const with '// &
         'par%rmin/rmax')
    call item(differs(par%density_zscale, d%density_zscale), 'density_zscale')
    call item(differs(par%density_rscale, d%density_rscale), 'density_rscale')
    call item(differs(par%density_powerlaw_index, d%density_powerlaw_index), &
         'density_powerlaw_index')
    call item(par%density_zprofile /= d%density_zprofile, 'density_zprofile')
    call item(par%spiral_m /= d%spiral_m, 'spiral_m')
    call item(par%spiral_pitch /= d%spiral_pitch, 'spiral_pitch')
    call item(par%spiral_amp /= d%spiral_amp, 'spiral_amp')
    call item(par%xyz_symmetry .neqv. d%xyz_symmetry, 'xyz_symmetry')
    call item(par%z_symmetry .neqv. d%z_symmetry, 'z_symmetry')
    call item(par%use_sed .neqv. d%use_sed, 'use_sed')
    call item(par%nlambda /= d%nlambda, 'nlambda')
    call item(par%lambda_min /= d%lambda_min, 'lambda_min')
    call item(par%lambda_max /= d%lambda_max, 'lambda_max')
    call item(par%lambda0 /= d%lambda0, 'lambda0')
    call item(par%kext_file /= d%kext_file, 'kext_file', &
         'the grain optics are par%dust_model, par%sed_kext and par%ion_dust_kext')
    call item(par%source_spectrum /= d%source_spectrum, 'source_spectrum', &
         'the source spectra are par%ion_spectrum, par%src_spectrum_file and '// &
         'par%ext_spectrum')
    call item(par%hgg /= d%hgg, 'hgg')
    call item(par%albedo /= d%albedo, 'albedo')

    if (nset > 0) call refuse_input('the input sets '//integer_text(nset)// &
       ' item(s) that MoCHII does not use; remove them:'//report)

  contains

    subroutine item(set, name, instead)
      logical,          intent(in)           :: set
      character(len=*), intent(in)           :: name
      character(len=*), intent(in), optional :: instead
      if (.not. set) return
      nset = nset + 1
      report = report//new_line('a')//'       par%'//name// &
         ' is a MoCafe dust-continuum parameter that MoCHII does not use'
      if (present(instead)) report = report//' ('//instead//')'
    end subroutine item

    !--- a real item against its default, where the default may be the NaN
    !--- "unset" marker (compared by class: -fp-model fast does not keep
    !--- NaN comparisons)
    logical function differs(a, b)
      real(kind=wp), intent(in) :: a, b
      if (ieee_class(b) == ieee_quiet_nan) then
         differs = ieee_class(a) /= ieee_quiet_nan
      else
         differs = ieee_class(a) == ieee_quiet_nan .or. a /= b
      end if
    end function differs

  end subroutine refuse_unused_mocafe_items

  !--- an integer as its shortest decimal text
  function integer_text(i) result(text)
    integer, intent(in) :: i
    character(len=:), allocatable :: text
    character(len=16) :: buf
    write(buf, '(i0)') i
    text = trim(buf)
  end function integer_text

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
