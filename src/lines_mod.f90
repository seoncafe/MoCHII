module lines_mod
! Line luminosities and optically thin emissivities on the written gas state.
! The selected radiative branches use populations from the complete atom.
  use define
  use octree_mod, only : leaf_half, leaf_cx, leaf_cy, leaf_cz
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHII, gas_xHeI, gas_xHeII, &
                            gas_xHeIII, gas_ne, gas_Te, gas_nleaf, &
                            gas_leaf_owner_first, gas_leaf_owner_stride
  use mpi
  use species_mod, only : n_elements, elem_name, elem_nstage, elem_abund, &
                          species_stage_fraction, species_write
  use nlevel_mod
  use sh95_mod
  use hydrogenic_recombination_mod, only : lyman_line_case
  use utility, only : get_base_name, abort_run
  use iofile_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only : iostat_end
  implicit none
  private
  public :: lines_write, line_selection_setup

  type line_request
     character(len=8) :: ion
     integer :: lower, upper
     real(kind=wp) :: wavelength
     integer :: index = 0
     integer :: row = 0
  end type line_request

  real(kind=wp), parameter :: WAVELENGTH_TOL = 0.0051_wp  ! vacuum A, enough for two decimals
  character(len=4), parameter :: RECOMB_ION(3) = ['h_1 ', 'he_2', 'he_1']
  character(len=4), parameter :: RECOMB_BLOCK(3) = ['h   ', 'heii', 'hei ']
  real(kind=wp), parameter :: BALMER_WL(3) = [6564.60_wp, 4862.65_wp, 4341.66_wp]
  character(len=8), parameter :: BALMER_LABEL(3) = ['Halpha_c', 'Hbeta_c ', 'Hgamma_c']
  type(line_request), allocatable :: requested(:)
  logical :: all_lines, file_lines

contains

  subroutine line_selection_setup()
    type(nlevel_atom_type) :: atom
    character(len=8) :: ion
    character(len=16) :: absent
    real(kind=wp), allocatable :: wavelengths(:)
    real(kind=wp) :: energy
    integer, allocatable :: lower(:), upper(:)
    integer :: ie, stage, k, series
    logical :: ok

    all_lines = par%line_selection == 'all'
    file_lines = index(par%line_selection, 'file:') == 1
    if (all_lines) return

    if (file_lines) then
       call line_list_read(trim(par%line_selection(6:)))
       if (any(requested%ion == 'h_1' .or. requested%ion == 'he_1' .or. requested%ion == 'he_2')) then
          call sh95_setup()
          do series = 1, size(RECOMB_ION)
             ion = RECOMB_ION(series)
             if (.not. any(requested%ion == ion)) cycle
             allocate(lower(sh95_nlines(series)), upper(sh95_nlines(series)), wavelengths(sh95_nlines(series)))
             do k = 1, size(lower)
                call sh95_levels(k, series, lower(k), upper(k))
                wavelengths(k) = sh95_wl(k, series)
             end do
             call resolve_transitions(ion, lower, upper, wavelengths)
             deallocate(lower, upper, wavelengths)
          end do
       end if
       if (par%h_coll_effects .and. any(requested%ion == 'hc_1')) then
          call nlevel_load(trim(par%atomic_dir)//'/nlevel_h_1.txt', atom, ok)
          if (ok) call resolve_transitions('hc_1', [2, 2, 2], [3, 4, 5], BALMER_WL)
       end if
    else
       call diagnostic_transitions()
       absent = ''
       if (par%abund_He <= 0.0_wp) absent = ' he'
       do k = 1, size(requested)
          ion = requested(k)%ion
          if (any([(elem_name(ie) == ion(:index(ion, '_')-1), ie = 1, n_elements)])) cycle
          if (index(absent, ' '//ion(:index(ion, '_')-1)) == 0) &
             absent = trim(absent)//' '//ion(:index(ion, '_')-1)
          requested(k)%index = -1
       end do
       if (mpar%p_rank == 0 .and. len_trim(absent) > 0) &
          write(*,'(a)') ' LINE: diagnostic elements absent from model:'//trim(absent)
    end if

    do ie = 1, n_elements
       do stage = 1, elem_nstage(ie)
          write(ion,'(a,"_",i0)') trim(elem_name(ie)), stage
          if (.not. any(requested%ion == ion)) cycle
          call nlevel_load(metal_excitation_data_file('nlevel', trim(elem_name(ie)), stage), atom, ok)
          if (.not. ok) cycle
          allocate(lower(atom%nrad), upper(atom%nrad), wavelengths(atom%nrad))
          do k = 1, atom%nrad
             call nlevel_transition(atom, k, lower(k), upper(k), energy)
             wavelengths(k) = nlevel_line_ident(atom, k)
          end do
          call resolve_transitions(ion, lower, upper, wavelengths)
          deallocate(lower, upper, wavelengths)
       end do
    end do

    do k = 1, size(requested)
       if (requested(k)%index == 0) call line_request_error(requested(k), 'ion or transition unavailable')
    end do
  end subroutine line_selection_setup

  subroutine diagnostic_transitions()
    ! These are vacuum wavelengths of the same transitions in the two models.
    ! The O II upper-level identifiers are reversed in the MOCASSIN atom.
    requested = [line_request('n_2', 2, 4, 6549.86_wp), line_request('n_2', 3, 4, 6585.27_wp), &
                 line_request('o_1', 1, 4, 6302.05_wp), line_request('o_2', 1, 3, 3729.84_wp), &
                 line_request('o_2', 1, 2, 3727.09_wp), line_request('o_3', 2, 4, 4960.29_wp), &
                 line_request('o_3', 3, 4, 5008.24_wp), line_request('s_2', 1, 2, 6732.67_wp), &
                 line_request('s_2', 1, 3, 6718.29_wp), line_request('s_3', 2, 4, 9070.05_wp), &
                 line_request('s_3', 3, 4, 9532.25_wp)]
    if (par%metal_excitation_model == 'mocassin_v2') then
       requested(3)%wavelength = 6302.03_wp
       requested(4)%upper = 2
       requested(4)%wavelength = 3729.88_wp
       requested(5)%upper = 3
       requested(8:11)%wavelength = [6732.69_wp, 6718.31_wp, 9072.11_wp, 9535.16_wp]
    end if
  end subroutine diagnostic_transitions

  subroutine line_list_read(path)
    character(len=*), intent(in) :: path
    character(len=1024) :: line
    character(len=64) :: ion, wavelength, extra
    character(len=32) :: row_text
    integer :: unit, ios, rows, row, n, comment

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) call abort_run('par%line_selection: cannot open '//path)
    rows = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios == iostat_end) exit
       if (ios /= 0) call abort_run('par%line_selection: cannot read '//path)
       rows = rows + 1
    end do
    rewind(unit)
    allocate(requested(rows))
    n = 0
    do row = 1, rows
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) call abort_run('par%line_selection: cannot reread '//path)
       comment = index(line, '#')
       if (comment > 0) line(comment:) = ''
       if (verify(line, ' '//achar(9)) == 0) cycle
       n = n + 1
       requested(n)%lower = 0
       requested(n)%upper = 0
       wavelength = ''
       ion = ''
       read(line,*,iostat=ios) ion, requested(n)%lower, requested(n)%upper, wavelength, extra
       if (ios /= iostat_end .or. len_trim(line) == len(line) .or. &
           scan(line, ',/*!''"') /= 0 .or. len_trim(ion) > len(requested(n)%ion) .or. &
           len_trim(ion) == 0 .or. requested(n)%lower < 1 .or. requested(n)%upper < 1 .or. &
           len_trim(wavelength) == 0 .or. verify(trim(wavelength), '0123456789.+-eEdD') /= 0) then
          write(row_text,'(i0)') row
          call abort_run('par%line_selection: expected ion lower upper vacuum_A in '//path//':'//trim(row_text))
       end if
       requested(n)%ion = ion
       requested(n)%row = row
       read(wavelength,*,iostat=ios) requested(n)%wavelength
       if (ios /= 0) call line_request_error(requested(n), 'invalid reference wavelength')
       if (requested(n)%wavelength <= 0.0_wp) call line_request_error(requested(n), 'nonpositive reference wavelength')
    end do
    close(unit)
    if (n == 0) call abort_run('par%line_selection: no transitions in '//path)
    requested = requested(:n)
  end subroutine line_list_read

  subroutine resolve_transitions(ion, lower, upper, wavelength)
    character(len=*), intent(in) :: ion
    integer, intent(in) :: lower(:), upper(:)
    real(kind=wp), intent(in) :: wavelength(:)
    integer :: i, k, match
    logical :: taken(size(lower))
    taken = .false.
    do i = 1, size(requested)
       if (requested(i)%ion /= ion) cycle
       match = 0
       do k = 1, size(lower)
          if (requested(i)%lower /= lower(k) .or. requested(i)%upper /= upper(k)) cycle
          if (abs(requested(i)%wavelength - wavelength(k)) > WAVELENGTH_TOL) cycle
          if (match /= 0) call line_request_error(requested(i), 'ambiguous transition')
          match = k
       end do
       if (match == 0) call line_request_error(requested(i), 'no matching levels and vacuum wavelength')
       if (taken(match)) call line_request_error(requested(i), 'duplicate transition request')
       requested(i)%index = match
       taken(match) = .true.
    end do
  end subroutine resolve_transitions

  subroutine line_request_error(line, reason)
    type(line_request), intent(in) :: line
    character(len=*), intent(in) :: reason
    character(len=PATH_LEN+32) :: location
    character(len=100) :: transition
    location = par%line_selection
    if (line%row > 0) write(location,'(a,":",i0)') trim(par%line_selection), line%row
    write(transition,'(a,2(1x,i0),1x,g0,1x,a)') trim(line%ion), line%lower, line%upper, &
         line%wavelength, 'vacuum A'
    call abort_run('par%line_selection: '//trim(location)//': '//trim(transition)//': '//reason)
  end subroutine line_request_error

  function selected_transitions(ion, n, hydrogen_helium) result(indices)
    character(len=*), intent(in) :: ion
    integer, intent(in) :: n
    logical, intent(in) :: hydrogen_helium
    integer, allocatable :: indices(:)
    integer :: j
    logical :: chosen(n)
    chosen = all_lines .or. (hydrogen_helium .and. .not. file_lines)
    if (.not. all_lines .and. (.not. hydrogen_helium .or. file_lines)) then
       do j = 1, size(requested)
          if (requested(j)%ion == ion .and. requested(j)%index > 0) chosen(requested(j)%index) = .true.
       end do
    end if
    indices = pack([(j, j = 1, n)], chosen)
  end function selected_transitions

  subroutine lines_write()
    use gas_rates_mod, only : run_stop_reason
    type(io_file_type) :: efile
    character(len=PATH_LEN) :: outname, emisname
    real(kind=wp) :: LHb
    integer :: unit, status, ierr, ie, ist

    status = 0
    LHb = 0.0_wp
    unit = 0
    if (mpar%p_rank == 0) then
       call sh95_setup()
       outname = trim(get_base_name(par%out_file))//'_lines.txt'
       open(newunit=unit, file=trim(outname), status='replace', iostat=status)
       if (status /= 0) call abort_run('cannot open '//trim(outname))
       if (par%emis_output) then
          emisname = trim(get_base_name(par%out_file))//'_emis'//trim(io_file_extension(par%file_format))
          call io_open_new(efile, emisname, status)
          call emis_write_state(efile, status)
       end if
       write(unit,'(a)') '# MoCHII line luminosities ('//trim(run_stop_reason)//')'
       write(unit,'(a)') '# H I recombination lines: Storey & Hummer (1995) '// &
            'case '//trim(lyman_line_case())//', bilinear (log T, log ne) interpolation'
       write(unit,'(a)') '# (SH95 H I/He II lines use case '// &
            trim(lyman_line_case())//'; He I Porter lines are case B)'
       write(unit,'(a)') '# metals: n-level solve leaf by leaf'
       call recombination_lines_write(efile, status, unit, LHb)
       if (par%h_coll_effects) call collisional_balmer_write(efile, status, unit, LHb)
       if (par%emis_output .and. all_lines .and. LHb <= 0.0_wp) &
          write(*,'(a)') ' LINE: zero Hbeta luminosity; no metal maps for line_selection=all.'
    end if
    call MPI_BCAST(LHb, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    do ie = 1, n_elements
       do ist = 1, elem_nstage(ie)
          call metal_lines_write(ie, ist, efile, status, unit, LHb)
       end do
    end do

    if (mpar%p_rank /= 0) return
    close(unit)
    write(*,'(2a)') ' LINE: line luminosities written to: ', trim(outname)
    if (par%emis_output) then
       call io_close(efile, status)
       if (status /= 0) call abort_run('failed to write emissivity file: '//trim(emisname))
       write(*,'(2a)') ' LINE: leaf emissivities + state written to: ', trim(emisname)
    end if
  end subroutine lines_write

  real(kind=wp) function hbeta_ratio(luminosity, LHb) result(ratio)
    real(kind=wp), intent(in) :: luminosity, LHb
    if (LHb > 0.0_wp) then
       ratio = luminosity / LHb
    else
       ratio = ieee_value(0.0_wp, ieee_quiet_nan)
    end if
  end function hbeta_ratio

  subroutine recombination_lines_write(efile, status, unit, LHb)
    type(io_file_type), intent(inout) :: efile
    integer, intent(inout) :: status
    integer, intent(in) :: unit
    real(kind=wp), intent(inout) :: LHb
    integer :: series, n, il, j, k, hb, stage
    integer, allocatable :: chosen(:), evaluate(:)
    real(kind=wp), allocatable :: maps(:,:), luminosity(:), wavelengths(:)
    real(kind=wp) :: ne, nH, T, density, volume, em
    character(len=2) :: el

    do series = 1, size(RECOMB_ION)
       if (series > 1 .and. .not. all_lines .and. .not. file_lines .and. par%abund_He <= 0.0_wp) cycle
       n = sh95_nlines(series)
       chosen = selected_transitions(RECOMB_ION(series), n, .true.)
       evaluate = chosen
       hb = 0
       if (series == 1) then
          do k = 1, n
             if (sh95_label(k) == 'Hbeta') hb = k
          end do
          if (hb == 0) call abort_run('Hbeta is missing from the H I recombination table')
          if (.not. any(chosen == hb)) evaluate = [chosen, hb]
       end if
       if (size(evaluate) == 0) cycle
       allocate(luminosity(n), source=0.0_wp)
       if (par%emis_output .and. size(chosen) > 0) allocate(maps(size(chosen),gas_nleaf), source=0.0_wp)
       do il = 1, gas_nleaf
          nH = gas_nH(il)
          if (nH <= 0.0_wp) cycle
          ne = gas_ne(il)
          T = gas_Te(il)
          volume = (2.0_wp*leaf_half(il)*par%distance2cm)**3
          select case (series)
          case (1)
             density = ne*nH*gas_xHII(il)
          case (2)
             if (gas_xHeIII(il) <= 1.0e-30_wp) cycle
             density = ne*nH*par%abund_He*gas_xHeIII(il)
          case (3)
             if (gas_xHeII(il) <= 1.0e-30_wp) cycle
             density = ne*nH*par%abund_He*gas_xHeII(il)
          end select
          do j = 1, size(evaluate)
             k = evaluate(j)
             em = density*sh95_emis(k, T, ne, series)
             luminosity(k) = luminosity(k) + em*volume
             if (allocated(maps) .and. j <= size(chosen)) maps(j, il) = em
          end do
       end do
       if (series == 1) then
          LHb = luminosity(hb)
          write(unit,'(a,es14.6,a)') '# L(Hbeta) = ', LHb, ' erg/s'
          write(unit,'(a)') '# elem stage  lambda[A]      L[erg/s]     L/L(Hbeta)'
       end if
       el = 'he'
       if (series == 1) el = 'h'
       stage = 1
       if (series == 2) stage = 2
       wavelengths = [(sh95_wl(chosen(j), series), j = 1, size(chosen))]
       do j = 1, size(chosen)
          k = chosen(j)
          write(unit,'(a4,i4,f12.2,2es15.5e3,2x,a)') trim(el), stage, wavelengths(j), &
               luminosity(k), hbeta_ratio(luminosity(k), LHb), trim(sh95_label(k, series))
       end do
       if (allocated(maps)) then
          call emissivity_block_write(efile, status, trim(RECOMB_BLOCK(series)), maps, wavelengths)
          deallocate(maps)
       end if
       deallocate(luminosity)
    end do
  end subroutine recombination_lines_write

  subroutine collisional_balmer_write(efile, status, unit, LHb)
    type(io_file_type), intent(inout) :: efile
    integer, intent(inout) :: status
    integer, intent(in) :: unit
    real(kind=wp), intent(in) :: LHb
    type(nlevel_atom_type) :: atom
    real(kind=wp), allocatable :: emc(:,:), emis(:)
    real(kind=wp) :: Lc(3), gemis(3), wlk, nHI, volume
    integer, allocatable :: chosen(:), blend(:)
    integer :: il, k, j
    logical :: ok

    chosen = selected_transitions('hc_1', size(BALMER_WL), .true.)
    if (size(chosen) == 0) return
    call nlevel_load(trim(par%atomic_dir)//'/nlevel_h_1.txt', atom, ok)
    if (.not. ok) then
       write(*,'(a)') ' LINE: nlevel_h_1.txt unavailable; skipping collisional Balmer lines.'
       return
    end if
    allocate(emis(atom%nrad), blend(atom%nrad))
    blend = 0
    do k = 1, atom%nrad
       wlk = nlevel_line_ident(atom, k)
       if (wlk >= 6520.0_wp .and. wlk <= 6600.0_wp) blend(k) = 1
       if (wlk >= 4840.0_wp .and. wlk <= 4880.0_wp) blend(k) = 2
       if (wlk >= 4320.0_wp .and. wlk <= 4360.0_wp) blend(k) = 3
    end do
    Lc = 0.0_wp
    if (par%emis_output) allocate(emc(size(chosen),gas_nleaf), source=0.0_wp)
    do il = 1, gas_nleaf
       nHI = gas_nH(il)*gas_xHI(il)
       if (nHI <= 1.0e-30_wp) cycle
       volume = (2.0_wp*leaf_half(il)*par%distance2cm)**3
       ! The optically thin atomic cascade supplies the existing blend groups.
       call nlevel_emissivities(atom, gas_Te(il), gas_ne(il), emis)
       gemis = 0.0_wp
       do k = 1, atom%nrad
          j = blend(k)
          if (j > 0) gemis(j) = gemis(j) + emis(k)
       end do
       Lc = Lc + nHI*gemis*volume
       if (allocated(emc)) emc(:, il) = nHI*gemis(chosen)
    end do
    write(unit,'(a)') '# hc: collisional H I Balmer component'// &
         ' (excitation of neutral H from the ground level; separate from H I recombination)'
    do j = 1, size(chosen)
       k = chosen(j)
       write(unit,'(a4,i4,f12.2,2es15.5e3,2x,a)') 'hc', 1, BALMER_WL(k), Lc(k), &
            hbeta_ratio(Lc(k), LHb), trim(BALMER_LABEL(k))
    end do
    if (allocated(emc)) call emissivity_block_write(efile, status, 'hc', emc, BALMER_WL(chosen))
  end subroutine collisional_balmer_write

  subroutine metal_lines_write(ie, ist, efile, status, unit, LHb)
    integer, intent(in) :: ie, ist, unit
    type(io_file_type), intent(inout) :: efile
    integer, intent(inout) :: status
    real(kind=wp), intent(in) :: LHb
    type(nlevel_atom_type) :: atom
    character(len=8) :: ion
    real(kind=wp), allocatable :: luminosity(:), wavelength(:), einstein(:), energy(:), pop(:), emis(:)
    real(kind=wp), allocatable :: upper_pop(:,:), ion_density(:), maps(:,:), packed(:,:), received(:,:)
    real(kind=wp), allocatable :: electron_rates(:,:), radiative_rates(:,:)
    real(kind=wp) :: nH, nHI, nion, T, ne, volume, ignored(1), rates_T
    integer, allocatable :: chosen(:), mapped(:), upper(:), slot(:), used_levels(:), level_slot(:)
    integer, allocatable :: counts(:), offsets(:)
    integer :: n, i, k, lower, il, rank_first, stride, nowned, owned, start, count_leaf, chunk_size, ierr
    integer :: rank, first_leaf, last_leaf, leaf_count
    integer, parameter :: EMIS_CHUNK_VALUES = 2**24
    logical :: ok, save_maps

    write(ion,'(a,"_",i0)') trim(elem_name(ie)), ist
    if (.not. all_lines) then
       if (.not. any(requested%ion == ion)) return
    end if
    call nlevel_load(metal_excitation_data_file('nlevel', trim(elem_name(ie)), ist), atom, ok)
    if (.not. ok) return
    chosen = selected_transitions(ion, atom%nrad, .false.)
    n = size(chosen)
    if (n == 0) return

    allocate(luminosity(n), wavelength(n), einstein(n), energy(n), upper(n), pop(atom%nlev), emis(n))
    luminosity = 0.0_wp
    einstein = atom%A(chosen)
    do i = 1, n
       k = chosen(i)
       call nlevel_transition(atom, k, lower, upper(i), energy(i))
       wavelength(i) = nlevel_line_ident(atom, k)
    end do

    rank_first = gas_leaf_owner_first()
    stride = gas_leaf_owner_stride()
    nowned = (gas_nleaf + stride - rank_first)/stride
    save_maps = par%emis_output .and. (.not. all_lines .or. LHb > 0.0_wp)
    if (save_maps) then
       ! Save the emitting-level populations and the ion density of each
       ! owned leaf for this ion only. In 'all' the integrated flux determines
       ! which branches to map, after the population solve has finished.
       allocate(level_slot(atom%nlev), source=0)
       do i = 1, n
          level_slot(upper(i)) = 1
       end do
       used_levels = pack([(i, i = 1, atom%nlev)], level_slot /= 0)
       do i = 1, size(used_levels)
          level_slot(used_levels(i)) = i
       end do
       slot = level_slot(upper)
       allocate(upper_pop(size(used_levels), nowned))
       allocate(ion_density(nowned), source=0.0_wp)
    end if
    allocate(electron_rates(atom%nlev,atom%nlev), radiative_rates(atom%nlev,atom%nlev))
    rates_T = -1.0_wp
    do owned = 1, nowned
       il = rank_first + (owned-1)*stride
       nH = gas_nH(il)
       if (nH <= 0.0_wp) cycle
       T = gas_Te(il)
       ne = gas_ne(il)
       nHI = 0.0_wp
       if (par%neutral_hydrogen_excitation) nHI = nH*gas_xHI(il)
       nion = elem_abund(ie)*nH*species_stage_fraction(ie, ist, il)
       if (nion <= 1.0e-30_wp) cycle
       volume = (2.0_wp*leaf_half(il)*par%distance2cm)**3
       ! The electron-collision and radiative rate matrices depend on T
       ! alone: rebuild them only when T differs from the last leaf's.  Runs
       ! of equal T are common (most neutral leaves sit exactly at te_min; a
       ! fixed-temperature model has one T), and the 196-level O II matrix is
       ! most of the cost of its solve.
       if (T /= rates_T) then
          call nlevel_cooling_prepT(atom, T, electron_rates, radiative_rates)
          rates_T = T
       end if
       call nlevel_line_populations(atom, T, ne, pop, nHI, electron_rates, radiative_rates)
       emis = pop(upper)*einstein*energy*nion
       luminosity = luminosity + emis*volume
       if (save_maps) then
          upper_pop(:,owned) = pop(used_levels)
          ion_density(owned) = nion
       end if
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE, luminosity, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (mpar%p_rank == 0) then
       do i = 1, n
          write(unit,'(a4,i4,f12.2,2es15.5e3)') trim(elem_name(ie)), ist, wavelength(i), &
               luminosity(i), hbeta_ratio(luminosity(i), LHb)
       end do
    end if

    if (.not. save_maps) return
    mapped = [(i, i = 1, n)]
    if (all_lines) mapped = pack(mapped, luminosity >= par%emis_line_floor*LHb .and. luminosity > 0.0_wp)
    n = size(mapped)
    if (n == 0) return
    chunk_size = max(1, min(gas_nleaf, EMIS_CHUNK_VALUES/n))
    allocate(packed(n, 1+(chunk_size-1)/stride), counts(stride), offsets(stride))
    if (mpar%p_rank == 0) allocate(maps(n,gas_nleaf), received(n,chunk_size))
    do start = 1, gas_nleaf, chunk_size
       count_leaf = min(chunk_size, gas_nleaf - start + 1)
       last_leaf = start + count_leaf - 1
       first_leaf = start + modulo(rank_first - start, stride)
       leaf_count = (last_leaf - first_leaf + stride)/stride
       packed(:,1:leaf_count) = 0.0_wp
       do il = first_leaf, last_leaf, stride
          owned = (il-rank_first)/stride + 1
          if (ion_density(owned) <= 0.0_wp) cycle
          i = (il-first_leaf)/stride + 1
          packed(:,i) = upper_pop(slot(mapped),owned)*einstein(mapped)*energy(mapped)*ion_density(owned)
       end do
       ! Pack just the leaves owned by each world rank. Restore the original
       ! leaf order at the writer, without sending zero copies from other ranks.
       if (mpar%p_rank == 0) then
          offsets(1) = 0
          do rank = 1, stride
             first_leaf = start + modulo(rank - start, stride)
             counts(rank) = n*((last_leaf - first_leaf + stride)/stride)
             if (rank < stride) offsets(rank+1) = offsets(rank) + counts(rank)
          end do
          call MPI_GATHERV(packed, n*leaf_count, MPI_DOUBLE_PRECISION, received, counts, offsets, &
                           MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
          do rank = 1, stride
             first_leaf = start + modulo(rank - start, stride)
             i = offsets(rank)/n
             do il = first_leaf, last_leaf, stride
                i = i + 1
                maps(:,il) = received(:,i)
             end do
          end do
       else
          call MPI_GATHERV(packed, n*leaf_count, MPI_DOUBLE_PRECISION, ignored, counts, offsets, &
                           MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       end if
    end do
    if (mpar%p_rank == 0) call emissivity_block_write(efile, status, trim(ion), maps, wavelength(mapped))
  end subroutine metal_lines_write

  subroutine emissivity_block_write(efile, status, ion, maps, wavelengths)
    type(io_file_type), intent(inout) :: efile
    integer, intent(inout) :: status
    character(len=*), intent(in) :: ion
    real(kind=wp), intent(in) :: maps(:,:), wavelengths(:)
    call io_append_image(efile, maps, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','emis_'//ion, 'emissivities [erg/s/cm^3]; sum(emis*V)=L',status)
    call io_append_image(efile, wavelengths, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','wl_'//ion, 'vacuum wavelengths [A]',status)
  end subroutine emissivity_block_write

  subroutine emis_write_state(efile, status)
    use iofile_mod
    implicit none
    type(io_file_type), intent(inout) :: efile
    integer, intent(inout) :: status
    real(kind=wp), allocatable :: tmp(:), lxyz(:,:)
    integer :: il

    if (status /= 0) return
    allocate(lxyz(gas_nleaf,3))
    do il = 1, gas_nleaf
       lxyz(il,1) = leaf_cx(il)
       lxyz(il,2) = leaf_cy(il)
       lxyz(il,3) = leaf_cz(il)
    end do
    call io_append_image(efile, lxyz, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','LeafXYZ', &
         'leaf center x,y,z (code units)',status)
    call io_put_keyword(efile,'DIST_CM', par%distance2cm, &
         'distance unit (cm)', status)
    deallocate(lxyz)

    allocate(tmp(gas_nleaf))
    do il = 1, gas_nleaf
       tmp(il) = 2.0_wp*leaf_half(il)        ! leaf edge length, code units
    end do
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','LeafSize', &
         'leaf edge length (code units); V = (size*DIST_CM)^3',status)
    tmp = gas_nH(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','n_H','H density [cm^-3]',status)
    tmp = gas_ne(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','n_e','electron density [cm^-3]',status)
    tmp = gas_Te(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','T_e','electron temperature [K]',status)
    tmp = gas_xHI(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HI','n_HI/n_H',status)
    tmp = gas_xHII(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HII','n_HII/n_H',status)
    tmp = gas_xHeI(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HeI','n_HeI/n_He',status)
    tmp = gas_xHeII(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HeII','n_HeII/n_He',status)
    tmp = gas_xHeIII(1:gas_nleaf)
    call io_append_image(efile, tmp, status, bitpix=-64)
    call io_put_keyword(efile,'EXTNAME','x_HeIII','n_HeIII/n_He',status)
    deallocate(tmp)

    !--- metal stage fractions (same blocks as the rates file).
    if (par%use_metals) call species_write(efile, status, include_rates=.false.)
  end subroutine emis_write_state

end module lines_mod
