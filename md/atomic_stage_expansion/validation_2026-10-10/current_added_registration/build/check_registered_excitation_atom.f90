program check_registered_excitation_atom
! Check the registered 401-by-161 cooling table against its direct atom solve.
  use mpi
  use define, only: wp, par, mpar
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_cooling_at
  use nlevel_cooling_mod, only: nlevel_cooling_add, nlevel_cooling_lam_se, &
       nlevel_cooling_apply, nlevel_cooling_report
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  type(nlevel_atom_type) :: atom
  real(kind=wp), parameter :: grid_t(5) = [100.0_wp, 1000.0_wp, &
       1.0e4_wp, 1.0e5_wp, 1.0e6_wp]
  real(kind=wp), parameter :: grid_ne(4) = [1.0e-6_wp, 1.0e2_wp, &
       1.0e6_wp, 1.0e10_wp]
  real(kind=wp), parameter :: middle_log_t(5) = [2.005_wp, 3.005_wp, &
       4.005_wp, 5.005_wp, 5.995_wp]
  real(kind=wp), parameter :: middle_log_ne(4) = [-5.95_wp, 2.05_wp, &
       6.05_wp, 9.95_wp]
  character(len=1024) :: directory, name, stage_arg, output, atom_file
  character(len=8) :: sample_kind
  real(kind=wp) :: t, ne, direct, table, applied, reference, table_floor
  real(kind=wp) :: expected_apply, table_error, apply_error
  real(kind=wp) :: worst_grid, worst_middle, worst_apply, elapsed
  integer(kind=8) :: tick_start, tick_end, tick_rate
  integer :: ierr, stage, ios, idx, i, j, pass_index, unit, failures, rows
  logical :: ok, normal, density_ratio

  call MPI_INIT(ierr)
  if (command_argument_count() /= 4) error stop &
       'usage: check_registered_excitation_atom.x ATOMIC_DIR ELEMENT STAGE OUTPUT'
  call get_command_argument(1, directory)
  call get_command_argument(2, name)
  call get_command_argument(3, stage_arg)
  call get_command_argument(4, output)
  read(stage_arg, *, iostat=ios) stage
  if (ios /= 0 .or. stage < 1) error stop 'Invalid ion stage'
  if (len_trim(directory) == 0 .or. len_trim(name) == 0 .or. &
      len_trim(output) == 0) error stop 'Empty argument'

  par%atomic_dir = trim(directory)
  par%fe_levels_full = .false.
  mpar%p_rank = 0
  write(atom_file, '(a,a,a,i0,a)') trim(directory)//'/nlevel_', &
       trim(name), '_', stage, '.txt'
  call nlevel_load(trim(atom_file), atom, ok)
  if (.not. ok) error stop 'Atom did not load'

  open(newunit=unit, file=trim(output), status='new', action='write', iostat=ios)
  if (ios /= 0) error stop 'Cannot open result file'
  write(unit, '(a)') '# Registered cooling table and direct n-level solve'
  write(unit, '(a)') '# kind T[K] ne[cm^-3] direct table applied direct_floor '// &
       'table_floor relative_table_error relative_apply_error'
  call system_clock(tick_start, tick_rate)
  idx = nlevel_cooling_add(trim(name), stage)
  call system_clock(tick_end)
  elapsed = real(tick_end-tick_start, wp)/real(tick_rate, wp)
  if (idx <= 0) then
     write(unit, '(a)') '# FAIL: ion was not registered'
     close(unit)
     error stop 'Ion was not registered'
  end if
  write(unit, '(a,i0)') '# nlev ', atom%nlev
  write(unit, '(a,i0)') '# nrad ', atom%nrad
  write(unit, '(a,i0)') '# nups ', atom%nups
  density_ratio = atom%density_ratio_cooling
  write(unit, '(a,l1)') '# density_ratio_cooling ', density_ratio
  if (density_ratio) then
     write(unit, '(a)') '# apply is compared with the direct selected-atom cooling'
  else
     write(unit, '(a)') '# apply is checked against the legacy partial-atom combination'
  end if
  write(unit, '(a,es24.16)') '# registration_wall_seconds ', elapsed

  failures = 0
  rows = 0
  worst_grid = 0.0_wp
  worst_middle = 0.0_wp
  worst_apply = 0.0_wp
  do pass_index = 1, 2
     if (pass_index == 1) then
        sample_kind = 'grid'
     else
        sample_kind = 'middle'
     end if
     do i = 1, size(grid_t)
        if (pass_index == 1) then
           t = grid_t(i)
        else
           t = 10.0_wp**middle_log_t(i)
        end if
        reference = nlevel_cooling_at(atom, t, grid_ne(1))
        table_floor = nlevel_cooling_lam_se(idx, t, grid_ne(1))
        do j = 1, size(grid_ne)
           if (pass_index == 1) then
              ne = grid_ne(j)
           else
              ne = 10.0_wp**middle_log_ne(j)
           end if
           direct = nlevel_cooling_at(atom, t, ne)
           table = nlevel_cooling_lam_se(idx, t, ne)
           applied = nlevel_cooling_apply(idx, reference, t, ne)
           if (density_ratio) then
              if (table_floor > 0.0_wp) then
                 expected_apply = reference*(table/table_floor)
              else
                 expected_apply = table
              end if
           else if (table_floor <= reference) then
              expected_apply = table + (reference-table_floor)
           else
              expected_apply = (reference/table_floor)*table
           end if
           normal = direct > 1.0e-200_wp
           if (normal) then
              table_error = abs(table/direct-1.0_wp)
              apply_error = abs(applied/direct-1.0_wp)
           else
              table_error = abs(table-direct)
              apply_error = abs(applied-direct)
           end if
           if (.not. all(ieee_is_finite([direct, table, applied, reference, &
                table_floor, expected_apply, table_error, apply_error]))) failures = failures + 1
           if (normal) then
              if (pass_index == 1) then
                 worst_grid = max(worst_grid, table_error)
                 if (table_error > 1.0e-10_wp) failures = failures + 1
              else
                 worst_middle = max(worst_middle, table_error)
                 if (table_error > 5.0e-3_wp) failures = failures + 1
              end if
              worst_apply = max(worst_apply, apply_error)
              if (density_ratio) then
                 if (pass_index == 1 .and. apply_error > 1.0e-10_wp) failures = failures + 1
                 if (pass_index == 2 .and. apply_error > 5.0e-3_wp) failures = failures + 1
              end if
           else
              if (table_error > 1.0e-200_wp) failures = failures + 1
              if (density_ratio .and. apply_error > 1.0e-200_wp) failures = failures + 1
           end if
           if (abs(applied-expected_apply) > &
               1.0e-12_wp*max(abs(expected_apply), 1.0e-200_wp)) failures = failures + 1
           write(unit, '(a8,1x,*(es24.16,1x))') sample_kind, t, ne, direct, &
                table, applied, reference, table_floor, table_error, apply_error
           rows = rows + 1
        end do
     end do
  end do
  write(unit, '(a,i0)') '# rows ', rows
  write(unit, '(a,es24.16)') '# worst_grid_relative_error ', worst_grid
  write(unit, '(a,es24.16)') '# worst_middle_relative_error ', worst_middle
  write(unit, '(a,es24.16)') '# worst_apply_relative_error ', worst_apply
  write(unit, '(a,i0)') '# failures ', failures
  if (failures == 0) then
     write(unit, '(a)') '# PASS'
  else
     write(unit, '(a)') '# FAIL'
  end if
  close(unit)
  call nlevel_cooling_report()
  if (failures /= 0 .or. rows /= 40) error stop 'Registered cooling check failed'
  call MPI_FINALIZE(ierr)
end program check_registered_excitation_atom
