program check_registry
  use mpi
  use define, only: wp, par, mpar
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_cooling_at
  use nlevel_cooling_mod, only: nlevel_cooling_add, nlevel_cooling_lam_se
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
       ieee_quiet_nan, ieee_positive_inf
  implicit none
  type(nlevel_atom_type) :: atom
  character(len=1024) :: directory, name, stage_text, mode, atom_file
  real(kind=wp) :: temperatures(6), densities(4), direct, table, t, ne
  integer :: ierr, stage, ios, idx, i, j, ntemp
  logical :: ok

  call MPI_INIT(ierr)
  if (command_argument_count() /= 4) error stop 'usage: DRIVER DIR NAME STAGE MODE'
  call get_command_argument(1, directory)
  call get_command_argument(2, name)
  call get_command_argument(3, stage_text)
  call get_command_argument(4, mode)
  read(stage_text, *, iostat=ios) stage
  if (ios /= 0) error stop 'invalid stage'
  par%atomic_dir = trim(directory)
  par%fe_levels_full = .false.
  mpar%p_rank = 0
  write(atom_file,'(a,a,a,i0,a)') trim(directory)//'/nlevel_', &
       trim(name), '_', stage, '.txt'
  call nlevel_load(trim(atom_file), atom, ok)
  if (.not. ok) error stop 'atom did not load'
  idx = nlevel_cooling_add(trim(name), stage)
  if (idx <= 0) error stop 'atom did not register'
  if (trim(mode) == 'reject') then
     table = nlevel_cooling_lam_se(idx, 158124.0_wp, 100.0_wp)
     print *, 'unexpected value', table
     error stop 'unsupported temperature was accepted'
  end if
  if (trim(mode) == 'reject_near') then
     table = nlevel_cooling_lam_se(idx, 999901.0_wp, 100.0_wp)
     print *, 'unexpected value', table
     error stop 'unsupported near-ceiling temperature was accepted'
  end if
  if (trim(mode) == 'nanT') then
     table = nlevel_cooling_lam_se(idx, ieee_value(1.0_wp, ieee_quiet_nan), 100.0_wp)
     error stop 'NaN temperature was accepted'
  end if
  if (trim(mode) == 'infT') then
     table = nlevel_cooling_lam_se(idx, ieee_value(1.0_wp, ieee_positive_inf), 100.0_wp)
     error stop 'infinite temperature was accepted'
  end if
  if (trim(mode) == 'zeroT') then
     table = nlevel_cooling_lam_se(idx, 0.0_wp, 100.0_wp)
     error stop 'zero temperature was accepted'
  end if
  if (trim(mode) == 'nanNe') then
     table = nlevel_cooling_lam_se(idx, 1000.0_wp, ieee_value(1.0_wp, ieee_quiet_nan))
     error stop 'NaN electron density was accepted'
  end if
  if (trim(mode) == 'infNe') then
     table = nlevel_cooling_lam_se(idx, 1000.0_wp, ieee_value(1.0_wp, ieee_positive_inf))
     error stop 'infinite electron density was accepted'
  end if
  if (trim(mode) == 'zeroNe') then
     table = nlevel_cooling_lam_se(idx, 1000.0_wp, 0.0_wp)
     error stop 'zero electron density was accepted'
  end if
  if (trim(mode) == 'argon') then
     temperatures = [100.0_wp, 1000.0_wp, 10000.0_wp, &
          154881.661891_wp, 156000.0_wp, 158000.0_wp]
     ntemp = 6
  else if (trim(mode) == 'toy') then
     temperatures = [100.0_wp, 1000.0_wp, 10000.0_wp, &
          154881.661891_wp, 156000.0_wp, 158123.0_wp]
     ntemp = 6
  else if (trim(mode) == 'unlimited') then
     temperatures = [100.0_wp, 1000.0_wp, 10000.0_wp, &
          100000.0_wp, 1000000.0_wp, 10000000.0_wp]
     ntemp = 6
  else if (trim(mode) == 'near') then
     temperatures = [100.0_wp, 1000.0_wp, 10000.0_wp, &
          977237.220956_wp, 990000.0_wp, 999900.0_wp]
     ntemp = 6
  else
     error stop 'invalid mode'
  end if
  densities = [1.0e-6_wp, 1.0e2_wp, 1.0e6_wp, 1.0e10_wp]
  do i = 1, ntemp
     t = temperatures(i)
     do j = 1, size(densities)
        ne = densities(j)
        table = nlevel_cooling_lam_se(idx, t, ne)
        if (trim(mode) == 'unlimited' .and. t > 1.0e6_wp) then
           direct = nlevel_cooling_at(atom, 1.0e6_wp, ne)
        else
           direct = nlevel_cooling_at(atom, t, ne)
        end if
        if (.not. ieee_is_finite(table) .or. .not. ieee_is_finite(direct)) &
             error stop 'nonfinite cooling'
        write(*,'(a,1x,a,1x,i0,4(1x,es25.16))') 'ROW', trim(mode), &
             i, t, ne, direct, table
     end do
  end do
  call MPI_FINALIZE(ierr)
end program check_registry
