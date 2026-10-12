program check_small_excitation_cooling
! Exercise absolute registered cooling, including positive subnormal values.
  use mpi
  use define, only: wp, par, mpar
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_cooling_at
  use nlevel_cooling_mod, only: nlevel_cooling_add, nlevel_cooling_lam_se, &
       nlevel_cooling_apply
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(nlevel_atom_type) :: atom
  character(len=1024) :: directory, name, stage_arg, output, atom_file
  real(kind=wp), parameter :: temperatures(4) = [100.0_wp, 10.0_wp**2.0005_wp, &
       10.0_wp**2.005_wp, 1000.0_wp]
  real(kind=wp), parameter :: densities(3) = [1.0e-6_wp, 10.0_wp**(-5.95_wp), 1.0e2_wp]
  real(kind=wp) :: t, ne, direct, table, reference, fit, applied, expected
  integer :: ierr, stage, ios, idx, i, j, k, unit
  logical :: ok

  call MPI_INIT(ierr)
  if (command_argument_count() /= 4) error stop &
       'usage: check_small_excitation_cooling.x ATOMIC_DIR ELEMENT STAGE OUTPUT'
  call get_command_argument(1, directory)
  call get_command_argument(2, name)
  call get_command_argument(3, stage_arg)
  call get_command_argument(4, output)
  read(stage_arg, *, iostat=ios) stage
  if (ios /= 0 .or. stage < 1) error stop 'Invalid ion stage'
  par%atomic_dir = trim(directory)
  par%fe_levels_full = .false.
  mpar%p_rank = 0
  write(atom_file, '(a,a,a,i0,a)') trim(directory)//'/nlevel_', &
       trim(name), '_', stage, '.txt'
  call nlevel_load(trim(atom_file), atom, ok)
  if (.not. ok) error stop 'Atom did not load'
  idx = nlevel_cooling_add(trim(name), stage)
  if (idx <= 0) error stop 'Cooling table was not registered'
  open(newunit=unit, file=trim(output), status='new', action='write', iostat=ios)
  if (ios /= 0) error stop 'Cannot create output'
  write(unit, '(a,l1)') '# density_ratio_cooling ', atom%density_ratio_cooling
  write(unit, '(a)') '# T ne fit_multiplier direct table reference fit applied'
  do i = 1, size(temperatures)
     t = temperatures(i)
     reference = nlevel_cooling_lam_se(idx, t, densities(1))
     do j = 1, size(densities)
        ne = densities(j)
        direct = nlevel_cooling_at(atom, t, ne)
        table = nlevel_cooling_lam_se(idx, t, ne)
        do k = 1, 4
           select case (k)
           case (1)
              fit = reference
           case (2)
              fit = 2.0_wp*reference
           case (3)
              fit = 0.5_wp*reference
           case (4)
              fit = 1.0e-25_wp
           end select
           applied = nlevel_cooling_apply(idx, fit, t, ne)
           if (.not. all(ieee_is_finite([direct, table, reference, fit, applied]))) &
                error stop 'Nonfinite cooling value'
           if (min(direct, table, reference, fit, applied) < 0.0_wp) &
                error stop 'Negative cooling value'
           if (direct > 0.0_wp) then
              if (table <= 0.0_wp) error stop 'Positive cooling was lost'
              if (abs(table/direct-1.0_wp) > 5.0e-3_wp) &
                   error stop 'Registered cooling differs from direct cooling'
           else if (table /= 0.0_wp) then
              error stop 'Exact zero cooling was not preserved'
           end if
           if (atom%density_ratio_cooling) then
              if (reference > 0.0_wp) then
                 expected = fit*(table/reference)
              else
                 expected = table
              end if
           else if (reference <= fit) then
              expected = table + (fit-reference)
           else
              expected = (fit/reference)*table
           end if
           if (expected > 0.0_wp) then
              if (abs(applied/expected-1.0_wp) > 1.0e-11_wp) &
                   error stop 'Cooling combination differs from expected value'
           else if (applied /= 0.0_wp) then
              error stop 'Exact zero combination was not preserved'
           end if
           ! Independent identities: a fit equal to the source reference
           ! recovers the absolute registered table at every density.
           if (k == 1) then
              if (table > 0.0_wp) then
                 if (abs(applied/table-1.0_wp) > 1.0e-11_wp) &
                      error stop 'Reference fit did not recover the absolute table'
              else if (applied /= 0.0_wp) then
                 error stop 'Reference fit did not preserve the zero table'
              end if
           end if
           ! With a positive reference, the reference-density table is the
           ! fit itself. For a true zero reference, the complete atom follows
           ! the zero source and the partial atom retains the fit residual.
           if (j == 1) then
              if (reference > 0.0_wp .or. .not. atom%density_ratio_cooling) then
                 if (fit > 0.0_wp) then
                    if (abs(applied/fit-1.0_wp) > 1.0e-11_wp) &
                         error stop 'Low-density fit was not recovered'
                 else if (applied /= 0.0_wp) then
                    error stop 'Low-density zero fit was not recovered'
                 end if
              else if (applied /= table) then
                 error stop 'Complete atom did not preserve true zero source'
              end if
           end if
           write(unit, '(2(i0,1x),*(es25.16e3,1x))') i, k, t, ne, direct, table, &
                reference, fit, applied
        end do
     end do
  end do
  close(unit)
  call MPI_FINALIZE(ierr)
end program check_small_excitation_cooling
