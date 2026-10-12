program check_si_extended
  use mpi
  use define, only: wp
  use cooling_mod, only: cooling_fit_type, cooling_load, cooling_eval
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(cooling_fit_type) :: fit
  character(len=1024) :: input, output
  real(kind=wp) :: temperature, cooling, minimum, maximum
  integer :: i, unit, ierr, zero_count
  call MPI_INIT(ierr)
  call get_command_argument(1, input)
  call get_command_argument(2, output)
  call cooling_load(trim(input), fit)
  minimum = huge(1.0_wp)
  maximum = 0.0_wp
  zero_count = 0
  do i = 0, 10000
     temperature = 10.0_wp**(9.0_wp*real(i,wp)/10000.0_wp)
     cooling = cooling_eval(fit, temperature)
     if (.not.ieee_is_finite(cooling) .or. cooling < 0.0_wp) error stop 'Invalid cooling'
     minimum = min(minimum, cooling)
     maximum = max(maximum, cooling)
     if (cooling == 0.0_wp) zero_count = zero_count + 1
  end do
  open(newunit=unit, file=trim(output), status='replace')
  write(unit,'(i0,1x,es26.16e3,1x,es26.16e3,1x,i0)') fit%n, minimum, maximum, zero_count
  close(unit)
  call MPI_FINALIZE(ierr)
end program check_si_extended
