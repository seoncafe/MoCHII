program check_cooling_fit
  use mpi
  use define, only: wp
  use cooling_mod, only: cooling_fit_type, cooling_load, cooling_eval
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(cooling_fit_type) :: fit
  character(len=1024) :: input, output
  real(kind=wp) :: temperature, coefficient
  integer :: i, unit, ierr
  call MPI_INIT(ierr)
  call get_command_argument(1, input)
  call get_command_argument(2, output)
  call cooling_load(trim(input), fit)
  open(newunit=unit, file=trim(output), status='replace')
  do i = 0, 800
     temperature = 10.0_wp**(3.0_wp+2.0_wp*real(i,wp)/800.0_wp)
     coefficient = cooling_eval(fit, temperature)
     if (.not. ieee_is_finite(coefficient) .or. coefficient < 0.0_wp) stop 1
     write(unit,'(2(es25.16e3,1x))') temperature, coefficient
  end do
  close(unit)
  call MPI_FINALIZE(ierr)
end program check_cooling_fit
