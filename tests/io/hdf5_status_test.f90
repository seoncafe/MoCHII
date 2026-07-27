program hdf5_status_test
  use, intrinsic :: iso_fortran_env, only : real64
  use hdf5io_mod
  use iofile_mod, only : io_detect_format, IO_FMT_HDF5, IO_FMT_FITS
  implicit none
  type(hdf5_state_type) :: file
  real(real64) :: values(3), restored(3)
  character(len=512) :: good_path, bad_path
  integer :: status

  call get_command_argument(1, good_path)
  call get_command_argument(2, bad_path)
  if (len_trim(good_path) == 0 .or. len_trim(bad_path) == 0) &
     error stop 'usage: hdf5_status_test GOOD_PATH BAD_PATH'

  if (io_detect_format('result.H5') /= IO_FMT_HDF5) &
     error stop 'uppercase .H5 was not detected'
  if (io_detect_format('result.HDF5') /= IO_FMT_HDF5) &
     error stop 'uppercase .HDF5 was not detected'
  if (io_detect_format('result.fits') /= IO_FMT_FITS) &
     error stop 'FITS detection regressed'

  values = [1.0_real64, 2.0_real64, 3.0_real64]

  status = 0
  call hdf5_open_new(file, trim(good_path), status)
  call hdf5_append_1D_real64(file, values, status)
  call hdf5_close(file, status)
  if (status /= 0) error stop 'normal HDF5 write failed'

  status = 0
  restored = 0.0_real64
  call hdf5_open_old(file, trim(good_path), status)
  call hdf5_read_1D_real64(file, restored, status)
  call hdf5_close(file, status)
  if (status /= 0) error stop 'HDF5 readback failed'
  if (any(restored /= values)) error stop 'HDF5 readback mismatch'

  status = 0
  call hdf5_open_new(file, trim(bad_path), status)
  if (status == 0) error stop 'HDF5 create failure was not reported'

  status = 0
  call hdf5_open_new(file, trim(good_path), status)
  if (status /= 0) error stop 'status-preservation setup failed'
  status = -777
  call hdf5_close(file, status)
  if (status /= -777) error stop 'cleanup overwrote the first error'

  print '(a)', 'HDF5 status regression: PASS'
end program hdf5_status_test
