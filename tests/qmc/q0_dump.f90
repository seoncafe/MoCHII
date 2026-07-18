program q0_dump
!---------------------------------------------------------------------------
! Standalone Q0 driver for qmc_mod: dump the first N unscrambled Sobol codes
! (for the scipy cross-check) and the scrambled uniforms (for the dyadic
! stratification / replicate checks).
!
! Build (from the repo root):
!   mpiifort -O2 -module src -c src/qmc_mod.f90 -o src/qmc_mod.o   (already built by make)
!   mpiifort -O2 -Isrc tests/qmc/q0_dump.f90 src/qmc_mod.o -o tests/qmc/q0_dump.x
! Run:
!   tests/qmc/q0_dump.x <N> <seed>
! Writes q0_raw.txt (i c1 c2 c3, unscrambled 32-bit codes) and
! q0_scr_<seed>.txt (i u1 u2 u3, scrambled uniforms).
!---------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use qmc_mod, only : qmc_setup, qmc_uniforms, qmc_uniforms_raw, QMC_NDIM_USED
  implicit none
  integer, parameter :: wp = real64
  real(wp), parameter :: TWO32 = 4294967296.0_wp
  integer(int64) :: i, n
  integer :: seed, un
  character(len=64) :: arg, fname
  real(wp) :: u(QMC_NDIM_USED)
  integer(int64) :: c(QMC_NDIM_USED)

  n = 1024
  seed = 12345
  if (command_argument_count() >= 1) then
     call get_command_argument(1, arg);  read(arg,*) n
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, arg);  read(arg,*) seed
  end if

  call qmc_setup(seed)

  !--- unscrambled integer codes (seed independent).
  open(newunit=un, file='tests/qmc/q0_raw.txt', status='replace', action='write')
  do i = 0, n-1
     call qmc_uniforms_raw(i, u)
     c = nint(u*TWO32 - 0.5_wp, int64)
     write(un,'(i10,3(1x,i12))') i, c(1), c(2), c(3)
  end do
  close(un)

  !--- scrambled uniforms for this seed.
  write(fname,'(a,i0,a)') 'tests/qmc/q0_scr_', seed, '.txt'
  open(newunit=un, file=trim(fname), status='replace', action='write')
  do i = 0, n-1
     call qmc_uniforms(i, u)
     write(un,'(i10,3(1x,es22.15))') i, u(1), u(2), u(3)
  end do
  close(un)

  write(*,'(a,i0,a,i0)') 'q0_dump: wrote N=', n, ' points, seed=', seed
end program q0_dump
