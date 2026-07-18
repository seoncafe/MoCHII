program q2_dump
!---------------------------------------------------------------------------
! Stage-2 Q0 driver for qmc_mod: dump the first N unscrambled Sobol codes for
! ALL 12 embedded dimensions (scipy d=12 cross-check) and the scrambled
! uniforms for BOTH streams (stellar / diffuse) - for the dyadic-balance and
! stream-decorrelation checks.
!
! Build (from the repo root, after `make` has built src/qmc_mod.o):
!   mpiifort -O2 -Isrc tests/qmc/q2_dump.f90 src/qmc_mod.o -o tests/qmc/q2_dump.x
! Run:
!   tests/qmc/q2_dump.x <N> <seed>
! Writes q2_raw.txt (i c1..c12), q2_scr_<seed>.txt (i u1..u12, stellar) and
! q2_dif_<seed>.txt (i u1..u12, diffuse stream).
!---------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use qmc_mod, only : qmc_setup, qmc_uniforms, qmc_uniforms_stream, &
                      qmc_uniforms_raw, QMC_MAXDIM, &
                      QMC_STREAM_STELLAR, QMC_STREAM_DIFFUSE
  implicit none
  integer, parameter :: wp = real64
  real(wp), parameter :: TWO32 = 4294967296.0_wp
  integer(int64) :: i, n
  integer :: seed, un
  character(len=64) :: arg, fname
  real(wp) :: u(QMC_MAXDIM)
  integer(int64) :: c(QMC_MAXDIM)

  n = 1024
  seed = 12345
  if (command_argument_count() >= 1) then
     call get_command_argument(1, arg);  read(arg,*) n
  end if
  if (command_argument_count() >= 2) then
     call get_command_argument(2, arg);  read(arg,*) seed
  end if

  call qmc_setup(seed)

  !--- unscrambled integer codes for all 12 dims (seed independent).
  open(newunit=un, file='tests/qmc/q2_raw.txt', status='replace', action='write')
  do i = 0, n-1
     call qmc_uniforms_raw(i, u)
     c = nint(u*TWO32 - 0.5_wp, int64)
     write(un,'(i10,12(1x,i12))') i, c
  end do
  close(un)

  !--- scrambled uniforms, stellar stream.
  write(fname,'(a,i0,a)') 'tests/qmc/q2_scr_', seed, '.txt'
  open(newunit=un, file=trim(fname), status='replace', action='write')
  do i = 0, n-1
     call qmc_uniforms_stream(i, u, QMC_STREAM_STELLAR)
     write(un,'(i10,12(1x,es22.15))') i, u
  end do
  close(un)

  !--- scrambled uniforms, diffuse stream.
  write(fname,'(a,i0,a)') 'tests/qmc/q2_dif_', seed, '.txt'
  open(newunit=un, file=trim(fname), status='replace', action='write')
  do i = 0, n-1
     call qmc_uniforms_stream(i, u, QMC_STREAM_DIFFUSE)
     write(un,'(i10,12(1x,es22.15))') i, u
  end do
  close(un)

  write(*,'(a,i0,a,i0)') 'q2_dump: wrote N=', n, ' points (12 dims), seed=', seed
end program q2_dump
