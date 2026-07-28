program energy_sampler_mpi
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use mpi
  use energy_sampler_mod, only : energy_sampler_type, build_planck_sampler, &
       sample_energy_cdf
  use qmc_mod, only : qmc_setup, qmc_uniforms
  implicit none

  integer, parameter :: wp = real64, ntotal = 4096
  type(energy_sampler_type) :: sampler
  real(wp), allocatable :: local_energy(:), all_energy(:)
  real(wp) :: u(3), thresholds(4)
  integer, allocatable :: counts(:), displs(:)
  integer :: ierr, rank, nproc, nlocal, first, i, un
  character(len=512) :: outfile

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  call get_command_argument(1, outfile)
  if (len_trim(outfile) == 0) call MPI_ABORT(MPI_COMM_WORLD, 2, ierr)

  allocate(counts(nproc), displs(nproc))
  do i = 0, nproc - 1
     counts(i+1) = ntotal/nproc + merge(1,0,i < mod(ntotal,nproc))
     displs(i+1) = i*(ntotal/nproc) + min(i,mod(ntotal,nproc))
  end do
  nlocal = counts(rank+1)
  first = displs(rank+1)
  allocate(local_energy(nlocal))
  if (rank == 0) then
     allocate(all_energy(ntotal))
  else
     allocate(all_energy(1))
  end if

  thresholds = [13.598_wp, 24.383_wp, 24.587_wp, 54.418_wp]
  call build_planck_sampler(sampler, 4.0e4_wp, 13.598_wp, 100.0_wp, &
       1.0e-8_wp, thresholds)
  call qmc_setup(20260728)
  do i = 1, nlocal
     call qmc_uniforms(int(first+i-1,int64), u)
     local_energy(i) = sample_energy_cdf(sampler,u(1))
  end do

  call MPI_GATHERV(local_energy, nlocal, MPI_DOUBLE_PRECISION, all_energy, &
       counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (rank == 0) then
     open(newunit=un, file=trim(outfile), access='stream', form='unformatted', &
          status='replace', action='write')
     write(un) all_energy
     close(un)
  end if
  call MPI_FINALIZE(ierr)
end program energy_sampler_mpi
