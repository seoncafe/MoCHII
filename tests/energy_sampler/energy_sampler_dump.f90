program energy_sampler_dump
  use, intrinsic :: iso_fortran_env, only : int64, real64
  use energy_sampler_mod, only : energy_sampler_type, build_tabulated_sampler, &
       build_tabulated_band_sampler, build_planck_sampler, sample_energy_cdf, &
       energy_cdf_value, sample_full_planck_bc
  use qmc_mod, only : qmc_setup, qmc_uniforms, qmc_uniforms_stream, &
       QMC_STREAM_DIFFUSE
  implicit none

  integer, parameter :: wp = real64
  integer, parameter :: nrandom = 40000, nbc = 200000
  integer, parameter :: nqmc = 32768, ntab = 4096
  type(energy_sampler_type) :: p20, p40, constant, linear, broken, clipped
  real(wp) :: thresholds(4), etab(4), ftab(4), u(9), ur(5)
  real(wp) :: e20, e40
  character(len=512) :: outdir, fname
  integer :: i, un, nseed
  integer, allocatable :: seed(:)

  call get_command_argument(1, outdir)
  if (len_trim(outdir) == 0) error stop 'usage: energy_sampler_dump.x OUTDIR'

  thresholds = [13.598_wp, 24.383_wp, 24.587_wp, 54.418_wp]
  call build_planck_sampler(p20, 2.0e4_wp, 13.598_wp, 100.0_wp, 1.0e-8_wp, thresholds)
  call build_planck_sampler(p40, 4.0e4_wp, 13.598_wp, 100.0_wp, 1.0e-8_wp, thresholds)

  etab = [13.598_wp, 24.383_wp, 24.587_wp, 100.0_wp]
  ftab = 1.0_wp
  call build_tabulated_sampler(constant, etab, ftab)
  ftab = 0.25_wp + 0.03_wp*(etab-etab(1))
  call build_tabulated_sampler(linear, etab, ftab)
  ftab = [0.2_wp, 3.0_wp, 0.4_wp, 0.05_wp]
  call build_tabulated_sampler(broken, etab, ftab)
  call build_tabulated_band_sampler(clipped, [10.0_wp,20.0_wp,110.0_wp], &
       [10.0_wp,20.0_wp,110.0_wp], 13.598_wp, 100.0_wp)
  if (abs(clipped%total_integral - 0.5_wp*(100.0_wp**2-13.598_wp**2)) > &
      2.0e-12_wp*clipped%total_integral) &
     error stop 'band-clipped tabulated integral failed'

  call random_seed(size=nseed)
  allocate(seed(nseed))
  do i = 1, nseed
     seed(i) = 7919 + 104729*i
  end do
  call random_seed(put=seed)
  deallocate(seed)
  call qmc_setup(20260728)

  fname = trim(outdir)//'/random.dat'
  open(newunit=un, file=trim(fname), status='replace', action='write')
  do i = 1, nrandom
     call random_number(ur)
     write(un,'(5(es24.16,1x))') sample_energy_cdf(p20,ur(1)), &
          sample_energy_cdf(p40,ur(2)), sample_full_planck_bc(2.0e4_wp,ur), &
          ur(1), ur(2)
  end do
  close(un)

  fname = trim(outdir)//'/full_planck.dat'
  open(newunit=un, file=trim(fname), status='replace', action='write')
  do i = 1, nbc
     call random_number(ur)
     write(un,'(es24.16)') sample_full_planck_bc(2.0e4_wp,ur)
  end do
  close(un)

  fname = trim(outdir)//'/qmc.dat'
  open(newunit=un, file=trim(fname), status='replace', action='write')
  do i = 0, nqmc-1
     call qmc_uniforms(int(i,int64), u)
     e20 = sample_energy_cdf(p20,u(1))
     e40 = sample_energy_cdf(p40,u(1))
     write(un,'(9(es24.16,1x))') u(1), e20, e40, u(2), &
          sample_energy_cdf(p40,u(2)), &
          merge(sample_energy_cdf(p20,u(2)),sample_energy_cdf(p40,u(2)),u(1)<0.3_wp), &
          merge(1.0_wp,2.0_wp,u(1)<0.3_wp), &
          energy_cdf_value(p20,e20), energy_cdf_value(p40,e40)
  end do
  close(un)

  fname = trim(outdir)//'/diffuse.dat'
  open(newunit=un, file=trim(fname), status='replace', action='write')
  do i = 0, nqmc-1
     call qmc_uniforms_stream(int(i,int64), u, QMC_STREAM_DIFFUSE)
     write(un,'(3(es24.16,1x))') u(8), sample_energy_cdf(p40,u(8)), u(9)
  end do
  close(un)

  fname = trim(outdir)//'/tabulated.dat'
  open(newunit=un, file=trim(fname), status='replace', action='write')
  do i = 1, ntab
     u(1) = (real(i,wp)-0.5_wp)/real(ntab,wp)
     if (abs(energy_cdf_value(clipped,sample_energy_cdf(clipped,u(1)))-u(1)) > &
         2.0e-14_wp) error stop 'band-clipped tabulated round trip failed'
     write(un,'(7(es24.16,1x))') u(1), sample_energy_cdf(constant,u(1)), &
          sample_energy_cdf(linear,u(1)), sample_energy_cdf(broken,u(1)), &
          energy_cdf_value(constant,sample_energy_cdf(constant,u(1))), &
          energy_cdf_value(linear,sample_energy_cdf(linear,u(1))), &
          energy_cdf_value(broken,sample_energy_cdf(broken,u(1)))
  end do
  close(un)
end program energy_sampler_dump
