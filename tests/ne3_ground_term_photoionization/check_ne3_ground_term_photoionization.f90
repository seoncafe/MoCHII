program check_ne3_ground_term_photoionization
!---------------------------------------------------------------------------
! Gate for the tabulated R-matrix ground-term cross section of Ne III
! (element_ne.txt PHTABLE, Nahar 2019 BPRM; generator
! tools/fitting/ground_term_photoionization_table.py) in the two places it
! enters: the transport and the metal ground continuum of the diffuse field
! (channel 6).  The references come from the published source alone
! (make_reference.py, which states how; reference.txt):
!
!  (a) species_cross_section at 4000 photon energies below E_join (1000 of
!      them in the first 2 eV, among the resonances) against the term
!      average of the three level tables interpolated each on its own grid.
!  (b) alpha_1(T) of channel 6 (milne_metal_alpha1, the node-stepped Milne
!      integral of the table at the common threshold ETH) against the Milne
!      relation on the three LEVEL tables, each at its own threshold, at
!      3e3-5e4 K.  The difference measures the common-threshold choice of
!      the table (<= 0.114 eV, weight 1/9 + 3/9).  Nahar has published no
!      Ne III recombination rate coefficients (make_reference.py), so this
!      Milne relation on her cross sections is the comparison available.
!  (c) The sampled channel-6 spectrum resolves the resonances: 2e6 packets
!      (stratified u) at 7.3e3, 1e4 and 1.37e4 K (between and on table rows)
!      histogrammed in 5 meV bins over the first 2 eV above threshold,
!      against the expected counts from the energy integral of each bin
!      (milne_metal_energy_integral, node-stepped); every bin with >= 100
!      expected counts within 5 Poisson sigma, and the bin-to-bin structure
!      reported (the largest/smallest bin ratio).  Sampling from an x-indexed
!      CDF interpolated between T rows fails this at 1e4 K by up to 435
!      sigma (30x the expected count 0.1 eV beside a resonance), which is
!      why the tabulated transitions sample from fixed energy blocks.
!
! Build and run: tests/ne3_ground_term_photoionization/run_check.sh
!---------------------------------------------------------------------------
  use define
  use species_mod, only : species_setup, species_cross_section, &
                          elem_name, elem_eth, species_transition, &
                          species_transition_count
  use milne_recomb_spectrum_mod, only : milne_setup, milne_metal_alpha1, &
                          milne_metal_sample_energy, milne_metal_energy_integral, &
                          milne_metal_band_coefficient
  use mpi
  implicit none
  real(kind=wp), parameter :: TOL_SIGMA = 1.0e-5_wp, TOL_ALPHA = 5.0e-3_wp
  integer,       parameter :: NS = 2000000, NBIN = 400
  real(kind=wp), parameter :: BIN = 5.0e-3_wp
  real(kind=wp), parameter :: TSAMP(3) = [7.3e3_wp, 1.0e4_wp, 1.37e4_wp]
  real(kind=wp) :: TS
  integer :: its
  character(len=256) :: line
  character(len=16)  :: key
  real(kind=wp), allocatable :: eref(:), sref(:), tref(:), aref(:)
  real(kind=wp) :: sig, heat, dev, worst, a1, eth, ptot, u, e, expct, obs, z, zmax
  real(kind=wp) :: emax_bin, emin_bin
  integer :: ierr, unit, n, na, i, k, ie, it, kne, nbig
  integer :: counts(NBIN)
  logical :: ok

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir = 'data/atomic'
  par%use_metals = .true.
  par%cooling_model = 'low_density'
  par%eion_min = 13.598_wp
  par%eion_max = 150.0_wp
  par%abund_Ne = 1.0e-4_wp
  call species_setup(1)
  call milne_setup()
  ok = .true.

  !--- the Ne III transition
  kne = 0
  do k = 1, species_transition_count()
     call species_transition(k, ie, it)
     if (trim(elem_name(ie)) == 'ne' .and. it == 3) kne = k
  end do
  if (kne == 0) then
     write(*,'(a)') 'FAIL: no Ne III transition';  error stop 1
  end if
  call species_transition(kne, ie, it)
  eth = elem_eth(ie, it)

  !--- reference file
  open(newunit=unit, file='tests/ne3_ground_term_photoionization/reference.txt', status='old')
  do
     read(unit,'(a)') line
     if (line(1:1) == '#') cycle
     read(line,*) key, na
     if (trim(key) == 'ALPHA1') exit
  end do
  allocate(tref(na), aref(na))
  do i = 1, na
     read(unit,*) tref(i), aref(i)
  end do
  read(unit,*) key, n
  allocate(eref(n), sref(n))
  do i = 1, n
     read(unit,*) eref(i), sref(i)
  end do
  close(unit)

  !--- (a)
  worst = 0.0_wp
  do i = 1, n
     call species_cross_section(ie, it, eref(i), sig, heat)
     dev = abs(sig/max(sref(i), 1.0e-40_wp) - 1.0_wp)
     worst = max(worst, dev)
     if (dev > TOL_SIGMA) then
        ok = .false.
        write(*,'(a,f14.8,2es14.6,es10.2)') '  FAIL sigma ', eref(i), sig, sref(i), dev
     end if
  end do
  write(*,'(a,i0,a,es10.3,a,es9.1,a)') '(a) sigma at ', n, ' energies vs the level tables: worst ', &
     worst, ' (tol ', TOL_SIGMA, ')'

  !--- (b)
  do i = 1, na
     a1  = milne_metal_alpha1(kne, tref(i))
     dev = a1/aref(i) - 1.0_wp
     write(*,'(a,f9.0,a,es12.5,a,es12.5,a,es10.2)') '(b) T =', tref(i), '  alpha_1 ', a1, &
        '  levels at own thresholds ', aref(i), '  dev ', dev
     if (abs(dev) > TOL_ALPHA) ok = .false.
  end do

  !--- (c)
  do its = 1, size(TSAMP)
  TS = TSAMP(its)
  counts = 0
  emin_bin = eth;  emax_bin = eth + BIN*NBIN
  do i = 1, NS
     u = (real(i, wp) - 0.5_wp)/real(NS, wp)
     e = milne_metal_sample_energy(kne, TS, u)
     k = int((e - emin_bin)/BIN) + 1
     if (k >= 1 .and. k <= NBIN) counts(k) = counts(k) + 1
  end do
  ptot = milne_metal_band_coefficient(kne, TS)
  zmax = 0.0_wp;  nbig = 0
  do k = 1, NBIN
     expct = real(NS, wp)*milne_metal_energy_integral(kne, TS, emin_bin + BIN*(k-1), &
             emin_bin + BIN*k)/ptot
     if (expct < 100.0_wp) cycle
     nbig = nbig + 1
     obs = real(counts(k), wp)
     z = abs(obs - expct)/sqrt(expct)
     zmax = max(zmax, z)
     if (z > 5.0_wp) then
        ok = .false.
        write(*,'(a,f10.4,2f12.1,f8.2)') '  FAIL bin ', emin_bin + BIN*(k-0.5_wp), obs, expct, z
     end if
  end do
  write(*,'(a,f8.0,a,i0,a,f6.2,a,f8.2,a)') '(c) sampled spectrum, 5 meV bins over 2 eV, T =', &
     TS, ' K: ', nbig, ' bins, largest |obs - exp|/sqrt(exp) = ', zmax, &
     '; largest/smallest bin count ', &
     real(maxval(counts),wp)/max(real(minval(counts, mask=counts > 0),wp), 1.0_wp), &
     ' (the resonances)'
  end do

  if (ok) then
     write(*,'(/,a)') 'GATE (Ne III tabulated ground-term cross section: transport, ' // &
        'alpha_1, sampled spectrum): PASS'
  else
     write(*,'(/,a)') 'GATE (Ne III tabulated ground-term cross section): FAIL'
     call MPI_FINALIZE(ierr);  error stop 1
  end if
  call MPI_FINALIZE(ierr)
end program check_ne3_ground_term_photoionization
