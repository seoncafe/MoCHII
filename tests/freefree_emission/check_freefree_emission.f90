program check_freefree_emission
!---------------------------------------------------------------------------
! Gate for the free-free emission of the code: the spectral coefficient of
! the written continuum and of the transported free-free packets
! (recomb_continuum_mod::gamma_ff, diffuse_mod channel 7) and the cooling
! (cooling_mod, species_mod::metal_freefree).
!
!  (a) Closure: the frequency integral of the spectral coefficient equals
!      the cooling coefficient,
!          Int_0^inf gamma_ff(nu, T, Z) 1e-40 dnu = FF_TOTAL Z^2 gbar_ff(T, Z) sqrt(T),
!      for Z = 1..6, 3e3-1e5 K (independent quadrature: 4000 log-spaced
!      points in u = h nu/kT over 1e-9..300 with the trapezoid rule in ln u).
!      The gas loses exactly the energy the spectrum carries.
!  (b) The coefficient's own exponent: gamma_ff / [FF_SPECTRAL Z^2 g_ff
!      T^(-1/2) exp(-h nu/kT)] = 1 with h nu = nu x RYD_ERG (the frequency
!      unit of the tables is the infinite-mass Rydberg, 13.6057 eV).
!  (c) Band coefficient of the diffuse channel against a direct quadrature
!      of gamma_ff over [eion_min, eion_max], at table nodes and between.
!  (d) The packets: equal-energy packets drawn by
!      freefree_band_sample_energy (stratified u) represent the photon rate
!      G <1/E>, which must equal Int_band gamma_ff/(h nu) dnu; all sampled
!      energies inside the band.
!  (e) Reported: the band share of the free-free emission, G/(FF_TOTAL Z^2
!      gbar sqrt(T)), the fraction the transported channel carries.
!
! Build and run: tests/freefree_emission/run_check.sh
!---------------------------------------------------------------------------
  use define
  use free_free_mod,        only : free_free_gaunt, FF_SPECTRAL, FF_TOTAL
  use recomb_continuum_mod, only : gamma_ff, freefree_band_setup, RYD_ERG, &
                                   freefree_band_coefficient, &
                                   freefree_band_sample_energy, FF_ZMAX
  use cooling_mod,          only : cooling_setup, gbar_ff
  use mpi
  implicit none
  real(kind=wp), parameter :: RYD_EV = 13.605693_wp
  real(kind=wp), parameter :: TT(7) = [3.0e3_wp, 5.0e3_wp, 8.0e3_wp, 1.0e4_wp, &
                                       1.37e4_wp, 2.0e4_wp, 5.0e4_wp]
  real(kind=wp), parameter :: TCL(6) = [3.0e3_wp, 7.0e3_wp, 1.0e4_wp, 2.0e4_wp, &
                                        5.0e4_wp, 1.0e5_wp]
  real(kind=wp), parameter :: TOL_CLOSURE = 1.0e-5_wp, TOL_EXP = 1.0e-12_wp
  real(kind=wp), parameter :: TOL_BAND = 2.0e-4_wp, TOL_RATE = 2.0e-4_wp
  integer, parameter :: NU_Q = 4000, NS = 100000, NB = 20000
  real(kind=wp) :: T, kT_eV, lu, u, du, nu, acc, gprev, gcur, uprev, ref, dev
  real(kind=wp) :: worst_cl, worst_exp, worst_band, worst_rate, g, gband, gdir
  real(kind=wp) :: rate, sinv, e, ea, eb, em, j0
  integer :: ierr, iz, it, i
  logical :: ok

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir = 'data/atomic'
  par%data_dir   = 'data'
  par%cooling_model = 'low_density'
  par%eion_min = 13.598_wp
  par%eion_max = 100.0_wp
  call cooling_setup()
  call freefree_band_setup()
  ok = .true.

  !--- (a) closure against the cooling coefficient
  worst_cl = 0.0_wp
  write(*,'(a)') '(a) Int gamma_ff dnu / (FF_TOTAL Z^2 gbar sqrt(T)) - 1'
  do iz = 1, FF_ZMAX
     do it = 1, size(TCL)
        T = TCL(it);  kT_eV = kboltz_cgs*T/ev2erg
        acc = 0.0_wp
        do i = 0, NU_Q
           lu = log(1.0e-9_wp) + (log(300.0_wp) - log(1.0e-9_wp))*real(i,wp)/real(NU_Q,wp)
           u  = exp(lu)
           nu = u*kT_eV/RYD_EV
           !--- Int gamma dnu = Int gamma (kT/h) u d ln u
           gcur = gamma_ff(nu, T, iz)*1.0e-40_wp*u*kboltz_cgs*T/h_planck_cgs
           if (i > 0) acc = acc + 0.5_wp*(gprev + gcur)*(lu - uprev)
           gprev = gcur;  uprev = lu
        end do
        ref = FF_TOTAL*real(iz*iz,wp)*gbar_ff(T, iz)*sqrt(T)
        dev = acc/ref - 1.0_wp
        worst_cl = max(worst_cl, abs(dev))
        write(*,'(a,i2,a,f8.0,a,es11.3)') '    Z =', iz, '  T =', T, '  dev =', dev
     end do
  end do
  if (worst_cl > TOL_CLOSURE) ok = .false.
  write(*,'(a,es10.3,a,es9.1,a)') '(a) worst closure ', worst_cl, ' (tol ', TOL_CLOSURE, ')'

  !--- (b) the exponent
  worst_exp = 0.0_wp
  do iz = 1, 2
     do it = 1, size(TT)
        T = TT(it)
        do i = 1, 20
           nu = 0.05_wp*real(i,wp)
           u  = nu*RYD_ERG/(kboltz_cgs*T)
           ref = FF_SPECTRAL*real(iz*iz,wp)*free_free_gaunt(real(iz,wp), T, u) &
                 *exp(-u)/sqrt(T)*1.0e40_wp
           if (ref <= 0.0_wp) cycle
           worst_exp = max(worst_exp, abs(gamma_ff(nu, T, iz)/ref - 1.0_wp))
        end do
     end do
  end do
  if (worst_exp > TOL_EXP) ok = .false.
  write(*,'(a,es10.3)') '(b) gamma_ff against FF_SPECTRAL Z^2 g exp(-h nu/kT)/sqrt(T): worst ', &
     worst_exp

  !--- (c), (d)
  worst_band = 0.0_wp;  worst_rate = 0.0_wp
  do iz = 1, FF_ZMAX
     do it = 1, size(TT)
        T = TT(it);  kT_eV = kboltz_cgs*T/ev2erg
        gband = freefree_band_coefficient(iz, T)
        !--- direct: midpoint slices in E over the band (up to 60 kT)
        gdir = 0.0_wp;  rate = 0.0_wp
        ea = par%eion_min
        eb = min(par%eion_max, par%eion_min + 60.0_wp*kT_eV)
        do i = 1, NB
           e  = ea + (eb - ea)*(real(i,wp) - 0.5_wp)/real(NB,wp)
           j0 = gamma_ff(e/RYD_EV, T, iz)*1.0e-40_wp*(eb - ea)/real(NB,wp) &
                *ev2erg/h_planck_cgs
           gdir = gdir + j0
           rate = rate + j0/(e*ev2erg)
        end do
        if (gdir <= 0.0_wp) cycle
        dev = abs(gband/gdir - 1.0_wp)
        worst_band = max(worst_band, dev)
        if (dev > TOL_BAND) write(*,'(a,i2,f9.0,2es14.6,es10.2)') '  FAIL band', &
           iz, T, gband, gdir, dev
        sinv = 0.0_wp
        do i = 1, NS
           e = freefree_band_sample_energy(iz, T, (real(i,wp) - 0.5_wp)/real(NS,wp))
           if (e < par%eion_min .or. e > par%eion_max) then
              ok = .false.
              write(*,'(a,i2,f9.0,es14.6)') '  FAIL sampled energy outside the band', iz, T, e
           end if
           sinv = sinv + 1.0_wp/(e*ev2erg)
        end do
        sinv = sinv/real(NS,wp)
        dev = abs(gband*sinv/rate - 1.0_wp)
        worst_rate = max(worst_rate, dev)
        if (dev > TOL_RATE) write(*,'(a,i2,f9.0,2es14.6,es10.2)') '  FAIL photon rate', &
           iz, T, gband*sinv, rate, dev
     end do
  end do
  if (worst_band > TOL_BAND .or. worst_rate > TOL_RATE) ok = .false.
  write(*,'(a,es10.3,a,es9.1,a)') '(c) band coefficient vs direct quadrature: worst ', &
     worst_band, ' (tol ', TOL_BAND, ')'
  write(*,'(a,es10.3,a,es9.1,a)') '(d) packets'' photon rate vs Int_band j/E: worst ', &
     worst_rate, ' (tol ', TOL_RATE, ')'

  !--- (e) band share of the free-free emission
  write(*,'(/,a)') '(e) share of the free-free emission inside [eion_min, eion_max]'
  do it = 1, size(TT)
     T = TT(it)
     write(*,'(a,f8.0,6es11.3)') '    T =', T, &
        (freefree_band_coefficient(iz, T)/(FF_TOTAL*real(iz*iz,wp)*gbar_ff(T, iz) &
         *sqrt(T)), iz = 1, FF_ZMAX)
  end do

  if (ok) then
     write(*,'(/,a)') 'GATE (free-free: closure, coefficient, band table, packets): PASS'
  else
     write(*,'(/,a)') 'GATE (free-free): FAIL'
     call MPI_FINALIZE(ierr);  error stop 1
  end if
  call MPI_FINALIZE(ierr)
end program check_freefree_emission
