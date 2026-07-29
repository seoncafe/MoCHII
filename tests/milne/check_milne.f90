program check_milne
!---------------------------------------------------------------------------
! Gate for milne_recomb_spectrum_mod.
!
! The sampler is checked against quadrature of the SAME Milne integrand it is
! built from, so this catches table, truncation, interpolation and inversion
! errors.  It does not re-derive the physics; milne_setup's alpha_1 report
! does that against independent recombination-rate data.
!
! Checks:
!   1. every sampled energy lies in [E_th, eion_max];
!   2. the sampled mean energy reproduces the quadrature energy-weighted mean
!      (the sampler draws from E n(E), the energy-luminosity density);
!   3. the sampled <1/E> equals 1/<E>_n, the identity that makes the emitted
!      photon rate equal n_e n_ion alpha_1 V for equal-energy packets;
!   4. selected quantiles match the quadrature CDF, so the shape is right and
!      not only the first moment;
!   5. the CDF is continuous across a temperature-table boundary.
!
! Build:  see tests/milne/run_check.sh
!---------------------------------------------------------------------------
  use define
  use photo_xsec, only : sigma_HI, sigma_HeI, sigma_HeII
  use milne_recomb_spectrum_mod, only : milne_setup, milne_sample_energy, &
                                        milne_energy_per_recomb, &
                                        MILNE_HI, MILNE_HEI, MILNE_HEII
  use mpi
  implicit none

  integer, parameter :: NS = 200000
  real(kind=wp), parameter :: TT(4) = [3.0e3_wp, 8.0e3_wp, 1.5e4_wp, 3.0e4_wp]
  character(len=6), parameter :: NAME(3) = &
     [character(len=6) :: 'H I', 'He I', 'He II']
  real(kind=wp) :: eth(3)
  integer :: ierr, ich, itt, i, nbad
  real(kind=wp) :: T, kT, u, E, s1, s2, sinv, emin, emax
  real(kind=wp) :: qn, qe, qmean_n, qmean_e, dev
  logical :: ok

  call MPI_INIT(ierr)
  par%eion_max = 100.0_wp
  par%diffuse_energy_model = 'milne'
  eth = [eth_HI, eth_HeI, eth_HeII]
  call milne_setup()
  ok = .true.

  do ich = 1, 3
     write(*,'(/,3a)') '=== ', trim(NAME(ich)), ' ground continuum ==='
     write(*,'(a8,5a13)') 'T [K]', '<E> samp', '<E> quad', &
                          '<E>_n quad', '<1/E>*<E>_n', 'q50 dev'
     do itt = 1, size(TT)
        T  = TT(itt)
        kT = kboltz_cgs*T/ev2erg
        call quad_moments(ich, T, qn, qe, qmean_n, qmean_e)

        s1 = 0.0_wp;  s2 = 0.0_wp;  sinv = 0.0_wp
        emin = huge(1.0_wp);  emax = -huge(1.0_wp)
        nbad = 0
        do i = 1, NS
           u = (real(i, wp) - 0.5_wp)/real(NS, wp)     ! stratified, low noise
           E = milne_sample_energy(ich, T, u)
           if (E < eth(ich) - 1.0e-9_wp .or. E > par%eion_max) nbad = nbad + 1
           s1   = s1   + E
           sinv = sinv + 1.0_wp/E
           emin = min(emin, E);  emax = max(emax, E)
        end do
        s1   = s1/real(NS, wp)
        sinv = sinv/real(NS, wp)

        !--- 1. range
        if (nbad > 0) then
           write(*,'(a,i0,a)') '  FAIL: ', nbad, ' samples outside the band'
           ok = .false.
        end if
        !--- 2. energy-weighted mean
        dev = abs(s1/qmean_e - 1.0_wp)
        if (dev > 2.0e-4_wp) then
           write(*,'(a,es10.3)') '  FAIL: sampled <E> deviates by ', dev
           ok = .false.
        end if
        !--- 3. photon-rate identity  <1/E>_{En} = 1/<E>_n
        dev = abs(sinv*qmean_n - 1.0_wp)
        if (dev > 5.0e-4_wp) then
           write(*,'(a,es10.3)') &
              '  FAIL: <1/E> <E>_n departs from 1 by ', dev
           ok = .false.
        end if
        !--- 4. median
        E = milne_sample_energy(ich, T, 0.5_wp)
        dev = abs(E - quantile_quad(ich, T, 0.5_wp))/kT
        if (dev > 5.0e-3_wp) then
           write(*,'(a,es10.3,a)') '  FAIL: median off by ', dev, ' kT'
           ok = .false.
        end if
        write(*,'(f8.0,4f13.6,es13.3)') T, s1, qmean_e, qmean_n, &
             sinv*qmean_n, dev

        !--- 5. the mean energy diffuse_build uses must be the number-weighted
        !--- mean of the same distribution (in band, so equal to <E>_n here)
        !--- 2e-5 admits the linear interpolation of <x> between temperature
        !--- nodes; what must hold tightly is the photon-rate identity above,
        !--- not the absolute accuracy of either mean on its own.
        dev = abs(milne_energy_per_recomb(ich, T)/qmean_n - 1.0_wp)
        if (dev > 2.0e-5_wp) then
           write(*,'(a,es10.3)') &
              '  FAIL: milne_energy_per_recomb vs quadrature <E>_n ', dev
           ok = .false.
        end if
     end do
  end do

  !--- 6. continuity in temperature across a table node
  write(*,'(/,a)') '=== continuity across a temperature node ==='
  do ich = 1, 3
     dev = 0.0_wp
     do i = 1, 200
        T = 1.0e4_wp*(1.0_wp + 1.0e-6_wp*real(i-100, wp))
        dev = max(dev, abs(milne_sample_energy(ich, T, 0.7_wp) &
                         - milne_sample_energy(ich, 1.0e4_wp, 0.7_wp)))
     end do
     write(*,'(3a,es10.3,a)') '  ', trim(NAME(ich)), &
        ': max |dE| over a 1e-4 relative T span = ', dev, ' eV'
     if (dev > 1.0e-3_wp) then
        write(*,'(a)') '  FAIL: discontinuous in temperature'
        ok = .false.
     end if
  end do

  write(*,'(/,a)') 'GATE (Milne ground-continuum sampler): ' // &
     merge('PASS', 'FAIL', ok)
  call MPI_FINALIZE(ierr)
  if (.not. ok) error stop 1

contains

  real(kind=wp) function sig(ich, E)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: E
    select case (ich)
    case (1);      sig = sigma_HI(E)
    case (2);      sig = sigma_HeI(E)
    case default;  sig = sigma_HeII(E)
    end select
  end function sig

  !--- fine-grid quadrature of the photon-number density n(E) and of E n(E)
  subroutine quad_moments(ich, T, qn, qe, mean_n, mean_e)
    integer,       intent(in)  :: ich
    real(kind=wp), intent(in)  :: T
    real(kind=wp), intent(out) :: qn, qe, mean_n, mean_e
    integer, parameter :: NQ = 400000
    real(kind=wp) :: kT, xhi, dx, x, E, f, q2
    integer :: i
    kT  = kboltz_cgs*T/ev2erg
    xhi = min(60.0_wp, (par%eion_max - eth(ich))/kT)
    dx  = xhi/real(NQ, wp)
    qn = 0.0_wp;  qe = 0.0_wp;  q2 = 0.0_wp
    do i = 1, NQ
       x = (real(i, wp) - 0.5_wp)*dx
       E = eth(ich) + x*kT
       f = E*E*sig(ich, E)*exp(-x)*dx
       qn = qn + f
       qe = qe + E*f
       q2 = q2 + E*E*f
    end do
    mean_n = qe/qn          ! photon-number-weighted mean energy
    mean_e = q2/qe          ! energy-weighted mean energy
  end subroutine quad_moments

  real(kind=wp) function quantile_quad(ich, T, p) result(Eq)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: T, p
    integer, parameter :: NQ = 400000
    real(kind=wp) :: kT, xhi, dx, x, E, f, tot, acc
    integer :: i
    kT  = kboltz_cgs*T/ev2erg
    xhi = min(60.0_wp, (par%eion_max - eth(ich))/kT)
    dx  = xhi/real(NQ, wp)
    tot = 0.0_wp
    do i = 1, NQ
       x = (real(i, wp) - 0.5_wp)*dx
       E = eth(ich) + x*kT
       tot = tot + E*(E*E*sig(ich, E)*exp(-x))*dx
    end do
    acc = 0.0_wp;  Eq = eth(ich)
    do i = 1, NQ
       x = (real(i, wp) - 0.5_wp)*dx
       E = eth(ich) + x*kT
       f = E*(E*E*sig(ich, E)*exp(-x))*dx
       acc = acc + f
       if (acc >= p*tot) then
          Eq = E;  return
       end if
    end do
  end function quantile_quad

end program check_milne
