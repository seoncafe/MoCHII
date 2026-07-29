program check_hei_2ph
!---------------------------------------------------------------------------
! Gate for hei_twophoton_mod: the He I 2^1S two-photon in-band spectrum.
!
! Checks, against independent quadrature of the same data/HeI2phot.dat table:
!   1. every sampled energy lies in [E_th(H I), 20.62 eV];
!   2. the sampled mean reproduces the ENERGY-weighted mean, because the
!      sampler draws from E n(E) for equal-energy packets;
!   3. <1/E> equals 1/<E>_n, the identity that makes the represented photon
!      rate equal the branch's decay rate;
!   4. quantiles match the quadrature CDF, so the shape is right and not only
!      the first moment;
!   5. the reported photons and energy per decay match quadrature;
!   6. the tabulated shape is symmetric about y = 1/2 and integrates to two
!      photons per decay, which is what makes the factor 2 in the module
!      correct rather than an arbitrary choice.
!---------------------------------------------------------------------------
  use define
  use hei_twophoton_mod, only : hei_2ph_setup, hei_2ph_sample_energy, &
                                hei_2ph_photons_per_decay, &
                                hei_2ph_energy_per_decay
  use mpi
  implicit none

  integer, parameter :: NS = 200000, NQ = 400000
  real(kind=wp), parameter :: E_TOT = 20.62_wp, A2G = 51.3_wp
  real(kind=wp) :: ty(256), ta(256)
  integer :: nt
  integer :: ierr, i, nbad
  real(kind=wp) :: u, E, s1, sinv, tot, qn, qe, q2, mean_n, mean_e, dev
  real(kind=wp) :: ya, yb, ym, ylo, yhi, dy, y
  logical :: ok

  call MPI_INIT(ierr)
  par%eion_max = 100.0_wp
  par%data_dir = '../../data'
  call hei_2ph_setup()
  call read_table()
  ok = .true.

  !--- 6. normalization and symmetry of the source table
  tot = 0.0_wp
  do i = 1, nt-1
     tot = tot + 0.5_wp*(ta(i) + ta(i+1))*(ty(i+1) - ty(i))
  end do
  write(*,'(a,f10.4,a,f8.5,a)') ' table integral = ', tot, &
     '  = ', tot/A2G, ' x A_2gamma  (expect 2 photons per decay)'
  if (abs(tot/A2G - 2.0_wp) > 0.02_wp) then
     write(*,'(a)') '  FAIL: table does not integrate to two photons per decay'
     ok = .false.
  end if
  dev = 0.0_wp
  do i = 1, nt
     dev = max(dev, abs(ta(i) - ta(nt+1-i))/max(maxval(ta), tinest))
  end do
  write(*,'(a,es10.3)') ' symmetry about y = 1/2, max relative asymmetry = ', dev
  if (dev > 5.0e-3_wp) then
     write(*,'(a)') '  FAIL: shape is not symmetric'
     ok = .false.
  end if

  !--- quadrature moments over the in-band window
  ylo = eth_HI/E_TOT;  yhi = 1.0_wp
  dy = (yhi - ylo)/real(NQ, wp)
  qn = 0.0_wp;  qe = 0.0_wp;  q2 = 0.0_wp
  do i = 1, NQ
     y = ylo + (real(i, wp) - 0.5_wp)*dy
     qn = qn + interp(y)*dy
     qe = qe + interp(y)*y*E_TOT*dy
     q2 = q2 + interp(y)*(y*E_TOT)**2*dy
  end do
  mean_n = qe/qn                 ! number-weighted mean energy
  mean_e = q2/qe                 ! energy-weighted mean energy

  !--- 5. reported per-decay quantities
  write(*,'(/,a)') ' per 2^1S decay:'
  write(*,'(a,f10.6,a,f10.6)') '   photons  module ', &
     hei_2ph_photons_per_decay(), '   quadrature ', 2.0_wp*qn/tot
  write(*,'(a,f10.6,a,f10.6)') '   energy   module ', &
     hei_2ph_energy_per_decay(), '   quadrature ', 2.0_wp*qe/tot
  dev = abs(hei_2ph_photons_per_decay()/(2.0_wp*qn/tot) - 1.0_wp)
  if (dev > 1.0e-4_wp) then
     write(*,'(a,es10.3)') '  FAIL: photons per decay off by ', ok_dev(dev)
     ok = .false.
  end if
  dev = abs(hei_2ph_energy_per_decay()/(2.0_wp*qe/tot) - 1.0_wp)
  if (dev > 1.0e-4_wp) then
     write(*,'(a,es10.3)') '  FAIL: energy per decay off by ', ok_dev(dev)
     ok = .false.
  end if

  !--- 1-3. sampled moments, stratified so the estimate is not noise-limited
  s1 = 0.0_wp;  sinv = 0.0_wp;  nbad = 0
  do i = 1, NS
     u = (real(i, wp) - 0.5_wp)/real(NS, wp)
     E = hei_2ph_sample_energy(u)
     if (E < eth_HI - 1.0e-9_wp .or. E > E_TOT + 1.0e-9_wp) nbad = nbad + 1
     s1 = s1 + E;  sinv = sinv + 1.0_wp/E
  end do
  s1 = s1/real(NS, wp);  sinv = sinv/real(NS, wp)
  write(*,'(/,a)') ' sampler:'
  write(*,'(a,i0)')            '   samples outside the band : ', nbad
  write(*,'(a,f10.6,a,f10.6)') '   <E>   sampled ', s1, '   quadrature ', mean_e
  write(*,'(a,f12.8)')         '   <1/E> <E>_n              ', sinv*mean_n
  if (nbad > 0) then
     write(*,'(a)') '  FAIL: samples outside [E_th, E_tot]'
     ok = .false.
  end if
  if (abs(s1/mean_e - 1.0_wp) > 2.0e-4_wp) then
     write(*,'(a)') '  FAIL: sampled mean is not the energy-weighted mean'
     ok = .false.
  end if
  if (abs(sinv*mean_n - 1.0_wp) > 5.0e-4_wp) then
     write(*,'(a)') '  FAIL: photon-rate identity <1/E> = 1/<E>_n broken'
     ok = .false.
  end if

  !--- 4. quantiles
  write(*,'(/,a)') ' quantiles [eV]:'
  do i = 1, 3
     u = 0.25_wp*real(i, wp)
     E = hei_2ph_sample_energy(u)
     dev = abs(E - quantile(u))
     write(*,'(a,f4.2,a,f10.5,a,f10.5,a,es9.2)') '   q', u, ' sampled ', E, &
        '  quadrature ', quantile(u), '  |d| ', dev
     if (dev > 5.0e-3_wp) then
        write(*,'(a)') '  FAIL: quantile mismatch'
        ok = .false.
     end if
  end do

  write(*,'(/,a)') 'GATE (He I 2^1S two-photon spectrum): ' // &
     merge('PASS', 'FAIL', ok)
  call MPI_FINALIZE(ierr)
  if (.not. ok) error stop 1

contains

  real(kind=wp) function ok_dev(d) result(r)
    real(kind=wp), intent(in) :: d
    r = d
  end function ok_dev

  subroutine read_table()
    integer :: unit, ios
    real(kind=wp) :: a, b
    open(newunit=unit, file=trim(par%data_dir)//'/HeI2phot.dat', status='old')
    nt = 0
    do
       read(unit,*,iostat=ios) a, b
       if (ios /= 0) exit
       nt = nt + 1;  ty(nt) = a;  ta(nt) = b
    end do
    close(unit)
  end subroutine read_table

  real(kind=wp) function interp(y) result(a)
    real(kind=wp), intent(in) :: y
    integer :: i
    a = 0.0_wp
    if (y <= ty(1) .or. y >= ty(nt)) return
    i = 1
    do while (i < nt-1 .and. ty(i+1) < y)
       i = i + 1
    end do
    a = ta(i) + (ta(i+1) - ta(i))*(y - ty(i))/(ty(i+1) - ty(i))
  end function interp

  !--- quantile of the ENERGY-weighted density, the one the sampler uses
  real(kind=wp) function quantile(p) result(Eq)
    real(kind=wp), intent(in) :: p
    real(kind=wp) :: acc, total, yy, ddy
    integer :: i
    ddy = (1.0_wp - eth_HI/E_TOT)/real(NQ, wp)
    total = 0.0_wp
    do i = 1, NQ
       yy = eth_HI/E_TOT + (real(i, wp) - 0.5_wp)*ddy
       total = total + interp(yy)*yy*ddy
    end do
    acc = 0.0_wp;  Eq = eth_HI
    do i = 1, NQ
       yy = eth_HI/E_TOT + (real(i, wp) - 0.5_wp)*ddy
       acc = acc + interp(yy)*yy*ddy
       if (acc >= p*total) then
          Eq = yy*E_TOT;  return
       end if
    end do
  end function quantile

end program check_hei_2ph
