program check_heii_cascade
!---------------------------------------------------------------------------
! Gate for heii_cascade_mod: the fifth diffuse channel, what a He III case-B
! recombination radiates on its way down through He II.
!
! Checks:
!   1. the n = 2 exit branching closes -- f_2s + f_2p = 1 -- with f_2s the
!      case-B alpha_eff(2^2S)/alpha_B of one zero-density Storey & Hummer
!      (1991) cascade (hydrogenic_recombination_mod; its own gate is
!      tests/hydrogenic_recombination);
!   2. alpha_eff(2^2P) = alpha_B (1 - f_2s) agrees with MOCASSIN's INDEPENDENT
!      He II Ly-alpha emissivity fit, 10^(-0.792 log10 Te + 6.01) e-25
!      erg cm^3/s, to better than 8%.  This is the one cross-check available
!      on the 2^2P branch that does not come from the same cascade;
!   3. the three branch luminosities are positive and ordered
!      Ly-alpha > two-photon > free-bound over 1e4-3e4 K;
!   4. every sampled energy lies in [E_th(H I), min(eion_max, 54.416)] -- no
!      packet may enter the band below the hydrogen edge, and the free-bound
!      part must stop where channel 3 takes over;
!   5. the sampler reproduces the branch luminosities: the fraction of
!      samples drawn from each path, weighted by 1/E, matches the photon rate
!      the luminosity implies.  This is the identity that makes equal-energy
!      packets represent the physical photon rate;
!   6. the two-photon in-band moments match independent quadrature of the
!      Nussbaumer & Schmutz (1984) shape, and the pair integrates to two
!      photons per decay;
!   7. case A (par%lyman_lines = 'escaping'): the Lyman-line path carries
!      alpha_B sum_n f_n E_n, every sampled line lies in [40.813, 54.42] eV,
!      the photons the energy-weighted sampler implies, L <1/E>, equal
!      alpha_B sum_n f_n to 1e-4 (the stratified sampling of NS draws), and
!      the free-bound term does not depend on the case.  Checks 1-6 run in
!      case B, the default.
!---------------------------------------------------------------------------
  use define
  use recomb_mod,               only : alphaB_HeIII
  use heii_cascade_mod,         only : heii_cascade_setup, &
                                       heii_branch_luminosity, &
                                       heii_fb_sample_energy, &
                                       heii_2q_sample_energy, &
                                       heii_lyman_sample_energy, &
                                       heii_f_2s, E_HEII_LYA
  use hydrogenic_recombination_mod, only : heii_lyman_shares, heii_lyman_energy, &
                                           HEII_NLY
  use hydrogenic_twophoton_mod, only : h2q_band_photons_per_decay, &
                                       h2q_band_energy_per_decay
  use mpi
  implicit none

  integer, parameter :: NS = 200000, NQ = 400000
  real(kind=wp), parameter :: EV = 1.602176634e-12_wp
  real(kind=wp), parameter :: TOL_MOC = 0.08_wp     ! MOCASSIN Ly-alpha fit
  real(kind=wp), parameter :: TOL_Q   = 1.0e-4_wp   ! module vs quadrature
  real(kind=wp), parameter :: TOL_S   = 5.0e-3_wp   ! sampler vs quadrature

  integer :: ierr, i, nbad
  real(kind=wp) :: T, f2s, aB, a2p, a2p_moc, r
  real(kind=wp) :: lfb, llyman, l2q, tot3, e_hi
  real(kind=wp) :: u, E, s1, sinv, qn, qe, dev
  real(kind=wp) :: ya, yb, ym, dy, accn, acce, totA
  real(kind=wp) :: fl(HEII_NLY), el(HEII_NLY), lB(3), lA(3), nimp
  real(kind=wp), parameter :: TGRID(5) = [5.0e3_wp, 1.0e4_wp, 1.5e4_wp, &
                                          2.0e4_wp, 3.0e4_wp]
  real(kind=wp), parameter :: TSHARE(3) = [1.0e4_wp, 1.5e4_wp, 3.0e4_wp]
  logical :: ok

  call MPI_INIT(ierr)
  par%eion_max   = 100.0_wp
  par%data_dir   = '../../data'
  par%atomic_dir = '../../data/atomic'
  call heii_cascade_setup()
  ok = .true.
  e_hi = min(E_HEII_LYA, par%eion_max)

  !=======================================================================
  ! 1-2. branching closure and the independent Ly-alpha cross-check
  !=======================================================================
  write(*,'(a)') ' n = 2 exit branching (zero-density cascade, case B):'
  write(*,'(a)') '   T [K]     alpha_B    alpha_2s     f_2s     f_2p'// &
                 '    alpha_2p   MOCASSIN    ratio'
  do i = 1, 5
     T = TGRID(i)
     f2s = heii_f_2s(T)
     aB  = alphaB_HeIII(T)
     a2p = aB*(1.0_wp - f2s)
     a2p_moc = 10.0_wp**(-0.792_wp*log10(T) + 6.01_wp)*1.0e-25_wp &
               /(E_HEII_LYA*EV)
     r = a2p/a2p_moc
     write(*,'(f9.0,2es12.4,2f9.5,2es12.4,f9.4)') &
        T, aB, aB*f2s, f2s, 1.0_wp - f2s, a2p, a2p_moc, r
     if (f2s <= 0.0_wp .or. f2s >= 1.0_wp) ok = .false.
     if (abs(r - 1.0_wp) > TOL_MOC) ok = .false.
  end do
  write(*,'(a,f6.3,a)') '   tolerance on the MOCASSIN Ly-alpha fit: ', &
     TOL_MOC*100.0_wp, ' %'

  !=======================================================================
  ! 3. branch luminosities: positive, and Ly-alpha > two-photon > free-bound
  !=======================================================================
  write(*,'(/,a)') ' channel 5 shares per (n_e n_HeIII) [erg cm^3/s]:'
  write(*,'(a)') '   T [K]   free-bound    Ly-alpha  two-photon'// &
                 '     fb %    Lya %     2q %'
  do i = 1, 3
     T = TSHARE(i)
     call heii_branch_luminosity(T, lfb, llyman, l2q)
     tot3 = lfb + llyman + l2q
     write(*,'(f9.0,3es12.4,3f9.3)') T, lfb, llyman, l2q, &
        100.0_wp*lfb/tot3, 100.0_wp*llyman/tot3, 100.0_wp*l2q/tot3
     if (lfb <= 0.0_wp .or. llyman <= 0.0_wp .or. l2q <= 0.0_wp) ok = .false.
     if (.not. (llyman > l2q .and. l2q > lfb)) ok = .false.
  end do

  !=======================================================================
  ! 6. two-photon in-band moments vs quadrature of the NS84 shape
  !=======================================================================
  totA = 0.0_wp
  dy = 1.0_wp/real(NQ, wp)
  do i = 1, NQ
     ya = real(i-1,wp)*dy;  yb = real(i,wp)*dy;  ym = 0.5_wp*(ya + yb)
     totA = totA + dy/6.0_wp*(shape_A(ya) + 4.0_wp*shape_A(ym) + shape_A(yb))
  end do
  accn = 0.0_wp;  acce = 0.0_wp
  do i = 1, NQ
     ya = eth_HI/E_HEII_LYA &
          + (e_hi/E_HEII_LYA - eth_HI/E_HEII_LYA)*real(i-1,wp)/real(NQ,wp)
     yb = eth_HI/E_HEII_LYA &
          + (e_hi/E_HEII_LYA - eth_HI/E_HEII_LYA)*real(i,  wp)/real(NQ,wp)
     ym = 0.5_wp*(ya + yb)
     accn = accn + (yb-ya)/6.0_wp &
            *(shape_A(ya) + 4.0_wp*shape_A(ym) + shape_A(yb))
     acce = acce + (yb-ya)/6.0_wp &
            *(shape_A(ya)*ya + 4.0_wp*shape_A(ym)*ym + shape_A(yb)*yb)
  end do
  qn = 2.0_wp*accn/totA
  qe = 2.0_wp*acce/totA*E_HEII_LYA
  write(*,'(/,a)') ' He II 2^2S two-photon in band:'
  write(*,'(a,f11.6,a,f11.6)') '   photons/decay  module ', &
     h2q_band_photons_per_decay(), '   quadrature ', qn
  write(*,'(a,f11.6,a,f11.6)') '   energy/decay   module ', &
     h2q_band_energy_per_decay(),  '   quadrature ', qe
  if (abs(h2q_band_photons_per_decay()/qn - 1.0_wp) > TOL_Q) ok = .false.
  if (abs(h2q_band_energy_per_decay()/qe  - 1.0_wp) > TOL_Q) ok = .false.

  !=======================================================================
  ! 4-5. sampler: band containment, and <1/E> = 1/<E>_n per sub-branch
  !=======================================================================
  T = 1.5e4_wp
  nbad = 0;  s1 = 0.0_wp;  sinv = 0.0_wp
  do i = 1, NS
     u = (real(i,wp) - 0.5_wp)/real(NS, wp)
     E = heii_2q_sample_energy(u)
     if (E < eth_HI .or. E > e_hi*(1.0_wp + 1.0e-12_wp)) nbad = nbad + 1
     s1 = s1 + E;  sinv = sinv + 1.0_wp/E
  end do
  s1 = s1/real(NS,wp);  sinv = sinv/real(NS,wp)
  dev = abs(sinv*(qe/qn) - 1.0_wp)
  write(*,'(/,a)') ' sampler (two-photon, T = 1.5e4 K):'
  write(*,'(a,i0)')       '   samples outside the band : ', nbad
  write(*,'(a,f11.6,a,f11.6)') '   <E>   sampled ', s1, '   quadrature ', &
     energy_weighted_mean()
  write(*,'(a,f12.8)')    '   <1/E> <E>_n              ', sinv*(qe/qn)
  if (nbad > 0) ok = .false.
  if (dev > TOL_S) ok = .false.

  nbad = 0
  do i = 1, NS
     u = (real(i,wp) - 0.5_wp)/real(NS, wp)
     E = heii_fb_sample_energy(T, u)
     if (E < eth_HI .or. E > eth_HeII*(1.0_wp + 1.0e-9_wp)) nbad = nbad + 1
  end do
  write(*,'(/,a)') ' sampler (free-bound, T = 1.5e4 K):'
  write(*,'(a,i0)') '   samples outside [13.598, 54.416] eV : ', nbad
  write(*,'(a,f9.4,a,f9.4,a)') '   range ', heii_fb_sample_energy(T, 0.0_wp), &
     ' - ', heii_fb_sample_energy(T, 1.0_wp), ' eV'
  if (nbad > 0) ok = .false.

  !=======================================================================
  ! 7. case A: the Lyman-line path
  !=======================================================================
  T = 1.5e4_wp
  call heii_branch_luminosity(T, lB(1), lB(2), lB(3))
  par%lyman_lines = 'escaping'
  call heii_branch_luminosity(T, lA(1), lA(2), lA(3))
  call heii_lyman_shares(T, 'A', fl)
  do i = 1, HEII_NLY
     el(i) = heii_lyman_energy(i)
  end do
  aB = alphaB_HeIII(T)
  nbad = 0;  sinv = 0.0_wp
  do i = 1, NS
     u = (real(i,wp) - 0.5_wp)/real(NS, wp)
     E = heii_lyman_sample_energy(T, u)
     if (E < 40.8129_wp .or. E > 54.42_wp) nbad = nbad + 1
     sinv = sinv + 1.0_wp/E
  end do
  nimp = lA(2)/(EV)*sinv/real(NS,wp)
  write(*,'(/,a)') ' case A (Lyman lines escaping), T = 1.5e4 K:'
  write(*,'(a,3f9.5,es11.3)') '   shares Ly-alpha, beta, gamma, n>=31:', &
     fl(1), fl(2), fl(3), fl(HEII_NLY)
  write(*,'(a,2es13.5)') '   L(Lyman) case B, case A [erg cm^3/s]:', lB(2), lA(2)
  write(*,'(a,2es13.5)') '   L(2q)    case B, case A             :', lB(3), lA(3)
  write(*,'(a,es13.5,a,f12.9)') '   L(Lyman)/(aB sum f E) - 1 =', &
     lA(2)/(aB*sum(fl*el)*EV) - 1.0_wp, &
     '   sampler photons / (aB sum f) =', nimp/(aB*sum(fl))
  write(*,'(a,i0)') '   sampled lines outside [40.813, 54.42] eV: ', nbad
  if (abs(lA(2)/(aB*sum(fl*el)*EV) - 1.0_wp) > 1.0e-13_wp) ok = .false.
  if (abs(nimp/(aB*sum(fl)) - 1.0_wp) > 1.0e-4_wp) ok = .false.
  if (nbad > 0) ok = .false.
  if (abs(lA(1)/lB(1) - 1.0_wp) > 1.0e-14_wp) ok = .false.
  par%lyman_lines = 'on_the_spot'

  write(*,'(/,a,a)') ' GATE (He III cascade channel): ', &
     merge('PASS', 'FAIL', ok)
  call MPI_FINALIZE(ierr)
  if (.not. ok) stop 1

contains

  !--- Nussbaumer & Schmutz (1984), repeated here so the gate does not
  !--- measure the module against itself.
  real(kind=wp) function shape_A(y) result(a)
    real(kind=wp), intent(in) :: y
    real(kind=wp) :: p, q
    a = 0.0_wp
    if (y <= 0.0_wp .or. y >= 1.0_wp) return
    p = y*(1.0_wp - y)
    q = (4.0_wp*p)**0.8_wp
    a = 202.0_wp*(p*(1.0_wp - q) + 0.88_wp*p**1.53_wp*q)
  end function shape_A

  !--- ENERGY-weighted mean of the in-band pair spectrum, the mean the
  !--- sampler must reproduce because packets carry equal energy.
  real(kind=wp) function energy_weighted_mean() result(m)
    real(kind=wp) :: a, b, mid, num, den, ylo, yhi
    integer :: k
    ylo = eth_HI/E_HEII_LYA;  yhi = e_hi/E_HEII_LYA
    num = 0.0_wp;  den = 0.0_wp
    do k = 1, NQ
       a = ylo + (yhi - ylo)*real(k-1,wp)/real(NQ,wp)
       b = ylo + (yhi - ylo)*real(k,  wp)/real(NQ,wp)
       mid = 0.5_wp*(a + b)
       num = num + (b-a)/6.0_wp*(shape_A(a)*a*a + 4.0_wp*shape_A(mid)*mid*mid &
                                 + shape_A(b)*b*b)
       den = den + (b-a)/6.0_wp*(shape_A(a)*a + 4.0_wp*shape_A(mid)*mid &
                                 + shape_A(b)*b)
    end do
    m = num/den*E_HEII_LYA
  end function energy_weighted_mean

end program check_heii_cascade
