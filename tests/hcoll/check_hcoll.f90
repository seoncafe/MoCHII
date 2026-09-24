program check_hcoll
!---------------------------------------------------------------------------
! Gate for the atomic-hydrogen impact rates of the C I, C II and O I ground
! terms (data/atomic/hcoll_<el>_<stage>.txt read by nlevel_mod):
!  (1) low-density limit: with n_e -> 0 and n_HI = 1e-4 cm^-3, every excitation
!      out of the ground level ends in a photon of the ground-term cascade, so
!      the photon rate of each fine-structure line reproduces the published
!      excitation rates -- Abrahamsson et al. (2007) Tables 1 and 3 at 100,
!      300, 1000 K; Barinovs et al. (2005) law for C II;
!  (2) detailed balance: at n_HI = 1e14 the ground-term populations are
!      Boltzmann at T, including above the 1000 K upper edge of the fits,
!      where the de-excitation rate is held and excitation follows T;
!  (3) the thermal increment nlevel_hcoll_cooling equals the emissivity
!      difference with and without n_HI.
! Also prints the old two-level [O I] 63 um / [C II] 158 um H-impact cooling
! that metal_cooling_H evaluated before 2026-09-24, for the record.
!---------------------------------------------------------------------------
  use define
  use nlevel_mod
  implicit none
  type(nlevel_atom_type) :: ao, ac, ac2
  logical :: ok1, ok2, ok3, fail
  real(kind=wp) :: emis(600), emis0(600), T, worst, dev
  real(kind=wp), parameter :: hc = 1.9864458571489287e-16_wp
  real(kind=wp), parameter :: Ts(3) = [100.0_wp, 300.0_wp, 1000.0_wp]
  real(kind=wp), parameter :: Tb(3) = [50.0_wp, 700.0_wp, 5000.0_wp]
  !--- O I Table 3: k(2->1), k(2->0), k(1->0) [1e-11]; C I Table 1 accurate
  real(kind=wp), parameter :: tabO(3,3) = reshape([2.20_wp, 0.243_wp, 5.38_wp, &
       15.83_wp, 3.24_wp, 19.65_wp, 43.75_wp, 11.01_wp, 32.03_wp], [3,3])
  real(kind=wp), parameter :: tabC(3,3) = reshape([29.75_wp, 23.83_wp, 29.79_wp, &
       42.90_wp, 51.70_wp, 54.79_wp, 64.15_wp, 88.64_wp, 90.61_wp], [3,3])
  !--- low enough that H collisions do not redistribute the excited levels
  !--- (n_HI k << A for the 370/609 um lines of C I)
  real(kind=wp), parameter :: NHI_LOW = 1.0e-4_wp
  character(len=*), parameter :: dir = '../../data/atomic/'
  integer :: it

  fail = .false.
  call nlevel_load(dir//'nlevel_o_1.txt', ao, ok1)
  call nlevel_load(dir//'nlevel_c_1.txt', ac, ok2)
  call nlevel_load(dir//'nlevel_c_2.txt', ac2, ok3)
  if (.not. (ok1 .and. ok2 .and. ok3 .and. nlevel_has_hcoll(ao) .and. &
      nlevel_has_hcoll(ac) .and. nlevel_has_hcoll(ac2))) then
     write(*,'(a)') 'FAIL: atoms or hydrogen rates not loaded'
     stop 1
  end if

  !--- (1) photon rates of the ground-term lines at n_e -> 0, divided by n_HI.
  !--- O I (1=3P2 2=3P1 3=3P0): photons(2->1) = k12 + k13 (3P0 cascades
  !--- through 3P1: A(3->1) = 1e-10 << A(3->2)), photons(3->2) = k13.
  !--- C I (1=3P0 2=3P1 3=3P2): photons(2->1) = k12 + k13 b32 with the
  !--- 3P2 branching to 3P1 dominant (A(3->1) is a 2e-14 quadrupole).
  worst = 0.0_wp
  do it = 1, 3
     T = Ts(it)
     call nlevel_emissivities(ao, T, 1.0e-12_wp, emis, NHI_LOW)
     dev = max(abs(line_rate(ao, emis, 2, 1)/NHI_LOW/(1e-11_wp*(tabO(1,it)+tabO(2,it))) - 1.0_wp), &
               abs(line_rate(ao, emis, 3, 2)/NHI_LOW/(1e-11_wp*tabO(2,it)) - 1.0_wp))
     worst = max(worst, dev)
     write(*,'(a,f7.0,a,2es12.4,a,es9.2)') ' O I  T=', T, '  photons 63/145 um ', &
        line_rate(ao, emis, 2, 1)/NHI_LOW, line_rate(ao, emis, 3, 2)/NHI_LOW, '  dev ', dev
     call nlevel_emissivities(ac, T, 1.0e-12_wp, emis, NHI_LOW)
     dev = abs(line_rate(ac, emis, 2, 1)/NHI_LOW/(1e-11_wp*(tabC(1,it)+tabC(2,it))) - 1.0_wp)
     worst = max(worst, dev)
     write(*,'(a,f7.0,a,es12.4,a,es9.2)') ' C I  T=', T, '  photons 609 um    ', &
        line_rate(ac, emis, 2, 1)/NHI_LOW, '  dev ', dev
     call nlevel_emissivities(ac2, T, 1.0e-12_wp, emis, NHI_LOW)
     dev = abs(line_rate(ac2, emis, 2, 1)/NHI_LOW / (2.0_wp*7.6e-10_wp*(T/100.0_wp)**0.14_wp &
               *exp(-(ac2%e_cm1(2))*hc/(kboltz_cgs*T))) - 1.0_wp)
     worst = max(worst, dev)
     write(*,'(a,f7.0,a,es12.4,a,es9.2)') ' C II T=', T, '  photons 158 um    ', &
        line_rate(ac2, emis, 2, 1)/NHI_LOW, '  dev ', dev
  end do
  write(*,'(a,es9.2)') ' (1) low-density line rates vs published: worst ', worst
  !--- tolerance: the eq. (1) fit reproduces its tables to 1.0%; the n_e =
  !--- 1e-12 electron term and the weak radiative branches add < 1e-4.
  if (worst > 1.2e-2_wp) fail = .true.

  !--- (2) Boltzmann at n_HI = 1e14
  worst = 0.0_wp
  do it = 1, 3
     T = Tb(it)
     worst = max(worst, boltzmann_dev(ao, T), boltzmann_dev(ac, T), boltzmann_dev(ac2, T))
  end do
  write(*,'(a,es9.2)') ' (2) Boltzmann ground term at n_HI = 1e14: worst ', worst
  if (worst > 1.0e-6_wp) fail = .true.

  !--- (3) increment = emissivity difference
  worst = 0.0_wp
  do it = 1, 3
     T = Ts(it)
     call nlevel_emissivities(ao, T, 10.0_wp, emis, 1.0e4_wp)
     call nlevel_emissivities(ao, T, 10.0_wp, emis0)
     dev = abs(nlevel_hcoll_cooling(ao, T, 10.0_wp, 1.0e4_wp) &
               / (sum(emis(1:ao%nrad)) - sum(emis0(1:ao%nrad))) - 1.0_wp)
     worst = max(worst, dev)
  end do
  write(*,'(a,es9.2)') ' (3) increment vs emissivity difference: worst ', worst
  if (worst > 1.0e-10_wp) fail = .true.

  !--- record: old two-level H-impact cooling [erg/s per ion] vs the new
  !--- increment at n_e = 1e-6 cm^-3 (hydrogen alone) for three n_HI.
  write(*,'(a)') ' record: H-impact cooling per ion, old two-level vs joint solve (n_e = 1e-6)'
  do it = 1, 3
     T = Ts(it)
     block
       real(kind=wp) :: nHI, qul, qlu, old
       integer :: jn
       do jn = 1, 3
          nHI = 10.0_wp**(2*jn)
          qul = 9.2e-11_wp*(T/100.0_wp)**0.67_wp;  qlu = 0.6_wp*qul*exp(-227.7_wp/T)
          old = nHI*qlu*227.7_wp*kboltz_cgs/(1.0_wp + nHI*(qul+qlu)/8.865e-5_wp)
          write(*,'(a,f6.0,a,es8.1,a,2es11.3,f8.3)') '   O I  T=', T, ' nHI=', nHI, &
             '  old/new ', old, nlevel_hcoll_cooling(ao, T, 1.0e-6_wp, nHI), &
             old/nlevel_hcoll_cooling(ao, T, 1.0e-6_wp, nHI)
          qul = 7.6e-10_wp*(T/100.0_wp)**0.14_wp;  qlu = 2.0_wp*qul*exp(-91.25_wp/T)
          old = nHI*qlu*91.25_wp*kboltz_cgs/(1.0_wp + nHI*(qul+qlu)/2.321e-6_wp)
          write(*,'(a,f6.0,a,es8.1,a,2es11.3,f8.3)') '   C II T=', T, ' nHI=', nHI, &
             '  old/new ', old, nlevel_hcoll_cooling(ac2, T, 1.0e-6_wp, nHI), &
             old/nlevel_hcoll_cooling(ac2, T, 1.0e-6_wp, nHI)
          write(*,'(a,f6.0,a,es8.1,a,11x,es11.3)') '   C I  T=', T, ' nHI=', nHI, &
             '  old/new ', nlevel_hcoll_cooling(ac, T, 1.0e-6_wp, nHI)
       end do
     end block
  end do

  if (fail) then
     write(*,'(a)') 'FAIL'
     stop 1
  end if
  write(*,'(a)') 'PASS'

contains

  !--- photon rate [s^-1 per ion] of the radiative row joining levels u -> l
  real(kind=wp) function line_rate(a, em, u, l) result(r)
    type(nlevel_atom_type), intent(in) :: a
    real(kind=wp), intent(in) :: em(:)
    integer, intent(in) :: u, l
    integer :: k
    r = 0.0_wp
    do k = 1, a%nrad
       if ((a%ru(k) == u .and. a%rl(k) == l) .or. (a%ru(k) == l .and. a%rl(k) == u)) &
          r = r + em(k)/((a%e_cm1(u) - a%e_cm1(l))*hc)
    end do
  end function line_rate

  !--- worst |n_i/n_1 / Boltzmann - 1| over the hydrogen-coupled levels,
  !--- read from the emissivities (n_u = emis/(A dE)).
  real(kind=wp) function boltzmann_dev(a, T) result(d)
    type(nlevel_atom_type), intent(in) :: a
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: em(600), pop(40), b
    integer :: k, u, l
    call nlevel_emissivities(a, T, 1.0e-12_wp, em, 1.0e14_wp)
    pop = -1.0_wp
    do k = 1, a%nrad
       u = a%ru(k);  l = a%rl(k)
       if (a%e_cm1(u) < a%e_cm1(l)) then
          u = a%rl(k);  l = a%ru(k)
       end if
       if (a%A(k) > 0.0_wp .and. a%e_cm1(u) > a%e_cm1(l)) &
          pop(u) = em(k)/(a%A(k)*(a%e_cm1(u) - a%e_cm1(l))*hc)
    end do
    !--- two hydrogen-coupled levels (C II): n_2 = b/(1+b) of the ion, the
    !--- excited terms being empty; three (C I, O I): n_3/n_2 = b_32.
    d = 0.0_wp
    if (a%nhc == 1) then
       b = (a%g(2)/a%g(1))*exp(-(a%e_cm1(2) - a%e_cm1(1))*hc/(kboltz_cgs*T))
       d = abs(pop(2)/(b/(1.0_wp + b)) - 1.0_wp)
    else
       b = (a%g(3)/a%g(2))*exp(-(a%e_cm1(3) - a%e_cm1(2))*hc/(kboltz_cgs*T))
       d = abs((pop(3)/pop(2))/b - 1.0_wp)
    end if
  end function boltzmann_dev

end program check_hcoll
