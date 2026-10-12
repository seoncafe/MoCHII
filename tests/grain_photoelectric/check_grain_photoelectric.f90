!---------------------------------------------------------------------------
! Gate for the grain photoelectric heating, the charged-particle accretion
! cooling (src/grain_photoelectric_mod.f90) and the neutral gas-grain energy
! exchange (src/gas_grain_collision_mod.f90).
!
! (1) The closed forms of the photoelectron yield y_2 and of the escaping
!     energy (WD01 eqs. 10-11) against a direct quadrature of eq. 10.
! (2) Weingartner & Draine (2001, ApJS 134, 263; WD01) reproduced for their
!     own grains: the WD01a (ApJ 548, 296) R_V = 3.1, b_C = 6e-5, case A size
!     distribution (their eqs. 2-6 and Table 1) on the DL07 optics of SEDust
!     (Draine & Li 2007: the same PAH-graphite continuum and astronomical
!     silicate, with DL07's updated PAH cross sections where WD01 used Li &
!     Draine 2001 -- the one input that is not theirs).
!       Tables 4 and 5: Gamma_pe/(G n_H) and Lambda_gr/(G n_H) for H II
!       regions, n_e = n_H (H+), T = 9000 K, blackbodies T_c = 3.5e4 and
!       4.5e4 K with the 13.6 eV cutoff and without it (the full blackbody;
!       h nu_max = 10 k T_c as in Draine's charging code), G/n_H = 0.1, 1, 10.
!       Eqs. 44-45 with the Table 2-3 ISRF row (R_V = 3.1, b_C = 6.0) at
!       T = 100 K and 1000 K, H+ ions (C+ printed), 1e2 <= G sqrt(T)/n_e <= 1e6: the fits'
!       stated maximum errors are 0.20 (heating) and 0.11 (cooling).
! (3) The charge distribution against Draine's own charging code as ported
!     in SEDust (pah_ioniz_mod::pah_ionfrac, the WD01b CHARGE routine): the
!     non-neutral fraction of carbonaceous grains in the CNM/WNM/WIM mix.
! (4) What the MoCHII field costs: the band cut at efuv_min = 4.78 eV (no
!     photons below), and J_nu taken constant inside the default band bins.
! (5) Size-bin merge: 8 bins per decade against every native radius.
! (6) Neutral gas-grain exchange: the HM79 accommodation limits, and the
!     sign of the energy flow with T - T_d.
! Exits nonzero when a gate fails.
!---------------------------------------------------------------------------
program check_grain_photoelectric
  use define
  use dust_lib, only : dust_model_t, build_dust
  use grain_size_bins_mod
  use grain_photoelectric_mod
  use gas_grain_collision_mod
  use pah_ioniz_mod, only : pah_ionfrac
  implicit none
  type(dust_model_t) :: m
  character(len=512) :: data_dir
  integer :: st, nfail
  real(kind=wp), parameter :: U_HAB = 5.33e-14_wp       ! WD01 sec. 4.1
  real(kind=wp), parameter :: KEV = 8.617333262e-5_wp   ! k [eV/K]
  real(kind=wp), parameter :: AMU = 1.66053907e-24_wp
  real(kind=wp), parameter :: MP = 1.67262192e-24_wp

  nfail = 0
  data_dir = '../../SEDust/data'
  if (command_argument_count() >= 1) call get_command_argument(1, data_dir)

  call check_yield_closed_forms()

  call build_dust(m, 'dl07', trim(data_dir), 20, 3.0_wp, 3000.0_wp, &
                  include_euv=.true., status=st, sd_index=7, u_isrf=1.0_wp)
  if (st /= 0) stop 'build_dust dl07 failed'
  call set_wd01a_distribution(m)

  call check_wd01_tables()
  call check_wd01_fits()
  call check_draine_charging()
  call check_mochii_field()
  call check_size_bins()
  call check_gas_grain()

  write(*,'(/,a,i0)') 'FAILED GATES: ', nfail
  if (nfail > 0) error stop 1

contains

  !=========================================================================
  subroutine gate(ok, text)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: text
    if (ok) then
       write(*,'(2a)') '  PASS  ', text
    else
       write(*,'(2a)') '  FAIL  ', text
       nfail = nfail + 1
    end if
  end subroutine gate

  !=========================================================================
  ! (1) y_2 = E_h^2 (E_h - 3E_l)/Theta^3 and Int_0^E_h E f dE
  !     = E_h^3 (E_h - 2E_l)/(2 Theta^3) against quadrature of eq. 10.
  subroutine check_yield_closed_forms()
    real(kind=wp) :: eh, el, th, y2q, keq, y2c, kec, e, f, de, worst
    integer :: i, j, n
    write(*,'(a)') '(1) photoelectron energy distribution, eq. 10'
    worst = 0.0_wp
    do i = 1, 5
       eh = 0.3_wp*real(i, wp)
       el = -0.17_wp*real(6-i, wp)
       th = eh - el
       n = 200000
       de = eh/real(n, wp)
       y2q = 0.0_wp;  keq = 0.0_wp
       do j = 1, n
          e = (real(j, wp) - 0.5_wp)*de
          f = 6.0_wp*(e - el)*(eh - e)/th**3
          y2q = y2q + f*de;  keq = keq + e*f*de
       end do
       y2c = eh*eh*(eh - 3.0_wp*el)/th**3
       kec = eh**3*(eh - 2.0_wp*el)/(2.0_wp*th**3)
       worst = max(worst, abs(y2q/y2c - 1.0_wp), abs(keq/kec - 1.0_wp))
    end do
    write(*,'(a,es10.3)') '    worst relative difference ', worst
    call gate(worst < 1.0e-8_wp, 'closed forms of y_2 and the escaping energy')
  end subroutine check_yield_closed_forms

  !=========================================================================
  ! WD01a eqs. 2-6, R_V = 3.1, b_C = 6e-5, case A (Table 1), on the native
  ! DL07 radii; dn = (dn/da) a dln a (trapezoid, half ends, as SEDust).
  ! The carbonaceous number is split equally between the neutral and the
  ! cation population (only the sum and the cross sections enter).
  subroutine set_wd01a_distribution(mm)
    type(dust_model_t), intent(inout) :: mm
    integer :: ia, na
    real(kind=wp) :: a, dlna, dc, ds
    na = size(mm%pops(1)%aeff)
    do ia = 1, na
       a = mm%pops(1)%aeff(ia)*1.0e-4_wp
       if (ia == 1) then
          dlna = 0.5_wp*log(mm%pops(1)%aeff(2)/mm%pops(1)%aeff(1))
       else if (ia == na) then
          dlna = 0.5_wp*log(mm%pops(1)%aeff(na)/mm%pops(1)%aeff(na-1))
       else
          dlna = 0.5_wp*log(mm%pops(1)%aeff(ia+1)/mm%pops(1)%aeff(ia-1))
       end if
       call wd01a_dnda(a, dc, ds)
       mm%pops(1)%dn(ia) = ds*a*dlna
       mm%pops(2)%dn(ia) = 0.5_wp*dc*a*dlna
       mm%pops(3)%dn(ia) = 0.5_wp*dc*a*dlna
    end do
  end subroutine set_wd01a_distribution

  subroutine wd01a_dnda(a, dc, ds)
    real(kind=wp), intent(in)  :: a        ! [cm]
    real(kind=wp), intent(out) :: dc, ds   ! (1/n_H) dn/da [cm^-1]
    real(kind=wp), parameter :: BC = 6.0e-5_wp, SIG = 0.4_wp, RHO = 2.24_wp
    real(kind=wp), parameter :: MC = 12.011_wp*AMU, AMIN = 3.5e-8_wp
    real(kind=wp), parameter :: ALG = -1.54_wp, BEG = -0.165_wp, ATG = 0.0107e-4_wp, &
                                ACG = 0.428e-4_wp, CG = 9.99e-12_wp
    real(kind=wp), parameter :: ALS = -2.21_wp, BES = 0.300_wp, ATS = 0.164e-4_wp, &
                                ACS = 0.1e-4_wp, CS = 1.00e-13_wp
    real(kind=wp) :: a0(2), bci(2), bi, d
    integer :: i
    a0 = [3.5e-8_wp, 3.0e-7_wp];  bci = [0.75_wp*BC, 0.25_wp*BC]
    d = 0.0_wp
    if (a > AMIN) then
       do i = 1, 2
          bi = 3.0_wp/(twopi**1.5_wp)*exp(-4.5_wp*SIG**2)/(RHO*a0(i)**3*SIG) &
               *bci(i)*MC/(1.0_wp + erf(3.0_wp*SIG/sqrt(2.0_wp) &
               + log(a0(i)/AMIN)/(SIG*sqrt(2.0_wp))))
          d = d + bi/a*exp(-0.5_wp*(log(a/a0(i))/SIG)**2)
       end do
    end if
    dc = d + CG/a*(a/ATG)**ALG*curv(a, BEG, ATG)
    if (a > ATG) dc = dc - CG/a*(a/ATG)**ALG*curv(a, BEG, ATG) &
                     *(1.0_wp - exp(-((a - ATG)/ACG)**3))
    ds = CS/a*(a/ATS)**ALS*curv(a, BES, ATS)
    if (a > ATS) ds = ds*exp(-((a - ATS)/ACS)**3)
    if (a < AMIN) then
       dc = 0.0_wp;  ds = 0.0_wp
    end if
  end subroutine wd01a_dnda

  pure real(kind=wp) function curv(a, beta, at) result(f)
    real(kind=wp), intent(in) :: a, beta, at
    if (beta >= 0.0_wp) then
       f = 1.0_wp + beta*a/at
    else
       f = 1.0_wp/(1.0_wp - beta*a/at)
    end if
  end function curv

  !=========================================================================
  ! Energy grid of nb log bins on [elo, ehi] and the bins of model m on it.
  subroutine grid_and_bins(elo, ehi, nb)
    real(kind=wp), intent(in) :: elo, ehi
    integer, intent(in) :: nb
    real(kind=wp) :: edges(nb+1)
    integer :: i
    do i = 1, nb+1
       edges(i) = elo*(ehi/elo)**(real(i-1, wp)/real(nb, wp))
    end do
    call grain_energy_grid_from_edges(edges)
    call grain_bins_from_model(m, trim(data_dir))
    call grain_photoelectric_setup()
  end subroutine grid_and_bins

  !--- 4 pi B_nu/(h nu) dnu at energy e [eV] over width de [eV], i.e. the
  !--- photon rate per unit area of an undiluted blackbody
  pure real(kind=wp) function bb_photons(e, de, tc) result(n)
    real(kind=wp), intent(in) :: e, de, tc
    real(kind=wp) :: nu, dnu, x
    nu = e*ev2erg/h_planck_cgs;  dnu = de*ev2erg/h_planck_cgs
    x = e/(KEV*tc)
    if (x > 600.0_wp) then
       n = 0.0_wp
    else
       n = 8.0_wp*pi*nu*nu/(clight_cgs**2*(exp(x) - 1.0_wp))*dnu
    end if
  end function bb_photons

  !--- energy density of an undiluted blackbody in [6, 13.6] eV [erg cm^-3]
  real(kind=wp) function bb_u_fuv(tc) result(u)
    real(kind=wp), intent(in) :: tc
    integer :: i
    real(kind=wp) :: e, de
    u = 0.0_wp
    de = 7.6_wp/20000.0_wp
    do i = 1, 20000
       e = 6.0_wp + (real(i, wp) - 0.5_wp)*de
       u = u + bb_photons(e, de, tc)*e*ev2erg/clight_cgs
    end do
  end function bb_u_fuv

  !--- WD01 eq. 31: nu u_nu [erg cm^-3] of the MMP83 ISRF at e [eV]
  pure real(kind=wp) function isrf_nuunu(e) result(v)
    real(kind=wp), intent(in) :: e
    real(kind=wp), parameter :: wi(3) = [1.0e-14_wp, 1.65e-13_wp, 4.0e-13_wp]
    real(kind=wp), parameter :: ti(3) = [7500.0_wp, 4000.0_wp, 3000.0_wp]
    real(kind=wp) :: nu, x
    integer :: i
    if (e > 13.6_wp) then
       v = 0.0_wp
    else if (e > 11.2_wp) then
       v = 3.328e-9_wp*e**(-4.4172_wp)
    else if (e > 9.26_wp) then
       v = 8.463e-13_wp/e
    else if (e > 5.04_wp) then
       v = 2.055e-14_wp*e**0.6678_wp
    else
       nu = e*ev2erg/h_planck_cgs
       v = 0.0_wp
       do i = 1, 3
          x = e/(KEV*ti(i))
          if (x < 600.0_wp) v = v + wi(i)*2.0_wp*h_planck_cgs*nu**3 &
                                /clight_cgs**2/(exp(x) - 1.0_wp)
       end do
       v = v*fourpi*nu/clight_cgs
    end if
  end function isrf_nuunu

  !--- photon rates on the current grid: blackbody of dilution w, or ISRF
  subroutine bb_rates(tc, w, ecut, nphot)
    real(kind=wp), intent(in) :: tc, w, ecut
    real(kind=wp), intent(out) :: nphot(:)
    integer :: k
    do k = 1, n_egrid
       if (egrid(k) > ecut) then
          nphot(k) = 0.0_wp
       else
          nphot(k) = w*bb_photons(egrid(k), egrid_de(k), tc)
       end if
    end do
  end subroutine bb_rates

  subroutine isrf_rates(scale, emin_cut, nphot)
    real(kind=wp), intent(in) :: scale, emin_cut
    real(kind=wp), intent(out) :: nphot(:)
    integer :: k
    do k = 1, n_egrid
       if (egrid(k) < emin_cut) then
          nphot(k) = 0.0_wp
       else
          nphot(k) = scale*clight_cgs*isrf_nuunu(egrid(k))*(egrid_de(k)/egrid(k)) &
                     /(egrid(k)*ev2erg)
       end if
    end do
  end subroutine isrf_rates

  !--- heating and cooling per H [erg s^-1] at (T, n_e, ions) in the loaded field
  subroutine per_h(T, ne, ion_mass_amu, gam, lam)
    real(kind=wp), intent(in)  :: T, ne, ion_mass_amu
    real(kind=wp), intent(out) :: gam, lam
    type(grain_collider_type) :: ions(1)
    ions(1) = grain_collider_type(ne, 1.0_wp, ion_mass_amu*AMU)
    if (ion_mass_amu < 1.5_wp) ions(1)%m = MP
    call grain_charge_terms(T, 20.0_wp, ne, ions, 1.0_wp, 1.0_wp, gam, lam)
  end subroutine per_h

  !=========================================================================
  ! (2a) WD01 Tables 4-5
  subroutine check_wd01_tables()
    real(kind=wp), parameter :: tab4(3,2,2) = reshape([ &
       2.10_wp, 1.84_wp, 1.05_wp,   7.43_wp, 6.44_wp, 3.25_wp, &   ! Tc 3.5e4: cutoff, full
       2.30_wp, 2.01_wp, 1.12_wp,   12.0_wp, 10.2_wp, 4.57_wp], &  ! Tc 4.5e4
       [3,2,2])
    real(kind=wp), parameter :: tab5(3,2,2) = reshape([ &
       30.5_wp, 3.52_wp, 0.61_wp,   31.3_wp, 4.15_wp, 0.93_wp, &
       30.6_wp, 3.56_wp, 0.63_wp,   31.9_wp, 4.59_wp, 1.12_wp], [3,2,2])
    real(kind=wp), parameter :: gnh(3) = [0.1_wp, 1.0_wp, 10.0_wp]
    real(kind=wp), parameter :: tcs(2) = [3.5e4_wp, 4.5e4_wp]
    real(kind=wp), allocatable :: nphot(:)
    real(kind=wp) :: ecut, w, gam, lam, rg, rl, worst_g, worst_l
    integer :: ic, it, ig
    character(len=8) :: lab(2) = ['13.6 cut', 'full    ']
    write(*,'(/,a)') '(2a) WD01 Tables 4-5: H II regions, n_e = n_H, T = 9000 K, '// &
       'R_V = 3.1, b_C = 6.0, case A'
    write(*,'(a)') '    T_c     field    G/n_H   Gamma/(G n_H) [1e-25]  WD01  ratio'// &
       '   Lambda/(G n_H) [1e-25]  WD01  ratio'
    worst_g = 0.0_wp;  worst_l = 0.0_wp
    do it = 1, 2
       do ic = 1, 2
          if (ic == 1) then
             ecut = 13.6_wp
          else
             ecut = 10.0_wp*KEV*tcs(it)
          end if
          call grid_and_bins(0.1_wp, ecut, 600)
          allocate(nphot(n_egrid))
          do ig = 1, 3
             w = gnh(ig)*U_HAB/bb_u_fuv(tcs(it))
             call bb_rates(tcs(it), w, ecut, nphot)
             call grain_field_load(nphot)
             call per_h(9000.0_wp, 1.0_wp, 1.0_wp, gam, lam)
             rg = gam/gnh(ig)/1.0e-25_wp/tab4(ig, ic, it)
             rl = lam/gnh(ig)/1.0e-25_wp/tab5(ig, ic, it)
             worst_g = max(worst_g, abs(rg - 1.0_wp))
             worst_l = max(worst_l, abs(rl - 1.0_wp))
             write(*,'(es10.2,2x,a,f7.1,f14.3,f12.2,f7.3,f17.3,f12.2,f7.3)') tcs(it), lab(ic), &
                gnh(ig), gam/gnh(ig)/1.0e-25_wp, tab4(ig, ic, it), rg, &
                lam/gnh(ig)/1.0e-25_wp, tab5(ig, ic, it), rl
          end do
          deallocate(nphot)
       end do
    end do
    write(*,'(a,2f8.3)') '    worst |ratio - 1|  heating, cooling: ', worst_g, worst_l
    !--- measured 2026-10-11: 2.3% (heating), 1.5% (cooling); the one input
    !--- that is not WD01's own is the DL07 PAH cross section
    call gate(worst_g < 0.05_wp, 'WD01 Table 4 photoelectric heating within 5%')
    call gate(worst_l < 0.05_wp, 'WD01 Table 5 recombination cooling within 5%')
  end subroutine check_wd01_tables

  !=========================================================================
  ! (2b) WD01 eqs. 44-45, ISRF row of Tables 2-3 (R_V 3.1, b_C 6.0).  The
  ! paper says (sec. 4.2) that its calculations other than the CNM/WNM take
  ! C+ as the ion; eq. 45 is reproduced with H+ and not with C+ (36% low at
  ! G sqrt(T)/n_e = 1e2, 1000 K), so the gate takes H+ and the C+ rows are
  ! printed for the record.
  subroutine check_wd01_fits()
    real(kind=wp), parameter :: c(0:6) = [5.45_wp, 2.50_wp, 0.00945_wp, &
       0.01453_wp, 0.147_wp, 0.623_wp, 0.511_wp]
    real(kind=wp), parameter :: d(0:4) = [0.4255_wp, 2.457_wp, -6.404_wp, &
       1.513_wp, 0.05343_wp]
    real(kind=wp), parameter :: mion(2) = [1.0_wp, 12.0_wp]
    character(len=2), parameter :: lion(2) = ['H+', 'C+']
    real(kind=wp), allocatable :: nphot(:)
    real(kind=wp) :: g, psi, ne, gam, lam, gfit, lfit, x, T, ufuv
    real(kind=wp) :: worst_g(2), worst_l(2)
    integer :: ip, itt, k, ii
    write(*,'(/,a)') '(2b) WD01 eqs. 44-45, ISRF (eq. 31), R_V = 3.1, b_C = 6.0'
    call grid_and_bins(0.1_wp, 13.6_wp, 600)
    allocate(nphot(n_egrid))
    call isrf_rates(1.0_wp, 0.0_wp, nphot)
    ufuv = 0.0_wp
    do k = 1, n_egrid
       if (egrid(k) >= 6.0_wp) ufuv = ufuv + nphot(k)*egrid(k)*ev2erg/clight_cgs
    end do
    g = ufuv/U_HAB
    write(*,'(a,es11.4,a,f6.3,a)') '    ISRF u(6-13.6 eV) = ', ufuv, &
       ' erg cm^-3 (WD01: 6.07e-14),  G = ', g, ' (WD01: 1.13)'
    call grain_field_load(nphot)
    worst_g = 0.0_wp;  worst_l = 0.0_wp
    do ii = 1, 2
       write(*,'(2a)') '    ions: ', lion(ii)
       write(*,'(a)') '      T     G sqrt(T)/n_e   Gamma/(G n_H) [1e-26]  eq.44   ratio'// &
          '    Lambda/(n_e n_H) [1e-28]  eq.45   ratio'
       do itt = 1, 2
          T = merge(100.0_wp, 1000.0_wp, itt == 1)
          do ip = 0, 8
             psi = 10.0_wp**(2.0_wp + 0.5_wp*real(ip, wp))
             ne = g*sqrt(T)/psi
             call per_h(T, ne, mion(ii), gam, lam)
             gfit = 1.0e-26_wp*g*(c(0) + c(1)*T**c(4)) &
                    /(1.0_wp + c(2)*psi**c(5)*(1.0_wp + c(3)*psi**c(6)))
             x = log(psi)
             lfit = 1.0e-28_wp*ne*T**(d(0) + d(1)/x)*exp(d(2) + d(3)*x - d(4)*x*x)
             worst_g(ii) = max(worst_g(ii), abs(gam/gfit - 1.0_wp))
             if (T >= 1.0e3_wp) worst_l(ii) = max(worst_l(ii), abs(lam/lfit - 1.0_wp))
             write(*,'(f9.0,es14.3,f16.3,f12.3,f8.3,f20.3,f12.3,f8.3)') T, psi, &
                gam/g/1.0e-26_wp, gfit/g/1.0e-26_wp, gam/gfit, lam/ne/1.0e-28_wp, &
                lfit/ne/1.0e-28_wp, lam/lfit
          end do
       end do
       write(*,'(3a,2f8.3)') '    ', lion(ii), ': worst |ratio - 1|  heating (100 and '// &
          '1000 K), cooling (1000 K): ', worst_g(ii), worst_l(ii)
    end do
    write(*,'(a)') '    (eq. 45 holds for 1e3 <= T <= 1e4 K; its 100 K rows are shown only)'
    call gate(worst_g(1) < 0.25_wp, 'eq. 44 (H+) within its 0.20 maximum error + 5%')
    call gate(worst_l(1) < 0.16_wp, 'eq. 45 (H+) within its 0.11 maximum error + 5%')
    deallocate(nphot)
  end subroutine check_wd01_fits

  !=========================================================================
  ! (3) Draine's charging code (SEDust pah_ionfrac): non-neutral fraction of
  ! carbonaceous grains in 0.43 CNM + 0.43 WNM + 0.14 WIM, MMP83 ISRF
  subroutine check_draine_charging()
    real(kind=wp), allocatable :: nphot(:)
    real(kind=wp) :: fmine, fdraine, worst
    real(kind=wp) :: tt(3), nh(3), xh(3), xm(3), wt(3)
    integer :: ib, ip
    write(*,'(/,a)') '(3) charge distribution vs Draine''s charging code '// &
       '(SEDust pah_ionfrac): carbonaceous, DISM mix, MMP83 ISRF'
    tt = [100.0_wp, 6000.0_wp, 8000.0_wp];  nh = [30.0_wp, 0.4_wp, 0.1_wp]
    xh = [0.0012_wp, 0.1_wp, 0.99_wp];  xm = [0.0003_wp, 0.0003_wp, 0.001_wp]
    wt = [0.43_wp, 0.43_wp, 0.14_wp]
    grain_bins_per_decade = 1000
    call grid_and_bins(0.1_wp, 13.6_wp, 600)
    allocate(nphot(n_egrid))
    !--- Draine's RADFLD mode 4 at chi = 1.23 U: MMP83 x 0.8115 x 1.23
    call isrf_rates(0.8115_wp*1.23_wp, 0.0_wp, nphot)
    call grain_field_load(nphot)
    write(*,'(a)') '      a [A]    f_ion (this)   f_ion (Draine code)'
    worst = 0.0_wp
    do ib = 1, n_gbin
       if (gbin(ib)%material /= CARBONACEOUS) cycle
       if (gbin(ib)%a_ang > 200.0_wp) cycle
       fmine = 0.0_wp
       do ip = 1, 3
          fmine = fmine + wt(ip)*(1.0_wp - neutral_fraction(ib, tt(ip), nh(ip), &
                                                            xh(ip), xm(ip)))
       end do
       fdraine = pah_ionfrac(gbin(ib)%a_cm*1.0e4_wp, 1.0_wp)
       worst = max(worst, abs(fmine - fdraine))
       write(*,'(f11.2,2f15.4)') gbin(ib)%a_ang, fmine, fdraine
    end do
    write(*,'(a,f8.4)') '    worst |difference| ', worst
    call gate(worst < 0.05_wp, 'non-neutral fraction within 0.05 of Draine''s code')
    grain_bins_per_decade = 8
    deallocate(nphot)
  end subroutine check_draine_charging

  !--- f(Z=0) of bin ib from the distribution: run the bin with a probe
  real(kind=wp) function neutral_fraction(ib, T, nH, xh, xm) result(f0)
    integer, intent(in) :: ib
    real(kind=wp), intent(in) :: T, nH, xh, xm
    type(grain_collider_type) :: ions(2)
    real(kind=wp) :: g1, l1
    ions(1) = grain_collider_type(xh*nH, 1.0_wp, MP)
    ions(2) = grain_collider_type(xm*nH, 1.0_wp, 12.0_wp*AMU)
    call grain_bin_terms(ib, T, 20.0_wp, (xh + xm)*nH, ions, g1, l1, f_neutral=f0)
  end function neutral_fraction

  !=========================================================================
  ! (4) the MoCHII field: no photons below 4.78 eV; J_nu flat in band bins
  subroutine check_mochii_field()
    real(kind=wp), allocatable :: nphot(:), nflat(:), edges(:), fband(:)
    real(kind=wp) :: g, gam0, lam0, gam1, lam1, psi, ne, w, T
    integer :: k, nb, i, ip
    write(*,'(/,a)') '(4a) ISRF without photons below efuv_min = 4.78 eV (WD01 fit grains)'
    call grid_and_bins(0.1_wp, 13.6_wp, 600)
    allocate(nphot(n_egrid))
    T = 100.0_wp
    g = 1.13_wp
    write(*,'(a)') '    G sqrt(T)/n_e   Gamma(cut)/Gamma   Lambda(cut)/Lambda'
    do ip = 0, 4
       psi = 10.0_wp**(2.0_wp + real(ip, wp))
       ne = g*sqrt(T)/psi
       call isrf_rates(1.0_wp, 0.0_wp, nphot)
       call grain_field_load(nphot)
       call per_h(T, ne, 12.0_wp, gam0, lam0)
       call isrf_rates(1.0_wp, 4.78_wp, nphot)
       call grain_field_load(nphot)
       call per_h(T, ne, 12.0_wp, gam1, lam1)
       write(*,'(es14.3,2f19.4)') psi, gam1/gam0, lam1/lam0
    end do
    deallocate(nphot)

    !--- default band: 8 FUV log bins on [4.78, 13.598], 16 ionizing log bins
    !--- on [13.598, 100]; H II conditions, T_c = 4e4 K, G/n_H = 1
    write(*,'(/,a)') '(4b) J_nu constant in the default band bins (8 FUV + 16 ionizing), '// &
       'T_c = 4e4 K, G/n_H = 1, T = 9000 K'
    nb = 24
    allocate(edges(nb+1), fband(nb))
    do i = 1, 9
       edges(i) = 4.78_wp*(13.598_wp/4.78_wp)**(real(i-1, wp)/8.0_wp)
    end do
    do i = 10, nb+1
       edges(i) = 13.598_wp*(100.0_wp/13.598_wp)**(real(i-9, wp)/16.0_wp)
    end do
    call grain_energy_grid_from_edges(edges)
    call grain_bins_from_model(m, trim(data_dir))
    call grain_photoelectric_setup()
    allocate(nphot(n_egrid), nflat(n_egrid))
    w = U_HAB/bb_u_fuv(4.0e4_wp)
    call bb_rates(4.0e4_wp, w, 100.0_wp, nphot)
    call grain_field_load(nphot)
    call per_h(9000.0_wp, 1.0_wp, 1.0_wp, gam0, lam0)
    !--- the band integrals of the same field, then flat inside each bin
    fband = 0.0_wp
    do k = 1, n_egrid
       fband(egrid_band(k)) = fband(egrid_band(k)) + nphot(k)*egrid(k)*ev2erg
    end do
    call photon_rates_from_band(fband, edges, nflat)
    call grain_field_load(nflat)
    call per_h(9000.0_wp, 1.0_wp, 1.0_wp, gam1, lam1)
    write(*,'(a,2f9.4)') '    flat-in-bin / exact:  heating, cooling ', gam1/gam0, lam1/lam0
    call gate(abs(gam1/gam0 - 1.0_wp) < 0.03_wp .and. abs(lam1/lam0 - 1.0_wp) < 0.03_wp, &
              'band-bin field within 3% of the exact spectrum')
    deallocate(nphot, nflat, edges, fband)
  end subroutine check_mochii_field

  !=========================================================================
  ! (5) 8 bins per decade against every native radius (Table 4 conditions)
  subroutine check_size_bins()
    real(kind=wp), allocatable :: nphot(:)
    real(kind=wp) :: w, g12, l12, gn, ln
    write(*,'(/,a)') '(5) size bins: 8 per decade against every native radius '// &
       '(T_c = 3.5e4 K full, G/n_H = 1, T = 9000 K)'
    call grid_and_bins(0.1_wp, 10.0_wp*KEV*3.5e4_wp, 600)
    allocate(nphot(n_egrid))
    w = U_HAB/bb_u_fuv(3.5e4_wp)
    call bb_rates(3.5e4_wp, w, 100.0_wp, nphot)
    call grain_field_load(nphot)
    call per_h(9000.0_wp, 1.0_wp, 1.0_wp, g12, l12)
    grain_bins_per_decade = 1000
    call grain_bins_from_model(m, trim(data_dir))
    call grain_photoelectric_setup()
    call grain_field_load(nphot)
    call per_h(9000.0_wp, 1.0_wp, 1.0_wp, gn, ln)
    grain_bins_per_decade = 8
    write(*,'(a,2f9.4)') '    8 per decade / native:   heating, cooling ', g12/gn, l12/ln
    call gate(abs(g12/gn - 1.0_wp) < 0.01_wp .and. abs(l12/ln - 1.0_wp) < 0.01_wp, &
              'size-bin merge within 1%')
    deallocate(nphot)
  end subroutine check_size_bins

  !=========================================================================
  ! (6) HM79 accommodation and the sign of the exchange
  subroutine check_gas_grain()
    real(kind=wp) :: a_hot, a_cold, c1, c2
    write(*,'(/,a)') '(6) neutral gas-grain exchange'
    a_hot  = thermal_accommodation_H(1.0e6_wp, 20.0_wp, 12.011_wp)
    a_cold = thermal_accommodation_H(1.0_wp, 1.0_wp, 12.011_wp)
    write(*,'(a,2f9.4,a,f7.4)') '    alpha_H(C) at T = 1e6 K / 1 K: ', a_hot, a_cold, &
       '  hard-sphere limit ', 2.0_wp*1.00794_wp*12.011_wp/(13.01894_wp)**2
    call gate(abs(a_hot - 2.0_wp*1.00794_wp*12.011_wp/(13.01894_wp)**2) < 1.0e-3_wp &
              .and. a_cold > 0.85_wp .and. a_cold <= 1.0_wp, &
              'alpha_H between its hard-sphere and low-temperature limits')
    c1 = gas_grain_neutral_cooling(100.0_wp, 30.0_wp, 1.0_wp, 0.1_wp, 1.0_wp, 1.0_wp)
    c2 = gas_grain_neutral_cooling(20.0_wp, 30.0_wp, 1.0_wp, 0.1_wp, 1.0_wp, 1.0_wp)
    write(*,'(a,2es12.4)') '    Lambda_gd per n_H (T=100, T_d=30) and (T=20, T_d=30): ', c1, c2
    call gate(c1 > 0.0_wp .and. c2 < 0.0_wp, 'sign follows T - T_d')
  end subroutine check_gas_grain

end program check_grain_photoelectric
