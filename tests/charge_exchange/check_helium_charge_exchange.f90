program check_helium_charge_exchange
!---------------------------------------------------------------------------
! Gate: He-H charge exchange in the H/He balance (par%helium_charge_exchange,
! ion_balance_mod::cx_HeII_H and companions).
!  (a) the rate coefficients are the published forms: the Stancil, Lepp &
!      Dalgarno (1998) reaction (19) radiative fit (Glover & Jappsen 2007
!      R26 until 2026-09-26) plus the non-radiative rate of Loreau et al.
!      (2014) by reciprocity -- the values below are
!      tools/fitting/make_heii_h_charge_exchange.py's direct Maxwell
!      averages at 3000, 8000, 1e4, 3e4 K (1e4 K a table node, the others
!      between nodes, so this also bounds the interpolation) -- and the KF96
!      He++ row;
!  (b) the two reverses are the detailed-balance partners of the forward
!      rates: He + H+ of the non-radiative He+ + H channel alone,
!      k12/k_nr = 4 exp(-10.989 eV/kT), and He+ + H+ of the He++ row,
!      k23/k32 = exp(-40.819 eV/kT);
!  (c) a solve_ion_cell call reaches a state that satisfies all three
!      balances WITH the exchange terms on both sides, in partially neutral
!      gas where they matter, and charge neutrality;
!  (d) the exchange lowers He+ and raises H+ there (sign of the effect), and
!      helium_charge_exchange = .false. removes it.
!---------------------------------------------------------------------------
  use define
  use recomb_mod
  use ion_balance_mod, only : solve_ion_cell, cx_HeII_H, cx_HeI_Hp, &
                              cx_HeIII_H, cx_HeII_Hp
  use mpi
  implicit none
  integer, parameter :: NST = 4
  real(kind=wp), parameter :: GH(NST) = [1.0e-10_wp, 3.0e-12_wp, 3.0e-13_wp, 3.0e-14_wp]
  real(kind=wp), parameter :: TOL = 1.0e-8_wp
  real(kind=wp) :: T, nH, x, h1, h2, ne, dev, worst, knr, ref, gHe1, gHe2
  real(kind=wp) :: aH, aHe2, aHe3, cH, cHe1, cHe2, k21, k12, k32, k23
  real(kind=wp) :: nHI, nHII, nHe, h3, xoff, h1off, h2off, neoff, eh, e1, e2, eq
  real(kind=wp) :: Ts(4) = [3000.0_wp, 8000.0_wp, 1.0e4_wp, 3.0e4_wp]
  !--- Loreau et al. (2014) by reciprocity at Ts (make_heii_h_charge_exchange.py)
  real(kind=wp), parameter :: KNR_REF(4) = [9.753526e-17_wp, 9.720719e-17_wp, &
                                            9.643531e-17_wp, 9.961932e-17_wp]
  integer :: is, it, ierr
  logical :: ok

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%use_metals = .false.
  par%recombination_continuum = 'on_the_spot'
  par%abund_He   = 0.1_wp
  par%ci_model   = 'voronov'
  par%recomb_model = 'detailed_balance'
  ok = .true.

  !--- (a) + (b)
  worst = 0.0_wp
  do it = 1, 4
     T = Ts(it)
     knr = KNR_REF(it)
     ref = 1.20e-15_wp*(T/300.0_wp)**0.25_wp + knr
     worst = max(worst, abs(cx_HeII_H(T)/ref - 1.0_wp))
     ref = 4.0_wp*knr*exp(-10.988953_wp*11604.518_wp/T)
     worst = max(worst, abs(cx_HeI_Hp(T)/ref - 1.0_wp))
     worst = max(worst, abs(cx_HeIII_H(T)/1.0e-14_wp - 1.0_wp))
     if (T >= 1.0e4_wp) worst = max(worst, abs(cx_HeII_Hp(T)/(1.0e-14_wp* &
                 exp(-40.819331_wp*11604.518_wp/T)) - 1.0_wp))
     write(*,'(a,f8.0,a,4es11.3)') '  T =', T, '  k21 k12 k32 k23 =', &
        cx_HeII_H(T), cx_HeI_Hp(T), cx_HeIII_H(T), cx_HeII_Hp(T)
  end do
  write(*,'(a,es9.2)') ' (a,b) published rate forms and both detailed balances: worst ', worst
  !--- 1e-3: the log-log interpolation of the 0.02 dex table between nodes
  !--- (the eV-to-K conversion rounding is 1e-5)
  if (worst > 1.0e-3_wp) ok = .false.

  !--- (c) + (d): partially neutral cells at 1e4 K, n_H = 100
  T = 1.0e4_wp;  nH = 100.0_wp;  nHe = par%abund_He*nH
  gHe1 = 5.0e-11_wp;  gHe2 = 1.0e-13_wp
  aH = alphaB_HII(T);  aHe2 = alphaB_HeII(T);  aHe3 = alphaB_HeIII(T)
  cH = ci_HI(T);  cHe1 = ci_HeI(T);  cHe2 = ci_HeII(T)
  k21 = cx_HeII_H(T);  k12 = cx_HeI_Hp(T);  k32 = cx_HeIII_H(T);  k23 = cx_HeII_Hp(T)
  worst = 0.0_wp
  do is = 1, NST
     par%helium_charge_exchange = .true.
     x = 0.5_wp;  h1 = 0.5_wp;  h2 = 0.4_wp
     call solve_ion_cell(GH(is), gHe1, gHe2, nH, T, .false., x, h1, h2, ne)
     h3 = max(0.0_wp, 1.0_wp - h1 - h2)
     nHI = nH*x;  nHII = nH*(1.0_wp - x)
     eh = abs(nHI*(GH(is) + cH*ne + nHe*(h2*k21 + h3*k32)) &
              - nHII*(aH*ne + nHe*(h1*k12 + h2*k23))) &
          / (nHI*(GH(is) + cH*ne + nHe*(h2*k21 + h3*k32)))
     e1 = abs(h1*(gHe1 + cHe1*ne + k12*nHII) - h2*(aHe2*ne + k21*nHI)) &
          / (h1*(gHe1 + cHe1*ne + k12*nHII))
     e2 = abs(h2*(gHe2 + cHe2*ne + k23*nHII) - h3*(aHe3*ne + k32*nHI)) &
          / (h2*(gHe2 + cHe2*ne + k23*nHII))
     eq = abs(ne - nH*((1.0_wp - x) + par%abund_He*(h2 + 2.0_wp*h3)))/ne
     worst = max(worst, eh, e1, e2, eq)
     par%helium_charge_exchange = .false.
     xoff = 0.5_wp;  h1off = 0.5_wp;  h2off = 0.4_wp
     call solve_ion_cell(GH(is), gHe1, gHe2, nH, T, .false., xoff, h1off, h2off, neoff)
     write(*,'(a,es9.2,a,f9.6,a,es10.3,a,es10.3,a,4es9.1)') '  Gamma_H =', GH(is), &
        '  x_HII =', 1.0_wp - x, '  x_HeII/off =', h2/h2off, &
        '  x_HII/off =', (1.0_wp - x)/(1.0_wp - xoff), '  resid', eh, e1, e2, eq
     if (.not. (h2 < h2off .and. x < xoff)) then
        write(*,'(a)') '  FAIL: the exchange must lower He+ and raise H+'
        ok = .false.
     end if
  end do
  write(*,'(a,es9.2)') ' (c) balances with the exchange on both sides: worst ', worst
  if (worst > TOL) ok = .false.

  if (ok) then
     write(*,'(a)') 'PASS'
  else
     write(*,'(a)') 'FAIL'
     call MPI_FINALIZE(ierr)
     stop 1
  end if
  call MPI_FINALIZE(ierr)
end program check_helium_charge_exchange
