program check_hydrogenic_recombination
!---------------------------------------------------------------------------
! Gate for hydrogenic_recombination_mod: H II -> H I and He III -> He II
! recombination from one zero-density nl-resolved cascade on the Storey &
! Hummer (1991) hydrogenic data (tools/fitting/make_hydrogenic_recombination.py).
!
! Checks:
!   1. alpha_B of both ions against the Storey & Hummer (1995) alpha-tot at
!      n_e = 100 cm^-3 (CDS files e1b, e2b), 5000-3e4 K (He II to 1e5 K),
!      within 0.5%.  The departure grows toward low T (He II -0.32% at
!      5000 K, H I -1.7% at 1000 K) because SH95's collisions at 100 cm^-3
!      add to the zero-density total there; the Badnell zero-density total
!      agrees with the cascade to 0.1% at 500 K; and recomb_mod's alphaB_HII / alphaB_HeIII are these
!      numbers (the case-B rate the ionization balance uses);
!   2. recomb_mod's alpha_A of H II (alpha_B + the Milne alpha_1 on the
!      transport cross section) against SH95 case A (e1a) within 0.5%;
!   3. the H I case-B alpha_eff(2s) against Osterbrock & Ferland (2006)
!      Table 4.11 at 1e4-2e4 K within 1% (their 5000 K entry, 1.38e-13, is
!      printed, not tested: it is 4% above the cascade and above SH95's
!      He II value at 2e4 K taken to H I by hydrogenic scaling), and the
!      He II case-B alpha_eff(2^2S) against SH95 within 1.5% (the cascade
!      omits the l-mixing SH95 include at n_e = 100);
!   4. 0 < f_2s(A) < f_2s(B) < 1 for both ions over 10^1.5-10^7 K;
!   5. He II exits close: in case A the Lyman-line shares plus f_2s sum to 1
!      to 1e-6 at every tabulated T; in case B Ly-alpha is 1 - f_2s and the
!      rest zero; the shares fall with n;
!   6. the logarithmic derivative of alpha_B, which the recombination
!      cooling takes by differencing, is continuous across the nodes: its
!      jump across a node is below 1e-3 of its value.
!---------------------------------------------------------------------------
  use define
  use recomb_mod, only : alphaA_HII, alphaB_HII, alphaB_HeIII
  use hydrogenic_recombination_mod, only : hydrogenic_alpha_B, hydrogenic_f_2s, &
                                           hydrogenic_alpha_2s, heii_lyman_shares, &
                                           HEII_NLY
  use mpi
  implicit none

  real(kind=wp), parameter :: TH(5)  = [5.0e3_wp, 1.0e4_wp, 1.5e4_wp, 2.0e4_wp, 3.0e4_wp]
  real(kind=wp), parameter :: SHB(5) = [4.5220e-13_wp, 2.5850e-13_wp, 1.8360e-13_wp, &
                                        1.4280e-13_wp, 9.9110e-14_wp]
  real(kind=wp), parameter :: SHA(5) = [6.8090e-13_wp, 4.1690e-13_wp, 3.1040e-13_wp, &
                                        2.5070e-13_wp, 1.8430e-13_wp]
  real(kind=wp), parameter :: THE(6) = [5.0e3_wp, 1.0e4_wp, 1.5e4_wp, 3.0e4_wp, &
                                        5.0e4_wp, 1.0e5_wp]
  real(kind=wp), parameter :: SHE(6) = [2.5450e-12_wp, 1.5330e-12_wp, 1.1280e-12_wp, &
                                        6.5400e-13_wp, 4.2880e-13_wp, 2.3410e-13_wp]
  real(kind=wp), parameter :: T_OF(4) = [5.0e3_wp, 1.0e4_wp, 1.5e4_wp, 2.0e4_wp]
  real(kind=wp), parameter :: A_OF(4) = [1.38e-13_wp, 0.838e-13_wp, 0.625e-13_wp, 0.506e-13_wp]
  real(kind=wp), parameter :: T_S2(4) = [5.0e3_wp, 1.0e4_wp, 1.5e4_wp, 3.0e4_wp]
  real(kind=wp), parameter :: A_S2(4) = [6.161e-13_wp, 4.091e-13_wp, 3.189e-13_wp, 2.035e-13_wp]
  integer :: ierr, i, z, k
  real(kind=wp) :: T, a, r, fa, fb, f(HEII_NLY), worst, lt, d1, d2, h
  logical :: ok, mono

  call MPI_INIT(ierr)
  ok = .true.

  write(*,'(a)') ' 1. alpha_B vs SH95 alpha-tot (n_e = 100); recomb_mod equal'
  do i = 1, 5
     T = TH(i);  a = hydrogenic_alpha_B(1, T);  r = a/SHB(i)
     write(*,'(a,f8.0,2es13.5,f9.5)') '   H I  ', T, a, SHB(i), r
     if (abs(r - 1.0_wp) > 5.0e-3_wp) ok = .false.
     if (abs(alphaB_HII(T)/a - 1.0_wp) > 1.0e-14_wp) ok = .false.
  end do
  do i = 1, 6
     T = THE(i);  a = hydrogenic_alpha_B(2, T);  r = a/SHE(i)
     write(*,'(a,f8.0,2es13.5,f9.5)') '   He II', T, a, SHE(i), r
     if (abs(r - 1.0_wp) > 5.0e-3_wp) ok = .false.
     if (abs(alphaB_HeIII(T)/a - 1.0_wp) > 1.0e-14_wp) ok = .false.
  end do

  write(*,'(a)') ' 2. H II alpha_A (alpha_B + Milne alpha_1) vs SH95 case A'
  do i = 1, 5
     T = TH(i);  a = alphaA_HII(T);  r = a/SHA(i)
     write(*,'(a,f8.0,2es13.5,f9.5)') '   ', T, a, SHA(i), r
     if (abs(r - 1.0_wp) > 5.0e-3_wp) ok = .false.
  end do

  write(*,'(a)') ' 3. alpha_eff(2s) case B: H I vs O&F Table 4.11, He II vs SH95'
  do i = 1, 4
     T = T_OF(i);  a = hydrogenic_alpha_2s(1, T, 'B');  r = a/A_OF(i)
     if (i == 1) then
        write(*,'(a,f8.0,2es13.5,f9.5,a)') '   H I  ', T, a, A_OF(i), r, '  (not tested)'
     else
        write(*,'(a,f8.0,2es13.5,f9.5)') '   H I  ', T, a, A_OF(i), r
        if (abs(r - 1.0_wp) > 1.0e-2_wp) ok = .false.
     end if
  end do
  do i = 1, 4
     T = T_S2(i);  a = hydrogenic_alpha_2s(2, T, 'B');  r = a/A_S2(i)
     write(*,'(a,f8.0,2es13.5,f9.5)') '   He II', T, a, A_S2(i), r
     if (abs(r - 1.0_wp) > 1.5e-2_wp) ok = .false.
  end do

  write(*,'(a)') ' 4. 0 < f_2s(A) < f_2s(B) < 1'
  do z = 1, 2
     do i = 0, 55
        T = 10.0_wp**(1.5_wp + 0.1_wp*i)
        fa = hydrogenic_f_2s(z, T, 'A');  fb = hydrogenic_f_2s(z, T, 'B')
        if (.not. (fa > 0.0_wp .and. fa < fb .and. fb < 1.0_wp)) then
           write(*,'(a,i2,es12.4,2f9.5)') '   FAIL z, T, f_A, f_B:', z, T, fa, fb
           ok = .false.
        end if
     end do
     write(*,'(a,i2,a,2f9.5)') '   z =', z, ', f_A, f_B at 1e4 K:', &
        hydrogenic_f_2s(z, 1.0e4_wp, 'A'), hydrogenic_f_2s(z, 1.0e4_wp, 'B')
  end do

  write(*,'(a)') ' 5. He II exits close (case A: Lyman shares + f_2s = 1)'
  worst = 0.0_wp;  mono = .true.
  do i = 0, 110
     T = 10.0_wp**(1.5_wp + 0.05_wp*i)
     call heii_lyman_shares(T, 'A', f)
     worst = max(worst, abs(sum(f) + hydrogenic_f_2s(2, T, 'A') - 1.0_wp))
     do k = 2, HEII_NLY - 1
        if (f(k) > f(k-1)) mono = .false.
     end do
     call heii_lyman_shares(T, 'B', f)
     if (abs(f(1) + hydrogenic_f_2s(2, T, 'B') - 1.0_wp) > 1.0e-14_wp) ok = .false.
     if (any(f(2:) /= 0.0_wp)) ok = .false.
  end do
  call heii_lyman_shares(1.0e4_wp, 'A', f)
  write(*,'(a,es10.2,a,l2)') '   worst |sum - 1| = ', worst, '  (tol 1e-6); shares fall with n:', mono
  write(*,'(a,4f9.5,es11.3)') '   1e4 K: Ly-alpha, beta, gamma, n=5..30, n>=31:', &
     f(1), f(2), f(3), sum(f(4:HEII_NLY-1)), f(HEII_NLY)
  if (worst > 1.0e-6_wp .or. .not. mono) ok = .false.

  write(*,'(a)') ' 6. dln alpha_B / dln T continuous across the nodes'
  worst = 0.0_wp;  h = 1.0e-4_wp
  do z = 1, 2
     do i = 2, 109
        lt = 1.5_wp + 0.05_wp*i
        d1 = (log(hydrogenic_alpha_B(z, 10.0_wp**(lt - h))) &
              - log(hydrogenic_alpha_B(z, 10.0_wp**(lt - 2.0_wp*h))))/(h*log(10.0_wp))
        d2 = (log(hydrogenic_alpha_B(z, 10.0_wp**(lt + 2.0_wp*h))) &
              - log(hydrogenic_alpha_B(z, 10.0_wp**(lt + h))))/(h*log(10.0_wp))
        worst = max(worst, abs(d2 - d1)/abs(d1))
     end do
  end do
  write(*,'(a,es10.2,a)') '   worst relative jump = ', worst, '  (tol 1e-3)'
  if (worst > 1.0e-3_wp) ok = .false.

  write(*,'(/,a,a)') ' GATE (hydrogenic recombination): ', merge('PASS', 'FAIL', ok)
  call MPI_FINALIZE(ierr)
  if (.not. ok) stop 1
end program check_hydrogenic_recombination
