program check_detailed_balance
!---------------------------------------------------------------------------
! The two directions of a charge-exchange pair are not independent.  For
!
!     X^0 + H^+  <->  X^+ + H^0 ,      dE = IP(X) - IP(H)
!
! detailed balance fixes the ratio of the rate coefficients at
!
!     k_i / k_r = [g(X+) g(H0)] / [g(X0) g(H+)] * exp(-dE/kT)
!
! with g the ground-term statistical weights.  MoCHII fits the two directions
! independently (Huang et al. 2023 for the first transition, Kingdon & Ferland
! 1996 above it), so nothing enforces this — it has to be checked.
!
! It is checked because it caught a real defect.  Huang et al. (2023, ApJ 951,
! 123) Table 4 is headed "Reactants".  Three species in it exchange with
! hydrogen and have a first ionization potential above hydrogen's — He, O and N
! — but only helium carries an explicit Boltzmann factor exp(-dE_IP/T) in its
! printed rate, so helium is the control: the table puts exp(-12.75/T4) on the
! "He + H+" row, matching |IP(He)-IP(H)|/k = 12.752 in T4 units exactly.
! For oxygen, the near-resonant case, it puts exp(-227/T) on the "O+ + H" row
! instead: the right magnitude — the printed 227 matches the 227.7 K defect to
! 0.3% — on the wrong side of the reaction.
! Read as written the oxygen pair failed this test by 1.47x; with the two
! directions exchanged it passes at 0.91, and the rates then match Cloudy c25's
! independent fits to 1-3% over 4000-20000 K.  Reinstating the transcription as
! printed makes this gate fail again.
!
! It caught a second one on 2026-09-24.  Nitrogen's two Huang rows are lnT
! polynomials under which the endothermic N + H+ falls with T (8.9e-20 at
! 1e4 K) and N+ + H sits 70-86x below Kingdon & Ferland (1996).  This gate
! used to screen out every pair whose smaller rate was below 1e-15 at 1e4 K as
! "too small to affect the balance" -- which let exactly that pair through:
! a rate made small by a defect was taken as evidence that it did not matter.
! There is no such screen now.  Every pair with data in both directions is
! tested, and the only exemption is by name, for a physical reason:
!
!  ONE-SIDED pairs, whose exothermic direction carries a radiative channel
!  that has no collisional reverse.  Detailed balance then bounds the
!  endothermic rate from above -- it is the reverse of part of the forward
!  rate -- so the test is k_exo / k_endo >= its detailed-balance value, i.e.
!  ratio >= 1 - TOL on the exothermic-over-endothermic side.  Carbon is one:
!  C + H+ -> C+ + H is dominated by radiative charge transfer (Stancil et al.
!  1998), and C+ + H -> C + H+ is the non-radiative channel alone.
!---------------------------------------------------------------------------
  use define
  use species_mod, only : species_setup, species_cx_coefficients, &
                          n_elements, elem_name, elem_nstage, elem_eth, &
                          elem_cx_form
  use mpi
  implicit none

  !--- ground-term statistical weights of the neutral and singly ionized
  !--- stage (NIST): C 3P/2P, N 4S/3P, O 3P/4S, Ne 1S/2P, Mg 1S/2S,
  !--- Si 3P/2P, S 3P/4S, Ar 1S/2P, Fe a5D/a6D, Cl 2P/3P, Ca 1S/2S.
  integer, parameter :: NEL = 11
  character(len=2), parameter :: ELN(NEL) = &
     ['c ','n ','o ','ne','mg','si','s ','ar','fe','cl','ca']
  integer, parameter :: GW0(NEL) = [9, 4, 9, 1, 1, 9, 9, 1, 25, 6, 1]
  integer, parameter :: GW1(NEL) = [6, 9, 4, 6, 2, 6, 4, 6, 30, 9, 2]

  real(kind=wp), parameter :: IP_H    = 13.5984_wp   ! eV
  real(kind=wp), parameter :: TOL     = 0.35_wp      ! |ratio-1| allowed
  !--- pairs whose exothermic direction carries radiative charge transfer,
  !--- tested one-sided (see the header), by name
  integer, parameter :: NRADIATIVE_EXO = 1
  character(len=2), parameter :: RADIATIVE_EXO(NRADIATIVE_EXO) = ['c ']
  !--- Tested from 3000 K.  The Kingdon & Ferland forms 6 and 7 clamp only
  !--- their fitted factor into t4 = [0.5, 5]; the Boltzmann factor of form 7
  !--- takes the actual T, so a KF96 pair that shares b, c, d keeps detailed
  !--- balance below 5000 K too.
  integer,       parameter :: NT = 5
  real(kind=wp), parameter :: TT(NT) = [3.0e3_wp, 5.0e3_wp, 8.0e3_wp, &
                                        1.0e4_wp, 2.0e4_wp]

  real(kind=wp) :: k_i, k_r, dE, req, dev, worst, dside, zp, kT_eV, kr0
  !--- S+ terms (NIST / CHIANTI observed, (2J+1)-weighted): 4S, 2D, 2P
  real(kind=wp), parameter :: SP_G(3) = [4.0_wp, 10.0_wp, 6.0_wp]
  real(kind=wp), parameter :: SP_E(3) = [0.0_wp, 1.8439_wp, 3.0446_wp]
  !--- Zhao et al. (2005) Table III, S(3P) + H+ -> S+(4S) + H at 1e4 K
  real(kind=wp), parameter :: ZHAO_4S_1E4 = 8.90e-16_wp
  character(len=2) :: nm
  integer :: ie, j, k, ierr, ntest, nnodata
  logical :: ok, one_sided, fail

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir = 'data/atomic'
  par%use_metals = .true.
  par%abund_C  = 2.2e-4_wp;  par%abund_N  = 4.0e-5_wp
  par%abund_O  = 3.3e-4_wp;  par%abund_Ne = 5.0e-5_wp
  par%abund_S  = 9.0e-6_wp;  par%abund_Ar = 3.0e-6_wp
  par%abund_Mg = 3.0e-6_wp;  par%abund_Fe = 3.0e-6_wp
  par%abund_Si = 3.0e-6_wp;  par%abund_Cl = 3.0e-7_wp
  par%abund_Ca = 2.0e-8_wp
  call species_setup(1)

  ok = .true.;  worst = 0.0_wp;  ntest = 0;  nnodata = 0
  write(*,'(/,a)') '=== detailed balance of the transition-1 charge-exchange pairs ==='
  write(*,'(a4,a10,a12,a12,a12,a12,a10)') 'el', 'T [K]', 'k_cxi', 'k_cxr', &
     'measured', 'required', 'ratio'

  do ie = 1, n_elements
     nm = elem_name(ie)
     k = 0
     do j = 1, NEL
        if (trim(ELN(j)) == trim(nm)) k = j
     end do
     if (k == 0) cycle
     if (elem_nstage(ie) < 2) cycle
     !--- dE from the registry's own threshold, so the gate cannot drift from
     !--- the data it is checking
     dE = (elem_eth(ie, 1) - IP_H) * ev2erg / kboltz_cgs
     one_sided = any(RADIATIVE_EXO == nm)

     !--- no data in one direction (form 0): nothing to test
     call species_cx_coefficients(ie, 1, TT(1), k_i, k_r)
     if (k_i <= 0.0_wp .or. k_r <= 0.0_wp) then
        nnodata = nnodata + 1
        write(*,'(a4,a)') nm, '   no charge-exchange data in one direction'
        cycle
     end if

     !--- resolved by product term (form 9, sulfur): the forward rate sums
     !--- the terms of X+, the reverse weights each by its actual
     !--- population.  Detailed balance of the TOTAL then holds only with
     !--- X+ in LTE, i.e. with the partition function of X+ in place of its
     !--- ground weight, at a density that thermalizes its terms; and at
     !--- n_e -> 0 the reverse is the ground term's channel alone.
     if (elem_cx_form(ie, 1) == 9) then
        if (trim(nm) /= 's') then
           write(*,'(a,a)') '  FAIL: form 9 checked for sulfur only: ', nm
           ok = .false.;  cycle
        end if
        do j = 1, NT
           kT_eV = kboltz_cgs*TT(j)/ev2erg
           zp = sum(SP_G*exp(-SP_E/kT_eV))
           call species_cx_coefficients(ie, 1, TT(j), k_i, k_r, 1.0e12_wp)
           req = (zp*2.0_wp)/(real(GW0(k), wp)*1.0_wp)*exp(-dE/TT(j))
           dev = k_i/k_r/req
           ntest = ntest + 1
           write(*,'(a4,f10.0,es12.3,es12.3,es12.3,es12.3,es10.3,a)') &
              nm, TT(j), k_i, k_r, k_i/k_r, req, dev, '  LTE (n_e = 1e12)'
           if (abs(dev - 1.0_wp) > 0.01_wp) then
              write(*,'(a,f8.0)') '  FAIL: sulfur LTE limit at T = ', TT(j)
              ok = .false.
           end if
        end do
        call species_cx_coefficients(ie, 1, 1.0e4_wp, k_i, k_r)
        kr0 = ZHAO_4S_1E4*real(GW0(k), wp)/(SP_G(1)*2.0_wp) &
              *exp(dE/1.0e4_wp)
        write(*,'(a,es12.4,a,es12.4)') '   s  n_e -> 0, 1e4 K: reverse ', k_r, &
           '   ground-term channel (Zhao 4S) by detailed balance ', kr0
        if (abs(k_r/kr0 - 1.0_wp) > 1.0e-3_wp) then
           write(*,'(a)') '  FAIL: sulfur reverse at n_e -> 0 is not the 4S channel'
           ok = .false.
        end if
        cycle
     end if

     do j = 1, NT
        call species_cx_coefficients(ie, 1, TT(j), k_i, k_r)
        req = (real(GW1(k), wp)*2.0_wp) / (real(GW0(k), wp)*1.0_wp) &
              * exp(-dE/TT(j))
        dev = k_i/k_r/req
        ntest = ntest + 1
        if (one_sided) then
           !--- exothermic over endothermic must be at least detailed balance
           dside = dev
           if (dE > 0.0_wp) dside = 1.0_wp/dev
           fail = dside < 1.0_wp - TOL
        else
           worst = max(worst, abs(dev - 1.0_wp))
           fail = abs(dev - 1.0_wp) > TOL
        end if
        if (j == 1 .or. fail .or. one_sided) &
           write(*,'(a4,f10.0,es12.3,es12.3,es12.3,es12.3,es10.3,a)') &
              nm, TT(j), k_i, k_r, k_i/k_r, req, dev, &
              merge('  one-sided', '           ', one_sided)
        if (fail) then
           write(*,'(a,a,a,f8.0,a,es10.3)') '  FAIL: ', trim(nm), ' at T = ', &
              TT(j), ' K, ratio to detailed balance = ', dev
           ok = .false.
        end if
     end do
  end do

  write(*,'(/,a,i0,a,i0,a)') ' tested ', ntest, ' (element, T) points; ', &
     nnodata, ' elements without data in both directions'
  write(*,'(a,f6.3,a,f4.2)') ' worst |ratio - 1| of the two-sided pairs = ', &
     worst, '   tolerance ', TOL
  if (ok) then
     write(*,'(a)') ' GATE (charge-exchange detailed balance): PASS'
  else
     write(*,'(a)') ' GATE (charge-exchange detailed balance): FAIL'
     call MPI_FINALIZE(ierr);  stop 1
  end if
  call MPI_FINALIZE(ierr)
end program check_detailed_balance
