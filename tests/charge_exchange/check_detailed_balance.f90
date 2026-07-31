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
! Nitrogen's two rows are numerical polynomial fits with no separated factor to
! misplace, and rates too small to matter, so they are screened out below.  For
! oxygen, the near-resonant case, it puts exp(-227/T) on the "O+ + H" row
! instead: the right magnitude — the printed 227 matches the 227.7 K defect to
! 0.3% — on the wrong side of the reaction.
! Read as written the oxygen pair failed this test by 1.47x; with the two
! directions exchanged it passes at 0.91, and the rates then match Cloudy c25's
! independent fits to 1-3% over 4000-20000 K.  Reinstating the transcription as
! printed makes this gate fail again.
!
! Pairs where the smaller of the two coefficients is below K_FLOOR cannot
! affect the ionization balance at nebular densities, and are reported but not
! failed: their fits are not constrained by anything this test could see.
!---------------------------------------------------------------------------
  use define
  use species_mod, only : species_setup, species_cx_coefficients, &
                          n_elements, elem_name, elem_nstage, elem_eth
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
  real(kind=wp), parameter :: K_FLOOR = 1.0e-15_wp   ! cm^3 s^-1
  real(kind=wp), parameter :: TOL     = 0.35_wp      ! |ratio-1| allowed
  !--- Tested from 5000 K up.  Below that the Kingdon & Ferland forms 6 and 7
  !--- are evaluated with t4 CLAMPED into [0.5, 5] by cx_rate, which is the
  !--- documented handling of their ~5e3-5e4 K validity range but necessarily
  !--- breaks detailed balance: the clamp freezes the Boltzmann factor of one
  !--- direction while the other keeps running.  Testing below 5000 K would
  !--- test the clamp, not the data.
  integer,       parameter :: NT = 4
  real(kind=wp), parameter :: TT(NT) = [5.0e3_wp, 8.0e3_wp, 1.0e4_wp, 2.0e4_wp]
  real(kind=wp), parameter :: T_SCREEN = 1.0e4_wp

  real(kind=wp) :: k_i, k_r, dE, req, dev, worst
  character(len=2) :: nm
  integer :: ie, j, k, ierr, ntest, nskip
  logical :: ok, screened

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

  ok = .true.;  worst = 0.0_wp;  ntest = 0;  nskip = 0
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

     !--- screen once per element, at a representative nebular temperature:
     !--- either the pair can affect the ionization balance or it cannot.
     call species_cx_coefficients(ie, 1, T_SCREEN, k_i, k_r)
     if (k_i <= 0.0_wp .or. k_r <= 0.0_wp) cycle
     screened = min(k_i, k_r) < K_FLOOR
     if (screened) then
        nskip = nskip + 1
        call species_cx_coefficients(ie, 1, T_SCREEN, k_i, k_r)
        req = (real(GW1(k), wp)*2.0_wp)/(real(GW0(k), wp)) * exp(-dE/T_SCREEN)
        write(*,'(a4,a10,es12.3,es12.3,es12.3,es12.3,f10.3,a)') nm, 'screened', &
           k_i, k_r, k_i/k_r, req, k_i/k_r/req, '   min k < floor'
        cycle
     end if

     do j = 1, NT
        call species_cx_coefficients(ie, 1, TT(j), k_i, k_r)
        if (k_i <= 0.0_wp .or. k_r <= 0.0_wp) cycle
        req = (real(GW1(k), wp)*2.0_wp) / (real(GW0(k), wp)*1.0_wp) &
              * exp(-dE/TT(j))
        dev = k_i/k_r/req
        ntest = ntest + 1
        worst = max(worst, abs(dev - 1.0_wp))
        if (j == 1 .or. abs(dev-1.0_wp) > TOL) &
           write(*,'(a4,f10.0,es12.3,es12.3,es12.3,es12.3,f10.3)') &
              nm, TT(j), k_i, k_r, k_i/k_r, req, dev
        if (abs(dev - 1.0_wp) > TOL) then
           write(*,'(a,a,a,f8.0,a,f6.3)') '  FAIL: ', trim(nm), ' at T = ', &
              TT(j), ' K, ratio to detailed balance = ', dev
           ok = .false.
        end if
     end do
  end do

  write(*,'(/,a,i0,a,i0,a)') ' tested ', ntest, ' (element, T) points; ', &
     nskip, ' screened below the rate floor'
  write(*,'(a,f6.3,a,f4.2)') ' worst |ratio - 1| = ', worst, &
     '   tolerance ', TOL
  if (ok) then
     write(*,'(a)') ' GATE (charge-exchange detailed balance): PASS'
  else
     write(*,'(a)') ' GATE (charge-exchange detailed balance): FAIL'
     call MPI_FINALIZE(ierr);  stop 1
  end if
  call MPI_FINALIZE(ierr)
end program check_detailed_balance
