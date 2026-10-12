module hei_metastable_mod
!---------------------------------------------------------------------------
! MoCHII: the He I 2^3S (triplet metastable) population.
!
! n_3 is a ~1e-7 trace of He, so its balance is LINEAR in n_3 and closes in
! each leaf as a quotient (hei_triplet_population):
!
!   n_3 = n_e [ alpha_3 n_HeII + (q_13 + q_feed) n_HeI ] / D
!   D   = A_31 + Phi_3 + n_e (q_31g + q_31a + q_31b + q_31s + C_3) + n_HI Q_31
!
!   alpha_3  recombination into the triplets, f_3S(T) alpha_B(He II): the
!            triplet share of the case-B cascade (hei_cascade_mod) times the
!            balance's own case-B coefficient (recomb_mod).  Every triplet
!            recombination ends in 2^3S; ground captures feed no triplet, so
!            it is alpha_B whatever par%recombination_continuum is.
!   q_13     e-impact 1^1S -> 2^3S
!   q_feed   e-impact 1^1S -> n^3L (every other triplet term n <= 5), each
!            atom cascading within the triplets to 2^3S
!   A_31     2^3S -> 1^1S radiative (1.272e-4 s^-1, 19.82 eV photon)
!   Phi_3    2^3S photoionization (threshold 4.768 eV; the cross-section fit
!            starts at 4.78 eV), a path estimator
!   q_31g    e-impact 2^3S -> 1^1S, the detailed-balance inverse of q_13
!   q_31a    e-impact 2^3S -> 2^1S;  q_31b  2^3S -> 2^1P
!   q_31s    e-impact 2^3S -> n^1L, n >= 3 (into the singlets)
!   C_3      e-impact ionization of 2^3S, Black (1981, MNRAS 197, 553) eq. 12
!   Q_31     Penning + associative ionization He(2^3S) + H0 (Garcia Munoz
!            2025, A&A 698, A199, Fig. 5)
! and q_33, e-impact 2^3S -> n^3L (2^3P and above), which returns to 2^3S
! by radiative cascade (10830 A for 2^3P): no population change, an energy
! loss.
!
! The electron-impact rates are q = beta Upsilon exp(-dE/kT)/(g_l sqrt T)
! with the Bray, Burgess, Fursa & Tully (2000, A&AS 146, 481) effective
! collision strengths (Tables 2 and 3, n <= 5), beta = 8.629132e-6 exactly,
! observed term energies; data/atomic/hei_metastable.txt carries them
! (tools/fitting/make_hei_metastable.py, which also documents the treatment
! below 10^3.75 K) and this module evaluates them.
!
! One population, four users.  hei_triplet_population is the one n_3 of
!  (1) the output diagnostic (par%hei_metastable, hei_metastable_run),
!  and, with explicitly enabled par%hei_metastable_coupling (default off),
!  (2) the ionization balance (ion_balance_mod::solve_ion_cell): C_3 and
!      Phi_3 return He to He+ -- a share u = (n_e C_3 + Phi_3)/D of every
!      triplet recombination and collisional triplet excitation is undone --
!      and Q_31 ionizes H;
!  (3) the thermal balance (thermal_mod::net_rate), through
!      hei_triplet_net_heating, whose energy ledger is derived below;
!  (4) the He I diffuse channel (diffuse_mod): n_3 A_31 photons at 19.82 eV,
!      and the collisional transfers into 2^1S / 2^1P.
!
! ENERGY LEDGER (hei_triplet_net_heating), energies above 1^1S:
!   E_3 = E(2^3S) = 19.8196 eV, IP_3 = IP(He) - E_3 = 4.7678 eV.
!   excitation 1^1S -> 2^3S, n^3L (rate n_e n_HeI q)    gas loses E(n^3L)
!       [the E(n^3L) - E_3 of a feed excitation leaves as triplet lines]
!   2^3S -> 1^1S by electrons (n_e n_3 q_31g)          gas gains E_3
!   2^3S -> 2^1S, 2^1P, n^1L, n^3L (n_e n_3 q)          gas loses dE
!       [the atom then radiates E(upper): two-photon, 584 A, cascades]
!   C_3 (n_e n_3 C_3)                                  gas loses IP_3
!   Phi_3 (n_3 Phi_3)                                  gas gains h nu - IP_3
!       [heating estimator Sum path_lum sigma_3 (1 - IP_3/E)]
!   Penning / associative (n_HI n_3 Q_31)              gas gains E_3 - IP(H)
!   A_31 (n_3 A_31)                                     no gas change
! so, per unit volume,
!   H_3 = n_3 [ n_e (q_31g E_3 - sum_{31a,31b,31s,33} q dE - C_3 IP_3)
!               + Gamma_3 + n_HI Q_31 (E_3 - IP(H)) ]
!         - n_e n_HeI sum_{13, feed} q E(n^3L) .
! It closes: with n_3 steady, recombination E_3 + absorbed photon energy
! = H_3 + emitted photon energy + ionization energy stored (IP(He) for C_3
! and Phi_3, IP(H) for Q_31), to rounding (gate tests/hei_metastable_coupling).
! What is NOT in it: excitation of the singlets from 1^1S (the open He I
! collisional-excitation cooling item), and the recombination kinetic
! energy, already in cooling_mod's beta on the balance's own alpha.
!
! Phi_3 is a path estimator accumulated during transport at each packet's
! own energy, Sum path_lum sigma_3(E)/(h nu), with its heating companion
! Sum path_lum sigma_3(E) (1 - IP_3/E) (hei_metastable_path_score, called
! from ion_score_mod), turned into rates per leaf at the reduction.  It
! covers whatever band is transported: without par%add_fuv the band starts
! at 13.6 eV and the [4.78, 13.6) eV part of the integral is missing, which
! for a 40 kK star is ~10x the part it counts.  So when this population is
! used, an input that leaves add_fuv and efuv_min out transports the FUV
! extension from 4.78 eV (setup.f90), launched as its own stratum
! (par%fuv_sampling = 'separate') so the ionizing band keeps its packets.  The 2^3S opacity is not
! carried in the transport (n_3 sigma_3 against n_HI sigma_H is ~1e-5 in the
! ionized gas), so the photons Phi_3 absorbs are not removed from the field.
!
! CAVEATS (see md/HEI_10833_PLAN.md sec. 8):
!  (i)   local equilibrium: the advection term (v df_3/dr) is absent - valid
!        in dense H II regions (2^3S destruction time short vs any flow time),
!        not a faithful reproduction of an outflow benchmark;
!  (ii)  10833 can be optically thick in dense metastable regions (a line-
!        transfer treatment of the absorption is out of scope here - this
!        module gives n_3, not tau);
!  (iii) the associative 10% of Q_31 (Garcia Munoz 2025) forms HeH+, which
!        MoCHII does not carry; it is counted as H0 -> H+ like the Penning
!        90% (the ion carries the hydrogen's charge) and deposits E_3 - IP(H),
!        omitting the HeH+ binding.
!  (iv)  the 2^1S / 2^1P populations are not solved: a collisional transfer
!        into them is taken to decay radiatively (A_2q = 51 s^-1 against
!        n_e q(2^1S -> 2^3S) ~ 2e-7 n_e s^-1), so no singlet atom returns
!        to 2^3S.
!---------------------------------------------------------------------------
  use define
  use mpi
  use memory_mod,   only : create_shared_mem, destroy_shared_mem, get_window, &
                          shared_host_rank, shared_host_size, shared_host_comm
  use octree_mod,   only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
  use photo_xsec,   only : sigma_vfky96
  implicit none
  private

  public :: hei_metastable_run, hei_metastable_load
  public :: hei_metastable_path_setup, hei_metastable_path_reset
  public :: hei_metastable_path_score, hei_metastable_path_reduce
  public :: hei_triplet_rates_type, hei_triplet_rates
  public :: hei_triplet_destruction, hei_triplet_population
  public :: hei_triplet_net_heating
  public :: hei_phi3_rate, hei_phi3_heating, hei_phi3_assign
  public :: hei_upsilon

  !--- ionization potentials [eV] (NIST ASD), as in ion_balance_mod
  real(kind=wp), parameter, public :: IP_HE_EV = 24.587387_wp
  real(kind=wp), parameter, public :: IP_H_EV  = 13.598434_wp
  !--- roles of the electron-impact rows
  integer, parameter, public :: ROLE_Q13 = 1, ROLE_FEED = 2, ROLE_Q31A = 3, &
                                ROLE_Q31B = 4, ROLE_SINGLET = 5, ROLE_TRIPLET = 6

  !--- the rate coefficients at one temperature [cm^3 s^-1]; cool13 and
  !--- cool31 the energy-weighted sums [erg cm^3 s^-1] of the excitations
  !--- from 1^1S (Q13 + FEED rows, sum q E(upper)) and from 2^3S
  !--- (Q31A, Q31B, SINGLET, TRIPLET rows, sum q dE).
  type :: hei_triplet_rates_type
     real(kind=wp) :: T      = 0.0_wp
     real(kind=wp) :: alpha3 = 0.0_wp
     real(kind=wp) :: q13    = 0.0_wp, qfeed = 0.0_wp, cool13 = 0.0_wp
     real(kind=wp) :: q31g   = 0.0_wp, q31a  = 0.0_wp, q31b   = 0.0_wp
     real(kind=wp) :: q31s   = 0.0_wp, q33   = 0.0_wp, cool31 = 0.0_wp
     real(kind=wp) :: c3     = 0.0_wp, pen   = 0.0_wp
     real(kind=wp) :: f1s    = 0.0_wp, f1p   = 0.0_wp   ! singlet recomb. branches
  end type hei_triplet_rates_type

  !--- coefficients read from data/atomic/hei_metastable.txt.
  !--- the validity ranges of the C_3 and Q_31 closed forms (C3TRANGE,
  !--- Q31TRANGE); zero until read, which the loader requires
  real(kind=wp) :: hm_c3_T(2) = 0.0_wp, hm_q31_T(2) = 0.0_wp
  real(kind=wp), public, protected :: hm_A31 = 1.272e-4_wp   ! [s^-1]
  real(kind=wp), public, protected :: hm_beta = 0.0_wp   ! [cm^3 s^-1 K^0.5]
  real(kind=wp), public, protected :: hm_E23S = 0.0_wp   ! E(2^3S) [eV]
  !--- electron-impact rows: node grid log10 T, log10 Upsilon at the nodes
  integer, parameter :: MAXROW = 64, MAXNODE = 32
  integer, public, protected :: hm_nrow = 0, hm_nnode = 0
  real(kind=wp), public, protected :: hm_logT(MAXNODE)
  real(kind=wp), public, protected :: hm_lups(MAXNODE, MAXROW)
  real(kind=wp), public, protected :: hm_dE(MAXROW), hm_gl(MAXROW)
  integer,       public, protected :: hm_role(MAXROW), hm_lower(MAXROW), &
                                      hm_upper(3, MAXROW)   ! n, 2S+1, l
  !--- CION3: A, p, T0 of C_3 = A T^p exp(-T0/T)
  real(kind=wp) :: hm_cion3(3) = 0.0_wp
  !--- PENNING: A, c, d1, d2, d3 (Garcia Munoz 2025)
  real(kind=wp) :: hm_pen(5)
  !--- sigma_3 VFKY96 wings [E_th,E_0,s_0,y_a,P,y_w,y_0,y_1] + bridge edges
  real(kind=wp) :: hm_sigA(8), hm_sigB(8), hm_E3, hm_E4
  logical :: hm_loaded = .false.
  !--- Phi_3 path estimators Sum(path_lum sigma_3(E)/(h nu)) and
  !--- Sum(path_lum sigma_3(E) (1 - IP_3/E)), rank local then reduced, and
  !--- the rates [s^-1] and heating [erg s^-1 per 2^3S atom] they give.
  real(kind=wp), allocatable :: phi3_path(:), heat3_path(:)
  !--- the rates are the same on every rank: node-shared (2026-10-12)
  real(kind=wp), pointer, contiguous :: phi3_leaf(:) => null(), heat3_leaf(:) => null()
  integer :: phi3_win(2) = MPI_WIN_NULL

  real(kind=wp), parameter :: kb_eV = kboltz_cgs / ev2erg    ! [eV/K]

contains

  !=========================================================================
  ! Public entry: evaluate n_3 leaf by leaf on the converged state and write
  ! '<base>_heimeta'.
  !=========================================================================
  subroutine hei_metastable_run()
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, &
                              gas_ne, gas_Te
    implicit none
    real(kind=wp), allocatable :: phi3(:), n2s3(:), a3(:), heat3(:)
    real(kind=wp) :: efloor, nH, yHe
    type(hei_triplet_rates_type) :: r
    integer :: il, nleaf

    call hei_metastable_load()

    !--- band-floor note: sigma_3 starts at hm_sigA(1) = 4.78 eV; if the band
    !--- starts above it the soft-UV part of the integral is missed (Phi_3
    !--- underestimated).  The default band starts there (setup.f90).
    efloor = par%eion_min
    if (par%add_fuv) efloor = par%efuv_min
    if (mpar%p_rank == 0 .and. efloor > hm_sigA(1)) &
       write(*,'(a,f7.3,a,f4.2,a)') ' HEIMETA: NOTE: band floor ', efloor, &
       ' eV > the 2^3S edge ', hm_sigA(1), ' eV -> Phi_3 misses the part '// &
       'below the floor; leave add_fuv and efuv_min out of the input (or set '// &
       'add_fuv with efuv_min at the edge) for the full 2^3S photoionization integral.'

    nleaf = amr_grid%nleaf
    allocate(phi3(nleaf), n2s3(nleaf), a3(nleaf), heat3(nleaf))

    yHe = par%abund_He
    do il = 1, nleaf
       phi3(il) = hei_phi3_rate(il)
       heat3(il) = hei_phi3_heating(il)
       nH = gas_nH(il)
       r  = hei_triplet_rates(gas_Te(il))
       a3(il) = r%alpha3
       n2s3(il) = hei_triplet_population(r, gas_ne(il), nH*gas_xHI(il), &
                     nH*yHe*gas_xHeI(il), nH*yHe*gas_xHeII(il), phi3(il))
    end do

    call hei_metastable_write(n2s3, phi3, a3, heat3, nleaf)
    deallocate(phi3, n2s3, a3, heat3)
  end subroutine hei_metastable_run

  !=========================================================================
  ! The rate set at T.  Electron-impact rows: log10 Upsilon linear in log10 T
  ! between the nodes, end values held outside them, Boltzmann factor at the
  ! actual T.  C_3 and Q_31 are each evaluated at T clamped to the range its
  ! closed form is stated for (C3TRANGE: Black 1981, 1e3-1e5 K; Q31TRANGE:
  ! Garcia Munoz 2025, <2% over 200-1e4 K), the end value held outside.
  ! alpha_3 = f_3S alpha_B(He II).
  !=========================================================================
  function hei_triplet_rates(T) result(r)
    use recomb_mod,      only : alphaB_HeII
    use hei_cascade_mod, only : hei_branch_fractions
    real(kind=wp), intent(in) :: T
    type(hei_triplet_rates_type) :: r
    real(kind=wp) :: x, w, ups, q, sqT, kT, Tc, lnT, f3s
    integer :: j, i

    if (.not. hm_loaded) call hei_metastable_load()
    r%T = T
    call hei_branch_fractions(T, f3s, r%f1s, r%f1p)
    r%alpha3 = f3s*alphaB_HeII(T)

    !--- node interval, shared by every row
    x = min(max(log10(T), hm_logT(1)), hm_logT(hm_nnode))
    j = 1
    do while (j < hm_nnode - 1 .and. x > hm_logT(j+1))
       j = j + 1
    end do
    w = (x - hm_logT(j))/(hm_logT(j+1) - hm_logT(j))
    sqT = sqrt(T)
    kT  = kb_eV*T
    do i = 1, hm_nrow
       ups = 10.0_wp**((1.0_wp - w)*hm_lups(j,i) + w*hm_lups(j+1,i))
       if (hm_role(i) == ROLE_Q13) &
          r%q31g = hm_beta*ups/(3.0_wp*sqT)   ! g_1/g_3 exp(+E_13/kT) q_13
       if (hm_dE(i)/kT < 700.0_wp) then
          q = hm_beta*ups*exp(-hm_dE(i)/kT)/(hm_gl(i)*sqT)
       else
          q = 0.0_wp                           ! Boltzmann underflow
       end if
       select case (hm_role(i))
       case (ROLE_Q13)
          r%q13 = r%q13 + q;      r%cool13 = r%cool13 + q*hm_dE(i)*ev2erg
       case (ROLE_FEED)
          r%qfeed = r%qfeed + q;  r%cool13 = r%cool13 + q*hm_dE(i)*ev2erg
       case (ROLE_Q31A)
          r%q31a = r%q31a + q;    r%cool31 = r%cool31 + q*hm_dE(i)*ev2erg
       case (ROLE_Q31B)
          r%q31b = r%q31b + q;    r%cool31 = r%cool31 + q*hm_dE(i)*ev2erg
       case (ROLE_SINGLET)
          r%q31s = r%q31s + q;    r%cool31 = r%cool31 + q*hm_dE(i)*ev2erg
       case (ROLE_TRIPLET)
          r%q33 = r%q33 + q;      r%cool31 = r%cool31 + q*hm_dE(i)*ev2erg
       end select
    end do

    Tc   = min(max(T, hm_c3_T(1)), hm_c3_T(2))
    r%c3 = hm_cion3(1)*Tc**hm_cion3(2)*exp(-hm_cion3(3)/Tc)
    Tc   = min(max(T, hm_q31_T(1)), hm_q31_T(2))
    lnT  = log(Tc)
    r%pen = hm_pen(1)*exp(hm_pen(2)/Tc + lnT*(hm_pen(3) + lnT*(hm_pen(4) &
            + lnT*hm_pen(5))))
  end function hei_triplet_rates

  !=========================================================================
  ! Effective collision strength of row i at T (the interpolation the rate
  ! set uses), for the gates.
  !=========================================================================
  real(kind=wp) function hei_upsilon(i, T) result(ups)
    integer,       intent(in) :: i
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: x, w
    integer :: j
    x = min(max(log10(T), hm_logT(1)), hm_logT(hm_nnode))
    j = 1
    do while (j < hm_nnode - 1 .and. x > hm_logT(j+1))
       j = j + 1
    end do
    w = (x - hm_logT(j))/(hm_logT(j+1) - hm_logT(j))
    ups = 10.0_wp**((1.0_wp - w)*hm_lups(j,i) + w*hm_lups(j+1,i))
  end function hei_upsilon

  !=========================================================================
  ! D [s^-1]: the total destruction rate of one 2^3S atom.
  !=========================================================================
  pure real(kind=wp) function hei_triplet_destruction(r, ne, nHI, phi3) result(D)
    type(hei_triplet_rates_type), intent(in) :: r
    real(kind=wp), intent(in) :: ne, nHI, phi3
    D = hm_A31 + phi3 + ne*(r%q31g + r%q31a + r%q31b + r%q31s + r%c3) &
        + nHI*r%pen
  end function hei_triplet_destruction

  !=========================================================================
  ! n_3 [cm^-3], the steady 2^3S density (D >= A_31 > 0 always).
  !=========================================================================
  pure real(kind=wp) function hei_triplet_population(r, ne, nHI, nHeI, nHeII, &
                                                    phi3) result(n3)
    type(hei_triplet_rates_type), intent(in) :: r
    real(kind=wp), intent(in) :: ne, nHI, nHeI, nHeII, phi3
    n3 = ne*(r%alpha3*nHeII + (r%q13 + r%qfeed)*nHeI) &
         / hei_triplet_destruction(r, ne, nHI, phi3)
  end function hei_triplet_population

  !=========================================================================
  ! Net gas heating by the 2^3S processes [erg cm^-3 s^-1] (the ledger in
  ! the module header); heat3 = Gamma_3, the photoionization heating of one
  ! 2^3S atom [erg s^-1].
  !=========================================================================
  pure real(kind=wp) function hei_triplet_net_heating(r, ne, nHI, nHeI, n3, &
                                                     heat3) result(h)
    type(hei_triplet_rates_type), intent(in) :: r
    real(kind=wp), intent(in) :: ne, nHI, nHeI, n3, heat3
    real(kind=wp) :: e3, ip3, eh
    e3  = hm_E23S*ev2erg
    ip3 = (IP_HE_EV - hm_E23S)*ev2erg
    eh  = (hm_E23S - IP_H_EV)*ev2erg
    h = n3*(ne*(r%q31g*e3 - r%cool31 - r%c3*ip3) + heat3 + nHI*r%pen*eh) &
        - ne*nHeI*r%cool13
  end function hei_triplet_net_heating

  !=========================================================================
  ! Phi_3 [s^-1] and its heating [erg s^-1 per 2^3S atom] of leaf il, from
  ! the reduced path estimators of the last transport pass (0 before any).
  !=========================================================================
  real(kind=wp) function hei_phi3_rate(il) result(p)
    integer, intent(in) :: il
    p = 0.0_wp
    if (associated(phi3_leaf)) p = phi3_leaf(il)
  end function hei_phi3_rate

  real(kind=wp) function hei_phi3_heating(il) result(h)
    integer, intent(in) :: il
    h = 0.0_wp
    if (associated(heat3_leaf)) h = heat3_leaf(il)
  end function hei_phi3_heating

  !--- set the leaf rates directly (gate programs, which run no transport,
  !--- and a restart, previous_state_mod::restore_hei_phi3).  Collective on
  !--- the node; every rank passes the same values, the node's rank 0
  !--- stores them.
  subroutine hei_phi3_assign(phi3, heat3)
    real(kind=wp), intent(in) :: phi3(:), heat3(:)
    call phi3_leaf_windows(size(phi3))
    if (shared_host_rank() == 0) then
       phi3_leaf  = phi3
       heat3_leaf = heat3
    end if
    call phi3_leaf_sync()
  end subroutine hei_phi3_assign

  !--- the node-shared leaf rates at n leaves (re)created when needed
  subroutine phi3_leaf_windows(n)
    integer, intent(in) :: n
    if (associated(phi3_leaf)) then
       if (size(phi3_leaf) /= n) call phi3_leaf_free()
    end if
    if (.not. associated(phi3_leaf)) then
       call create_shared_mem(phi3_leaf,  [n]);  call get_window(phi3_leaf,  phi3_win(1))
       call create_shared_mem(heat3_leaf, [n]);  call get_window(heat3_leaf, phi3_win(2))
    end if
  end subroutine phi3_leaf_windows

  subroutine phi3_leaf_free()
    if (.not. associated(phi3_leaf)) return
    call destroy_shared_mem(phi3_leaf)
    call destroy_shared_mem(heat3_leaf)
    phi3_win = MPI_WIN_NULL
  end subroutine phi3_leaf_free

  subroutine phi3_leaf_sync()
    integer :: ierr
    call MPI_WIN_FENCE(0, phi3_win(1), ierr)
    call MPI_WIN_FENCE(0, phi3_win(2), ierr)
    call MPI_BARRIER(shared_host_comm(), ierr)
  end subroutine phi3_leaf_sync

  !=========================================================================
  ! Phi_3 path estimators: allocated with the other transport estimators
  ! (ion_score_setup, also on re-refinement), zeroed at the start of each
  ! transport pass, scored for every ionizing- or FUV-band segment at the
  ! packet's own energy, and reduced with the rates.
  !=========================================================================
  subroutine hei_metastable_path_setup(nleaf)
    integer, intent(in) :: nleaf
    call hei_metastable_load()
    if (allocated(phi3_path)) deallocate(phi3_path, heat3_path)
    call phi3_leaf_free()
    allocate(phi3_path(nleaf), heat3_path(nleaf))
    phi3_path  = 0.0_wp
    heat3_path = 0.0_wp
  end subroutine hei_metastable_path_setup

  subroutine hei_metastable_path_reset()
    if (allocated(phi3_path)) then
       phi3_path  = 0.0_wp
       heat3_path = 0.0_wp
    end if
  end subroutine hei_metastable_path_reset

  subroutine hei_metastable_path_score(il, E, path_lum)
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: E, path_lum
    real(kind=wp) :: s
    if (E <= hm_sigA(1)) return                       ! below 4.78 eV
    s = path_lum*hei_sigma3(E)
    phi3_path(il)  = phi3_path(il)  + s/(E*ev2erg)
    heat3_path(il) = heat3_path(il) + s*max(1.0_wp - (IP_HE_EV - hm_E23S)/E, 0.0_wp)
  end subroutine hei_metastable_path_score

  !--- reduce, then turn Sum(path_lum coeff) into Int 4 pi J_nu coeff dnu
  !--- per leaf, as for Gamma_i.
  subroutine hei_metastable_path_reduce()
    integer :: ierr, il, n
    real(kind=wp) :: fac
    if (.not. allocated(phi3_path)) return
    n = size(phi3_path)
    call MPI_ALLREDUCE(MPI_IN_PLACE, phi3_path, n, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, heat3_path, n, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    call phi3_leaf_windows(n)
    do il = shared_host_rank() + 1, n, shared_host_size()
       fac = 1.0_wp/((2.0_wp*leaf_half(il))**3*par%distance2cm**2)
       phi3_leaf(il)  = phi3_path(il)*fac
       heat3_leaf(il) = heat3_path(il)*fac
    end do
    call phi3_leaf_sync()
  end subroutine hei_metastable_path_reduce

  !=========================================================================
  ! sigma_3(E) [cm^2]: two VFKY96 wings joined by a log-linear bridge between
  ! the bridge edges hm_E3 (Cooper minimum) and hm_E4 (bump).  Below threshold
  ! (4.78 eV) wing A returns 0 by construction.
  !=========================================================================
  real(kind=wp) function hei_sigma3(E) result(sig)
    real(kind=wp), intent(in) :: E
    real(kind=wp) :: logE, lx3, lx4, y3, y4, mb
    if (E <= 0.0_wp) then
       sig = 0.0_wp;  return
    end if
    logE = log10(E)
    lx3  = log10(hm_E3);  lx4 = log10(hm_E4)
    if (logE <= lx3) then
       sig = sigma_vfky96(E, hm_sigA(1), hm_sigA(2), hm_sigA(3), hm_sigA(4), &
                          hm_sigA(5), hm_sigA(6), hm_sigA(7), hm_sigA(8))
    else if (logE >= lx4) then
       sig = sigma_vfky96(E, hm_sigB(1), hm_sigB(2), hm_sigB(3), hm_sigB(4), &
                          hm_sigB(5), hm_sigB(6), hm_sigB(7), hm_sigB(8))
    else
       y3 = log10(sigma_vfky96(hm_E3, hm_sigA(1), hm_sigA(2), hm_sigA(3), &
                  hm_sigA(4), hm_sigA(5), hm_sigA(6), hm_sigA(7), hm_sigA(8)))
       y4 = log10(sigma_vfky96(hm_E4, hm_sigB(1), hm_sigB(2), hm_sigB(3), &
                  hm_sigB(4), hm_sigB(5), hm_sigB(6), hm_sigB(7), hm_sigB(8)))
       mb  = (y4 - y3)/(lx4 - lx3)
       sig = 10.0_wp**(mb*(logE - lx3) + y3)
    end if
  end function hei_sigma3

  !=========================================================================
  ! Parse data/atomic/hei_metastable.txt (keyword lines; '#' comments).
  !=========================================================================
  subroutine hei_metastable_load()
    use mpi
    implicit none
    character(len=PATH_LEN) :: fname
    character(len=512) :: line
    character(len=16)  :: key, rolename
    real(kind=wp) :: ups(MAXNODE), gl, dE
    integer :: unit, ios, ierr, lower, n, mult, l
    logical :: old_layout

    if (hm_loaded) return
    fname = trim(par%atomic_dir)//'/hei_metastable.txt'
    open(newunit=unit, file=trim(fname), status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(3a)') 'ERROR: cannot open ', &
          trim(fname), ' (par%atomic_dir)'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    old_layout = .false.
    hm_nrow = 0;  hm_nnode = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) key
       select case (trim(key))
       case ('C3TRANGE');  read(line,*) key, hm_c3_T
       case ('Q31TRANGE'); read(line,*) key, hm_q31_T
       case ('TRANGE');    old_layout = .true.   ! one range for both forms
       case ('A31');       read(line,*) key, hm_A31
       case ('BETA');      read(line,*) key, hm_beta
       case ('E23S');      read(line,*) key, hm_E23S
       case ('UPSGRID')
          read(line,*) key, hm_nnode
          if (hm_nnode < 2 .or. hm_nnode > MAXNODE) then
             old_layout = .true.;  cycle
          end if
          read(line,*) key, hm_nnode, hm_logT(1:hm_nnode)
       case ('UPSILON')
          if (hm_nnode < 2 .or. hm_nrow >= MAXROW) then
             old_layout = .true.;  cycle
          end if
          read(line,*) key, rolename, lower, n, mult, l, gl, dE, ups(1:hm_nnode)
          hm_nrow = hm_nrow + 1
          select case (trim(rolename))
          case ('Q13');     hm_role(hm_nrow) = ROLE_Q13
          case ('FEED');    hm_role(hm_nrow) = ROLE_FEED
          case ('Q31A');    hm_role(hm_nrow) = ROLE_Q31A
          case ('Q31B');    hm_role(hm_nrow) = ROLE_Q31B
          case ('SINGLET'); hm_role(hm_nrow) = ROLE_SINGLET
          case ('TRIPLET'); hm_role(hm_nrow) = ROLE_TRIPLET
          case default;     old_layout = .true.
          end select
          hm_lower(hm_nrow) = lower
          hm_upper(:,hm_nrow) = [n, mult, l]
          hm_gl(hm_nrow) = gl
          hm_dE(hm_nrow) = dE
          hm_lups(1:hm_nnode,hm_nrow) = log10(ups(1:hm_nnode))
       case ('CION3');     read(line,*) key, hm_cion3
       case ('PENNING');   read(line,*) key, hm_pen
       case ('SIGMA3A');   read(line,*) key, hm_sigA
       case ('SIGMA3B');   read(line,*) key, hm_sigB
       case ('SIGBRIDGE'); read(line,*) key, hm_E3, hm_E4
       case ('LINE10833'); ! Stage 2 line constants, read by the Python tools
       case ('ALPHA3', 'ALPHA1', 'QEXC13', 'QEXC31A', 'QEXC31B', 'Q10833')
          old_layout = .true.          ! the pre-2026-09-25 file
       end select
    end do
    close(unit)
    !--- the recombination source is alpha_B(He II) f_3S since 2026-09-25 and
    !--- the electron-impact rates are the Bray et al. (2000) rows; a file
    !--- still carrying ALPHA3 / QEXC rows, or missing a row set, would be
    !--- read silently wrong.
    if (old_layout .or. hm_nrow == 0 .or. hm_beta <= 0.0_wp .or. &
        hm_E23S <= 0.0_wp .or. hm_cion3(1) <= 0.0_wp .or. &
        hm_c3_T(1) <= 0.0_wp .or. hm_q31_T(1) <= 0.0_wp .or. &
        count(hm_role(1:hm_nrow) == ROLE_Q13) /= 1) then
       if (mpar%p_rank == 0) write(*,'(3a)') 'ERROR: ', trim(fname), &
          ' is not in the current UPSILON/BETA/E23S layout; regenerate it'// &
          ' with tools/fitting/make_hei_metastable.py'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    !--- the default FUV floor (par%efuv_min when the input leaves it out) is
    !--- define.f90's E_HEI_23S_SIGMA_EDGE, set before this file is read; it
    !--- has to be the edge sigma_3 is evaluated with, or the default band
    !--- would start above (or below) the photoionization it exists to cover.
    if (abs(hm_sigA(1) - E_HEI_23S_SIGMA_EDGE) > 1.0e-9_wp*E_HEI_23S_SIGMA_EDGE) then
       if (mpar%p_rank == 0) write(*,'(a,f10.6,a,f10.6,a)') &
          'ERROR: the SIGMA3A edge ', hm_sigA(1), ' eV differs from '// &
          'E_HEI_23S_SIGMA_EDGE = ', E_HEI_23S_SIGMA_EDGE, &
          ' eV in define.f90 (the default FUV floor); make them agree.'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    hm_loaded = .true.
    if (mpar%p_rank == 0) write(*,'(a,i3,2a)') &
       ' HEIMETA: 2^3S coefficients (', hm_nrow, ' Bray et al. 2000 rows) '// &
       'loaded from ', trim(fname)
  end subroutine hei_metastable_load

  !=========================================================================
  ! Write '<base>_heimeta' (rank 0): n_2s3, Phi_3, alpha_3, Gamma_3, n_H,
  ! n_e, T_e + leaf geometry.
  !=========================================================================
  subroutine hei_metastable_write(n2s3, phi3, a3, heat3, nleaf)
    use iofile_mod
    use utility,       only : get_base_name, fatal_error
    use gas_state_mod, only : gas_nH, gas_ne, gas_Te
    implicit none
    integer,       intent(in) :: nleaf
    real(kind=wp), intent(in) :: n2s3(nleaf), phi3(nleaf), a3(nleaf), heat3(nleaf)
    type(io_file_type) :: file
    character(len=PATH_LEN) :: filename
    real(kind=wp), allocatable :: lxyz(:,:), tmp(:)
    integer :: status, il

    if (mpar%p_rank /= 0) return
    status = 0
    filename = trim(get_base_name(par%out_file))//'_heimeta'// &
               trim(io_file_extension(par%file_format))
    call io_open_new(file, trim(filename), status)

    call io_append_image(file, n2s3, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','n_2s3','He I 2^3S density [cm^-3]',status)
    call io_put_keyword(file,'DIST_CM',par%distance2cm,'distance unit (cm)',status)
    call io_put_keyword(file,'A31',hm_A31,'2^3S->1^1S radiative rate [s^-1]',status)

    call io_append_image(file, phi3, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Phi_3','2^3S photoionization rate [s^-1]',status)

    call io_append_image(file, a3, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','alpha_3', &
         'recombination into 2^3S, f_3S alpha_B(He II) [cm^3 s^-1]',status)

    call io_append_image(file, heat3, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','Gamma_3', &
         '2^3S photoionization heating [erg s^-1 per atom]',status)

    allocate(tmp(nleaf))
    tmp = gas_nH(1:nleaf)
    call io_append_image(file, tmp, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','n_H','hydrogen density [cm^-3]',status)
    tmp = gas_ne(1:nleaf)
    call io_append_image(file, tmp, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','n_e','electron density [cm^-3]',status)
    tmp = gas_Te(1:nleaf)
    call io_append_image(file, tmp, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','T_e','electron temperature [K]',status)
    deallocate(tmp)

    allocate(lxyz(nleaf,3))
    do il = 1, nleaf
       lxyz(il,1) = leaf_cx(il);  lxyz(il,2) = leaf_cy(il);  lxyz(il,3) = leaf_cz(il)
    end do
    call io_append_image(file, lxyz, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','LeafXYZ','leaf center x,y,z (code units)',status)
    deallocate(lxyz)

    allocate(tmp(nleaf))
    do il = 1, nleaf
       tmp(il) = 2.0_wp*leaf_half(il)                  ! full width, code units
    end do
    call io_append_image(file, tmp, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','LeafSize', &
         'leaf full width (code units); V = (size*DIST_CM)^3',status)
    deallocate(tmp)

    call io_close(file, status)
    if (status /= 0) call fatal_error('failed to write He I metastable file: '//trim(filename))
    write(*,'(2a)') ' HEIMETA: 2^3S diagnostic written to: ', trim(filename)
  end subroutine hei_metastable_write

end module hei_metastable_mod
