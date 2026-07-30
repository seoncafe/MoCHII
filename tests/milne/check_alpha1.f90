program check_alpha1
!---------------------------------------------------------------------------
! Gate for the ground-level direct-capture rate alpha_1(T) of recomb_mod.
!
! alpha_1 is what turns the Badnell total alpha_A into case B, and MoCHII takes
! it from the Milne relation applied to the transport photoionization cross
! sections, fitted to the rr_mao form by tools/fitting/fit_alpha1_milne.py.
! This gate re-derives the same integral here, from photo_xsec directly and
! with its own quadrature, so it tests the coefficients rather than trusting
! either the fitting script or milne_recomb_spectrum_mod's internal copy.
!
!   alpha_1(T) = (g_n/2 g_+) (2/sqrt(pi)) (kT)^(-3/2) sqrt(2/m_e)/(m_e c^2)
!                Int E^2 sigma_pi(E) exp[-(E-E_th)/kT] dE
!
! Checks:
!   1. alpha_A - alpha_B (which is alpha_1) reproduces that quadrature over
!      1e3-1e5 K, the fitted range, to TOL_FIT;
!   2. alpha_1 < alpha_A everywhere, i.e. case B stays positive;
!   3. alpha_1 at 1e4 K agrees with Cloudy c25.00's level-resolved Milne value
!      to TOL_CLOUDY.  Cloudy derives H-like/He-like recombination the same
!      way from its own cross sections (source/iso_radiative_recomb.cpp); the
!      numbers below are row 0 (ground level) of its data/h_iso_recomb.dat and
!      data/he_iso_recomb.dat, whose last row reproduces the known alpha_A.
!
! Check 3 is the one that fails if alpha_1 is taken from a recombination data
! set unrelated to sigma_pi: the Mao & Kaastra (2016) level fit used before is
! 2.4% low on He I.
!
! Build:  see tests/milne/run_check.sh
!---------------------------------------------------------------------------
  use define
  use photo_xsec, only : sigma_HI, sigma_HeI, sigma_HeII
  use recomb_mod, only : alphaA_HII, alphaB_HII, alphaA_HeII, alphaB_HeII, &
                         alphaA_HeIII, alphaB_HeIII
  use mpi
  implicit none

  integer,       parameter :: NCH = 3
  integer,       parameter :: NT  = 81
  real(kind=wp), parameter :: TLO = 1.0e3_wp, THI = 1.0e5_wp
  !--- upper limit of the x = (E-E_th)/kT integration, the XMAX of
  !--- milne_recomb_spectrum_mod: exp(-60) of the near-threshold emissivity.
  real(kind=wp), parameter :: XMAX = 60.0_wp
  !--- g_n/(2 g_+): H I(1s) 2/(2*1), He I(1^1S) 1/(2*2), He II(1s) 2/(2*1)
  real(kind=wp), parameter :: GFAC(NCH) = [1.0_wp, 0.25_wp, 1.0_wp]
  !--- 5x the worst fit residual of fit_alpha1_milne.py (0.019%, on H I)
  real(kind=wp), parameter :: TOL_FIT    = 1.0e-3_wp
  real(kind=wp), parameter :: TOL_CLOUDY = 1.0e-2_wp
  !--- Cloudy c25.00 ground-level (row 0) values at 1e4 K [cm^3 s^-1]
  real(kind=wp), parameter :: A1_CLOUDY(NCH) = &
     [1.5840e-13_wp, 1.6032e-13_wp, 6.5145e-13_wp]
  character(len=6), parameter :: NAME(NCH) = &
     [character(len=6) :: 'H I', 'He I', 'He II']

  real(kind=wp) :: eth(NCH)
  real(kind=wp) :: T, a1_code, a1_quad, aA, dev, worst, T_worst, ratio_max
  integer :: ierr, ich, it
  logical :: ok

  call MPI_INIT(ierr)
  par%recomb_model = 'badnell_mao'
  eth = [eth_HI, eth_HeI, eth_HeII]
  ok  = .true.

  !--- 1 + 2: the fitted coefficients against the Milne quadrature, and the
  !--- positivity of case B, over the fitted temperature range.
  write(*,'(a)') '=== alpha_1 vs the Milne quadrature (1e3-1e5 K) ==='
  write(*,'(a8,4a14)') 'channel', 'worst dev', 'at T [K]', &
                       'max a1/aA', 'tolerance'
  do ich = 1, NCH
     worst = 0.0_wp;  T_worst = 0.0_wp;  ratio_max = 0.0_wp
     do it = 1, NT
        T = 10.0_wp**(log10(TLO) + (log10(THI) - log10(TLO)) &
                      *real(it-1, wp)/real(NT-1, wp))
        a1_code = alpha1_from_rates(ich, T)
        a1_quad = alpha1_milne_quad(ich, T)
        aA      = alphaA_of(ich, T)
        dev     = abs(a1_code/a1_quad - 1.0_wp)
        if (dev > worst) then
           worst = dev;  T_worst = T
        end if
        ratio_max = max(ratio_max, a1_code/aA)
        if (a1_code >= aA) then
           write(*,'(3a,f9.0,a)') '  FAIL: ', trim(NAME(ich)), &
              ' alpha_1 >= alpha_A at T =', T, ' K (case B would be negative)'
           ok = .false.
        end if
     end do
     write(*,'(a8,es14.3,f14.0,f14.5,es14.3)') NAME(ich), worst, T_worst, &
        ratio_max, TOL_FIT
     if (worst > TOL_FIT) then
        write(*,'(3a,es10.3)') '  FAIL: ', trim(NAME(ich)), &
           ' alpha_1 departs from the Milne quadrature by ', worst
        ok = .false.
     end if
  end do

  !--- 3: the independent reference at 1e4 K.
  write(*,'(/,a)') '=== alpha_1 at 1e4 K vs Cloudy c25.00 ==='
  write(*,'(a8,4a14)') 'channel', 'MoCHII', 'Milne quad', 'Cloudy c25', &
                       'MoCHII/C25'
  do ich = 1, NCH
     a1_code = alpha1_from_rates(ich, 1.0e4_wp)
     a1_quad = alpha1_milne_quad(ich, 1.0e4_wp)
     write(*,'(a8,3es14.5,f14.5)') NAME(ich), a1_code, a1_quad, &
        A1_CLOUDY(ich), a1_code/A1_CLOUDY(ich)
     dev = abs(a1_code/A1_CLOUDY(ich) - 1.0_wp)
     if (dev > TOL_CLOUDY) then
        write(*,'(3a,f7.3,a)') '  FAIL: ', trim(NAME(ich)), &
           ' alpha_1 is ', 100.0_wp*dev, '% from Cloudy c25.00'
        ok = .false.
     end if
  end do

  write(*,'(/,a)') 'GATE (ground-level alpha_1 from the Milne relation): ' // &
     merge('PASS', 'FAIL', ok)
  call MPI_FINALIZE(ierr)
  if (.not. ok) error stop 1

contains

  !--- alpha_1 as the code uses it: the difference the case-B rate is built
  !--- from, so this exercises alpha1_* through its only caller.
  real(kind=wp) function alpha1_from_rates(ich, T) result(a1)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: T
    select case (ich)
    case (1);      a1 = alphaA_HII(T)   - alphaB_HII(T)
    case (2);      a1 = alphaA_HeII(T)  - alphaB_HeII(T)
    case default;  a1 = alphaA_HeIII(T) - alphaB_HeIII(T)
    end select
  end function alpha1_from_rates

  real(kind=wp) function alphaA_of(ich, T) result(aA)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: T
    select case (ich)
    case (1);      aA = alphaA_HII(T)
    case (2);      aA = alphaA_HeII(T)
    case default;  aA = alphaA_HeIII(T)
    end select
  end function alphaA_of

  real(kind=wp) function sig(ich, E)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: E
    select case (ich)
    case (1);      sig = sigma_HI(E)
    case (2);      sig = sigma_HeI(E)
    case default;  sig = sigma_HeII(E)
    end select
  end function sig

  !=========================================================================
  ! alpha_1(T) [cm^3/s] from the Milne relation, integrated here.  Composite
  ! Simpson in x = (E - E_th)/kT on segments that resolve the near-threshold
  ! peak of E^2 sigma_pi(E) exp(-x).
  !=========================================================================
  real(kind=wp) function alpha1_milne_quad(ich, T_K) result(a1)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: T_K
    integer,       parameter :: NSEG = 6, NPAN = 2000
    real(kind=wp), parameter :: XEDGE(NSEG+1) = &
       [0.0_wp, 0.25_wp, 1.0_wp, 4.0_wp, 12.0_wp, 30.0_wp, XMAX]
    real(kind=wp) :: kT_eV, kT_erg, h, xa, xm, xb, acc
    integer :: is, ip
    kT_eV  = kboltz_cgs*T_K/ev2erg
    kT_erg = kboltz_cgs*T_K
    acc = 0.0_wp
    do is = 1, NSEG
       h = (XEDGE(is+1) - XEDGE(is))/real(NPAN, wp)
       do ip = 1, NPAN
          xa = XEDGE(is) + h*real(ip-1, wp)
          xb = xa + h
          xm = 0.5_wp*(xa + xb)
          acc = acc + h/6.0_wp*( shape_x(ich, xa, kT_eV) &
                               + 4.0_wp*shape_x(ich, xm, kT_eV) &
                               + shape_x(ich, xb, kT_eV) )
       end do
    end do
    !--- acc is Int E^2 sigma exp(-x) dx with E in eV; -> erg^3 cm^2
    acc = acc*kT_eV*ev2erg**3
    a1  = GFAC(ich)*(2.0_wp/sqrt(pi))*kT_erg**(-1.5_wp) &
          *sqrt(2.0_wp/me_cgs)/(me_cgs*clight_cgs**2)*acc
  end function alpha1_milne_quad

  real(kind=wp) function shape_x(ich, x, kT_eV) result(f)
    integer,       intent(in) :: ich
    real(kind=wp), intent(in) :: x, kT_eV
    real(kind=wp) :: E
    E = eth(ich) + x*kT_eV
    f = E*E*sig(ich, E)*exp(-x)
  end function shape_x

end program check_alpha1
