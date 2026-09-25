module hei_metastable_mod
!---------------------------------------------------------------------------
! MoCHII: He I 2^3S (triplet metastable) diagnostic.
!
! On the converged state this evaluates the local 2^3S population n_3 from the
! Oklopcic & Hirata (2018, ApJ 855 L11 = OH18) balance, which is LINEAR in n_3
! (a closed-form quotient per leaf, no iteration):
!
!   n_3 = [ alpha_3 n_e n_HeII + q_13 n_e n_HeI ]
!         / [ A_31 + Phi_3 + n_e (q_31a + q_31b + q_31g + C_3) + n_HI Q_31 ]
!
!   alpha_3     He II -> 2^3S recombination (triplet channel)
!   q_13        e-impact 1^1S -> 2^3S
!   A_31        2^3S -> 1^1S radiative (1.272e-4 s^-1)
!   Phi_3       2^3S photoionization rate (threshold 4.78 eV)
!   q_31a,q_31b e-impact 2^3S -> 2^1S and 2^3S -> 2^1P
!   q_31g       e-impact 2^3S -> 1^1S, the inverse of q_13 by detailed
!               balance, q_31g = q_13 (g_1/g_3) exp(E_13/kT), g_1/g_3 = 1/3
!   C_3         e-impact ionization of 2^3S, Black (1981, MNRAS 197, 553)
!               eq. (12), 8.38e-10 T^(1/2) exp(-55338/T), 1e3-1e5 K
!   Q_31        Penning + associative ionization He(2^3S) + H0
!
! The three excitation rates are q = beta Upsilon(T) exp(-E/kT)/(g_l sqrt T)
! with the exact Maxwellian constant beta = h^2/[(2 pi m_e)^(3/2) k^(1/2)] =
! 8.629132e-6 and the observed level energies (CHIANTI he_1.elvlc); OH18's
! 2.1e-8 sqrt(13.6 eV/kT) is the same constant 3.3% low.
!
! q_31g is not in the OH18 balance, but a rate whose forward direction is
! carried must carry its inverse, or the two levels cannot reach the
! Boltzmann ratio in equilibrium.  Measured against q_31a + q_31b it is 5.6%
! at 1e4 K, 7.2% at 8000 K, 15% at 5000 K and 57% at 3000 K (added
! 2026-09-24).  C_3, also absent from OH18, is 1.0% of q_31a + q_31b at
! 1e4 K, 6.8% at 1.5e4 K and 20% at 2e4 K.
!
! Coefficients are read from data/atomic/hei_metastable.txt (closed-form fits
! transcribed from the EXHALE triplet implementation; the Fortran here only
! evaluates the fits).  Phi_3 is a path estimator accumulated during
! transport at each packet's own energy, Sum path_lum sigma_3(E)/(h nu)
! (hei_metastable_path_score, called from ion_score_mod), the same form as
! gas_rates_mod's Gamma_i; it covers whatever band is transported (add the
! FUV extension down to 4.78 eV for the full integral).
!
! Runs once at output time (par%hei_metastable): a leaf loop that reads the
! converged gas state, no feedback into the ionization / thermal balance (n_3
! is a ~1e-7 trace of He, its energy budget negligible).  Output section file
! '<base>_heimeta' carries n_2s3, Phi_3, n_H, n_e, T_e + leaf geometry.
!
! CAVEATS (see md/HEI_10833_PLAN.md sec. 8):
!  (i)   local equilibrium: the advection term (v df_3/dr) is absent - valid
!        in dense H II regions (2^3S destruction time short vs any flow time),
!        not a faithful reproduction of an outflow benchmark;
!  (ii)  10833 can be optically thick in dense metastable regions (a line-
!        transfer treatment of the absorption is out of scope here - this
!        module gives n_3, not tau);
!  (iii) the n_HI Q_31 sink is the neutral/PDR channel (Penning ionization
!        with H0), included from the start.
!---------------------------------------------------------------------------
  use define
  use octree_mod,   only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
  use photo_xsec,   only : sigma_vfky96
  implicit none
  private

  public :: hei_metastable_run
  public :: hei_metastable_path_setup, hei_metastable_path_reset
  public :: hei_metastable_path_score, hei_metastable_path_reduce

  !--- coefficients read from data/atomic/hei_metastable.txt.
  real(kind=wp) :: hm_Tmin = 1.0e2_wp, hm_Tmax = 3.0e4_wp
  real(kind=wp) :: hm_A31  = 1.272e-4_wp
  real(kind=wp) :: hm_alpha3(2) = [2.10e-13_wp, -0.778_wp]   ! C, p
  !--- QEXC: beta, Eexc [eV], g_lower, a, b, c, d
  real(kind=wp) :: hm_q13(7), hm_q31a(7), hm_q31b(7)
  !--- CION3: A, p, T0 of C_3 = A T^p exp(-T0/T)
  real(kind=wp) :: hm_cion3(3) = 0.0_wp
  !--- PENNING: A, c, d1, d2, d3 (Garcia Munoz 2025)
  real(kind=wp) :: hm_pen(5)
  !--- sigma_3 VFKY96 wings [E_th,E_0,s_0,y_a,P,y_w,y_0,y_1] + bridge edges
  real(kind=wp) :: hm_sigA(8), hm_sigB(8), hm_E3, hm_E4
  logical :: hm_loaded = .false.
  !--- Phi_3 path estimator Sum(path_lum sigma_3(E)/(h nu)), rank local then
  !--- reduced, over the leaves.
  real(kind=wp), allocatable :: phi3_path(:)

  real(kind=wp), parameter :: kb_eV = kboltz_cgs / ev2erg    ! [eV/K]

contains

  !=========================================================================
  ! Public entry: read coefficients (once), turn the Phi_3 path estimator of
  ! the last transport pass into a rate, evaluate n_3 leaf by leaf, and write
  ! '<base>_heimeta'.
  !=========================================================================
  subroutine hei_metastable_run()
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, &
                              gas_ne, gas_Te
    implicit none
    real(kind=wp), allocatable :: phi3(:), n2s3(:)
    real(kind=wp) :: efloor, nH, yHe, ne, Te, Tc
    real(kind=wp) :: nHI, nHeI, nHeII, a3, q13, q31a, q31b, q31g, c3, Q31
    real(kind=wp) :: src, snk
    integer :: il, nleaf

    call hei_metastable_load()

    !--- band-floor note: Phi_3 threshold is 4.78 eV; if the band starts above
    !--- it the soft-UV part of the integral is missed (Phi_3 underestimated).
    efloor = par%eion_min
    if (par%add_fuv) efloor = par%efuv_min
    if (mpar%p_rank == 0 .and. efloor > 4.78_wp) &
       write(*,'(a,f7.3,a)') ' HEIMETA: NOTE: band floor ', efloor, &
       ' eV > 4.78 eV -> Phi_3 misses [4.78, floor); set add_fuv with '// &
       'efuv_min=4.78 for the full 2^3S photoionization integral.'

    nleaf = amr_grid%nleaf
    allocate(phi3(nleaf), n2s3(nleaf))

    !--- Phi_3(il) = Int 4 pi J_nu sigma_3 dnu / (h nu), from the reduced
    !--- path estimator of the last transport pass.
    call hei_phi3(phi3, nleaf)

    yHe = par%abund_He
    do il = 1, nleaf
       nH  = gas_nH(il);   ne = gas_ne(il);   Te = gas_Te(il)
       Tc  = min(max(Te, hm_Tmin), hm_Tmax)            ! clamp to the fit range
       nHI   = nH * gas_xHI(il)
       nHeI  = nH * yHe * gas_xHeI(il)
       nHeII = nH * yHe * gas_xHeII(il)
       a3    = hei_alpha_power(Tc, hm_alpha3)
       q13   = hei_qexc(Tc, hm_q13)
       q31a  = hei_qexc(Tc, hm_q31a)
       q31b  = hei_qexc(Tc, hm_q31b)
       q31g  = hei_qdeexc_ground(Tc)
       c3    = hm_cion3(1)*Tc**hm_cion3(2)*exp(-hm_cion3(3)/Tc)
       Q31   = hei_penning(Tc)
       src   = a3*ne*nHeII + q13*ne*nHeI
       snk   = hm_A31 + phi3(il) + ne*(q31a + q31b + q31g + c3) + nHI*Q31
       n2s3(il) = src / snk                            ! snk >= A31 > 0 always
    end do

    call hei_metastable_write(n2s3, phi3, nleaf)
    deallocate(phi3, n2s3)
  end subroutine hei_metastable_run

  !=========================================================================
  ! Phi_3 [s^-1] per leaf from the reduced path estimator: the factor turns
  ! Sum(path_lum coeff) into Int 4 pi J_nu coeff dnu, as for Gamma_i.
  !=========================================================================
  subroutine hei_phi3(phi3, nleaf)
    implicit none
    integer,       intent(in)  :: nleaf
    real(kind=wp), intent(out) :: phi3(nleaf)
    integer :: il

    if (.not. allocated(phi3_path)) &
       error stop 'Phi_3 path estimator is not allocated'
    do il = 1, nleaf
       phi3(il) = phi3_path(il) &
                  / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
    end do
  end subroutine hei_phi3

  !=========================================================================
  ! Phi_3 path estimator: allocated with the other transport estimators
  ! (ion_score_setup, also on re-refinement), zeroed at the start of each
  ! transport pass, scored for every ionizing- or FUV-band segment at the
  ! packet's own energy, and reduced with the rates.
  !=========================================================================
  subroutine hei_metastable_path_setup(nleaf)
    integer, intent(in) :: nleaf
    call hei_metastable_load()
    if (allocated(phi3_path)) deallocate(phi3_path)
    allocate(phi3_path(nleaf))
    phi3_path = 0.0_wp
  end subroutine hei_metastable_path_setup

  subroutine hei_metastable_path_reset()
    if (allocated(phi3_path)) phi3_path = 0.0_wp
  end subroutine hei_metastable_path_reset

  subroutine hei_metastable_path_score(il, E, path_lum)
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: E, path_lum
    if (E <= hm_sigA(1)) return                       ! below 4.78 eV
    phi3_path(il) = phi3_path(il) + path_lum*hei_sigma3(E)/(E*ev2erg)
  end subroutine hei_metastable_path_score

  subroutine hei_metastable_path_reduce()
    use mpi
    integer :: ierr
    if (allocated(phi3_path)) &
       call MPI_ALLREDUCE(MPI_IN_PLACE, phi3_path, size(phi3_path), &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine hei_metastable_path_reduce

  !=========================================================================
  ! alpha(T) = C (T/1e4)^p  [cm^3 s^-1].
  !=========================================================================
  real(kind=wp) function hei_alpha_power(T, cp) result(a)
    real(kind=wp), intent(in) :: T, cp(2)
    a = cp(1) * (T/1.0e4_wp)**cp(2)
  end function hei_alpha_power

  !=========================================================================
  ! q(T) = beta Upsilon(T) exp(-Eexc/kT) / (g_l sqrt T) [cm^3 s^-1],
  ! Upsilon = a exp(bT) + c exp(dT).  p(1:7) = beta, Eexc, g_l, a, b, c, d.
  !=========================================================================
  real(kind=wp) function hei_qexc(T, p) result(q)
    real(kind=wp), intent(in) :: T, p(7)
    real(kind=wp) :: ups, arg
    ups = p(4)*exp(p(5)*T) + p(6)*exp(p(7)*T)
    arg = -p(2)/(kb_eV*T)
    if (arg > -700.0_wp) then
       q = p(1)*ups*exp(arg)/(p(3)*sqrt(T))
    else
       q = 0.0_wp                                       ! Boltzmann underflow
    end if
  end function hei_qexc

  !=========================================================================
  ! q_31g(T) [cm^3 s^-1]: e-impact 2^3S -> 1^1S, the inverse of q_13 by
  ! detailed balance, q_31g = q_13 (g_1/g_3) exp(E_13/kT) with g_1 = 1,
  ! g_3 = 3.  Written with the Boltzmann factor of q_13 removed rather than
  ! multiplied back, so it does not underflow at low T.
  !=========================================================================
  real(kind=wp) function hei_qdeexc_ground(T) result(q)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: ups
    ups = hm_q13(4)*exp(hm_q13(5)*T) + hm_q13(6)*exp(hm_q13(7)*T)
    q   = hm_q13(1)*ups/(3.0_wp*sqrt(T))
  end function hei_qdeexc_ground

  !=========================================================================
  ! Q_31(T) [cm^3 s^-1]: Penning + associative ionization He(2^3S) + H0,
  ! Garcia Munoz (2025, A&A 698 A199, Fig. 5), the Maxwellian average of the
  ! Movre & Meyer (1997) cross sections:
  !   Q = A exp(c/T + d1 lnT + d2 (lnT)^2 + d3 (lnT)^3),
  ! hm_pen = A, c, d1, d2, d3.  <2% from 200 to 1e4 K; smooth and positive
  ! beyond, and T is clamped to TRANGE by the caller.
  !=========================================================================
  real(kind=wp) function hei_penning(T) result(Q)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lnT
    lnT = log(T)
    Q = hm_pen(1)*exp(hm_pen(2)/T + lnT*(hm_pen(3) + lnT*(hm_pen(4) + lnT*hm_pen(5))))
  end function hei_penning

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
    character(len=256) :: line
    character(len=16)  :: key
    integer :: unit, ios, ierr

    if (hm_loaded) return
    fname = trim(par%atomic_dir)//'/hei_metastable.txt'
    open(newunit=unit, file=trim(fname), status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(3a)') 'ERROR: cannot open ', &
          trim(fname), ' (par%atomic_dir)'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) key
       select case (trim(key))
       case ('TRANGE');    read(line,*) key, hm_Tmin, hm_Tmax
       case ('A31');       read(line,*) key, hm_A31
       case ('ALPHA3');    read(line,*) key, hm_alpha3
       case ('ALPHA1');    ! singlet recombination; not part of the n_3 balance
       case ('QEXC13');    read(line,*) key, hm_q13
       case ('QEXC31A');   read(line,*) key, hm_q31a
       case ('QEXC31B');   read(line,*) key, hm_q31b
       case ('CION3');     read(line,*) key, hm_cion3
       case ('PENNING');   read(line,*) key, hm_pen
       case ('SIGMA3A');   read(line,*) key, hm_sigA
       case ('SIGMA3B');   read(line,*) key, hm_sigB
       case ('SIGBRIDGE'); read(line,*) key, hm_E3, hm_E4
       case ('Q10833');    ! Stage 2 emission source; recorded, not used here
       end select
    end do
    close(unit)
    !--- the QEXC rows carry (beta, Eexc, g_l, ...) since 2026-09-24; a file in
    !--- the earlier (prefac, sconst, Eexc, div, ...) layout would be read
    !--- silently wrong, and a missing CION3 would drop the ionization sink.
    if (hm_q13(3) /= 1.0_wp .or. hm_q31a(3) /= 3.0_wp .or. &
        hm_q31b(3) /= 3.0_wp .or. hm_cion3(1) <= 0.0_wp) then
       if (mpar%p_rank == 0) write(*,'(3a)') 'ERROR: ', trim(fname), &
          ' is not in the current QEXC/CION3 layout; regenerate it with'// &
          ' tools/fitting/make_hei_metastable.py'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    hm_loaded = .true.
    if (mpar%p_rank == 0) write(*,'(2a)') &
       ' HEIMETA: 2^3S coefficients loaded from ', trim(fname)
  end subroutine hei_metastable_load

  !=========================================================================
  ! Write '<base>_heimeta' (rank 0): n_2s3, Phi_3, n_e, T_e + leaf geometry.
  !=========================================================================
  subroutine hei_metastable_write(n2s3, phi3, nleaf)
    use iofile_mod
    use utility,       only : get_base_name, fatal_error
    use gas_state_mod, only : gas_nH, gas_ne, gas_Te
    implicit none
    integer,       intent(in) :: nleaf
    real(kind=wp), intent(in) :: n2s3(nleaf), phi3(nleaf)
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
