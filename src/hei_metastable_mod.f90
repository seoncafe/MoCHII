module hei_metastable_mod
!---------------------------------------------------------------------------
! MoCHII: He I 2^3S (triplet metastable) diagnostic.
!
! On the converged state this evaluates the local 2^3S population n_3 from the
! Oklopcic & Hirata (2018, ApJ 855 L11 = OH18) balance, which is LINEAR in n_3
! (a closed-form quotient per leaf, no iteration):
!
!   n_3 = [ alpha_3 n_e n_HeII + q_13 n_e n_HeI ]
!         / [ A_31 + Phi_3 + n_e (q_31a + q_31b) + n_HI Q_31 ]
!
!   alpha_3     He II -> 2^3S recombination (triplet channel)
!   q_13        e-impact 1^1S -> 2^3S
!   A_31        2^3S -> 1^1S radiative (1.272e-4 s^-1)
!   Phi_3       2^3S photoionization rate (threshold 4.78 eV; full-band tally)
!   q_31a,q_31b e-impact 2^3S -> 2^1S and 2^3S -> 2^1P
!   Q_31        Penning + associative ionization He(2^3S) + H0
!
! Coefficients are read from data/atomic/hei_metastable.txt (closed-form fits
! transcribed from the EXHALE triplet implementation; the Fortran here only
! evaluates the fits).  Phi_3 uses the reduced ionizing-band tally jt_ion, the
! same integral pattern as gas_rates_mod's Gamma_i.
!
! Runs once at output time (par%hei_metastable): a leaf loop that reads the
! converged gas state, no feedback into the ionization / thermal balance (n_3
! is a ~1e-7 trace of He, its energy budget negligible).  Output section file
! '<base>_heimeta' carries n_2s3, Phi_3, n_H, n_e, T_e + leaf geometry.
!
! CAVEATS (see docs/HEI_10833_PLAN.md sec. 8):
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
  use jtally_mod,   only : jt_ion
  use ion_band_mod, only : ion_e, nnu_band
  use photo_xsec,   only : sigma_vfky96
  implicit none
  private

  public :: hei_metastable_run

  !--- coefficients read from data/atomic/hei_metastable.txt.
  real(kind=wp) :: hm_Tmin = 1.0e2_wp, hm_Tmax = 3.0e4_wp
  real(kind=wp) :: hm_A31  = 1.272e-4_wp
  real(kind=wp) :: hm_alpha3(2) = [2.10e-13_wp, -0.778_wp]   ! C, p
  real(kind=wp) :: hm_alpha1(2) = [1.54e-13_wp, -0.486_wp]
  !--- QEXC: prefac, sconst, Eexc, div, a, b, c, d
  real(kind=wp) :: hm_q13(8), hm_q31a(8), hm_q31b(8)
  !--- PENNING: Tsplit, A_lo, p_lo, A_hi, p_hi
  real(kind=wp) :: hm_pen(5)
  !--- sigma_3 VFKY96 wings [E_th,E_0,s_0,y_a,P,y_w,y_0,y_1] + bridge edges
  real(kind=wp) :: hm_sigA(8), hm_sigB(8), hm_E3, hm_E4
  logical :: hm_loaded = .false.

  real(kind=wp), parameter :: kb_eV = kboltz_cgs / ev2erg    ! [eV/K]

contains

  !=========================================================================
  ! Public entry: read coefficients (once), integrate Phi_3 from the band
  ! tally, evaluate n_3 leaf by leaf, and write '<base>_heimeta'.
  !=========================================================================
  subroutine hei_metastable_run()
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, &
                              gas_ne, gas_Te
    implicit none
    real(kind=wp), allocatable :: phi3(:), n2s3(:)
    real(kind=wp) :: efloor, nH, yHe, ne, Te, Tc
    real(kind=wp) :: nHI, nHeI, nHeII, a3, q13, q31a, q31b, Q31, src, snk
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

    !--- Phi_3(il) = Sum_nu 4 pi J_nu sigma_3(E) dnu / (h nu), the same
    !--- reduced-tally integral pattern as gas_rates_mod's Gamma_i.  jt_ion is
    !--- already ALLREDUCEd before output; sub-threshold bins give sigma_3 = 0.
    call hei_phi3(phi3, nleaf)

    yHe = par%He_abund
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
       Q31   = hei_penning(Tc)
       src   = a3*ne*nHeII + q13*ne*nHeI
       snk   = hm_A31 + phi3(il) + ne*(q31a + q31b) + nHI*Q31
       n2s3(il) = src / snk                            ! snk >= A31 > 0 always
    end do

    call hei_metastable_write(n2s3, phi3, nleaf)
    deallocate(phi3, n2s3)
  end subroutine hei_metastable_run

  !=========================================================================
  ! Phi_3 [s^-1] per leaf from the reduced ionizing-band tally.
  !=========================================================================
  subroutine hei_phi3(phi3, nleaf)
    implicit none
    integer,       intent(in)  :: nleaf
    real(kind=wp), intent(out) :: phi3(nleaf)
    real(kind=wp) :: sig3(nnu_band), vol, fac, fJ, hnu
    integer :: il, inu

    !--- sigma_3 depends only on the (fixed) band -> evaluate once per bin.
    do inu = 1, nnu_band
       sig3(inu) = hei_sigma3(ion_e(inu))
    end do
    phi3 = 0.0_wp
    do il = 1, nleaf
       vol = (2.0_wp*leaf_half(il))**3
       fac = 1.0_wp / (vol * par%distance2cm**2)       ! -> 4 pi J_nu dnu per bin
       do inu = 1, nnu_band
          fJ = jt_ion(inu, il) * fac                   ! 4 pi J dnu [erg/s/cm^2]
          if (fJ <= 0.0_wp) cycle
          hnu = ion_e(inu) * ev2erg
          phi3(il) = phi3(il) + fJ * sig3(inu) / hnu
       end do
    end do
  end subroutine hei_phi3

  !=========================================================================
  ! alpha(T) = C (T/1e4)^p  [cm^3 s^-1].
  !=========================================================================
  real(kind=wp) function hei_alpha_power(T, cp) result(a)
    real(kind=wp), intent(in) :: T, cp(2)
    a = cp(1) * (T/1.0e4_wp)**cp(2)
  end function hei_alpha_power

  !=========================================================================
  ! q(T) = prefac sqrt(sconst/kT) exp(-Eexc/kT) (a exp(bT) + c exp(dT)) / div
  ! [cm^3 s^-1].  p(1:8) = prefac, sconst, Eexc, div, a, b, c, d.
  !=========================================================================
  real(kind=wp) function hei_qexc(T, p) result(q)
    real(kind=wp), intent(in) :: T, p(8)
    real(kind=wp) :: kT, ups, arg
    kT  = kb_eV * T
    ups = p(5)*exp(p(6)*T) + p(7)*exp(p(8)*T)
    arg = -p(3)/kT
    if (arg > -700.0_wp) then
       q = p(1)*sqrt(p(2)/kT)*exp(arg)*ups/p(4)
    else
       q = 0.0_wp                                       ! Boltzmann underflow
    end if
  end function hei_qexc

  !=========================================================================
  ! Q_31(T) [cm^3 s^-1]: Q = A_lo (300/T)^p_lo for T <= Tsplit, else
  ! A_hi (300/T)^p_hi.  hm_pen = Tsplit, A_lo, p_lo, A_hi, p_hi.
  !=========================================================================
  real(kind=wp) function hei_penning(T) result(Q)
    real(kind=wp), intent(in) :: T
    if (T <= hm_pen(1)) then
       Q = hm_pen(2) * (300.0_wp/T)**hm_pen(3)
    else
       Q = hm_pen(4) * (300.0_wp/T)**hm_pen(5)
    end if
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
    character(len=256) :: fname, line
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
       case ('ALPHA1');    read(line,*) key, hm_alpha1
       case ('QEXC13');    read(line,*) key, hm_q13
       case ('QEXC31A');   read(line,*) key, hm_q31a
       case ('QEXC31B');   read(line,*) key, hm_q31b
       case ('PENNING');   read(line,*) key, hm_pen
       case ('SIGMA3A');   read(line,*) key, hm_sigA
       case ('SIGMA3B');   read(line,*) key, hm_sigB
       case ('SIGBRIDGE'); read(line,*) key, hm_E3, hm_E4
       case ('Q10833');    ! Stage 2 emission source; recorded, not used here
       end select
    end do
    close(unit)
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
    character(len=192) :: filename
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
