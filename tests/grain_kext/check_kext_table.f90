program check_kext_table
!---------------------------------------------------------------------------
! Gate: the dust extinction MoCHII transports is the extinction of the grains
! MoCHII lets radiate.
!
! grain_model_mod builds one SEDust model and hands gas_opacity_mod its
! (lambda, albedo, <cos>, C_ext) curve through dust_extinction.  That routine
! no longer integrates over grain size at call time: it serves the PRECOMPUTED
! size-integrated curve of a data/kext_*.dat table -- the file the builder was
! given through kext_path -- interpolated onto the model's own wavelength grid.
! The size integral over the model as built is a separate routine,
! size_integrated_extinction, and it is what wrote those tables.  A transport
! host must not call it -- SEDust's own rt_example says so, and for MoCHII's
! reasons too: it re-does the triple sum over populations, sizes and
! wavelengths, it returns the number already on disk, and taking some
! wavelengths from it and others from the table would make the run's own optics
! inconsistent.  Writing a table, and checking one, are the uses it exists for,
! and this gate is the second.
!
! So the two must agree, and where they do not the physics says how much.
! Four quantities are compared, each the way the transport reads it: C_ext
! (optical depth) and C_abs (the heated fraction) relative to themselves,
! because they set rates; the albedo (the chance an interaction scatters) and
! albedo*<cos> (the forward displacement a scattering carries) absolutely,
! because they are probabilities the photon uses as they stand.
!
!   (a) lambda >= 0.0912 um.  The model grid here IS the table grid (both come
!       from the T-matrix Q table, and Zubko's grid is its own DustEM table),
!       so no interpolation happens and the served value is the table value.
!       Any departure beyond the 13 significant digits the tables are written
!       with means the table does not belong to this model -- a different size
!       distribution, a different Q table, a different ZDA config.  That is the
!       one failure this design can hide, because a mismatched table still
!       loads, still covers the band, and still looks like a plausible
!       extinction curve.  TOL_TABLE_NODE.
!
!   (b) lambda < 0.0912 um.  MoCHII asks SEDust to extend the grid down to
!       1.23984/eion_max um so the ionizing band is resolved, and the extended
!       nodes are NOT the table's nodes, so the log(lambda)-log(C)
!       interpolation enters.  That error is real but bounded and is not a
!       defect.  TOL_EUV_INTERP.
!
!   (c) The two integrals the transport actually uses.  s_abs(E) = C_abs/C_ext
!       (lambda_ref) is the band absorption MoCHII carries per unit reference
!       extinction, and C_ext(lambda_ref = 0.55 um) is the normalization every
!       leaf's grey rhokap is built on.  The (b) departures alternate in sign
!       along the grid, so they cancel here: these two must agree far more
!       tightly than any single wavelength does.  TOL_SABS_BAND, TOL_CREF.
!
! The cases are the grids MoCHII actually builds: astrodust extended to
! par%eion_max = 100 eV and to 150 eV, astrodust unextended (an emission-only
! run, par%ion_add_dust = .false.), DL07 extended to 100 eV, and Zubko, which
! needs no extension.  Everything else is at the par% defaults of define.f90.
!
! Build and run:  see tests/grain_kext/run_check.sh  (must run in SEDust/sed/,
! the directory all these paths are relative to)
!---------------------------------------------------------------------------
  use constants, only : wp
  use dust_lib,  only : dust_model_t, build_astrodust, build_dl07, build_zubko, &
                        dust_extinction, size_integrated_extinction, &
                        dust_nlam, dust_lambda
  implicit none

  !--- par% defaults of src/define.f90, as grain_model_mod passes them
  character(len=*), parameter :: QTABLE   = &
     '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
  character(len=*), parameter :: SIZEDIST = '../data/release/size_distribution.dat'
  character(len=*), parameter :: ZUBKO_CONFIG = '../data/zubko/ZDA_BARE_GR_S_Config.dat'
  character(len=*), parameter :: ZUBKO_DIR    = '../data/zubko/'
  !--- SEDust's standard table for each model: what a run gets when
  !--- par%sed_kext is blank.  Named here so the gate exercises the same files.
  character(len=*), parameter :: KEXT_ASTRODUST = '../data/kext_astrodust_MW_euv.dat'
  character(len=*), parameter :: KEXT_DL07      = '../data/kext_dl07_MW_euv.dat'
  character(len=*), parameter :: KEXT_ZUBKO     = '../data/kext_zubko_BARE_GR_S.dat'
  integer,       parameter :: NT      = 200
  real(kind=wp), parameter :: T_LO    = 2.7_wp, T_HI = 5.0e3_wp
  integer,       parameter :: SDINDEX = 7
  real(kind=wp), parameter :: U_ISRF  = 1.0_wp
  real(kind=wp), parameter :: LAM_REF = 0.55_wp        ! par%lambda_ref [um]

  !--- 0.0912 um = 13.598 eV: the short end of the T-matrix Q table, hence the
  !--- boundary between "the model grid is the table grid" and "MoCHII extended
  !--- the grid below the table's nodes".
  real(kind=wp), parameter :: LAM_QTAB = 0.0912_wp
  !--- transported band [eV] the band-averaged s_abs of check (c) runs over
  real(kind=wp), parameter :: E_BAND_LO = 13.598_wp, E_BAND_HI = 100.0_wp

  !--- (a) node-coincident: the served value is a table read, so the only
  !--- difference allowed is the 13-digit precision of the table file itself --
  !--- measured 4.7-4.9e-13 on the cross sections and 2.7-4.2e-13 on the albedo
  !--- pair, all five cases.  Three orders of margin on that still leaves eight
  !--- orders between here and the smallest table mismatch a changed size
  !--- distribution could produce.
  real(kind=wp), parameter :: TOL_TABLE_NODE     = 1.0e-10_wp
  real(kind=wp), parameter :: TOL_TABLE_NODE_DIM = 1.0e-10_wp
  !--- (b) EUV extension: log-log interpolation between table nodes.  Measured
  !--- worst single wavelength 1.2e-2 (astrodust on the 100 eV grid), 8.2e-3
  !--- (150 eV), 1.1e-3 (DL07), 4.8e-13 (Zubko, whose grid is the table's own),
  !--- and 3.3e-3 / 9.7e-3 / 2.2e-4 / 2.8e-13 on the albedo pair.  3e-2 is ~3x
  !--- the worst: loose enough that grid refinement or an eion_max change does
  !--- not trip it, tight enough that a wrong table (which is O(1) here as well
  !--- as at (a)) still does.  The albedo and <cos> are O(0.1-0.5) across this
  !--- band, so the same number is the same fractional statement about them.
  real(kind=wp), parameter :: TOL_EUV_INTERP     = 3.0e-2_wp
  real(kind=wp), parameter :: TOL_EUV_INTERP_DIM = 3.0e-2_wp
  !--- (c) the transported integrals.  Measured on the band-averaged s_abs:
  !--- 4.1e-6 (astrodust 100 eV), 2.0e-5 (150 eV), 3.9e-6 (DL07), 2.9e-14
  !--- (Zubko) -- two to three orders below the single wavelengths they are made
  !--- of, which is the cancellation this check exists to hold.  On C_ext at
  !--- 0.55 um: 2.8e-14 / 7.0e-15 / 1.4e-14, pure rounding, since 0.55 um sits
  !--- in the node-coincident region.  The tolerances sit 5x and five orders
  !--- above those, and both must stay far below TOL_EUV_INTERP: an integral
  !--- that drifted up to the single-wavelength bound would mean the departures
  !--- had stopped cancelling and started accumulating.
  real(kind=wp), parameter :: TOL_SABS_BAND = 1.0e-4_wp
  real(kind=wp), parameter :: TOL_CREF      = 1.0e-9_wp

  type(dust_model_t) :: m
  integer :: st
  logical :: ok

  ok = .true.

  write(*,'(a)') '=== table-served extinction vs the size integral that wrote the table ==='
  write(*,'(a)') 'worst single-wavelength departure on each side of 0.0912 um -- relative'
  write(*,'(a)') 'for the cross sections (C_ext, C_abs), absolute for the scattering pair'
  write(*,'(a)') '(albedo, albedo*<cos>) -- then the two band integrals the transport uses.'
  write(*,'(/,a24,4a12,2a12)') 'case', 'C >=.0912', 'a,ag>=.0912', &
     'C <.0912', 'a,ag<.0912', 'band s_abs', 'C_ext(ref)'

  !--- astrodust, EUV-extended to par%eion_max = 100 eV
  call build_astrodust(m, QTABLE, SIZEDIST, NT, T_LO, T_HI, status=st, &
                       lam_min=1.23984_wp/100.0_wp, kext_path=KEXT_ASTRODUST)
  call compare_served_with_integral('astrodust, EUV 100 eV', m, st, ok)

  !--- astrodust, EUV-extended to 150 eV (a harder spectrum's band)
  call build_astrodust(m, QTABLE, SIZEDIST, NT, T_LO, T_HI, status=st, &
                       lam_min=1.23984_wp/150.0_wp, kext_path=KEXT_ASTRODUST)
  call compare_served_with_integral('astrodust, EUV 150 eV', m, st, ok)

  !--- astrodust with no extension: par%ion_add_dust = .false., the native
  !--- T-matrix grid, entirely inside the table's nodes.
  call build_astrodust(m, QTABLE, SIZEDIST, NT, T_LO, T_HI, status=st, &
                       kext_path=KEXT_ASTRODUST)
  call compare_served_with_integral('astrodust, no extension', m, st, ok)

  !--- DL07, EUV-extended to 100 eV
  call build_dl07(m, QTABLE, SIZEDIST, SDINDEX, U_ISRF, NT, T_LO, T_HI, status=st, &
                  lam_min=1.23984_wp/100.0_wp, kext_path=KEXT_DL07)
  call compare_served_with_integral('dl07, EUV 100 eV', m, st, ok)

  !--- Zubko: its own DustEM tables already span the ionizing band, so the model
  !--- grid is the table grid everywhere and nothing is interpolated.
  call build_zubko(m, ZUBKO_CONFIG, ZUBKO_DIR, NT, T_LO, T_HI, status=st, &
                   kext_path=KEXT_ZUBKO)
  call compare_served_with_integral('zubko', m, st, ok)

  write(*,'(/,a,4(a,es8.1))') 'tolerances:', '  node ', TOL_TABLE_NODE, &
     ' / ', TOL_TABLE_NODE_DIM, '   EUV ', TOL_EUV_INTERP, ' / ', TOL_EUV_INTERP_DIM
  write(*,'(a,2(a,es8.1))') '           ', ' band s_abs ', TOL_SABS_BAND, &
     '   C_ext(ref) ', TOL_CREF
  write(*,'(/,a)') 'GATE (dust extinction table = the model''s own size integral): ' // &
     merge('PASS', 'FAIL', ok)
  if (.not. ok) error stop 1

contains

  !=========================================================================
  !--- One model: compare what dust_extinction serves with what the size
  !--- integral over the same model gives, split at 0.0912 um, plus the two
  !--- band integrals the transport is normalized on.
  !=========================================================================
  subroutine compare_served_with_integral(label, m, st_build, ok)
    character(len=*),   intent(in)    :: label
    type(dust_model_t), intent(in)    :: m
    integer,            intent(in)    :: st_build
    logical,            intent(inout) :: ok
    real(kind=wp), allocatable :: lam(:)
    real(kind=wp), allocatable :: Ce_t(:), Ca_t(:), Cs_t(:), g_t(:)
    real(kind=wp), allocatable :: Ce_i(:), Ca_i(:), Cs_i(:), g_i(:)
    real(kind=wp) :: dev_node, dev_euv, lam_node, lam_euv
    real(kind=wp) :: dim_node, dim_euv, lam_dnode, lam_deuv
    real(kind=wp) :: cref_t, cref_i, dev_cref
    real(kind=wp) :: sabs_t, sabs_i, dev_sabs, lam_lo, lam_hi
    real(kind=wp) :: alb_t, alb_i
    integer :: n, i, nb, st_t, st_i
    logical :: have_band

    if (st_build /= 0) then
       write(*,'(a24,a,i0,a)') label, '   BUILD FAILED (status=', st_build, ')'
       ok = .false.
       return
    end if

    n = dust_nlam(m)
    allocate(lam(n), Ce_t(n), Ca_t(n), Cs_t(n), g_t(n), &
                     Ce_i(n), Ca_i(n), Cs_i(n), g_i(n))
    lam = dust_lambda(m)
    call dust_extinction(m, Ce_t, Ca_t, Cs_t, gbar=g_t, status=st_t)
    call size_integrated_extinction(m, Ce_i, Ca_i, Cs_i, gbar=g_i, status=st_i)
    if (st_t /= 0 .or. st_i /= 0) then
       write(*,'(a24,a,i0,a,i0,a)') label, '   EXTINCTION FAILED (served=', st_t, &
          ', integral=', st_i, ')'
       ok = .false.
       deallocate(lam, Ce_t, Ca_t, Cs_t, g_t, Ce_i, Ca_i, Cs_i, g_i)
       return
    end if

    !--- (a) + (b): worst single-wavelength departure on each side of the
    !--- T-matrix short end, over all four transported quantities.
    !--- The cross sections C_ext (optical depth) and C_abs (the heated
    !--- fraction) are compared RELATIVE to themselves: they set rates, and a
    !--- fractional error in them is a fractional error in the transport.
    !--- The scattering pair enters the transport as probabilities, so it is
    !--- compared ABSOLUTELY, and as the transport weights it: the albedo is the
    !--- chance that an interaction scatters, and <cos> only ever acts through a
    !--- scattering event, so what the photon sees is albedo and albedo*<cos>.
    !--- That weighting is physics, not leniency, and it is what keeps the far
    !--- infrared honest: there C_sca falls to 1e-10 of C_ext, its own relative
    !--- accuracy collapses (Mie's small-x scattering series is
    !--- cancellation-limited, and the shipped DL07 table and the size integral
    !--- part company at the 4e-8 level there), yet the scattering it describes
    !--- cannot move a photon.  In the optical, where the albedo is 0.6-0.7,
    !--- the same measure is as tight as a direct comparison.
    dev_node = 0.0_wp;  lam_node = 0.0_wp;  dim_node = 0.0_wp;  lam_dnode = 0.0_wp
    dev_euv  = 0.0_wp;  lam_euv  = 0.0_wp;  dim_euv  = 0.0_wp;  lam_deuv  = 0.0_wp
    do i = 1, n
       alb_t = 0.0_wp;  if (Ce_t(i) > 0.0_wp) alb_t = Cs_t(i)/Ce_t(i)
       alb_i = 0.0_wp;  if (Ce_i(i) > 0.0_wp) alb_i = Cs_i(i)/Ce_i(i)
       if (lam(i) >= LAM_QTAB) then
          call worst_rel(lam(i), Ce_t(i), Ce_i(i), dev_node, lam_node)
          call worst_rel(lam(i), Ca_t(i), Ca_i(i), dev_node, lam_node)
          call worst_abs(lam(i), alb_t,         alb_i,         dim_node, lam_dnode)
          call worst_abs(lam(i), alb_t*g_t(i),  alb_i*g_i(i),  dim_node, lam_dnode)
       else
          call worst_rel(lam(i), Ce_t(i), Ce_i(i), dev_euv, lam_euv)
          call worst_rel(lam(i), Ca_t(i), Ca_i(i), dev_euv, lam_euv)
          call worst_abs(lam(i), alb_t,         alb_i,         dim_euv, lam_deuv)
          call worst_abs(lam(i), alb_t*g_t(i),  alb_i*g_i(i),  dim_euv, lam_deuv)
       end if
    end do

    !--- (c) the reference extinction gas_opacity_mod normalizes on, taken the
    !--- same way it takes it: log-log interpolation at par%lambda_ref.
    cref_t = interp_log(lam, Ce_t, n, LAM_REF)
    cref_i = interp_log(lam, Ce_i, n, LAM_REF)
    dev_cref = abs(cref_t/cref_i - 1.0_wp)

    !--- (c) band-averaged absorption per unit reference extinction.  This is
    !--- MoCHII's ion_dust_sabs = (1 - albedo) C_ext/C_ext(ref) = C_abs/C_ext(ref),
    !--- averaged over the model wavelengths inside the transported band, each
    !--- side normalized by its OWN C_ext(ref) as the code would.
    lam_hi = 1.23984_wp/E_BAND_LO;  lam_lo = 1.23984_wp/E_BAND_HI
    sabs_t = 0.0_wp;  sabs_i = 0.0_wp;  nb = 0
    do i = 1, n
       if (lam(i) >= lam_lo .and. lam(i) <= lam_hi) then
          sabs_t = sabs_t + Ca_t(i)/cref_t
          sabs_i = sabs_i + Ca_i(i)/cref_i
          nb = nb + 1
       end if
    end do
    !--- an unextended grid touches the band only at its 0.0912 um endpoint,
    !--- so there is no band average to form there.
    have_band = (nb >= 2 .and. sabs_i > 0.0_wp)
    dev_sabs = 0.0_wp
    if (have_band) dev_sabs = abs(sabs_t/sabs_i - 1.0_wp)

    if (have_band) then
       write(*,'(a24,6es12.2)') label, dev_node, dim_node, dev_euv, dim_euv, &
          dev_sabs, dev_cref
    else
       write(*,'(a24,4es12.2,a12,es12.2)') label, dev_node, dim_node, dev_euv, dim_euv, &
          '  n/a', dev_cref
    end if

    if (dev_node > TOL_TABLE_NODE) then
       write(*,'(2a,es10.3,a,es11.4,a)') '  FAIL: ', &
          'the served cross sections are not the table at lam >= 0.0912 um: ', &
          dev_node, ' at ', lam_node, ' um (is the table this model''s?)'
       ok = .false.
    end if
    if (dim_node > TOL_TABLE_NODE_DIM) then
       write(*,'(2a,es10.3,a,es11.4,a)') '  FAIL: ', &
          'the served scattering pair is not the table at lam >= 0.0912 um: ', &
          dim_node, ' at ', lam_dnode, ' um (is the table this model''s?)'
       ok = .false.
    end if
    if (dev_euv > TOL_EUV_INTERP) then
       write(*,'(2a,es10.3,a,es11.4,a)') '  FAIL: ', &
          'EUV interpolation departs by ', dev_euv, ' at ', lam_euv, ' um'
       ok = .false.
    end if
    if (dim_euv > TOL_EUV_INTERP_DIM) then
       write(*,'(2a,es10.3,a,es11.4,a)') '  FAIL: ', &
          'EUV scattering-pair interpolation departs by ', dim_euv, ' at ', &
          lam_deuv, ' um'
       ok = .false.
    end if
    if (have_band .and. dev_sabs > TOL_SABS_BAND) then
       write(*,'(a,es10.3)') &
          '  FAIL: band-averaged s_abs (13.598-100 eV) departs by ', dev_sabs
       ok = .false.
    end if
    if (dev_cref > TOL_CREF) then
       write(*,'(a,es10.3)') '  FAIL: C_ext at lambda_ref = 0.55 um departs by ', dev_cref
       ok = .false.
    end if

    deallocate(lam, Ce_t, Ca_t, Cs_t, g_t, Ce_i, Ca_i, Cs_i, g_i)
  end subroutine compare_served_with_integral

  !--- Keep the worst RELATIVE departure and where it is (the cross sections).
  subroutine worst_rel(lam_i, a, b, dev, lam_worst)
    real(kind=wp), intent(in)    :: lam_i, a, b
    real(kind=wp), intent(inout) :: dev, lam_worst
    real(kind=wp) :: d, scale
    scale = max(abs(a), abs(b))
    if (scale <= 0.0_wp) return
    d = abs(a - b)/scale
    if (d > dev) then
       dev = d;  lam_worst = lam_i
    end if
  end subroutine worst_rel

  !--- Keep the worst ABSOLUTE departure and where it is (albedo and <cos>,
  !--- which the transport uses as probabilities, not as rates).
  subroutine worst_abs(lam_i, a, b, dev, lam_worst)
    real(kind=wp), intent(in)    :: lam_i, a, b
    real(kind=wp), intent(inout) :: dev, lam_worst
    real(kind=wp) :: d
    d = abs(a - b)
    if (d > dev) then
       dev = d;  lam_worst = lam_i
    end if
  end subroutine worst_abs

  !--- log-log interpolation on an ascending grid, the same expression
  !--- gas_opacity_mod uses to evaluate C_ext at par%lambda_ref.
  real(kind=wp) function interp_log(x, y, n, x0) result(y0)
    integer,       intent(in) :: n
    real(kind=wp), intent(in) :: x(n), y(n), x0
    integer :: i
    real(kind=wp) :: w
    if (x0 <= x(1)) then
       y0 = y(1);  return
    else if (x0 >= x(n)) then
       y0 = y(n);  return
    end if
    do i = 1, n-1
       if (x0 < x(i+1)) exit
    end do
    w  = log(x0/x(i))/log(x(i+1)/x(i))
    y0 = exp(log(y(i))*(1.0_wp-w) + log(y(i+1))*w)
  end function interp_log

end program check_kext_table
