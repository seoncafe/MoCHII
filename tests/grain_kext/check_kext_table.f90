program check_kext_table
!---------------------------------------------------------------------------
! Gate: the dust extinction MoCHII transports is the extinction of the grains
! MoCHII lets radiate.
!
! grain_model_mod builds one SEDust model -- one build_dust call, one data
! directory, the whole stored wavelength axis (include_euv) and no lam_min --
! and hands gas_opacity_mod its (lambda, albedo, <cos>, C_ext) curve through
! dust_extinction.  That routine does not integrate over grain size at call
! time: it serves the PRECOMPUTED size-integrated curve the builder loaded
! (/kext of the model's own product, or the file par%sed_kext names),
! interpolated onto the model's own wavelength grid.  The size integral over
! the model as built is a separate routine, size_integrated_extinction, and it
! is what wrote those curves.  A transport host must not call it -- SEDust's
! own rt_example says so, and for MoCHII's reasons too: it re-does the triple
! sum over populations, sizes and wavelengths, it returns the number already on
! disk, and taking some wavelengths from it and others from the stored curve
! would make the run's own optics inconsistent.  Writing a curve, and checking
! one, are the uses it exists for, and this gate is the second.
!
! So the two must agree, and where they do not the physics says how much.
! Four quantities are compared, each the way the transport reads it: C_ext
! (optical depth) and C_abs (the heated fraction) relative to themselves,
! because they set rates; the albedo (the chance an interaction scatters) and
! albedo*<cos> (the forward displacement a scattering carries) absolutely,
! because they are probabilities the photon uses as they stand.
!
!   (a) lambda >= 0.0912 um.  The model grid here IS the curve's grid -- both
!       are read from the same product -- so no interpolation happens and the
!       served value is the stored value.  Any departure beyond what separates
!       the two halves of that product (see the tolerances below) means the
!       curve does not belong to this model: a different size distribution,
!       different optics, a different ZDA config.  That is the one failure this
!       design can hide, because a mismatched curve still loads, still covers
!       the band, and still looks like a plausible extinction curve.
!
!   (b) lambda < 0.0912 um, the ionizing band itself.  Nothing changes there:
!       include_euv takes the whole axis, the stored curve spans that same
!       whole axis, and no lam_min extends the grid past it.  The check is
!       therefore the same statement as (a), kept separate because the two
!       halves sit in different numerical regimes (C_sca and C_abs are orders
!       apart across the seam).  Should a case with a genuine extension be
!       added (it needs par%eion_max above 12398.4 eV, and a library built
!       WITH_TMATRIX, or astrodust would refuse), the log(lambda)-log(C)
!       interpolation returns here and the tolerance has to go back to the
!       ~3e-2 the extended nodes used to cost.
!
!   (c) The two integrals the transport actually uses.  s_abs(E) = C_abs/C_ext
!       (lambda_ref) is the band absorption MoCHII carries per unit reference
!       extinction, and C_ext(lambda_ref = 0.55 um) is the normalization every
!       leaf's grey rhokap is built on.  TOL_SABS_BAND, and TOL_KEXT_NODE for
!       C_ext(lambda_ref).
!
! The cases are the three models MoCHII can build, each with the arguments
! grain_model_mod passes -- there is no lam_min case to vary any more, and
! include_euv is fixed .true. there, so one case per model is the whole set.
! Measured grids: astrodust NLAM = 1762 over 1.0000e-04 to 3.9811e+04 um, DL07
! NLAM = 1823 over 6.2054e-05 to 3.9811e+04 um (its OWN axis now, no longer the
! astrodust one it used to borrow), Zubko NLAM = 1201 over 1.0000e-03 to
! 1.0000e+04 um.  Everything else is at the par% defaults of define.f90.
!
! Build and run:  see tests/grain_kext/run_check.sh, which passes the data
! directory as argv(1).  Runs from any working directory.
!---------------------------------------------------------------------------
  use constants, only : wp
  use dust_lib,  only : dust_model_t, build_dust, &
                        dust_extinction, size_integrated_extinction, &
                        dust_nlam, dust_lambda
  implicit none

  !--- par% defaults of src/define.f90, as grain_model_mod passes them.
  !--- par%sed_data_dir comes in on the command line, the way a run gets it
  !--- from read_input: build_dust resolves every file the model is made of
  !--- inside it, dielectric functions included, so the gate needs no
  !--- particular working directory either.  Absent, the SEDust sed/ default.
  character(len=512) :: DATA_DIR = '../data'               ! par%sed_data_dir
  integer,       parameter :: NT      = 200
  real(kind=wp), parameter :: T_LO    = 2.7_wp, T_HI = 5.0e3_wp
  integer,       parameter :: SDINDEX = 7
  real(kind=wp), parameter :: U_ISRF  = 1.0_wp
  real(kind=wp), parameter :: LAM_REF = 0.55_wp        ! par%lambda_ref [um]

  !--- 0.0912 um = 13.598 eV: the Lyman limit, where SEDust cuts the stored
  !--- axis when include_euv is not asked for, hence the boundary the two
  !--- numerical regimes are reported on either side of.
  real(kind=wp), parameter :: LAM_QTAB = 0.0912_wp
  !--- transported band [eV] the band-averaged s_abs of check (c) runs over
  real(kind=wp), parameter :: E_BAND_LO = 13.598_wp, E_BAND_HI = 100.0_wp

  !--- What separates the served curve from the integral, model by model, is how
  !--- that model's product was written -- not how the gate reads it.  Every
  !--- model wavelength is a node of the stored curve, so nothing is
  !--- interpolated on either side of 0.0912 um; what is left is whether /kext
  !--- and /qtable of one sedust_<model>.h5 were computed from the SAME numbers.
  !---
  !--- They are, for all three, and one tolerance therefore covers them.  It was
  !--- two: DL07's and Zubko's /kext used to be written from the text q_*.dat
  !--- optics (7 significant digits) while /qtable of the same file carried full
  !--- double precision, which put those two models at the 7th digit and needed
  !--- 1e-6 where astrodust sat at 1e-9.  Both halves of every product now come
  !--- from one set of numbers.  Measured here (worst single wavelength, cross
  !--- sections / albedo pair, lam >= 0.0912 um then the ionizing band):
  !---   astrodust  0.0 / 0.0            3.3e-16 / 1.1e-16
  !---   dl07       2.7e-15 / 6.7e-16    4.2e-16 / 3.3e-16
  !---   zubko      5.8e-16 / 2.2e-16    6.3e-16 / 2.2e-16
  !--- and on C_ext(0.55 um): 0.0 (astrodust), 0.0 (DL07), 7.2e-15 (Zubko).
  !--- That is rounding on a double, so 1e-12 leaves three orders of headroom
  !--- and still catches a curve belonging to other grains, which is O(1).
  real(kind=wp), parameter :: TOL_KEXT_NODE = 1.0e-12_wp
  !--- (c) the band-averaged s_abs.  Measured 0.0 (astrodust), 0.0 (DL07),
  !--- 6.9e-15 (Zubko) -- an integral over 137 wavelengths, so the single-node
  !--- departures above average down rather than add.  What it has to catch is
  !--- an integral that has stopped agreeing, so it stays above the node
  !--- tolerance, by four orders; 1e-4 was eleven orders above the measurement
  !--- and no longer watched anything.
  real(kind=wp), parameter :: TOL_SABS_BAND = 1.0e-8_wp

  type(dust_model_t) :: m
  integer :: st
  logical :: ok

  ok = .true.
  if (command_argument_count() >= 1) call get_command_argument(1, DATA_DIR)

  write(*,'(2a)') 'par%sed_data_dir = ', trim(DATA_DIR)
  write(*,'(a)') '=== served extinction vs the size integral that wrote it ==='
  write(*,'(a)') 'worst single-wavelength departure on each side of 0.0912 um -- relative'
  write(*,'(a)') 'for the cross sections (C_ext, C_abs), absolute for the scattering pair'
  write(*,'(a)') '(albedo, albedo*<cos>) -- then the two band integrals the transport uses.'
  write(*,'(/,a22,a6,2a11,4a11,2a11)') 'case', 'NLAM', 'lam(1)', 'lam(N)', &
     'C >=.0912', 'a,ag>=.0912', 'C <.0912', 'a,ag<.0912', 'band s_abs', 'C_ext(ref)'

  !--- The three models grain_model_mod can build, with its arguments: the whole
  !--- stored axis (include_euv), no lam_min, and the model's own /kext (a blank
  !--- par%sed_kext, so kext_path is left out here too).
  call build_dust(m, 'astrodust', trim(DATA_DIR), NT, T_LO, T_HI, include_euv=.true., &
                  status=st, sd_index=SDINDEX, u_isrf=U_ISRF)
  call compare_served_with_integral('astrodust', m, st, TOL_KEXT_NODE, ok)

  call build_dust(m, 'dl07', trim(DATA_DIR), NT, T_LO, T_HI, include_euv=.true., &
                  status=st, sd_index=SDINDEX, u_isrf=U_ISRF)
  call compare_served_with_integral('dl07', m, st, TOL_KEXT_NODE, ok)

  call build_dust(m, 'zubko', trim(DATA_DIR), NT, T_LO, T_HI, include_euv=.true., &
                  status=st, sd_index=SDINDEX, u_isrf=U_ISRF)
  call compare_served_with_integral('zubko', m, st, TOL_KEXT_NODE, ok)

  write(*,'(/,a,a,es8.1)') 'tolerances:', &
     '  single wavelength (/kext vs its own size integral) ', TOL_KEXT_NODE
  write(*,'(a,es8.1)') '            band s_abs ', TOL_SABS_BAND
  write(*,'(/,a)') 'GATE (dust extinction curve = the model''s own size integral): ' // &
     merge('PASS', 'FAIL', ok)
  if (.not. ok) error stop 1

contains

  !=========================================================================
  !--- One model: compare what dust_extinction serves with what the size
  !--- integral over the same model gives, split at 0.0912 um, plus the two
  !--- band integrals the transport is normalized on.
  !=========================================================================
  subroutine compare_served_with_integral(label, m, st_build, tol, ok)
    character(len=*),   intent(in)    :: label
    type(dust_model_t), intent(in)    :: m
    integer,            intent(in)    :: st_build
    !--- How far apart this model's /kext and /qtable were written; see where
    !--- the two constants are declared for which model gets which and why.
    !--- The cross sections are held to it relatively, the albedo pair
    !--- absolutely -- both O(0.1-1) quantities, so it is the same statement.
    real(kind=wp),      intent(in)    :: tol
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
       write(*,'(a22,a,i0,a)') label, '   BUILD FAILED (status=', st_build, ')'
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
       write(*,'(a22,a,i0,a,i0,a)') label, '   EXTINCTION FAILED (served=', st_t, &
          ', integral=', st_i, ')'
       ok = .false.
       deallocate(lam, Ce_t, Ca_t, Cs_t, g_t, Ce_i, Ca_i, Cs_i, g_i)
       return
    end if

    !--- (a) + (b): worst single-wavelength departure on each side of the Lyman
    !--- limit, over all four transported quantities.
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
    !--- cancellation-limited, and the shipped DL07 curve and the size integral
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
    !--- a grid cut at the Lyman limit touches the band only at its 0.0912 um
    !--- endpoint, so there is no band average to form there.
    have_band = (nb >= 2 .and. sabs_i > 0.0_wp)
    dev_sabs = 0.0_wp
    if (have_band) dev_sabs = abs(sabs_t/sabs_i - 1.0_wp)

    if (have_band) then
       write(*,'(a22,i6,2es11.4,6es11.2)') label, n, lam(1), lam(n), &
          dev_node, dim_node, dev_euv, dim_euv, dev_sabs, dev_cref
    else
       write(*,'(a22,i6,2es11.4,4es11.2,a11,es11.2)') label, n, lam(1), lam(n), &
          dev_node, dim_node, dev_euv, dim_euv, '  n/a', dev_cref
    end if

    if (dev_node > tol) then
       write(*,'(2a,es10.3,a,es11.4,a)') '  FAIL: ', &
          'the served cross sections are not the stored curve at lam >= 0.0912 um: ', &
          dev_node, ' at ', lam_node, ' um (is the curve this model''s?)'
       ok = .false.
    end if
    if (dim_node > tol) then
       write(*,'(2a,es10.3,a,es11.4,a)') '  FAIL: ', &
          'the served scattering pair is not the stored curve at lam >= 0.0912 um: ', &
          dim_node, ' at ', lam_dnode, ' um (is the curve this model''s?)'
       ok = .false.
    end if
    if (dev_euv > tol) then
       write(*,'(2a,es10.3,a,es11.4,a)') '  FAIL: ', &
          'the ionizing band departs by ', dev_euv, ' at ', lam_euv, ' um'
       ok = .false.
    end if
    if (dim_euv > tol) then
       write(*,'(2a,es10.3,a,es11.4,a)') '  FAIL: ', &
          'the ionizing-band scattering pair departs by ', dim_euv, ' at ', &
          lam_deuv, ' um'
       ok = .false.
    end if
    if (have_band .and. dev_sabs > TOL_SABS_BAND) then
       write(*,'(a,es10.3)') &
          '  FAIL: band-averaged s_abs (13.598-100 eV) departs by ', dev_sabs
       ok = .false.
    end if
    if (dev_cref > tol) then
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
