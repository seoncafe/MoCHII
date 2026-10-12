module recomb_continuum_mod
!---------------------------------------------------------------------------
! MoCHII: free-bound + free-free continuum coefficients (port of the MOCASSIN
! emission_mod::fb_ff tables and expressions).
!
! Per (n_e n_ion), in units of 1e-40 erg cm^3 s^-1 Hz^-1:
!  - BELOW each species' first ionization edge: the tabulated coefficients
!    data/gamma{HI,HeI,HeII}.dat (21 T x paired-edge nu points; log-log T
!    interpolation, linear-in-nu interpolation of log gamma between edge
!    pairs, as MOCASSIN fb_ff), plus the ranges the tables leave out (the
!    gaps around each edge, the stretch up to the ground edge, and below the
!    first node), which returned 0 until 2026-09-26; see gamma_tab.
!    VERIFIED (2026-07-12, tools check vs a hydrogenic Kramers Milne sum):
!    the tables are free-bound PLUS free-free combined — do not add ff again
!    inside the table range.  The Ercolano & Storey (2006) machine-readable
!    swap remains a drop-in at this interface.
!  - ABOVE the table range: Milne-relation free-bound from the ground
!    photoionization cross sections (photo_xsec: exact hydrogenic for H I
!    and He II, VFKY96 for He I) + free-free with the Gaunt factor of
!    free_free_mod (Z=1 for H II/He II, Z=2 for He III).
!
! These coefficients have three consumers that must not drift apart: the
! written nebular continuum (nebcont_mod), the He II n>=2 free-bound
! component of the diffuse ionizing field (heii_cascade_mod), and the
! free-free packets of the diffuse field (diffuse_mod channel 7, through the
! band tables freefree_band_setup builds from gamma_ff).  They live here so
! that all read one coefficient.
!---------------------------------------------------------------------------
  use define
  use free_free_mod, only : free_free_gaunt, FF_SPECTRAL
  implicit none
  private

  public :: recomb_continuum_setup, recomb_continuum_is_ready
  public :: gamma_tab, gamma_milne, gamma_ff, gamma_tab_nodes
  public :: RYD_ERG, HCRYD_K, MILNE_C, NU_HI, NU_HEI, NU_HEII
  public :: freefree_band_setup, freefree_band_coefficient
  public :: freefree_band_sample_energy, freefree_band_is_ready
  public :: FF_ZMAX

  !--- tabulated continuum coefficients (log10 of 1e-40 erg cm^3/s/Hz)
  type gamma_table
     integer :: nT = 0, nnu = 0
     real(kind=wp), allocatable :: tk(:), nu(:), lg(:,:)   ! (nnu, nT)
  end type gamma_table
  type(gamma_table) :: gtab(3)          ! 1=HI, 2=HeI, 3=HeII
  logical, save :: ready = .false.

  !--- Free-free emission inside the transported band (diffuse_mod channel
  !--- 7), per net ion charge Z = 1..FF_ZMAX: energy-luminosity CDFs in
  !--- x = (E - E_lo)/kT and the band coefficients, on a log-T grid.  The
  !--- layout is that of the Milne ground-continuum tables
  !--- (milne_recomb_spectrum_mod): NX nodes clustered towards x = 0, cut
  !--- at x = 60.
  integer,       parameter :: FF_ZMAX = MAX_ION_STAGES - 1
  integer,       parameter :: FF_NT = 192, FF_NX = 512
  real(kind=wp), parameter :: FF_TLO = 1.0e2_wp, FF_THI = 1.0e6_wp
  real(kind=wp), parameter :: FF_XMAX = 60.0_wp
  real(kind=wp), save :: ff_lgT(FF_NT), ff_xg(FF_NX)
  real(kind=wp), save, allocatable :: ff_cdf(:,:,:)      ! (FF_NX, FF_NT, FF_ZMAX)
  real(kind=wp), save :: ff_lncoef(FF_NT, FF_ZMAX)       ! ln G_Z(T), see below
  real(kind=wp), save :: ff_elo = 0.0_wp, ff_ehi = 0.0_wp ! band edges [eV]
  logical,       save :: ff_ready = .false.

  real(kind=wp), parameter :: RYD_ERG   = 2.1798723611035e-11_wp
  !--- Ryd/k [K] for the frequencies of these tables, which are in units of
  !--- the infinite-mass Rydberg (RYD_ERG; the H edge sits at 0.999466 of
  !--- it, and MILNE_C carries R_inf^2).  It was MOCASSIN's 157807 K, close
  !--- to the H ionization energy over k (157802 K) rather than R_inf/k
  !--- (157887.5 K), so exp(-h nu/kT) was too large by exp(5.1e-4 h nu/kT):
  !--- the free-free coefficient +0.8% at 1 Ryd and 1e4 K, +3.3% at 4 Ryd.
  real(kind=wp), parameter :: HCRYD_K   = RYD_ERG/kboltz_cgs
  real(kind=wp), parameter :: MILNE_C   = 4.9874105e-6_wp   ! (h^2/2 pi me k)^1.5 Ryd^2 [cm K^1.5]
  real(kind=wp), parameter :: NU_HI     = 0.999466_wp       ! edges [Ryd]
  real(kind=wp), parameter :: NU_HEI    = 1.80804_wp
  real(kind=wp), parameter :: NU_HEII   = 4.001067_wp

contains

  !=========================================================================
  logical function recomb_continuum_is_ready() result(ok)
    ok = ready
  end function recomb_continuum_is_ready

  !=========================================================================
  ! Read the three tables once.  Idempotent, so any consumer may call it.
  !=========================================================================
  subroutine recomb_continuum_setup()
    use utility, only : fatal_error
    implicit none
    character(len=16) :: fn(3)
    integer :: k, i, unit, ios

    if (ready) return
    fn = [character(len=16) :: 'gammaHI.dat', 'gammaHeI.dat', 'gammaHeII.dat']
    do k = 1, 3
       open(newunit=unit, file=trim(par%data_dir)//'/'//trim(fn(k)), &
            status='old', iostat=ios)
       if (ios /= 0) call fatal_error('cannot open ' // &
          trim(par%data_dir)//'/'//trim(fn(k)) // ' (par%data_dir)')
       read(unit,*) gtab(k)%nT, gtab(k)%nnu
       allocate(gtab(k)%tk(gtab(k)%nT), gtab(k)%nu(gtab(k)%nnu), &
                gtab(k)%lg(gtab(k)%nnu, gtab(k)%nT))
       read(unit,*) gtab(k)%tk
       do i = 1, gtab(k)%nnu
          read(unit,*) gtab(k)%nu(i), gtab(k)%lg(i,:)
       end do
       close(unit)
       gtab(k)%lg = log10(max(gtab(k)%lg, 1.0e-30_wp))
    end do
    ready = .true.
  end subroutine recomb_continuum_setup

  !=========================================================================
  ! fb+ff coefficient of recombining species k (1 = H+, 2 = He+, 3 = He2+)
  ! at (nu, T) BELOW its ground edge (0 at and above it, where the caller
  ! takes the Milne ground continuum plus free-free).  Inside the table:
  ! linear log-gamma interpolation within an edge pair, log-log linear in T
  ! (clamped) -- the fb_ff scheme.  Three ranges the table does not cover,
  ! all of which used to return 0, i.e. no emission at all:
  !  - the gaps between pairs.  Each pair runs from just above one edge to
  !    just below the next, so the nodes (rounded to 1e-4 Ryd) leave a gap of
  !    about 2e-4 Ryd with the edge inside it: 10 in the H+ table, 7 in He+,
  !    1 in He2+, and five of the written-continuum frequencies fell in one.
  !    The edge is placed at the gap's midpoint, within 1e-4 Ryd of the true
  !    edge (the precision of the nodes themselves; for H at 1/121 Ryd it is
  !    3.5e-5 Ryd off), and each side extends its own pair.
  !  - from the last node to the ground edge: the last pair is extended
  !    (H+ 0.9994-0.99947, He+ 1.7999-1.8080, He2+ 3.9996-4.0011 Ryd).
  !  - below the first node.  He+ (table from 0.0911 Ryd, 1.0 micron): the
  !    levels whose edges lie there (n >= 4) are nearly hydrogenic, and a
  !    He+ capture into them weighs g_n/g_ion = 4n^2/2, the same 2n^2 as for
  !    H+, so the coefficient per ion is the H+ one; at the join the two
  !    tables agree to 0.4% (1e4 K) and 2.2% (5000 K).  Below the H+ table
  !    (0.0069 Ryd) and for He2+ below its table (0.0278 Ryd), free-free plus
  !    the captures into every level whose edge lies below h nu, summed in
  !    the Kramers approximation (bound-free Gaunt factor 1, the sum over n
  !    as an integral): with n1 the lowest such level and chi_n1 its binding
  !    energy, sum_{n >= n1} fb = [exp(chi_n1/kT) - 1] ff / g_ff
  !    (high_level_capture).  Just below the table's first node, 1e4 K, it
  !    is -0.3% (H+, 0.00689 against 0.0069 Ryd) and +3.6% (He2+, 0.02779
  !    against 0.0278 Ryd) off the table; free-free alone would be -6.6% and
  !    -30%.
  !=========================================================================
  recursive real(kind=wp) function gamma_tab(k, nu, T) result(g)
    implicit none
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: nu, T
    real(kind=wp), parameter :: NU_GROUND(3) = [NU_HI, NU_HEI, NU_HEII]
    integer :: i, ilo, n

    g = 0.0_wp
    n = gtab(k)%nnu
    if (nu >= NU_GROUND(k)) return
    if (nu < gtab(k)%nu(1)) then
       select case (k)
       case (1)
          g = high_level_capture(nu, T, 1)
       case (2)
          g = gamma_tab(1, nu, T)
       case (3)
          g = high_level_capture(nu, T, 2)
       end select
       return
    end if
    !--- the pair whose span, extended halfway across the gaps on either
    !--- side (and to the ground edge above the last node), contains nu
    ilo = n - 1
    do i = 1, n/2 - 1
       if (nu < 0.5_wp*(gtab(k)%nu(2*i) + gtab(k)%nu(2*i+1))) then
          ilo = 2*i - 1
          exit
       end if
    end do
    g = pair_value(k, ilo, nu, T)
  end function gamma_tab

  !--- Free-free plus Kramers free-bound into all levels n >= n1 of a
  !--- hydrogenic ion of charge Z, n1 the lowest level with Z^2/n^2 <= nu
  !--- [Ryd]: sum_n (2 chi_n/(n kT)) exp(chi_n/kT) taken as an integral over
  !--- u = chi_n/kT gives exp(chi_n1/kT) - 1 times the Kramers free-free,
  !--- i.e. that factor over the free-free Gaunt factor times gamma_ff.
  real(kind=wp) function high_level_capture(nu, T, Z) result(g)
    real(kind=wp), intent(in) :: nu, T
    integer,       intent(in) :: Z
    real(kind=wp) :: ff, gff, chi1
    integer :: n1
    ff  = gamma_ff(nu, T, Z)
    gff = ff/(FF_SPECTRAL*real(Z*Z,wp)*exp(-nu*HCRYD_K/T)/sqrt(T)*1.0e40_wp)
    n1  = ceiling(real(Z,wp)/sqrt(nu) - 1.0e-9_wp)
    chi1 = real(Z*Z,wp)/real(n1*n1,wp)
    g = ff*(1.0_wp + (exp(chi1*HCRYD_K/T) - 1.0_wp)/gff)
  end function high_level_capture

  !--- log gamma linear in nu through the pair (ilo, ilo+1), log-log in T
  !--- (clamped); used inside the pair and, extended, across a gap.
  real(kind=wp) function pair_value(k, ilo, nu, T) result(g)
    integer,       intent(in) :: k, ilo
    real(kind=wp), intent(in) :: nu, T
    integer :: iT
    real(kind=wp) :: w, lg1, lg2, lg
    if (T <= gtab(k)%tk(1)) then
       iT = 1;  w = 0.0_wp
    else if (T >= gtab(k)%tk(gtab(k)%nT)) then
       iT = gtab(k)%nT - 1;  w = 1.0_wp
    else
       do iT = 1, gtab(k)%nT - 1
          if (T < gtab(k)%tk(iT+1)) exit
       end do
       w = log10(T/gtab(k)%tk(iT)) / log10(gtab(k)%tk(iT+1)/gtab(k)%tk(iT))
    end if
    lg1 = gtab(k)%lg(ilo,  iT)*(1.0_wp - w) + gtab(k)%lg(ilo,  iT+1)*w
    lg2 = gtab(k)%lg(ilo+1,iT)*(1.0_wp - w) + gtab(k)%lg(ilo+1,iT+1)*w
    lg  = lg2 + (lg1 - lg2)*(gtab(k)%nu(ilo+1) - nu) &
          / (gtab(k)%nu(ilo+1) - gtab(k)%nu(ilo))
    g = 10.0_wp**lg
  end function pair_value

  !--- The frequencies [Ryd] where gamma_tab of table k is not smooth: the
  !--- table nodes and the edges gamma_tab places at the gap midpoints.  A
  !--- quadrature of it must break there.
  function gamma_tab_nodes(k) result(nu)
    integer, intent(in) :: k
    real(kind=wp), allocatable :: nu(:)
    integer :: i, n
    n = gtab(k)%nnu
    allocate(nu(n + n/2 - 1))
    nu(1:n) = gtab(k)%nu
    do i = 1, n/2 - 1
       nu(n+i) = 0.5_wp*(gtab(k)%nu(2*i) + gtab(k)%nu(2*i+1))
    end do
  end function gamma_tab_nodes

  !=========================================================================
  ! Milne free-bound to the ground level above the edge (fb_ff expression):
  ! gamma = 4 pi sigma statW hcRyd MILNE_C nu^3/T^1.5 exp(-(nu-nu0) Ryd/kT)
  ! in 1e-40 units.  statW = g_ground(recombined)/(2 g_ground(ion)).
  !=========================================================================
  real(kind=wp) function gamma_milne(nu, nu0, sigma_cm2, statW, T) result(g)
    implicit none
    real(kind=wp), intent(in) :: nu, nu0, sigma_cm2, statW, T
    g = 0.0_wp
    if (nu < nu0) return
    g = fourpi*sigma_cm2*statW*RYD_ERG*MILNE_C*nu**3/(T*sqrt(T)) &
        * exp(-(nu - nu0)*HCRYD_K/T) * 1.0e40_wp
  end function gamma_milne

  !=========================================================================
  ! Free-free coefficient [1e-40 erg cm^3/s/Hz], 4 pi j_nu / (n_e n_i):
  ! gamma_ff = FF_SPECTRAL Z^2 g_ff exp(-h nu/kT)/sqrt(T) * 1e40, with
  ! FF_SPECTRAL and g_ff from free_free_mod (until 2026-10-11 the rounded
  ! MOCASSIN 6.8391e-38 and the Hummer fit at every u, also above its
  ! stated u <= 10^1.5).
  !=========================================================================
  real(kind=wp) function gamma_ff(nu, T, Z) result(g)
    implicit none
    real(kind=wp), intent(in) :: nu, T
    integer,       intent(in) :: Z
    real(kind=wp) :: u
    u = nu*HCRYD_K/T
    g = FF_SPECTRAL*real(Z*Z,wp)*free_free_gaunt(real(Z,wp), T, u) &
        * exp(-u)/sqrt(T) * 1.0e40_wp
  end function gamma_ff

  !=========================================================================
  ! FREE-FREE EMISSION AS A DIFFUSE IONIZING SOURCE (diffuse_mod channel 7).
  !
  ! The free-free photons of the ionized gas above the band's lower edge
  ! E_lo = par%eion_min (13.598 eV by default) up to par%eion_max ionize H0,
  ! He0, He+ and the metals like any other band photon.  The emissivity is
  ! gamma_ff -- the coefficient of the written continuum -- so in the band
  !   L_Z = n_e n_Z V Int_{E_lo}^{E_hi} gamma_ff(nu, T, Z) 1e-40 dnu
  !       = n_e n_Z V G_Z(T) exp(-E_lo/kT),
  !   G_Z(T) = FF_SPECTRAL Z^2 T^(-1/2) (kT/h) Int_0^{x_hi} g_ff(u) e^(-x) dx,
  ! u = E_lo/kT + x, x_hi = min(60, (E_hi - E_lo)/kT).  The Boltzmann factor
  ! of the edge is applied at the leaf, so G_Z is smooth in T and is
  ! interpolated as ln G_Z by a Catmull-Rom cubic in log T.  Packet energies follow the
  ! energy emissivity itself (a packet carries energy; the photon rate it
  ! represents is L <1/E>, the physical one), from the CDF of
  ! g_ff(u) e^(-x), interpolated in log T and inverted once, continuous in E.
  !
  ! ENERGY.  The free-free cooling (cooling_mod, species_mod::metal_freefree,
  ! FF_TOTAL Z^2 gbar_ff) removes the whole emitted spectrum from the gas;
  ! the band part transported here is not removed again: absorbed, it heats
  ! the gas or ionizes, as stellar photons do, and what escapes leaves.
  ! Without it the band part simply escaped.  Its share of the free-free
  ! emission is ~ exp(-E_lo/kT): 1.4e-7 at 1e4 K, 3.7e-4 at 2e4 K.
  !
  ! VALIDITY.  Non-relativistic free-free (kT << m_e c^2) of a fully
  ! stripped-equivalent net charge Z (hydrogenic Gaunt factor); the
  ! Gaunt factor rule is free_free_mod's.
  !=========================================================================
  subroutine freefree_band_setup()
    implicit none
    real(kind=wp) :: t, T_K, kT_eV, x_hi, xa, xb, xm, fa, fm, fb, seg, acc, ulo
    integer :: iz, it, ix

    if (ff_ready) return
    ff_elo = par%eion_min
    ff_ehi = par%eion_max
    if (.not. allocated(ff_cdf)) allocate(ff_cdf(FF_NX, FF_NT, FF_ZMAX))
    do ix = 1, FF_NX
       t = real(ix-1, wp)/real(FF_NX-1, wp)
       ff_xg(ix) = FF_XMAX*t*t
    end do
    do it = 1, FF_NT
       ff_lgT(it) = log10(FF_TLO) + (log10(FF_THI) - log10(FF_TLO)) &
                    *real(it-1, wp)/real(FF_NT-1, wp)
    end do
    do iz = 1, FF_ZMAX
       do it = 1, FF_NT
          T_K   = 10.0_wp**ff_lgT(it)
          kT_eV = kboltz_cgs*T_K/ev2erg
          ulo   = ff_elo/kT_eV
          x_hi  = min(FF_XMAX, (ff_ehi - ff_elo)/kT_eV)
          acc = 0.0_wp
          ff_cdf(1, it, iz) = 0.0_wp
          fb = freefree_shape(iz, T_K, ulo, 0.0_wp)
          do ix = 1, FF_NX-1
             xa = ff_xg(ix)
             xb = min(ff_xg(ix+1), x_hi)
             if (xb > xa) then
                xm = 0.5_wp*(xa + xb)
                fa = fb                ! the previous segment's upper end
                fm = freefree_shape(iz, T_K, ulo, xm)
                fb = freefree_shape(iz, T_K, ulo, xb)
                seg = (xb - xa)/6.0_wp*(fa + 4.0_wp*fm + fb)
                acc = acc + seg
             end if
             ff_cdf(ix+1, it, iz) = acc
          end do
          if (acc > 0.0_wp) then
             ff_cdf(:, it, iz) = ff_cdf(:, it, iz)/acc
             ff_lncoef(it, iz) = log(FF_SPECTRAL*real(iz*iz,wp)/sqrt(T_K) &
                                     *kboltz_cgs*T_K/h_planck_cgs*acc)
          else
             ff_cdf(:, it, iz) = 1.0_wp
             ff_lncoef(it, iz) = -huge(1.0_wp)
          end if
          ff_cdf(FF_NX, it, iz) = 1.0_wp
       end do
    end do
    ff_ready = .true.
  end subroutine freefree_band_setup

  logical function freefree_band_is_ready() result(ok)
    ok = ff_ready
  end function freefree_band_is_ready

  !--- g_ff(u) e^(-x) at u = ulo + x, the free-free energy emissivity in the
  !--- band up to the factor FF_SPECTRAL Z^2 T^(-1/2) exp(-E_lo/kT)
  real(kind=wp) function freefree_shape(iz, T_K, ulo, x) result(f)
    integer,       intent(in) :: iz
    real(kind=wp), intent(in) :: T_K, ulo, x
    f = free_free_gaunt(real(iz,wp), T_K, ulo + x)*exp(-x)
  end function freefree_shape

  !--- temperature bracket on the log-T grid, clamped
  subroutine freefree_tbracket(T_K, it, w)
    real(kind=wp), intent(in)  :: T_K
    integer,       intent(out) :: it
    real(kind=wp), intent(out) :: w
    real(kind=wp) :: s
    s = (log10(max(T_K, tinest)) - ff_lgT(1))/(ff_lgT(2) - ff_lgT(1))
    s = min(max(s, 0.0_wp), real(FF_NT-1, wp))
    it = min(int(s) + 1, FF_NT-1)
    w  = s - real(it-1, wp)
  end subroutine freefree_tbracket

  !=========================================================================
  ! Band free-free luminosity per n_e n_Z [erg cm^3 s^-1] at T:
  ! G_Z(T) exp(-E_lo/kT) (0 for Z outside 1..FF_ZMAX).
  !=========================================================================
  real(kind=wp) function freefree_band_coefficient(iz, T_K) result(c)
    use mathlib, only : catmull_rom_uniform
    integer,       intent(in) :: iz
    real(kind=wp), intent(in) :: T_K
    real(kind=wp) :: w, lc
    integer :: it
    c = 0.0_wp
    if (iz < 1 .or. iz > FF_ZMAX) return
    call freefree_tbracket(T_K, it, w)
    lc = catmull_rom_uniform(ff_lncoef(:, iz), it, w)
    if (lc < -1.0e300_wp) return
    c  = exp(lc - ff_elo*ev2erg/(kboltz_cgs*T_K))
  end function freefree_band_coefficient

  !=========================================================================
  ! Photon energy [eV] of one free-free packet of net charge iz at T from
  ! the uniform u: the CDF interpolated in log T, inverted once, linear in
  ! x within the bracketing nodes.
  !=========================================================================
  real(kind=wp) function freefree_band_sample_energy(iz, T_K, u) result(E)
    integer,       intent(in) :: iz
    real(kind=wp), intent(in) :: T_K, u
    real(kind=wp) :: w, kT_eV, uu, ca, cb, cm, x
    integer :: it, lo, hi, mid
    call freefree_tbracket(T_K, it, w)
    kT_eV = kboltz_cgs*T_K/ev2erg
    uu = min(max(u, 0.0_wp), 1.0_wp)
    lo = 1;  hi = FF_NX
    ca = (1.0_wp - w)*ff_cdf(lo, it, iz) + w*ff_cdf(lo, it+1, iz)
    cb = (1.0_wp - w)*ff_cdf(hi, it, iz) + w*ff_cdf(hi, it+1, iz)
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       cm  = (1.0_wp - w)*ff_cdf(mid, it, iz) + w*ff_cdf(mid, it+1, iz)
       if (uu > cm) then
          lo = mid;  ca = cm
       else
          hi = mid;  cb = cm
       end if
    end do
    if (cb > ca) then
       x = ff_xg(lo) + (ff_xg(hi) - ff_xg(lo))*(uu - ca)/(cb - ca)
    else
       x = ff_xg(lo)
    end if
    E = min(ff_elo + x*kT_eV, ff_ehi*(1.0_wp - 1.0e-12_wp))
  end function freefree_band_sample_energy

end module recomb_continuum_mod
