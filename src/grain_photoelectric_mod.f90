module grain_photoelectric_mod
!---------------------------------------------------------------------------
! MoCHII: grain charge distribution, photoelectric heating of the gas and
! cooling of the gas by charged-particle accretion onto grains, after
! Weingartner & Draine (2001, ApJS 134, 263; WD01), for the grains of the run
! (grain_size_bins_mod) in the radiation field the transport has in
! each leaf.
!
! Physics (WD01 equation numbers):
!   thresholds   IP_V(Z) eq. 2; EA(Z) eq. 4 (carbonaceous), eq. 5 (silicate,
!                E_bg = 5 eV); E_min(Z<0) eq. 7; h nu_pet eq. 6;
!                h nu_pdt eq. 18;
!   yield        Y = y_2 min(y_0 y_1, 1), eqs. 8-17: Theta eq. 9, parabolic
!                photoelectron distribution eq. 10, y_2 eq. 11, y_1 eq. 13
!                (l_e = 10 A), y_0 eq. 16 (carbonaceous) / eq. 17 (silicate);
!   photodetachment  sigma_pdt eq. 19 with Delta E = 3 eV, f_pdt = 0.5, the
!                constant evaluated from eq. 19 (1.2197e-17 cm^2; eq. 20
!                prints it rounded to 1.2e-17), electron energy
!                h nu - h nu_pdt + E_min;
!   charge range Z_max eq. 22 (h nu_max = top of the transported band),
!                Z_min eqs. 23-24;
!   charging     f(Z+1)/f(Z) = [J_pe(Z) + J_ion(Z)]/J_e(Z+1), eq. 21;
!                J_pe eq. 25; accretion J_i eq. 26 with the Draine & Sutin
!                (1987, ApJ 320, 803; DS87) reduced rates J~ (their eqs.
!                3.3-3.5, theta_nu by eq. 2.4a);
!   sticking     s_e eqs. 28-30 (N_C = 468 (a/10 A)^3 for both materials),
!                s_i = 1 (sec. 3.2);
!   heating      Gamma_pe = Sum_Z f(Z) [Gamma'_pe,v(Z) + Gamma'_pd(Z)],
!                eqs. 38-40, with the photoelectron energy at infinity
!                averaged over eq. 10 between 0 (or E_min) and E_high;
!   cooling      Lambda_gr eq. 42: kT Lambda~ (DS87 eqs. 3.8-3.10) for each
!                colliding species, plus f(Z_min) J_e(Z_min) EA(Z_min).
!
! Colliders: electrons; H+, He+ (charge 1, 4 m_H), He2+ (charge 2), and the
! metal ions, as one species of charge 1 and mass 12 m_H (C+, which WD01 and
! Draine's charging code take for the metal ions), with the density that
! closes charge neutrality, n_e - n(H+) - n(He+) - 2 n(He2+), when positive.
!
! Two readings of WD01 that the paper leaves open, stated here:
!   (a) eq. 42's second term returns EA(Z_min) when a grain at Z_min captures
!       an electron and the state Z_min - 1 autoionizes.  With s_e(Z_min) = 0
!       (eq. 29) that capture never happens in eq. 21, so the term would
!       vanish; it is evaluated with the sticking of the neutral grain,
!       s_e(Z=0) (eq. 28), as the only reading in which the term exists, and
!       the energy the arriving electrons carry in is counted with it.
!   (b) an ion that reaches a grain is neutralized and returns to the gas as
!       an atom accommodated to the grain, carrying 2 k T_d per atom
!       (Hollenbach & McKee 1979, sec. IIIe: charged particles accommodate
!       with alpha of order unity).  eq. 42 counts the energy the ions bring
!       in; the 2 k T_d they take back is subtracted from Lambda_gr here.
!       T_d is the leaf's equilibrium mixture temperature (dust_temp_mod).
!
! Radiation field: the band tally jt_ion of the leaf, the bin-integrated
! 4 pi J_nu of every transported packet over [efuv_min, eion_max] (stellar,
! diffuse and dust-scattered), taken constant in nu inside a band bin; the
! energy grid of grain_size_bins_mod resolves the yield inside the
! bins.  Photons below par%efuv_min are not transported, so their
! photodetachment of negative grains and the photoemission of large
! carbonaceous grains between W = 4.4 eV and efuv_min are not counted;
! tests/grain_photoelectric measures that cut on the WD01 ISRF.
!
! Cost: the photoemission rates of a leaf depend on the field only, so they
! are computed once for each leaf and charge and kept while the thermal
! bisection varies T, n_e and the ions (photo_cache).  For grains whose charge
! window spans many states the rates are evaluated every CHARGE_STRIDE
! charges, with CHARGE_STRIDE e^2/a <= DE_STRIDE = 0.05 eV, and interpolated
! linearly between (they vary on the scale of an eV of threshold shift).
!
! Default on with dust (par%grain_pe; setup.f90): photoelectric emission and
! electron capture by grains are physically present wherever grains meet a
! radiation field and a plasma; WD01 find the photoelectric heating
! comparable to the photoionization heating in H II regions when
! G/n_H >~ 5 cm^3 (their sec. 7.1).
!---------------------------------------------------------------------------
  use define
  use grain_size_bins_mod, only : grain_bin_type, CARBONACEOUS, SILICATE, &
                                          gbin, n_gbin, egrid, n_egrid
  implicit none
  private

  public :: grain_photoelectric_setup, grain_field_load, grain_charge_terms
  public :: grain_bin_terms, grain_collider_type, grain_photoelectric_leaf
  public :: grain_bin_charge_range, grain_field_forget

  !--- one colliding species: density [cm^-3], charge [e], mass [g]
  type :: grain_collider_type
     real(kind=wp) :: n = 0.0_wp, q = 1.0_wp, m = 0.0_wp
  end type grain_collider_type

  real(kind=wp), parameter :: E2_EV_ANG = 14.39964548_wp      ! e^2 [eV Angstrom]
  real(kind=wp), parameter :: E_ESU = 4.80320471e-10_wp       ! e [esu]
  real(kind=wp), parameter :: AMU = 1.66053907e-24_wp         ! [g]
  real(kind=wp), parameter :: M_PROTON = 1.67262192e-24_wp    ! [g]
  real(kind=wp), parameter :: W_CARB = 4.4_wp, W_SIL = 8.0_wp, EBG_SIL = 5.0_wp
  real(kind=wp), parameter :: L_ESCAPE_ANG = 10.0_wp
  real(kind=wp), parameter :: DE_PDT = 3.0_wp, F_PDT = 0.5_wp
  real(kind=wp), parameter :: DE_STRIDE = 0.05_wp
  !--- f(Z) is followed out from its mode until it falls below this
  real(kind=wp), parameter :: LN_F_FLOOR = -20.0_wp
  real(kind=wp) :: sigma_pdt0 = 0.0_wp   ! eq. 19 constant [cm^2]
  real(kind=wp) :: e_top = 0.0_wp        ! h nu_max [eV]

  !--- per bin: charge limits, stride, and the photo rates of the loaded field
  type :: photo_cache_type
     integer :: zmin = 0, zmax = 0, stride = 1, nnode = 0
     integer :: zmode_last = 0      ! the last mode found: the search starts there
     integer,       allocatable :: stamp(:)
     real(kind=wp), allocatable :: jpe(:), gpe(:)
  end type photo_cache_type
  type(photo_cache_type), allocatable :: cache(:)
  real(kind=wp), allocatable :: nphot_loaded(:)
  integer :: field_stamp = 0
  integer :: leaf_loaded = -1
  logical :: ready = .false.

contains

  !=========================================================================
  ! Charge limits and caches for the bins of grain_size_bins_mod.
  !=========================================================================
  subroutine grain_photoelectric_setup()
    integer :: ib, zmin, zmax, s
    real(kind=wp) :: e2a

    sigma_pdt0 = twopi*E_ESU**2*h_planck_cgs*F_PDT &
                 /(3.0_wp*me_cgs*clight_cgs*DE_PDT*ev2erg)
    e_top = egrid(n_egrid)
    if (allocated(cache)) deallocate(cache)
    allocate(cache(n_gbin))
    do ib = 1, n_gbin
       call grain_bin_charge_range(gbin(ib), zmin, zmax)
       e2a = E2_EV_ANG/gbin(ib)%a_ang
       s = max(1, int(DE_STRIDE/e2a))
       cache(ib)%zmin = zmin;  cache(ib)%zmax = zmax;  cache(ib)%stride = s
       cache(ib)%nnode = (zmax - zmin)/s + 2
       allocate(cache(ib)%stamp(0:cache(ib)%nnode), cache(ib)%jpe(0:cache(ib)%nnode), &
                cache(ib)%gpe(0:cache(ib)%nnode))
       cache(ib)%stamp = -1
    end do
    if (allocated(nphot_loaded)) deallocate(nphot_loaded)
    allocate(nphot_loaded(n_egrid))
    nphot_loaded = 0.0_wp
    field_stamp = 0
    leaf_loaded = -1
    ready = .true.
  end subroutine grain_photoelectric_setup

  !--- int[x] of WD01: the greatest integer less than x
  pure integer function int_below(x) result(n)
    real(kind=wp), intent(in) :: x
    n = ceiling(x) - 1
  end function int_below

  !=========================================================================
  ! Z_min (eqs. 23-24) and Z_max (eq. 22, h nu_max = top of the energy grid)
  !=========================================================================
  subroutine grain_bin_charge_range(b, zmin, zmax)
    type(grain_bin_type), intent(in) :: b
    integer, intent(out) :: zmin, zmax
    real(kind=wp) :: a, uait, w
    a = b%a_ang
    if (b%material == CARBONACEOUS) then
       uait = -(3.9_wp + 0.12_wp*a + 2.0_wp/a)
       w = W_CARB
    else
       uait = -(2.5_wp + 0.07_wp*a + 8.0_wp/a)
       w = W_SIL
    end if
    zmin = int_below(uait/E2_EV_ANG*a) + 1
    zmax = int_below(((egrid(n_egrid) - w)/E2_EV_ANG*a + 0.5_wp - 0.3_wp/a) &
                     /(1.0_wp + 0.3_wp/a))
    zmax = max(zmax, 0)
  end subroutine grain_bin_charge_range

  !=========================================================================
  ! Load a radiation field: photon rates on the energy grid [cm^-2 s^-1].
  ! Every cached photoemission rate is invalidated.
  !=========================================================================
  subroutine grain_field_load(nphot)
    real(kind=wp), intent(in) :: nphot(:)
    nphot_loaded = nphot
    field_stamp = field_stamp + 1
    leaf_loaded = -1
  end subroutine grain_field_load

  !--- the band tally changes between iterations and leaves are renumbered by
  !--- re-refinement: the next leaf loads its field afresh
  subroutine grain_field_forget()
    leaf_loaded = -1
  end subroutine grain_field_forget

  !=========================================================================
  ! The field of leaf il from the band tally, loaded once for each leaf
  ! (photon_rates_from_band); then the heating and cooling at the state.
  ! Returns [erg cm^-3 s^-1]: gamma = photoelectric + photodetachment heating,
  ! lambda = accretion cooling net of the energy the neutralized ions return.
  !=========================================================================
  subroutine grain_photoelectric_leaf(il, T, Td, nH, ne, nHII, nHeII, nHeIII, &
                                      xHI, xHII, gamma, lambda)
    use jtally_mod,   only : jt_ion
    use ion_band_mod, only : ion_eedge, nnu_band
    use octree_mod,   only : leaf_half
    use grain_size_bins_mod, only : photon_rates_from_band, grain_leaf_scale
    integer,       intent(in)  :: il
    real(kind=wp), intent(in)  :: T, Td, nH, ne, nHII, nHeII, nHeIII, xHI, xHII
    real(kind=wp), intent(out) :: gamma, lambda
    real(kind=wp) :: s_other, s_pah, fac, nphot(n_egrid), fband(nnu_band)
    type(grain_collider_type) :: ions(4)

    gamma = 0.0_wp;  lambda = 0.0_wp
    call grain_leaf_scale(il, nH, xHI, xHII, s_other, s_pah)
    if (s_other <= 0.0_wp .and. s_pah <= 0.0_wp) return
    if (il /= leaf_loaded) then
       fac = 1.0_wp/((2.0_wp*leaf_half(il))**3*par%distance2cm**2)
       fband = jt_ion(1:nnu_band, il)*fac
       call photon_rates_from_band(fband, ion_eedge(1:nnu_band+1), nphot)
       call grain_field_load(nphot)
       leaf_loaded = il
    end if
    ions(1) = grain_collider_type(nHII, 1.0_wp, M_PROTON)
    ions(2) = grain_collider_type(nHeII, 1.0_wp, 4.002602_wp*AMU)
    ions(3) = grain_collider_type(nHeIII, 2.0_wp, 4.002602_wp*AMU)
    ions(4) = grain_collider_type(max(ne - nHII - nHeII - 2.0_wp*nHeIII, 0.0_wp), &
                                  1.0_wp, 12.0_wp*AMU)
    call grain_charge_terms(T, Td, ne, ions, nH*s_other, nH*s_pah, gamma, lambda)
  end subroutine grain_photoelectric_leaf

  !=========================================================================
  ! Sum over the bins of the loaded field: [erg cm^-3 s^-1] with dens_other
  ! (dens_pah) the H density times the abundance scale of the other (PAH)
  ! bins.
  !=========================================================================
  subroutine grain_charge_terms(T, Td, ne, ions, dens_other, dens_pah, gamma, lambda)
    real(kind=wp), intent(in) :: T, Td, ne, dens_other, dens_pah
    type(grain_collider_type), intent(in) :: ions(:)
    real(kind=wp), intent(out) :: gamma, lambda
    real(kind=wp) :: g1, l1, d
    integer :: ib
    gamma = 0.0_wp;  lambda = 0.0_wp
    do ib = 1, n_gbin
       if (gbin(ib)%pah) then
          d = dens_pah
       else
          d = dens_other
       end if
       if (d <= 0.0_wp) cycle
       call grain_bin_terms(ib, T, Td, ne, ions, g1, l1)
       gamma  = gamma  + d*gbin(ib)%n_per_H*g1
       lambda = lambda + d*gbin(ib)%n_per_H*l1
    end do
  end subroutine grain_charge_terms

  !=========================================================================
  ! One grain of bin ib at (T, n_e, ions) in the loaded field: the charge
  ! distribution, then the heating g1 and cooling l1 per grain [erg s^-1].
  ! Optional: zmean = <Z>, f_neutral = f(Z=0).
  !=========================================================================
  subroutine grain_bin_terms(ib, T, Td, ne, ions, g1, l1, zmean, f_neutral)
    integer,       intent(in)  :: ib
    real(kind=wp), intent(in)  :: T, Td, ne
    type(grain_collider_type), intent(in) :: ions(:)
    real(kind=wp), intent(out) :: g1, l1
    real(kind=wp), intent(out), optional :: zmean, f_neutral
    real(kind=wp) :: a, area, kT, tau_e, ve, se0, sep, ncarb
    real(kind=wp) :: vi(size(ions))
    real(kind=wp) :: lnf, fsum, gsum, lsum, zsum, jpe, gpe, up, down, fzero
    integer :: zmin, zmax, lo, hi, mid, zmode, z, i, step, s
    real(kind=wp) :: kap

    associate (c => cache(ib))
    zmin = c%zmin;  zmax = c%zmax
    a = gbin(ib)%a_cm
    area = pi*a*a
    kT = kboltz_cgs*T
    tau_e = a*kT/E_ESU**2
    ve = sqrt(8.0_wp*kT/(pi*me_cgs))
    do i = 1, size(ions)
       vi(i) = sqrt(8.0_wp*kT/(pi*ions(i)%m))
    end do
    ncarb = 468.0_wp*(gbin(ib)%a_ang/10.0_wp)**3
    sep = 0.5_wp*(1.0_wp - exp(-gbin(ib)%a_ang/L_ESCAPE_ANG))
    se0 = sep/(1.0_wp + exp(min(20.0_wp - ncarb, 700.0_wp)))

    !--- mode: the largest Z with f(Z) > f(Z-1), i.e. up(Z-1) > down(Z);
    !--- up(Z) falls and down(Z) rises with Z, so the test is monotone.  The
    !--- bracket is found by doubling steps from the last mode of this bin
    !--- (the thermal bisection moves it little), then bisected.
    z = min(max(c%zmode_last, zmin), zmax)
    if (rising(z)) then
       lo = z;  s = 1
       do
          hi = min(lo + s, zmax)
          if (hi == lo) exit
          if (.not. rising(hi)) then
             hi = hi - 1;  exit
          end if
          lo = hi;  s = 2*s
       end do
    else
       hi = z - 1;  s = 1
       do
          lo = max(hi - s, zmin)
          if (lo == zmin) exit
          if (rising(lo)) exit
          hi = lo - 1;  s = 2*s
       end do
    end if
    do while (hi > lo)
       mid = (lo + hi + 1)/2
       if (rising(mid)) then
          lo = mid
       else
          hi = mid - 1
       end if
    end do
    zmode = lo
    c%zmode_last = zmode

    !--- f by eq. 21 outward from the mode, in logs.  Where the distribution
    !--- spans many charges (large grains) it is sampled every `step`
    !--- charges, each sample standing for `step` states and ln f advanced by
    !--- step ln r at the middle of the step; step is at most the photo-rate
    !--- stride and at most half the width sigma = kappa^(-1/2),
    !--- kappa = -d ln r/dZ at the mode, so the error is O((step/sigma)^2).
    step = 1
    if (c%stride > 1 .and. zmode - c%stride > zmin .and. zmode + c%stride < zmax) then
       kap = (log_ratio(zmode - c%stride) - log_ratio(zmode + c%stride)) &
             /real(2*c%stride, wp)
       if (kap > 0.0_wp) step = max(1, min(c%stride, int(0.5_wp/sqrt(kap))))
    end if
    fsum = 0.0_wp;  gsum = 0.0_wp;  lsum = 0.0_wp;  zsum = 0.0_wp;  fzero = 0.0_wp
    call accumulate(zmode, 0.0_wp, real(step, wp))
    lnf = 0.0_wp
    z = zmode
    do while (z < zmax)
       s = min(step, zmax - z)
       up = rate_up(z + (s-1)/2);  down = rate_down(z + (s-1)/2 + 1)
       if (up <= 0.0_wp) exit
       lnf = lnf + real(s, wp)*log(up/down)
       if (lnf < LN_F_FLOOR) exit
       z = z + s
       call accumulate(z, lnf, real(s, wp))
    end do
    lnf = 0.0_wp
    z = zmode
    do while (z > zmin)
       s = min(step, z - zmin)
       down = rate_down(z - (s+1)/2 + 1);  up = rate_up(z - (s+1)/2)
       if (down <= 0.0_wp) exit
       lnf = lnf + real(s, wp)*log(down/up)
       if (lnf < LN_F_FLOOR) exit
       z = z - s
       call accumulate(z, lnf, real(s, wp))
    end do
    g1 = gsum/fsum
    l1 = lsum/fsum
    if (present(zmean)) zmean = zsum/fsum
    if (present(f_neutral)) f_neutral = fzero/fsum
    end associate

  contains

    !--- f(Z) = exp(lnf) relative to the mode: heating and cooling sums
    subroutine accumulate(zz, lnfz, weight)
      integer,       intent(in) :: zz
      real(kind=wp), intent(in) :: lnfz, weight
      real(kind=wp) :: ff, cool, ion_flux
      integer :: k
      ff = weight*exp(lnfz)
      call photo_rates(ib, zz, jpe, gpe)
      !--- kinetic energy brought in by the electrons and ions that stick
      cool = ne*electron_sticking(zz)*ve*area*kT &
             *reduced_cooling(tau_e, -real(zz, wp))
      ion_flux = 0.0_wp
      do k = 1, size(ions)
         if (ions(k)%n <= 0.0_wp) cycle
         cool = cool + ions(k)%n*vi(k)*area*kT &
                *reduced_cooling(tau_e/ions(k)%q**2, real(zz, wp)/ions(k)%q)
         ion_flux = ion_flux + ions(k)%n*vi(k)*area &
                *reduced_rate(tau_e/ions(k)%q**2, real(zz, wp)/ions(k)%q)
      end do
      !--- (b) the neutralized ions return as atoms at 2 k T_d
      cool = cool - ion_flux*2.0_wp*kboltz_cgs*Td
      !--- (a) eq. 42's second term at Z_min
      if (zz == zmin .and. zmin < 0) then
         cool = cool + ne*se0*ve*area*( kT*reduced_cooling(tau_e, -real(zz, wp)) &
                + reduced_rate(tau_e, -real(zz, wp))*electron_affinity(zz)*ev2erg )
      end if
      fsum = fsum + ff
      gsum = gsum + ff*gpe
      lsum = lsum + ff*cool
      zsum = zsum + ff*real(zz, wp)
      if (zz == 0) fzero = ff
    end subroutine accumulate

    !--- f(Z) > f(Z-1) (Z_min counts as rising: nothing lies below it)
    logical function rising(zz) result(r)
      integer, intent(in) :: zz
      if (zz <= zmin) then
         r = .true.
      else
         r = rate_up(zz-1) > rate_down(zz)
      end if
    end function rising

    !--- ln [f(Z+1)/f(Z)]
    real(kind=wp) function log_ratio(zz) result(lr)
      integer, intent(in) :: zz
      real(kind=wp) :: u, d
      u = rate_up(zz);  d = rate_down(zz+1)
      if (u > 0.0_wp .and. d > 0.0_wp) then
         lr = log(u/d)
      else
         lr = 0.0_wp
      end if
    end function log_ratio

    !--- rate out of Z upward: photoemission + photodetachment + ion capture
    real(kind=wp) function rate_up(zz) result(r)
      integer, intent(in) :: zz
      integer :: k
      call photo_rates(ib, zz, jpe, gpe)
      r = jpe
      do k = 1, size(ions)
         if (ions(k)%n <= 0.0_wp) cycle
         r = r + ions(k)%n*vi(k)*area &
             *reduced_rate(tau_e/ions(k)%q**2, real(zz, wp)/ions(k)%q)
      end do
    end function rate_up

    !--- rate into Z from Z+1... i.e. electron capture by a grain of charge zz
    real(kind=wp) function rate_down(zz) result(r)
      integer, intent(in) :: zz
      r = ne*electron_sticking(zz)*ve*area*reduced_rate(tau_e, -real(zz, wp))
    end function rate_down

    !--- WD01 eqs. 28-30
    real(kind=wp) function electron_sticking(zz) result(s)
      integer, intent(in) :: zz
      if (zz > 0) then
         s = sep
      else if (zz > zmin) then
         s = se0
      else
         s = 0.0_wp
      end if
    end function electron_sticking

    !--- EA(Z) [eV], eqs. 4-5
    real(kind=wp) function electron_affinity(zz) result(ea)
      integer, intent(in) :: zz
      real(kind=wp) :: e2a
      e2a = E2_EV_ANG/gbin(ib)%a_ang
      if (gbin(ib)%material == CARBONACEOUS) then
         ea = W_CARB + (real(zz, wp) - 0.5_wp)*e2a &
              - e2a*4.0_wp/(gbin(ib)%a_ang + 7.0_wp)
      else
         ea = W_SIL - EBG_SIL + (real(zz, wp) - 0.5_wp)*e2a
      end if
    end function electron_affinity

  end subroutine grain_bin_terms

  !=========================================================================
  ! Photoemission rate J_pe [s^-1] and heating Gamma'_pe + Gamma'_pd
  ! [erg s^-1] of one grain of bin ib at charge z in the loaded field; nodes
  ! every stride charges, cached, linear in z between nodes.
  !=========================================================================
  subroutine photo_rates(ib, z, jpe, gpe)
    integer,       intent(in)  :: ib, z
    real(kind=wp), intent(out) :: jpe, gpe
    integer :: j, s, z0
    real(kind=wp) :: w, j0, g0, j1, g1
    associate (c => cache(ib))
    s = c%stride
    if (s == 1) then
       call node(z - c%zmin, z, jpe, gpe)
       return
    end if
    j = (z - c%zmin)/s
    z0 = c%zmin + j*s
    if (z0 == z) then
       call node(j, z0, jpe, gpe)
       return
    end if
    call node(j, z0, j0, g0)
    call node(j+1, z0 + s, j1, g1)
    w = real(z - z0, wp)/real(s, wp)
    jpe = (1.0_wp - w)*j0 + w*j1
    gpe = (1.0_wp - w)*g0 + w*g1
    end associate
  contains
    subroutine node(jn, zn, jv, gv)
      integer,       intent(in)  :: jn, zn
      real(kind=wp), intent(out) :: jv, gv
      associate (c => cache(ib))
      if (c%stamp(jn) /= field_stamp) then
         call photoemission(gbin(ib), zn, nphot_loaded, c%jpe(jn), c%gpe(jn))
         c%stamp(jn) = field_stamp
      end if
      jv = c%jpe(jn);  gv = c%gpe(jn)
      end associate
    end subroutine node
  end subroutine photo_rates

  !=========================================================================
  ! WD01 eqs. 2-20, 25, 39-40 for one grain at charge z in the photon rates
  ! nphot [cm^-2 s^-1] on the energy grid.
  !=========================================================================
  subroutine photoemission(b, z, nphot, jpe, gpe)
    type(grain_bin_type), intent(in) :: b
    integer,              intent(in) :: z
    real(kind=wp),        intent(in) :: nphot(:)
    real(kind=wp),        intent(out) :: jpe, gpe
    real(kind=wp) :: a, e2a, w, ipv, emin, epet, epdt, ea1, phi, x
    real(kind=wp) :: e, theta, t, y0, eh, el, y2, ymin, ke, cabs, sig, epd
    real(kind=wp) :: zz
    integer :: k, k0
    logical :: carb

    a = b%a_ang
    zz = real(z, wp)
    e2a = E2_EV_ANG/a
    carb = b%material == CARBONACEOUS
    w = merge(W_CARB, W_SIL, carb)
    ipv = w + (zz + 0.5_wp)*e2a + (zz + 2.0_wp)*e2a*0.3_wp/a
    if (z < -1) then
       emin = -(zz + 1.0_wp)*e2a/(1.0_wp + (27.0_wp/a)**0.75_wp)
    else
       emin = 0.0_wp
    end if
    if (z >= -1) then
       epet = ipv
    else
       epet = max(ipv + emin, 0.0_wp)
    end if
    epdt = huge(1.0_wp)
    if (z < 0) then
       if (carb) then
          ea1 = W_CARB + (zz + 0.5_wp)*e2a - e2a*4.0_wp/(a + 7.0_wp)
       else
          ea1 = W_SIL - EBG_SIL + (zz + 0.5_wp)*e2a
       end if
       epdt = max(ea1 + emin, 0.0_wp)
    end if
    phi = (zz + 1.0_wp)*e2a
    jpe = 0.0_wp;  gpe = 0.0_wp
    !--- first grid energy above the lower threshold
    k0 = first_above(min(epet, epdt))
    do k = k0, n_egrid
       e = egrid(k)
       if (nphot(k) <= 0.0_wp) cycle
       if (e > epet) then
          if (z >= 0) then
             theta = e - epet + phi
          else
             theta = e - epet
          end if
          t = theta/w
          if (carb) then
             t = t*t*t*t*t
             y0 = 9.0e-3_wp*t/(1.0_wp + 3.7e-2_wp*t)
          else
             y0 = 0.5_wp*t/(1.0_wp + 5.0_wp*t)
          end if
          ymin = min(y0*b%y1(k), 1.0_wp)
          if (z == 0) then
             cabs = b%cabs_neutral(k)
          else
             cabs = b%cabs_charged(k)
          end if
          if (z >= 0) then
             !--- eq. 10 between E_low = -phi and E_high, escape above 0:
             !--- y_2 = E_h^2 (E_h - 3 E_l)/Theta^3 and the escaping
             !--- energy Int_0^E_h E f dE = E_h^3 (E_h - 2 E_l)/(2 Theta^3)
             eh = e - epet;  el = -phi
             y2 = eh*eh*(eh - 3.0_wp*el)/theta**3
             ke = eh**3*(eh - 2.0_wp*el)/(2.0_wp*theta**3)
             jpe = jpe + nphot(k)*cabs*y2*ymin
             gpe = gpe + nphot(k)*cabs*ymin*ke
          else
             !--- every attempting electron escapes; mean energy midway
             !--- between E_min and E_min + h nu - h nu_pet
             jpe = jpe + nphot(k)*cabs*ymin
             gpe = gpe + nphot(k)*cabs*ymin*(emin + 0.5_wp*(e - epet))
          end if
       end if
       if (z < 0 .and. e > epdt) then
          x = (e - epdt)/DE_PDT
          sig = sigma_pdt0*abs(zz)*x/(1.0_wp + x*x/3.0_wp)**2
          epd = e - epdt + emin
          jpe = jpe + nphot(k)*sig
          gpe = gpe + nphot(k)*sig*epd
       end if
    end do
    gpe = gpe*ev2erg
  end subroutine photoemission

  pure integer function first_above(e0) result(k)
    real(kind=wp), intent(in) :: e0
    integer :: lo, hi, mid
    if (egrid(1) > e0) then
       k = 1;  return
    end if
    lo = 1;  hi = n_egrid + 1
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       if (egrid(mid) > e0) then
          hi = mid
       else
          lo = mid
       end if
    end do
    k = hi
  end function first_above

  !=========================================================================
  ! DS87 reduced collision rate J~(tau, nu), nu = Z e/q: eqs. 3.3-3.5, with
  ! theta_nu = nu/(1 + nu^(-1/2)) (eq. 2.4a)
  !=========================================================================
  pure real(kind=wp) function reduced_rate(tau, nu) result(j)
    real(kind=wp), intent(in) :: tau, nu
    real(kind=wp) :: th
    if (nu == 0.0_wp) then
       j = 1.0_wp + sqrt(pi/(2.0_wp*tau))
    else if (nu < 0.0_wp) then
       j = (1.0_wp - nu/tau)*(1.0_wp + sqrt(2.0_wp/(tau - 2.0_wp*nu)))
    else
       th = nu/(1.0_wp + 1.0_wp/sqrt(nu))
       if (th/tau > 700.0_wp) then
          j = 0.0_wp
       else
          j = (1.0_wp + 1.0_wp/sqrt(4.0_wp*tau + 3.0_wp*nu))**2*exp(-th/tau)
       end if
    end if
  end function reduced_rate

  !=========================================================================
  ! DS87 reduced cooling Lambda~(tau, nu): eqs. 3.8-3.10
  !=========================================================================
  pure real(kind=wp) function reduced_cooling(tau, nu) result(l)
    real(kind=wp), intent(in) :: tau, nu
    real(kind=wp) :: th
    if (nu == 0.0_wp) then
       l = 2.0_wp + 1.5_wp*sqrt(pi/(2.0_wp*tau))
    else if (nu < 0.0_wp) then
       l = (2.0_wp - nu/tau)*(1.0_wp + 1.0_wp/sqrt(tau - nu))
    else
       th = nu/(1.0_wp + 1.0_wp/sqrt(nu))
       if (th/tau > 700.0_wp) then
          l = 0.0_wp
       else
          l = (2.0_wp + nu/tau)*(1.0_wp + 1.0_wp/sqrt(1.5_wp/tau + 3.0_wp*nu)) &
              *exp(-th/tau)
       end if
    end if
  end function reduced_cooling

end module grain_photoelectric_mod
