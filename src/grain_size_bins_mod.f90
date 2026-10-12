module grain_size_bins_mod
!---------------------------------------------------------------------------
! MoCHII: the run's grains as discrete size bins, for the gas-grain physics
! (grain_photoelectric_mod: charge distribution, photoelectric heating and
! charged-particle accretion cooling; gas_grain_collision_mod: neutral-atom
! energy exchange).
!
! The grains are those of the SEDust model the run carries (par%dust_model,
! grain_model_mod::dmodel): every population's size grid, number per H in each
! size bin (dn), and absorption cross section C_abs(lambda, a), so the grains
! that heat the gas are the grains that absorb and reemit the field.
!
! Material.  Weingartner & Draine (2001, ApJS 134, 263; WD01) distinguish two
! photoemission materials:
!   carbonaceous (W = 4.4 eV, eqs. 2, 4, 16, 23): the SEDust populations 'pah'
!     and 'gra' (the PAH-graphite continuum of DL07 and astrodust, MRN and
!     Zubko graphite, Zubko PAH);
!   silicate (W = 8 eV, E_bg = 5 eV, eqs. 2, 5, 17, 23): the populations
!     'sil' -- the silicates of DL07, MRN and Zubko, and the astrodust grains
!     of Hensley & Draine (2023).
! ASTRODUST IS GIVEN THE SILICATE PHOTOEMISSION PARAMETERS.  No photoelectric
! threshold or yield has been measured or calculated for the astrodust
! composite (silicate with a carbonaceous and Fe-bearing minority);
! amorphous silicate is its dominant component, so WD01's silicate W, E_bg and
! y_0 are taken, with the astrodust dielectric function for the photon
! attenuation length.  Validity: an assumption, not a published model; the
! measured share of the heating that the astrodust grains carry is printed at
! setup by grain_photoelectric_mod, and in the test cases the PAHs carry most
! of the heating.
!
! The photon attenuation length l_a of the yield enhancement y_1 (WD01 eqs.
! 13-15) comes from the dielectric function the model's optics were computed
! on: graphite (D03, 1/3 parallel + 2/3 perpendicular, eq. 15) for every
! carbonaceous grain -- WD01: "we use graphite dielectric functions for all
! carbonaceous grains" -- D03 astronomical silicate for 'sil', and the DH21
! astrodust function for astrodust.
!
! Size bins.  The native grids (astrodust 169 radii) are merged into bins of
! grain_bins_per_decade per decade in a; a bin carries its number per H, its
! area-weighted radius (Sum dn a^2 / Sum dn)^(1/2), and its number-weighted
! cross sections.  tests/grain_photoelectric measures what the merge costs.
! The two charge-state populations of the astrodust and DL07 PAH-graphite
! grains (neutral and cation, same radii) are one grain population for the
! charge distribution: their numbers are summed, and the neutral cross section
! is used at Z = 0 and the cation one at Z /= 0, as WD01 (footnote 5) do.
!
! Energy grid.  The radiation field a leaf has is the band tally jt_ion
! (bin-integrated 4 pi J_nu over the transported band, the same input SEDust
! takes); inside a band bin J_nu is taken constant in nu, and the bin is cut
! into sub-intervals no wider than max(DE_SUB_ABS, DE_SUB_REL E) so that the
! yield, which rises steeply above each threshold, is resolved continuously.
!
! Abundance in a leaf.  The leaf grey opacity is rhokap = n_dust C_ext(ref)
! per code length, so the grains per H in a leaf are the model's dn times
! (rhokap/distance2cm) / (n_H C_ext,model(lambda_ref)), C_ext,model the
! model's own extinction at par%lambda_ref (with par%ion_dust_kext the grains
! are the model's, scaled to that file's extinction at lambda_ref).  The PAH
! bins are always weighted by the PAH survival (x_HI + f_ion_pah x_HII), the
! rule sedust_mod applies to the PAH emission channel under par%sed_pah_live,
! divided, under laursen09_live, by the dust survival already in rhokap.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  integer, parameter, public :: CARBONACEOUS = 1, SILICATE = 2

  type, public :: grain_bin_type
     integer       :: material = CARBONACEOUS
     logical       :: pah = .false.     ! part of the model's PAH channel
     real(kind=wp) :: a_cm  = 0.0_wp    ! radius [cm]
     real(kind=wp) :: a_ang = 0.0_wp    ! radius [Angstrom]
     real(kind=wp) :: n_per_H = 0.0_wp  ! grains per H atom of the model
     !--- on the energy grid egrid(:): absorption cross section [cm^2] of the
     !--- neutral and of the charged grain, and WD01 eq. 13 y_1
     real(kind=wp), allocatable :: cabs_neutral(:), cabs_charged(:), y1(:)
  end type grain_bin_type

  public :: grain_size_bins_setup, grain_bins_from_model
  public :: grain_energy_grid_from_edges, grain_leaf_scale, photon_rates_from_band
  public :: grain_bins_ready

  type(grain_bin_type), allocatable, public :: gbin(:)
  integer, public :: n_gbin = 0
  !--- energy grid [eV]: sub-interval centers, widths, and the band bin each
  !--- lies in
  real(kind=wp), allocatable, public :: egrid(:), egrid_de(:)
  integer,       allocatable, public :: egrid_band(:)
  integer, public :: n_egrid = 0
  !--- the model's C_ext(lambda_ref) per H [cm^2/H]
  real(kind=wp), public :: cext_ref_model = 0.0_wp
  logical :: grain_bins_ready = .false.

  !--- size bins per decade of radius (the size-bin convergence of
  !--- tests/grain_photoelectric sets it)
  integer, public :: grain_bins_per_decade = 8
  real(kind=wp), parameter :: DE_SUB_ABS = 0.2_wp, DE_SUB_REL = 0.02_wp
  integer,       parameter :: MAX_SUB = 64
  real(kind=wp), parameter :: HC_EV_CM = 1.23984198e-4_wp   ! h c [eV cm]
  real(kind=wp), parameter :: L_ESCAPE = 1.0e-7_wp          ! l_e = 10 A, WD01 sec. 2.3.1

contains

  !=========================================================================
  ! Setup for a MoCHII run: the energy grid on the band bins and the bins of
  ! the grain model named by par%dust_model.
  !=========================================================================
  subroutine grain_size_bins_setup()
    use grain_model_mod, only : build_grain_model, dmodel, grain_extinction_table
    use ion_band_mod,    only : ion_eedge, nnu_band
    implicit none
    real(kind=wp), allocatable :: lam(:), alb(:), cosg(:), cext(:)
    integer :: n, i
    real(kind=wp) :: w

    if (grain_bins_ready) return
    call build_grain_model()
    call grain_energy_grid_from_edges(ion_eedge(1:nnu_band+1))
    call grain_bins_from_model(dmodel, trim(par%sed_data_dir))
    !--- the model's own extinction at lambda_ref, the normalization of rhokap
    call grain_extinction_table(lam, alb, cosg, cext, n)
    i = 1
    do while (i < n - 1 .and. lam(i+1) < par%lambda_ref)
       i = i + 1
    end do
    w = log(par%lambda_ref/lam(i))/log(lam(i+1)/lam(i))
    cext_ref_model = exp((1.0_wp - w)*log(cext(i)) + w*log(cext(i+1)))
    deallocate(lam, alb, cosg, cext)
    grain_bins_ready = .true.
    if (mpar%p_rank == 0) then
       write(*,'(a,i0,3a,i0,a)') ' GRAIN: ', n_gbin, ' size bins of the ''', &
          trim(par%dust_model), ''' grains on ', n_egrid, &
          ' energies of the transported band'
       write(*,'(a,es11.4,a,i0,a,i0,a)') ' GRAIN: model C_ext/H at lambda_ref = ', &
          cext_ref_model, ' cm^2; carbonaceous bins ', &
          count(gbin(:)%material == CARBONACEOUS), ', silicate bins ', &
          count(gbin(:)%material == SILICATE), &
          ' (astrodust takes the WD01 silicate photoemission parameters)'
       if (len_trim(par%ion_dust_kext) > 0) write(*,'(3a)') &
          ' GRAIN: par%ion_dust_kext sets the extinction; the gas-grain '// &
          'physics uses the ''', trim(par%dust_model), &
          ''' grains scaled to it at lambda_ref'
    end if
  end subroutine grain_size_bins_setup

  !=========================================================================
  ! Energy grid from ascending band edges [eV]: each band bin is cut into
  ! equal sub-intervals no wider than max(DE_SUB_ABS, DE_SUB_REL E).
  !=========================================================================
  subroutine grain_energy_grid_from_edges(edges)
    real(kind=wp), intent(in) :: edges(:)
    integer :: nb, b, j, k, nsub(size(edges)-1)
    real(kind=wp) :: de

    nb = size(edges) - 1
    do b = 1, nb
       de = edges(b+1) - edges(b)
       nsub(b) = min(MAX_SUB, max(1, ceiling(de/max(DE_SUB_ABS, &
                 DE_SUB_REL*0.5_wp*(edges(b) + edges(b+1))) - 1.0e-9_wp)))
    end do
    n_egrid = sum(nsub)
    if (allocated(egrid)) deallocate(egrid, egrid_de, egrid_band)
    allocate(egrid(n_egrid), egrid_de(n_egrid), egrid_band(n_egrid))
    k = 0
    do b = 1, nb
       de = (edges(b+1) - edges(b))/real(nsub(b), wp)
       do j = 1, nsub(b)
          k = k + 1
          egrid(k) = edges(b) + (real(j, wp) - 0.5_wp)*de
          egrid_de(k) = de
          egrid_band(k) = b
       end do
    end do
  end subroutine grain_energy_grid_from_edges

  !=========================================================================
  ! Photon rates per unit area on the energy grid [photons cm^-2 s^-1],
  ! n_k = 4 pi J_nu dnu_k / (h nu_k), from the band-integrated field
  ! fband(b) = Int_bin 4 pi J_nu dnu [erg cm^-2 s^-1], J_nu constant in nu
  ! inside a band bin.
  !=========================================================================
  subroutine photon_rates_from_band(fband, edges, nphot)
    real(kind=wp), intent(in)  :: fband(:), edges(:)
    real(kind=wp), intent(out) :: nphot(:)
    integer :: k, b
    do k = 1, n_egrid
       b = egrid_band(k)
       nphot(k) = fband(b)*egrid_de(k)/(edges(b+1) - edges(b)) &
                  /(egrid(k)*ev2erg)
    end do
  end subroutine photon_rates_from_band

  !=========================================================================
  ! Size bins of a SEDust model on the current energy grid.  data_dir is the
  ! SEDust data root (par%sed_data_dir), where the dielectric functions are.
  !=========================================================================
  subroutine grain_bins_from_model(m, data_dir)
    use dust_lib,  only : dust_model_t
    use utility,   only : fatal_error
    implicit none
    type(dust_model_t), intent(in) :: m
    character(len=*),   intent(in) :: data_dir
    real(kind=wp), allocatable :: lam(:), imn_gra(:), imn_sil(:), imn_ad(:)
    type(grain_bin_type), allocatable :: tmp(:)
    logical, allocatable :: used(:)
    integer :: ip, jp, npop, nbin
    logical :: paired

    npop = size(m%pops)
    !--- the model wavelength grid, ascending
    allocate(lam(m%NLAM))
    lam = m%lam
    !--- Im(m) on the energy grid for the three dielectric functions
    allocate(imn_gra(n_egrid), imn_sil(n_egrid), imn_ad(n_egrid))
    call graphite_imaginary_index(data_dir, imn_gra)
    call imaginary_index(trim(data_dir)//'/dielectric/index_silD03', imn_sil)
    if (trim(m%name) == 'astrodust') then
       call imaginary_index(trim(data_dir)// &
            '/dielectric/index_DH21Ad_P0.20_0.00_1.400', imn_ad)
    else
       imn_ad = imn_sil
    end if

    allocate(used(npop))
    used = .false.
    allocate(tmp(0))
    do ip = 1, npop
       if (used(ip)) cycle
       used(ip) = .true.
       !--- a PAH population followed by another on the same radii is its
       !--- cation counterpart (astrodust, DL07: neutral then cation)
       paired = .false.
       jp = 0
       if (trim(m%pops(ip)%grain_type) == 'pah') then
          do jp = ip+1, npop
             if (used(jp)) cycle
             if (trim(m%pops(jp)%grain_type) /= 'pah') cycle
             if (size(m%pops(jp)%aeff) /= size(m%pops(ip)%aeff)) cycle
             if (all(m%pops(jp)%aeff == m%pops(ip)%aeff)) then
                paired = .true.;  used(jp) = .true.
                exit
             end if
          end do
       end if
       select case (trim(m%pops(ip)%grain_type))
       case ('pah', 'gra')
          if (paired) then
             call merge_population(m%pops(ip)%aeff, m%pops(ip)%dn + m%pops(jp)%dn, &
                  m%pops(ip)%dn, m%pops(ip)%Cabs, m%pops(jp)%dn, m%pops(jp)%Cabs, &
                  lam, CARBONACEOUS, trim(m%pops(ip)%grain_type) == 'pah', &
                  imn_gra, tmp)
          else
             call merge_population(m%pops(ip)%aeff, m%pops(ip)%dn, &
                  m%pops(ip)%dn, m%pops(ip)%Cabs, m%pops(ip)%dn, m%pops(ip)%Cabs, &
                  lam, CARBONACEOUS, trim(m%pops(ip)%grain_type) == 'pah', &
                  imn_gra, tmp)
          end if
       case ('sil')
          call merge_population(m%pops(ip)%aeff, m%pops(ip)%dn, &
               m%pops(ip)%dn, m%pops(ip)%Cabs, m%pops(ip)%dn, m%pops(ip)%Cabs, &
               lam, SILICATE, .false., imn_ad, tmp)
       case default
          call fatal_error('grain_size_bins_mod: population type '''// &
               trim(m%pops(ip)%grain_type)//''' has no photoemission material')
       end select
    end do
    nbin = size(tmp)
    if (allocated(gbin)) deallocate(gbin)
    call move_alloc(tmp, gbin)
    n_gbin = nbin
  end subroutine grain_bins_from_model

  !=========================================================================
  ! Merge one population (radii aeff [um], number per H dn_tot) into bins of
  ! grain_bins_per_decade per decade; the neutral and the charged cross sections are
  ! the dn_neu- and dn_chg-weighted means of cabs_neu and cabs_chg [cm^2].
  !=========================================================================
  subroutine merge_population(aeff, dn_tot, dn_neu, cabs_neu, dn_chg, cabs_chg, &
                              lam, material, is_pah, imn, bins)
    real(kind=wp), intent(in) :: aeff(:), dn_tot(:), dn_neu(:), dn_chg(:)
    real(kind=wp), intent(in) :: cabs_neu(:,:), cabs_chg(:,:), lam(:)
    integer,       intent(in) :: material
    logical,       intent(in) :: is_pah
    real(kind=wp), intent(in) :: imn(:)
    type(grain_bin_type), allocatable, intent(inout) :: bins(:)
    type(grain_bin_type), allocatable :: more(:)
    real(kind=wp), allocatable :: cn(:), cc(:)
    real(kind=wp) :: lo, hi, la, wsum, a2sum, wneu, wchg, lnlo
    integer :: na, nb, ib, ia, k, nnew, ibin
    integer, allocatable :: owner(:)

    na = size(aeff)
    if (.not. any(dn_tot > 0.0_wp)) return
    lo = minval(aeff, mask=dn_tot > 0.0_wp)
    hi = maxval(aeff, mask=dn_tot > 0.0_wp)
    lnlo = log10(lo)
    nb = max(1, ceiling((log10(hi) - lnlo)*real(grain_bins_per_decade, wp) + 1.0e-9_wp))
    allocate(owner(na))
    do ia = 1, na
       owner(ia) = min(nb, 1 + int((log10(aeff(ia)) - lnlo)*real(grain_bins_per_decade, wp)))
    end do
    allocate(more(nb), cn(n_egrid), cc(n_egrid))
    nnew = 0
    do ib = 1, nb
       wsum = 0.0_wp;  a2sum = 0.0_wp;  wneu = 0.0_wp;  wchg = 0.0_wp
       cn = 0.0_wp;  cc = 0.0_wp
       do ia = 1, na
          if (owner(ia) /= ib .or. dn_tot(ia) <= 0.0_wp) cycle
          wsum  = wsum + dn_tot(ia)
          a2sum = a2sum + dn_tot(ia)*aeff(ia)**2
          if (dn_neu(ia) > 0.0_wp) then
             wneu = wneu + dn_neu(ia)
             do k = 1, n_egrid
                cn(k) = cn(k) + dn_neu(ia)*cabs_at(lam, cabs_neu(:,ia), egrid(k))
             end do
          end if
          if (dn_chg(ia) > 0.0_wp) then
             wchg = wchg + dn_chg(ia)
             do k = 1, n_egrid
                cc(k) = cc(k) + dn_chg(ia)*cabs_at(lam, cabs_chg(:,ia), egrid(k))
             end do
          end if
       end do
       if (wsum <= 0.0_wp) cycle
       nnew = nnew + 1
       ibin = nnew
       more(ibin)%material = material
       more(ibin)%pah = is_pah
       more(ibin)%n_per_H = wsum
       more(ibin)%a_cm = sqrt(a2sum/wsum)*1.0e-4_wp
       more(ibin)%a_ang = more(ibin)%a_cm*1.0e8_wp
       allocate(more(ibin)%cabs_neutral(n_egrid), more(ibin)%cabs_charged(n_egrid), &
                more(ibin)%y1(n_egrid))
       if (wneu > 0.0_wp) then
          more(ibin)%cabs_neutral = cn/wneu
       else
          more(ibin)%cabs_neutral = cc/wchg
       end if
       if (wchg > 0.0_wp) then
          more(ibin)%cabs_charged = cc/wchg
       else
          more(ibin)%cabs_charged = more(ibin)%cabs_neutral
       end if
       do k = 1, n_egrid
          la = HC_EV_CM/egrid(k)/(fourpi*max(imn(k), 1.0e-30_wp))
          more(ibin)%y1(k) = yield_enhancement(more(ibin)%a_cm, la)
       end do
    end do
    bins = [bins, more(1:nnew)]
  end subroutine merge_population

  !=========================================================================
  ! WD01 eq. 13: y_1 = (beta/alpha)^2 g(alpha)/g(beta),
  ! g(x) = x^2 - 2x + 2 - 2 exp(-x), beta = a/l_a, alpha = a/l_a + a/l_e;
  ! g by its series where the closed form cancels.
  !=========================================================================
  pure real(kind=wp) function yield_enhancement(a, la) result(y1)
    real(kind=wp), intent(in) :: a, la
    real(kind=wp) :: alpha, beta
    beta  = a/la
    alpha = beta + a/L_ESCAPE
    y1 = (beta/alpha)**2*gfun(alpha)/gfun(beta)
  contains
    pure real(kind=wp) function gfun(x) result(g)
      real(kind=wp), intent(in) :: x
      if (x < 1.0e-2_wp) then
         g = x**3/3.0_wp - x**4/12.0_wp + x**5/60.0_wp - x**6/360.0_wp
      else
         g = x*x - 2.0_wp*x + 2.0_wp - 2.0_wp*exp(-x)
      end if
    end function gfun
  end function yield_enhancement

  !=========================================================================
  ! C_abs at photon energy E [eV] from a population column on the model's
  ! ascending wavelength grid [um]: log-log in lambda.
  !=========================================================================
  real(kind=wp) function cabs_at(lam, c, e) result(v)
    real(kind=wp), intent(in) :: lam(:), c(:), e
    real(kind=wp) :: l, w
    integer :: lo, hi, mid, n
    n = size(lam)
    l = 1.0e4_wp*HC_EV_CM/e
    if (l <= lam(1)) then
       v = c(1);  return
    else if (l >= lam(n)) then
       v = c(n);  return
    end if
    lo = 1;  hi = n
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       if (lam(mid) <= l) then
          lo = mid
       else
          hi = mid
       end if
    end do
    w = log(l/lam(lo))/log(lam(hi)/lam(lo))
    v = exp((1.0_wp - w)*log(max(c(lo), tinest)) + w*log(max(c(hi), tinest)))
  end function cabs_at

  !=========================================================================
  ! Im(m) on the energy grid from a SEDust index file (two header lines,
  ! then E [eV], Re(m)-1, Im(m), ...), log-log in E.
  !=========================================================================
  subroutine imaginary_index(path, imn)
    use utility, only : fatal_error
    character(len=*), intent(in)  :: path
    real(kind=wp),    intent(out) :: imn(:)
    real(kind=wp), allocatable :: e(:), v(:)
    real(kind=wp) :: x(3), w
    character(len=512) :: line
    integer :: unit, ios, n, k, lo, hi, mid

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal_error('cannot open the dielectric function '//path)
    read(unit,'(a)') line;  read(unit,'(a)') line
    n = 0
    do
       read(unit,*,iostat=ios) x
       if (ios /= 0) exit
       n = n + 1
    end do
    rewind(unit)
    read(unit,'(a)') line;  read(unit,'(a)') line
    allocate(e(n), v(n))
    do k = 1, n
       read(unit,*) x
       e(k) = x(1);  v(k) = x(3)
    end do
    close(unit)
    if (e(n) < e(1)) then
       e = e(n:1:-1);  v = v(n:1:-1)
    end if
    do k = 1, n_egrid
       if (egrid(k) <= e(1)) then
          imn(k) = v(1)
       else if (egrid(k) >= e(n)) then
          imn(k) = v(n)
       else
          lo = 1;  hi = n
          do while (hi - lo > 1)
             mid = (lo + hi)/2
             if (e(mid) <= egrid(k)) then
                lo = mid
             else
                hi = mid
             end if
          end do
          w = log(egrid(k)/e(lo))/log(e(hi)/e(lo))
          imn(k) = exp((1.0_wp - w)*log(max(v(lo), tinest)) &
                       + w*log(max(v(hi), tinest)))
       end if
    end do
  end subroutine imaginary_index

  !--- graphite, WD01 eq. 15: Im(m)_eff = (2/3) Im(m_perp) + (1/3) Im(m_par)
  subroutine graphite_imaginary_index(data_dir, imn)
    character(len=*), intent(in)  :: data_dir
    real(kind=wp),    intent(out) :: imn(:)
    real(kind=wp) :: perp(size(imn)), para(size(imn))
    call imaginary_index(trim(data_dir)//'/dielectric/index_CpeD03', perp)
    call imaginary_index(trim(data_dir)//'/dielectric/index_CpaD03', para)
    imn = (2.0_wp*perp + para)/3.0_wp
  end subroutine graphite_imaginary_index

  !=========================================================================
  ! Grains per H in leaf il relative to the model: scale_other for every bin,
  ! scale_pah for the PAH bins (see the header).
  !=========================================================================
  subroutine grain_leaf_scale(il, nH, xHI, xHII, scale_other, scale_pah)
    use octree_mod, only : amr_grid
    integer,       intent(in)  :: il
    real(kind=wp), intent(in)  :: nH, xHI, xHII
    real(kind=wp), intent(out) :: scale_other, scale_pah
    real(kind=wp) :: wd

    scale_other = 0.0_wp;  scale_pah = 0.0_wp
    if (nH <= 0.0_wp .or. cext_ref_model <= 0.0_wp) return
    scale_other = (amr_grid%rhokap(il)/par%distance2cm)/(nH*cext_ref_model)
    !--- PAH survival applies always (author's decision, 2026-10-11): PAHs
    !--- are destroyed inside H II regions (their emission is absent from
    !--- ionized gas; WD01 Sec. 7.1 makes the same caveat), and with them
    !--- present in ionized gas the photoelectric heating rose by up to
    !--- 0.7 of the photoionization heating (pdr_L5, +1584 K).  The PAH
    !--- bins carry (x_HI + f_ion_pah x_HII), f_ion_pah = 0 by default,
    !--- independently of par%sed_pah_live (the emission switch).
    wd = 1.0_wp
    if (trim(par%dust_density_law) == 'laursen09_live') &
       wd = max(xHI + par%f_ion_dust*xHII, tinest)
    scale_pah = scale_other*(xHI + par%f_ion_pah*xHII)/wd
  end subroutine grain_leaf_scale

end module grain_size_bins_mod
