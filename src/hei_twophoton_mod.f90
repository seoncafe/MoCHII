module hei_twophoton_mod
!---------------------------------------------------------------------------
! MoCHII: He I 2^1S -> 1^1S two-photon continuum.
!
! The 2^1S term cannot decay to the ground state by a single electric-dipole
! photon; it emits a photon PAIR whose energies sum to E_tot = 20.62 eV.  The
! spectral distribution A(y) of one photon carrying the fraction
! y = E/E_tot is tabulated in data/HeI2phot.dat from Drake, Victor & Dalgarno
! (1969), symmetric about y = 1/2 and normalized so that
!
!     Int_0^1 A(y) dy = 2 A_2gamma ,     A_2gamma = 51.3 s^-1
!
! i.e. two photons per decay.  The table integrates to 1.9945 A_2gamma, which
! is the accuracy of the tabulation itself.
!
! Only the part above the H I edge matters for the diffuse ionizing field, and
! the split is not a free constant: with y_th = 13.598/20.62 = 0.65946,
!
!     photons per decay above the edge = 2 Int_{y_th}^1 A dy / Int_0^1 A dy
!                                      = 0.5564
!     in-band energy per decay         = 8.9589 eV
!     mean in-band photon energy       = 16.100 eV .
!
! diffuse_mod previously carried 0.56 and sampled the in-band photons FLAT in
! [13.598, 20.62] eV, whose mean is 17.109 eV.  The photon count was right to
! 0.6%, but a flat distribution has no memory of the peak at y = 1/2 and so
! overestimates the mean in-band energy by 6%, and the radiated in-band energy
! per decay by 6.5%.  The same A(y) table has been driving the He I
! two-photon component of nebcont_mod all along, so the nebular continuum
! MoCHII wrote and the packets it transported disagreed about one process.
!
! Sampling weight: MoCHII emits equal-energy packets whose number follows an
! ENERGY luminosity, so the energies are drawn from E n(E) and the luminosity
! uses the number-weighted mean, exactly as in milne_recomb_spectrum_mod.  The
! represented photon rate is then the physical one.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: hei_2ph_setup, hei_2ph_sample_energy
  public :: hei_2ph_photons_per_decay, hei_2ph_energy_per_decay
  public :: hei_2ph_is_ready

  !--- 2^1S - 1^1S transition energy [eV]
  real(kind=wp), parameter :: E_TOT = 20.62_wp
  integer,       parameter :: MAX_PTS = 128    ! table capacity
  integer,       parameter :: NCDF  = 1024   ! inversion grid above the edge

  integer,       save :: n_pts = 0
  real(kind=wp), save :: tab_y(MAX_PTS), tab_A(MAX_PTS)
  real(kind=wp), save :: cdf_e(NCDF), grid_y(NCDF)   ! energy-weighted CDF
  real(kind=wp), save :: n_per_decay = 0.0_wp        ! in-band photons per decay
  real(kind=wp), save :: e_per_decay = 0.0_wp        ! in-band energy per decay [eV]
  real(kind=wp), save :: y_lo = 0.0_wp, y_hi = 1.0_wp
  logical,       save :: ready = .false.

contains

  !=========================================================================
  logical function hei_2ph_is_ready() result(ok)
    ok = ready
  end function hei_2ph_is_ready

  !--- H-ionizing photons emitted per 2^1S decay (the quantity the fixed
  !--- constant 0.56 stood for).
  real(kind=wp) function hei_2ph_photons_per_decay() result(n)
    n = n_per_decay
  end function hei_2ph_photons_per_decay

  !--- In-band energy radiated per 2^1S decay [eV].  Multiplying by the 2^1S
  !--- branch rate gives that branch's contribution to the channel luminosity,
  !--- so it must be the mean of the same distribution the sampler draws from.
  real(kind=wp) function hei_2ph_energy_per_decay() result(e)
    e = e_per_decay
  end function hei_2ph_energy_per_decay

  !=========================================================================
  ! A(y) by linear interpolation of the tabulated shape; zero outside [0,1].
  !=========================================================================
  real(kind=wp) function shape_A(y) result(a)
    real(kind=wp), intent(in) :: y
    integer :: i
    a = 0.0_wp
    if (y <= tab_y(1) .or. y >= tab_y(n_pts)) return
    i = 1
    do while (i < n_pts-1 .and. tab_y(i+1) < y)
       i = i + 1
    end do
    a = tab_A(i) + (tab_A(i+1) - tab_A(i)) &
        *(y - tab_y(i))/(tab_y(i+1) - tab_y(i))
  end function shape_A

  !=========================================================================
  ! Read the shape and build the in-band moments and the sampling CDF.  Call
  ! once, after par%data_dir and par%eion_max are fixed.
  !=========================================================================
  subroutine hei_2ph_setup()
    use utility, only : fatal_error
    implicit none
    real(kind=wp) :: ya, yb, ym, fa, fb, fm, tot, acc_n, acc_e, seg
    real(kind=wp) :: e_hi
    integer :: unit, ios, i

    open(newunit=unit, file=trim(par%data_dir)//'/HeI2phot.dat', &
         status='old', iostat=ios)
    if (ios /= 0) call fatal_error('cannot open ' // &
       trim(par%data_dir)//'/HeI2phot.dat (par%data_dir)')
    n_pts = 0
    do
       read(unit,*,iostat=ios) ya, fa
       if (ios /= 0) exit
       n_pts = n_pts + 1
       if (n_pts > MAX_PTS) call fatal_error('HeI2phot.dat has more rows than ' // &
          'hei_twophoton_mod expects')
       tab_y(n_pts) = ya;  tab_A(n_pts) = fa
    end do
    close(unit)
    if (n_pts < 4) call fatal_error('HeI2phot.dat is too short to interpolate')

    !--- total over the full pair spectrum, the normalization of the shape
    tot = 0.0_wp
    do i = 1, n_pts-1
       ya = tab_y(i);  yb = tab_y(i+1);  ym = 0.5_wp*(ya + yb)
       tot = tot + (yb - ya)/6.0_wp &
             *(tab_A(i) + 4.0_wp*shape_A(ym) + tab_A(i+1))
    end do
    if (tot <= 0.0_wp) call fatal_error('HeI2phot.dat integrates to zero')

    !--- the ionizing window, truncated at the band edge if that is lower
    e_hi = min(E_TOT, par%eion_max)
    y_lo = eth_HI/E_TOT
    y_hi = e_hi/E_TOT
    if (y_hi <= y_lo) then
       n_per_decay = 0.0_wp;  e_per_decay = 0.0_wp
       cdf_e = 1.0_wp;  ready = .true.
       return
    end if

    !--- moments over the window; the factor two is the pair, and dividing by
    !--- tot turns the tabulated shape into a probability per photon
    acc_n = 0.0_wp;  acc_e = 0.0_wp
    do i = 1, NCDF-1
       ya = y_lo + (y_hi - y_lo)*real(i-1, wp)/real(NCDF-1, wp)
       yb = y_lo + (y_hi - y_lo)*real(i,   wp)/real(NCDF-1, wp)
       ym = 0.5_wp*(ya + yb)
       fa = shape_A(ya);  fm = shape_A(ym);  fb = shape_A(yb)
       acc_n = acc_n + (yb - ya)/6.0_wp*(fa + 4.0_wp*fm + fb)
       acc_e = acc_e + (yb - ya)/6.0_wp &
               *(fa*ya + 4.0_wp*fm*ym + fb*yb)*E_TOT
    end do
    n_per_decay = 2.0_wp*acc_n/tot
    e_per_decay = 2.0_wp*acc_e/tot

    !--- sampling CDF of the ENERGY-luminosity density y A(y)
    grid_y(1) = y_lo;  cdf_e(1) = 0.0_wp
    seg = 0.0_wp
    do i = 1, NCDF-1
       ya = y_lo + (y_hi - y_lo)*real(i-1, wp)/real(NCDF-1, wp)
       yb = y_lo + (y_hi - y_lo)*real(i,   wp)/real(NCDF-1, wp)
       ym = 0.5_wp*(ya + yb)
       fa = shape_A(ya)*ya;  fm = shape_A(ym)*ym;  fb = shape_A(yb)*yb
       seg = seg + (yb - ya)/6.0_wp*(fa + 4.0_wp*fm + fb)
       grid_y(i+1) = yb
       cdf_e(i+1)  = seg
    end do
    if (seg > 0.0_wp) then
       cdf_e = cdf_e/seg
    else
       cdf_e = 1.0_wp
    end if
    cdf_e(NCDF) = 1.0_wp
    ready = .true.

    if (mpar%p_rank == 0) then
       write(*,'(a,i0,a)') ' HE2PH: He I 2^1S two-photon shape loaded (', &
          n_pts, ' points, Drake et al. 1969)'
       write(*,'(a,f8.5,a,f9.5,a,f9.5,a)') &
          '        in-band: ', n_per_decay, ' photons/decay, ', &
          e_per_decay, ' eV/decay, mean ', &
          e_per_decay/max(n_per_decay, tinest), ' eV'
    end if
  end subroutine hei_2ph_setup

  !=========================================================================
  ! Energy [eV] of one H-ionizing two-photon packet, from the uniform u.
  ! Drawn from the energy-luminosity density, so that the emitted photon rate
  ! matches the branch rate when the packet count comes from the luminosity.
  !=========================================================================
  real(kind=wp) function hei_2ph_sample_energy(u) result(E)
    real(kind=wp), intent(in) :: u
    real(kind=wp) :: uu, y
    integer :: lo, hi, mid
    uu = min(max(u, 0.0_wp), 1.0_wp)
    lo = 1;  hi = NCDF
    do while (hi - lo > 1)
       mid = (lo + hi)/2
       if (uu > cdf_e(mid)) then
          lo = mid
       else
          hi = mid
       end if
    end do
    if (cdf_e(hi) > cdf_e(lo)) then
       y = grid_y(lo) + (grid_y(hi) - grid_y(lo)) &
           *(uu - cdf_e(lo))/(cdf_e(hi) - cdf_e(lo))
    else
       y = grid_y(lo)
    end if
    E = y*E_TOT
  end function hei_2ph_sample_energy

end module hei_twophoton_mod
