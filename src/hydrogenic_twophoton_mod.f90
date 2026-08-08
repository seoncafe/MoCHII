module hydrogenic_twophoton_mod
!---------------------------------------------------------------------------
! MoCHII: H-sequence 2s -> 1s two-photon continuum (H I, He II).
!
! The 2s level cannot decay to the ground state by a single electric-dipole
! photon; it emits a photon PAIR whose energies sum to E_tot (10.199 eV for
! H I, 40.813 eV for He II).  The spectral shape of one photon carrying the
! fraction y = E/E_tot is the Nussbaumer & Schmutz (1984) fit
!
!     A(y) = 202 [ y(1-y)(1 - (4y(1-y))^0.8)
!                  + 0.88 (y(1-y))^1.53 (4y(1-y))^0.8 ] ,
!
! symmetric about y = 1/2.  It is normalized here on its OWN integral, so
!
!     Int_0^1 n(y) dy = 2 ,     n(y) = 2 A(y) / Int A dy ,
!
! two photons per decay by construction.  The port carried instead a quoted
! total A_2q -- 8.2249 s^-1 for H I, 8.226 x 64 for He II -- which reproduces
! the pair count only to 0.014% and 0.05%.  A pair is two photons; that is a
! property of the decay, not a constant the shape may disagree with.
!
! Z dependence: checked against CHIANTI's hseq_2photon.dat psi tables (Goldman
! & Drake 1981) for Z = 1 and Z = 2, the normalized photon distributions agree
! to 0.42% at worst (at y = 0.98, where the tabulation ends) and to 0.2% over
! the bulk; the Z = 1 and Z = 2 CHIANTI tables themselves differ by < 0.05%.
! One shape therefore serves both ions.
!
! Band state (h2q_band_*): the He II pair spans 40.813 eV, so a large part of
! it photoionizes hydrogen and has to be transported.  h2q_band_setup builds
! the in-band moments and the sampling CDF for one E_tot; H I never needs it,
! since its pair tops out at 10.199 eV, below the Lyman limit.
!
! Sampling weight: MoCHII emits equal-energy packets whose number follows an
! ENERGY luminosity, so the energies are drawn from E n(E) and the luminosity
! uses the number-weighted mean, as in hei_twophoton_mod and
! milne_recomb_spectrum_mod.  The represented photon rate is then the
! physical one.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: h2q_photons_per_y
  public :: h2q_band_setup, h2q_band_sample_energy
  public :: h2q_band_photons_per_decay, h2q_band_energy_per_decay
  public :: h2q_band_is_ready

  integer, parameter :: NCDF = 1024        ! inversion grid over the band

  !--- Int_0^1 A(y) dy, evaluated once by h2q_norm (Simpson, NINT panels)
  integer,       parameter :: NINT  = 4096
  real(kind=wp), save :: tot_A = 0.0_wp
  logical,       save :: norm_done = .false.

  !--- band state (one E_tot: He II)
  real(kind=wp), save :: band_Etot = 0.0_wp
  real(kind=wp), save :: cdf_e(NCDF), grid_y(NCDF)
  real(kind=wp), save :: n_per_decay = 0.0_wp
  real(kind=wp), save :: e_per_decay = 0.0_wp
  logical,       save :: band_ready = .false.

contains

  !=========================================================================
  ! Nussbaumer & Schmutz (1984) shape; zero outside [0,1].
  !=========================================================================
  elemental real(kind=wp) function shape_A(y) result(a)
    real(kind=wp), intent(in) :: y
    real(kind=wp) :: p, q
    a = 0.0_wp
    if (y <= 0.0_wp .or. y >= 1.0_wp) return
    p = y*(1.0_wp - y)
    q = (4.0_wp*p)**0.8_wp
    a = 202.0_wp*(p*(1.0_wp - q) + 0.88_wp*p**1.53_wp*q)
  end function shape_A

  !--- Int_0^1 A dy by composite Simpson; the endpoints are zero.
  subroutine h2q_norm()
    real(kind=wp) :: h, s, y
    integer :: i
    if (norm_done) return
    h = 1.0_wp/real(NINT, wp)
    s = 0.0_wp
    do i = 1, NINT
       y = (real(i,wp) - 0.5_wp)*h
       s = s + h/6.0_wp*(shape_A(real(i-1,wp)*h) + 4.0_wp*shape_A(y) &
                         + shape_A(real(i,wp)*h))
    end do
    tot_A = s
    norm_done = .true.
  end subroutine h2q_norm

  !--- Photons per unit y per decay; Int_0^1 dy = 2 exactly.
  real(kind=wp) function h2q_photons_per_y(y) result(p)
    real(kind=wp), intent(in) :: y
    call h2q_norm()
    p = 2.0_wp*shape_A(y)/tot_A
  end function h2q_photons_per_y

  !=========================================================================
  logical function h2q_band_is_ready() result(ok)
    ok = band_ready
  end function h2q_band_is_ready

  !--- H-ionizing photons emitted per 2s decay.
  real(kind=wp) function h2q_band_photons_per_decay() result(n)
    n = n_per_decay
  end function h2q_band_photons_per_decay

  !--- In-band energy radiated per 2s decay [eV].  Multiplying by the branch
  !--- rate gives that branch's contribution to the channel luminosity, so it
  !--- must be the mean of the same distribution the sampler draws from.
  real(kind=wp) function h2q_band_energy_per_decay() result(e)
    e = e_per_decay
  end function h2q_band_energy_per_decay

  !=========================================================================
  ! Build the in-band moments and the sampling CDF for a pair of total energy
  ! E_tot over the window [E_lo, E_hi].  The window is the transported band:
  ! E_lo is the hydrogen threshold, E_hi the band ceiling (par%eion_max) or
  ! E_tot, whichever is lower.
  !=========================================================================
  subroutine h2q_band_setup(E_tot, E_lo, E_hi)
    real(kind=wp), intent(in) :: E_tot, E_lo, E_hi
    real(kind=wp) :: y_lo, y_hi, ya, yb, ym, fa, fm, fb, acc_n, acc_e, seg
    integer :: i

    call h2q_norm()
    band_Etot = E_tot
    y_lo = E_lo/E_tot
    y_hi = min(E_hi, E_tot)/E_tot
    if (y_hi <= y_lo) then
       n_per_decay = 0.0_wp;  e_per_decay = 0.0_wp
       grid_y = y_lo;  cdf_e = 1.0_wp;  band_ready = .true.
       return
    end if

    !--- moments over the window; the factor two is the pair, and dividing by
    !--- tot_A turns the shape into a probability per photon
    acc_n = 0.0_wp;  acc_e = 0.0_wp
    do i = 1, NCDF-1
       ya = y_lo + (y_hi - y_lo)*real(i-1, wp)/real(NCDF-1, wp)
       yb = y_lo + (y_hi - y_lo)*real(i,   wp)/real(NCDF-1, wp)
       ym = 0.5_wp*(ya + yb)
       fa = shape_A(ya);  fm = shape_A(ym);  fb = shape_A(yb)
       acc_n = acc_n + (yb - ya)/6.0_wp*(fa + 4.0_wp*fm + fb)
       acc_e = acc_e + (yb - ya)/6.0_wp &
               *(fa*ya + 4.0_wp*fm*ym + fb*yb)*E_tot
    end do
    n_per_decay = 2.0_wp*acc_n/tot_A
    e_per_decay = 2.0_wp*acc_e/tot_A

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
    band_ready = .true.
  end subroutine h2q_band_setup

  !=========================================================================
  ! Energy [eV] of one in-band two-photon packet, from the uniform u.
  !=========================================================================
  real(kind=wp) function h2q_band_sample_energy(u) result(E)
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
    E = y*band_Etot
  end function h2q_band_sample_energy

end module hydrogenic_twophoton_mod
