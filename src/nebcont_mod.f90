module nebcont_mod
!---------------------------------------------------------------------------
! MoCHII: nebular continuum emission (port of MOCASSIN emission_mod::fb_ff
! + twoPhoton/HeI2photSub).
!
! The free-bound + free-free coefficients and their provenance live in
! recomb_continuum_mod, which this module and heii_cascade_mod share so that
! the written continuum and the transported He II free-bound packets cannot
! come from different tables.  What is added here:
!  - Two-photon: H I 2s and He II 2s (Nussbaumer & Schmutz 1984 spectral
!    shape, hydrogenic_twophoton_mod; Pengelly 1964 alpha_eff(2s) with the
!    Osterbrock q(2^2S-2^2P) collisional suppression for H I, Storey &
!    Hummer 1995 alpha_eff for He II), He I 2^1S (Drake, Victor & Dalgarno
!    1969 shape from data/HeI2phot.dat, owned by hei_twophoton_mod;
!    Benjamin, Skillman & Smits 1999 alpha_eff).
!  - The output grid and the grid-integrated write.
!
! Output (nebcont_write): '<base>_nebcont.txt' — the grid-integrated
! continuum L_nu on a log frequency grid, total and per component.
!---------------------------------------------------------------------------
  use define
  use recomb_continuum_mod, only : recomb_continuum_setup, gamma_tab, &
                                   gamma_milne, gamma_ff, &
                                   HCRYD_K, NU_HI, NU_HEI, NU_HEII
  implicit none
  private

  public :: nebcont_setup, nebcont_write

  !--- output frequency grid [Ryd]
  integer, parameter :: NNU_OUT = 400
  real(kind=wp) :: nu_out(NNU_OUT)

contains

  !=========================================================================
  subroutine nebcont_setup()
    use hei_twophoton_mod, only : hei_2ph_setup, hei_2ph_is_ready
    use heii_cascade_mod,  only : heii_cascade_setup, heii_cascade_is_ready
    implicit none
    integer :: i

    call recomb_continuum_setup()

    !--- He I 2^1S two-photon shape.  hei_twophoton_mod owns HeI2phot.dat --
    !--- one reader, one normalization -- and main.f90 has already built it
    !--- when the He I cascade emits into the diffuse field; build it here
    !--- otherwise, since the written continuum needs the same shape.
    if (.not. hei_2ph_is_ready()) call hei_2ph_setup()
    !--- likewise the SH95 alpha_eff(2^2S) behind the He II two-photon term
    if (.not. heii_cascade_is_ready()) call heii_cascade_setup()

    !--- output grid 0.008 - 4.2 Ryd (lambda ~ 11.4 um - 217 A)
    do i = 1, NNU_OUT
       nu_out(i) = exp(log(0.008_wp) + (log(4.2_wp) - log(0.008_wp)) &
                   *real(i-1,wp)/real(NNU_OUT-1,wp))
    end do
    if (mpar%p_rank == 0) write(*,'(a)') ' NEBC: continuum tables loaded'
  end subroutine nebcont_setup

  !=========================================================================
  ! Two-photon coefficients (1e-40 erg cm^3/s/Hz), MOCASSIN formulas.
  !=========================================================================
  real(kind=wp) function gamma_2q_HI(nu, T, ne) result(g)
    use hydrogenic_twophoton_mod, only : h2q_photons_per_y
    implicit none
    real(kind=wp), intent(in) :: nu, T, ne
    real(kind=wp) :: y, gNu, a2s, q, fac
    real(kind=wp), parameter :: nu0 = 0.7496_wp
    g = 0.0_wp
    y = nu/nu0
    if (y >= 1.0_wp) return
    a2s = 0.8368_wp*(T*1.0e-4_wp)**(-0.723_wp)      ! [1e-13 cm^3/s], Pengelly 64
    gNu = h_planck_cgs*y*h2q_photons_per_y(y)       ! [erg/Hz]
    if (T <= 1.0e4_wp) then
       q = 5.31e-4_wp
    else if (T >= 2.0e4_wp) then
       q = 4.71e-4_wp
    else
       fac = log10(4.71e-4_wp/5.31e-4_wp)/log10(2.0_wp)
       q = 10.0_wp**(log10(5.31e-4_wp) + fac*log10(T/1.0e4_wp))
    end if
    g = a2s*1.0e-13_wp*gNu/(1.0_wp + ne*q/8.23_wp) * 1.0e40_wp
  end function gamma_2q_HI

  !--- He II two-photon.  alpha_eff(2^2S) comes from heii_cascade_mod, the
  !--- same SH95 table that sets how much of the He III cascade leaves through
  !--- 2^2S in the diffuse field, so the written continuum and the transported
  !--- packets describe one nebula.  The shape is normalized on its own
  !--- integral (two photons per decay); the port divided it by a quoted
  !--- 8.226 x 64 s^-1, which reproduced the pair count to 0.05%.
  real(kind=wp) function gamma_2q_HeII(nu, T) result(g)
    use hydrogenic_twophoton_mod, only : h2q_photons_per_y
    use heii_cascade_mod,         only : heii_alpha_2s
    implicit none
    real(kind=wp), intent(in) :: nu, T
    real(kind=wp) :: y, gNu, a2s
    real(kind=wp), parameter :: nu0 = 3.00_wp
    g = 0.0_wp
    y = nu/nu0
    if (y >= 1.0_wp) return
    a2s = heii_alpha_2s(T)                          ! [cm^3/s], SH95 n_e = 100
    gNu = h_planck_cgs*y*h2q_photons_per_y(y)       ! [erg/Hz]
    g = a2s*gNu * 1.0e40_wp
  end function gamma_2q_HeII

  !--- He I two-photon continuum.  The emitting population is 2^1S, fed both
  !--- directly by recombination and by the trapped 584 A line that is
  !--- collisionally transferred 2^1P -> 2^1S rather than photoionizing
  !--- hydrogen.  That second path is the same branching diffuse_mod applies
  !--- (hei_584_ionizes_H), so the transported packets and this written
  !--- continuum describe one nebula:
  !---
  !---     alpha_2q = alpha(2^1S) + [1 - p_584] alpha(2^1P) .
  !---
  !--- The port took MOCASSIN's constant 6.23e-14 (T4)^-0.827 here, which is
  !--- the singlet TOTAL -- 2.06e-14 + 4.17e-14 = 6.23e-14 exactly, refitted
  !--- as one power law -- under its stated assumption that every singlet ends
  !--- in 2^1S.  That is the p_584 -> 0 limit, and it disagreed with
  !--- diffuse_mod's p_584 -> 1 limit about the fate of the same decay.
  !---
  !--- The shape comes from hei_twophoton_mod, normalized on the table's own
  !--- integral so that the decay emits exactly two photons.  The port divided
  !--- the same table by a separately quoted A = 51.3 s^-1, but the tabulation
  !--- carries A = 50.93 (it is Drake, Victor & Dalgarno 1969 psi times that
  !--- A value, as CHIANTI's heseq_2photon.dat shows with A = 50.94), so the
  !--- written continuum radiated 1.9945 photons per decay while the packets
  !--- transported 2.  A pair is two photons; the normalization is not a
  !--- constant the shape can disagree with.
  real(kind=wp) function gamma_2q_HeI(nu, T, xHI, xHeI) result(g)
    use hei_cascade_mod,   only : hei_584_ionizes_H
    use hei_twophoton_mod, only : hei_2ph_photons_per_y
    implicit none
    real(kind=wp), intent(in) :: nu, T, xHI, xHeI
    real(kind=wp) :: y, a2s
    real(kind=wp), parameter :: nu0 = 1.514_wp
    g = 0.0_wp
    y = nu/nu0
    if (y >= 1.0_wp) return
    !--- BSS99 Table 1 + JSM, in units of 1e-14 cm^3 s^-1
    a2s = 2.06_wp*(T/1.0e4_wp)**(-0.451_wp) &
          + (1.0_wp - hei_584_ionizes_H(T, xHI, xHeI)) &
            *4.17_wp*(T/1.0e4_wp)**(-0.695_wp)
    !--- result in 1e-40 units: a2s [1e-14] x (h x 1e26 = 0.66262) x y n(y)
    g = a2s*0.66262_wp*y*hei_2ph_photons_per_y(y)
  end function gamma_2q_HeI

  !=========================================================================
  ! Grid-integrated continuum spectrum: '<base>_nebcont.txt'.
  !=========================================================================
  subroutine nebcont_write()
    use octree_mod, only : amr_grid, leaf_half
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, &
                              gas_Te, gas_nleaf
    use photo_xsec,    only : sigma_HI, sigma_HeI, sigma_HeII
    use utility,       only : get_base_name
    implicit none
    real(kind=wp) :: Lfb(NNU_OUT), Lff(NNU_OUT), L2q(NNU_OUT)
    real(kind=wp) :: T, ne, nH, vol, nHII, nHeII_n, nHeIII, xHeIII
    real(kind=wp) :: nu, e_eV, gH, gHe1, gHe2, ff1, ff2
    character(len=192) :: outname
    integer :: il, i, unit

    if (mpar%p_rank /= 0) return
    Lfb = 0.0_wp;  Lff = 0.0_wp;  L2q = 0.0_wp

    do il = 1, gas_nleaf
       nH = gas_nH(il)
       if (nH <= 0.0_wp) cycle
       T  = gas_Te(il);  ne = gas_ne(il)
       vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
       xHeIII  = max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il))
       nHII    = nH*(1.0_wp - gas_xHI(il))
       nHeII_n = nH*par%abund_He*gas_xHeII(il)
       nHeIII  = nH*par%abund_He*xHeIII
       do i = 1, NNU_OUT
          nu = nu_out(i)
          e_eV = nu*13.605693_wp
          !--- table region: fb+ff combined; above: Milne fb + explicit ff
          gH   = gamma_tab(1, nu, T)
          gHe1 = gamma_tab(2, nu, T)
          gHe2 = gamma_tab(3, nu, T)
          if (nu >= NU_HI)   gH   = gamma_milne(nu, NU_HI,  sigma_HI(e_eV), &
                                    1.0_wp, T)
          if (nu >= NU_HEI)  gHe1 = gamma_milne(nu, NU_HEI, sigma_HeI(e_eV), &
                                    0.25_wp, T)
          if (nu >= NU_HEII) gHe2 = gamma_milne(nu, NU_HEII, sigma_HeII(e_eV), &
                                    1.0_wp, T)
          if (nu >= NU_HI) then
             ff1 = gamma_ff(nu, T, 1)
             gH = gH + ff1
             if (nu >= NU_HEI)  gHe1 = gHe1 + ff1
             if (nu >= NU_HEII) gHe2 = gHe2 + gamma_ff(nu, T, 2) - ff1
          end if
          Lfb(i) = Lfb(i) + ne*(nHII*gH + nHeII_n*gHe1 + nHeIII*gHe2) &
                   *1.0e-40_wp*vol
          !--- explicit ff component for the report (table region estimate)
          if (nu < NU_HI) then
             ff1 = gamma_ff(nu, T, 1)
             ff2 = gamma_ff(nu, T, 2)
             Lff(i) = Lff(i) + ne*((nHII + nHeII_n)*ff1 + nHeIII*ff2) &
                      *1.0e-40_wp*vol
          end if
          L2q(i) = L2q(i) + ne*( nHII*gamma_2q_HI(nu, T, ne) &
                   + nHeII_n*gamma_2q_HeI(nu, T, gas_xHI(il), gas_xHeI(il)) &
                   + nHeIII*gamma_2q_HeII(nu, T) )*1.0e-40_wp*vol
       end do
    end do

    outname = trim(get_base_name(par%out_file))//'_nebcont.txt'
    open(newunit=unit, file=trim(outname), status='replace')
    write(unit,'(a)') '# MoCHII nebular continuum (grid-integrated), '// &
       'fb_ff port (see nebcont_mod.f90 header for provenance)'
    write(unit,'(a)') '# col 1: nu [Ryd]; 2: lambda [um]; '// &
       '3: L_nu total = fb+ff(+Milne) + two-photon [erg/s/Hz]; '// &
       '4: ff estimate below the H edge; 5: two-photon'
    do i = 1, NNU_OUT
       write(unit,'(f10.5,f12.5,3es14.5)') nu_out(i), &
          0.0911267_wp/nu_out(i), Lfb(i) + L2q(i), Lff(i), L2q(i)
    end do
    close(unit)
    write(*,'(2a)') ' NEBC: nebular continuum written to: ', trim(outname)
  end subroutine nebcont_write

end module nebcont_mod
