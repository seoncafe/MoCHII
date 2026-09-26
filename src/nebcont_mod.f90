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
!    shape, hydrogenic_twophoton_mod), fed at alpha_B f_2s with the share
!    f_2s of one zero-density Storey & Hummer (1991) cascade in the case of
!    the recombination lines (hydrogenic_recombination_mod); the
!    collisional 2s -> 2p suppression for H I (par%twophoton_quenching,
!    h2s_quench_denominator); He II 2s is not quenched, its A_2s = 526.6 s^-1
!    putting the critical density near 1e6 cm^-3.  He I 2^1S (Drake, Victor
!    & Dalgarno 1969 shape from data/HeI2phot.dat, owned by
!    hei_twophoton_mod), fed at the 2^1S decay rate the He I cascade exits
!    give (hei_recombination_photons_mod::hei_twophoton_decays).
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

  !--- H(2s) -> H(2p) collisional rate coefficients by collider (e, p, He+,
  !--- He2+) on a log T grid, and A_2s (data/atomic/h_2s_2p_mixing.txt)
  integer, parameter :: NMIX_MAX = 128
  integer :: n_mix = 0
  real(kind=wp) :: mix_logT(NMIX_MAX), mix_logq(4, NMIX_MAX), a_2s = 0.0_wp

contains

  !=========================================================================
  subroutine nebcont_setup()
    use hei_twophoton_mod, only : hei_2ph_setup, hei_2ph_is_ready
    implicit none
    integer :: i

    call recomb_continuum_setup()

    !--- He I 2^1S two-photon shape.  hei_twophoton_mod owns HeI2phot.dat --
    !--- one reader, one normalization -- and main.f90 has already built it
    !--- when the He I cascade emits into the diffuse field; build it here
    !--- otherwise, since the written continuum needs the same shape.
    if (.not. hei_2ph_is_ready()) call hei_2ph_setup()
    !--- the 2s -> 2p rate coefficients behind the H I quenching
    if (trim(par%twophoton_quenching) == 'colliders') call h2s_mixing_load()

    !--- output grid 0.008 - 4.2 Ryd (lambda ~ 11.4 um - 217 A)
    do i = 1, NNU_OUT
       nu_out(i) = exp(log(0.008_wp) + (log(4.2_wp) - log(0.008_wp)) &
                   *real(i-1,wp)/real(NNU_OUT-1,wp))
    end do
    if (mpar%p_rank == 0) write(*,'(a)') ' NEBC: continuum tables loaded'
  end subroutine nebcont_setup

  !=========================================================================
  ! Read data/atomic/h_2s_2p_mixing.txt (tools/fitting/
  ! make_h2s_collisional_mixing.py): A_2S and the table T, q_e, q_p, q_He+,
  ! q_He2+ [cm^3 s^-1].
  !=========================================================================
  subroutine h2s_mixing_load()
    use utility, only : fatal_error
    character(len=256) :: line
    character(len=16)  :: key
    real(kind=wp) :: row(5)
    integer :: unit, ios, i
    open(newunit=unit, file=trim(par%atomic_dir)//'/h_2s_2p_mixing.txt', &
         status='old', iostat=ios)
    if (ios /= 0) call fatal_error('cannot open '//trim(par%atomic_dir)// &
       '/h_2s_2p_mixing.txt (par%twophoton_quenching = ''colliders'')')
    n_mix = 0;  a_2s = 0.0_wp
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
       read(line,*) key
       if (trim(key) == 'A_2S') then
          read(line,*) key, a_2s
       else if (trim(key) == 'TABLE') then
          read(line,*) key, n_mix
          if (n_mix < 2 .or. n_mix > NMIX_MAX) &
             call fatal_error('h_2s_2p_mixing.txt: TABLE size out of range')
          do i = 1, n_mix
             read(unit,*) row
             mix_logT(i)    = log(row(1))
             mix_logq(:, i) = log(row(2:5))
          end do
       end if
    end do
    close(unit)
    if (n_mix < 2 .or. a_2s <= 0.0_wp) &
       call fatal_error('h_2s_2p_mixing.txt is missing A_2S or TABLE')
  end subroutine h2s_mixing_load

  !=========================================================================
  ! Denominator d of the surviving H I two-photon fraction 1/d at a leaf:
  ! the 2s level decays by two photons at A_2s and is moved to 2p (then
  ! Ly-alpha) by collisions.
  !  'colliders': d = 1 + (n_e q_e + n_p q_p + n(He+) q_He+ + n(He2+) q_He2+)/A_2s,
  !     q by log-log interpolation in T, held at the table ends; the table is
  !     the low-density limit (the Debye cut-off takes over only above
  !     n_e ~ 2e12 cm^-3).
  !  'lumped_ne': d = 1 + n_e q/8.23, q = 5.31e-4 at T <= 1e4 K, 4.71e-4 at
  !     T >= 2e4 K, power law between (the Osterbrock Table 4.10 totals as
  !     MOCASSIN carries them).  The totals are the proton + electron sums,
  !     multiplied here by n_e alone.
  !=========================================================================
  real(kind=wp) function h2s_quench_denominator(T, ne, nHII, nHeII, nHeIII) &
                                                                  result(d)
    real(kind=wp), intent(in) :: T, ne, nHII, nHeII, nHeIII
    real(kind=wp) :: q, fac, lt, w, qk(4)
    integer :: i
    if (trim(par%twophoton_quenching) == 'lumped_ne') then
       if (T <= 1.0e4_wp) then
          q = 5.31e-4_wp
       else if (T >= 2.0e4_wp) then
          q = 4.71e-4_wp
       else
          fac = log10(4.71e-4_wp/5.31e-4_wp)/log10(2.0_wp)
          q = 10.0_wp**(log10(5.31e-4_wp) + fac*log10(T/1.0e4_wp))
       end if
       d = 1.0_wp + ne*q/8.23_wp
       return
    end if
    lt = min(max(log(T), mix_logT(1)), mix_logT(n_mix))
    i = 1
    do while (i < n_mix - 1 .and. mix_logT(i+1) < lt)
       i = i + 1
    end do
    w  = (lt - mix_logT(i))/(mix_logT(i+1) - mix_logT(i))
    qk = exp((1.0_wp - w)*mix_logq(:, i) + w*mix_logq(:, i+1))
    d  = 1.0_wp + (ne*qk(1) + nHII*qk(2) + nHeII*qk(3) + nHeIII*qk(4))/a_2s
  end function h2s_quench_denominator

  !=========================================================================
  ! Two-photon coefficients (1e-40 erg cm^3/s/Hz), MOCASSIN formulas.
  !=========================================================================
  !--- H I two-photon; dq is h2s_quench_denominator at the leaf.  The 2s
  !--- population is fed at alpha_B f_2s: MoCHII's own alpha_B, the rate the
  !--- ionization balance uses, times the share of the n >= 2 captures that
  !--- end in 2s (hydrogenic_recombination_mod), in the case of the
  !--- recombination lines (lyman_line_case) -- case A when the Lyman lines
  !--- escape (par%lyman_lines = 'escaping'), case B otherwise.  This took
  !--- the case-B Pengelly (1964) power law 0.8368e-13 T4^-0.723 whatever the
  !--- lines did, so a case-A run wrote case-A Balmer lines beside a
  !--- two-photon continuum 2.63x its own case at 1e4 K.
  real(kind=wp) function gamma_2q_HI(nu, T, dq, lyman_case) result(g)
    use hydrogenic_twophoton_mod,     only : h2q_photons_per_y
    use hydrogenic_recombination_mod, only : hydrogenic_f_2s
    use recomb_mod,                   only : alphaB_HII
    implicit none
    real(kind=wp),    intent(in) :: nu, T, dq
    character(len=1), intent(in) :: lyman_case
    real(kind=wp) :: y, gNu, a2s
    real(kind=wp), parameter :: nu0 = 0.7496_wp
    g = 0.0_wp
    y = nu/nu0
    if (y >= 1.0_wp) return
    a2s = alphaB_HII(T)*hydrogenic_f_2s(1, T, lyman_case)   ! [cm^3/s]
    gNu = h_planck_cgs*y*h2q_photons_per_y(y)               ! [erg/Hz]
    g = a2s*gNu/dq * 1.0e40_wp
  end function gamma_2q_HI

  !--- He II two-photon, fed at alphaB_HeIII f_2s in the case of the
  !--- recombination lines, as for H I.  In case B this is the rate the He III
  !--- cascade of the diffuse field sends through 2^2S (heii_cascade_mod, the
  !--- same alpha_B and the same f_2s), so the written continuum and the
  !--- transported packets describe one nebula.  It used to take SH95's
  !--- alpha_eff(2^2S) at n_e = 100 on four anchors, 5000-30000 K, held at the
  !--- end values outside them, and case B only.  The shape is normalized on
  !--- its own integral (two photons per decay).
  real(kind=wp) function gamma_2q_HeII(nu, T, lyman_case) result(g)
    use hydrogenic_twophoton_mod,     only : h2q_photons_per_y
    use hydrogenic_recombination_mod, only : hydrogenic_f_2s
    use recomb_mod,                   only : alphaB_HeIII
    implicit none
    real(kind=wp),    intent(in) :: nu, T
    character(len=1), intent(in) :: lyman_case
    real(kind=wp) :: y, gNu, a2s
    real(kind=wp), parameter :: nu0 = 3.00_wp
    g = 0.0_wp
    y = nu/nu0
    if (y >= 1.0_wp) return
    a2s = alphaB_HeIII(T)*hydrogenic_f_2s(2, T, lyman_case) ! [cm^3/s]
    gNu = h_planck_cgs*y*h2q_photons_per_y(y)               ! [erg/Hz]
    g = a2s*gNu * 1.0e40_wp
  end function gamma_2q_HeII

  !--- He I two-photon continuum, per 2^1S decay rate r2q [cm^-3 s^-1]
  !--- (hei_recombination_photons_mod::hei_twophoton_decays): the population
  !--- is fed by recombination, by the trapped 584 A line that decays to
  !--- 2^1S (the 2.06 micron 2^1P -> 2^1S branch) rather than photoionizing
  !--- hydrogen, and, with the 2^3S population coupled, by its collisional
  !--- transfers into the singlets.  One routine gives those exits to the
  !--- diffuse channel 4, the on-the-spot sources and this written
  !--- continuum, so the three describe one nebula.  It used to take
  !--- alpha(2^1S) + [1 - p_584] alpha(2^1P) from the BSS99 + JSM fits in
  !--- absolute units, while the cascade normalizes the same fits by their
  !--- sum and scales by alpha_B(He II): the two differed by that sum,
  !--- 0.91 alpha_B at 8000 K and 0.985 at 1e4 K, and the coupled transfers
  !--- were missing (x1.09 at n_e = 100, x3.5 at 1e4 in the 2^1S energy).
  !--- The earlier port took MOCASSIN's 6.23e-14 (T4)^-0.827, the singlet
  !--- TOTAL (every singlet ends in 2^1S, the p_584 -> 0 limit).
  !---
  !--- The shape comes from hei_twophoton_mod, normalized on the table's own
  !--- integral so that the decay emits exactly two photons.  The port divided
  !--- the same table by a separately quoted A = 51.3 s^-1, but the tabulation
  !--- carries A = 50.93 (it is Drake, Victor & Dalgarno 1969 psi times that
  !--- A value, as CHIANTI's heseq_2photon.dat shows with A = 50.94), so the
  !--- written continuum radiated 1.9945 photons per decay while the packets
  !--- transported 2.  A pair is two photons; the normalization is not a
  !--- constant the shape can disagree with.  Result: r2q h y n(y) [erg/s/Hz
  !--- per cm^3], in the 1e-40 units of the other emissivities times n_e n_i.
  real(kind=wp) function gamma_2q_HeI(nu, r2q) result(g)
    use hei_twophoton_mod, only : hei_2ph_photons_per_y
    implicit none
    real(kind=wp), intent(in) :: nu, r2q
    real(kind=wp) :: y
    real(kind=wp), parameter :: nu0 = 1.514_wp
    g = 0.0_wp
    y = nu/nu0
    if (y >= 1.0_wp) return
    g = r2q*h_planck_cgs*y*hei_2ph_photons_per_y(y)*1.0e40_wp
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
    use hei_recombination_photons_mod, only : hei_twophoton_decays
    use hydrogenic_recombination_mod, only : lyman_line_case
    implicit none
    real(kind=wp) :: Lfb(NNU_OUT), Lff(NNU_OUT), L2q(NNU_OUT)
    real(kind=wp) :: T, ne, nH, vol, nHII, nHeII_n, nHeIII, xHeIII, dq, r2q_HeI
    real(kind=wp) :: nu, e_eV, gH, gHe1, gHe2, ff1, ff2
    character(len=PATH_LEN) :: outname
    character(len=1) :: lcase
    integer :: il, i, unit

    if (mpar%p_rank /= 0) return
    Lfb = 0.0_wp;  Lff = 0.0_wp;  L2q = 0.0_wp
    lcase = lyman_line_case()

    do il = 1, gas_nleaf
       nH = gas_nH(il)
       if (nH <= 0.0_wp) cycle
       T  = gas_Te(il);  ne = gas_ne(il)
       vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
       xHeIII  = max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il))
       nHII    = nH*(1.0_wp - gas_xHI(il))
       nHeII_n = nH*par%abund_He*gas_xHeII(il)
       nHeIII  = nH*par%abund_He*xHeIII
       dq      = h2s_quench_denominator(T, ne, nHII, nHeII_n, nHeIII)
       r2q_HeI = 0.0_wp
       if (par%abund_He > 0.0_wp) &
          r2q_HeI = hei_twophoton_decays(il, T, nH, ne, gas_xHI(il), &
                                         gas_xHeI(il), gas_xHeII(il))
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
          L2q(i) = L2q(i) + ( ne*( nHII*gamma_2q_HI(nu, T, dq, lcase) &
                   + nHeIII*gamma_2q_HeII(nu, T, lcase) ) &
                   + gamma_2q_HeI(nu, r2q_HeI) )*1.0e-40_wp*vol
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
