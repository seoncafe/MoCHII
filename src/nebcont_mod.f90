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
  use recomb_continuum_mod, only : recomb_continuum_setup, gamma_milne, &
                                   gamma_ff, NU_HI, NU_HEI, NU_HEII
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
    use utility, only : abort_run
    character(len=256) :: line
    character(len=16)  :: key
    real(kind=wp) :: row(5)
    integer :: unit, ios, i
    open(newunit=unit, file=trim(par%atomic_dir)//'/h_2s_2p_mixing.txt', &
         status='old', iostat=ios)
    if (ios /= 0) call abort_run('cannot open '//trim(par%atomic_dir)// &
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
             call abort_run('h_2s_2p_mixing.txt: TABLE size out of range')
          do i = 1, n_mix
             read(unit,*) row
             mix_logT(i)    = log(row(1))
             mix_logq(:, i) = log(row(2:5))
          end do
       end if
    end do
    close(unit)
    if (n_mix < 2 .or. a_2s <= 0.0_wp) &
       call abort_run('h_2s_2p_mixing.txt is missing A_2S or TABLE')
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
  ! Two-photon continua.  Each is a decay rate per unit volume times the
  ! spectral energy one decay radiates, h nu_pair y n(y) [erg/Hz] with
  ! y = nu/nu_pair and n(y) the photons per unit y, normalized on the
  ! shape's own integral so that a decay emits exactly two photons.  The
  ! spectral energy depends on the frequency alone and the rate on the leaf
  ! alone, so nebcont_write forms the first once on the output grid and the
  ! second once for each leaf.
  !  - H I 2s and He II 2s: Nussbaumer & Schmutz (1984) shape
  !    (hydrogenic_twophoton_mod, one evaluation for both).  The 2s level is
  !    fed at alpha_B f_2s (twophoton_2s_feed): MoCHII's own alpha_B, the
  !    rate the ionization balance uses, times the share of the n >= 2
  !    captures that end in 2s (hydrogenic_recombination_mod), in the case
  !    of the recombination lines (lyman_line_case) -- case A when the
  !    Lyman lines escape (par%lyman_lines = 'escaping'), case B otherwise.
  !    For He II in case B this is the rate the He III cascade of the
  !    diffuse field sends through 2^2S (heii_cascade_mod, the same alpha_B
  !    and the same f_2s), so the written continuum and the transported
  !    packets describe one nebula.  H I is divided by the collisional
  !    2s -> 2p suppression h2s_quench_denominator; He II 2s is not
  !    quenched (A_2s = 526.6 s^-1 puts its critical density near 1e6 cm^-3).
  !  - He I 2^1S: Drake, Victor & Dalgarno (1969) shape (hei_twophoton_mod,
  !    which owns data/HeI2phot.dat), at the 2^1S decay rate of
  !    hei_recombination_photons_mod::hei_twophoton_decays: the population
  !    is fed by recombination, by the trapped 584 A line that decays to
  !    2^1S (the 2.06 micron 2^1P -> 2^1S branch) rather than photoionizing
  !    hydrogen, and, with the 2^3S population coupled, by its collisional
  !    transfers into the singlets.  One routine gives those exits to the
  !    diffuse channel 4, the on-the-spot sources and this written
  !    continuum, so the three describe one nebula.  The table carries
  !    A = 50.93 s^-1 (Drake, Victor & Dalgarno 1969 psi times that A value;
  !    CHIANTI's heseq_2photon.dat gives A = 50.94); normalizing on its
  !    integral, not on a separately quoted A, keeps the pair at two
  !    photons, as in the transported packets.
  !=========================================================================
  !--- spectral energy of one hydrogenic 2s decay [erg/Hz]; nu, nu_pair in
  !--- Ryd (H I 0.7496, He II 3.00); 0 at and above nu_pair
  real(kind=wp) function hydrogenic_twophoton_energy(nu, nu_pair) result(e)
    use hydrogenic_twophoton_mod, only : h2q_photons_per_y
    implicit none
    real(kind=wp), intent(in) :: nu, nu_pair
    real(kind=wp) :: y
    e = 0.0_wp
    y = nu/nu_pair
    if (y >= 1.0_wp) return
    e = h_planck_cgs*y*h2q_photons_per_y(y)
  end function hydrogenic_twophoton_energy

  !--- spectral energy of one He I 2^1S decay [erg/Hz]; 0 at and above
  !--- the pair energy 1.514 Ryd
  real(kind=wp) function hei_twophoton_energy(nu) result(e)
    use hei_twophoton_mod, only : hei_2ph_photons_per_y
    implicit none
    real(kind=wp), intent(in) :: nu
    real(kind=wp) :: y
    real(kind=wp), parameter :: nu_pair = 1.514_wp
    e = 0.0_wp
    y = nu/nu_pair
    if (y >= 1.0_wp) return
    e = h_planck_cgs*y*hei_2ph_photons_per_y(y)
  end function hei_twophoton_energy

  !--- the rate coefficient alpha_B f_2s [cm^3/s] that feeds 2s of the
  !--- hydrogenic ion of nuclear charge z (1: H I from H II, 2: He II from
  !--- He III), in the case lyman_case of the recombination lines
  real(kind=wp) function twophoton_2s_feed(z, T, lyman_case) result(a)
    use hydrogenic_recombination_mod, only : hydrogenic_f_2s
    use recomb_mod,                   only : alphaB_HII, alphaB_HeIII
    implicit none
    integer,          intent(in) :: z
    real(kind=wp),    intent(in) :: T
    character(len=1), intent(in) :: lyman_case
    if (z == 1) then
       a = alphaB_HII(T)*hydrogenic_f_2s(1, T, lyman_case)
    else
       a = alphaB_HeIII(T)*hydrogenic_f_2s(2, T, lyman_case)
    end if
  end function twophoton_2s_feed

  !=========================================================================
  ! Grid-integrated continuum spectrum: '<base>_nebcont.txt'.  Collective:
  ! every rank sums the leaves it owns in the thermal solve
  ! (gas_state_mod::gas_leaf_owner_first/stride, cyclic over
  ! MPI_COMM_WORLD; the gas state, the octree and the Phi_3 leaf rates are
  ! node-shared), the three spectra are summed over the ranks, and rank 0
  ! writes the file.  Everything that depends on the frequency alone is
  ! evaluated once on the output grid before the leaf loop.
  !=========================================================================
  subroutine nebcont_write()
    use mpi
    use octree_mod,    only : leaf_half
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHII, gas_xHeI, gas_xHeII, &
                              gas_xHeIII, gas_ne, gas_Te, gas_nleaf, &
                              gas_leaf_owner_first, gas_leaf_owner_stride
    use photo_xsec,    only : sigma_HI, sigma_HeI, sigma_HeII
    use utility,       only : get_base_name
    use recomb_continuum_mod, only : gamma_tab_pair, gamma_tab_spectrum
    use hei_recombination_photons_mod, only : hei_twophoton_decays
    use hydrogenic_recombination_mod,  only : lyman_line_case
    implicit none
    integer, parameter :: L_FB = 1, L_FF = 2, L_2Q = 3
    real(kind=wp), parameter :: RYD_EV = 13.605693_wp      ! eV/Ryd
    !--- L_nu [erg/s/Hz]: fb+ff(+Milne) (L_FB), the free-free estimate below
    !--- the H edge (L_FF), two-photon (L_2Q)
    real(kind=wp) :: Lnu(NNU_OUT, 3)
    !--- frequency-only factors on nu_out: ground photoionization cross
    !--- sections [cm^2] (0 below each edge), the table pair of each
    !--- frequency, and the two-photon spectral energy of one decay [erg/Hz]
    real(kind=wp) :: sigma_HI_out(NNU_OUT), sigma_HeI_out(NNU_OUT), &
                     sigma_HeII_out(NNU_OUT)
    integer       :: tab_pair_out(NNU_OUT, 3)
    real(kind=wp) :: e2q_HI_out(NNU_OUT), e2q_HeII_out(NNU_OUT), &
                     e2q_HeI_out(NNU_OUT)
    logical       :: above_HI(NNU_OUT), above_HeI(NNU_OUT), above_HeII(NNU_OUT)
    !--- the table coefficients of the leaf on nu_out (1e-40 erg cm^3/s/Hz)
    real(kind=wp) :: gtab_out(NNU_OUT, 3)
    real(kind=wp) :: T, ne, nH, vol, nHII, nHeII_n, nHeIII, dq, r2q_HeI
    real(kind=wp) :: feed_2s_HI, feed_2s_HeII, g2q_HI, g2q_HeII, g2q_HeI
    real(kind=wp) :: nu, gH, gHe1, gHe2, ff1, ff2
    character(len=PATH_LEN) :: outname
    character(len=1) :: lcase
    integer :: il, i, k, unit, ierr

    lcase = lyman_line_case()
    do i = 1, NNU_OUT
       nu = nu_out(i)
       above_HI(i)   = nu >= NU_HI
       above_HeI(i)  = nu >= NU_HEI
       above_HeII(i) = nu >= NU_HEII
       sigma_HI_out(i) = 0.0_wp;  sigma_HeI_out(i) = 0.0_wp
       sigma_HeII_out(i) = 0.0_wp
       if (above_HI(i))   sigma_HI_out(i)   = sigma_HI(nu*RYD_EV)
       if (above_HeI(i))  sigma_HeI_out(i)  = sigma_HeI(nu*RYD_EV)
       if (above_HeII(i)) sigma_HeII_out(i) = sigma_HeII(nu*RYD_EV)
       do k = 1, 3
          tab_pair_out(i, k) = gamma_tab_pair(k, nu)
       end do
       e2q_HI_out(i)   = hydrogenic_twophoton_energy(nu, 0.7496_wp)
       e2q_HeII_out(i) = hydrogenic_twophoton_energy(nu, 3.00_wp)
       e2q_HeI_out(i)  = hei_twophoton_energy(nu)
    end do

    Lnu = 0.0_wp
    do il = gas_leaf_owner_first(), gas_nleaf, gas_leaf_owner_stride()
       nH = gas_nH(il)
       if (nH <= 0.0_wp) cycle
       T  = gas_Te(il);  ne = gas_ne(il)
       vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
       nHII    = nH*gas_xHII(il)
       nHeII_n = nH*par%abund_He*gas_xHeII(il)
       nHeIII  = nH*par%abund_He*gas_xHeIII(il)
       dq      = h2s_quench_denominator(T, ne, nHII, nHeII_n, nHeIII)
       r2q_HeI = 0.0_wp
       if (par%abund_He > 0.0_wp) &
          r2q_HeI = hei_twophoton_decays(il, T, nH, ne, gas_xHI(il), &
                                         gas_xHeI(il), gas_xHeII(il))
       feed_2s_HI   = twophoton_2s_feed(1, T, lcase)
       feed_2s_HeII = twophoton_2s_feed(2, T, lcase)
       do k = 1, 3
          call gamma_tab_spectrum(k, nu_out, tab_pair_out(:, k), T, &
                                  gtab_out(:, k))
       end do
       do i = 1, NNU_OUT
          nu = nu_out(i)
          !--- table region: fb+ff combined; above: Milne fb + explicit ff
          gH   = gtab_out(i, 1)
          gHe1 = gtab_out(i, 2)
          gHe2 = gtab_out(i, 3)
          if (above_HI(i))   gH   = gamma_milne(nu, NU_HI,  sigma_HI_out(i), &
                                    1.0_wp, T)
          if (above_HeI(i))  gHe1 = gamma_milne(nu, NU_HEI, sigma_HeI_out(i), &
                                    0.25_wp, T)
          if (above_HeII(i)) gHe2 = gamma_milne(nu, NU_HEII, sigma_HeII_out(i), &
                                    1.0_wp, T)
          if (above_HI(i)) then
             ff1 = gamma_ff(nu, T, 1)
             gH = gH + ff1
             if (above_HeI(i))  gHe1 = gHe1 + ff1
             if (above_HeII(i)) gHe2 = gHe2 + gamma_ff(nu, T, 2) - ff1
          else
             !--- explicit ff component for the report (table region estimate)
             ff1 = gamma_ff(nu, T, 1)
             ff2 = gamma_ff(nu, T, 2)
             Lnu(i, L_FF) = Lnu(i, L_FF) + ne*((nHII + nHeII_n)*ff1 + nHeIII*ff2) &
                            *1.0e-40_wp*vol
          end if
          Lnu(i, L_FB) = Lnu(i, L_FB) + ne*(nHII*gH + nHeII_n*gHe1 + nHeIII*gHe2) &
                         *1.0e-40_wp*vol
          !--- two-photon coefficients in the 1e-40 units of the others
          g2q_HI   = feed_2s_HI*e2q_HI_out(i)/dq * 1.0e40_wp
          g2q_HeII = feed_2s_HeII*e2q_HeII_out(i) * 1.0e40_wp
          g2q_HeI  = r2q_HeI*e2q_HeI_out(i)*1.0e40_wp
          Lnu(i, L_2Q) = Lnu(i, L_2Q) + ( ne*( nHII*g2q_HI + nHeIII*g2q_HeII ) &
                         + g2q_HeI )*1.0e-40_wp*vol
       end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE, Lnu, size(Lnu), MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (mpar%p_rank /= 0) return

    outname = trim(get_base_name(par%out_file))//'_nebcont.txt'
    open(newunit=unit, file=trim(outname), status='replace')
    write(unit,'(a)') '# MoCHII nebular continuum (grid-integrated), '// &
       'fb_ff port (see nebcont_mod.f90 header for provenance)'
    write(unit,'(a)') '# col 1: nu [Ryd]; 2: lambda [um]; '// &
       '3: L_nu total = fb+ff(+Milne) + two-photon [erg/s/Hz]; '// &
       '4: ff estimate below the H edge; 5: two-photon'
    do i = 1, NNU_OUT
       write(unit,'(f10.5,f12.5,3es14.5)') nu_out(i), &
          0.0911267_wp/nu_out(i), Lnu(i, L_FB) + Lnu(i, L_2Q), Lnu(i, L_FF), &
          Lnu(i, L_2Q)
    end do
    close(unit)
    write(*,'(2a)') ' NEBC: nebular continuum written to: ', trim(outname)
  end subroutine nebcont_write

end module nebcont_mod
