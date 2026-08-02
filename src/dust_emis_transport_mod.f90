module dust_emis_transport_mod
!---------------------------------------------------------------------------
! MoCHII: radiative transfer of the thermal infrared the grains reradiate.
!
! Grains absorb Heat_dust [erg/s/cm^3] of the transported EUV/FUV field and
! reradiate it.  Post-processing that heating into a spectrum (dust_temp_mod,
! sedust_mod) answers what is emitted; it does not answer where the emitted
! photons go.  This module launches them and transports them, so the infrared
! reabsorbed by grains elsewhere in the nebula returns to Heat_dust and the
! dust self-heating loop closes, as in MoCafe's Lucy iteration
! (MoCafe_v2.00/src/lucy_mod.f90 + dustemis_mod.f90).
!
! Transport is the SHARED packet transport (transport_ion_packet).  Every
! cross section is evaluated at the packet's exact photon energy, and in the
! infrared the gas ones vanish identically (the lowest threshold in the metal
! registry is Ca I at 6.113 eV, H I is at 13.598 eV), so an infrared packet is
! carried by grain absorption and grain scattering alone with no separate
! transport path.  That is not assumed: at setup the emission grid is scanned
! with the cross sections themselves and the spectrum is cut at the first
! wavelength that clears both the transported band's floor and every gas
! threshold, so no packet of this band can photoionize.  The ESTIMATORS are
! switched by dust_ir_band_mod's transported_band, so no infrared segment
! enters the ionizing-band mean intensity, photoionization, photoheating or
! Habing-G0 estimators.
!
! Emitting leaf and packet count.  Storing a wavelength CDF for every leaf
! would be leaf x wavelength (2.6e5 x 1.3e3 doubles = several GB), so nothing
! of the kind is kept: each rank walks ITS OWN leaves (the same strided
! partition sedust_mod uses), builds one leaf's spectrum, launches that leaf's
! packets immediately and discards the spectrum.  The packet count per leaf
! comes from the GLOBAL cumulative emission luminosity by systematic
! (stratified) sampling with one shared random offset: leaf il takes the
! packets k whose position (k-1+xi)/N falls in [C(il-1), C(il)]/L_tot.  The
! counts sum to exactly N, every packet carries the same luminosity
! L_tot/N, so the launched luminosity is exactly L_tot, and each rank
! evaluates its own leaves' counts from the cumulative array alone.
!
! Wavelength sampling is CONTINUOUS.  The leaf spectrum is a piecewise-linear
! density in lambda on the emission grid (the same representation the written
! SEDs integrate by trapezoid), so after the bin is drawn from the cumulative,
! the wavelength inside it is drawn from that linear density in closed form
! (sample_linear_spectrum_bin).  No packet is ever emitted at a bin center.
!
! Iteration.  Pass k launches the emission of the current Heat_dust, tallies
! what is reabsorbed, and sets Heat_dust = Heat_euv + Heat_ir; T_dust follows.
! par%dust_niter caps the passes and par%dust_tol is the relative change of
! the total emitted luminosity that stops them.  An H II region is optically
! thin to its own infrared (tau_IR ~ 1e-3), so the reabsorbed fraction f is
! small and successive passes change the total by f, f^2, ...: two passes
! reach 1e-6 there.  Optically thick models need more, which is what the
! reported per-pass change measures.
!
! Documented approximation: the reabsorbed infrared enters the EMITTED POWER
! of each leaf (exactly, through Heat_dust) but not the local mean intensity
! Jsed that sets the SEDust spectral SHAPE, because the ionizing-band tally
! jt_ion does not extend below par%eion_min and an infrared J tally would be
! the leaf x wavelength array excluded above.  The shape error is therefore
! of order the reabsorbed fraction itself, valid while tau_IR << 1; the
! emitted POWER carries no such error at any optical depth.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: dust_emis_transport_setup, transport_dust_emission

  !--- emission wavelength grid [um], ascending: the SEDust model grid when
  !--- par%dust_sed selects the stochastic/PAH spectrum, otherwise the
  !--- equilibrium modified-blackbody grid built here.
  real(kind=wp), allocatable :: lam_emis(:)
  !--- grain absorption cross section relative to the reference extinction,
  !--- (1-albedo) C_ext(lambda)/C_ext(lambda_ref), on that grid.  Only the
  !--- equilibrium spectrum needs it; SEDust carries its own optics.
  real(kind=wp), allocatable :: sabs_emis(:)
  integer :: nl_emis = 0
  logical :: from_sedust = .false.
  !--- leading grid nodes whose emission is set to zero because they lie in
  !--- the transported EUV/FUV band rather than in the thermal infrared, and
  !--- the shortest wavelength a packet can then be drawn at.  Grain thermal
  !--- emission there vanishes for any temperature below grain sublimation
  !--- (~2000 K puts 1e-15 of B_lambda shortward of 0.2 um), so the cut
  !--- removes nothing physical; it makes the band separation exact.
  integer       :: nzero_emis = 0
  real(kind=wp) :: lam_hardest = 0.0_wp
  !--- emission fraction the cut discarded, reported once.
  real(kind=wp) :: frac_discarded = 0.0_wp
  logical       :: discard_reported = .false.

  !--- equilibrium modified-blackbody grid: 0.1 - 3000 um.  A grain at
  !--- T_d < 3000 K puts less than 1e-6 of its emission shortward of 0.1 um
  !--- (Wien peak 0.1 um needs 2.9e4 K) and the Rayleigh-Jeans tail beyond
  !--- 3000 um carries less than 1e-4 for T_d > 3 K, so the grid truncation
  !--- is negligible against the reabsorbed fraction it feeds.
  integer,       parameter :: NLAM_BB = 300
  real(kind=wp), parameter :: LAM_BB_MIN = 0.1_wp, LAM_BB_MAX = 3.0e3_wp

contains

  !=========================================================================
  subroutine dust_emis_transport_setup()
    use mpi
    use octree_mod,      only : amr_grid
    use gas_opacity_mod, only : ion_dust_sabs_at
    use dust_ir_band_mod,only : dust_ir_band_setup
    use sedust_mod,      only : sedust_setup, sedust_lambda_grid
    implicit none
    real(kind=wp) :: E_floor, lam_floor, lam_start
    integer :: k, kpass, ierr

    from_sedust = par%dust_sed
    if (allocated(lam_emis)) deallocate(lam_emis)
    if (allocated(sabs_emis)) deallocate(sabs_emis)

    !--- the thermal infrared is the radiation below the transported band's
    !--- floor; at and above it the packets would be band photons.
    E_floor = par%eion_min
    if (par%add_fuv) E_floor = par%efuv_min
    lam_floor = 1.23984_wp/E_floor

    if (from_sedust) then
       call sedust_setup()
       lam_emis = sedust_lambda_grid()
       nl_emis  = size(lam_emis)
    else
       lam_start = max(LAM_BB_MIN, lam_floor)
       nl_emis   = NLAM_BB
       allocate(lam_emis(nl_emis), sabs_emis(nl_emis))
       do k = 1, nl_emis
          lam_emis(k) = lam_start*(LAM_BB_MAX/lam_start) &
                        **(real(k-1,wp)/real(nl_emis-1,wp))
          sabs_emis(k) = ion_dust_sabs_at(1.23984_wp/lam_emis(k))
       end do
    end if

    !--- Shortest wavelength a packet of this band may be drawn at.  It has to
    !--- clear BOTH the transported band's floor and every gas photoionization
    !--- threshold: an infrared packet with a nonzero gas cross section would
    !--- be attenuated by gas while the infrared estimators book all of its
    !--- extinction as grain absorption, i.e. energy that ionized an atom would
    !--- come back as grain heating.  The thresholds are not hard-wired; the
    !--- grid is scanned with the cross sections themselves, so a registry
    !--- change cannot invalidate the cut.
    kpass = 0
    do k = 1, nl_emis
       if (lam_emis(k) < lam_floor) cycle
       if (gas_absorption_at(1.23984_wp/lam_emis(k)) <= 0.0_wp) then
          kpass = k
          exit
       end if
    end do
    if (kpass == 0) then
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: no wavelength of the dust emission grid is free of both '// &
          'the transported band and every gas photoionization threshold.'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    !--- zeroing the nodes up to and INCLUDING kpass leaves lam_emis(kpass) as
    !--- the shortest wavelength the linear ramp out of the zero can reach.
    nzero_emis  = 0
    if (kpass > 1) nzero_emis = kpass
    lam_hardest = lam_emis(kpass)

    call dust_ir_band_setup(amr_grid%nleaf, lam_emis)

    if (mpar%p_rank == 0) then
       write(*,'(a)') ' DUSTRT: dust thermal emission is transported'
       if (from_sedust) then
          write(*,'(a,i0,a,2es11.3)') &
             ' DUSTRT: emission spectrum = SEDust stochastic + PAH, ', &
             nl_emis, ' wavelengths [um] = ', lam_emis(1), lam_emis(nl_emis)
       else
          write(*,'(a,i0,a,2es11.3)') &
             ' DUSTRT: emission spectrum = equilibrium modified blackbody, ', &
             nl_emis, ' wavelengths [um] = ', lam_emis(1), lam_emis(nl_emis)
       end if
       write(*,'(a,es11.3,a,f8.3,a)') &
          ' DUSTRT: emission cut below ', lam_hardest, ' um (', &
          1.23984_wp/lam_hardest, ' eV: band floor and gas thresholds cleared)'
    end if
  end subroutine dust_emis_transport_setup

  !=========================================================================
  ! Total gas photoionization cross section at photon energy E [eV]: H I,
  ! He I, He II and every active metal stage.  Zero below all thresholds.
  !=========================================================================
  real(kind=wp) function gas_absorption_at(energy) result(sgas)
    use photo_xsec, only : sigma_HI, sigma_HeI, sigma_HeII
    implicit none
    real(kind=wp), intent(in) :: energy
    integer :: nmetal

    sgas = sigma_HI(energy) + sigma_HeI(energy) + sigma_HeII(energy)
    if (par%use_metals) then
       block
         use species_mod, only : species_packet_cross_sections
         real(kind=wp) :: smetal(MAX_METAL_TRANSITIONS)
         call species_packet_cross_sections(energy, smetal, nmetal)
         if (nmetal > 0) sgas = sgas + sum(smetal(1:nmetal))
       end block
    end if
  end function gas_absorption_at

  !=========================================================================
  ! Transport the dust thermal emission and close the self-heating loop.
  ! heat_dust enters as the EUV/FUV grain heating and leaves as the total,
  ! EUV/FUV plus reabsorbed infrared.  t_dust is left consistent with it.
  !=========================================================================
  subroutine transport_dust_emission(heat_dust)
    use mpi
    use random,           only : rand_number
    use octree_mod,       only : amr_grid, leaf_half
    use jtally_mod,       only : slab_tally_on
    use dust_temp_mod,    only : dust_temp_compute
    use dust_ir_band_mod, only : transported_band, BAND_IONIZING, BAND_DUST_IR, &
                                 dust_ir_band_reset, dust_ir_band_reduce, &
                                 dust_ir_heating, dust_ir_energy_budget
    implicit none
    real(kind=wp), intent(inout) :: heat_dust(:)
    real(kind=wp), allocatable :: heat_euv(:), heat_ir(:), lcum(:), pdf(:), cdf(:)
    real(kind=wp) :: L_tot, Lpacket, vol, xi, esum, L_abs_path
    real(kind=wp) :: emitted, extinct, escaped, closure, L_new, L_prev, drel
    integer(kind=int64) :: npacket, k_lo, k_hi, ip
    integer :: il, iter, ierr, nleaf, niter_done
    logical :: slab_saved, converged

    nleaf = amr_grid%nleaf
    allocate(heat_euv(nleaf), heat_ir(nleaf), lcum(0:nleaf))
    allocate(pdf(nl_emis), cdf(nl_emis))
    heat_euv = heat_dust(1:nleaf)

    npacket = int(max(par%dust_no_photons, 1.0_wp), int64)

    !--- the escaping infrared must not be booked as an ionizing-band packet
    !--- leaving the plane-parallel slab boundary.
    slab_saved   = slab_tally_on
    slab_tally_on = .false.
    transported_band = BAND_DUST_IR

    L_prev = 0.0_wp
    converged = .false.
    niter_done = 0
    do iter = 1, max(par%dust_niter, 1)
       niter_done = iter
       !--- cumulative emission luminosity over ALL leaves (identical on every
       !--- rank: heat_dust is the reduced, global array).
       lcum(0) = 0.0_wp
       do il = 1, nleaf
          vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
          lcum(il) = lcum(il-1) + max(heat_dust(il), 0.0_wp)*vol
       end do
       L_tot = lcum(nleaf)
       !--- the change measured after pass 1 is against the emission BEFORE any
       !--- infrared was transported, i.e. the reabsorbed fraction itself.
       if (iter == 1) L_prev = L_tot
       if (.not. (L_tot > 0.0_wp)) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             ' DUSTRT: no dust emission (Heat_dust = 0); nothing transported.'
          exit
       end if
       Lpacket = L_tot/real(npacket, wp)

       !--- one shared stratification offset per pass, so the systematic
       !--- allocation is unbiased and identical on every rank.
       if (mpar%p_rank == 0) xi = rand_number()
       call MPI_BCAST(xi, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       xi = min(max(xi, 1.0e-12_wp), 1.0_wp - 1.0e-12_wp)

       call dust_ir_band_reset()
       do il = mpar%p_rank+1, nleaf, mpar%nproc
          !--- packets k whose stratified position (k-1+xi)/N falls in this
          !--- leaf's share [C(il-1), C(il)]/L_tot of the cumulative emission.
          !--- The endpoint expressions are shared with the neighbouring
          !--- leaves, so the counts telescope to exactly N over the grid.
          k_lo = ceiling(real(npacket,wp)*(lcum(il-1)/L_tot) - xi, int64) + 1_int64
          k_hi = ceiling(real(npacket,wp)*(lcum(il  )/L_tot) - xi, int64)
          if (k_hi < k_lo) cycle
          call leaf_emission_spectrum(il, pdf, esum)
          if (.not. (esum > 0.0_wp)) cycle
          call emission_cumulative(pdf, cdf)
          do ip = k_lo, k_hi
             call launch_and_transport(il, cdf, pdf, Lpacket)
          end do
       end do

       call dust_ir_band_reduce()
       call dust_ir_heating(heat_ir)
       heat_dust(1:nleaf) = heat_euv + heat_ir
       call dust_temp_compute(heat_dust)

       call dust_ir_energy_budget(emitted, extinct, escaped)
       closure = abs(extinct + escaped - emitted)/max(emitted, tiny(1.0_wp))
       !--- L_abs_path is the SAME absorbed power seen by the path estimator
       !--- that feeds heat_ir, so it is the physical one; without dust
       !--- scattering it equals the extinction ledger of the packets exactly
       !--- (Int kappa_abs exp(-tau) dl = 1 - exp(-tau_edge) cell by cell),
       !--- with scattering it is smaller by the scattered-out energy.
       L_abs_path = 0.0_wp
       L_new      = 0.0_wp
       do il = 1, nleaf
          vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
          L_abs_path = L_abs_path + heat_ir(il)*vol
          L_new      = L_new + heat_dust(il)*vol
       end do
       drel = abs(L_new - L_prev)/max(L_new, tinest)
       if (mpar%p_rank == 0) then
          write(*,'(a,i0,a,i0,a,es13.6)') ' DUSTRT: pass ', iter, ': ', &
             npacket, ' packets, L_emit = ', L_tot
          write(*,'(a,3es14.6)') &
             ' DUSTRT: IR energy emitted/extinguished/escaped = ', &
             emitted, extinct, escaped
          write(*,'(a,es11.3,a,es11.3)') &
             ' DUSTRT: IR energy closure = ', closure, &
             ',  extinguished fraction = ', extinct/max(emitted, tinest)
          write(*,'(a,es13.6,a,es11.3)') &
             ' DUSTRT: grain-absorbed IR (path estimator) = ', L_abs_path, &
             ' erg/s;  vs ledger, rel. diff = ', &
             abs(L_abs_path - extinct)/max(emitted, tinest)
          write(*,'(a,es11.3)') &
             ' DUSTRT: relative change of the total emitted luminosity = ', drel
       end if
       if (closure > 5.0e-13_wp) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: the transported infrared does not conserve energy.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       !--- a leaf with Heat_dust > 0 whose spectrum comes back empty would
       !--- drop its packets, and the launched luminosity would fall short of
       !--- the emission it is supposed to carry.
       if (abs(emitted - L_tot) > 1.0e-12_wp*L_tot) then
          if (mpar%p_rank == 0) write(*,'(a,es11.3)') &
             ' DUSTRT: WARNING: launched luminosity misses the emission by ', &
             (emitted - L_tot)/L_tot
       end if
       L_prev = L_new
       if (drel < par%dust_tol) then
          converged = .true.
          exit
       end if
    end do

    transported_band = BAND_IONIZING
    slab_tally_on    = slab_saved

    !--- how much emission the band cut removed, worst leaf over all ranks.
    if (nzero_emis > 0 .and. .not. discard_reported) then
       call MPI_ALLREDUCE(MPI_IN_PLACE, frac_discarded, 1, &
                          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
       discard_reported = .true.
       if (mpar%p_rank == 0 .and. frac_discarded > 1.0e-8_wp) &
          write(*,'(a,es11.3)') ' DUSTRT: WARNING: the band cut discarded up '// &
             'to this fraction of a leaf''s emission: ', frac_discarded
    end if

    if (mpar%p_rank == 0) then
       if (converged) then
          write(*,'(a,i0,a)') ' DUSTRT: dust self-heating converged after ', &
             niter_done, ' pass(es)'
       else
          write(*,'(a,i0,a)') ' DUSTRT: WARNING: dust self-heating still '// &
             'changing after ', niter_done, ' pass(es); raise par%dust_niter'
       end if
    end if
    call dust_ir_write_escaped()

    deallocate(heat_euv, heat_ir, lcum, pdf, cdf)
  end subroutine transport_dust_emission

  !=========================================================================
  ! Emission spectrum of one leaf as a density in lambda on lam_emis,
  ! pdf(k) proportional to L_lambda, with esum = Int pdf dlambda.
  !=========================================================================
  subroutine leaf_emission_spectrum(il, pdf, esum)
    use sedust_mod, only : sedust_leaf_spectrum
    implicit none
    integer,       intent(in)  :: il
    real(kind=wp), intent(out) :: pdf(:), esum
    real(kind=wp) :: esum_full
    if (from_sedust) then
       call sedust_leaf_spectrum(il, pdf, esum)
    else
       call equilibrium_leaf_spectrum(il, pdf, esum)
    end if
    if (nzero_emis > 0 .and. esum > 0.0_wp) then
       esum_full = esum
       pdf(1:nzero_emis) = 0.0_wp
       esum = trapezoid(lam_emis, pdf)
       frac_discarded = max(frac_discarded, (esum_full - esum)/esum_full)
    end if
  end subroutine leaf_emission_spectrum

  !=========================================================================
  ! Equilibrium modified blackbody of one leaf:
  !   L_lambda  proportional to  s_abs(lambda) B_lambda(T_dust),
  ! the same emission dust_temp_mod's dust_ir_write integrates, with s_abs
  ! the grain absorption the transport itself uses.
  !=========================================================================
  subroutine equilibrium_leaf_spectrum(il, pdf, esum)
    use dust_temp_mod, only : t_dust
    implicit none
    integer,       intent(in)  :: il
    real(kind=wp), intent(out) :: pdf(:), esum
    real(kind=wp) :: T
    integer :: k

    pdf = 0.0_wp;  esum = 0.0_wp
    if (.not. allocated(t_dust)) return
    T = t_dust(il)
    if (.not. (T > 0.0_wp)) return
    do k = 1, nl_emis
       pdf(k) = sabs_emis(k)*planck_lambda(lam_emis(k), T)
    end do
    esum = trapezoid(lam_emis, pdf)
  end subroutine equilibrium_leaf_spectrum

  !=========================================================================
  ! B_lambda(T) up to a constant factor, lambda in um.
  !=========================================================================
  pure real(kind=wp) function planck_lambda(lam_um, T) result(b)
    implicit none
    real(kind=wp), intent(in) :: lam_um, T
    real(kind=wp) :: lam_cm, x
    lam_cm = lam_um*1.0e-4_wp
    x = h_planck_cgs*clight_cgs/(lam_cm*kboltz_cgs*T)
    if (x > 700.0_wp) then
       b = 0.0_wp
    else
       b = 1.0_wp/(lam_cm**5*(exp(x) - 1.0_wp))
    end if
  end function planck_lambda

  !=========================================================================
  pure real(kind=wp) function trapezoid(x, y) result(s)
    implicit none
    real(kind=wp), intent(in) :: x(:), y(:)
    integer :: k
    s = 0.0_wp
    do k = 1, size(x)-1
       s = s + 0.5_wp*(y(k) + y(k+1))*(x(k+1) - x(k))
    end do
  end function trapezoid

  !=========================================================================
  ! Cumulative emission of a leaf, cdf(k) = Int_lam(1)^lam(k) pdf dlambda,
  ! by the SAME trapezoid the written SEDs integrate: the sampled spectrum is
  ! then the piecewise-linear density those integrals define.
  !=========================================================================
  subroutine emission_cumulative(pdf, cdf)
    implicit none
    real(kind=wp), intent(in)  :: pdf(:)
    real(kind=wp), intent(out) :: cdf(:)
    integer :: k
    cdf(1) = 0.0_wp
    do k = 1, nl_emis-1
       cdf(k+1) = cdf(k) + 0.5_wp*(pdf(k) + pdf(k+1)) &
                  *(lam_emis(k+1) - lam_emis(k))
    end do
  end subroutine emission_cumulative

  !=========================================================================
  ! Wavelength inside one bin from a density that varies LINEARLY between
  ! p0 at lam0 and p1 at lam1.  Inverting
  !   F(t) = [p0 t + (p1-p0) t^2/2] / [(p0+p1)/2],  t = (lam-lam0)/(lam1-lam0)
  ! gives t = 2c/(p0 + sqrt(p0^2 + 2 d c)) with d = p1-p0, c = u (p0+p1)/2,
  ! the numerically stable root (t -> u as d -> 0, t -> sqrt(u) for p0 = 0,
  ! t -> 1-sqrt(1-u) for p1 = 0).
  !=========================================================================
  pure real(kind=wp) function sample_linear_spectrum_bin(lam0, lam1, p0, p1, u) &
       result(lam)
    implicit none
    real(kind=wp), intent(in) :: lam0, lam1, p0, p1, u
    real(kind=wp) :: uu, d, twoc, t, disc
    uu = min(max(u, 0.0_wp), 1.0_wp)
    if (p0 + p1 <= 0.0_wp) then
       lam = lam0 + uu*(lam1 - lam0)
       return
    end if
    d    = p1 - p0
    twoc = uu*(p0 + p1)
    disc = p0*p0 + d*twoc
    t    = twoc/(p0 + sqrt(max(disc, 0.0_wp)))
    t    = min(max(t, 0.0_wp), 1.0_wp)
    lam  = lam0 + t*(lam1 - lam0)
  end function sample_linear_spectrum_bin

  !=========================================================================
  ! One emitted packet: uniform position in the leaf, isotropic direction,
  ! continuous wavelength from the leaf spectrum, then the shared transport.
  !=========================================================================
  subroutine launch_and_transport(il, cdf, pdf, Lpacket)
    use random,           only : rand_number
    use octree_mod,       only : leaf_cx, leaf_cy, leaf_cz, leaf_half
    use raytrace_amr_mod, only : transport_ion_packet
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: cdf(:), pdf(:), Lpacket
    type(photon_type) :: photon
    real(kind=wp) :: u, ulo, ubin, cost, sint, phi, half, lam
    integer :: lo, hi, mid, k

    half = leaf_half(il)
    photon%x = leaf_cx(il) + (2.0_wp*rand_number() - 1.0_wp)*half
    photon%y = leaf_cy(il) + (2.0_wp*rand_number() - 1.0_wp)*half
    photon%z = leaf_cz(il) + (2.0_wp*rand_number() - 1.0_wp)*half

    cost = 2.0_wp*rand_number() - 1.0_wp
    sint = sqrt(1.0_wp - cost*cost)
    phi  = twopi*rand_number()
    photon%kx = sint*cos(phi)
    photon%ky = sint*sin(phi)
    photon%kz = cost

    !--- wavelength: bin from the cumulative, then continuously inside it.
    u  = rand_number()*cdf(nl_emis)
    lo = 1;  hi = nl_emis-1
    do while (lo < hi)
       mid = (lo + hi)/2
       if (u <= cdf(mid+1)) then
          hi = mid
       else
          lo = mid + 1
       end if
    end do
    k    = lo
    ulo  = cdf(k)
    ubin = 0.0_wp
    if (cdf(k+1) > ulo) ubin = (u - ulo)/(cdf(k+1) - ulo)
    lam = sample_linear_spectrum_bin(lam_emis(k), lam_emis(k+1), &
                                     pdf(k), pdf(k+1), ubin)

    photon%lambda    = lam
    photon%energy_eV = 1.23984_wp/lam
    photon%inu       = 1          ! no ionizing-band bin: the infrared
                                  ! estimators never address jt_ion
    photon%il        = k
    photon%icell_amr = il
    photon%nscatt    = 0
    photon%wgt       = 1.0_wp
    photon%inside    = .true.
    photon%Lpacket   = Lpacket
    photon%s_ext     = 1.0_wp
    photon%from_external = .false.
    photon%id        = 0_int64
    photon%istream   = 0

    call transport_ion_packet(photon)
  end subroutine launch_and_transport

  !=========================================================================
  ! Escaping infrared spectrum -> '<base>_dustesc.txt'.  The emitted
  ! counterparts ('<base>_dustir.txt', '<base>_dustsed.txt') are untouched:
  ! this one is what leaves the grid after transport, so the two differ by
  ! the reabsorbed fraction and, in an optically thick model, by the
  ! reprocessing that fraction drives.
  !=========================================================================
  subroutine dust_ir_write_escaped()
    use dust_ir_band_mod, only : dust_ir_escaped_spectrum, dust_ir_energy_budget
    use utility,          only : get_base_name
    implicit none
    real(kind=wp), allocatable :: lam(:), Llam(:), dlam(:)
    real(kind=wp) :: emitted, extinct, escaped, esum
    character(len=192) :: outname
    integer :: k, unit

    if (mpar%p_rank /= 0) return
    call dust_ir_escaped_spectrum(lam, Llam, dlam)
    call dust_ir_energy_budget(emitted, extinct, escaped)
    esum = sum(Llam*dlam)

    outname = trim(get_base_name(par%out_file))//'_dustesc.txt'
    open(newunit=unit, file=trim(outname), status='replace')
    write(unit,'(a)') '# MoCHII escaping dust thermal emission '// &
       '(transported infrared packets that left the grid)'
    if (from_sedust) then
       write(unit,'(a)') '# emission spectrum: SEDust stochastic + PAH'
    else
       write(unit,'(a)') '# emission spectrum: equilibrium modified blackbody'
    end if
    write(unit,'(a,es14.6,a,es14.6,a,es14.6)') '# L_emit = ', emitted, &
       ' erg/s;  L_abs = ', extinct, ' erg/s;  L_esc = ', escaped
    write(unit,'(a,es14.6,a)') &
       '# sum L_lambda dlambda over the bins below = ', esum, &
       ' erg/s (equals L_esc by construction)'
    write(unit,'(a)') '# lambda[um]   L_lambda[erg/s/um]'
    do k = 1, size(lam)
       write(unit,'(es13.5,es15.6)') lam(k), Llam(k)
    end do
    close(unit)
    write(*,'(2a)') ' DUSTRT: escaping IR spectrum written to: ', trim(outname)
    deallocate(lam, Llam, dlam)
  end subroutine dust_ir_write_escaped

end module dust_emis_transport_mod
