program check_metal_recombination_continuum
!---------------------------------------------------------------------------
! Gate for the metal ground continua of the diffuse field (diffuse_mod
! channel 6, milne_recomb_spectrum_mod::milne_metal_*), on every transition
! of the 11-element, 77-stage registry.
!
!  (a) The cross section the Milne relation uses is phfit2's outer shell:
!      species_outer_shell_cross_section against Verner's phfit2.f called
!      for that shell alone (PH2 below the first inner edge, PH1 above), on
!      a log grid 5-1000 eV, edges +-1% skipped (NIST vs PH1 thresholds).
!      Below the last node E_join of a PHTABLE (Ne III since 2026-10-11) the
!      reference is that table, read from the element file here and
!      interpolated linearly in eps = E - E_th, at every node and interval
!      midpoint; phfit2 at and above E_join.  The integrals (b)-(d) then run
!      over the resonances the table resolves (node-stepped quadrature in
!      milne_recomb_spectrum_mod); tests/ne3_ground_term_photoionization
!      checks the table against its source and the sampled spectrum.
!  (b) The band coefficient the leaf luminosity is built from equals the
!      energy the captures radiate inside the band,
!          P(T) = Int_{E_lo}^{eion_max} E d alpha_1(E),
!      by an independent quadrature (milne_metal_energy_integral: step
!      0.02 kT from E_th, against the 512-node table interpolated in T), at
!      table nodes and between them; both with the capture cap
!      min(1, alpha_rec/alpha_1) the channel carries.
!  (c) Energy bookkeeping.  Over all E the captures radiate
!          Int E d alpha_1 = alpha_1 IP + beta_1,
!          beta_1 = kT [(3/2) alpha_1 + d alpha_1/d ln T],
!      the kinetic energy taken from the gas (inside the beta of the total
!      rate that species_mod::metal_ionization_cooling removes) plus the
!      ionization potential.  Checked with alpha_1(T) from the direct
!      quadrature and d alpha_1/d ln T by a centered difference -- a test of
!      the identity on the code's own integrand, independent of (b).
!  (d) The photon rate the packets represent: equal-energy packets drawn
!      from milne_metal_sample_energy (stratified u) carry P <1/E>, which
!      must equal the in-band photon rate Int_band d alpha_1.
!  (f) Prefactor and weights on a case with an exact answer: C VI(1s) is
!      hydrogenic, so alpha_1(Z = 6, T) = 3 alpha_1(Z = 2, T/9), against
!      He II's alpha_A - alpha_B (recomb_mod, the Milne fit on the exact
!      hydrogenic cross section); within 2% (accuracy of the VY95 1s fit).
!  (e) Reported, not gated: alpha_1/alpha_RR at 1e4 K (alpha_1 is a part
!      of the radiative rate), the share of the Milne integral above the
!      first excited term of the ion (where the outer-shell fit counts
!      excited-term photoionization as ground recombination), and the bound
!      n_e n_+ (alpha_tot - alpha_1)(IP + kT) on the ionizing photons of
!      the captures into excited levels that the diffuse field does not
!      carry, relative to the ground continuum.
!
! Build and run: tests/metal_recombination_continuum/run_check.sh
!---------------------------------------------------------------------------
  use define
  use species_mod, only : species_setup, species_outer_shell_cross_section, &
                          n_elements, elem_name, elem_nstage, elem_eth, &
                          species_transition, species_transition_count, &
                          species_radiative_recombination, &
                          species_recombination_coefficient
  use milne_recomb_spectrum_mod, only : milne_setup, milne_metal_count, &
                          milne_metal_band_coefficient, milne_metal_sample_energy, &
                          milne_metal_alpha1, milne_metal_energy_integral, &
                          milne_metal_excited_share, milne_metal_capture_cap
  use recomb_mod, only : alphaA_HeIII, alphaB_HeIII
  use mpi
  implicit none
  integer, parameter :: NEL = 11
  character(len=2), parameter :: ELN(NEL) = &
     ['c ','n ','o ','ne','s ','ar','mg','fe','si','cl','ca']
  integer, parameter :: ZEL(NEL) = [6,7,8,10,16,18,12,26,14,17,20]
  integer, parameter :: NE_GRID = 400, NS = 100000
  real(kind=wp), parameter :: MB = 1.0e-18_wp
  real(kind=wp), parameter :: TOL_SIGMA = 2.0e-5_wp, EDGE_SKIP = 1.0e-2_wp
  real(kind=wp), parameter :: TOL_BAND = 2.0e-4_wp, TOL_IDENTITY = 1.0e-5_wp
  real(kind=wp), parameter :: TOL_RATE = 2.0e-4_wp
  !--- table nodes (log T grid 1e2-1e6 K, 192 nodes) and points between them
  real(kind=wp), parameter :: TT(8) = [3.0e3_wp, 5.0e3_wp, 7.3e3_wp, 1.0e4_wp, &
                                       1.37e4_wp, 2.0e4_wp, 3.1e4_wp, 5.0e4_wp]
  real(kind=wp), parameter :: DLNT = 1.0e-3_wp
  real :: e4, s4, eth4(7)
  real :: ph1
  integer :: ninn, ntot
  common/ph1/ph1(6,30,30,7)
  common/ninn/ninn(30)
  common/ntot/ntot(30)
  real(kind=wp) :: E, sig, sref, dev, worst_sig, worst_band, worst_id, worst_rate
  real(kind=wp) :: T, kT_eV, p_tab, p_dir, a1, ap, am, beta1, ptot, sinv, u, ea
  real(kind=wp) :: rate_dir, arr, atot, share, bound, ip
  integer, parameter :: MAXTAB = 200000
  real(kind=wp), parameter :: TOL_TABLE = 1.0e-8_wp
  real(kind=wp), allocatable :: teps(:), tsig(:)
  real(kind=wp) :: ejoin
  integer :: ntab
  integer :: ierr, k, ie, it, nz, ne, is, is_outer, j, itt, i, ntest, nskip
  logical :: ok, near
  character(len=2) :: nm
  character(len=8) :: ion

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir = 'data/atomic'
  par%use_metals = .true.
  !--- no (T, n_e) line-cooling tables: the gate does not use the cooling,
  !--- and the full CHIANTI tables take hours on one process
  par%cooling_model = 'low_density'
  par%eion_min = 13.598_wp
  par%eion_max = 100.0_wp
  par%abund_C  = 1.0e-4_wp;  par%abund_N  = 1.0e-4_wp
  par%abund_O  = 1.0e-4_wp;  par%abund_Ne = 1.0e-4_wp
  par%abund_S  = 1.0e-5_wp;  par%abund_Ar = 1.0e-5_wp
  par%abund_Mg = 1.0e-5_wp;  par%abund_Fe = 1.0e-5_wp
  par%abund_Si = 1.0e-5_wp;  par%abund_Cl = 1.0e-7_wp
  par%abund_Ca = 1.0e-8_wp
  call species_setup(1)
  call milne_setup()
  allocate(teps(MAXTAB), tsig(MAXTAB))
  ok = .true.
  if (milne_metal_count() /= species_transition_count()) then
     write(*,'(a)') 'FAIL: metal channel count differs from the transition count'
     ok = .false.
  end if
  write(*,'(a,i0)') 'metal transitions: ', milne_metal_count()

  !--- (a) outer-shell cross section against phfit2
  worst_sig = 0.0_wp;  ntest = 0;  nskip = 0
  do k = 1, milne_metal_count()
     call species_transition(k, ie, it)
     nm = elem_name(ie)
     nz = 0
     do j = 1, NEL
        if (trim(ELN(j)) == trim(nm)) nz = ZEL(j)
     end do
     ne = nz - it + 1
     do is = 1, 7
        eth4(is) = ph1(1, nz, ne, is)
     end do
     is_outer = 0
     do is = 1, 7
        if (eth4(is) <= 0.0) cycle
        if (is_outer == 0) then
           is_outer = is
        else if (abs(eth4(is) - elem_eth(ie,it)) < &
                 abs(eth4(is_outer) - elem_eth(ie,it))) then
           is_outer = is
        end if
     end do
     call read_phtable(trim(nm), it, teps, tsig, ntab)
     ejoin = 0.0_wp
     if (ntab > 0) then
        ejoin = elem_eth(ie,it) + teps(ntab)
        do j = 1, 2*ntab - 2
           i = (j + 1)/2
           u = 0.5_wp*real(mod(j + 1, 2), wp)       ! 0 at a node, 0.5 between
           E = elem_eth(ie,it) + teps(i) + u*(teps(i+1) - teps(i))
           sref = tsig(i) + u*(tsig(i+1) - tsig(i))
           sig = species_outer_shell_cross_section(ie, it, E)
           dev = abs(sig - sref)/max(sref, 1.0e-30_wp)
           ntest = ntest + 1
           worst_sig = max(worst_sig, dev)
           if (dev > TOL_TABLE) then
              ok = .false.
              write(*,'(a,a2,i2,f12.6,2es14.6,es10.2)') '  FAIL table sigma ', nm, it, E, sig, sref, dev
           end if
        end do
     end if
     do j = 0, NE_GRID
        E = 5.0_wp*(1000.0_wp/5.0_wp)**(real(j,wp)/NE_GRID)
        if (E < ejoin) cycle                 ! on the table: checked above
        near = .false.
        do is = 1, 7
           if (eth4(is) > 0.0 .and. abs(E/eth4(is) - 1.0_wp) < EDGE_SKIP) near = .true.
        end do
        if (abs(E/elem_eth(ie,it) - 1.0_wp) < EDGE_SKIP) near = .true.
        if (near) then;  nskip = nskip + 1;  cycle;  end if
        sig = species_outer_shell_cross_section(ie, it, E)
        e4 = real(E)
        call phfit2(nz, ne, is_outer, e4, s4)
        sref = real(s4, wp)*MB
        if (sref <= 0.0_wp .and. sig <= 0.0_wp) cycle
        dev = abs(sig - sref)/max(sref, 1.0e-30_wp)
        ntest = ntest + 1
        worst_sig = max(worst_sig, dev)
        if (dev > TOL_SIGMA) then
           ok = .false.
           write(*,'(a,a2,i2,f10.3,2es14.6,es10.2)') '  FAIL outer sigma ', nm, it, E, sig, sref, dev
        end if
     end do
  end do
  write(*,'(a,i0,a,i0,a,es10.3)') '(a) outer shell vs phfit2: ', ntest, &
     ' points (', nskip, ' near edges skipped), worst ', worst_sig

  !--- (b), (c), (d)
  worst_band = 0.0_wp;  worst_id = 0.0_wp;  worst_rate = 0.0_wp
  do k = 1, milne_metal_count()
     call species_transition(k, ie, it)
     do itt = 1, size(TT)
        T = TT(itt)
        kT_eV = kboltz_cgs*T/ev2erg
        !--- (b)
        p_tab = milne_metal_band_coefficient(k, T)
        p_dir = milne_metal_energy_integral(k, T, max(elem_eth(ie,it), par%eion_min), &
                                            par%eion_max) &
                *min(1.0_wp, species_recombination_coefficient(ie, it, T) &
                             /milne_metal_alpha1(k, T))
        if (p_dir > 0.0_wp .and. p_dir > 1.0e-300_wp) then
           dev = abs(p_tab/p_dir - 1.0_wp)
           worst_band = max(worst_band, dev)
           if (dev > TOL_BAND) then
              ok = .false.
              write(*,'(a,a2,i2,f9.0,2es14.6,es10.2)') '  FAIL band coefficient ', &
                 elem_name(ie), it, T, p_tab, p_dir, dev
           end if
        else if (p_tab > 0.0_wp) then
           ok = .false.
           write(*,'(a,a2,i2)') '  FAIL band coefficient nonzero where the band has none ', &
              elem_name(ie), it
        end if
        !--- (c) the identity on all E
        a1 = milne_metal_alpha1(k, T)
        ap = milne_metal_alpha1(k, T*exp(DLNT))
        am = milne_metal_alpha1(k, T*exp(-DLNT))
        beta1 = kboltz_cgs*T*(1.5_wp*a1 + (ap - am)/(2.0_wp*DLNT))
        ptot  = milne_metal_energy_integral(k, T, elem_eth(ie,it), huge(1.0_wp))
        ip = elem_eth(ie,it)*ev2erg
        dev = abs((a1*ip + beta1)/ptot - 1.0_wp)
        worst_id = max(worst_id, dev)
        if (dev > TOL_IDENTITY) then
           ok = .false.
           write(*,'(a,a2,i2,f9.0,2es14.6,es10.2)') '  FAIL alpha_1 IP + beta_1 ', &
              elem_name(ie), it, T, a1*ip + beta1, ptot, dev
        end if
        !--- (d) photon rate of equal-energy packets
        if (p_tab > 0.0_wp) then
           sinv = 0.0_wp
           do i = 1, NS
              u  = (real(i,wp) - 0.5_wp)/real(NS,wp)
              ea = milne_metal_sample_energy(k, T, u)
              if (ea < max(elem_eth(ie,it), par%eion_min)*(1.0_wp - 1.0e-12_wp) .or. &
                  ea > par%eion_max) then
                 ok = .false.
                 write(*,'(a,a2,i2,f9.0,es14.6)') '  FAIL sampled energy outside the band ', &
                    elem_name(ie), it, T, ea
              end if
              sinv = sinv + 1.0_wp/(ea*ev2erg)
           end do
           sinv = sinv/real(NS,wp)
           rate_dir = band_photon_rate(k, ie, it, T) &
                      *min(1.0_wp, species_recombination_coefficient(ie, it, T) &
                                   /milne_metal_alpha1(k, T))
           dev = abs(p_tab*sinv/rate_dir - 1.0_wp)
           worst_rate = max(worst_rate, dev)
           if (dev > TOL_RATE) then
              ok = .false.
              write(*,'(a,a2,i2,f9.0,2es14.6,es10.2)') '  FAIL represented photon rate ', &
                 elem_name(ie), it, T, p_tab*sinv, rate_dir, dev
           end if
        end if
     end do
  end do
  write(*,'(a,es10.3,a,es9.1,a)') '(b) band coefficient vs direct quadrature: worst ', &
     worst_band, ' (tol ', TOL_BAND, ')'
  write(*,'(a,es10.3,a,es9.1,a)') '(c) Int E d alpha_1 = alpha_1 IP + beta_1: worst ', &
     worst_id, ' (tol ', TOL_IDENTITY, ')'
  write(*,'(a,es10.3,a,es9.1,a)') '(d) packets'' photon rate vs Int_band d alpha_1: worst ', &
     worst_rate, ' (tol ', TOL_RATE, ')'

  !--- (f) hydrogenic C VI
  do k = 1, milne_metal_count()
     call species_transition(k, ie, it)
     if (trim(elem_name(ie)) /= 'c' .or. it /= 6) cycle
     do itt = 1, size(TT)
        T = TT(itt)
        a1 = milne_metal_alpha1(k, T)
        ap = 3.0_wp*(alphaA_HeIII(T/9.0_wp) - alphaB_HeIII(T/9.0_wp))
        dev = a1/ap - 1.0_wp
        write(*,'(a,f9.0,a,es12.4,a,es12.4,a,es10.2)') '(f) C VI T =', T, &
           '  Milne ', a1, '  hydrogenic ', ap, '  dev ', dev
        if (abs(dev) > 2.0e-2_wp) ok = .false.
     end do
  end do

  !--- (e) report at 1e4 K
  T = 1.0e4_wp;  kT_eV = kboltz_cgs*T/ev2erg
  write(*,'(/,a)') '(e) at 1e4 K: ion  alpha_1  alpha_RR  alpha_rec(RR+DR)  ' // &
     'alpha_1/alpha_rec  cap  excited-term share  bound(excited captures)/P_ground'
  do k = 1, milne_metal_count()
     call species_transition(k, ie, it)
     a1   = milne_metal_alpha1(k, T)
     arr  = species_radiative_recombination(ie, it, T)
     atot = species_recombination_coefficient(ie, it, T)
     share = milne_metal_excited_share(k, T)
     ptot = milne_metal_energy_integral(k, T, elem_eth(ie,it), huge(1.0_wp))
     bound = 0.0_wp
     if (elem_eth(ie,it) > eth_HI) bound = max(atot - a1, 0.0_wp) &
        *(elem_eth(ie,it) + kT_eV)*ev2erg/max(ptot, tinest)
     write(ion,'(a,i0)') trim(elem_name(ie)), it
     write(*,'(3x,a8,3es11.3,f9.4,f8.4,f10.4,es11.3)') ion, a1, arr, atot, &
        a1/max(atot,tinest), milne_metal_capture_cap(k, T), share, bound
  end do

  if (ok) then
     write(*,'(/,a)') 'GATE (metal ground continua: cross section, band energy, ' // &
        'alpha_1 IP + beta_1, photon rate): PASS'
  else
     write(*,'(/,a)') 'GATE (metal ground continua): FAIL'
     call MPI_FINALIZE(ierr);  error stop 1
  end if
  call MPI_FINALIZE(ierr)

contains

  !--- the PHTABLE rows (eps [eV], sigma [cm^2]) of transition it of element
  !--- el, read from data/atomic/element_<el>.txt; n = 0 when it has none
  subroutine read_phtable(el, it, eps, sig, n)
    character(len=*), intent(in)  :: el
    integer,          intent(in)  :: it
    real(kind=wp),    intent(out) :: eps(:), sig(:)
    integer,          intent(out) :: n
    character(len=512) :: line
    character(len=32)  :: key
    integer :: unit, ios, cur, k
    n = 0;  cur = 0
    open(newunit=unit, file='data/atomic/element_'//trim(el)//'.txt', status='old')
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) key
       if (trim(key) == 'TRANSITION') read(line,*) key, cur
       if (trim(key) == 'PHTABLE' .and. cur == it) then
          read(line,*) key, n
          do k = 1, n
             read(unit,*) eps(k), sig(k)
          end do
          sig(1:n) = sig(1:n)*MB
          exit
       end if
    end do
    close(unit)
  end subroutine read_phtable

  !--- Int_band d alpha_1, the in-band photon rate per n_e n_+: the direct
  !--- energy quadrature over narrow slices, each divided by its central
  !--- photon energy.
  real(kind=wp) function band_photon_rate(k, ie, it, T) result(r)
    integer,       intent(in) :: k, ie, it
    real(kind=wp), intent(in) :: T
    integer, parameter :: NSL = 20000
    real(kind=wp) :: elo, ehi, ea, eb, em, pa
    integer :: i
    elo = max(elem_eth(ie,it), par%eion_min)
    ehi = min(par%eion_max, elem_eth(ie,it) + 60.0_wp*kboltz_cgs*T/ev2erg)
    r = 0.0_wp
    if (ehi <= elo) return
    !--- slices narrow enough that 1/E is constant to 1e-8 inside each
    do i = 1, NSL
       ea = elo + (ehi - elo)*real(i-1,wp)/real(NSL,wp)
       eb = elo + (ehi - elo)*real(i,wp)/real(NSL,wp)
       em = 0.5_wp*(ea + eb)
       pa = milne_metal_energy_integral(k, T, ea, eb)
       r  = r + pa/(em*ev2erg)
    end do
  end function band_photon_rate

end program check_metal_recombination_continuum
