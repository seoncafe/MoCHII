program check_photo_shells
!---------------------------------------------------------------------------
! Metal photoionization cross sections against Verner's own phfit2.f.
!
! The VFKY96 outer-shell fit (phfit2's PH2 table) is the cross section of an
! ion only BELOW its first inner-shell edge; above it phfit2 sums the Verner
! & Yakovlev (1995) fits of every occupied subshell.  MoCHII used to carry
! the PH2 fit alone and extrapolate it through the Mg 2p (55-66 eV), Ca 3p
! (34-41 eV), Fe 3p (66-87 eV) and Cl 3s (~25 eV) edges, all inside the
! transported band -- Mg II at 85 eV came out 0.050 Mb against 6.68 Mb.
! species_cross_section now follows phfit2 exactly, and this gate holds it
! to the subroutine the data came from, summed over shells 1-7, for every
! registry ion on a log grid over 5-1000 eV.
!
! Energies within EDGE_SKIP of any threshold are skipped: MoCHII switches
! the outer shell on at the NIST ionization potential (the ETH row) where
! phfit2 uses its PH1 threshold, which differ by up to ~0.1 eV.
! phfit2 is single precision, hence TOL.
!
! Where an element file carries a PHTABLE (Ne III since 2026-10-11: the
! Nahar 2019 R-matrix ground-term cross section, which replaces the VFKY96
! fit below its last node E_join), the reference below E_join is that table,
! read here from the element file by this program and interpolated linearly
! in the photoelectron energy, with heat sigma (1 - E_th/E); phfit2 is the
! reference at and above E_join.  The table is also checked at every node
! and at the midpoint of every interval, and its last node against phfit2's
! PH2 fit (continuity at E_join; the PH1 threshold is not involved there).
!
! The photoheating coefficient is checked against its definition,
! sum_s sigma_s (1 - E_th,s/E), with the same phfit2 shells; the outer shell
! (PH1 threshold nearest the ionization potential) takes the NIST ionization
! potential of the ETH row as MoCHII does, the inner shells their PH1
! thresholds.
!---------------------------------------------------------------------------
  use define
  use species_mod, only : species_setup, species_cross_section, &
                          n_elements, elem_name, elem_nstage, elem_eth
  use mpi
  implicit none
  integer, parameter :: NEL = 11
  character(len=2), parameter :: ELN(NEL) = &
     ['c ','n ','o ','ne','s ','ar','mg','fe','si','cl','ca']
  integer, parameter :: ZEL(NEL) = [6,7,8,10,16,18,12,26,14,17,20]
  integer, parameter :: NE_GRID = 400
  real(kind=wp), parameter :: TOL = 2.0e-5_wp, EDGE_SKIP = 1.0e-2_wp
  real(kind=wp) :: E, sig, heat, sref, href, dev, worst, worsth
  real :: e4, s4, eth4(7)
  real(kind=wp), parameter :: MB = 1.0e-18_wp
  integer :: ie, it, k, is, nz, ne, ierr, j, ntest, nskip, ncontrib, is_outer
  real(kind=wp) :: eth_shell
  logical :: ok, near
  character(len=2) :: nm
  integer, parameter :: MAXTAB = 200000
  !--- the table is interpolated in eps = E - E_th, recovered from E ~ 70 eV:
  !--- a relative rounding of ~1e-16 in E is ~1e-14 in eps and, on the steep
  !--- side of a resonance (1.36 meV nodes), up to ~1e-10 in sigma
  real(kind=wp), parameter :: TOL_TABLE = 1.0e-8_wp
  real(kind=wp), allocatable :: teps(:), tsig(:)
  integer :: ntab
  real(kind=wp) :: ejoin
  real :: ph1
  integer :: ninn, ntot
  common/ph1/ph1(6,30,30,7)
  common/ninn/ninn(30)
  common/ntot/ntot(30)

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir = 'data/atomic'
  par%use_metals = .true.
  par%abund_C  = 1.0e-4_wp;  par%abund_N  = 1.0e-4_wp
  par%abund_O  = 1.0e-4_wp;  par%abund_Ne = 1.0e-4_wp
  par%abund_S  = 1.0e-5_wp;  par%abund_Ar = 1.0e-5_wp
  par%abund_Mg = 1.0e-5_wp;  par%abund_Fe = 1.0e-5_wp
  par%abund_Si = 1.0e-5_wp;  par%abund_Cl = 1.0e-7_wp
  par%abund_Ca = 1.0e-8_wp
  !--- no (T, n_e) line-cooling tables: the gate does not use the cooling
  par%cooling_model = 'low_density'
  call species_setup(1)
  allocate(teps(MAXTAB), tsig(MAXTAB))

  ok = .true.;  worst = 0.0_wp;  worsth = 0.0_wp;  ntest = 0;  nskip = 0
  do ie = 1, n_elements
     nm = elem_name(ie)
     nz = 0
     do k = 1, NEL
        if (trim(ELN(k)) == trim(nm)) nz = ZEL(k)
     end do
     do it = 1, elem_nstage(ie) - 1
        ne = nz - it + 1
        do is = 1, 7
           eth4(is) = ph1(1, nz, ne, is)
        end do
        !--- the outer shell is the one whose PH1 threshold is the ionization
        !--- potential (phfit2 carries empty shells too, e.g. Ca II 3d at
        !--- 1.0 eV, so the lowest threshold is not it)
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
           call check_table_nodes(ie, it, nm, ok, worst, worsth, ntest)
           !--- continuity at E_join: the last node against phfit2's PH2 fit
           e4 = real(ejoin)
           call phfit2(nz, ne, is_outer, e4, s4)
           dev = abs(tsig(ntab)/(real(s4,wp)*MB) - 1.0_wp)
           write(*,'(a,a2,i2,a,f10.4,a,es10.2)') 'PHTABLE ', nm, it, ': E_join = ', &
              ejoin, ' eV, last node / phfit2 PH2 - 1 = ', dev
           if (dev > TOL) ok = .false.
        end if
        do j = 0, NE_GRID
           E = 5.0_wp*(1000.0_wp/5.0_wp)**(real(j,wp)/NE_GRID)
           if (E < ejoin) cycle             ! on the table: checked above
           near = .false.
           do is = 1, 7
              if (eth4(is) > 0.0 .and. abs(E/eth4(is) - 1.0_wp) < EDGE_SKIP) near = .true.
           end do
           if (near) then;  nskip = nskip + 1;  cycle;  end if
           call species_cross_section(ie, it, E, sig, heat)
           e4 = real(E);  sref = 0.0_wp;  href = 0.0_wp;  ncontrib = 0
           do is = 1, 7
              call phfit2(nz, ne, is, e4, s4)
              sref = sref + real(s4,wp)*MB
              if (s4 > 0.0) then
                 ncontrib = ncontrib + 1
                 eth_shell = real(eth4(is),wp)
                 if (is == is_outer) eth_shell = elem_eth(ie,it)
                 href = href + real(s4,wp)*MB*max(1.0_wp - eth_shell/E, 0.0_wp)
              end if
           end do
           !--- one contributing shell = phfit2's PH2 branch (or the outer
           !--- shell alone), where MoCHII takes the photoelectron energy
           !--- from the NIST ionization potential of its ETH row.
           if (ncontrib == 1) href = sref*max(1.0_wp - elem_eth(ie,it)/E, 0.0_wp)
           if (sref <= 0.0_wp .and. sig <= 0.0_wp) cycle
           dev = abs(sig - sref)/max(sref, 1.0e-30_wp)
           ntest = ntest + 1
           worst = max(worst, dev)
           if (dev > TOL) then
              ok = .false.
              write(*,'(a,a2,i2,f10.3,2es14.6,es10.2)') '  FAIL sigma ', nm, it, E, sig, sref, dev
           end if
           dev = abs(heat - href)/max(href, 1.0e-30_wp)
           worsth = max(worsth, dev)
           if (dev > TOL) then
              ok = .false.
              write(*,'(a,a2,i2,f10.3,2es14.6,es10.2)') '  FAIL heat  ', nm, it, E, heat, href, dev
           end if
        end do
     end do
  end do
  write(*,'(a,i0,a,i0,a)') 'compared ', ntest, ' (ion, energy) points, skipped ', nskip, ' near edges'
  write(*,'(a,es10.3,a,es10.3)') 'worst |sigma/phfit2 - 1| = ', worst, '   worst heating = ', worsth
  if (ok) then
     write(*,'(a)') 'GATE (metal photoionization = phfit2 shell sum): PASS'
  else
     write(*,'(a)') 'GATE (metal photoionization = phfit2 shell sum): FAIL'
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

  !--- species_cross_section on the table: at every node and at the midpoint
  !--- of every interval, against linear interpolation of the rows
  subroutine check_table_nodes(ie, it, nm, ok, worst, worsth, ntest)
    integer,          intent(in)    :: ie, it
    character(len=2), intent(in)    :: nm
    logical,          intent(inout) :: ok
    real(kind=wp),    intent(inout) :: worst, worsth
    integer,          intent(inout) :: ntest
    real(kind=wp) :: E, sig, heat, sref, href, dev
    integer :: k, h
    do k = 1, ntab - 1
       do h = 0, 1
          E = elem_eth(ie,it) + teps(k) + 0.5_wp*real(h,wp)*(teps(k+1) - teps(k))
          sref = tsig(k) + 0.5_wp*real(h,wp)*(tsig(k+1) - tsig(k))
          href = sref*max(1.0_wp - elem_eth(ie,it)/E, 0.0_wp)
          call species_cross_section(ie, it, E, sig, heat)
          ntest = ntest + 1
          dev = abs(sig - sref)/max(sref, 1.0e-30_wp)
          if (sref <= 0.0_wp) dev = abs(sig)/1.0e-30_wp
          worst = max(worst, dev)
          if (dev > TOL_TABLE) then
             ok = .false.
             write(*,'(a,a2,i2,f12.6,2es14.6,es10.2)') '  FAIL table sigma ', nm, it, E, sig, sref, dev
          end if
          if (href > 0.0_wp) then
             dev = abs(heat - href)/href
             worsth = max(worsth, dev)
             if (dev > TOL_TABLE) then
                ok = .false.
                write(*,'(a,a2,i2,f12.6,2es14.6,es10.2)') '  FAIL table heat  ', nm, it, E, heat, href, dev
             end if
          end if
       end do
    end do
  end subroutine check_table_nodes

end program check_photo_shells
