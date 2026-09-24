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
! The photoheating coefficient is checked against its definition,
! sum_s sigma_s (1 - E_th,s/E), with the same phfit2 shells.
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
  integer :: ie, it, k, is, nz, ne, ierr, j, ntest, nskip, ncontrib
  logical :: ok, near
  character(len=2) :: nm
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
  call species_setup(1)

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
        do j = 0, NE_GRID
           E = 5.0_wp*(1000.0_wp/5.0_wp)**(real(j,wp)/NE_GRID)
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
                 href = href + real(s4,wp)*MB* &
                        max(1.0_wp - real(eth4(is),wp)/E, 0.0_wp)
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
end program check_photo_shells
