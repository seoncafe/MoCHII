module ion_band_mod
!---------------------------------------------------------------------------
! MoCHII: ionizing frequency band + ionizing-photon generation.
!
! The band covers photon energies [par%eion_min, par%eion_max] eV with
! par%nnu_ion log-spaced bins (10-20 bins reach percent-level rate
! integrals).  The SED grid covers dust
! wavelengths and is untouched; the ionizing band carries its own grids,
! source sampling, opacity, and J tally.
!
! Source spectrum in the band: par%ion_spectrum (2-column file: E [eV],
! L_E [arb per eV]) or, when empty, a Planck function B_nu(par%tstar).
! Bin luminosities are integrated on a 32-point sub-grid per bin and
! normalized so that sum(ion_lum) = par%luminosity [erg/s]: par%luminosity
! is the luminosity OF THE BAND.  Packets sample bins from the luminosity
! CDF and carry Lpacket = par%luminosity / nphotons.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: ion_setup, gen_ion_photon
  public :: ion_e, ion_de, ion_nu, ion_dnu, ion_lum, ion_Ltot

  real(kind=wp), allocatable :: ion_e(:)     ! bin center [eV]
  real(kind=wp), allocatable :: ion_de(:)    ! bin width  [eV]
  real(kind=wp), allocatable :: ion_nu(:)    ! bin center [Hz]
  real(kind=wp), allocatable :: ion_dnu(:)   ! bin width  [Hz]
  real(kind=wp), allocatable :: ion_lum(:)   ! bin luminosity [erg/s]
  real(kind=wp), allocatable :: ion_cdf(:)   ! sampling CDF over bins
  !--- total band luminosity carried by the packets: par%luminosity
  !--- ([eion_min, eion_max]) plus, with par%add_fuv, the FUV part of the
  !--- same source spectrum.
  real(kind=wp), protected :: ion_Ltot = 0.0_wp

contains

  !=========================================================================
  subroutine ion_setup(metal_eth, nmetal)
    use mpi
    implicit none
    !--- metal_eth(1:nmetal): active metal photoionization thresholds [eV],
    !--- gathered by the caller from the species registry.  Only consulted
    !--- when par%ion_align_edges is on; absent for the H/He-only path.
    real(kind=wp), intent(in), optional :: metal_eth(:)
    integer,       intent(in), optional :: nmetal
    integer, parameter :: NSUB = 32
    real(kind=wp), allocatable :: eedge(:), meth(:)
    real(kind=wp) :: lo, hi, e1, e2, es, fsum, lion
    integer :: nnu, nfuv, i, k, ierr, nmet

    !--- FUV extension: nnu_fuv extra log bins on [efuv_min, eion_min]
    !--- BELOW the nnu_ion ionizing bins (par%nnu_ion becomes the total
    !--- bin count consumed by every band array).
    nfuv = 0
    if (par%add_fuv) then
       if (par%efuv_min >= par%eion_min) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: add_fuv requires par%efuv_min < par%eion_min.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       nfuv = par%nnu_fuv
       par%nnu_ion = par%nnu_ion + nfuv
    end if
    nnu = par%nnu_ion
    allocate(eedge(nnu+1), ion_e(nnu), ion_de(nnu), ion_nu(nnu), &
             ion_dnu(nnu), ion_lum(nnu), ion_cdf(nnu))

    lo = log(par%eion_min);  hi = log(par%eion_max)
    if (nfuv > 0) then
       do i = 1, nfuv                   ! FUV segment [efuv_min, eion_min)
          eedge(i) = exp(log(par%efuv_min) + (lo - log(par%efuv_min)) &
                     *real(i-1,wp)/real(nfuv,wp))
       end do
    end if
    if (par%ion_align_edges) then
       !--- ionizing segment [eion_min, eion_max]: bin edges pinned to the
       !--- ionization thresholds so no bin straddles one.
       nmet = 0
       if (present(nmetal))   nmet = nmetal
       if (nmet > 0 .and. present(metal_eth)) then
          allocate(meth(nmet));  meth = metal_eth(1:nmet)
       else
          allocate(meth(1));  meth = 0.0_wp;  nmet = 0
       end if
       call build_aligned_edges(eedge, nfuv, nnu, meth, nmet)
       deallocate(meth)
    else
       do i = nfuv+1, nnu+1             ! ionizing segment [eion_min, eion_max]
          eedge(i) = exp(lo + (hi - lo)*real(i-1-nfuv,wp)/real(nnu-nfuv,wp))
       end do
    end if
    do i = 1, nnu
       ion_e(i)   = sqrt(eedge(i)*eedge(i+1))          ! geometric bin center
       ion_de(i)  = eedge(i+1) - eedge(i)
       ion_nu(i)  = ion_e(i)  * ev2erg / h_planck_cgs
       ion_dnu(i) = ion_de(i) * ev2erg / h_planck_cgs
    end do

    !--- bin luminosities from the source spectrum (trapezoid on NSUB points).
    if (len_trim(par%ion_spectrum) > 0) then
       call bin_lum_from_file(eedge, nnu, NSUB)
    else
       if (par%tstar <= 0.0_wp) then
          if (mpar%p_rank == 0) write(*,'(a)') &
             'ERROR: use_ion_band requires par%tstar > 0 or par%ion_spectrum.'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       do i = 1, nnu
          e1 = eedge(i);  e2 = eedge(i+1)
          fsum = 0.0_wp
          do k = 1, NSUB
             es = e1 + (e2 - e1)*(real(k,wp) - 0.5_wp)/real(NSUB,wp)
             fsum = fsum + planck_nu(es, par%tstar)
          end do
          ion_lum(i) = fsum * (e2 - e1)/real(NSUB,wp)
       end do
    end if

    !--- normalize and build the sampling CDF.  par%luminosity is the
    !--- luminosity of the IONIZING segment [eion_min, eion_max]; with
    !--- add_fuv the FUV segment of the same spectrum rides on top, so
    !--- the ionizing photon budget (Q_H) is preserved without manual
    !--- band-ratio rescaling.
    lion = sum(ion_lum(nfuv+1:nnu))
    ion_lum  = ion_lum / lion * par%luminosity
    ion_Ltot = sum(ion_lum)
    ion_cdf(1) = ion_lum(1)
    do i = 2, nnu
       ion_cdf(i) = ion_cdf(i-1) + ion_lum(i)
    end do
    ion_cdf = ion_cdf / ion_cdf(nnu)

    if (mpar%p_rank == 0) then
       if (nfuv > 0) then
          write(*,'(a,i4,a,i4,a,f8.3,a,f8.3,a,f8.3,a)') ' ION: band: ', &
             nnu-nfuv, ' ionizing +', nfuv, ' FUV bins, ', par%efuv_min, &
             ' /', par%eion_min, ' -', par%eion_max, ' eV'
       else
          write(*,'(a,i4,a,f8.3,a,f8.3,a)') ' ION: band: ', nnu, ' bins, ', &
             par%eion_min, ' - ', par%eion_max, ' eV'
       end if
       if (len_trim(par%ion_spectrum) > 0) then
          write(*,'(2a)')        ' ION: spectrum file = ', trim(par%ion_spectrum)
       else
          write(*,'(a,es12.4)')  ' ION: Planck spectrum, tstar [K] = ', par%tstar
       end if
       write(*,'(a,es12.4)')     ' ION: ionizing luminosity [erg/s] = ', par%luminosity
       if (nfuv > 0) write(*,'(a,es12.4,a,f7.4)') &
          ' ION: band total with FUV        = ', ion_Ltot, &
          ',  L_FUV/L_ion = ', (ion_Ltot - par%luminosity)/par%luminosity
    end if
    deallocate(eedge)
  end subroutine ion_setup

  !=========================================================================
  ! Threshold-aligned ionizing bin edges.  Fills eedge(nfuv+1 .. nnu+1)
  ! (the ionizing segment [eion_min, eion_max]) so that every ionization
  ! threshold strictly inside the band falls ON a bin edge.  The ionizing
  ! segment carries nion = nnu - nfuv bins, split among the sub-segments
  ! between consecutive thresholds; each sub-segment is log-uniform.
  !=========================================================================
  subroutine build_aligned_edges(eedge, nfuv, nnu, metal_eth, nmetal)
    use mpi
    implicit none
    real(kind=wp), intent(inout) :: eedge(:)
    integer,       intent(in)    :: nfuv, nnu, nmetal
    real(kind=wp), intent(in)    :: metal_eth(:)
    real(kind=wp) :: cand(256), thr(256), points(256), w(64)
    real(kind=wp) :: tmp, lo, hi, wsum, frac
    integer :: ncand, nthr, nseg, nion, i, j, ipt, ib, r
    integer :: nb(64), off, best, ierr, ntot, itmp, ndrop
    integer :: cpri(256), tpri(256)
    logical :: dup
    real(kind=wp), parameter :: MERGE_TOL = 0.01_wp   ! |ln(T2/T1)| merge

    lo = par%eion_min;  hi = par%eion_max
    nion = nnu - nfuv

    !--- 1) collect candidate thresholds STRICTLY inside (eion_min, eion_max):
    !---    H I / He I / He II plus the active metal photoionization edges
    !---    (gathered by the caller from the species registry).  H/He carry
    !---    priority 1 so that a near-coincident metal edge does not displace
    !---    them in the merge (e.g. C II 24.383 vs He I 24.587 within 1%).
    ncand = 0
    call add_cand(cand, cpri, ncand, eth_HI,   1, lo, hi)
    call add_cand(cand, cpri, ncand, eth_HeI,  1, lo, hi)
    call add_cand(cand, cpri, ncand, eth_HeII, 1, lo, hi)
    do i = 1, nmetal
       if (metal_eth(i) > 0.0_wp) &
          call add_cand(cand, cpri, ncand, metal_eth(i), 0, lo, hi)
    end do

    !--- 2) sort ascending (insertion sort; ncand is small; carry priority)
    do i = 2, ncand
       tmp = cand(i);  itmp = cpri(i);  j = i - 1
       do while (j >= 1)
          if (cand(j) <= tmp) exit
          cand(j+1) = cand(j);  cpri(j+1) = cpri(j);  j = j - 1
       end do
       cand(j+1) = tmp;  cpri(j+1) = itmp
    end do

    !--- 2b) merge near-coincident thresholds.  Within a merged cluster keep
    !---    the higher-priority member (H/He over metal); ties keep the first.
    nthr = 0
    do i = 1, ncand
       dup = .false.
       if (nthr >= 1) then
          if (abs(log(cand(i)/thr(nthr))) < MERGE_TOL) dup = .true.
       end if
       if (dup) then
          if (cpri(i) > tpri(nthr)) then     ! promote H/He onto the edge
             thr(nthr)  = cand(i)
             tpri(nthr) = cpri(i)
          end if
       else
          nthr = nthr + 1
          thr(nthr)  = cand(i)
          tpri(nthr) = cpri(i)
       end if
    end do

    !--- 2c) graceful degradation: if the thresholds outnumber the ionizing
    !---    bins (nseg = nthr+1 > nion), drop the lowest-priority (metal)
    !---    thresholds - keeping H/He always - until they fit, removing the
    !---    most redundant metal first (smallest log-gap to a neighbor; the
    !---    band edges lo/hi are the outer neighbors).  So the aligned band is
    !---    safe at any nnu_ion: it resolves the big He edges plus as many
    !---    metal edges as the bins allow.
    ndrop = 0
    do while (nthr + 1 > nion)
       best = 0;  frac = 1.0e30_wp
       do i = 1, nthr
          if (tpri(i) /= 0) cycle            ! never drop H/He (priority 1)
          if (i == 1) then;  tmp = log(thr(i)/lo)
          else;              tmp = log(thr(i)/thr(i-1));  end if
          if (i == nthr) then;  tmp = min(tmp, log(hi/thr(i)))
          else;                 tmp = min(tmp, log(thr(i+1)/thr(i)));  end if
          if (tmp < frac) then;  frac = tmp;  best = i;  end if
       end do
       if (best == 0) exit                   ! only H/He remain - cannot reduce
       do i = best, nthr-1
          thr(i) = thr(i+1);  tpri(i) = tpri(i+1)
       end do
       nthr = nthr - 1;  ndrop = ndrop + 1
    end do
    if (ndrop > 0 .and. mpar%p_rank == 0) write(*,'(a,i0,a,i0,a)') &
       ' ION: align_edges: nnu_ion=', nion, ' too small for all edges; dropped ', &
       ndrop, ' metal threshold(s) (H/He kept; raise nnu_ion to align all).'

    !--- 3) split points: eion_min, interior thresholds, eion_max
    points(1) = lo
    do i = 1, nthr
       points(i+1) = thr(i)
    end do
    points(nthr+2) = hi
    nseg = nthr + 1

    !--- 5) need at least one bin per sub-segment.  After degradation this
    !---    only trips for an absurd nnu_ion below the H/He edge count.
    if (nion < nseg) then
       if (mpar%p_rank == 0) write(*,'(a,i0,a,i0,a,i0,a)') &
          'ERROR: ion_align_edges: nnu_ion=', nion, ' cannot fit ', nthr, &
          ' H/He thresholds; raise nnu_ion to >= ', nseg, '.'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if

    !--- 4) distribute nion bins over the nseg sub-segments by log width,
    !---    largest-remainder (Hamilton) with a floor of 1 per sub-segment.
    wsum = 0.0_wp
    do j = 1, nseg
       w(j) = log(points(j+1)) - log(points(j))
       wsum = wsum + w(j)
    end do
    ntot = 0
    do j = 1, nseg
       nb(j) = max(1, int(real(nion,wp)*w(j)/wsum))   ! floor(ideal), >= 1
       ntot  = ntot + nb(j)
    end do
    r = nion - ntot
    do while (r > 0)                 ! hand out by largest fractional remainder
       best = 1;  frac = -1.0_wp
       do j = 1, nseg
          tmp = real(nion,wp)*w(j)/wsum - real(nb(j),wp)
          if (tmp > frac) then
             frac = tmp;  best = j
          end if
       end do
       nb(best) = nb(best) + 1;  r = r - 1
    end do
    do while (r < 0)                 ! too many (floor-of-1 excess): remove
       best = 0;  frac = -1.0e30_wp
       do j = 1, nseg
          if (nb(j) <= 1) cycle
          tmp = real(nb(j),wp) - real(nion,wp)*w(j)/wsum   ! excess over ideal
          if (tmp > frac) then
             frac = tmp;  best = j
          end if
       end do
       if (best == 0) exit           ! cannot reduce further (all at 1)
       nb(best) = nb(best) - 1;  r = r + 1
    end do

    !--- 6) fill eedge inside each sub-segment (log-uniform); sub-segment
    !---    boundaries are the split points = thresholds.
    off = nfuv + 1                   ! eedge index of points(1) = eion_min
    eedge(off) = points(1)
    ipt = off
    do j = 1, nseg
       do ib = 1, nb(j)
          ipt = ipt + 1
          eedge(ipt) = exp(log(points(j)) + (log(points(j+1)) - log(points(j))) &
                       *real(ib,wp)/real(nb(j),wp))
       end do
    end do

    !--- 7) rank-0 log
    if (mpar%p_rank == 0) then
       write(*,'(a,i0,a,i0,a)') ' ION: align_edges on: ', nthr, &
          ' interior thresholds, ', nseg, ' sub-segments (bins each:'
       write(*,'(a,20(1x,i0))') '      ', (nb(j), j=1,nseg)
    end if
  end subroutine build_aligned_edges

  !=========================================================================
  ! Append E (with priority pri) to the candidate list if strictly inside
  ! (lo, hi).  Priority resolves near-coincident merges (H/He = 1 > metal 0).
  !=========================================================================
  subroutine add_cand(cand, cpri, ncand, E, pri, lo, hi)
    implicit none
    real(kind=wp), intent(inout) :: cand(:)
    integer,       intent(inout) :: cpri(:)
    integer,       intent(inout) :: ncand
    real(kind=wp), intent(in)    :: E, lo, hi
    integer,       intent(in)    :: pri
    if (E > lo .and. E < hi) then
       ncand = ncand + 1
       cand(ncand) = E
       cpri(ncand) = pri
    end if
  end subroutine add_cand

  !=========================================================================
  ! Planck B_nu at photon energy E [eV] and temperature T [K], arbitrary
  ! normalization per unit E (the constant prefactor cancels in ion_lum).
  !=========================================================================
  real(kind=wp) function planck_nu(E, T) result(B)
    real(kind=wp), intent(in) :: E, T
    real(kind=wp) :: x
    x = E * ev2erg / (kboltz_cgs * T)
    if (x < 700.0_wp) then
       B = E**3 / (exp(x) - 1.0_wp)
    else
       B = E**3 * exp(-x)          ! Wien tail, avoids overflow
    end if
  end function planck_nu

  !=========================================================================
  ! 2-column spectrum file (E [eV], L_E [arb per eV]; '#' comments), linear
  ! interpolation onto the sub-grid; zero outside the tabulated range.
  !=========================================================================
  subroutine bin_lum_from_file(eedge, nnu, nsub)
    use mpi
    implicit none
    real(kind=wp), intent(in) :: eedge(:)
    integer,       intent(in) :: nnu, nsub
    real(kind=wp), allocatable :: etab(:), ftab(:)
    character(len=256) :: line
    real(kind=wp) :: e1, e2, es, fsum, ei, fi
    integer :: unit, ios, ntab, i, k, j, ierr

    ntab = 0
    open(newunit=unit, file=trim(par%ion_spectrum), status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(2a)') 'ERROR: cannot open ', trim(par%ion_spectrum)
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       ntab = ntab + 1
    end do
    allocate(etab(ntab), ftab(ntab))
    rewind(unit)
    i = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       i = i + 1
       read(line,*) etab(i), ftab(i)
    end do
    close(unit)

    do i = 1, nnu
       e1 = eedge(i);  e2 = eedge(i+1)
       fsum = 0.0_wp
       do k = 1, nsub
          es = e1 + (e2 - e1)*(real(k,wp) - 0.5_wp)/real(nsub,wp)
          fi = 0.0_wp
          if (es >= etab(1) .and. es <= etab(ntab)) then
             do j = 1, ntab-1
                if (es <= etab(j+1)) then
                   ei = (es - etab(j)) / (etab(j+1) - etab(j))
                   fi = ftab(j)*(1.0_wp - ei) + ftab(j+1)*ei
                   exit
                end if
             end do
          end if
          fsum = fsum + fi
       end do
       ion_lum(i) = fsum * (e2 - e1)/real(nsub,wp)
    end do
    deallocate(etab, ftab)
  end subroutine bin_lum_from_file

  !=========================================================================
  ! Ionizing source packet: point source at par%xs/ys/zs_point, isotropic
  ! direction, frequency bin sampled from the band-luminosity CDF.
  !=========================================================================
  subroutine gen_ion_photon(photon)
    use random,     only : rand_number
    use octree_mod, only : amr_find_leaf
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp) :: cost, sint, phi, u
    integer :: i

    photon%x = par%xs_point
    photon%y = par%ys_point
    photon%z = par%zs_point

    cost = 2.0_wp*rand_number() - 1.0_wp
    sint = sqrt(1.0_wp - cost*cost)
    phi  = twopi*rand_number()
    photon%kx = sint*cos(phi)
    photon%ky = sint*sin(phi)
    photon%kz = cost

    u = rand_number()
    photon%inu = par%nnu_ion
    do i = 1, par%nnu_ion
       if (u <= ion_cdf(i)) then
          photon%inu = i
          exit
       end if
    end do

    photon%wgt     = 1.0_wp
    photon%Lpacket = ion_Ltot / real(par%nphotons, wp)
    photon%nscatt  = 0
    photon%inside  = .true.
    photon%icell_amr = amr_find_leaf(photon%x, photon%y, photon%z)
  end subroutine gen_ion_photon

end module ion_band_mod
