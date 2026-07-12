module ion_band_mod
!---------------------------------------------------------------------------
! MoCHII: ionizing frequency band + ionizing-photon generation (new, G0).
!
! The band covers photon energies [par%eion_min, par%eion_max] eV with
! par%nnu_ion log-spaced bins (10-20 bins reach percent-level rate
! integrals; docs/PLAN.md section 2.1).  The SED grid covers dust
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
  subroutine ion_setup()
    use mpi
    implicit none
    integer, parameter :: NSUB = 32
    real(kind=wp), allocatable :: eedge(:)
    real(kind=wp) :: lo, hi, e1, e2, es, fsum, lion
    integer :: nnu, nfuv, i, k, ierr

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
    do i = nfuv+1, nnu+1                ! ionizing segment [eion_min, eion_max]
       eedge(i) = exp(lo + (hi - lo)*real(i-1-nfuv,wp)/real(nnu-nfuv,wp))
    end do
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
