module sedust_mod
!---------------------------------------------------------------------------
! MoCHII: stochastic dust emission via the SEDust library (adapted from
! MoCafe_v2.00/src/dustemis_mod.f90; libsedust.a + .mod copied to
! SEDust_lib/, provenance 2026-07-12).
!
! Per leaf, the local mean intensity built from the transported band tally
! is handed to SEDust (equilibrium + stochastically heated grains + PAHs;
! astrodust / DL07 / Zubko models) for the emission SPECTRAL SHAPE; the
! absolute luminosity is the locally absorbed power Heat_dust*V
! (radiative equilibrium — exact by construction, sidestepping the
! SEDust normalization convention, as in MoCafe).  Because HII regions
! are optically thin to the re-emitted IR, no Lucy re-iteration is
! needed and the grid-integrated SED is accumulated directly.
!
! Field mapping: the band covers lambda = 1.24/E_max ... 1.24/E_min um;
! SEDust's optics grid starts at ~0.0912 um, so the EUV part of J below
! its grid minimum is deposited energy-conservingly into the shortest
! SEDust bin (the stochastic spike hardness is slightly underestimated
! there — documented approximation).
!
! SEDust reads its dielectric tables relative to its sed/ directory:
! par%sed_workdir must point at a SEDust sed/ tree (e.g. the MoCafe copy,
! read-only) and par%sed_qtable / par%sed_sizedist at the optics tables.
!---------------------------------------------------------------------------
  use define
  use dust_lib, only : dust_model_t, build_astrodust, build_dl07, &
                       build_zubko, dust_emission, dust_nlam, dust_lambda
  implicit none
  private

  public :: sedust_setup, sedust_compute_write

  type(dust_model_t)         :: dmodel
  integer                    :: nl_sed = 0
  real(kind=wp), allocatable :: lam_sed(:)

contains

  !=========================================================================
  subroutine sedust_setup()
    use mpi
    use ifport, only : chdir, getcwd
    implicit none
    integer :: ierr, cstat
    character(len=512) :: cwd_save

    cstat = getcwd(cwd_save)
    if (len_trim(par%sed_workdir) > 0) then
       cstat = chdir(trim(par%sed_workdir))
       if (cstat /= 0 .and. mpar%p_rank == 0) write(*,'(3a)') &
          'WARNING: could not chdir to par%sed_workdir = ''', &
          trim(par%sed_workdir), ''''
    end if

    select case (trim(par%dust_model_sed))
    case ('astrodust')
       call build_astrodust(dmodel, trim(par%sed_qtable), &
            trim(par%sed_sizedist), par%sed_NT, par%sed_Tlo, par%sed_Thi)
    case ('dl07')
       call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
            par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
            par%sed_NT, par%sed_Tlo, par%sed_Thi)
    case ('zubko')
       call build_zubko(dmodel, trim(par%sed_zubko_config), &
            trim(par%sed_zubko_dir), par%sed_NT, par%sed_Tlo, par%sed_Thi)
    case default
       cstat = chdir(trim(cwd_save))
       if (mpar%p_rank == 0) write(*,'(3a)') 'ERROR: dust_model_sed = ''', &
          trim(par%dust_model_sed), ''' (astrodust/dl07/zubko).'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end select
    cstat = chdir(trim(cwd_save))

    nl_sed  = dust_nlam(dmodel)
    lam_sed = dust_lambda(dmodel)
    if (mpar%p_rank == 0) then
       write(*,'(a)')    ' SEDU: SEDust dust-emission model ready'
       write(*,'(2a)')   ' SEDU: model  = ', trim(dmodel%name)
       write(*,'(a,i6,a,2es11.3)') ' SEDU: grid   = ', nl_sed, &
          ' points, lambda range [um] = ', lam_sed(1), lam_sed(nl_sed)
    end if
  end subroutine sedust_setup

  !=========================================================================
  ! Grid-integrated stochastic dust SED.  heat_dust [erg/s/cm^3] sets each
  ! leaf's absolute emission; SEDust sets the shape from the local field.
  !=========================================================================
  subroutine sedust_compute_write(heat_dust)
    use mpi
    use octree_mod,      only : amr_grid
    use jtally_mod,      only : jt_ion
    use ion_band_mod,    only : ion_e, ion_dnu
    use utility,         only : get_base_name, is_finite
    implicit none
    real(kind=wp), intent(in) :: heat_dust(:)
    real(kind=wp), allocatable :: Ltot(:), Jsed(:), lamI(:), pdf(:)
    real(kind=wp), allocatable :: lam_b(:), Jlam_b(:)
    real(kind=wp), allocatable :: em_band(:,:)
    real(kind=wp) :: band_wl(8), vol, Labs, esum, dlam1, extra
    character(len=192) :: outname
    integer :: il, ic, k, b, unit, ierr, ndone, nmine, nband, ib

    allocate(Ltot(nl_sed), Jsed(nl_sed), lamI(nl_sed), pdf(nl_sed))
    allocate(lam_b(par%nnu_ion), Jlam_b(par%nnu_ion))
    Ltot = 0.0_wp
    !--- dust-band leaf emissivities (par%dust_emis_bands, um).
    nband = count(is_finite(par%dust_emis_bands))
    if (nband > 0) then
       band_wl(1:nband) = pack(par%dust_emis_bands, &
                               is_finite(par%dust_emis_bands))
       allocate(em_band(nband, amr_grid%nleaf))
       em_band = 0.0_wp
    end if
    do b = 1, par%nnu_ion
       lam_b(b) = 1.23984_wp/ion_e(b)          ! um (descending in b)
    end do
    dlam1 = lam_sed(2) - lam_sed(1)

    nmine = (amr_grid%nleaf - mpar%p_rank + mpar%nproc - 1)/mpar%nproc
    ndone = 0
    if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
       ' SEDU: solving stochastic emission for ~', nmine, ' leaves each rank'

    do il = mpar%p_rank+1, amr_grid%nleaf, mpar%nproc
       ndone = ndone + 1
       if (mpar%p_rank == 0 .and. nmine >= 500 .and. &
           mod(ndone, nmine/5) == 0) write(*,'(a,i0,a,i0)') &
          ' SEDU: rank0 ', ndone, '/', nmine
       if (heat_dust(il) <= 0.0_wp) cycle
       ic  = amr_grid%icell_of_leaf(il)
       vol = (2.0_wp*amr_grid%ch(ic)*par%distance2cm)**3
       Labs = heat_dust(il)*vol
       !--- J_lambda [SI, W/m^2/sr/m] of this leaf on the band bins.
       do b = 1, par%nnu_ion
          !--- J_nu = jt/(4 pi V_code d_cm^2 dnu); J_lambda = J_nu c/lambda^2
          Jlam_b(b) = jt_ion(b,il)/(fourpi*(vol/par%distance2cm**3) &
                      *par%distance2cm**2*ion_dnu(b)) &
                      * (2.99792458e14_wp/lam_b(b)**2)      ! cgs per um
          Jlam_b(b) = Jlam_b(b)*1.0e3_wp                     ! -> SI per m
       end do
       !--- map onto the SEDust grid: interpolate inside the overlap,
       !--- deposit the below-grid EUV energy into the first bin.
       Jsed = 0.0_wp
       extra = 0.0_wp
       do b = 1, par%nnu_ion
          if (lam_b(b) < lam_sed(1)) then
             !--- energy flux of this bin [SI]: J_lambda dlambda with
             !--- dlambda = lambda^2 dnu / c
             extra = extra + Jlam_b(b)*(lam_b(b)**2*ion_dnu(b) &
                     /2.99792458e14_wp)*1.0e-6_wp            ! m
          end if
       end do
       do k = 1, nl_sed
          if (lam_sed(k) > lam_b(1)) exit                    ! beyond band max
          if (lam_sed(k) < lam_b(par%nnu_ion)) cycle         ! below band min handled above
          Jsed(k) = interp_band(lam_b, Jlam_b, lam_sed(k))
       end do
       Jsed(1) = Jsed(1) + extra/(dlam1*1.0e-6_wp)
       !--- SEDust emission shape, normalized to the absorbed power.
       call dust_emission(dmodel, Jsed, lamI)
       do k = 1, nl_sed
          pdf(k) = max(lamI(k), 0.0_wp)/lam_sed(k)
       end do
       esum = 0.0_wp
       do k = 1, nl_sed-1
          esum = esum + 0.5_wp*(pdf(k) + pdf(k+1))*(lam_sed(k+1) - lam_sed(k))
       end do
       if (esum <= 0.0_wp) cycle
       do k = 1, nl_sed
          Ltot(k) = Ltot(k) + Labs*pdf(k)/esum               ! L_lambda [erg/s/um]
       end do
       !--- band emissivities of this leaf [erg/s/cm^3/um].
       do ib = 1, nband
          em_band(ib, il) = Labs*pdf_interp(band_wl(ib))/esum/vol
       end do
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, Ltot, nl_sed, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    if (nband > 0) then
       call MPI_ALLREDUCE(MPI_IN_PLACE, em_band, size(em_band), &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       call dustemis_write(em_band, band_wl, nband)
       deallocate(em_band)
    end if

    if (mpar%p_rank == 0) then
       esum = 0.0_wp
       do k = 1, nl_sed-1
          esum = esum + 0.5_wp*(Ltot(k) + Ltot(k+1))*(lam_sed(k+1) - lam_sed(k))
       end do
       outname = trim(get_base_name(par%out_file))//'_dustsed.txt'
       open(newunit=unit, file=trim(outname), status='replace')
       write(unit,'(a)') '# MoCHII stochastic dust emission (SEDust '// &
          'shape x locally absorbed power; optically thin in the IR)'
       write(unit,'(3a)') '# model: ', trim(dmodel%name), &
          '; EUV J below the SEDust grid minimum deposited in the first bin'
       write(unit,'(a,es14.6,a)') '# integral L_lambda dlambda = ', esum, &
          ' erg/s (equals sum Heat_dust V by construction)'
       write(unit,'(a)') '# lambda[um]   L_lambda[erg/s/um]'
       do k = 1, nl_sed
          write(unit,'(es13.5,es15.6)') lam_sed(k), Ltot(k)
       end do
       close(unit)
       write(*,'(2a)') ' SEDU: dust SED written to: ', trim(outname)
       write(*,'(a,es12.4,a)') ' SEDU: total dust luminosity = ', esum, ' erg/s'
    end if
    deallocate(Ltot, Jsed, lamI, pdf, lam_b, Jlam_b)

  contains

    !--- linear interpolation of the current leaf's emission pdf at
    !--- wavelength w [um] (0 outside the SEDust grid).
    real(kind=wp) function pdf_interp(w) result(p)
      real(kind=wp), intent(in) :: w
      real(kind=wp) :: f
      integer :: j
      p = 0.0_wp
      if (w <= lam_sed(1) .or. w >= lam_sed(nl_sed)) return
      do j = 1, nl_sed-1
         if (w < lam_sed(j+1)) exit
      end do
      f = (w - lam_sed(j))/(lam_sed(j+1) - lam_sed(j))
      p = pdf(j)*(1.0_wp - f) + pdf(j+1)*f
    end function pdf_interp
  end subroutine sedust_compute_write

  !=========================================================================
  ! Dust-band leaf emissivities -> '<base>_dustemis' (HDF5/FITS): the
  ! same block layout the Python EmisData reader consumes (emis_dust +
  ! wl_dust + LeafXYZ + LeafSize).  The IR is optically thin, so a
  ! flux-conserving column map of these blocks IS the dust-band image.
  !=========================================================================
  subroutine dustemis_write(em_band, band_wl, nband)
    use octree_mod, only : amr_grid
    use iofile_mod
    use utility,    only : get_base_name
    implicit none
    real(kind=wp), intent(in) :: em_band(:,:), band_wl(:)
    integer,       intent(in) :: nband
    type(io_file_type) :: file
    character(len=192) :: outname
    real(kind=wp), allocatable :: lxyz(:,:), tmp(:)
    real(kind=wp) :: wl_A(nband)
    integer :: il, ic, status

    if (mpar%p_rank /= 0) return
    status = 0
    outname = trim(get_base_name(par%out_file))//'_dustemis'// &
              trim(io_file_extension(par%file_format))
    call io_open_new(file, trim(outname), status)

    allocate(lxyz(amr_grid%nleaf,3), tmp(amr_grid%nleaf))
    do il = 1, amr_grid%nleaf
       ic = amr_grid%icell_of_leaf(il)
       lxyz(il,1) = amr_grid%cx(ic)
       lxyz(il,2) = amr_grid%cy(ic)
       lxyz(il,3) = amr_grid%cz(ic)
       tmp(il)    = 2.0_wp*amr_grid%ch(ic)
    end do
    call io_append_image(file, lxyz, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','LeafXYZ', &
         'leaf center x,y,z (code units)',status)
    call io_put_keyword(file,'DIST_CM', par%distance2cm, &
         'distance unit (cm)', status)
    call io_append_image(file, tmp, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','LeafSize', &
         'leaf edge length (code units)',status)
    call io_append_image(file, em_band(1:nband,:), status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','emis_dust', &
         'dust band emissivity [erg/s/cm^3/um]',status)
    wl_A = band_wl(1:nband)*1.0e4_wp
    call io_append_image(file, wl_A, status, bitpix=-64)
    call io_put_keyword(file,'EXTNAME','wl_dust', &
         'band wavelengths [A]',status)
    call io_close(file, status)
    deallocate(lxyz, tmp)
    write(*,'(2a)') ' SEDU: dust band emissivities written to: ', trim(outname)
  end subroutine dustemis_write

  !=========================================================================
  ! Interpolate J_lambda from the (descending-lambda) band bins.
  !=========================================================================
  real(kind=wp) function interp_band(lam_b, J_b, l0) result(j)
    implicit none
    real(kind=wp), intent(in) :: lam_b(:), J_b(:), l0
    integer :: b, n
    real(kind=wp) :: w
    n = size(lam_b)
    j = 0.0_wp
    if (l0 > lam_b(1) .or. l0 < lam_b(n)) return
    do b = 1, n-1
       if (l0 <= lam_b(b) .and. l0 >= lam_b(b+1)) then
          w = (lam_b(b) - l0)/(lam_b(b) - lam_b(b+1))
          j = J_b(b)*(1.0_wp - w) + J_b(b+1)*w
          return
       end if
    end do
  end function interp_band

end module sedust_mod
