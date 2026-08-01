module sedust_mod
!---------------------------------------------------------------------------
! MoCHII: stochastic dust emission via the SEDust dust_lib layer, built from
! source under SEDust/sed (build_lib.sh); a self-contained copy of the
! MoCafe v2.00 SEDust tree.  Adapted from MoCafe_v2.00/src/dustemis_mod.f90.
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
! Field mapping: the band covers lambda = 1.24/E_max ... 1.24/E_min um.
! With par%ion_add_dust the shared grain model (grain_model_mod) is built
! down to 1.24/eion_max, so the whole EUV band falls on the model grid and
! the grains absorb it with their own cross sections.  In an emission-only
! run (ion_add_dust off) the model keeps its native grid, which starts at
! ~0.0912 um, and the EUV part of J below that minimum is deposited
! energy-conservingly into the shortest SEDust bin (the stochastic spike
! hardness is slightly underestimated there — documented approximation).
!
! SEDust reads its dielectric tables relative to its sed/ directory:
! par%sed_workdir must point at a SEDust sed/ tree (e.g. the MoCafe copy,
! read-only) and par%sed_qtable / par%sed_sizedist at the optics tables.
!---------------------------------------------------------------------------
  use define
  use dust_lib,        only : dust_model_t, dust_emission, dust_nlam, dust_lambda
  !--- the grain population model is shared with gas_opacity_mod / dust_temp_mod
  !--- (grain_model_mod): the absorbed power and the reemitted spectrum then
  !--- refer to the SAME grains.  gas_opacity_setup builds it first when the
  !--- ionizing band carries dust, so build_grain_model here is a no-op reuse.
  use grain_model_mod, only : dmodel, build_grain_model
  implicit none
  private

  public :: sedust_setup, sedust_compute_write
  public :: sedust_leaf_spectrum, sedust_lambda_grid

  integer                    :: nl_sed = 0
  real(kind=wp), allocatable :: lam_sed(:)
  !--- index of the model's 'PAH' output channel (0 = none; the
  !--- par%sed_pah_live weighting needs it).
  integer                    :: ipah_chan = 0
  !--- scratch of one leaf's solve, allocated once: the band-bin wavelengths
  !--- and J_lambda, the same field mapped onto the SEDust grid, the emitted
  !--- lambda I_lambda and its per-channel split.
  real(kind=wp), allocatable :: lam_b(:), Jlam_b(:), Jsed(:), lamI(:)
  real(kind=wp), allocatable :: lamI_ch(:,:)
  real(kind=wp)              :: dlam1 = 0.0_wp

contains

  !=========================================================================
  subroutine sedust_setup()
    use ion_band_mod, only : ion_e, nnu_band
    implicit none
    integer :: b

    !--- idempotent: the emission transport builds the model before the
    !--- output pass asks for it again.
    if (nl_sed > 0) return

    !--- build (or reuse) the shared SEDust grain model.
    call build_grain_model()

    nl_sed  = dust_nlam(dmodel)
    lam_sed = dust_lambda(dmodel)

    allocate(lam_b(nnu_band), Jlam_b(nnu_band))
    allocate(Jsed(nl_sed), lamI(nl_sed))
    do b = 1, nnu_band
       lam_b(b) = 1.23984_wp/ion_e(b)          ! um (descending in b)
    end do
    dlam1 = lam_sed(2) - lam_sed(1)

    !--- locate the 'PAH' output channel for par%sed_pah_live.
    ipah_chan = 0
    if (par%sed_pah_live) then
       block
         integer :: ic
         do ic = 1, dmodel%n_channel
            if (trim(dmodel%channel_name(ic)) == 'PAH') ipah_chan = ic
         end do
         if (ipah_chan == 0 .and. mpar%p_rank == 0) write(*,'(a)') &
            ' SEDU: WARNING: sed_pah_live needs a model with a ''PAH'' '// &
            'channel (astrodust); weighting disabled.'
       end block
    end if

    if (mpar%p_rank == 0) then
       write(*,'(a)')    ' SEDU: SEDust dust-emission model ready'
       write(*,'(2a)')   ' SEDU: model  = ', trim(dmodel%name)
       write(*,'(a,i6,a,2es11.3)') ' SEDU: grid   = ', nl_sed, &
          ' points, lambda range [um] = ', lam_sed(1), lam_sed(nl_sed)
       if (ipah_chan > 0) write(*,'(a,i2,a)') &
          ' SEDU: live PAH weighting on channel', ipah_chan, ' (''PAH'')'
    end if
  end subroutine sedust_setup

  !=========================================================================
  ! The SEDust model's wavelength grid [um].
  !=========================================================================
  function sedust_lambda_grid() result(lam)
    implicit none
    real(kind=wp), allocatable :: lam(:)
    lam = lam_sed
  end function sedust_lambda_grid

  !=========================================================================
  ! Emission spectrum of ONE leaf: the local mean intensity built from the
  ! transported band tally is handed to SEDust, which returns lambda
  ! I_lambda; pdf(k) = lamI(k)/lam_sed(k) is then proportional to L_lambda
  ! and esum = Int pdf dlambda is its trapezoid integral, so the leaf's
  ! absolutely normalized emission is L_lambda(k) = Labs*pdf(k)/esum with
  ! Labs = Heat_dust*V.  esum <= 0 means this leaf emits nothing.
  !
  ! Single source of the SEDust leaf spectrum: the grid-integrated SED output
  ! and the emitted-photon launch both come through here.
  !=========================================================================
  subroutine sedust_leaf_spectrum(il, pdf, esum)
    use octree_mod,    only : leaf_half
    use jtally_mod,    only : jt_ion
    use ion_band_mod,  only : ion_dnu, nnu_band
    use gas_state_mod, only : gas_xHI
    implicit none
    integer,       intent(in)  :: il
    real(kind=wp), intent(out) :: pdf(:), esum
    real(kind=wp) :: vol, extra, xh, wd, wpah
    integer :: b, k

    pdf  = 0.0_wp
    esum = 0.0_wp
    vol  = (2.0_wp*leaf_half(il)*par%distance2cm)**3
    !--- J_lambda [SI, W/m^2/sr/m] of this leaf on the band bins.
    do b = 1, nnu_band
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
    do b = 1, nnu_band
       if (lam_b(b) < lam_sed(1)) then
          !--- energy flux of this bin [SI]: J_lambda dlambda with
          !--- dlambda = lambda^2 dnu / c
          extra = extra + Jlam_b(b)*(lam_b(b)**2*ion_dnu(b) &
                  /2.99792458e14_wp)*1.0e-6_wp            ! m
       end if
    end do
    do k = 1, nl_sed
       if (lam_sed(k) > lam_b(1)) exit                    ! beyond band max
       if (lam_sed(k) < lam_b(nnu_band)) cycle         ! below band min handled above
       Jsed(k) = interp_band(lam_b, Jlam_b, lam_sed(k))
    end do
    Jsed(1) = Jsed(1) + extra/(dlam1*1.0e-6_wp)
    !--- SEDust emission shape, normalized to the absorbed power.
    !--- With sed_pah_live the PAH channel is weighted by the leaf's
    !--- PAH survival (xHI + f_ion_pah xHII), divided under
    !--- laursen09_live by the dust survival already in rhokap.
    if (ipah_chan > 0) then
       if (.not. allocated(lamI_ch)) &
          allocate(lamI_ch(nl_sed, dmodel%n_channel))
       call dust_emission(dmodel, Jsed, lamI, lamI_ch)
       xh   = gas_xHI(il)
       wpah = xh + par%f_ion_pah*(1.0_wp - xh)
       wd   = 1.0_wp
       if (trim(par%dust_model) == 'laursen09_live') &
          wd = max(xh + par%f_ion_dust*(1.0_wp - xh), tinest)
       lamI = lamI + (wpah/wd - 1.0_wp)*lamI_ch(:, ipah_chan)
    else
       call dust_emission(dmodel, Jsed, lamI)
    end if
    do k = 1, nl_sed
       pdf(k) = max(lamI(k), 0.0_wp)/lam_sed(k)
    end do
    esum = 0.0_wp
    do k = 1, nl_sed-1
       esum = esum + 0.5_wp*(pdf(k) + pdf(k+1))*(lam_sed(k+1) - lam_sed(k))
    end do
  end subroutine sedust_leaf_spectrum

  !=========================================================================
  ! Grid-integrated stochastic dust SED.  heat_dust [erg/s/cm^3] sets each
  ! leaf's absolute emission; SEDust sets the shape from the local field.
  !=========================================================================
  subroutine sedust_compute_write(heat_dust)
    use mpi
    use octree_mod, only : amr_grid, leaf_half
    use utility,         only : get_base_name, is_finite
    implicit none
    real(kind=wp), intent(in) :: heat_dust(:)
    real(kind=wp), allocatable :: Ltot(:), pdf(:)
    real(kind=wp) :: band_wl(8), vol, Labs, esum
    real(kind=wp), allocatable :: em_band(:,:)
    character(len=192) :: outname
    integer :: il, ic, k, unit, ierr, ndone, nmine, nband, ib

    allocate(Ltot(nl_sed), pdf(nl_sed))
    Ltot = 0.0_wp
    !--- dust-band leaf emissivities (par%dust_emis_bands, um).
    nband = count(is_finite(par%dust_emis_bands))
    if (nband > 0) then
       band_wl(1:nband) = pack(par%dust_emis_bands, &
                               is_finite(par%dust_emis_bands))
       allocate(em_band(nband, amr_grid%nleaf))
       em_band = 0.0_wp
    end if

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
       vol = (2.0_wp*leaf_half(il)*par%distance2cm)**3
       Labs = heat_dust(il)*vol
       call sedust_leaf_spectrum(il, pdf, esum)
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
    deallocate(Ltot, pdf)

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
    use octree_mod, only : amr_grid, leaf_half, leaf_cx, leaf_cy, leaf_cz
    use iofile_mod
    use utility,    only : get_base_name, fatal_error
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
       lxyz(il,1) = leaf_cx(il)
       lxyz(il,2) = leaf_cy(il)
       lxyz(il,3) = leaf_cz(il)
       tmp(il)    = 2.0_wp*leaf_half(il)
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
    if (status /= 0) call fatal_error('failed to write dust-emissivity file: '//trim(outname))
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
