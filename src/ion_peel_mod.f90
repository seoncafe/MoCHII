module ion_peel_mod
!---------------------------------------------------------------------------
! MoCHII: peel-off imaging of the ionizing/FUV band (adapted from
! MoCafe_v2.00/src/observer_mod.f90 + peelingoff_mod.f90, 2026-07-12).
!
! With par%ion_peel, every emitted packet (stellar AND diffuse) peels a
! direct contribution toward each observer, and — when dust scattering is
! sampled (par%ion_dust_scatter) — every interaction peels a scattered
! contribution with the Henyey-Greenstein phase function of the packet's
! bin.  The optical depth to the box edge along the observer direction is
! integrated through kap_ion at the packet's bin (extinction: gas + dust
! absorption + dust scattering), so the images carry the same physics as
! the transport.
!
! Image planes (TAN projection, MoCafe convention): direct and scattered
! surface-brightness images [erg/s/cm^2/sr at the observer], split into
! the ionizing (E >= par%eion_min) and FUV channels when par%add_fuv.
! Observers are defined exactly as in MoCafe: par%obsx/y/z (positions or
! directions) or par%alpha/beta/gamma angle lists, par%distance,
! par%nxim/nyim, par%dxim/dyim (auto-sized from par%rmax when unset).
! Output: '<base>_image' (HDF5/FITS) — direct_ion/scatt_ion (+ _fuv)
! blocks for each observer (suffix _obs<k> when par%nobs > 1).
!
! Images are plain per-rank arrays reduced onto rank 0 at output time
! (small next to the leaf arrays).
!---------------------------------------------------------------------------
  use define
  use utility, only : is_finite
  implicit none
  private

  public :: ion_peel_setup, ion_peel_direct, ion_peel_scatter, &
            ion_peel_write

  !--- (nxim, nyim, channel, nobs); channel 1 = ionizing, 2 = FUV
  real(kind=wp), allocatable :: img_dir(:,:,:,:), img_sca(:,:,:,:)
  integer :: nchan = 1

contains

  !=========================================================================
  ! Observer geometry (trimmed MoCafe observer_create: no scan/sed/stokes
  ! image blocks; the rotation-matrix and image-plane conventions are
  ! identical).
  !=========================================================================
  subroutine ion_peel_setup()
    use mpi
    implicit none
    real(kind=wp), allocatable :: cosa(:), cosb(:), cosg(:)
    real(kind=wp), allocatable :: sina(:), sinb(:), sing(:)
    real(kind=wp) :: dist_scale
    integer :: i, ierr

    if (any(is_finite(par%phase_angle)))       par%alpha = -par%phase_angle
    if (any(is_finite(par%inclination_angle))) par%beta  = -par%inclination_angle
    if (any(is_finite(par%position_angle)))    par%gamma = -par%position_angle
    where (is_finite(par%beta)  .and. .not. is_finite(par%alpha)) par%alpha = 0.0_wp
    where (is_finite(par%alpha) .and. .not. is_finite(par%beta))  par%beta  = 0.0_wp

    !--- default: one observer on the +z axis, far away.
    if (.not. (is_finite(par%alpha(1)) .and. is_finite(par%beta(1))) .and. &
        .not. (is_finite(par%obsx(1))  .and. is_finite(par%obsy(1)) .and. &
               is_finite(par%obsz(1)))) then
       if (.not. is_finite(par%distance)) &
          par%distance = maxval([par%xmax, par%ymax, par%zmax])*1.0e3_wp
       par%obsx(1) = 0.0_wp;  par%obsy(1) = 0.0_wp;  par%obsz(1) = par%distance
       par%alpha(1) = 0.0_wp; par%beta(1) = 0.0_wp
    end if

    if (is_finite(par%alpha(1)) .and. is_finite(par%beta(1))) then
       par%nobs = count(is_finite(par%alpha) .and. is_finite(par%beta))
       allocate(cosa(par%nobs), cosb(par%nobs), cosg(par%nobs), &
                sina(par%nobs), sinb(par%nobs), sing(par%nobs))
       if (.not. allocated(observer)) allocate(observer(par%nobs))
       if (.not. is_finite(par%distance)) &
          par%distance = maxval([par%xmax, par%ymax, par%zmax])*1.0e3_wp
       do i = 1, par%nobs
          if (.not. is_finite(par%gamma(i))) then
             if (par%beta(i) > 0.0_wp .and. par%beta(i) <= 90.0_wp) then
                par%gamma(i) = 90.0_wp
             else if (par%beta(i) > 90.0_wp) then
                par%gamma(i) = -90.0_wp
             else
                par%gamma(i) = 0.0_wp
             end if
          end if
          cosa(i) = cos(par%alpha(i)*deg2rad); sina(i) = sin(par%alpha(i)*deg2rad)
          cosb(i) = cos(par%beta(i)*deg2rad);  sinb(i) = sin(par%beta(i)*deg2rad)
          cosg(i) = cos(par%gamma(i)*deg2rad); sing(i) = sin(par%gamma(i)*deg2rad)
          observer(i)%x = par%distance*cosa(i)*sinb(i)
          observer(i)%y = par%distance*sina(i)*sinb(i)
          observer(i)%z = par%distance*cosb(i)
          observer(i)%alpha = par%alpha(i);  observer(i)%beta = par%beta(i)
          observer(i)%gamma = par%gamma(i);  observer(i)%distance = par%distance
       end do
    else
       par%nobs = count(is_finite(par%obsx) .and. is_finite(par%obsy) .and. &
                        is_finite(par%obsz))
       allocate(cosa(par%nobs), cosb(par%nobs), cosg(par%nobs), &
                sina(par%nobs), sinb(par%nobs), sing(par%nobs))
       if (.not. allocated(observer)) allocate(observer(par%nobs))
       if (.not. is_finite(par%distance)) then
          par%distance = sqrt(par%obsx(1)**2 + par%obsy(1)**2 + par%obsz(1)**2)
          if (par%distance < 10.0_wp*maxval([par%xmax, par%ymax, par%zmax])) &
             par%distance = maxval([par%xmax, par%ymax, par%zmax])*1.0e3_wp
       end if
       do i = 1, par%nobs
          if (.not. is_finite(par%gamma(i))) par%gamma(i) = 0.0_wp
          dist_scale = par%distance &
                       /sqrt(par%obsx(i)**2 + par%obsy(i)**2 + par%obsz(i)**2)
          if (dist_scale > 1.001_wp) then
             observer(i)%x = par%obsx(i)*dist_scale
             observer(i)%y = par%obsy(i)*dist_scale
             observer(i)%z = par%obsz(i)*dist_scale
          else
             observer(i)%x = par%obsx(i)
             observer(i)%y = par%obsy(i)
             observer(i)%z = par%obsz(i)
          end if
          cosb(i) = observer(i)%z/par%distance
          if (abs(cosb(i) - 1.0_wp) < eps) cosb(i) =  1.0_wp
          if (abs(cosb(i) + 1.0_wp) < eps) cosb(i) = -1.0_wp
          sinb(i) = sqrt(1.0_wp - cosb(i)**2)
          par%beta(i) = atan2(sinb(i), cosb(i))*rad2deg
          cosg(i) = cos(par%gamma(i)*deg2rad); sing(i) = sin(par%gamma(i)*deg2rad)
          if (sinb(i) == 0.0_wp) then
             cosa(i) = 1.0_wp;  sina(i) = 0.0_wp;  par%alpha(i) = 0.0_wp
          else
             par%alpha(i) = atan2(observer(i)%y, observer(i)%x)
             cosa(i) = cos(par%alpha(i));  sina(i) = sin(par%alpha(i))
             par%alpha(i) = par%alpha(i)*rad2deg
          end if
          observer(i)%alpha = par%alpha(i);  observer(i)%beta = par%beta(i)
          observer(i)%gamma = par%gamma(i);  observer(i)%distance = par%distance
       end do
    end if

    !--- rotation matrix, grid -> observer frame.
    do i = 1, par%nobs
       observer(i)%rmatrix(1,1) =  cosa(i)*cosb(i)*cosg(i) - sina(i)*sing(i)
       observer(i)%rmatrix(1,2) =  sina(i)*cosb(i)*cosg(i) + cosa(i)*sing(i)
       observer(i)%rmatrix(1,3) = -sinb(i)*cosg(i)
       observer(i)%rmatrix(2,1) = -cosa(i)*cosb(i)*sing(i) - sina(i)*cosg(i)
       observer(i)%rmatrix(2,2) = -sina(i)*cosb(i)*sing(i) + cosa(i)*cosg(i)
       observer(i)%rmatrix(2,3) =  sinb(i)*sing(i)
       observer(i)%rmatrix(3,1) =  cosa(i)*sinb(i)
       observer(i)%rmatrix(3,2) =  sina(i)*sinb(i)
       observer(i)%rmatrix(3,3) =  cosb(i)
    end do
    deallocate(cosa, cosb, cosg, sina, sinb, sing)

    !--- image plane: auto pixel scale from the sphere radius (par%rmax
    !--- is the half box for the AMR grid).
    if (.not. (is_finite(par%dxim) .and. is_finite(par%dyim))) then
       par%dxim = atan2(par%rmax*sqrt(3.0_wp), par%distance) &
                  /(par%nxim/2.0_wp)*rad2deg
       par%dyim = atan2(par%rmax*sqrt(3.0_wp), par%distance) &
                  /(par%nyim/2.0_wp)*rad2deg
    end if
    do i = 1, par%nobs
       observer(i)%nxim = par%nxim;  observer(i)%nyim = par%nyim
       observer(i)%dxim = par%dxim;  observer(i)%dyim = par%dyim
       observer(i)%steradian_pix = par%dxim*par%dyim*deg2rad**2
    end do

    !--- channel axis: band-integrated (1 or 2 channels) or, with
    !--- par%peel_bins, one image per band bin.
    if (par%peel_bins) then
       nchan = par%nnu_ion
    else
       nchan = merge(2, 1, par%add_fuv)
    end if
    allocate(img_dir(par%nxim, par%nyim, nchan, par%nobs), &
             img_sca(par%nxim, par%nyim, nchan, par%nobs))
    img_dir = 0.0_wp
    img_sca = 0.0_wp

    if (mpar%p_rank == 0) then
       write(*,'(a,i3,a,i5,a,i5,a)') ' PEEL: ', par%nobs, &
          ' observer(s), image ', par%nxim, ' x', par%nyim, ' pixels'
       write(*,'(a,2es12.4,a,es12.4)') ' PEEL: dxim, dyim [deg] = ', &
          par%dxim, par%dyim, ',  distance = ', par%distance
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine ion_peel_setup

  !=========================================================================
  ! channel of a band bin: 1 = ionizing, 2 = FUV (below par%eion_min).
  !=========================================================================
  integer function peel_channel(inu) result(ch)
    use ion_band_mod, only : ion_e
    integer, intent(in) :: inu
    if (par%peel_bins) then
       ch = inu
       return
    end if
    ch = 1
    if (nchan == 2 .and. ion_e(inu) < par%eion_min) ch = 2
  end function peel_channel

  !=========================================================================
  ! Direct peel at emission (stellar and diffuse packets): the isotropic
  ! source contributes Lpacket wgt/(4 pi r^2) e^-tau to its pixel.
  !=========================================================================
  subroutine ion_peel_direct(photon)
    use raytrace_amr_mod, only : raytrace_ion_tau_only_amr
    implicit none
    type(photon_type), intent(in) :: photon
    type(photon_type) :: pobs
    real(kind=wp) :: r2, r, vdet(3), tau, contrib
    integer :: ix, iy, i, k, ch

    ch = peel_channel(photon%inu)
    do k = 1, par%nobs
       pobs = photon
       pobs%kx = observer(k)%x - photon%x
       pobs%ky = observer(k)%y - photon%y
       pobs%kz = observer(k)%z - photon%z
       r2 = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
       r  = sqrt(r2)
       pobs%kx = pobs%kx/r;  pobs%ky = pobs%ky/r;  pobs%kz = pobs%kz/r
       do i = 1, 3
          vdet(i) = observer(k)%rmatrix(i,1)*pobs%kx &
                  + observer(k)%rmatrix(i,2)*pobs%ky &
                  + observer(k)%rmatrix(i,3)*pobs%kz
       end do
       ix = floor(atan2(-vdet(1), vdet(3))*rad2deg/observer(k)%dxim &
            + observer(k)%nxim/2.0_wp) + 1
       iy = floor(atan2(-vdet(2), vdet(3))*rad2deg/observer(k)%dyim &
            + observer(k)%nyim/2.0_wp) + 1
       if (ix < 1 .or. ix > observer(k)%nxim .or. &
           iy < 1 .or. iy > observer(k)%nyim) cycle
       call raytrace_ion_tau_only_amr(pobs, tau)
       if (tau > 500.0_wp) cycle      ! e^-tau = 0 to any precision
       contrib = photon%Lpacket*photon%wgt/(fourpi*r2)*exp(-tau)
       img_dir(ix, iy, ch, k) = img_dir(ix, iy, ch, k) + contrib
    end do
  end subroutine ion_peel_direct

  !=========================================================================
  ! Scattered peel at a dust interaction: called AFTER the albedo weight
  ! is applied and BEFORE the new direction is sampled (photon%k is the
  ! incident direction).  HG phase function with the bin's g.
  !=========================================================================
  subroutine ion_peel_scatter(photon)
    use raytrace_amr_mod, only : raytrace_ion_tau_only_amr
    use gas_opacity_mod,  only : ion_dust_g
    implicit none
    type(photon_type), intent(in) :: photon
    type(photon_type) :: pobs
    real(kind=wp) :: r2, r, vdet(3), tau, cosa, gg, peel, contrib
    integer :: ix, iy, i, k, ch

    ch = peel_channel(photon%inu)
    !--- the D03 table carries <cos> = 1.000 in the EUV; g = 1 with a
    !--- perfectly forward peel direction gives 0/0 in the HG phase
    !--- function — clamp (one NaN contribution poisons the image sum).
    gg = min(ion_dust_g(photon%inu), 0.9999_wp)
    do k = 1, par%nobs
       pobs = photon
       pobs%kx = observer(k)%x - photon%x
       pobs%ky = observer(k)%y - photon%y
       pobs%kz = observer(k)%z - photon%z
       r2 = pobs%kx**2 + pobs%ky**2 + pobs%kz**2
       r  = sqrt(r2)
       pobs%kx = pobs%kx/r;  pobs%ky = pobs%ky/r;  pobs%kz = pobs%kz/r
       do i = 1, 3
          vdet(i) = observer(k)%rmatrix(i,1)*pobs%kx &
                  + observer(k)%rmatrix(i,2)*pobs%ky &
                  + observer(k)%rmatrix(i,3)*pobs%kz
       end do
       ix = floor(atan2(-vdet(1), vdet(3))*rad2deg/observer(k)%dxim &
            + observer(k)%nxim/2.0_wp) + 1
       iy = floor(atan2(-vdet(2), vdet(3))*rad2deg/observer(k)%dyim &
            + observer(k)%nyim/2.0_wp) + 1
       if (ix < 1 .or. ix > observer(k)%nxim .or. &
           iy < 1 .or. iy > observer(k)%nyim) cycle
       call raytrace_ion_tau_only_amr(pobs, tau)
       !--- e^-tau = 0 to any precision; skipping also avoids a
       !--- -fp-model fast underflow pathology (the -O0 -check build is
       !--- clean, the fast build produced Inf/NaN pixels through the
       !--- extreme-dynamic-range exp/accumulate chain).
       if (tau > 500.0_wp) cycle
       cosa = photon%kx*pobs%kx + photon%ky*pobs%ky + photon%kz*pobs%kz
       peel = (1.0_wp - gg*gg) &
              /max((1.0_wp + gg*gg) - 2.0_wp*gg*cosa, tinest)**1.5_wp/fourpi
       contrib = photon%Lpacket*photon%wgt*peel/r2*exp(-tau)
       !--- defensive: a non-finite contribution would poison the whole
       !--- image sum.  Bit test on the exponent field — an ordinary
       !--- NaN comparison is folded away under -fp-model fast.
       if (nonfinite(contrib)) then
          write(*,'(a,i4,7es13.5,4l2)') ' PEEL: non-finite contrib:', &
             photon%inu, gg, cosa, tau, photon%wgt, photon%Lpacket, r2, &
             peel, nonfinite(tau), nonfinite(photon%wgt), &
             nonfinite(peel), nonfinite(cosa)
          cycle
       end if
       img_sca(ix, iy, ch, k) = img_sca(ix, iy, ch, k) + contrib
       !--- event trap: the pixel must stay finite if every added
       !--- contrib is finite; firing here with a finite contrib means
       !--- the += chain itself (or another writer) is at fault.
       if (nonfinite(img_sca(ix, iy, ch, k))) &
          write(*,'(a,2i5,2es13.5)') ' PEEL: pixel went non-finite at', &
             ix, iy, contrib, img_sca(ix, iy, ch, k)
    end do
  end subroutine ion_peel_scatter

  !=========================================================================
  ! .true. for NaN or Inf: exponent bits all set.  Immune to -fp-model
  ! fast folding (pure integer arithmetic).
  !=========================================================================
  logical function nonfinite(x)
    use iso_fortran_env, only : int64
    real(kind=wp), intent(in) :: x
    integer(int64) :: bits
    bits = transfer(x, 0_int64)
    nonfinite = iand(ishft(bits, -52), 2047_int64) == 2047_int64
  end function nonfinite

  !=========================================================================
  ! Reduce and write the images: surface brightness at the observer,
  ! [erg/s/cm^2/sr] (contributions are erg/s per code-length^2; divide by
  ! distance2cm^2 and by the pixel solid angle).
  !=========================================================================
  subroutine ion_peel_write()
    use mpi
    use iofile_mod
    use utility, only : get_base_name
    implicit none
    type(io_file_type) :: file
    character(len=192) :: outname
    character(len=64)  :: extname, suffix
    real(kind=wp) :: fac
    integer :: k, ch, ierr, status
    character(len=8), parameter :: chname(2) = ['ion     ', 'fuv     ']

    !--- NaN bisection: count non-finite pixels rank-locally BEFORE the
    !--- reduce and globally after (the guarded contrib line never
    !--- fires, so the poison enters elsewhere).
    block
      integer :: i1, i2, i3, i4, nbad
      nbad = 0
      do i4 = 1, size(img_sca,4);  do i3 = 1, size(img_sca,3)
      do i2 = 1, size(img_sca,2);  do i1 = 1, size(img_sca,1)
         if (nonfinite(img_sca(i1,i2,i3,i4))) then
            nbad = nbad + 1
            if (nbad <= 3) write(*,'(a,i4,a,4i5,es13.5)') &
               ' PEEL: rank', mpar%p_rank, ' local non-finite at', &
               i1, i2, i3, i4, img_sca(i1,i2,i3,i4)
         end if
      end do;  end do;  end do;  end do
      if (nbad > 0) write(*,'(a,i4,a,i6)') ' PEEL: rank', mpar%p_rank, &
         ' local non-finite count =', nbad
    end block
    call MPI_ALLREDUCE(MPI_IN_PLACE, img_dir, size(img_dir), &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, img_sca, size(img_sca), &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (mpar%p_rank /= 0) return

    status = 0
    outname = trim(get_base_name(par%out_file))//'_image'// &
              trim(io_file_extension(par%file_format))
    call io_open_new(file, trim(outname), status)
    do k = 1, par%nobs
       fac = 1.0_wp/(par%distance2cm**2 * observer(k)%steradian_pix)
       suffix = ''
       if (par%nobs > 1) write(suffix,'(a,i0)') '_obs', k

       if (par%peel_bins) then
          !--- bin-resolved cubes + the band-integrated channel images
          !--- derived from them.
          block
            use ion_band_mod, only : ion_e
            real(kind=wp), allocatable :: img2(:,:)
            integer :: b, ich, nch2
            call io_append_image(file, img_dir(:,:,:,k)*fac, status, bitpix=-64)
            write(extname,'(2a)') 'direct_cube', trim(suffix)
            call io_put_keyword(file,'EXTNAME',trim(extname), &
                 'direct I(x,y,bin) [erg/s/cm^2/sr]',status)
            call io_put_keyword(file,'DXIM', par%dxim, 'pixel [deg]', status)
            call io_put_keyword(file,'DIST', observer(k)%distance, &
                 'observer distance (code units)', status)
            call io_append_image(file, img_sca(:,:,:,k)*fac, status, bitpix=-64)
            write(extname,'(2a)') 'scatt_cube', trim(suffix)
            call io_put_keyword(file,'EXTNAME',trim(extname), &
                 'scattered I(x,y,bin) [erg/s/cm^2/sr]',status)
            allocate(img2(par%nxim, par%nyim))
            nch2 = merge(2, 1, par%add_fuv)
            do ich = 1, nch2
               img2 = 0.0_wp
               do b = 1, par%nnu_ion
                  if (ich == 1 .and. ion_e(b) <  par%eion_min) cycle
                  if (ich == 2 .and. ion_e(b) >= par%eion_min) cycle
                  img2 = img2 + img_dir(:,:,b,k)
               end do
               call io_append_image(file, img2*fac, status, bitpix=-64)
               write(extname,'(4a)') 'direct_', trim(chname(ich)), trim(suffix)
               call io_put_keyword(file,'EXTNAME',trim(extname), &
                    'direct surface brightness [erg/s/cm^2/sr]',status)
               img2 = 0.0_wp
               do b = 1, par%nnu_ion
                  if (ich == 1 .and. ion_e(b) <  par%eion_min) cycle
                  if (ich == 2 .and. ion_e(b) >= par%eion_min) cycle
                  img2 = img2 + img_sca(:,:,b,k)
               end do
               call io_append_image(file, img2*fac, status, bitpix=-64)
               write(extname,'(4a)') 'scatt_', trim(chname(ich)), trim(suffix)
               call io_put_keyword(file,'EXTNAME',trim(extname), &
                    'scattered surface brightness [erg/s/cm^2/sr]',status)
            end do
            deallocate(img2)
          end block
       else
          do ch = 1, nchan
             call io_append_image(file, img_dir(:,:,ch,k)*fac, status, bitpix=-64)
             write(extname,'(4a)') 'direct_', trim(chname(ch)), trim(suffix)
             call io_put_keyword(file,'EXTNAME',trim(extname), &
                  'direct surface brightness [erg/s/cm^2/sr]',status)
             call io_put_keyword(file,'DXIM',   par%dxim,  'pixel [deg]', status)
             call io_put_keyword(file,'DIST',   observer(k)%distance, &
                  'observer distance (code units)', status)
             call io_put_keyword(file,'ALPHA',  observer(k)%alpha, '[deg]', status)
             call io_put_keyword(file,'BETA',   observer(k)%beta,  '[deg]', status)
             call io_put_keyword(file,'GAMMA',  observer(k)%gamma, '[deg]', status)
             call io_append_image(file, img_sca(:,:,ch,k)*fac, status, bitpix=-64)
             write(extname,'(4a)') 'scatt_', trim(chname(ch)), trim(suffix)
             call io_put_keyword(file,'EXTNAME',trim(extname), &
                  'scattered surface brightness [erg/s/cm^2/sr]',status)
          end do
       end if
    end do
    call io_close(file, status)
    write(*,'(2a)') ' PEEL: images written to: ', trim(outname)
    write(*,'(a,es12.4,a,es12.4)') ' PEEL: sum(direct) = ', sum(img_dir), &
       ',  sum(scatt) = ', sum(img_sca)
  end subroutine ion_peel_write

end module ion_peel_mod
