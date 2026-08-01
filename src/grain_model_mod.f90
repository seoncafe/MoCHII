module grain_model_mod
!---------------------------------------------------------------------------
! MoCHII: the run's grain population model (SEDust dust_model_t), shared as
! the SINGLE source of the dust optics.  One object supplies all three
! halves of the dust physics from the same size distribution:
!   extinction / absorption  gas_opacity_mod   (ionizing-band transport)
!   equilibrium T_dust        dust_temp_mod    (grain heating <-> IR)
!   stochastic IR emission    sedust_mod       (reemitted spectrum)
! Taking C_ext, albedo and <cos> from here rather than from a separate kext
! file (par%ion_dust_kext) is what makes the absorbed power and the reemitted
! spectrum refer to the SAME grains: the energy a cell absorbs is set by
! C_abs of the size distribution that then radiates it away.  Naming the
! model (par%dust_model_sed) is therefore enough to fix the dust physics.
!
! EUV extension: the T-matrix Q table bottoms out at 0.0912 um (13.6 eV), so
! a photoionization run that transports shortward of that (par%ion_add_dust,
! band up to par%eion_max) needs the model grid extended down to
! 1.23984/eion_max um.  Only the astrodust and DL07 models can be extended
! (their optics are dielectric-function Mie throughout); the Zubko model's
! DustEM Q tables ARE its definition and cannot be extended, so a Zubko run
! that carries dust in the ionizing band must supply par%ion_dust_kext.
!
! Built once for each run; both callers go through build_grain_model, which
! returns immediately after the first call.
!---------------------------------------------------------------------------
  use define
  use dust_lib, only : dust_model_t, build_astrodust, build_dl07, build_zubko, &
                       dust_extinction, dust_nlam, dust_lambda
  implicit none
  private

  type(dust_model_t), save, public :: dmodel
  logical,            save, public :: grain_model_ready = .false.

  public :: build_grain_model, grain_extinction_table

contains
  !=========================================================================
  !--- Build the model named by par%dust_model_sed.  Every rank builds its
  !--- own copy from the same files, so the result is identical without a
  !--- broadcast.  When par%ion_add_dust, the astrodust/DL07 grid is extended
  !--- down to the shortest transported wavelength (1.23984/eion_max um) so
  !--- the ionizing band is resolved; without it (emission-only run) the grid
  !--- is the native table, bit-identical to the pre-shared inline build.
  subroutine build_grain_model()
    use mpi
    use ifport, only : chdir, getcwd
    implicit none
    integer :: ierr, cstat, st
    character(len=512) :: cwd_save
    real(kind=wp)      :: lam_min_band
    logical            :: want_euv

    if (grain_model_ready) return
    st = 0

    !--- transport / T_dust in the ionizing band need the EUV extension.
    want_euv     = par%ion_add_dust
    lam_min_band = 1.23984_wp/par%eion_max      ! shortest transported wavelength [um]

    !--- SEDust reads its dielectric tables via paths hard-coded relative to
    !--- its sed/ directory ('../data/dielectric/...'), so build the model from
    !--- par%sed_workdir (a SEDust sed/ tree) and restore the working directory.
    !--- par%sed_qtable / par%sed_sizedist are given relative to that directory
    !--- (or absolute).  chdir/getcwd are the Intel IFPORT integer functions
    !--- (return 0 on success).
    cstat = getcwd(cwd_save)
    if (len_trim(par%sed_workdir) > 0) then
       cstat = chdir(trim(par%sed_workdir))
       if (cstat /= 0 .and. mpar%p_rank == 0) write(*,'(3a)') &
          'WARNING: could not chdir to par%sed_workdir = ''', trim(par%sed_workdir), &
          ''' (SEDust dielectric files may not be found).'
    end if

    select case (trim(par%dust_model_sed))
    case ('astrodust')
       if (want_euv) then
          call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
               par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, lam_min=lam_min_band)
       else
          call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
               par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
       end if
    case ('dl07')
       if (want_euv) then
          call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
               par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
               par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, lam_min=lam_min_band)
       else
          call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
               par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
               par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
       end if
    case ('zubko')
       !--- The Zubko Q tables are the model definition and cannot be extended
       !--- into the EUV, so a run that transports dust in the ionizing band
       !--- must supply the extinction through par%ion_dust_kext instead.
       if (want_euv .and. mpar%p_rank == 0) write(*,'(a)') &
          ' WARNING: the zubko model has no EUV extension; '// &
          'set par%ion_dust_kext to cover the ionizing band.'
       call build_zubko(dmodel, trim(par%sed_zubko_config), trim(par%sed_zubko_dir), &
            par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
    case default
       cstat = chdir(trim(cwd_save))
       if (mpar%p_rank == 0) write(*,'(3a)') &
          'ERROR: par%dust_model_sed = ''', trim(par%dust_model_sed), &
          ''' unknown (use ''astrodust'', ''dl07'', or ''zubko'').'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end select

    cstat = chdir(trim(cwd_save))

    !--- SEDust reports a missing or malformed input through status instead of
    !--- stopping itself, so a bad par%sed_* path fails cleanly on every rank
    !--- here rather than aborting mid-build.
    if (st /= 0) then
       if (mpar%p_rank == 0) write(*,'(3a,i0,a)') &
          'ERROR: SEDust failed to build the ''', trim(par%dust_model_sed), &
          ''' dust model (status=', st, '); check the par%sed_* input paths.'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if

    grain_model_ready = .true.
  end subroutine build_grain_model

  !=========================================================================
  !--- Size-integrated extinction of the model, on the model's own wavelength
  !--- grid, in the (lambda, albedo, <cos>, C_ext) form the transport reads:
  !---
  !---   C_ext(lambda) = sum_pop sum_a dn(a) [C_abs(lambda,a) + C_sca(lambda,a)]
  !---   albedo        = C_sca / C_ext
  !---   <cos>         = sum dn C_sca g / sum dn C_sca
  !---
  !--- Cross sections are per H atom [cm^2/H], the same normalization as the
  !--- kext table read from a file, so gas_opacity_mod / dust_temp_mod divide
  !--- by C_ext(lambda_ref) into the same dimensionless shapes either way.
  !--- Populations without scattering optics (the PAHs) enter through
  !--- absorption only, so the albedo and <cos> fall where they dominate.
  subroutine grain_extinction_table(lam, alb, cosg, cext, n)
    use mpi
    implicit none
    real(kind=wp), allocatable, intent(out) :: lam(:), alb(:), cosg(:), cext(:)
    integer,                    intent(out) :: n
    real(kind=wp), allocatable :: cabs(:), csca(:)
    integer :: ierr, st, il

    call build_grain_model()

    n = dust_nlam(dmodel)
    if (n < 2) then
       if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
          'ERROR: the dust model has ', n, ' wavelength points; at least 2 are needed.'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if

    allocate(lam(n), alb(n), cosg(n), cext(n), cabs(n), csca(n))
    lam = dust_lambda(dmodel)
    call dust_extinction(dmodel, cext, cabs, csca, gbar=cosg, status=st)
    if (st /= 0) then
       if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
          'ERROR: dust_extinction failed (status=', st, ').'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if

    !--- Every model SEDust builds carries scattering optics, so a C_sca that
    !--- vanishes at every wavelength means the optics failed to load rather
    !--- than that the dust does not scatter.  An albedo forced to zero would
    !--- turn grains whose true albedo is ~0.5-0.7 in the optical into a purely
    !--- absorbing medium, so the run stops rather than transport that silently.
    if (maxval(csca) <= 0.0_wp) then
       if (mpar%p_rank == 0) then
          write(*,'(3a)') 'ERROR: the ''', trim(par%dust_model_sed), &
             ''' model returned no scattering at any wavelength;'
          write(*,'(a)')  '       check the par%sed_* optics paths.'
       end if
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if

    do il = 1, n
       if (cext(il) > 0.0_wp) then
          alb(il) = csca(il)/cext(il)
       else
          alb(il) = 0.0_wp
       end if
    end do
    deallocate(cabs, csca)
  end subroutine grain_extinction_table

end module grain_model_mod
