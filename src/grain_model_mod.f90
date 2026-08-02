module grain_model_mod
!---------------------------------------------------------------------------
! MoCHII: the run's grain population model (SEDust dust_model_t), shared as
! the SINGLE source of the dust optics.  One object supplies all three parts
! of the dust physics from the same size distribution:
!   extinction / absorption  gas_opacity_mod   (ionizing-band transport)
!   equilibrium T_dust        dust_temp_mod    (grain heating <-> IR)
!   stochastic IR emission    sedust_mod       (reemitted spectrum)
! Taking C_ext, albedo and <cos> from here rather than from a separate kext
! file (par%ion_dust_kext) is what makes the absorbed power and the reemitted
! spectrum refer to the SAME grains: the energy a cell absorbs is set by
! C_abs of the size distribution that then radiates it away.  Naming the
! model (par%dust_model) is therefore enough to fix the dust physics.
!
! EUV extension: the astrodust and DL07 models take their wavelength grid from
! the T-matrix Q table, which bottoms out at 0.0912 um (13.6 eV), so a
! photoionization run that transports shortward of that (par%ion_add_dust, band
! up to par%eion_max) asks SEDust for the grid extended down to
! 1.23984/eion_max um.  Both can be extended: DL07 is dielectric-function Mie
! throughout, and astrodust computes the extended band from the DH21 dielectric
! function.  The Zubko model needs no extension -- its own DustEM Q tables span
! 0.001-10000 um (1.24e-4 eV to 1.24 keV), so build_zubko already returns a grid
! that reaches well into the EUV, which is why it takes no lam_min.
!
! Whether the grid actually spans the transported band is checked, with the
! numbers, in gas_opacity_mod at setup -- for every model alike.  A model that
! falls short stops the run there and names par%ion_dust_kext as the way out.
!
! Where the transport optics come from: SEDust serves them from the model's
! PRECOMPUTED size-integrated extinction table (data/kext_*.dat), loaded by the
! builder and interpolated onto the model grid by dust_extinction.  The size
! integral itself is no longer done at call time; see grain_extinction_table
! below for what the two differ by and for what the table presumes about the
! par%sed_* inputs.  Which table is a run's choice: par%sed_kext blank (the
! default) takes SEDust's own standard table for the named model, so those file
! names live in SEDust alone; naming a file makes it mandatory.
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
  !--- Build the model named by par%dust_model.  Every rank builds its
  !--- own copy from the same files, so the result is identical without a
  !--- broadcast.  When par%ion_add_dust, the astrodust/DL07 grid is extended
  !--- down to the shortest transported wavelength (1.23984/eion_max um) so
  !--- the ionizing band is resolved; without it (emission-only run) the grid
  !--- is the native table, bit-identical to the pre-shared inline build.
  subroutine build_grain_model()
    use mpi
    use utility, only : fatal_error
    use ifport,  only : chdir, getcwd
    implicit none
    integer :: ierr, cstat, st, st_kext
    character(len=512) :: cwd_save, msg
    real(kind=wp)      :: lam_min_band
    logical            :: want_euv, kext_named

    if (grain_model_ready) return
    st = 0
    !--- the build status that means "the named extinction table could not be
    !--- read": 5 for build_astrodust / build_dl07, 6 for build_zubko.
    st_kext = 0

    !--- A blank par%sed_kext leaves the extinction table to SEDust's own default
    !--- for the model, so the standard table names stay in one place (SEDust).
    !--- A named file is instead mandatory: the build fails if it cannot be read.
    kext_named = len_trim(par%sed_kext) > 0

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

    select case (trim(par%dust_model))
    case ('astrodust')
       st_kext = 5
       if (want_euv) then
          if (kext_named) then
             call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                  par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, lam_min=lam_min_band, &
                  kext_path=trim(par%sed_kext))
          else
             call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                  par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, lam_min=lam_min_band)
          end if
       else
          if (kext_named) then
             call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                  par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, &
                  kext_path=trim(par%sed_kext))
          else
             call build_astrodust(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                  par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
          end if
       end if
    case ('dl07')
       st_kext = 5
       if (want_euv) then
          if (kext_named) then
             call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                  par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
                  par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, lam_min=lam_min_band, &
                  kext_path=trim(par%sed_kext))
          else
             call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                  par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
                  par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, lam_min=lam_min_band)
          end if
       else
          if (kext_named) then
             call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                  par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
                  par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, &
                  kext_path=trim(par%sed_kext))
          else
             call build_dl07(dmodel, trim(par%sed_qtable), trim(par%sed_sizedist), &
                  par%sed_dl07_sdindex, par%sed_dl07_uisrf, &
                  par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
          end if
       end if
    case ('zubko')
       !--- No lam_min here, and none is needed: the Zubko DustEM Q tables span
       !--- 0.001-10000 um (1.24e-4 eV to 1.24 keV), and build_zubko takes the
       !--- model grid straight from them, so the ionizing band is already on
       !--- the grid.  gas_opacity_mod checks that for real at setup.
       st_kext = 6
       if (kext_named) then
          call build_zubko(dmodel, trim(par%sed_zubko_config), trim(par%sed_zubko_dir), &
               par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st, &
               kext_path=trim(par%sed_kext))
       else
          call build_zubko(dmodel, trim(par%sed_zubko_config), trim(par%sed_zubko_dir), &
               par%sed_NT, par%sed_Tlo, par%sed_Thi, status=st)
       end if
    case default
       cstat = chdir(trim(cwd_save))
       if (mpar%p_rank == 0) write(*,'(3a)') &
          'ERROR: par%dust_model = ''', trim(par%dust_model), &
          ''' unknown (use ''astrodust'', ''dl07'', or ''zubko'').'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end select

    cstat = chdir(trim(cwd_save))

    !--- SEDust reports a missing or malformed input through status instead of
    !--- stopping itself, so a bad par%sed_* path fails cleanly on every rank
    !--- here rather than aborting mid-build.
    !--- Every rank builds the same model from the same files, so every rank
    !--- reaches the same status and every rank aborts.  A rank other than 0
    !--- therefore tears the job down while rank 0 is still writing, and a
    !--- rank-0-only message is lost -- the user sees nothing but MPI_Abort.
    !--- utility::fatal_error exists for exactly this: it writes from whichever
    !--- rank reached the failure and flushes before aborting.  The rank-0 lines
    !--- below are the extra detail, and the one-line cause goes through it.
    if (st /= 0) then
       if (kext_named .and. st == st_kext) then
          !--- the one status a named par%sed_kext can produce: the transport
          !--- extinction curve this run asked for is not readable.  A blank
          !--- par%sed_kext cannot land here -- SEDust only OFFERS its default
          !--- table, and a missing one surfaces later as dust_extinction
          !--- status 2 instead.
          write(msg,'(5a)') 'par%sed_kext could not be read: ''', &
             trim(par%sed_kext), ''' (relative to par%sed_workdir = ''', &
             trim(par%sed_workdir), ''')'
       else
          write(msg,'(3a,i0,a)') 'SEDust could not build the ''', &
             trim(par%dust_model), ''' dust model (status=', st, &
             '); check the par%sed_* input paths'
       end if
       call fatal_error(trim(msg))
    end if

    grain_model_ready = .true.
  end subroutine build_grain_model

  !=========================================================================
  !--- The model's size-integrated extinction curve, on the model's own
  !--- wavelength grid, in the (lambda, albedo, <cos>, C_ext) form the transport
  !--- reads.  The curve is not integrated here: dust_extinction returns the
  !--- PRECOMPUTED table the builder loaded (par%sed_kext, or SEDust's standard
  !--- table for the model), interpolated onto
  !--- the model grid -- log(lambda)-log(C) for the cross sections, linear in
  !--- log(lambda) for <cos>, and the table value itself wherever a model
  !--- wavelength coincides with a table wavelength.
  !---
  !--- The size integral
  !---
  !---   C_ext(lambda) = sum_pop sum_a dn(a) [C_abs(lambda,a) + C_sca(lambda,a)]
  !---   albedo        = C_sca / C_ext
  !---   <cos>         = sum dn C_sca g / sum dn C_sca
  !---
  !--- is the calculation that WROTE that table, not the calculation this
  !--- routine performs.  It is still reachable as SEDust's
  !--- size_integrated_extinction for a host that wants to redo it.
  !---
  !--- How far the two are apart, as tests/grain_kext measures it -- worst
  !--- single wavelength, relative on C_ext and C_abs, absolute on the pair
  !--- (albedo, albedo*<cos>) the transport uses as probabilities:
  !---   * lambda >= 0.0912 um the model grid IS the table grid, so the table
  !---     value comes back unchanged: 4.9e-13 on the cross sections, 4.2e-13 on
  !---     the albedo pair.
  !---   * below 0.0912 um, on the EUV extension MoCHII asks for, the extended
  !---     nodes are not the table's nodes, so the log-log interpolation enters:
  !---     1.16e-2 for astrodust on the eion_max=100 eV grid (8.20e-3 at 150 eV),
  !---     1.13e-3 for DL07, and 4.82e-13 for Zubko, whose grid is the table's
  !---     own; 3.3e-3 / 9.7e-3 / 2.2e-4 / 2.8e-13 on the albedo pair.
  !---   * those departures alternate in sign and cancel in the integrals that
  !---     matter: the band-averaged s_abs over 13.598-100 eV agrees to 4.1e-6
  !---     (astrodust 100 eV), 2.0e-5 (150 eV), 3.9e-6 (DL07), and C_ext at
  !---     lambda_ref = 0.55 um to 2.8e-14.
  !---
  !--- What the table presumes: it was computed from the model's DEFAULT inputs.
  !--- Changing par%sed_sizedist / par%sed_qtable / par%sed_dl07_sdindex /
  !--- par%sed_zubko_config moves the grains that emit without moving the curve
  !--- that extinguishes, and the two then no longer describe the same dust.
  !--- tests/grain_kext catches exactly that: it compares the served curve with
  !--- the size integral over the model as built, so a table belonging to other
  !--- grains shows up as a large departure at lambda >= 0.0912 um.
  !---
  !--- Cross sections are per H atom [cm^2/H], the same normalization as the
  !--- kext table read from a file, so gas_opacity_mod / dust_temp_mod divide
  !--- by C_ext(lambda_ref) into the same dimensionless shapes either way.
  !--- Populations without scattering optics (the PAHs) enter through
  !--- absorption only, so the albedo and <cos> fall where they dominate.
  subroutine grain_extinction_table(lam, alb, cosg, cext, n)
    use mpi
    use utility, only : fatal_error
    implicit none
    real(kind=wp), allocatable, intent(out) :: lam(:), alb(:), cosg(:), cext(:)
    integer,                    intent(out) :: n
    real(kind=wp), allocatable :: cabs(:), csca(:)
    character(len=512) :: msg
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
       select case (st)
       case (2)
          !--- SEDust's standard table for this model could not be read.  It is
          !--- only an OFFER, so a missing one does not fail the build -- an
          !--- emission-only driver still runs -- and it surfaces here, the first
          !--- time transport asks for optics.  A named par%sed_kext fails
          !--- earlier, in build_grain_model.
          write(msg,'(5a)') 'the ''', trim(par%dust_model), &
             ''' model carries no extinction table: SEDust''s standard table '// &
             'for it is missing under par%sed_workdir = ''', &
             trim(par%sed_workdir), &
             ''' (run SEDust/populate_data.sh, or name one with par%sed_kext)'
       case (3)
          write(msg,'(a,es10.3,a,es10.3,3a)') &
             'the model grid [um] = [', lam(1), ', ', lam(n), &
             '] reaches outside the extinction table ''', &
             trim(dmodel%kext_path), &
             ''' -- no extrapolation: name a wider par%sed_kext, or raise '// &
             '1.23984/par%eion_max'
       case default
          write(msg,'(a,i0,a)') 'dust_extinction failed (status=', st, ')'
       end select
       call fatal_error(trim(msg))
    end if

    !--- Every model SEDust builds carries scattering optics, so a C_sca that
    !--- vanishes at every wavelength means the optics failed to load rather
    !--- than that the dust does not scatter.  An albedo forced to zero would
    !--- turn grains whose true albedo is ~0.5-0.7 in the optical into a purely
    !--- absorbing medium, so the run stops rather than transport that silently.
    if (maxval(csca) <= 0.0_wp) then
       if (mpar%p_rank == 0) then
          write(*,'(3a)') 'ERROR: the ''', trim(par%dust_model), &
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
