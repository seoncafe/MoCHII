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
! Where the optics come from: ONE SEDust entry point, build_dust, and ONE
! directory.  par%sed_data_dir/<par%dust_model>/sedust_<model>.h5 carries the
! model's wavelength axis, its cross-section tables and its size-integrated
! extinction curve (/kext) together, so the wavelengths the transport runs on
! and the curve it reads its cross sections off come from the same file and
! cannot be paired wrongly.  What is not optics sits beside that directory: the
! size distribution (release/size_distribution.dat), the dielectric functions
! and PAH cross sections the optics were computed on (dielectric/), and, for
! Zubko, the ZDA component definition.  Naming the model is therefore all a run
! does; only par%sed_kext overrides part of it, with an extinction curve of its
! own.  Nothing outside par%sed_data_dir is read, and nothing is resolved
! against the working directory, so the run needs no particular one.
!
! include_euv = .true., always, and not tied to par%ion_add_dust.  It selects
! the whole of the model's stored wavelength axis instead of its non-ionizing
! part, and MoCHII is a photoionization host: the grid is then the same object
! in every run, so a dust-emission-only run heats the same grains, on the same
! wavelengths, as a run that transports the ionizing band.
!
! No lam_min is passed, and none is needed.  (a) With include_euv the astrodust
! axis already reaches 1.0e-4 um (12398.4 eV), so lam_min = 1.23984/par%eion_max
! would do nothing at any par%eion_max below that.  (b) lam_min is a COVERAGE
! requirement rather than a grid: a model that cannot reach the asked-for
! wavelength refuses to build (Zubko's grid IS its optics table, so it refuses
! at 1.0e-3 um), and refusing is not what a host wants from a floor its own
! band never reaches.  (c) Whether the grid really spans the transported band
! is not something this module has to presume: gas_opacity_mod checks it, with
! the numbers, at setup, for every model alike, and a model that falls short
! stops the run there and names par%ion_dust_kext as the way out.
!
! The extinction the transport reads is the PRECOMPUTED size-integrated curve,
! loaded by the builder and interpolated onto the model grid by
! dust_extinction; the size integral itself is not done at call time.  See
! grain_extinction_table below for what the two differ by and for what the
! curve presumes about the par%sed_* inputs.
!
! Built once for each run; both callers go through build_grain_model, which
! returns immediately after the first call.
!---------------------------------------------------------------------------
  use define
  use dust_lib, only : dust_model_t, build_dust, &
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
  !--- broadcast.  One call for every model: build_dust resolves the optics
  !--- product, the size distribution and the ZDA config from par%sed_data_dir
  !--- and takes the whole stored wavelength axis (include_euv), so the only
  !--- choice left here is whether par%sed_kext replaces the extinction curve.
  subroutine build_grain_model()
    use utility, only : fatal_error
    implicit none
    integer :: st
    character(len=512) :: msg
    character(len=256) :: sed_msg
    logical            :: kext_named

    if (grain_model_ready) return
    st = 0

    !--- A blank par%sed_kext leaves the curve to /kext of the model's own
    !--- product, so the standard file names stay in one place (SEDust); a named
    !--- file replaces it.  Either way the curve is required and the build fails
    !--- without it, so both cases reach status 5 below.
    kext_named = len_trim(par%sed_kext) > 0

    !--- sd_index / u_isrf are DL07's alone (the WD01 size distribution and the
    !--- field its PAH ionization balance is computed at); the other models
    !--- ignore them, so they are passed unconditionally.
    if (kext_named) then
       call build_dust(dmodel, trim(par%dust_model), trim(par%sed_data_dir), &
            par%sed_NT, par%sed_Tlo, par%sed_Thi, include_euv=.true., status=st, &
            message=sed_msg, &
            sd_index=par%sed_dl07_sdindex, u_isrf=par%sed_dl07_uisrf, &
            kext_path=trim(par%sed_kext))
    else
       call build_dust(dmodel, trim(par%dust_model), trim(par%sed_data_dir), &
            par%sed_NT, par%sed_Tlo, par%sed_Thi, include_euv=.true., status=st, &
            message=sed_msg, &
            sd_index=par%sed_dl07_sdindex, u_isrf=par%sed_dl07_uisrf)
    end if

    !--- What st says.  build_dust points the whole library at par%sed_data_dir
    !--- for the length of the build, so EVERY file the model is made of is
    !--- resolved from there and every failure to read one comes back as a
    !--- status in one vocabulary, the same numbers whichever model is named:
    !--- 1 optics table, 2 size distribution, 3 dielectric function, 5 the
    !--- extinction curve, 6 calorimetry, 9 the model definition, and 90 / 91
    !--- for a model name MoCHII cannot offer.  A directory that does not exist
    !--- is a status too, not a runtime abort, so there is nothing left for this
    !--- module to guard with a working directory of its own.
    !--- Every rank builds the same model from the same files, so every rank
    !--- reaches the same status and every rank aborts.  A rank other than 0
    !--- therefore tears the job down while rank 0 is still writing, and a
    !--- rank-0-only message is lost -- the user sees nothing but MPI_Abort.
    !--- utility::fatal_error exists for exactly this: it writes from whichever
    !--- rank reached the failure and flushes before aborting.
    if (st /= 0) then
       if (st == 90) then
          !--- build_dust knows four models and MoCHII offers three of them.
          write(msg,'(3a)') 'par%dust_model = ''', trim(par%dust_model), &
             ''' unknown (use ''astrodust'', ''dl07'', or ''zubko'')'
       else if (st == 91) then
          !--- SEDust's fourth model is defined by a descriptor file, and MoCHII
          !--- has no input that names one.
          write(msg,'(a)') 'par%dust_model = ''from_files'' needs a model '// &
             'descriptor, which MoCHII has no input for (use ''astrodust'', '// &
             '''dl07'', or ''zubko'')'
       else if (st == 5) then
          !--- the extinction curve the transport reads its cross sections off
          !--- is not readable.
          if (kext_named) then
             write(msg,'(3a)') 'par%sed_kext could not be read: ''', &
                trim(par%sed_kext), ''''
          else
             write(msg,'(5a)') 'the ''', trim(par%dust_model), &
                ''' model carries no extinction curve: /kext of its product '// &
                'under par%sed_data_dir = ''', trim(par%sed_data_dir), &
                ''' could not be read (run SEDust/populate_data.sh, or name a '// &
                'curve with par%sed_kext)'
          end if
       else
          !--- Everything else is a file under par%sed_data_dir, and SEDust
          !--- names which one in words.
          write(msg,'(5a,i0,3a)') 'SEDust could not build the ''', &
             trim(par%dust_model), ''' dust model from par%sed_data_dir = ''', &
             trim(par%sed_data_dir), ''' (status=', st, '): ', trim(sed_msg), &
             ' -- run SEDust/populate_data.sh if the tree is incomplete'
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
  !--- (albedo, albedo*<cos>) the transport uses as probabilities.  The curve
  !--- and the model grid come out of one product, so every wavelength MoCHII
  !--- builds is a node of the curve; nothing is interpolated, and what is left
  !--- is rounding on a double:
  !---   * lambda >= 0.0912 um: 0.0 / 0.0 (astrodust), 2.7e-15 / 6.7e-16
  !---     (DL07), 5.8e-16 / 2.2e-16 (Zubko).
  !---   * lambda < 0.0912 um, the ionizing band itself: 3.3e-16 / 1.1e-16
  !---     (astrodust), 4.2e-16 / 3.3e-16 (DL07), 6.3e-16 / 2.2e-16 (Zubko).
  !---   * the integrals: the band-averaged s_abs over 13.598-100 eV agrees to
  !---     0.0 (astrodust), 0.0 (DL07), 6.9e-15 (Zubko), and C_ext at
  !---     lambda_ref = 0.55 um to 0.0 / 0.0 / 7.2e-15.
  !---
  !--- What the table presumes: it was computed from the model's DEFAULT inputs.
  !--- Changing par%sed_dl07_sdindex / par%sed_dl07_uisrf, or pointing
  !--- par%sed_kext at a curve computed for other grains, moves the grains that
  !--- emit without moving the curve that extinguishes (or the reverse), and the
  !--- two then no longer describe the same dust.
  !--- tests/grain_kext catches exactly that: it compares the served curve with
  !--- the size integral over the model as built, so a table belonging to other
  !--- grains shows up as a large departure, at any wavelength.
  !---
  !--- Cross sections are per H atom [cm^2/H], the same normalization as the
  !--- kext table read from a file, so gas_opacity_mod / dust_temp_mod divide
  !--- by C_ext(lambda_ref) into the same dimensionless shapes either way.
  !--- A population stored without scattering optics enters through absorption
  !--- only, so the albedo and <cos> fall where it dominates.  That is the PAHs
  !--- of astrodust and DL07, whose products carry Q_abs alone for pah_neu and
  !--- pah_ion; Zubko's default optics are the Camps et al. (2015) tables, whose
  !--- PAH component scatters as the graphite sphere of the same radius.
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
          !--- No extinction curve was loaded at all.  build_grain_model asks
          !--- build_dust for one in every run -- named (par%sed_kext) or the
          !--- model's own /kext -- and a curve that fails to load fails the
          !--- build there, so this cannot be reached through that path; it is
          !--- kept because dust_extinction documents the status.
          write(msg,'(5a)') 'the ''', trim(par%dust_model), &
             ''' model carries no extinction curve: none was loaded from '// &
             'par%sed_data_dir = ''', trim(par%sed_data_dir), &
             ''' (run SEDust/populate_data.sh, or name one with par%sed_kext)'
       case (3)
          !--- The curve does not cover the model grid.  Both come from the same
          !--- product unless par%sed_kext names another file, which is then the
          !--- one thing to look at.
          write(msg,'(a,es10.3,a,es10.3,3a)') &
             'the model grid [um] = [', lam(1), ', ', lam(n), &
             '] reaches outside the extinction curve ''', &
             trim(dmodel%kext_path), &
             ''' -- no extrapolation: leave par%sed_kext blank to take the '// &
             'model''s own curve, or name one that spans the grid'
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
