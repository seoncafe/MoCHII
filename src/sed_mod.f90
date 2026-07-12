! MoCHII: copied from MoCafe_v2.00/src/sed_mod.f90 (2026-07-11)
module sed_mod
!--- Multi-wavelength (SED) infrastructure for MoCafe v2.00 (Stage 1).
!--- Provides: a log-spaced wavelength grid; wavelength-dependent dust
!--- properties C_ext(lambda), albedo(lambda), g(lambda) read from an
!--- extinction table (e.g. SEDust calc_kext_astrodust.x output); and a
!--- stellar source spectrum (Planck or 2-column file) sampled per photon
!--- with an alias table.
!---
!--- Transport of a photon at wavelength bin il uses the grey-rescaling
!--- factor sed_sext(il) = C_ext(il)/C_ext(lambda_ref) (photon%s_ext):
!--- the grid rhokap is the opacity at the reference wavelength, and every
!--- optical depth is scaled by s_ext (same mathematics as the Jonsson 2006
!--- tau scan), so the raytrace routines (car/clump/amr) need no changes.
  use define
  use random,  only : random_alias_setup, rand_alias_choise
  implicit none
  public

  integer :: sed_nlam = 0
  real(kind=wp), allocatable :: sed_wave(:)     ! bin centers [um]
  real(kind=wp), allocatable :: sed_dwave(:)    ! bin widths  [um]
  real(kind=wp), allocatable :: sed_cext(:)     ! C_ext/H [cm^2/H] at bin centers
  real(kind=wp), allocatable :: sed_albedo(:)   ! dust albedo(lambda)
  real(kind=wp), allocatable :: sed_hgg(:)      ! asymmetry g(lambda)
  real(kind=wp), allocatable :: sed_sext(:)     ! C_ext(lambda)/C_ext(lambda_ref)
  real(kind=wp), allocatable :: sed_lum(:)      ! source luminosity fraction per bin (sum = 1)
  real(kind=wp), allocatable :: sed_src_pdf(:)  ! alias probability table
  integer,       allocatable :: sed_src_alias(:)
  real(kind=wp) :: sed_cext_ref = 0.0_wp        ! C_ext/H at par%lambda_ref

contains
  !---------------------------------------------------------------
  subroutine setup_sed()
  use mpi
  implicit none
  ! local variables
  real(kind=wp), allocatable :: tb_lam(:), tb_alb(:), tb_cos(:), tb_cext(:)
  real(kind=wp), allocatable :: sp_lam(:), sp_lum(:)
  real(kind=wp), allocatable :: edge(:)
  real(kind=wp) :: dlnlam, lum_sum
  integer       :: ntab, nsp, il, ierr

  if (par%nlambda < 2) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: par%nlambda must be >= 2 in SED mode.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (.not. (par%lambda_max > par%lambda_min .and. par%lambda_min > 0.0_wp)) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: require 0 < lambda_min < lambda_max in SED mode.'
     call MPI_FINALIZE(ierr);  stop
  endif
  if (len_trim(par%kext_file) == 0) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: SED mode requires par%kext_file (lambda, albedo, <cos>, C_ext/H table).'
     call MPI_FINALIZE(ierr);  stop
  endif

  !--- log-spaced wavelength grid (bin edges and geometric bin centers).
  sed_nlam = par%nlambda
  allocate(edge(sed_nlam+1))
  allocate(sed_wave(sed_nlam), sed_dwave(sed_nlam))
  allocate(sed_cext(sed_nlam), sed_albedo(sed_nlam), sed_hgg(sed_nlam), sed_sext(sed_nlam))
  allocate(sed_lum(sed_nlam), sed_src_pdf(sed_nlam), sed_src_alias(sed_nlam))
  dlnlam = log(par%lambda_max/par%lambda_min)/sed_nlam
  do il = 1, sed_nlam+1
     edge(il) = par%lambda_min * exp((il-1)*dlnlam)
  enddo
  do il = 1, sed_nlam
     sed_wave(il)  = sqrt(edge(il)*edge(il+1))
     sed_dwave(il) = edge(il+1) - edge(il)
  enddo

  !--- read the dust extinction table (rank 0) and broadcast.
  if (mpar%p_rank == 0) call read_kext_table(trim(par%kext_file), tb_lam, tb_alb, tb_cos, tb_cext, ntab)
  call MPI_BCAST(ntab, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (mpar%p_rank /= 0) allocate(tb_lam(ntab), tb_alb(ntab), tb_cos(ntab), tb_cext(ntab))
  call MPI_BCAST(tb_lam,  ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tb_alb,  ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tb_cos,  ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(tb_cext, ntab, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  !--- interpolate onto the RT wavelength grid: C_ext log-log, albedo/g
  !--- linear in ln(lambda); clamp outside the table range.
  do il = 1, sed_nlam
     sed_cext(il)   = exp(interp_clamped(log(tb_lam), log(tb_cext), log(sed_wave(il))))
     sed_albedo(il) = interp_clamped(log(tb_lam), tb_alb, log(sed_wave(il)))
     sed_hgg(il)    = interp_clamped(log(tb_lam), tb_cos, log(sed_wave(il)))
  enddo
  sed_cext_ref  = exp(interp_clamped(log(tb_lam), log(tb_cext), log(par%lambda_ref)))
  sed_sext(:)   = sed_cext(:)/sed_cext_ref
  !--- the grid opacity (rhokap) is defined at the reference wavelength.
  par%cext_dust = sed_cext_ref
  par%albedo    = interp_clamped(log(tb_lam), tb_alb, log(par%lambda_ref))
  par%hgg       = interp_clamped(log(tb_lam), tb_cos, log(par%lambda_ref))
  par%lambda0   = par%lambda_ref

  !--- source spectrum.  With multiple source components (par%nsource > 1) the
  !--- spectra for each source are set up in sources_mod, so the global single-source
  !--- spectrum is optional here: build a flat placeholder and skip the checks.
  if (par%nsource > 1 .and. len_trim(par%source_spectrum) == 0 .and. par%tstar <= 0.0_wp) then
     sed_lum(:) = sed_dwave(:)
  else if (len_trim(par%source_spectrum) > 0) then
     if (mpar%p_rank == 0) call read_spectrum_file(trim(par%source_spectrum), sp_lam, sp_lum, nsp)
     call MPI_BCAST(nsp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     if (mpar%p_rank /= 0) allocate(sp_lam(nsp), sp_lum(nsp))
     call MPI_BCAST(sp_lam, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(sp_lum, nsp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     do il = 1, sed_nlam
        !--- outside the tabulated spectrum: zero luminosity (no clamping).
        if (sed_wave(il) < sp_lam(1) .or. sed_wave(il) > sp_lam(nsp)) then
           sed_lum(il) = 0.0_wp
        else
           sed_lum(il) = interp_clamped(log(sp_lam), sp_lum, log(sed_wave(il))) * sed_dwave(il)
        endif
     enddo
     deallocate(sp_lam, sp_lum)
  else if (par%tstar > 0.0_wp) then
     do il = 1, sed_nlam
        sed_lum(il) = planck_shape(sed_wave(il), par%tstar) * sed_dwave(il)
     enddo
  else
     if (mpar%p_rank == 0) write(*,'(a)') &
        'ERROR: SED mode requires a source spectrum: par%source_spectrum (file) or par%tstar (Planck, K).'
     call MPI_FINALIZE(ierr);  stop
  endif

  lum_sum = sum(sed_lum)
  if (.not. (lum_sum > 0.0_wp)) then
     if (mpar%p_rank == 0) write(*,'(a)') 'ERROR: source spectrum has zero luminosity on the wavelength grid.'
     call MPI_FINALIZE(ierr);  stop
  endif
  sed_lum(:) = sed_lum(:)/lum_sum

  !--- alias table for luminosity-weighted wavelength-bin sampling.
  sed_src_pdf(:) = sed_lum(:)
  call random_alias_setup(sed_src_pdf, sed_src_alias)

  if (mpar%p_rank == 0) then
     write(*,'(a)')          '--- SED (multi-wavelength) mode ---'
     write(*,'(a,i6)')       'N wavelength bins         : ', sed_nlam
     write(*,'(a,2es12.4)')  'lambda_min, lambda_max(um): ', par%lambda_min, par%lambda_max
     write(*,'(a,es12.4)')   'reference lambda (um)     : ', par%lambda_ref
     write(*,'(a,es12.4)')   'C_ext/H at lambda_ref     : ', sed_cext_ref
     write(*,'(a,2f8.4)')    'albedo, g at lambda_ref   : ', par%albedo, par%hgg
     if (par%tstar > 0.0_wp .and. len_trim(par%source_spectrum) == 0) then
        write(*,'(a,es12.4)') 'source: Planck T_star (K) : ', par%tstar
     else
        write(*,'(2a)')       'source spectrum file      : ', trim(par%source_spectrum)
     endif
  endif

  deallocate(tb_lam, tb_alb, tb_cos, tb_cext, edge)
  end subroutine setup_sed

  !---------------------------------------------------------------
  !--- sample a wavelength-bin index from the source spectrum.
  function sample_sed_lambda() result(il)
  implicit none
  integer :: il
  il = rand_alias_choise(sed_src_pdf, sed_src_alias)
  end function sample_sed_lambda

  !---------------------------------------------------------------
  !--- Planck function B_lambda (arbitrary normalization), lambda in um.
  pure function planck_shape(lam_um, T) result(b)
  implicit none
  real(kind=wp), intent(in) :: lam_um, T
  real(kind=wp) :: b, x
  real(kind=wp), parameter :: hc_over_k = 1.43877687750393e4_wp  ! [um K]
  x = hc_over_k/(lam_um*T)
  if (x > 700.0_wp) then
     b = 0.0_wp
  else
     b = 1.0_wp/(lam_um**5 * (exp(x) - 1.0_wp))
  endif
  end function planck_shape

  !---------------------------------------------------------------
  !--- linear interpolation with clamping at the table ends.
  pure function interp_clamped(x, y, xnew) result(ynew)
  implicit none
  real(kind=wp), intent(in) :: x(:), y(:), xnew
  real(kind=wp) :: ynew, t
  integer :: n, i
  n = size(x)
  if (xnew <= x(1)) then
     ynew = y(1)
  else if (xnew >= x(n)) then
     ynew = y(n)
  else
     do i = 2, n
        if (xnew <= x(i)) exit
     enddo
     t    = (xnew - x(i-1))/(x(i) - x(i-1))
     ynew = y(i-1) + t*(y(i) - y(i-1))
  endif
  end function interp_clamped

  !---------------------------------------------------------------
  !--- read the extinction table: comment lines start with '#';
  !--- columns: lambda[um]  albedo  <cos>  C_ext/H [cm^2/H]  (extra columns ignored).
  subroutine read_kext_table(fname, lam, alb, cosg, cext, n)
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=wp), allocatable, intent(out) :: lam(:), alb(:), cosg(:), cext(:)
  integer, intent(out) :: n
  character(len=512) :: line
  integer :: unit, ios, i
  open(newunit=unit, file=fname, status='old', action='read')
  n = 0
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
     n = n + 1
  enddo
  if (n < 2) then
     write(*,'(3a)') 'ERROR: extinction table ', trim(fname), ' has fewer than 2 data rows.'
     stop
  endif
  allocate(lam(n), alb(n), cosg(n), cext(n))
  rewind(unit)
  i = 0
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
     i = i + 1
     read(line,*) lam(i), alb(i), cosg(i), cext(i)
  enddo
  close(unit)
  end subroutine read_kext_table

  !---------------------------------------------------------------
  !--- read a 2-column spectrum file: lambda[um], L_lambda[arbitrary].
  subroutine read_spectrum_file(fname, lam, lum, n)
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=wp), allocatable, intent(out) :: lam(:), lum(:)
  integer, intent(out) :: n
  character(len=512) :: line
  integer :: unit, ios, i
  open(newunit=unit, file=fname, status='old', action='read')
  n = 0
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
     n = n + 1
  enddo
  if (n < 2) then
     write(*,'(3a)') 'ERROR: spectrum file ', trim(fname), ' has fewer than 2 data rows.'
     stop
  endif
  allocate(lam(n), lum(n))
  rewind(unit)
  i = 0
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) exit
     line = adjustl(line)
     if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
     i = i + 1
     read(line,*) lam(i), lum(i)
  enddo
  close(unit)
  end subroutine read_spectrum_file

end module sed_mod
