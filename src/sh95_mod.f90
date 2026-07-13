module sh95_mod
!---------------------------------------------------------------------------
! MoCHII: case-B recombination line emissivities on a (T, n_e) grid:
! H I (ion 1, sh95_hi_caseB.txt) and He II (ion 2, sh95_heii_caseB.txt)
! from Storey & Hummer (1995, MNRAS 272, 41) via
! tools/fitting/make_sh95_lines.py; He I (ion 3, hei_porter_caseB.txt)
! from Porter et al. (2012, 2013) via tools/fitting/make_hei_lines.py.
! Bilinear interpolation in (log T, log n_e), clamped at the grid edges.
! sh95_emis returns 4 pi j / (n_e n_ion) [erg cm^3 s^-1] with n_ion =
! n_p (H I), n_HeIII (He II) or n_HeII (He I).  The ion argument is
! optional and defaults to 1 (H I) for the original callers.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: sh95_setup, sh95_nlines, sh95_label, sh95_wl, sh95_emis

  integer, parameter :: MAXL = 12

  type sh95_tab_type
     integer :: nT = 0, nD = 0, nlin = 0
     real(kind=wp), allocatable :: tgrid(:), dgrid(:), etab(:,:,:)
     character(len=12) :: labels(MAXL)
     real(kind=wp)     :: wls(MAXL)
  end type sh95_tab_type
  type(sh95_tab_type) :: tab(3)

  !--- warn once (not each of the millions of leaf x line evaluations) when
  !--- T exceeds the 30000 K table maximum: the H line emissivity is then
  !--- frozen at the grid edge, which under the default te_max = 50000 K can
  !--- affect hot cells.  Low-T and low/high-ne clamps are benign.
  logical :: sh95_warned_thi = .false.

contains

  integer function ion_of(ion)
    integer, intent(in), optional :: ion
    ion_of = 1
    if (present(ion)) ion_of = ion
  end function ion_of

  integer function sh95_nlines(ion)
    integer, intent(in), optional :: ion
    sh95_nlines = tab(ion_of(ion))%nlin
  end function sh95_nlines

  character(len=12) function sh95_label(k, ion)
    integer, intent(in) :: k
    integer, intent(in), optional :: ion
    sh95_label = tab(ion_of(ion))%labels(k)
  end function sh95_label

  real(kind=wp) function sh95_wl(k, ion)
    integer, intent(in) :: k
    integer, intent(in), optional :: ion
    sh95_wl = tab(ion_of(ion))%wls(k)
  end function sh95_wl

  !=========================================================================
  subroutine sh95_setup()
    use mpi
    implicit none
    call sh95_load(trim(par%atomic_dir)//'/sh95_hi_caseB.txt', tab(1), &
                   'H I', required=.true.)
    call sh95_load(trim(par%atomic_dir)//'/sh95_heii_caseB.txt', tab(2), &
                   'He II', required=.false.)
    call sh95_load(trim(par%atomic_dir)//'/hei_porter_caseB.txt', tab(3), &
                   'He I', required=.false.)
  end subroutine sh95_setup

  !=========================================================================
  subroutine sh95_load(fname, t, name, required)
    use mpi
    implicit none
    character(len=*),    intent(in)    :: fname, name
    type(sh95_tab_type), intent(inout) :: t
    logical,             intent(in)    :: required
    character(len=256) :: line
    character(len=16)  :: key
    integer :: unit, ios, ia, nu_, nl_, ierr

    if (t%nlin > 0) return                     ! already loaded
    open(newunit=unit, file=fname, status='old', iostat=ios)
    if (ios /= 0) then
       if (required) then
          if (mpar%p_rank == 0) write(*,'(3a)') 'ERROR: cannot open ', &
             trim(fname), ' (par%atomic_dir)'
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
       return
    end if
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) key
       select case (trim(key))
       case ('GRID')
          read(line,*) key, t%nT, t%nD
          allocate(t%tgrid(t%nT), t%dgrid(t%nD), t%etab(t%nT, t%nD, MAXL))
       case ('T')
          read(line,*) key, t%tgrid
       case ('NE')
          read(line,*) key, t%dgrid
       case ('LINE')
          t%nlin = t%nlin + 1
          read(line,*) key, t%labels(t%nlin), nu_, nl_, t%wls(t%nlin)
          do ia = 1, t%nT
             read(unit,*) t%etab(ia, :, t%nlin)
          end do
       end select
    end do
    close(unit)
    if (mpar%p_rank == 0) write(*,'(4a,i0,a)') &
       ' SH95: ', trim(name), ' case-B emissivities loaded', &
       ' (', t%nlin, ' lines)'
  end subroutine sh95_load

  !=========================================================================
  real(kind=wp) function sh95_emis(k, T, ne, ion) result(e)
    implicit none
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: T, ne
    integer,       intent(in), optional :: ion
    real(kind=wp) :: lt, ld, wT, wD
    integer :: iT, iD, io

    io = ion_of(ion)
    associate(nT => tab(io)%nT, nD => tab(io)%nD, &
              tgrid => tab(io)%tgrid, dgrid => tab(io)%dgrid, &
              etab => tab(io)%etab)
    if (T > tgrid(nT) .and. .not. sh95_warned_thi) then
       sh95_warned_thi = .true.
       if (mpar%p_rank == 0) write(*,'(a,f0.0,a,f0.0,a)') &
          ' SH95: WARNING - T = ', T, ' K exceeds the table max ', tgrid(nT), &
          ' K; H/He II line emissivity frozen at the grid edge (shown once).'
    end if
    lt = min(max(T,  tgrid(1)), tgrid(nT))
    ld = min(max(ne, dgrid(1)), dgrid(nD))
    do iT = 1, nT-1
       if (lt < tgrid(iT+1)) exit
    end do
    do iD = 1, nD-1
       if (ld < dgrid(iD+1)) exit
    end do
    wT = log(lt/tgrid(iT)) / log(tgrid(iT+1)/tgrid(iT))
    wD = log(ld/dgrid(iD)) / log(dgrid(iD+1)/dgrid(iD))
    e = exp( log(etab(iT,  iD,  k))*(1.0_wp-wT)*(1.0_wp-wD) &
           + log(etab(iT+1,iD,  k))*wT*(1.0_wp-wD) &
           + log(etab(iT,  iD+1,k))*(1.0_wp-wT)*wD &
           + log(etab(iT+1,iD+1,k))*wT*wD )
    end associate
  end function sh95_emis

end module sh95_mod
