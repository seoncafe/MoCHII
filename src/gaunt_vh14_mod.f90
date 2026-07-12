module gaunt_vh14_mod
!---------------------------------------------------------------------------
! MoCHII: van Hoof et al. (2014, MNRAS 444, 420) thermally averaged
! free-free Gaunt factors.  Reads data/gauntff_vh14.dat (downloaded from
! the published table set, 2026-07-12): g_ff on an (log10 gam2, log10 u)
! grid with gam2 = Z^2 Ry/kT and u = h nu/kT; 81 x 146 points, both axes
! starting at (-6, -16) with step 0.2 dex.  Bilinear interpolation,
! clamped at the grid edges.  Selected by par%gaunt_vh14 at the same
! call sites as the ported Hummer (1988) getGauntFF (the PLAN item-3
! swap); the Hummer path remains the default so the recorded gates
! reproduce.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: gaunt_vh14_setup, gauntff_vh14

  real(kind=wp), parameter :: RY_OVER_K = 157807.4_wp   ! Ry/k_B [K]
  integer  :: ng = 0, nu_ = 0
  real(kind=wp) :: lg0, lu0, dstep
  real(kind=wp), allocatable :: gtab(:,:)               ! (ng, nu_)

contains

  !=========================================================================
  subroutine gaunt_vh14_setup()
    use mpi
    implicit none
    character(len=256) :: line
    integer :: unit, ios, ierr, iu, magic

    if (ng > 0) return
    open(newunit=unit, file='../../data/gauntff_vh14.dat', status='old', &
         iostat=ios)
    if (ios /= 0) open(newunit=unit, file='data/gauntff_vh14.dat', &
                       status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: par%gaunt_vh14 needs data/gauntff_vh14.dat'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    !--- header: comment lines, then magic, ng nu, lg0, lu0, step (each
    !--- followed by a trailing comment).
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) magic
       exit
    end do
    call read_noncomment(unit, line);  read(line,*) ng, nu_
    call read_noncomment(unit, line);  read(line,*) lg0
    call read_noncomment(unit, line);  read(line,*) lu0
    call read_noncomment(unit, line);  read(line,*) dstep
    allocate(gtab(ng, nu_))
    do iu = 1, nu_
       call read_values(unit, gtab(:, iu))
    end do
    close(unit)
    if (mpar%p_rank == 0) write(*,'(a,i4,a,i4,a)') &
       ' GFF: van Hoof et al. (2014) Gaunt table loaded (', ng, ' x', &
       nu_, ' grid)'
  end subroutine gaunt_vh14_setup

  subroutine read_noncomment(unit, line)
    integer,          intent(in)  :: unit
    character(len=*), intent(out) :: line
    integer :: ios
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) return
       line = adjustl(line)
       if (len_trim(line) > 0 .and. line(1:1) /= '#') return
    end do
  end subroutine read_noncomment

  !--- read ng values spanning multiple lines, skipping comments.
  subroutine read_values(unit, row)
    integer,       intent(in)  :: unit
    real(kind=wp), intent(out) :: row(:)
    character(len=4096) :: line
    integer :: got, ios, ntok
    real(kind=wp) :: tmp(256)
    got = 0
    do while (got < size(row))
       call read_noncomment(unit, line)
       ntok = count_tokens(line)
       read(line,*,iostat=ios) tmp(1:ntok)
       row(got+1:got+ntok) = tmp(1:ntok)
       got = got + ntok
    end do
  end subroutine read_values

  integer function count_tokens(line) result(n)
    character(len=*), intent(in) :: line
    logical :: insp
    integer :: i
    n = 0;  insp = .true.
    do i = 1, len_trim(line)
       if (line(i:i) == ' ') then
          insp = .true.
       else if (insp) then
          n = n + 1;  insp = .false.
       end if
    end do
  end function count_tokens

  !=========================================================================
  ! g_ff at charge Zc, temperature T [K], photon frequency nu_ryd [Ryd].
  !=========================================================================
  real(kind=wp) function gauntff_vh14(Zc, T, nu_ryd) result(g)
    implicit none
    real(kind=wp), intent(in) :: Zc, T, nu_ryd
    real(kind=wp) :: lg, lu, xg, xu, wg, wu
    integer :: ig, iu

    lg = log10(Zc*Zc*RY_OVER_K/T)
    lu = log10(max(nu_ryd, tinest)*RY_OVER_K/T)
    xg = min(max((lg - lg0)/dstep, 0.0_wp), real(ng-1, wp) - 1.0e-9_wp)
    xu = min(max((lu - lu0)/dstep, 0.0_wp), real(nu_-1, wp) - 1.0e-9_wp)
    ig = int(xg) + 1;  wg = xg - real(ig-1, wp)
    iu = int(xu) + 1;  wu = xu - real(iu-1, wp)
    g = gtab(ig,  iu)*(1.0_wp-wg)*(1.0_wp-wu) &
      + gtab(ig+1,iu)*wg*(1.0_wp-wu) &
      + gtab(ig,  iu+1)*(1.0_wp-wg)*wu &
      + gtab(ig+1,iu+1)*wg*wu
  end function gauntff_vh14

end module gaunt_vh14_mod
