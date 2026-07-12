module sh95_mod
!---------------------------------------------------------------------------
! MoCHII: Storey & Hummer (1995, MNRAS 272, 41) H I case-B recombination
! line emissivities.  Reads data/atomic/sh95_hi_caseB.txt (produced by
! tools/fitting/make_sh95_lines.py from the SH95 e1bx.d data): the
! principal lines on the (T, n_e) grid; bilinear interpolation in
! (log T, log n_e), clamped at the grid edges.
! sh95_emis returns 4 pi j / (n_e n_p) [erg cm^3 s^-1].
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: sh95_setup, sh95_nlines, sh95_label, sh95_wl, sh95_emis

  integer, parameter :: MAXL = 12
  integer :: nT = 0, nD = 0, nlin = 0
  real(kind=wp), allocatable :: tgrid(:), dgrid(:), etab(:,:,:)  ! (nT,nD,nlin)
  character(len=12) :: labels(MAXL)
  real(kind=wp)     :: wls(MAXL)

contains

  integer function sh95_nlines();  sh95_nlines = nlin;  end function
  character(len=12) function sh95_label(k);  integer, intent(in) :: k
    sh95_label = labels(k);  end function
  real(kind=wp) function sh95_wl(k);  integer, intent(in) :: k
    sh95_wl = wls(k);  end function

  !=========================================================================
  subroutine sh95_setup()
    use mpi
    implicit none
    character(len=256) :: line
    character(len=16)  :: key
    integer :: unit, ios, k, ia, nu_, nl_, ierr

    open(newunit=unit, file=trim(par%atomic_dir)//'/sh95_hi_caseB.txt', &
         status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: cannot open sh95_hi_caseB.txt (par%atomic_dir)'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    nlin = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) key
       select case (trim(key))
       case ('GRID')
          read(line,*) key, nT, nD
          allocate(tgrid(nT), dgrid(nD), etab(nT, nD, MAXL))
       case ('T')
          read(line,*) key, tgrid
       case ('NE')
          read(line,*) key, dgrid
       case ('LINE')
          nlin = nlin + 1
          read(line,*) key, labels(nlin), nu_, nl_, wls(nlin)
          do ia = 1, nT
             read(unit,*) etab(ia, :, nlin)
          end do
       end select
    end do
    close(unit)
    if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
       ' SH95: H I case-B emissivities loaded (', nlin, ' lines)'
  end subroutine sh95_setup

  !=========================================================================
  real(kind=wp) function sh95_emis(k, T, ne) result(e)
    implicit none
    integer,       intent(in) :: k
    real(kind=wp), intent(in) :: T, ne
    real(kind=wp) :: lt, ld, wT, wD
    integer :: iT, iD

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
  end function sh95_emis

end module sh95_mod
