module nlevel_mod
!---------------------------------------------------------------------------
! MoCHII: generic n-level atom solver (Tier 2, output-time diagnostics).
!
! Reads data/atomic/nlevel_<el>_<stage>.txt (fit_nlevel_tier2.py):
!   NLEV  index  g  E[cm^-1]
!   NRAD  l u A_ul
!   NUPS  l u type C dE[Ry] st_lo st_hi logf ncheb c_0..c_{n-1} maxerr
! and evaluates EXACTLY the recipe in those headers (the reference
! implementation is tools/fitting/verify_nlevel_pyneb.py):
!   et = kT/dE ;  st = 1 - ln(C)/ln(et+C)  (type 1,4)
!                 st = et/(et+C)           (type 2,3,5,6)
!   clip st to [st_lo, st_hi];  u = 2(st-st_lo)/(st_hi-st_lo) - 1
!   y = Chebyshev series;  logf=1 -> y = 10^y
!   descale: y ln(et+e) [1]; y [2]; y/(et+1) [3]; y ln(et+C) [4];
!            y/et [5]; 10^y [6]
! Statistical equilibrium at (T_e, n_e): collisional rates from Upsilon
! (de-excitation q_ul = 8.629e-6 Ups/(g_u sqrt T), excitation by detailed
! balance) + radiative A; closure sum(n_i) = 1 (Gaussian elimination).
! Line emissivity: j_ul = n_u A_ul dE_ul [erg/s per ion].
!
! Evaluated per leaf only at output time on the converged state
! (docs/PLAN.md section 8, Tier 2).
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: nlevel_load, nlevel_emissivities, nlevel_atom_type
  public :: nlevel_nlines, nlevel_line_ident

  integer, parameter :: MAXLEV = 16, MAXTR = 100

  type nlevel_atom_type
     logical :: loaded = .false.
     integer :: nlev = 0
     real(kind=wp) :: g(MAXLEV) = 0.0_wp, e_cm1(MAXLEV) = 0.0_wp
     integer :: nrad = 0
     integer :: rl(MAXTR) = 0, ru(MAXTR) = 0
     real(kind=wp) :: A(MAXTR) = 0.0_wp
     integer :: nups = 0
     integer :: ul_l(MAXTR) = 0, ul_u(MAXTR) = 0, ttype(MAXTR) = 0
     integer :: logf(MAXTR) = 0, nc(MAXTR) = 0
     real(kind=wp) :: cups(MAXTR) = 0.0_wp, dery(MAXTR) = 0.0_wp
     real(kind=wp) :: st_lo(MAXTR) = 0.0_wp, st_hi(MAXTR) = 0.0_wp
     real(kind=wp) :: coef(16, MAXTR) = 0.0_wp
  end type nlevel_atom_type

  real(kind=wp), parameter :: COLL_PREF = 8.629e-6_wp
  real(kind=wp), parameter :: KB_OVER_RY = kboltz_cgs/2.1798723611035e-11_wp
  real(kind=wp), parameter :: HC_ERG_CM = 1.9864458571489287e-16_wp

contains

  !=========================================================================
  subroutine nlevel_load(fname, atom, ok)
    implicit none
    character(len=*),       intent(in)  :: fname
    type(nlevel_atom_type), intent(out) :: atom
    logical,                intent(out) :: ok
    character(len=1024) :: line
    character(len=16)   :: key
    integer :: unit, ios, i, k, n
    real(kind=wp) :: err

    ok = .false.
    open(newunit=unit, file=fname, status='old', iostat=ios)
    if (ios /= 0) return
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) key, n
       select case (trim(key))
       case ('NLEV')
          atom%nlev = n
          do k = 1, n
             read(unit,*) i, atom%g(k), atom%e_cm1(k)
          end do
       case ('NRAD')
          atom%nrad = n
          do k = 1, n
             read(unit,*) atom%rl(k), atom%ru(k), atom%A(k)
          end do
       case ('NUPS')
          atom%nups = n
          do k = 1, n
             read(unit,'(a)') line
             read(line,*) atom%ul_l(k), atom%ul_u(k), atom%ttype(k), &
                atom%cups(k), atom%dery(k), atom%st_lo(k), atom%st_hi(k), &
                atom%logf(k), atom%nc(k)
             read(line,*) atom%ul_l(k), atom%ul_u(k), atom%ttype(k), &
                atom%cups(k), atom%dery(k), atom%st_lo(k), atom%st_hi(k), &
                atom%logf(k), atom%nc(k), atom%coef(1:atom%nc(k),k)
          end do
       end select
    end do
    close(unit)
    atom%loaded = atom%nlev > 0 .and. atom%nups > 0
    ok = atom%loaded
  end subroutine nlevel_load

  !=========================================================================
  ! Upsilon(T) from the stored Chebyshev fit (the header recipe).
  !=========================================================================
  real(kind=wp) function ups_eval(atom, k, T) result(ups)
    implicit none
    type(nlevel_atom_type), intent(in) :: atom
    integer,                intent(in) :: k
    real(kind=wp),          intent(in) :: T
    real(kind=wp) :: et, st, u, y, b0, b1, b2
    integer :: j

    et = KB_OVER_RY*T/atom%dery(k)
    if (atom%ttype(k) == 1 .or. atom%ttype(k) == 4) then
       st = 1.0_wp - log(atom%cups(k))/log(et + atom%cups(k))
    else
       st = et/(et + atom%cups(k))
    end if
    st = min(max(st, atom%st_lo(k)), atom%st_hi(k))
    if (atom%st_hi(k) > atom%st_lo(k)) then
       u = 2.0_wp*(st - atom%st_lo(k))/(atom%st_hi(k) - atom%st_lo(k)) - 1.0_wp
    else
       u = 0.0_wp
    end if
    !--- Clenshaw evaluation of the Chebyshev series
    b0 = 0.0_wp;  b1 = 0.0_wp
    do j = atom%nc(k), 2, -1
       b2 = b1;  b1 = b0
       b0 = 2.0_wp*u*b1 - b2 + atom%coef(j,k)
    end do
    y = u*b0 - b1 + atom%coef(1,k)
    if (atom%logf(k) == 1) y = 10.0_wp**y
    select case (atom%ttype(k))
    case (1);  ups = y*log(et + exp(1.0_wp))
    case (2);  ups = y
    case (3);  ups = y/(et + 1.0_wp)
    case (4);  ups = y*log(et + atom%cups(k))
    case (5);  ups = y/et
    case default;  ups = 10.0_wp**y
    end select
    ups = max(ups, 0.0_wp)
  end function ups_eval

  !=========================================================================
  ! Level populations and line emissivities at (T, ne).
  ! emis(k) [erg/s per ion] for each NRAD transition k.
  !=========================================================================
  subroutine nlevel_emissivities(atom, T, ne, emis)
    implicit none
    type(nlevel_atom_type), intent(in)  :: atom
    real(kind=wp),          intent(in)  :: T, ne
    real(kind=wp),          intent(out) :: emis(MAXTR)
    real(kind=wp) :: R(MAXLEV,MAXLEV), M(MAXLEV,MAXLEV+1), pop(MAXLEV)
    real(kind=wp) :: q_ul, q_lu, de_erg, ups, fac
    integer :: k, l, u_, i, j, n, ip

    n = atom%nlev
    R = 0.0_wp
    do k = 1, atom%nups
       l = atom%ul_l(k);  u_ = atom%ul_u(k)
       de_erg = (atom%e_cm1(u_) - atom%e_cm1(l))*HC_ERG_CM
       ups  = ups_eval(atom, k, T)
       q_ul = COLL_PREF/(atom%g(u_)*sqrt(T))*ups
       q_lu = (atom%g(u_)/atom%g(l))*q_ul*exp(-de_erg/(kboltz_cgs*T))
       R(u_,l) = R(u_,l) + ne*q_lu
       R(l,u_) = R(l,u_) + ne*q_ul
    end do
    do k = 1, atom%nrad
       R(atom%rl(k), atom%ru(k)) = R(atom%rl(k), atom%ru(k)) + atom%A(k)
    end do

    !--- M n = b with rate matrix (dn/dt = 0) + closure row.
    M = 0.0_wp
    do i = 1, n
       do j = 1, n
          if (i /= j) then
             M(i,j) = R(i,j)              ! rate j -> i
             M(i,i) = M(i,i) - R(j,i)     ! loss from i
          end if
       end do
    end do
    M(1,1:n) = 1.0_wp
    M(1,n+1) = 1.0_wp
    !--- Gaussian elimination with partial pivoting.
    do i = 1, n-1
       ip = i
       do j = i+1, n
          if (abs(M(j,i)) > abs(M(ip,i))) ip = j
       end do
       if (ip /= i) then
          do j = i, n+1
             fac = M(i,j);  M(i,j) = M(ip,j);  M(ip,j) = fac
          end do
       end if
       do j = i+1, n
          fac = M(j,i)/M(i,i)
          M(j,i:n+1) = M(j,i:n+1) - fac*M(i,i:n+1)
       end do
    end do
    do i = n, 1, -1
       pop(i) = M(i,n+1)
       do j = i+1, n
          pop(i) = pop(i) - M(i,j)*pop(j)
       end do
       pop(i) = pop(i)/M(i,i)
    end do

    do k = 1, atom%nrad
       de_erg = (atom%e_cm1(atom%ru(k)) - atom%e_cm1(atom%rl(k)))*HC_ERG_CM
       emis(k) = pop(atom%ru(k))*atom%A(k)*de_erg
    end do
  end subroutine nlevel_emissivities

  !=========================================================================
  integer function nlevel_nlines(atom)
    type(nlevel_atom_type), intent(in) :: atom
    nlevel_nlines = atom%nrad
  end function nlevel_nlines

  !=========================================================================
  ! Line identification: vacuum wavelength [Angstrom] of NRAD transition k.
  !=========================================================================
  real(kind=wp) function nlevel_line_ident(atom, k) result(wl_A)
    type(nlevel_atom_type), intent(in) :: atom
    integer,                intent(in) :: k
    wl_A = 1.0e8_wp/(atom%e_cm1(atom%ru(k)) - atom%e_cm1(atom%rl(k)))
  end function nlevel_line_ident

end module nlevel_mod
