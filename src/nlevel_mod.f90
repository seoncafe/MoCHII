module nlevel_mod
!---------------------------------------------------------------------------
! MoCHII: generic n-level atom solver (output-time diagnostics).
!
! Reads data/atomic/nlevel_<el>_<stage>.txt (fit_nlevel.py):
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
! Atomic-hydrogen impact (ground-term fine structure of C I, C II, O I):
! nlevel_load also reads data/atomic/hcoll_<el>_<stage>.txt beside the atom
! file when it exists (make_hcoll_fine_structure.py).  Given n_HI, those
! rates join the electron rates in the SAME statistical equilibrium, so the
! populations saturate against both colliders together:
!   R = n_e q^e + n_HI q^H + A.
! Row: lo hi form Tmin Tmax c_1..c_10.  form 1 is Abrahamsson, Krems &
! Dalgarno (2007) eq. (1), the EXCITATION rate lo -> hi,
!   ln(k/1e-11) = sum_{i=0..8} c_(i+2) T^(-i c_1);
! form 2 is a DE-excitation power law q_hi->lo = c_1 (T/100)^c_2.  The
! de-excitation rate is taken at T clamped to [Tmin, Tmax] (the range of the
! calculation; no extrapolation) and excitation follows by detailed balance
! at the actual T with this atom's level energies.
! Every transition is oriented by ENERGY (emitting_level), never by the order
! of its two index columns -- CHIANTI does not order them consistently.
!
! Evaluated per leaf only at output time on the converged state.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: nlevel_load, nlevel_emissivities, nlevel_atom_type
  public :: nlevel_nlines, nlevel_line_ident
  !--- cooling-table siblings (nlevel_cooling_mod): the collisional line
  !--- cooling Lambda_SE(T, ne) = sum_k n_u A dE, split by n_e-dependence so
  !--- the Upsilon work is done once for each T across the whole n_e column.
  public :: nlevel_cooling_prepT, nlevel_cooling_solve, nlevel_cooling_at
  !--- the level populations themselves at (Rcoll, Rrad, ne): the term
  !--- populations of a charge-exchange product ion (species_mod, form 9)
  public :: nlevel_populations
  public :: MAXLEV
  !--- atomic-hydrogen impact: the increment of the line cooling over the
  !--- electron-only populations (species_mod::metal_cooling_H).
  public :: nlevel_has_hcoll, nlevel_hcoll_cooling

  integer, parameter :: MAXLEV = 40, MAXTR = 600, MAXHC = 10

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
     !--- atomic-hydrogen impact rates (hcoll_<el>_<stage>.txt); nhc = 0
     !--- when the atom has none.
     integer :: nhc = 0
     integer :: hc_lo(MAXHC) = 0, hc_hi(MAXHC) = 0, hc_form(MAXHC) = 0
     real(kind=wp) :: hc_tmin(MAXHC) = 0.0_wp, hc_tmax(MAXHC) = 0.0_wp
     real(kind=wp) :: hc_c(10, MAXHC) = 0.0_wp
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
          if (n > MAXLEV) call nlevel_toobig(fname, 'NLEV', n, MAXLEV)
          atom%nlev = n
          do k = 1, n
             read(unit,*) i, atom%g(k), atom%e_cm1(k)
          end do
       case ('NRAD')
          if (n > MAXTR) call nlevel_toobig(fname, 'NRAD', n, MAXTR)
          atom%nrad = n
          do k = 1, n
             read(unit,*) atom%rl(k), atom%ru(k), atom%A(k)
          end do
       case ('NUPS')
          if (n > MAXTR) call nlevel_toobig(fname, 'NUPS', n, MAXTR)
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
    if (ok) call hcoll_load(fname, atom)
  end subroutine nlevel_load

  !=========================================================================
  ! Atomic-hydrogen impact rates from the sibling file: the basename
  ! nlevel_<el>_<stage>[_suffix].txt with 'nlevel_' replaced by 'hcoll_'.
  ! No file -> nhc = 0 (electron impact only).  A malformed file stops the
  ! run: a silently dropped collider would change the populations.
  !=========================================================================
  subroutine hcoll_load(fname, atom)
    implicit none
    character(len=*),       intent(in)    :: fname
    type(nlevel_atom_type), intent(inout) :: atom
    character(len=1024) :: hname, line
    character(len=16)   :: key
    integer :: unit, ios, k, n, ib
    logical :: exists

    ib = index(fname, 'nlevel_', back=.true.)
    if (ib == 0) return
    hname = fname(1:ib-1)//'hcoll_'//fname(ib+7:)
    inquire(file=trim(hname), exist=exists)
    if (.not. exists) return
    open(newunit=unit, file=trim(hname), status='old', iostat=ios)
    if (ios /= 0) call hcoll_bad(hname, 'cannot be opened')
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) key, n
       if (trim(key) /= 'NHCOL') cycle
       if (n > MAXHC) call nlevel_toobig(hname, 'NHCOL', n, MAXHC)
       do k = 1, n
          read(unit,*,iostat=ios) atom%hc_lo(k), atom%hc_hi(k), &
             atom%hc_form(k), atom%hc_tmin(k), atom%hc_tmax(k), atom%hc_c(:,k)
          if (ios /= 0) call hcoll_bad(hname, 'has a short NHCOL row')
          if (min(atom%hc_lo(k), atom%hc_hi(k)) < 1 .or. &
              max(atom%hc_lo(k), atom%hc_hi(k)) > atom%nlev .or. &
              atom%e_cm1(atom%hc_hi(k)) <= atom%e_cm1(atom%hc_lo(k))) &
             call hcoll_bad(hname, 'names a level pair that is not '// &
                'lower -> upper in '//trim(fname))
          if (atom%hc_form(k) /= 1 .and. atom%hc_form(k) /= 2) &
             call hcoll_bad(hname, 'has an unknown form')
       end do
       atom%nhc = n
    end do
    close(unit)
  end subroutine hcoll_load

  subroutine hcoll_bad(hname, what)
    character(len=*), intent(in) :: hname, what
    write(*,'(a)') 'ERROR: '//trim(hname)//' '//what//'.'
    error stop
  end subroutine hcoll_bad

  logical function nlevel_has_hcoll(atom)
    type(nlevel_atom_type), intent(in) :: atom
    nlevel_has_hcoll = atom%nhc > 0
  end function nlevel_has_hcoll

  !=========================================================================
  ! Atomic-hydrogen impact rates without the n_HI factor, added into Rh
  ! (Rh(i,j) = rate j -> i, the convention of the electron matrices).
  !=========================================================================
  subroutine hcoll_rate_matrix(atom, T, Rh)
    implicit none
    type(nlevel_atom_type), intent(in)    :: atom
    real(kind=wp),          intent(in)    :: T
    real(kind=wp),          intent(inout) :: Rh(MAXLEV,MAXLEV)
    real(kind=wp) :: Tc, de_erg, q_ul, q_lu, s
    integer :: k, l, u_, i

    do k = 1, atom%nhc
       l  = atom%hc_lo(k);  u_ = atom%hc_hi(k)
       de_erg = (atom%e_cm1(u_) - atom%e_cm1(l))*HC_ERG_CM
       Tc = min(max(T, atom%hc_tmin(k)), atom%hc_tmax(k))
       select case (atom%hc_form(k))
       case (1)   ! Abrahamsson et al. (2007) eq. (1): excitation at Tc
          s = 0.0_wp
          do i = 9, 1, -1
             s = s*Tc**(-atom%hc_c(1,k)) + atom%hc_c(i+1,k)
          end do
          q_lu = 1.0e-11_wp*exp(s)
          q_ul = q_lu*(atom%g(l)/atom%g(u_))*exp(de_erg/(kboltz_cgs*Tc))
       case default   ! de-excitation power law
          q_ul = atom%hc_c(1,k)*(Tc/100.0_wp)**atom%hc_c(2,k)
       end select
       q_lu = (atom%g(u_)/atom%g(l))*q_ul*exp(-de_erg/(kboltz_cgs*T))
       Rh(u_,l) = Rh(u_,l) + q_lu    ! excitation    l -> u
       Rh(l,u_) = Rh(l,u_) + q_ul    ! de-excitation u -> l
    end do
  end subroutine hcoll_rate_matrix

  !=========================================================================
  subroutine nlevel_toobig(fname, key, n, nmax)
    character(len=*), intent(in) :: fname, key
    integer,          intent(in) :: n, nmax
    write(*,'(a,i0,a,i0,a)') 'ERROR: '//trim(fname)//' '//trim(key)//' = ', &
       n, ' exceeds the nlevel_mod limit ', nmax, ' (raise MAXLEV/MAXTR).'
    error stop
  end subroutine nlevel_toobig

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
  ! Orient a transition by ENERGY.  Spontaneous emission goes downward, so the
  ! emitting level of a radiative row is the one with the higher term energy;
  ! the collisional de-excitation rate q_ul = 8.629e-6 Ups/(g_u sqrt T) is the
  ! rate out of that same level, and the Boltzmann factor of the excitation
  ! rate needs dE >= 0.  The two index columns of the CHIANTI .wgfa/.scups
  ! rows are NOT consistently ordered (O II lists 2D3/2 - 2D5/2, i.e. the
  ! upper level first; N I lists 2D5/2 - 2D3/2, i.e. the lower level first),
  ! so the pair is ordered here instead of being assumed.
  ! Exactly degenerate levels (dE = 0) radiate no energy: the pair is then
  ! ordered by the index arguments so the result stays deterministic, dE = 0
  ! leaves the Boltzmann factor at unity, and the line emissivity vanishes.
  !=========================================================================
  subroutine emitting_level(atom, i1, i2, lo, hi, de_erg)
    implicit none
    type(nlevel_atom_type), intent(in)  :: atom
    integer,                intent(in)  :: i1, i2
    integer,                intent(out) :: lo, hi
    real(kind=wp),          intent(out) :: de_erg
    if (atom%e_cm1(i2) >= atom%e_cm1(i1)) then
       lo = i1;  hi = i2
    else
       lo = i2;  hi = i1
    end if
    de_erg = (atom%e_cm1(hi) - atom%e_cm1(lo))*HC_ERG_CM
  end subroutine emitting_level

  !=========================================================================
  ! Level populations and line emissivities at (T, ne[, nHI]).
  ! emis(k) [erg/s per ion] for each NRAD transition k.  nHI adds the
  ! atomic-hydrogen impact rates when the atom carries them.
  !=========================================================================
  subroutine nlevel_emissivities(atom, T, ne, emis, nHI)
    implicit none
    type(nlevel_atom_type), intent(in)  :: atom
    real(kind=wp),          intent(in)  :: T, ne
    real(kind=wp),          intent(out) :: emis(MAXTR)
    real(kind=wp), intent(in), optional :: nHI
    real(kind=wp) :: Rcoll(MAXLEV,MAXLEV), Rrad(MAXLEV,MAXLEV)
    real(kind=wp) :: Rh(MAXLEV,MAXLEV), R(MAXLEV,MAXLEV), pop(MAXLEV)
    real(kind=wp) :: de_erg
    integer :: k, l, u_

    call nlevel_cooling_prepT(atom, T, Rcoll, Rrad)
    R = ne*Rcoll + Rrad
    if (present(nHI) .and. atom%nhc > 0) then
       Rh = 0.0_wp
       call hcoll_rate_matrix(atom, T, Rh)
       R = R + nHI*Rh
    end if
    call level_populations(atom%nlev, R, pop)
    do k = 1, atom%nrad
       call emitting_level(atom, atom%rl(k), atom%ru(k), l, u_, de_erg)
       emis(k) = pop(u_)*atom%A(k)*de_erg
    end do
  end subroutine nlevel_emissivities

  !=========================================================================
  ! Statistical equilibrium: dn/dt = 0 for the rate matrix R (R(i,j) = rate
  ! j -> i [s^-1]) with the closure sum(n_i) = 1 replacing the first row,
  ! by Gaussian elimination with partial pivoting.
  !=========================================================================
  subroutine level_populations(n, R, pop)
    implicit none
    integer,       intent(in)  :: n
    real(kind=wp), intent(in)  :: R(MAXLEV,MAXLEV)
    real(kind=wp), intent(out) :: pop(MAXLEV)
    real(kind=wp) :: M(MAXLEV,MAXLEV+1), fac
    integer :: i, j, ip

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
    pop = 0.0_wp
    do i = n, 1, -1
       pop(i) = M(i,n+1)
       do j = i+1, n
          pop(i) = pop(i) - M(i,j)*pop(j)
       end do
       pop(i) = pop(i)/M(i,i)
    end do
  end subroutine level_populations

  !=========================================================================
  ! Collisional line cooling COEFFICIENT
  !   Lambda_SE(T, ne) = (sum_k n_u A_ul dE_ul) / ne   [erg cm^3 s^-1],
  ! summed over the NRAD radiative transitions -- the same quantity as
  ! cooling_eval, so the volume cooling is ne n_ion Lambda_SE and the n_e -> 0
  ! limit is the low-density fit.  It decreases with n_e (the populations
  ! saturate), which is the density suppression.  Split into two steps so the
  ! n_e-independent rate matrices (and the Upsilon evaluations inside them)
  ! are built once for each T and reused for the whole n_e column:
  !   nlevel_cooling_prepT -> Rcoll (collisional rates without the n_e factor)
  !                           and Rrad (radiative A-values), matching the rate
  !                           assembly of nlevel_emissivities;
  !   nlevel_cooling_solve -> R = ne*Rcoll + Rrad, statistical equilibrium,
  !                           and the cooling sum.
  ! nlevel_cooling_at combines both (no cache) as the direct reference.
  !=========================================================================
  subroutine nlevel_cooling_prepT(atom, T, Rcoll, Rrad)
    implicit none
    type(nlevel_atom_type), intent(in)  :: atom
    real(kind=wp),          intent(in)  :: T
    real(kind=wp),          intent(out) :: Rcoll(MAXLEV,MAXLEV)
    real(kind=wp),          intent(out) :: Rrad(MAXLEV,MAXLEV)
    real(kind=wp) :: q_ul, q_lu, de_erg, ups
    integer :: k, l, u_

    Rcoll = 0.0_wp
    Rrad  = 0.0_wp
    do k = 1, atom%nups
       call emitting_level(atom, atom%ul_l(k), atom%ul_u(k), l, u_, de_erg)
       ups  = ups_eval(atom, k, T)
       q_ul = COLL_PREF/(atom%g(u_)*sqrt(T))*ups
       q_lu = (atom%g(u_)/atom%g(l))*q_ul*exp(-de_erg/(kboltz_cgs*T))
       Rcoll(u_,l) = Rcoll(u_,l) + q_lu   ! excitation    l -> u
       Rcoll(l,u_) = Rcoll(l,u_) + q_ul   ! de-excitation u -> l
    end do
    do k = 1, atom%nrad
       call emitting_level(atom, atom%rl(k), atom%ru(k), l, u_, de_erg)
       Rrad(l,u_) = Rrad(l,u_) + atom%A(k)   ! radiative decay u -> l
    end do
  end subroutine nlevel_cooling_prepT

  !=========================================================================
  real(kind=wp) function nlevel_cooling_solve(atom, Rcoll, Rrad, ne) result(lam)
    implicit none
    type(nlevel_atom_type), intent(in) :: atom
    real(kind=wp),          intent(in) :: Rcoll(MAXLEV,MAXLEV)
    real(kind=wp),          intent(in) :: Rrad(MAXLEV,MAXLEV)
    real(kind=wp),          intent(in) :: ne
    real(kind=wp) :: R(MAXLEV,MAXLEV), pop(MAXLEV)
    real(kind=wp) :: de_erg
    integer :: k, lo, hi

    R = ne*Rcoll + Rrad
    call level_populations(atom%nlev, R, pop)

    lam = 0.0_wp
    do k = 1, atom%nrad
       call emitting_level(atom, atom%rl(k), atom%ru(k), lo, hi, de_erg)
       lam = lam + pop(hi)*atom%A(k)*de_erg
    end do
    !--- cooling coefficient: divide the cooling of one ion by n_e so the volume
    !--- cooling is ne n_ion Lambda_SE (matching cooling_eval).  ne > 0 always
    !--- here (the table floor is 1e-6 cm^-3).
    lam = lam/ne
  end function nlevel_cooling_solve

  !=========================================================================
  ! Level populations normalized to the ion (sum = 1) at the rate matrices of
  ! nlevel_cooling_prepT and electron density ne, electrons only.
  !=========================================================================
  subroutine nlevel_populations(atom, Rcoll, Rrad, ne, pop)
    implicit none
    type(nlevel_atom_type), intent(in)  :: atom
    real(kind=wp),          intent(in)  :: Rcoll(MAXLEV,MAXLEV)
    real(kind=wp),          intent(in)  :: Rrad(MAXLEV,MAXLEV)
    real(kind=wp),          intent(in)  :: ne
    real(kind=wp),          intent(out) :: pop(MAXLEV)
    real(kind=wp) :: R(MAXLEV,MAXLEV)
    R = ne*Rcoll + Rrad
    call level_populations(atom%nlev, R, pop)
  end subroutine nlevel_populations

  !=========================================================================
  real(kind=wp) function nlevel_cooling_at(atom, T, ne) result(lam)
    implicit none
    type(nlevel_atom_type), intent(in) :: atom
    real(kind=wp),          intent(in) :: T, ne
    real(kind=wp) :: Rcoll(MAXLEV,MAXLEV), Rrad(MAXLEV,MAXLEV)
    call nlevel_cooling_prepT(atom, T, Rcoll, Rrad)
    lam = nlevel_cooling_solve(atom, Rcoll, Rrad, ne)
  end function nlevel_cooling_at

  !=========================================================================
  ! Atomic-hydrogen impact increment of the line cooling of one ion
  ! [erg/s per ion]:
  !   dP = sum_k [n_u(e + H) - n_u(e)] A_ul dE_ul,
  ! the radiated power with both colliders in one statistical equilibrium
  ! minus that with electrons alone at the same (T, n_e).  The thermal
  ! balance carries the electron part from the (T, n_e) table
  ! (nlevel_cooling_mod), so adding n_ion dP makes the total the joint
  ! solve: the two colliders saturate the populations together instead of
  ! each against its own critical density.  0 when the atom has no hydrogen
  ! rates or nHI <= 0.
  !=========================================================================
  real(kind=wp) function nlevel_hcoll_cooling(atom, T, ne, nHI) result(dP)
    implicit none
    type(nlevel_atom_type), intent(in) :: atom
    real(kind=wp),          intent(in) :: T, ne, nHI
    real(kind=wp) :: Rcoll(MAXLEV,MAXLEV), Rrad(MAXLEV,MAXLEV)
    real(kind=wp) :: Rh(MAXLEV,MAXLEV), R(MAXLEV,MAXLEV)
    real(kind=wp) :: pe(MAXLEV), peh(MAXLEV), de_erg
    integer :: k, lo, hi

    dP = 0.0_wp
    if (atom%nhc == 0 .or. nHI <= 0.0_wp) return
    call nlevel_cooling_prepT(atom, T, Rcoll, Rrad)
    Rh = 0.0_wp
    call hcoll_rate_matrix(atom, T, Rh)
    R = ne*Rcoll + Rrad
    call level_populations(atom%nlev, R, pe)
    R = R + nHI*Rh
    call level_populations(atom%nlev, R, peh)
    do k = 1, atom%nrad
       call emitting_level(atom, atom%rl(k), atom%ru(k), lo, hi, de_erg)
       dP = dP + (peh(hi) - pe(hi))*atom%A(k)*de_erg
    end do
  end function nlevel_hcoll_cooling

  !=========================================================================
  integer function nlevel_nlines(atom)
    type(nlevel_atom_type), intent(in) :: atom
    nlevel_nlines = atom%nrad
  end function nlevel_nlines

  !=========================================================================
  ! Line identification: vacuum wavelength [Angstrom] of NRAD transition k,
  ! from |E_u - E_l| so the two index columns may come in either order.
  ! Exactly degenerate levels emit no photon; wl_A = 0 flags that row (its
  ! emissivity is 0 as well) instead of dividing by zero.
  !=========================================================================
  real(kind=wp) function nlevel_line_ident(atom, k) result(wl_A)
    type(nlevel_atom_type), intent(in) :: atom
    integer,                intent(in) :: k
    real(kind=wp) :: dcm
    dcm = abs(atom%e_cm1(atom%ru(k)) - atom%e_cm1(atom%rl(k)))
    if (dcm > 0.0_wp) then
       wl_A = 1.0e8_wp/dcm
    else
       wl_A = 0.0_wp
    end if
  end function nlevel_line_ident

end module nlevel_mod
