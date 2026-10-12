module nlevel_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
!---------------------------------------------------------------------------
! MoCHII: electron excitation, statistical equilibrium, and line emissivity.
!
! Reads nlevel_<el>_<stage>.txt atoms (data/atomic/excitation/<model>/ for
! metals, data/atomic/ for H I) in two forms: Chebyshev fits of Upsilon (the
! legacy form, still used by the H I atom) and NLEVEL_VERSION 2 source-node
! atoms. Both forms use NLEV, NRAD, and NUPS.
! Version 2 adds NODES before NUPS, with BT_SOURCE rows storing original
! scaled nodes and TEMPERATURE_TABLE rows storing physical T and Upsilon.
! An optional BT_TMIN_K record holds source Upsilon constant below its
! stated lower temperature; Maxwellian rates still use the actual T.
! TEMPERATURE_TABLE rows hold Upsilon at their end nodes outside the table
! (nlevel_upsilon); an UPSILON_CONSTANT_UPPER_K record after NUPS, which
! once bounded that hold, is checked and otherwise ignored.
! The legacy fit evaluates the recipe in its file header:
!   et = kT/dE ;  st = 1 - ln(C)/ln(et+C)  (type 1,4)
!                 st = et/(et+C)           (type 2,3,5,6)
!   clip st to [st_lo, st_hi];  u = 2(st-st_lo)/(st_hi-st_lo) - 1
!   y = Chebyshev series;  logf=1 -> y = 10^y
!   descale: y ln(et+e) [1]; y [2]; y/(et+1) [3]; y ln(et+C) [4];
!            y/et [5]; 10^y [6]
! Version 2 evaluates source nodes with not-a-knot cubic or PCHIP in the
! scaled coordinate. Temperature tables use PCHIP in log T and log Upsilon.
! Statistical equilibrium at (T_e, n_e): collisional rates from Upsilon
! (de-excitation q_ul = 8.629e-6 Ups/(g_u sqrt T), excitation by detailed
! balance) + radiative A; closure sum(n_i) = 1 (LU with partial pivoting
! and iterative refinement, level_populations).
! Line emissivity: j_ul = n_u A_ul dE_ul [erg/s per ion].
!
! Atomic-hydrogen impact (ground-term fine structure of C I, C II, O I):
! nlevel_load also reads hcoll_<el>_<stage>.txt beside the selected atom
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
! Used for registered electron cooling, local hydrogen-impact thermal terms,
! charge-transfer product populations, and output line diagnostics.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: nlevel_load, nlevel_emissivities, nlevel_atom_type
  public :: nlevel_nlines, nlevel_line_ident
  public :: metal_excitation_data_file
  !--- cooling-table siblings (nlevel_cooling_mod): the collisional line
  !--- cooling Lambda_SE(T, ne) = sum_k n_u A dE, split by n_e-dependence so
  !--- the Upsilon work is done once for each T across the whole n_e column.
  public :: nlevel_cooling_prepT, nlevel_cooling_solve, nlevel_cooling_at
  !--- the level populations themselves at (Rcoll, Rrad, ne): the term
  !--- populations of a charge-exchange product ion (species_mod, form 9)
  public :: nlevel_populations
  public :: nlevel_upsilon
  !--- atomic-hydrogen impact: the increment of the line cooling over the
  !--- electron-only populations (species_mod::metal_cooling_H).
  public :: nlevel_has_hcoll, nlevel_hcoll_cooling

  type nlevel_atom_type
     logical :: loaded = .false.
     logical :: density_ratio_cooling = .false.
     real(kind=wp) :: bt_tmin_K = 0.0_wp
     integer :: nlev = 0
     real(kind=wp), allocatable :: g(:), e_cm1(:)
     integer :: nrad = 0
     integer, allocatable :: rl(:), ru(:)
     real(kind=wp), allocatable :: A(:)
     integer :: nups = 0
     integer, allocatable :: ul_l(:), ul_u(:), ttype(:), logf(:), nc(:)
     integer, allocatable :: source_model(:), interpolation(:), offset(:)
     real(kind=wp), allocatable :: cups(:), dery(:), st_lo(:), st_hi(:)
     real(kind=wp), allocatable :: coef(:,:)
     real(kind=wp), allocatable :: node_x(:), node_y(:)
     real(kind=wp), allocatable :: poly_a(:), poly_b(:), poly_c(:), poly_d(:)
     !--- atomic-hydrogen impact rates (hcoll_<el>_<stage>.txt); nhc = 0
     !--- when the atom has none.
     integer :: nhc = 0
     integer, allocatable :: hc_lo(:), hc_hi(:), hc_form(:)
     real(kind=wp), allocatable :: hc_tmin(:), hc_tmax(:), hc_c(:,:)
  end type nlevel_atom_type

  integer, parameter :: LEGACY_FIT=0, BT_SOURCE=1, TEMPERATURE_TABLE=2
  integer, parameter :: CUBIC=1, PCHIP=2, LOG_PCHIP_T=3

  real(kind=wp), parameter :: COLL_PREF = 8.629e-6_wp
  real(kind=wp), parameter :: KB_OVER_RY = kboltz_cgs/2.1798723611035e-11_wp
  real(kind=wp), parameter :: HC_ERG_CM = 1.9864458571489287e-16_wp

contains
  function metal_excitation_data_file(quantity, element, stage) result(filename)
    character(len=*), intent(in) :: quantity, element
    integer, intent(in) :: stage
    character(len=PATH_LEN) :: filename, directory
    logical :: exists
    select case(trim(par%metal_excitation_model))
    case('chianti', 'mocassin_v2')
    case default
       error stop 'metal_excitation_model must be chianti or mocassin_v2'
    end select
    if (quantity /= 'nlevel' .and. quantity /= 'cooling') &
       error stop 'Unknown metal excitation data quantity'
    directory=trim(par%atomic_dir)//'/excitation/'//trim(par%metal_excitation_model)
    inquire(file=trim(directory)//'/manifest.json',exist=exists)
    if (.not. exists) error stop 'Requested metal excitation model directory is missing'
    write(filename,'(a,a,a,a,i0,a)') trim(directory)//'/',trim(quantity),'_', &
         trim(element)//'_',stage,'.txt'
  end function metal_excitation_data_file



  subroutine nlevel_load(fname, atom, ok)
    implicit none
    character(len=*),       intent(in)  :: fname
    type(nlevel_atom_type), intent(out) :: atom
    logical,                intent(out) :: ok
    character(len=1024) :: line
    character(len=32) :: key, model, interpolation_name, extra
    real(kind=wp) :: upper_K
    integer :: unit, ios, i, k, n, version, total_nodes, next_node, j
    logical :: first_record, have_nodes, have_density_ratio, have_bt_tmin
    logical :: have_constant_upper
    logical, allocatable :: pair_seen(:,:)

    ok = .false.
    open(newunit=unit, file=fname, status='old', iostat=ios)
    if (ios /= 0) return
    version = 1
    first_record = .true.
    have_nodes = .false.
    have_density_ratio = .false.
    have_bt_tmin = .false.
    have_constant_upper = .false.
    total_nodes = 0
    next_node = 1
    do
       call atom_record(unit, fname, line, ios)
       if (ios < 0) exit
       read(line,*,iostat=ios) key
       if (ios /= 0) call atom_bad(fname, 'invalid section header')
       if (trim(key) == 'UPSILON_CONSTANT_UPPER_K') then
          if (version /= 2 .or. have_constant_upper .or. have_density_ratio .or. &
              .not. allocated(atom%source_model)) call atom_bad(fname, 'invalid constant upper placement')
          if (atom%nups < 1 .or. any(atom%source_model /= TEMPERATURE_TABLE)) &
              call atom_bad(fname, 'constant upper requires temperature table rows')
          ! The record bounded the temperature to which the final-node
          ! Upsilon was held; nlevel_upsilon now holds it at every higher
          ! temperature, so the bound is read and checked but not stored.
          read(line,*,iostat=ios) key, upper_K
          if (ios /= 0) call atom_bad(fname, 'invalid constant upper temperature')
          read(line,*,iostat=ios) key, upper_K, extra
          if (ios == 0) call atom_bad(fname, 'extra constant upper value')
          if (ios > 0) call atom_bad(fname, 'invalid constant upper record')
          if (.not. ieee_is_finite(upper_K) .or. upper_K <= 0.0_wp) &
              call atom_bad(fname, 'invalid constant upper temperature')
          have_constant_upper = .true.
          first_record = .false.
          cycle
       end if
       if (trim(key) == 'BT_TMIN_K') then
          if (version /= 2 .or. have_bt_tmin .or. have_density_ratio .or. &
              .not. allocated(atom%source_model)) call atom_bad(fname, 'invalid BT_TMIN_K placement')
          if (count(atom%source_model == BT_SOURCE) == 0) &
              call atom_bad(fname, 'BT_TMIN_K requires a BT source row')
          read(line,*,iostat=ios) key, atom%bt_tmin_K
          if (ios /= 0) call atom_bad(fname, 'invalid BT_TMIN_K value')
          read(line,*,iostat=ios) key, atom%bt_tmin_K, extra
          if (ios == 0) call atom_bad(fname, 'extra BT_TMIN_K value')
          if (ios > 0) call atom_bad(fname, 'invalid BT_TMIN_K record')
          if (.not. ieee_is_finite(atom%bt_tmin_K)) &
              call atom_bad(fname, 'nonfinite BT_TMIN_K')
          if (atom%bt_tmin_K <= 0.0_wp) &
              call atom_bad(fname, 'nonpositive BT_TMIN_K')
          have_bt_tmin = .true.
          first_record = .false.
          cycle
       end if
       read(line,*,iostat=ios) key, n
       if (ios /= 0) call atom_bad(fname, 'invalid section header')
       if (first_record) then
          if (trim(key) == 'NLEVEL_VERSION') then
             if (n /= 2) call atom_bad(fname, 'unsupported format version')
             version = 2
             first_record = .false.
             cycle
          end if
          first_record = .false.
       end if
       select case (trim(key))
       case ('DENSITY_RATIO')
          if (version /= 2 .or. have_density_ratio .or. (n /= 0 .and. n /= 1)) &
             call atom_bad(fname, 'invalid DENSITY_RATIO')
          atom%density_ratio_cooling = n == 1
          have_density_ratio = .true.
       case ('NLEV')
          if (allocated(atom%g) .or. n < 1) call atom_bad(fname, 'invalid NLEV')
          atom%nlev = n
          allocate(atom%g(n), atom%e_cm1(n))
          do k = 1, n
             call atom_record(unit, fname, line, ios)
             if (ios < 0) call atom_bad(fname, 'truncated level rows')
             read(line,*,iostat=ios) i, atom%g(k), atom%e_cm1(k)
             if (ios /= 0) call atom_bad(fname, 'invalid level row')
             if (i /= k .or. .not. ieee_is_finite(atom%g(k)) .or. &
                 .not. ieee_is_finite(atom%e_cm1(k))) call atom_bad(fname, 'invalid level index or value')
             if (atom%g(k) <= 0.0_wp) call atom_bad(fname, 'nonpositive statistical weight')
          end do
       case ('NRAD')
          if (.not. allocated(atom%g) .or. allocated(atom%rl) .or. n < 0) &
             call atom_bad(fname, 'invalid NRAD')
          atom%nrad = n
          allocate(atom%rl(n), atom%ru(n), atom%A(n))
          do k = 1, n
             call atom_record(unit, fname, line, ios)
             if (ios < 0) call atom_bad(fname, 'truncated radiative rows')
             read(line,*,iostat=ios) atom%rl(k), atom%ru(k), atom%A(k)
             if (ios /= 0) call atom_bad(fname, 'invalid radiative row')
             if (.not. valid_pair(atom, atom%rl(k), atom%ru(k)) .or. &
                 .not. ieee_is_finite(atom%A(k))) call atom_bad(fname, 'invalid radiative pair or rate')
             if (atom%A(k) < 0.0_wp) call atom_bad(fname, 'negative radiative rate')
          end do
       case ('NODES')
          if (version /= 2 .or. have_nodes .or. .not. allocated(atom%rl) .or. n < 0) &
             call atom_bad(fname, 'invalid NODES')
          total_nodes = n
          have_nodes = .true.
          allocate(atom%node_x(n), atom%node_y(n), atom%poly_a(n), &
             atom%poly_b(n), atom%poly_c(n), atom%poly_d(n))
       case ('NUPS')
          if (.not. allocated(atom%rl) .or. allocated(atom%ul_l) .or. n < 1) &
             call atom_bad(fname, 'invalid NUPS')
          if (version == 2 .and. .not. have_nodes) call atom_bad(fname, 'missing NODES')
          atom%nups = n
          allocate(atom%ul_l(n), atom%ul_u(n), atom%ttype(n), atom%cups(n), &
             atom%dery(n), atom%st_lo(n), atom%st_hi(n), atom%logf(n), atom%nc(n), &
             atom%source_model(n), atom%interpolation(n), atom%offset(n))
          atom%source_model = LEGACY_FIT
          atom%interpolation = 0
          atom%offset = 0
          atom%nc = 0
          atom%logf = 0
          if (version == 1) allocate(atom%coef(16,n))
          if (version == 2) then
             allocate(pair_seen(atom%nlev,atom%nlev))
             pair_seen = .false.
          end if
          do k = 1, n
             call atom_record(unit, fname, line, ios)
             if (ios < 0) call atom_bad(fname, 'truncated collision rows')
             if (version == 1) then
                read(line,*,iostat=ios) atom%ul_l(k), atom%ul_u(k), atom%ttype(k), &
                   atom%cups(k), atom%dery(k), atom%st_lo(k), atom%st_hi(k), &
                   atom%logf(k), atom%nc(k)
                if (ios /= 0) call atom_bad(fname, 'invalid legacy collision row')
                if (atom%nc(k) < 1 .or. atom%nc(k) > 16) call atom_bad(fname, 'invalid coefficient count')
                read(line,*,iostat=ios) atom%ul_l(k), atom%ul_u(k), atom%ttype(k), &
                   atom%cups(k), atom%dery(k), atom%st_lo(k), atom%st_hi(k), &
                   atom%logf(k), atom%nc(k), atom%coef(1:atom%nc(k),k)
                if (ios /= 0 .or. .not. all(ieee_is_finite(atom%coef(1:atom%nc(k),k)))) &
                   call atom_bad(fname, 'invalid legacy coefficients')
                if (atom%logf(k) /= 0 .and. atom%logf(k) /= 1) call atom_bad(fname, 'invalid log flag')
                if (.not. ieee_is_finite(atom%st_lo(k)) .or. .not. ieee_is_finite(atom%st_hi(k))) &
                   call atom_bad(fname, 'invalid source interval')
                if (atom%st_lo(k) > atom%st_hi(k)) call atom_bad(fname, 'reversed source interval')
             else
                read(line,*,iostat=ios) atom%ul_l(k), atom%ul_u(k), model, atom%ttype(k), &
                   atom%cups(k), atom%dery(k), interpolation_name, atom%nc(k)
                if (ios /= 0) call atom_bad(fname, 'invalid source collision row')
                if (.not. ieee_is_finite(atom%cups(k)) .or. &
                    .not. ieee_is_finite(atom%dery(k))) call atom_bad(fname, 'nonfinite source parameters')
                select case (trim(model))
                case ('BT_SOURCE'); atom%source_model(k) = BT_SOURCE
                case ('TEMPERATURE_TABLE'); atom%source_model(k) = TEMPERATURE_TABLE
                case default; call atom_bad(fname, 'unknown collision source model')
                end select
                select case (trim(interpolation_name))
                case ('CUBIC'); atom%interpolation(k) = CUBIC
                case ('PCHIP'); atom%interpolation(k) = PCHIP
                case ('LOG_PCHIP_T'); atom%interpolation(k) = LOG_PCHIP_T
                case default; call atom_bad(fname, 'unknown interpolation')
                end select
                if (atom%nc(k) < 2 .or. atom%nc(k) > total_nodes-next_node+1) &
                   call atom_bad(fname, 'invalid source node count')
                atom%offset(k) = next_node
                if (atom%source_model(k) == BT_SOURCE) then
                   if (atom%interpolation(k) == LOG_PCHIP_T .or. &
                       atom%ttype(k) < 1 .or. atom%ttype(k) > 6) &
                      call atom_bad(fname, 'invalid BT source rule')
                else
                   if (atom%interpolation(k) /= LOG_PCHIP_T .or. &
                       atom%ttype(k) /= 0 .or. atom%cups(k) /= 0.0_wp .or. &
                       atom%dery(k) /= 0.0_wp) call atom_bad(fname, 'invalid temperature table rule')
                end if
                do j = next_node, next_node+atom%nc(k)-1
                   call atom_record(unit, fname, line, ios)
                   if (ios < 0) call atom_bad(fname, 'truncated source nodes')
                   read(line,*,iostat=ios) atom%node_x(j), atom%node_y(j)
                   if (ios /= 0) call atom_bad(fname, 'invalid source node')
                   if (.not. ieee_is_finite(atom%node_x(j)) .or. &
                       .not. ieee_is_finite(atom%node_y(j))) call atom_bad(fname, 'nonfinite source node')
                   if (atom%source_model(k) == TEMPERATURE_TABLE) then
                      if (atom%node_x(j) <= 0.0_wp .or. atom%node_y(j) <= 0.0_wp) &
                         call atom_bad(fname, 'nonpositive temperature table node')
                      atom%node_x(j) = log(atom%node_x(j))
                      atom%node_y(j) = log(atom%node_y(j))
                   else if (atom%ttype(k) /= 6 .and. atom%node_y(j) < 0.0_wp) then
                      call atom_bad(fname, 'negative BT source ordinate')
                   end if
                   if (j > next_node) then
                      if (atom%node_x(j) <= atom%node_x(j-1)) call atom_bad(fname, 'unordered source nodes')
                   end if
                end do
                call source_polynomial(atom, k, fname)
                next_node = next_node + atom%nc(k)
             end if
             if (.not. valid_pair(atom, atom%ul_l(k), atom%ul_u(k))) &
                call atom_bad(fname, 'invalid collision level pair')
             if (version == 2) then
                i=min(atom%ul_l(k),atom%ul_u(k))
                j=max(atom%ul_l(k),atom%ul_u(k))
                if (pair_seen(i,j)) call atom_bad(fname, 'duplicate collision level pair')
                pair_seen(i,j)=.true.
             end if
             if (atom%source_model(k) /= TEMPERATURE_TABLE) then
                if (atom%ttype(k) < 1 .or. atom%ttype(k) > 6 .or. &
                    .not. ieee_is_finite(atom%cups(k)) .or. &
                    .not. ieee_is_finite(atom%dery(k))) call atom_bad(fname, 'invalid BT parameters')
                if (atom%cups(k) <= 0.0_wp .or. atom%dery(k) <= 0.0_wp) &
                   call atom_bad(fname, 'nonpositive BT parameters')
             end if
          end do
       case default
          call atom_bad(fname, 'unknown section')
       end select
    end do
    close(unit)
    if (version == 2 .and. next_node-1 /= total_nodes) call atom_bad(fname, 'NODES count mismatch')
    if (.not. allocated(atom%g) .or. .not. allocated(atom%rl) .or. &
        .not. allocated(atom%ul_l)) call atom_bad(fname, 'missing required section')
    atom%loaded = .true.
    ok = .true.
    call hcoll_load(fname, atom)
  end subroutine nlevel_load

  subroutine atom_record(unit, fname, line, ios)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: fname
    character(len=*), intent(out) :: line
    integer, intent(out) :: ios
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) then
          if (ios > 0) call atom_bad(fname, 'cannot read a record')
          return
       end if
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       return
    end do
  end subroutine atom_record

  subroutine atom_bad(fname, what)
    character(len=*), intent(in) :: fname, what
    write(*,'(a)') 'ERROR: '//trim(fname)//' '//what//'.'
    error stop
  end subroutine atom_bad

  logical function valid_pair(atom, i, j)
    type(nlevel_atom_type), intent(in) :: atom
    integer, intent(in) :: i, j
    valid_pair = i >= 1 .and. j >= 1 .and. i <= atom%nlev .and. &
       j <= atom%nlev .and. i /= j
  end function valid_pair

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
       call atom_record(unit, hname, line, ios)
       if (ios < 0) exit
       read(line,*,iostat=ios) key, n
       if (ios /= 0) call hcoll_bad(hname, 'has an invalid header')
       if (trim(key) /= 'NHCOL') cycle
       if (allocated(atom%hc_lo) .or. n < 0) call hcoll_bad(hname, 'has an invalid NHCOL count')
       allocate(atom%hc_lo(n), atom%hc_hi(n), atom%hc_form(n), &
          atom%hc_tmin(n), atom%hc_tmax(n), atom%hc_c(10,n))
       do k = 1, n
          call atom_record(unit, hname, line, ios)
          if (ios < 0) call hcoll_bad(hname, 'has truncated NHCOL rows')
          read(line,*,iostat=ios) atom%hc_lo(k), atom%hc_hi(k), &
             atom%hc_form(k), atom%hc_tmin(k), atom%hc_tmax(k), atom%hc_c(:,k)
          if (ios /= 0) call hcoll_bad(hname, 'has a short NHCOL row')
          if (.not. valid_pair(atom, atom%hc_lo(k), atom%hc_hi(k))) &
             call hcoll_bad(hname, 'has an invalid level pair')
          if (atom%e_cm1(atom%hc_hi(k)) <= atom%e_cm1(atom%hc_lo(k))) &
             call hcoll_bad(hname, 'names a level pair that is not '// &
                'lower -> upper in '//trim(fname))
          if (atom%hc_form(k) /= 1 .and. atom%hc_form(k) /= 2) &
             call hcoll_bad(hname, 'has an unknown form')
          if (.not. ieee_is_finite(atom%hc_tmin(k)) .or. &
              .not. ieee_is_finite(atom%hc_tmax(k)) .or. &
              .not. all(ieee_is_finite(atom%hc_c(:,k)))) &
             call hcoll_bad(hname, 'has a nonfinite rate parameter')
          if (atom%hc_tmin(k) <= 0.0_wp .or. atom%hc_tmax(k) < atom%hc_tmin(k)) &
             call hcoll_bad(hname, 'has an invalid temperature interval')
          if (atom%hc_form(k) == 2 .and. atom%hc_c(1,k) < 0.0_wp) &
             call hcoll_bad(hname, 'has a negative rate coefficient')
       end do
       atom%nhc = n
    end do
    close(unit)
    if (.not. allocated(atom%hc_lo)) call hcoll_bad(hname, 'is missing NHCOL')
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
    real(kind=wp),          intent(inout) :: Rh(:,:)
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
  ! Construct the source piecewise cubic once when the atom is read.
  !=========================================================================
  subroutine source_polynomial(atom, k, fname)
    type(nlevel_atom_type), intent(inout) :: atom
    integer, intent(in) :: k
    character(len=*), intent(in) :: fname
    real(kind=wp), allocatable :: x(:), y(:), slope(:), secant(:), h(:), m(:,:), rhs(:)
    real(kind=wp) :: w1, w2, fac, deriv, root1, root2, disc, z
    integer :: n, o, i, j, ip

    n = atom%nc(k)
    o = atom%offset(k)
    allocate(x(n), y(n), slope(n), secant(n-1), h(n-1))
    x = atom%node_x(o:o+n-1)
    y = atom%node_y(o:o+n-1)
    h = x(2:n)-x(1:n-1)
    secant = (y(2:n)-y(1:n-1))/h
    if (atom%interpolation(k) == CUBIC) then
       if (n == 2) then
          slope = secant(1)
       else if (n == 3) then
          slope(1) = secant(1) - h(1)*(secant(2)-secant(1))/(h(1)+h(2))
          slope(2) = secant(1) + h(1)*(secant(2)-secant(1))/(h(1)+h(2))
          slope(3) = secant(2) + h(2)*(secant(2)-secant(1))/(h(1)+h(2))
       else
          allocate(m(n,n), rhs(n))
          m = 0.0_wp
          rhs = 0.0_wp
          m(1,1) = h(2); m(1,2) = -(h(1)+h(2)); m(1,3) = h(1)
          m(n,n-2) = h(n-1); m(n,n-1) = -(h(n-2)+h(n-1)); m(n,n) = h(n-2)
          do i = 2,n-1
             m(i,i-1) = h(i-1)
             m(i,i) = 2.0_wp*(h(i-1)+h(i))
             m(i,i+1) = h(i)
             rhs(i) = 6.0_wp*(secant(i)-secant(i-1))
          end do
          do i = 1,n-1
             ip = maxloc(abs(m(i:n,i)),dim=1)+i-1
             if (m(ip,i) == 0.0_wp) call atom_bad(fname, 'singular source interpolation')
             if (ip /= i) then
                do j = i,n
                   fac=m(i,j); m(i,j)=m(ip,j); m(ip,j)=fac
                end do
                fac=rhs(i); rhs(i)=rhs(ip); rhs(ip)=fac
             end if
             do j=i+1,n
                fac=m(j,i)/m(i,i)
                m(j,i:n)=m(j,i:n)-fac*m(i,i:n)
                rhs(j)=rhs(j)-fac*rhs(i)
             end do
          end do
          if (m(n,n) == 0.0_wp) call atom_bad(fname, 'singular source interpolation')
          rhs(n)=rhs(n)/m(n,n)
          do i=n-1,1,-1
             rhs(i)=(rhs(i)-dot_product(m(i,i+1:n),rhs(i+1:n)))/m(i,i)
          end do
          slope(1)=secant(1)-h(1)*(2.0_wp*rhs(1)+rhs(2))/6.0_wp
          do i=2,n
             slope(i)=secant(i-1)+h(i-1)*(rhs(i-1)+2.0_wp*rhs(i))/6.0_wp
          end do
       end if
    else
       ! PCHIP: weighted harmonic interior slopes and shape-preserving ends.
       if (n == 2) then
          slope=secant(1)
       else
          do i=2,n-1
             if (secant(i-1)*secant(i) <= 0.0_wp) then
                slope(i)=0.0_wp
             else
                w1=2.0_wp*h(i)+h(i-1)
                w2=h(i)+2.0_wp*h(i-1)
                slope(i)=(w1+w2)/(w1/secant(i-1)+w2/secant(i))
             end if
          end do
          slope(1)=pchip_end(h(1),h(2),secant(1),secant(2))
          slope(n)=pchip_end(h(n-1),h(n-2),secant(n-1),secant(n-2))
       end if
    end if
    do i=1,n-1
       j=o+i-1
       atom%poly_a(j)=y(i)
       atom%poly_b(j)=slope(i)
       atom%poly_c(j)=(3.0_wp*secant(i)-2.0_wp*slope(i)-slope(i+1))/h(i)
       atom%poly_d(j)=(slope(i)+slope(i+1)-2.0_wp*secant(i))/h(i)**2
       if (atom%source_model(k) == BT_SOURCE .and. atom%ttype(k) /= 6 .and. &
           atom%interpolation(k) == CUBIC) then
          ! A nonnegative source cubic must have nonnegative interior extrema.
          if (atom%poly_d(j) /= 0.0_wp) then
             disc=4.0_wp*atom%poly_c(j)**2-12.0_wp*atom%poly_d(j)*atom%poly_b(j)
             if (disc >= 0.0_wp) then
                root1=(-2.0_wp*atom%poly_c(j)-sqrt(disc))/(6.0_wp*atom%poly_d(j))
                root2=(-2.0_wp*atom%poly_c(j)+sqrt(disc))/(6.0_wp*atom%poly_d(j))
                if (root1 > 0.0_wp .and. root1 < h(i)) then
                   z=atom%poly_a(j)+root1*(atom%poly_b(j)+root1*(atom%poly_c(j)+root1*atom%poly_d(j)))
                   if (z < 0.0_wp) call atom_bad(fname, 'negative cubic source extremum')
                end if
                if (root2 > 0.0_wp .and. root2 < h(i)) then
                   z=atom%poly_a(j)+root2*(atom%poly_b(j)+root2*(atom%poly_c(j)+root2*atom%poly_d(j)))
                   if (z < 0.0_wp) call atom_bad(fname, 'negative cubic source extremum')
                end if
             end if
          else if (atom%poly_c(j) /= 0.0_wp) then
             deriv=-atom%poly_b(j)/(2.0_wp*atom%poly_c(j))
             if (deriv > 0.0_wp .and. deriv < h(i)) then
                z=atom%poly_a(j)+deriv*(atom%poly_b(j)+deriv*atom%poly_c(j))
                if (z < 0.0_wp) call atom_bad(fname, 'negative cubic source extremum')
             end if
          end if
       end if
    end do
    if (.not. all(ieee_is_finite(atom%poly_a(o:o+n-2))) .or. &
        .not. all(ieee_is_finite(atom%poly_b(o:o+n-2))) .or. &
        .not. all(ieee_is_finite(atom%poly_c(o:o+n-2))) .or. &
        .not. all(ieee_is_finite(atom%poly_d(o:o+n-2)))) &
       call atom_bad(fname, 'nonfinite source interpolation')
  end subroutine source_polynomial

  real(kind=wp) function pchip_end(h1,h2,d1,d2) result(m)
    real(kind=wp), intent(in) :: h1,h2,d1,d2
    m=((2.0_wp*h1+h2)*d1-h1*d2)/(h1+h2)
    if (m*d1 <= 0.0_wp) then
       m=0.0_wp
    else if (d1*d2 < 0.0_wp .and. abs(m) > 3.0_wp*abs(d1)) then
       m=3.0_wp*d1
    end if
  end function pchip_end

  real(kind=wp) function source_value(atom,k,x) result(y)
    type(nlevel_atom_type), intent(in) :: atom
    integer, intent(in) :: k
    real(kind=wp), intent(in) :: x
    integer :: lo,hi,mid,o,n,j
    real(kind=wp) :: dx
    o=atom%offset(k); n=atom%nc(k)
    if (x <= atom%node_x(o)) then
       j=o
    else if (x >= atom%node_x(o+n-1)) then
       j=o+n-2
    else
       lo=o; hi=o+n-1
       do while (hi-lo > 1)
          mid=(lo+hi)/2
          if (x < atom%node_x(mid)) then
             hi=mid
          else
             lo=mid
          end if
       end do
       j=lo
    end if
    dx=x-atom%node_x(j)
    y=atom%poly_a(j)+dx*(atom%poly_b(j)+dx*(atom%poly_c(j)+dx*atom%poly_d(j)))
  end function source_value

  ! Upsilon(T) from either the legacy Chebyshev fit or source nodes.
  !=========================================================================
  real(kind=wp) function nlevel_upsilon(atom, k, T) result(ups)
    implicit none
    type(nlevel_atom_type), intent(in) :: atom
    integer,                intent(in) :: k
    real(kind=wp),          intent(in) :: T
    real(kind=wp) :: et, st, u, y, b0, b1, b2, T_source
    integer :: j

    if (.not. ieee_is_finite(T)) error stop 'nonfinite excitation temperature'
    if (T <= 0.0_wp) error stop 'nonpositive excitation temperature'
    if (k < 1 .or. k > atom%nups) error stop 'invalid collision row index'
    if (atom%source_model(k) == TEMPERATURE_TABLE) then
       ! Beyond the tabulated temperatures Upsilon is held at the end node
       ! (both ends); the Maxwellian rate and the Boltzmann factor take the
       ! actual T.  Valid while Upsilon is slowly varying there: forbidden
       ! and intercombination rows tend to a constant or fall slowly at
       ! high T; an allowed row grows ~ ln T, which this underestimates
       ! (the Burgess-Tully rows above keep that asymptote).  Setup prints
       ! the atoms whose tables end below te_max (nlevel_cooling_add).
       st=log(T)
       st=min(st,atom%node_x(atom%offset(k)+atom%nc(k)-1))
       if (st <= atom%node_x(atom%offset(k))) then
          ups=exp(atom%node_y(atom%offset(k)))
       else
          ups=exp(source_value(atom,k,st))
       end if
       if (.not. ieee_is_finite(ups) .or. ups <= 0.0_wp) &
          error stop 'invalid temperature table strength'
       return
    end if
    T_source = T
    if (atom%source_model(k) == BT_SOURCE) T_source = max(T, atom%bt_tmin_K)
    et = KB_OVER_RY*T_source/atom%dery(k)
    if (atom%ttype(k) == 1 .or. atom%ttype(k) == 4) then
       st = 1.0_wp - log(atom%cups(k))/log(et + atom%cups(k))
    else
       st = et/(et + atom%cups(k))
    end if
    if (atom%source_model(k) == BT_SOURCE) then
       st=min(max(st,atom%node_x(atom%offset(k))), &
          atom%node_x(atom%offset(k)+atom%nc(k)-1))
       y=source_value(atom,k,st)
    else
       st = min(max(st, atom%st_lo(k)), atom%st_hi(k))
       if (atom%st_hi(k) > atom%st_lo(k)) then
          u = 2.0_wp*(st - atom%st_lo(k))/(atom%st_hi(k) - atom%st_lo(k)) - 1.0_wp
       else
          u = 0.0_wp
       end if
       b0 = 0.0_wp;  b1 = 0.0_wp
       do j = atom%nc(k), 2, -1
          b2 = b1;  b1 = b0
          b0 = 2.0_wp*u*b1 - b2 + atom%coef(j,k)
       end do
       y = u*b0 - b1 + atom%coef(1,k)
       if (atom%logf(k) == 1) y = 10.0_wp**y
    end if
    select case (atom%ttype(k))
    case (1);  ups = y*log(et + exp(1.0_wp))
    case (2);  ups = y
    case (3);  ups = y/(et + 1.0_wp)
    case (4);  ups = y*log(et + atom%cups(k))
    case (5);  ups = y/et
    case default;  ups = 10.0_wp**y
    end select
    if (.not. ieee_is_finite(ups) .or. ups < 0.0_wp) &
       error stop 'invalid effective collision strength'
  end function nlevel_upsilon

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
    real(kind=wp),          intent(out) :: emis(:)
    real(kind=wp), intent(in), optional :: nHI
    real(kind=wp), allocatable :: Rcoll(:,:), Rrad(:,:)
    real(kind=wp), allocatable :: Rh(:,:), R(:,:), pop(:)
    real(kind=wp) :: de_erg
    integer :: k, l, u_

    if (size(emis) < atom%nrad) error stop 'short emissivity array'
    allocate(Rcoll(atom%nlev,atom%nlev), Rrad(atom%nlev,atom%nlev), &
       Rh(atom%nlev,atom%nlev), R(atom%nlev,atom%nlev), pop(atom%nlev))
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
  ! j -> i [s^-1]) with the closure sum(n_i) = 1 replacing the first row;
  ! rows 2..n are scaled to unit maximum, then LU with partial pivoting
  ! (LAPACK dgetrf/dgetrs, blocked, which the full CHIANTI atoms -- up to 322
  ! levels -- need) followed by one step of iterative refinement.
  !
  ! The refinement is required, not cosmetic. LU with partial pivoting is
  ! normwise backward stable, but the populations span tens of decades, so
  ! a residual small against the largest flux of the atom can still be large
  ! against the fluxes of a weakly populated level. Measured with dgesv
  ! alone (MOCASSIN v2 O I at 1.58e5 K, n_e = 1e2; O III at 6.31e4 K,
  ! n_e = 1e6; condition estimates 6.5e13 and 3.2e12): level balances off by
  ! 1.1e-10 and 7.6e-10 of the largest flux into or out of that level, the
  ! populations by up to 1.2e-10 / 6.5e-10. One refinement step in working
  ! precision (r = b - M x, M dx = r with the same factors, x + dx; Skeel
  ! 1980) restores componentwise backward stability: every level of the 125
  ! atoms of both models then balances to <= 1.4e-11 of its own throughput
  ! (tests/ion_stages/check_source_excitation_atom.f90, 41 temperatures x 4
  ! densities), the same as LAPACK dgerfs with its ferr estimate, which was
  ! measured to add 38% to a Fe II (142-level) table build against 19% for
  ! the single step.
  !
  ! A quantity below the smallest normal number (tiny = 2.2e-308) is zero
  ! under flush-to-zero (the production -fp-model fast=2), so each of the n
  ! eliminated terms of a population is exact only to tiny. Where the fluxes
  ! of a level are near that limit (CHIANTI C VI, O VII below ~1e4 K, largest
  ! flux ~1e-300 s^-1) its balance holds only to (n+1) tiny times its rates,
  ! the floor LAPACK's refinement applies too; such populations contribute
  ! nothing representable to the line emission.
  !=========================================================================
  subroutine level_populations(n, R, pop)
    implicit none
    integer,       intent(in)  :: n
    real(kind=wp), intent(in)  :: R(:,:)
    real(kind=wp), intent(out) :: pop(:)
    real(kind=wp), allocatable :: M(:,:), LU(:,:), rhs(:,:), x(:,:), correction(:,:)
    integer, allocatable :: ipiv(:)
    real(kind=wp) :: scale
    integer :: i, j, info
    external :: dgetrf, dgetrs, dgemv

    if (size(R,1) < n .or. size(R,2) < n .or. size(pop) < n) &
       error stop 'short population matrix or vector'
    if (.not. all(ieee_is_finite(R(1:n,1:n)))) error stop 'nonfinite transition rate'
    allocate(M(n,n), LU(n,n), rhs(n,1), x(n,1), correction(n,1), ipiv(n))
    M = 0.0_wp
    do j = 1, n
       do i = 1, n
          if (i /= j) then
             M(i,j) = R(i,j)              ! rate j -> i
             M(j,j) = M(j,j) - R(i,j)     ! loss from j
          end if
       end do
    end do
    do i=2,n
       scale=maxval(abs(M(i,1:n)))
       if (scale == 0.0_wp .or. .not. ieee_is_finite(scale)) &
          error stop 'singular statistical equilibrium row'
       M(i,1:n)=M(i,1:n)/scale
    end do
    M(1,1:n) = 1.0_wp
    rhs = 0.0_wp
    rhs(1,1) = 1.0_wp
    LU = M
    call dgetrf(n, n, LU, n, ipiv, info)
    if (info /= 0) error stop 'singular statistical equilibrium pivot'
    x = rhs
    call dgetrs('N', n, 1, LU, n, ipiv, x, n, info)
    if (info /= 0) error stop 'statistical equilibrium substitution failed'
    correction = rhs
    call dgemv('N', n, n, -1.0_wp, M, n, x(:,1), 1, 1.0_wp, correction(:,1), 1)
    call dgetrs('N', n, 1, LU, n, ipiv, correction, n, info)
    if (info /= 0) error stop 'statistical equilibrium refinement failed'
    x = x + correction
    pop = 0.0_wp
    pop(1:n) = x(1:n,1)
    if (.not. all(ieee_is_finite(pop(1:n)))) error stop 'nonfinite level population'
    if (minval(pop(1:n)) < -1.0e-12_wp) error stop 'negative level population'
    if (abs(sum(pop(1:n))-1.0_wp) > 1.0e-9_wp) &
       error stop 'unnormalized level populations'
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
    real(kind=wp),          intent(out) :: Rcoll(:,:)
    real(kind=wp),          intent(out) :: Rrad(:,:)
    real(kind=wp) :: q_ul, q_lu, de_erg, ups
    integer :: k, l, u_

    if (size(Rcoll,1) < atom%nlev .or. size(Rcoll,2) < atom%nlev .or. &
        size(Rrad,1) < atom%nlev .or. size(Rrad,2) < atom%nlev) &
       error stop 'short rate matrix'
    if (.not. ieee_is_finite(T)) error stop 'nonfinite excitation temperature'
    if (T <= 0.0_wp) error stop 'nonpositive excitation temperature'
    Rcoll = 0.0_wp
    Rrad  = 0.0_wp
    do k = 1, atom%nups
       call emitting_level(atom, atom%ul_l(k), atom%ul_u(k), l, u_, de_erg)
       ups  = nlevel_upsilon(atom, k, T)
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
    real(kind=wp),          intent(in) :: Rcoll(:,:)
    real(kind=wp),          intent(in) :: Rrad(:,:)
    real(kind=wp),          intent(in) :: ne
    real(kind=wp), allocatable :: R(:,:), pop(:)
    real(kind=wp) :: de_erg
    integer :: k, lo, hi

    if (.not. ieee_is_finite(ne)) error stop 'nonfinite electron density'
    if (ne <= 0.0_wp) error stop 'nonpositive electron density'
    if (size(Rcoll,1) < atom%nlev .or. size(Rcoll,2) < atom%nlev .or. &
        size(Rrad,1) < atom%nlev .or. size(Rrad,2) < atom%nlev) &
       error stop 'short rate matrix'
    allocate(R(atom%nlev,atom%nlev),pop(atom%nlev))
    R = ne*Rcoll(1:atom%nlev,1:atom%nlev) + Rrad(1:atom%nlev,1:atom%nlev)
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
    real(kind=wp),          intent(in)  :: Rcoll(:,:)
    real(kind=wp),          intent(in)  :: Rrad(:,:)
    real(kind=wp),          intent(in)  :: ne
    real(kind=wp),          intent(out) :: pop(:)
    real(kind=wp), allocatable :: R(:,:)
    if (size(Rcoll,1) < atom%nlev .or. size(Rcoll,2) < atom%nlev .or. &
        size(Rrad,1) < atom%nlev .or. size(Rrad,2) < atom%nlev) &
       error stop 'short rate matrix'
    allocate(R(atom%nlev,atom%nlev))
    R = ne*Rcoll(1:atom%nlev,1:atom%nlev) + Rrad(1:atom%nlev,1:atom%nlev)
    call level_populations(atom%nlev, R, pop)
  end subroutine nlevel_populations

  !=========================================================================
  real(kind=wp) function nlevel_cooling_at(atom, T, ne) result(lam)
    implicit none
    type(nlevel_atom_type), intent(in) :: atom
    real(kind=wp),          intent(in) :: T, ne
    real(kind=wp), allocatable :: Rcoll(:,:), Rrad(:,:)
    allocate(Rcoll(atom%nlev,atom%nlev),Rrad(atom%nlev,atom%nlev))
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
    real(kind=wp), allocatable :: Rcoll(:,:), Rrad(:,:)
    real(kind=wp), allocatable :: Rh(:,:), R(:,:)
    real(kind=wp), allocatable :: pe(:), peh(:)
    real(kind=wp) :: de_erg
    integer :: k, lo, hi

    dP = 0.0_wp
    if (atom%nhc == 0 .or. nHI <= 0.0_wp) return
    allocate(Rcoll(atom%nlev,atom%nlev), Rrad(atom%nlev,atom%nlev), &
       Rh(atom%nlev,atom%nlev), R(atom%nlev,atom%nlev), &
       pe(atom%nlev), peh(atom%nlev))
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
