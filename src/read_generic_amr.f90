! MoCHII: copied from MoCafe_v2.00/src/read_generic_amr.f90 (2026-07-11)
module read_generic_amr_mod
!---------------------------------------------------------------------------
! Reads generic AMR input files (FITS, HDF5, or text) and returns flat arrays
! of leaf-cell data.  Dust-only: only the geometry + dust-relevant fields are
! read.  Mandatory leaf fields are read by COLUMN NAME (with an index fallback
! for legacy 9-column files):
!
!   x, y, z      leaf center coordinates [code units]
!   level        AMR level
!   nH           total hydrogen number density [cm^-3]   (alias: dens, gasDen)
!
! Optional leaf fields (read by name when present):
!   metallicity  mass fraction Z
!   xHI          neutral hydrogen fraction
!   ndust        dust pseudo-number density [cm^-3]
!
! Any T / vx / vy / vz columns in the file are IGNORED (velocity and
! temperature are out of scope; see AMR_CLUMPS_PLAN.md Part A).  RAMSES /
! Illustris snapshots are converted to this generic format by the Python
! converters, not read here.  Reads on every rank (identical read-only data);
! amr_build_tree then shares the octree across a node via MPI-3 shared memory.
!---------------------------------------------------------------------------
  use define
  use iofile_mod
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  private

  public :: generic_amr_read

contains

  subroutine generic_amr_read(filename, xleaf, yleaf, zleaf, leaf_level, &
      nH_cgs, nleaf, boxlen_phys, metallicity, xHI, ndust, &
      origin_x, origin_y, origin_z)
    character(len=*), intent(in)       :: filename
    real(wp), allocatable, intent(out) :: xleaf(:), yleaf(:), zleaf(:)
    integer,  allocatable, intent(out) :: leaf_level(:)
    real(wp), allocatable, intent(out) :: nH_cgs(:)
    integer,               intent(out) :: nleaf
    real(wp),              intent(out) :: boxlen_phys
    real(wp), allocatable, intent(out), optional :: metallicity(:), xHI(:), ndust(:)
    real(wp),              intent(out), optional :: origin_x, origin_y, origin_z

    real(wp) :: x, y, z, nH, t6, t7, t8, t9, bl
    integer  :: lv, n, unit, ios
    character(len=:), allocatable :: resolved

    resolved = io_resolve_filename(filename)

    if (is_binary_amr_file(resolved)) then
      call generic_amr_read_binary(resolved, xleaf, yleaf, zleaf, leaf_level, &
          nH_cgs, nleaf, boxlen_phys, metallicity, xHI, ndust, &
          origin_x, origin_y, origin_z)
      return
    end if

    ! --- Text format: header "nleaf boxlen", rows x y z level nH [T vx vy vz ...] ---
    open(newunit=unit, file=resolved, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(6,'(a,a)') 'generic_amr_read: cannot open file: ', resolved
      stop 'generic_amr_read: cannot open file'
    end if
    read(unit, *) n, bl
    nleaf       = n
    boxlen_phys = bl
    allocate(xleaf(n), yleaf(n), zleaf(n), leaf_level(n), nH_cgs(n))
    do n = 1, nleaf
      ! tolerate a 5-column (dust) or legacy 9-column row; extra cols ignored.
      read(unit, *, iostat=ios) x, y, z, lv, nH, t6, t7, t8, t9
      if (ios /= 0) then
        read(unit, *, iostat=ios) x, y, z, lv, nH
        if (ios /= 0) exit
      end if
      xleaf(n) = x;  yleaf(n) = y;  zleaf(n) = z
      leaf_level(n) = lv;  nH_cgs(n) = nH
    end do
    close(unit)
    if (present(origin_x)) origin_x = -0.5_wp * boxlen_phys
    if (present(origin_y)) origin_y = -0.5_wp * boxlen_phys
    if (present(origin_z)) origin_z = -0.5_wp * boxlen_phys
  end subroutine generic_amr_read

  subroutine generic_amr_read_binary(filename, xleaf, yleaf, zleaf, leaf_level, &
      nH_cgs, nleaf, boxlen_phys, metallicity, xHI, ndust, &
      origin_x_out, origin_y_out, origin_z_out)
    character(len=*), intent(in)       :: filename
    real(wp), allocatable, intent(out) :: xleaf(:), yleaf(:), zleaf(:)
    integer,  allocatable, intent(out) :: leaf_level(:)
    real(wp), allocatable, intent(out) :: nH_cgs(:)
    integer,               intent(out) :: nleaf
    real(wp),              intent(out) :: boxlen_phys
    real(wp), allocatable, intent(out), optional :: metallicity(:), xHI(:), ndust(:)
    real(wp),              intent(out), optional :: origin_x_out, origin_y_out, origin_z_out

    type(io_file_type) :: iofh
    integer        :: status
    integer(int32) :: nleaf_i4
    real(wp) :: origin_x, origin_y, origin_z
    integer  :: status_ox, status_oy, status_oz

    status = 0
    call io_open_old(iofh, trim(filename), status)
    if (status /= 0) stop 'generic_amr_read_binary: cannot open file'
    call io_move_to_next_section(iofh, status)
    if (status /= 0) stop 'generic_amr_read_binary: cannot move to binary table HDU'

    call io_get_keyword(iofh, 'NAXIS2', nleaf_i4, status)
    if (status /= 0) stop 'generic_amr_read_binary: cannot read NAXIS2'
    boxlen_phys = 0.0_wp
    call io_get_keyword(iofh, 'BOXLEN', boxlen_phys, status)
    if (status /= 0) stop 'generic_amr_read_binary: cannot read BOXLEN'

    origin_x = -0.5_wp * boxlen_phys;  origin_y = origin_x;  origin_z = origin_x
    status_ox = 0;  status_oy = 0;  status_oz = 0
    call io_get_keyword(iofh, 'ORIGINX', origin_x, status_ox)
    call io_get_keyword(iofh, 'ORIGINY', origin_y, status_oy)
    call io_get_keyword(iofh, 'ORIGINZ', origin_z, status_oz)
    if (status_ox /= 0) origin_x = -0.5_wp * boxlen_phys
    if (status_oy /= 0) origin_y = -0.5_wp * boxlen_phys
    if (status_oz /= 0) origin_z = -0.5_wp * boxlen_phys

    nleaf = int(nleaf_i4)
    allocate(xleaf(nleaf), yleaf(nleaf), zleaf(nleaf), leaf_level(nleaf), nH_cgs(nleaf))

    !--- mandatory columns by name (with index fallback for legacy files)
    call read_real_col(iofh, 'x',            1, nleaf, xleaf)
    call read_real_col(iofh, 'y',            2, nleaf, yleaf)
    call read_real_col(iofh, 'z',            3, nleaf, zleaf)
    call read_int_col (iofh, 'level',        4, nleaf, leaf_level)
    call read_real_col(iofh, 'nH,dens,gasDen', 5, nleaf, nH_cgs)

    if (present(origin_x_out)) origin_x_out = origin_x
    if (present(origin_y_out)) origin_y_out = origin_y
    if (present(origin_z_out)) origin_z_out = origin_z

    !--- optional columns by name
    call try_read_optional_column(iofh, nleaf, 'metallicity', metallicity)
    call try_read_optional_column(iofh, nleaf, 'xHI',         xHI)
    call try_read_optional_column(iofh, nleaf, 'ndust',       ndust)

    call io_close(iofh, status)
  end subroutine generic_amr_read_binary

  !--- Read a mandatory real column by a comma-separated list of candidate
  !    names; fall back to the given 1-based column index; abort if neither.
  subroutine read_real_col(iofh, names, idx_fallback, nleaf, arr)
    type(io_file_type), intent(inout) :: iofh
    character(len=*),   intent(in)    :: names
    integer,            intent(in)    :: idx_fallback, nleaf
    real(wp),           intent(out)   :: arr(nleaf)
    integer :: colnum, status, ks, ke, lstr
    character(len=32) :: cand
    logical :: found
    found = .false.
    lstr = len_trim(names);  ks = 1
    do while (ks <= lstr .and. .not. found)
      ke = index(names(ks:lstr), ',')
      if (ke == 0) then
        cand = adjustl(names(ks:lstr));  ks = lstr + 1
      else
        cand = adjustl(names(ks:ks+ke-2));  ks = ks + ke
      end if
      if (len_trim(cand) == 0) cycle
      status = 0
      call io_get_column_number(iofh, trim(cand), colnum, status)
      if (status == 0) then
        call io_read_table_column(iofh, colnum, arr, status)
        if (status == 0) found = .true.
      end if
    end do
    if (.not. found) then
      status = 0
      call io_read_table_column(iofh, idx_fallback, arr, status)
      if (status /= 0) then
        write(6,'(3a,i0)') 'generic_amr_read: missing mandatory column {', &
             trim(names), '} and index fallback ', idx_fallback
        stop 'generic_amr_read: missing mandatory column'
      end if
    end if
  end subroutine read_real_col

  subroutine read_int_col(iofh, name, idx_fallback, nleaf, arr)
    type(io_file_type), intent(inout) :: iofh
    character(len=*),   intent(in)    :: name
    integer,            intent(in)    :: idx_fallback, nleaf
    integer,            intent(out)   :: arr(nleaf)
    integer(int32), allocatable :: tmp(:)
    integer :: colnum, status
    allocate(tmp(nleaf))
    status = 0
    call io_get_column_number(iofh, trim(name), colnum, status)
    if (status /= 0) colnum = idx_fallback
    status = 0
    call io_read_table_column(iofh, colnum, tmp, status)
    if (status /= 0) then
      write(6,'(3a)') 'generic_amr_read: cannot read mandatory column ', trim(name), '.'
      stop 'generic_amr_read: missing level column'
    end if
    arr(:) = int(tmp(:))
    deallocate(tmp)
  end subroutine read_int_col

  subroutine try_read_optional_column(iofh, nleaf, colname, arr)
    type(io_file_type), intent(inout) :: iofh
    integer,            intent(in)    :: nleaf
    character(len=*),   intent(in)    :: colname
    real(wp), allocatable, intent(out), optional :: arr(:)
    integer :: colnum, status
    if (.not. present(arr)) return
    status = 0
    call io_get_column_number(iofh, colname, colnum, status)
    if (status /= 0 .or. colnum <= 0) return
    allocate(arr(nleaf))
    status = 0
    call io_read_table_column(iofh, colnum, arr, status)
    if (status /= 0) deallocate(arr)
  end subroutine try_read_optional_column

  logical function is_binary_amr_file(filename)
    character(len=*), intent(in) :: filename
    character(len=:), allocatable :: name
    integer :: n
    name = trim(filename);  n = len(name)
    is_binary_amr_file = .false.
    if (n >= 5) then
      if (name(n-4:n) == '.fits' .or. name(n-4:n) == '.FITS') is_binary_amr_file = .true.
      if (name(n-4:n) == '.hdf5' .or. name(n-4:n) == '.HDF5') is_binary_amr_file = .true.
    end if
    if (n >= 8) then
      if (name(n-7:n) == '.fits.gz' .or. name(n-7:n) == '.FITS.GZ') is_binary_amr_file = .true.
    end if
    if (n >= 3) then
      if (name(n-2:n) == '.h5' .or. name(n-2:n) == '.H5') is_binary_amr_file = .true.
    end if
  end function is_binary_amr_file

end module read_generic_amr_mod
