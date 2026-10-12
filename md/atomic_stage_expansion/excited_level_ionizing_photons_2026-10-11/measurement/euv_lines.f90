! Number of EUV lines (h nu >= 13.598 eV) needed to carry 1 - 1e-4 of the
! collisional EUV emissivity at each node of a coarse (T, n_e) grid; prints
! the largest count over the nodes and the size of the union over nodes.
program euv_lines
  use define, only: wp
  use nlevel_mod
  implicit none
  type(nlevel_atom_type) :: atom
  real(kind=wp), allocatable :: rc(:,:), rr(:,:), pop(:), em(:)
  integer, allocatable :: idx(:)
  logical, allocatable :: used(:)
  real(kind=wp) :: t, ne, e_ev, tot, acc
  real(kind=wp), parameter :: HC = 1.9864458571489287e-16_wp, EV = 1.602176634e-12_wp
  character(len=1024) :: atom_file
  integer :: it, ie, k, hi, m, nmax, neuv
  logical :: ok
  call get_command_argument(1, atom_file)
  call nlevel_load(trim(atom_file), atom, ok)
  allocate(rc(atom%nlev,atom%nlev), rr(atom%nlev,atom%nlev), pop(atom%nlev), em(atom%nrad), idx(atom%nrad), used(atom%nrad))
  used = .false.; nmax = 0; neuv = 0
  do k = 1, atom%nrad
     if (abs(atom%e_cm1(atom%ru(k)) - atom%e_cm1(atom%rl(k)))*HC/EV >= 13.598_wp) neuv = neuv + 1
  end do
  do it = 0, 15
     t = 10.0_wp**(3.5_wp + 0.1_wp*it)
     call nlevel_cooling_prepT(atom, t, rc, rr)
     do ie = 0, 4
        ne = 10.0_wp**ie
        call nlevel_populations(atom, rc, rr, ne, pop)
        em = 0
        do k = 1, atom%nrad
           e_ev = abs(atom%e_cm1(atom%ru(k)) - atom%e_cm1(atom%rl(k)))*HC/EV
           hi = atom%ru(k); if (atom%e_cm1(atom%rl(k)) > atom%e_cm1(hi)) hi = atom%rl(k)
           if (e_ev >= 13.598_wp) em(k) = max(pop(hi),0.0_wp)*atom%A(k)*e_ev
        end do
        tot = sum(em); if (tot <= 0) cycle
        call sort_desc(em, idx, atom%nrad)
        acc = 0; m = 0
        do while (acc < (1.0_wp-1.0e-4_wp)*tot .and. m < atom%nrad)
           m = m + 1; acc = acc + em(idx(m)); used(idx(m)) = .true.
        end do
        nmax = max(nmax, m)
     end do
  end do
  write(*,'(a,4i8)') trim(atom_file(index(atom_file,'/',back=.true.)+1:)), neuv, nmax, count(used), atom%nlev
contains
  subroutine sort_desc(a, ix, n)
    real(kind=wp), intent(in) :: a(:)
    integer, intent(out) :: ix(:)
    integer, intent(in) :: n
    integer :: i, j, tmp
    ix = [(i, i=1,n)]
    ! simple shell sort descending
    do i = 2, n
       tmp = ix(i); j = i - 1
       do while (j >= 1)
          if (a(ix(j)) >= a(tmp)) exit
          ix(j+1) = ix(j); j = j - 1
       end do
       ix(j+1) = tmp
    end do
  end subroutine sort_desc
end program euv_lines
