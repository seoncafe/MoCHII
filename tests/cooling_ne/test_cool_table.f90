program test_cool_table
!---------------------------------------------------------------------------
! Standalone gate for the (T, n_e) line-cooling table (nlevel_cooling_mod):
!   G-a  interpolation: table Lambda_SE (nlevel_cooling_lam_se) vs the direct
!        n-level solve (nlevel_cooling_at) at 200 random off-grid (T, n_e)
!        points for each ion.  Criterion: worst <= 0.5%, median <= 0.1%.
!   G-c  low-density identity: nlevel_cooling_apply(idx, 1, T, ne_floor) = 1
!        for a positive source reference; exact source zeros follow the
!        complete or partial atom rule.
!   G-e  the runtime complete/partial cooling combination against the direct
!        n-level solve and the actual low-density cooling fit.
! Links only define + nlevel_mod + nlevel_cooling_mod (no MPI/HDF5).
!---------------------------------------------------------------------------
  use define
  use nlevel_mod,         only : nlevel_atom_type, nlevel_load, nlevel_cooling_at, metal_excitation_data_file
  use nlevel_cooling_mod, only : nlevel_cooling_add, nlevel_cooling_lam_se, &
                                 nlevel_cooling_apply
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  integer, parameter :: NION = 9, NPT = 200, NTC = 40
  character(len=8) :: names(NION)
  integer :: stages(NION), idx(NION)
  type(nlevel_atom_type) :: atom
  character(len=PATH_LEN) :: fname
  logical :: ok
  integer :: k, ip, seedn, nfail, total_fail
  integer, allocatable :: seed(:)
  real(kind=wp) :: r1, r2, T, ne, lt, lne
  real(kind=wp) :: rel, relmax, relarr(NPT), relmed
  real(kind=wp) :: lam_tab, lam_dir, Tworst, neworst
  real(kind=wp) :: worst_max, worst_med, gc_worst, dev
  real(kind=wp) :: lam_fit, lam_ref, lam_exp
  real(kind=wp) :: ref_worst, fit_worst, direct_worst, table_worst, expected_worst
  real(kind=wp) :: fitA(64), fitT(64)
  integer :: nfit

  par%atomic_dir     = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/data/atomic'
  !--- the default excitation model unless a second argument names one
  if (command_argument_count() > 2) error stop 'usage: test_cool_table [ATOMIC_DIR [MODEL]]'
  if (command_argument_count() >= 1) call get_command_argument(1, par%atomic_dir)
  if (command_argument_count() == 2) call get_command_argument(2, par%metal_excitation_model)
  mpar%p_rank        = 0

  !--- C I and N I close the neutral-gas gap: their FIR/optical lines
  !--- ([C I] 609.1/370.4um, [N I] 5198/5200) have low critical densities.
  names  = [character(len=8) :: 'o', 'o', 'n', 's', 'ne', 'c', 'fe', 'c', 'n']
  stages = [3, 2, 2, 2, 3, 2, 3, 1, 1]

  write(*,'(a)') '=== register ions (nlevel_cooling_add) ==='
  total_fail = 0
  do k = 1, NION
     idx(k) = nlevel_cooling_add(trim(names(k)), stages(k))
     if (idx(k) <= 0) then
        write(*,'(a,a,i0)') ' Missing n-level atom: ', trim(names(k)), stages(k)
        total_fail = total_fail + 1
     end if
  end do

  !--- fixed seed for reproducibility.
  call random_seed(size=seedn)
  allocate(seed(seedn))
  seed = 20260725
  call random_seed(put=seed)

  write(*,'(a)') ''
  write(*,'(a)') '=== G-a: table interp vs direct solve (200 off-grid pts/ion) ==='
  write(*,'(a)') '  ion stage        max_rel       median_rel  idx nlev'// &
     '   T_worst  ne_worst'
  worst_max = 0.0_wp;  worst_med = 0.0_wp;  nfail = 0
  do k = 1, NION
     if (idx(k) <= 0) cycle
     ! Read the same selected atom used to construct this cooling table.
     fname=metal_excitation_data_file('nlevel',trim(names(k)),stages(k))
     call nlevel_load(trim(fname), atom, ok)
     if (.not. ok) error stop 'Cannot reload n-level atom'
     relmax = 0.0_wp;  Tworst = 0.0_wp;  neworst = 0.0_wp
     do ip = 1, NPT
        call random_number(r1);  call random_number(r2)
        lt  = 3.0_wp + 2.0_wp*r1           ! log10 T  in [3, 5]  (1e3..1e5 K)
        lne = -2.0_wp + 9.0_wp*r2          ! log10 ne in [-2, 7] (1e-2..1e7)
        T   = 10.0_wp**lt
        ne  = 10.0_wp**lne
        lam_tab = nlevel_cooling_lam_se(idx(k), T, ne)
        lam_dir = nlevel_cooling_at(atom, T, ne)
        if (.not. all(ieee_is_finite([lam_tab, lam_dir])) .or. &
            min(lam_tab, lam_dir) < 0.0_wp) error stop 'Invalid cooling value'
        if (lam_dir > 0.0_wp) then
           rel = abs(lam_tab/lam_dir - 1.0_wp)
        else
           rel = 0.0_wp
           if (lam_tab /= 0.0_wp) rel = huge(1.0_wp)
        end if
        relarr(ip) = rel
        if (rel > relmax) then
           relmax = rel;  Tworst = T;  neworst = ne
        end if
     end do
     relmed = median(relarr, NPT)
     write(*,'(2x,a4,i5,2es16.5,i5,i5,2es10.2)') trim(names(k)), stages(k), &
        relmax, relmed, idx(k), atom%nlev, Tworst, neworst
     worst_max = max(worst_max, relmax)
     worst_med = max(worst_med, relmed)
     if (relmax > 5.0e-3_wp .or. relmed > 1.0e-3_wp) nfail = nfail + 1
  end do
  write(*,'(a,es12.4,a,es12.4)') ' WORST across ions:  max = ', worst_max, &
     '   median = ', worst_med
  if (nfail == 0) then
     write(*,'(a)') ' G-a PASS (every ion: max <= 0.5%, median <= 0.1%)'
  else
     write(*,'(a,i0,a)') ' G-a FAIL: ', nfail, ' ion(s) exceed the criterion'
  end if
  total_fail = total_fail + nfail

  write(*,'(a)') ''
  write(*,'(a)') '=== G-c: low-density identity and exact-zero rule ==='
  gc_worst = 0.0_wp
  do k = 1, NION
     if (idx(k) <= 0) cycle
     fname=metal_excitation_data_file('nlevel',trim(names(k)),stages(k))
     call nlevel_load(trim(fname), atom, ok)
     if (.not. ok) error stop 'Cannot reload n-level atom'
     do ip = 1, NTC
        lt = 3.0_wp + 2.0_wp*real(ip-1, wp)/real(NTC-1, wp)   ! 1e3..1e5 K
        T  = 10.0_wp**lt
        lam_ref = nlevel_cooling_lam_se(idx(k), T, 1.0e-6_wp)
        lam_tab = nlevel_cooling_apply(idx(k), 1.0_wp, T, 1.0e-6_wp)
        if (.not. all(ieee_is_finite([lam_ref, lam_tab])) .or. &
            min(lam_ref, lam_tab) < 0.0_wp) error stop 'Invalid low-density cooling'
        if (lam_ref > 0.0_wp) then
           dev = abs(lam_tab - 1.0_wp)
        else if (lam_ref == 0.0_wp) then
           dev = abs(lam_tab - 1.0_wp)
           if (atom%density_ratio_cooling) dev = abs(lam_tab)
        else
           error stop 'Negative low-density reference'
        end if
        gc_worst = max(gc_worst, dev)
     end do
  end do
  write(*,'(a,es12.4)') ' worst |apply - 1| at ne_floor = ', gc_worst
  if (gc_worst <= 1.0e-9_wp) then
     write(*,'(a)') ' G-c PASS (identity to ~1e-10)'
  else
     write(*,'(a)') ' G-c FAIL'
     total_fail = total_fail + 1
  end if

  !--- G-e: compare the runtime cooling combination with the same physical
  !--- rule evaluated using an independent direct n-level solve. A converged
  !--- atom applies the direct density ratio. A partial atom retains the
  !--- low-density residual when the fit covers more cooling than the atom;
  !--- otherwise it scales the direct density-dependent cooling. No bound
  !--- Lambda <= lam_fit is assumed: density-dependent enhancement is possible.
  write(*,'(a)') ''
  write(*,'(a)') '=== G-e: complete/partial combiner vs direct solve ==='
  write(*,'(a)') '  ion stage        max_rel       median_rel'// &
     '       T@worst    ne@worst'
  worst_max = 0.0_wp;  worst_med = 0.0_wp;  nfail = 0
  do k = 1, NION
     if (idx(k) <= 0) cycle
     call load_fit(trim(names(k)), stages(k), fitA, fitT, nfit, ok)
     if (.not. ok) then
        write(*,'(2x,a4,i5,a)') trim(names(k)), stages(k), '   NO FIT'
        nfail = nfail + 1
        cycle
     end if
     fname=metal_excitation_data_file('nlevel',trim(names(k)),stages(k))
     call nlevel_load(trim(fname), atom, ok)
     if (.not. ok) error stop 'Cannot reload n-level atom'
     relmax = 0.0_wp;  Tworst = 0.0_wp;  neworst = 0.0_wp
     do ip = 1, NPT
        call random_number(r1);  call random_number(r2)
        lt  = 3.0_wp + 2.0_wp*r1
        lne = -2.0_wp + 9.0_wp*r2
        T   = 10.0_wp**lt
        ne  = 10.0_wp**lne
        lam_fit  = fiteval(fitA, fitT, nfit, T)
        lam_ref  = nlevel_cooling_at(atom, T, 1.0e-6_wp)
        lam_dir  = nlevel_cooling_at(atom, T, ne)
        if (.not. all(ieee_is_finite([lam_fit, lam_ref, lam_dir])) .or. &
            lam_fit <= 0.0_wp .or. min(lam_ref, lam_dir) < 0.0_wp) &
            error stop 'Invalid direct cooling or fit'
        if (atom%density_ratio_cooling) then
           if (lam_ref > 0.0_wp) then
              lam_exp = lam_fit*(lam_dir/lam_ref)
           else
              lam_exp = lam_dir
           end if
        else if (lam_ref <= lam_fit) then
           lam_exp = lam_dir + (lam_fit - lam_ref)
        else
           lam_exp = (lam_fit/lam_ref)*lam_dir
        end if
        lam_tab  = nlevel_cooling_apply(idx(k), lam_fit, T, ne)
        if (.not. all(ieee_is_finite([lam_exp, lam_tab])) .or. &
            min(lam_exp, lam_tab) < 0.0_wp) error stop 'Invalid combined cooling'
        if (lam_exp > 0.0_wp) then
           rel = abs(lam_tab/lam_exp-1.0_wp)
        else
           rel = 0.0_wp
           if (lam_tab /= 0.0_wp) rel = huge(1.0_wp)
        end if
        relarr(ip) = rel
        if (rel > relmax) then
           relmax = rel;  Tworst = T;  neworst = ne
           ref_worst = lam_ref;  fit_worst = lam_fit
           direct_worst = lam_dir;  table_worst = lam_tab
           expected_worst = lam_exp
        end if
     end do
     relmed = median(relarr, NPT)
     write(*,'(2x,a4,i5,2es16.5,2es11.2)') trim(names(k)), stages(k), &
        relmax, relmed, Tworst, neworst
     if (relmax > 5.0e-3_wp) then
        write(*,'(a,5es13.4)') '   ref fit direct table expected: ', &
           ref_worst, fit_worst, direct_worst, table_worst, expected_worst
     end if
     worst_max = max(worst_max, relmax)
     worst_med = max(worst_med, relmed)
     if (relmax > 5.0e-3_wp .or. relmed > 1.0e-3_wp) nfail = nfail + 1
  end do
  write(*,'(a,es12.4,a,es12.4)') ' WORST across ions:  max = ', &
     worst_max, '   median = ', worst_med
  if (nfail == 0) then
     write(*,'(a)') ' G-e PASS (every ion: max <= 0.5%, median <= 0.1%)'
  else
     write(*,'(a,i0,a)') ' G-e FAIL: ', nfail, ' ion(s) exceed the criterion'
  end if
  total_fail = total_fail + nfail
  if (total_fail /= 0) error stop 'Cooling-table gate failed'

contains

  !--- read a cooling fit file (cooling_<name>_<stage>.txt): '#' comments,
  !--- then nterm, then nterm rows of A_i T_i.  Same layout as cooling_load.
  subroutine load_fit(name, stage, A, Ti, nt, ok)
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: stage
    real(kind=wp),    intent(out) :: A(:), Ti(:)
    integer,          intent(out) :: nt
    logical,          intent(out) :: ok
    character(len=PATH_LEN) :: fn
    character(len=256) :: line
    integer :: u, ios, i
    fn = metal_excitation_data_file('cooling', trim(name), stage)
    ok = .false.;  nt = 0
    open(newunit=u, file=trim(fn), status='old', iostat=ios)
    if (ios /= 0) return
    do
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) then;  close(u);  return;  end if
       line = adjustl(line)
       if (len_trim(line) > 0 .and. line(1:1) /= '#') exit
    end do
    read(line,*,iostat=ios) nt
    if (ios /= 0) then;  close(u);  return;  end if
    if (nt < 1 .or. nt > size(A) .or. nt > size(Ti)) then
       close(u)
       return
    end if
    do i = 1, nt
       read(u,*,iostat=ios) A(i), Ti(i)
       if (ios /= 0) then;  close(u);  return;  end if
    end do
    close(u)
    ok = .true.
  end subroutine load_fit

  real(kind=wp) function fiteval(A, Ti, nt, T) result(lam)
    real(kind=wp), intent(in) :: A(:), Ti(:), T
    integer,       intent(in) :: nt
    integer :: i
    lam = 0.0_wp
    do i = 1, nt
       lam = lam + A(i)*exp(-Ti(i)/T)
    end do
    lam = lam/sqrt(T)
  end function fiteval

  real(kind=wp) function median(a, n) result(m)
    integer, intent(in) :: n
    real(kind=wp), intent(in) :: a(n)
    real(kind=wp) :: b(n), tmp
    integer :: i, j
    b = a
    do i = 1, n-1
       do j = i+1, n
          if (b(j) < b(i)) then
             tmp = b(i);  b(i) = b(j);  b(j) = tmp
          end if
       end do
    end do
    if (mod(n,2) == 0) then
       m = 0.5_wp*(b(n/2) + b(n/2+1))
    else
       m = b((n+1)/2)
    end if
  end function median

end program test_cool_table
