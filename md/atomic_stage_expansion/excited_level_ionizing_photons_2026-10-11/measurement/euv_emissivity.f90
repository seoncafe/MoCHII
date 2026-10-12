! Collisionally excited line emission above the hydrogen edge for one atom,
! on a (T, n_e) grid, from the MoCHII statistical-equilibrium solve.
! Output rows: log10 T, log10 n_e, eps_ion [erg cm^3 s^-1] = sum over lines
! with h nu >= 13.598 eV of pop_u A h nu / n_e, lambda_all (all lines, same norm),
! and the same EUV sum split into 13.598-24.587, 24.587-54.418, >54.418 eV.
program euv_emissivity
  use define, only: wp
  use nlevel_mod
  implicit none
  type(nlevel_atom_type) :: atom
  real(kind=wp), allocatable :: rc(:,:), rr(:,:), pop(:)
  real(kind=wp) :: t, ne, e_ev, p, eps, all, b(3)
  real(kind=wp), parameter :: HC = 1.9864458571489287e-16_wp, EV = 1.602176634e-12_wp
  character(len=1024) :: atom_file
  integer :: it, ie, k, hi
  logical :: ok
  call get_command_argument(1, atom_file)
  call nlevel_load(trim(atom_file), atom, ok)
  if (.not. ok) error stop 'Atom did not load'
  allocate(rc(atom%nlev,atom%nlev), rr(atom%nlev,atom%nlev), pop(atom%nlev))
  do it = 0, 64
     t = 10.0_wp**(3.4_wp + 0.025_wp*it)
     call nlevel_cooling_prepT(atom, t, rc, rr)
     do ie = 0, 4
        ne = 10.0_wp**ie
        call nlevel_populations(atom, rc, rr, ne, pop)
        eps = 0; all = 0; b = 0
        do k = 1, atom%nrad
           e_ev = abs(atom%e_cm1(atom%ru(k)) - atom%e_cm1(atom%rl(k)))*HC/EV
           hi = atom%ru(k)
           if (atom%e_cm1(atom%rl(k)) > atom%e_cm1(hi)) hi = atom%rl(k)
           p = pop(hi)*atom%A(k)*e_ev*EV/ne
           all = all + p
           if (e_ev >= 13.598_wp) then
              eps = eps + p
              if (e_ev < 24.587_wp) then; b(1)=b(1)+p
              else if (e_ev < 54.418_wp) then; b(2)=b(2)+p
              else; b(3)=b(3)+p; end if
           end if
        end do
        write(*,'(2f8.4,6es16.6e3)') log10(t), log10(ne), eps, all, b
     end do
  end do
end program euv_emissivity
