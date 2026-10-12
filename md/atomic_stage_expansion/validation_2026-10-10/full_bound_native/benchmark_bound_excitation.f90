program benchmark_bound_excitation
use define, only: wp
use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_cooling_prepT, nlevel_cooling_solve
implicit none
type(nlevel_atom_type) :: atom
real(wp), allocatable :: rc(:,:), rr(:,:)
real(wp) :: start, finish, cooling
character(1024) :: path
integer :: i
logical :: ok
call get_command_argument(1,path)
call nlevel_load(trim(path),atom,ok)
if (.not.ok) error stop 'cannot load atom'
allocate(rc(atom%nlev,atom%nlev),rr(atom%nlev,atom%nlev))
call nlevel_cooling_prepT(atom,1e4_wp,rc,rr)
call cpu_time(start)
do i=1,100
 cooling=nlevel_cooling_solve(atom,rc,rr,1e4_wp)
end do
call cpu_time(finish)
write(*,'(i0,1x,2(es26.16e3,1x))') atom%nlev,cooling,finish-start
end program benchmark_bound_excitation
