program check_bound_excitation_atom
! Stationary electron excitation and spontaneous emission; no photoexcitation.
! The cooling coefficient of the level atom and its fitted low-density
! cooling file on 801 temperatures, 1e3-1e5 K, and four densities; the level
! populations are written unless a fourth argument 'cooling' asks for the
! cooling columns alone.
  use mpi
  use define, only: wp
  use cooling_mod, only: cooling_fit_type, cooling_load, cooling_eval
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, &
       nlevel_cooling_prepT, nlevel_cooling_solve, nlevel_populations
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(cooling_fit_type) :: fit
  type(nlevel_atom_type) :: atom
  real(kind=wp) :: t, ne, lam
  real(kind=wp), allocatable :: pop(:), rc(:,:), rr(:,:)
  real(kind=wp), parameter :: density(4)=[1e-6_wp,100.0_wp,1e6_wp,1e10_wp]
  integer :: i, j, ierr, unit
  logical :: ok
  character(len=1024) :: cooling_file, atom_file, output, columns
  call MPI_INIT(ierr)
  call get_command_argument(1,cooling_file)
  call get_command_argument(2,atom_file)
  call get_command_argument(3,output)
  call get_command_argument(4,columns)
  call cooling_load(trim(cooling_file),fit)
  call nlevel_load(trim(atom_file),atom,ok)
  if (.not.ok) error stop 'Excitation atom failed to load'
  allocate(pop(atom%nlev), rc(atom%nlev,atom%nlev), rr(atom%nlev,atom%nlev))
  open(newunit=unit,file=trim(output),status='new')
  do i=0,800
     t=10.0_wp**(3.0_wp+2.0_wp*real(i,wp)/800.0_wp)
     call nlevel_cooling_prepT(atom,t,rc,rr)
     do j=1,4
        ne=density(j)
        lam=nlevel_cooling_solve(atom,rc,rr,ne)
        call nlevel_populations(atom,rc,rr,ne,pop)
        if (.not.ieee_is_finite(lam).or.lam<0.or.any(.not.ieee_is_finite(pop))) &
             error stop 'Invalid excitation cooling or populations'
        if (minval(pop(1:atom%nlev))<0.or.abs(sum(pop(1:atom%nlev))-1.0_wp)>2e-13_wp) &
             error stop 'Excitation population normalization'
        if (trim(columns)=='cooling') then
           write(unit,'(4(es26.16e3,1x))') t,ne,cooling_eval(fit,t),lam
        else
           write(unit,'(*(es26.16e3,1x))') t,ne,cooling_eval(fit,t),lam,pop(1:atom%nlev)
        end if
     end do
  end do
  close(unit)
  call MPI_FINALIZE(ierr)
end program check_bound_excitation_atom
