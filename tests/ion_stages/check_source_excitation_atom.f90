program check_source_excitation_atom
! Direct source strengths (mode 1), or (mode 0) the level populations: their
! statistical equilibrium, sum(R(i,j) n_j) = n_i sum(R(j,i)) with
! R = n_e q + A for every level i, and the emitted line power
! sum n_u A_ul dE_ul they give, against both public consumers
! (nlevel_cooling_solve, nlevel_emissivities).
! Balance criterion, level by level:
!   |gains - losses| <= 1e-10 max(gains, losses)
!                       + (n+1) tiny (sum_j R(i,j) + sum_j R(j,i)).
! The first term is relative to the level's own throughput, not to the
! largest flux of the atom: populations span tens of decades, and a weakly
! populated level must balance as well as the ground level. The second is
! the underflow floor: under flush-to-zero (the production -fp-model fast=2)
! every quantity below tiny = 2.2e-308 is zero, so each of the n eliminated
! terms of a population is exact only to tiny, and a level whose fluxes are
! near that limit (CHIANTI C VI and O VII below ~1e4 K, largest flux
! ~1e-300 s^-1; level populations ~tiny at 100-800 K in Fe II, O III) cannot
! balance more finely. It is the floor LAPACK's refinement (dgerfs, SAFE1 =
! (n+1) safmin) itself applies to a row of the scaled system. The reported
! measure is |gains - losses|/(max(gains, losses) + 1e10 x the second term),
! which the criterion bounds by 1e-10.
  use mpi
  use define, only: wp
  use nlevel_mod, only: nlevel_atom_type, nlevel_load, nlevel_upsilon, &
       nlevel_cooling_prepT, nlevel_populations, nlevel_cooling_solve, nlevel_emissivities
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(nlevel_atom_type) :: atom
  real(kind=wp), allocatable :: temperature(:), rc(:,:), rr(:,:), pop(:), emis(:), rate(:,:)
  real(kind=wp), parameter :: density(4)=[1e-6_wp,100.0_wp,1e6_wp,1e10_wp]
  real(kind=wp), parameter :: hc=1.9864458571489287e-16_wp
  real(kind=wp) :: t, ne, lam, actual, gap, residual, worst_balance, worst_t, worst_ne, worst_flux
  real(kind=wp) :: gain, loss, throughput, quantum, measure
  character(len=1024) :: atom_file, grid_file, output
  integer :: i, j, k, count, unit, ierr, lo, hi, mode, worst_level
  logical :: ok
  call MPI_INIT(ierr)
  call get_command_argument(1,atom_file)
  call get_command_argument(2,grid_file)
  call get_command_argument(3,output)
  call nlevel_load(trim(atom_file),atom,ok)
  if (.not.ok) error stop 'Atom did not load'
  open(newunit=unit,file=trim(grid_file),status='old')
  read(unit,*) mode,count
  allocate(temperature(count))
  read(unit,*) temperature
  close(unit)
  open(newunit=unit,file=trim(output),status='new')
  if (mode == 1) then
     do i=1,count
        t=temperature(i)
        write(unit,'(*(es26.16e3,1x))') t,(nlevel_upsilon(atom,k,t),k=1,atom%nups)
     end do
  else
     allocate(rc(atom%nlev,atom%nlev),rr(atom%nlev,atom%nlev), &
              pop(atom%nlev),emis(atom%nrad),rate(atom%nlev,atom%nlev))
     worst_balance=0.0_wp; worst_t=0.0_wp; worst_ne=0.0_wp; worst_level=0; worst_flux=0.0_wp
     do i=1,count
        t=temperature(i)
        call nlevel_cooling_prepT(atom,t,rc,rr)
        do j=1,4
           ne=density(j)
           call nlevel_populations(atom,rc,rr,ne,pop)
           if (any(.not.ieee_is_finite(pop)).or.minval(pop)<-1e-13_wp &
               .or.abs(sum(pop)-1.0_wp)>1e-12_wp) error stop 'Invalid populations'
           rate=ne*rc+rr
           do k=1,atom%nlev
              rate(k,k)=0.0_wp
           end do
           do k=1,atom%nlev
              gain=dot_product(rate(k,:),pop)
              loss=pop(k)*sum(rate(:,k))
              residual=abs(gain-loss)
              throughput=max(gain,loss)
              quantum=real(atom%nlev+1,wp)*tiny(1.0_wp)*(sum(rate(k,:))+sum(rate(:,k)))
              measure=residual/max(throughput+1.0e10_wp*quantum,tiny(1.0_wp))
              if (measure>worst_balance) then
                 worst_balance=measure
                 worst_t=t; worst_ne=ne; worst_level=k; worst_flux=throughput
              end if
           end do
           lam=0.0_wp
           do k=1,atom%nrad
              lo=atom%rl(k);hi=atom%ru(k)
              gap=(atom%e_cm1(hi)-atom%e_cm1(lo))*hc
              lam=lam+pop(hi)*atom%A(k)*gap/ne
           end do
           ! Exercise both downstream public consumers across the profile.
           if (i==1.or.i==count.or.mod(i-1,80)==0) then
              actual=nlevel_cooling_solve(atom,rc,rr,ne)
              if (abs(actual-lam)>1e-12_wp*max(abs(lam),tiny(lam))) &
                   error stop 'Cooling consumer disagrees with populations'
              call nlevel_emissivities(atom,t,ne,emis)
              if (abs(sum(emis)/ne-lam)>1e-12_wp*max(abs(lam),tiny(lam))) &
                   error stop 'Line consumer disagrees with populations'
           end if
           if (.not.ieee_is_finite(lam).or.lam<0) error stop 'Invalid cooling'
           write(unit,'(*(es26.16e3,1x))') t,ne,lam,pop
        end do
     end do
  end if
  close(unit)
  if (mode /= 1) then
     write(*,'(a,es12.4,a,es12.4,a,es12.4,a,i0,a,es12.4e3)') 'BALANCE worst relative residual ', &
          worst_balance, ' at T ', worst_t, ' n_e ', worst_ne, ' level ', worst_level, &
          ' level throughput [s^-1] ', worst_flux
     if (worst_balance>1e-10_wp) error stop 'Populations not stationary'
  end if
  call MPI_FINALIZE(ierr)
end program check_source_excitation_atom
