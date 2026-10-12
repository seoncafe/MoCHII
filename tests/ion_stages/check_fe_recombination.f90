program check_fe_recombination
! Recombination coefficients of the Fe stages read back from the stationary
! stage balance: with photoionization (gamma) and recombination alone --
! electron-impact ionization switched off and no hydrogen, so no charge
! transfer -- n(q+1)/n(q) = gamma/(n_e alpha_q) for each adjacent pair, so
! alpha_q = gamma f(q)/(n_e f(q+1)).  A coefficient measurement on the
! production Fe file, not an Fe plasma model.
  use define
  use species_mod
  use mpi
  implicit none
  integer :: ierr, j, q
  real(kind=wp) :: t, alpha, nem, direct_ne, cxi, cxr
  real(kind=wp) :: frac(MAX_ION_STAGES), gamma(MAX_ION_STAGES-1,1)
  character(len=PATH_LEN) :: directory, output
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD; mpar%h_rank=0
  call get_command_argument(1,directory)
  call get_command_argument(2,output)
  par%atomic_dir=trim(directory); par%cooling_model='low_density'
  par%metal_collisional_ionization=.false.
  par%abund_C=0; par%abund_N=0; par%abund_O=0; par%abund_Ne=0
  par%abund_S=0; par%abund_Ar=0; par%abund_Mg=0; par%abund_Si=0
  par%abund_Cl=0; par%abund_Ca=0; par%abund_Fe=1e-4_wp
  call species_setup(1)
  if(n_elements/=1.or.elem_nstage(1)/=7) error stop 'Fe rate fixture roster'
  gamma=1e-8_wp
  call species_set_photoionization_rates(1,gamma)
  open(unit=71,file=trim(output),status='new')
  do q=1,6
     do j=0,300
        t=real(q*q,wp)*20.0_wp*(5e5_wp)**(real(j,wp)/300.0_wp)
        call species_fractions(1,1,t,10.0_wp,0.0_wp,0.0_wp,frac)
        if(any(frac<=0).or.abs(sum(frac)-1.0_wp)>2e-14_wp) error stop 'Fe fractions'
        alpha=gamma(q,1)*frac(q)/(10.0_wp*frac(q+1))
        call species_ne_prepare(1,t)
        call species_metal_sums(100.0_wp,10.0_wp,0.0_wp,0.0_wp,nem,cxi,cxr)
        direct_ne=0
        block
          integer :: k
          do k=1,7
             direct_ne=direct_ne+100.0_wp*elem_abund(1)*real(k-1,wp)*frac(k)
          end do
        end block
        if(abs(nem-direct_ne)>2e-14_wp) error stop 'Fe direct/cache charge'
        write(71,'(i4,2es27.17e3)') q,t,alpha
     end do
  end do
  close(71)
  print *, 'PASS: 1806 RR+DR values and direct/cache charge'
  call MPI_FINALIZE(ierr)
end program check_fe_recombination
