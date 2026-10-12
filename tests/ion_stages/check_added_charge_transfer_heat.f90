program check_added_charge_transfer_heat
! Charge conservation and kinetic-energy accounting of the metal-H charge
! exchange at fixed rates and densities (not a transported nebular model):
! the cached sums against the reaction-by-reaction products of the stage
! fractions, the electron closure, and metal_charge_exchange_heating against
! the sum of rate x kinetic energy of every reaction.
  use define
  use species_mod
  use mpi
  implicit none
  integer :: ierr, ie, it, ns, state, nstage, ntransition, expected_stages
  real(kind=wp) :: temperature, nh, ne, nhi, nhii, ki, kr, di, dr, split
  real(kind=wp) :: frac(MAX_ION_STAGES), gamma(MAX_ION_STAGES-1,1)
  real(kind=wp) :: incoming(MAX_ION_STAGES,MAX_METAL_ELEMENTS)
  real(kind=wp) :: outgoing(MAX_ION_STAGES,MAX_METAL_ELEMENTS)
  real(kind=wp) :: nem, sum_in, sum_out, direct_ne, direct_in, direct_out
  real(kind=wp) :: direct_heat, computed_heat, limited_ion_rate, nx, gross_heat
  real(kind=wp), parameter :: temperatures(3)=[5000.0_wp,10000.0_wp,30000.0_wp]
  character(len=PATH_LEN) :: directory, output, expected_text, option_file
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD; mpar%h_rank=0
  call get_command_argument(1,directory)
  call get_command_argument(2,output)
  call get_command_argument(3,expected_text)
  call get_command_argument(4,option_file)
  if (len_trim(option_file) > 0) then
     block
       use setup_mod, only: read_parameters_with_implicit_defaults
       call read_parameters_with_implicit_defaults(trim(option_file))
     end block
  end if
  read(expected_text,*) expected_stages
  par%atomic_dir=trim(directory); par%cooling_model='low_density'
  par%use_metals=.true.; par%charge_exchange='two_way'
  par%abund_C=1e-4_wp; par%abund_N=1e-4_wp; par%abund_O=1e-4_wp
  par%abund_Ne=1e-4_wp; par%abund_S=1e-4_wp; par%abund_Ar=1e-4_wp
  par%abund_Mg=1e-4_wp; par%abund_Fe=1e-4_wp; par%abund_Si=1e-4_wp
  par%abund_Cl=1e-4_wp; par%abund_Ca=1e-4_wp
  call species_setup(1)
  nstage=0; ntransition=0; gamma=2e-8_wp
  do ie=1,n_elements
     ns=elem_nstage(ie)
     nstage=nstage+ns; ntransition=ntransition+ns-1
     call species_set_photoionization_rates(ie,gamma(1:ns-1,:))
  end do
  if (nstage/=expected_stages .or. ntransition/=expected_stages-11) &
       error stop 'staged roster mismatch'
  open(unit=71,file=trim(output),status='replace')
  do state=1,3
     temperature=temperatures(state); nh=100; ne=10; nhi=50; nhii=50
     call species_ne_prepare(1,temperature)
     call species_metal_sums(nh,ne,nhi,nhii,nem,sum_in,sum_out,incoming,outgoing)
     direct_ne=0; direct_in=0; direct_out=0; direct_heat=0; gross_heat=0
     do ie=1,n_elements
        ns=elem_nstage(ie); nx=elem_abund(ie)*nh
        call species_fractions(ie,1,temperature,ne,nhi,nhii,frac)
        if (abs(sum(frac(1:ns))-1.0_wp)>1e-14_wp .or. any(frac(1:ns)<0)) &
             error stop 'invalid staged fractions'
        do it=1,ns
           direct_ne=direct_ne+nx*real(it-1,wp)*frac(it)
        end do
        do it=1,ns-1
           call species_cx_coefficients(ie,it,temperature,ki,kr,ne)
           call species_cx_kinetic_energy(ie,it,di,dr,split,temperature,ne)
           direct_in=direct_in+ki*nx*frac(it)
           direct_out=direct_out+kr*nx*frac(it+1)
           if (abs(incoming(it,ie)-ki*nx*frac(it))>1e-24_wp .or. &
               abs(outgoing(it,ie)-kr*nx*frac(it+1))>1e-24_wp) &
                error stop 'reaction-resolved charge conservation'
           limited_ion_rate=ki
           if (split>0) limited_ion_rate=min(ki,split*kr*exp(min(di*ev2erg/ &
                                               (kboltz_cgs*temperature),600.0_wp)))
           direct_heat=direct_heat+nx*(nhi*frac(it+1)*kr*dr+nhii*frac(it)*limited_ion_rate*di)
           gross_heat=gross_heat+nx*(abs(nhi*frac(it+1)*kr*dr)+abs(nhii*frac(it)*limited_ion_rate*di))
           write(71,'(a2,2i4,7es27.17e3)') elem_name(ie),it,state,temperature,ki,kr,di,dr,split,frac(it+1)
        end do
     end do
     computed_heat=metal_charge_exchange_heating(1,temperature,nh,ne,nhi,nhii)
     if (abs(nem-direct_ne)>1e-14_wp) error stop 'staged electron closure'
     if (abs(sum_in/direct_in-1)>1e-13_wp .or. abs(sum_out/direct_out-1)>1e-13_wp) &
          error stop 'summed charge conservation'
     if (abs(computed_heat/ev2erg-direct_heat)/gross_heat>1e-13_wp) &
          error stop 'staged thermal energy accounting'
     write(*,'(a,f9.1,a,3es12.4)') 'T=',temperature,' closure deviations: ', &
          abs(nem-direct_ne),abs(sum_in/direct_in-1),abs(sum_out/direct_out-1)
  end do
  close(71)
  print *, 'PASS: staged ions=',nstage,', transitions=',ntransition, &
       ', direct/cached charge and mean-heat accounting'
  call MPI_FINALIZE(ierr)
end program check_added_charge_transfer_heat
