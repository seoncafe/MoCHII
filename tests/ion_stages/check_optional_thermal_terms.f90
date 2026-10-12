program check_optional_thermal_terms
! thermal_mod::net_rate against the sum of its consumers, assembled here
! term by term, for metal-only, helium-only and combined gas, with the
! switches par%charge_exchange_heating, par%metal_freefree and
! par%neutral_hydrogen_excitation as the input sets them.
  use define
  use mpi
  use setup_mod, only: read_parameters_with_implicit_defaults
  use thermal_mod, only: net_rate
  use gas_state_mod, only: gas_nH
  use gas_rates_mod, only: gamma_HI, gamma_HeI, gamma_HeII, heat_HI, heat_HeI, heat_HeII
  use cooling_mod, only: cooling_setup, cooling_total
  use species_mod, only: species_setup, metal_cooling, metal_heating, &
       metal_ionization_cooling, metal_freefree, metal_charge_exchange_heating, &
       metal_cooling_H
  use ion_balance_mod, only: helium_charge_exchange_heating
  implicit none
  character(len=PATH_LEN) :: input
  real(kind=wp), parameter :: temperatures(3)=[5000.0_wp,10000.0_wp,20000.0_wp]
  real(kind=wp), parameter :: radiation(3)=[1e-14_wp,1e-12_wp,1e-10_wp]
  real(kind=wp) :: xHI,xHII,xHeI,xHeII,xHeIII,ne,net,base,cx,ff,expected,T,nH
  integer :: ierr,i,j,k,count
  logical :: metals
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD
  ! Compiled defaults: both kinetic-energy terms and the H-impact excitation on.
  if (.not.par%charge_exchange_heating.or..not.par%metal_freefree &
      .or..not.par%neutral_hydrogen_excitation) stop 1
  call get_command_argument(1,input)
  call read_parameters_with_implicit_defaults(trim(input))
  par%atomic_dir='data/atomic'
  par%data_dir='data'
  par%cooling_model='low_density'
  par%hei_metastable_coupling=.false.
  par%hei_on_the_spot=.false.;par%heii_on_the_spot=.false.
  par%grain_pe=.false.;par%use_sec_ion=.false.
  par%metal_heat=.true.;par%metal_ne=.true.
  par%helium_charge_exchange=.true.
  par%abund_C=2.2e-4_wp;par%abund_O=3.3e-4_wp
  par%abund_N=0;par%abund_Ne=0;par%abund_S=0;par%abund_Ar=0
  par%abund_Mg=0;par%abund_Fe=0;par%abund_Si=0;par%abund_Cl=0;par%abund_Ca=0
  call species_setup(1)
  call cooling_setup()
  allocate(gas_nH(1),gamma_HI(1),gamma_HeI(1),gamma_HeII(1))
  allocate(heat_HI(1),heat_HeI(1),heat_HeII(1))
  nH=100;gas_nH=nH;count=0
  ! Metal-only, helium-only, and their combined kinetic-energy consumers.
  do k=1,3
    metals=k/=2;par%use_metals=metals
    par%abund_He=0
    if(k/=1) par%abund_He=0.1_wp
    do i=1,3
      T=temperatures(i)
      do j=1,3
        gamma_HI=radiation(j);gamma_HeI=0.2_wp*radiation(j);gamma_HeII=0
        heat_HI=gamma_HI*5*ev2erg;heat_HeI=gamma_HeI*5*ev2erg;heat_HeII=0
        xHI=0.5_wp;xHII=0.5_wp;xHeI=0.5_wp;xHeII=0.5_wp;xHeIII=0
        call net_rate(1,T,.false.,xHI,xHII,xHeI,xHeII,xHeIII,ne,net)
        base=nH*(xHI*heat_HI(1)+par%abund_He*xHeI*heat_HeI(1)) &
             -cooling_total(T,nH,ne,xHI,xHII,xHeI,xHeII,xHeIII)
        ff=0;cx=0
        if(metals) then
          base=base-metal_cooling(1,T,nH,ne,nH*xHI,nH*xHII) &
               +metal_heating(1,T,nH,ne,nH*xHI,nH*xHII) &
               -metal_ionization_cooling(1,T,nH,ne,nH*xHI,nH*xHII)
          if(par%neutral_hydrogen_excitation) &
               base=base-metal_cooling_H(1,T,nH,ne,nH*xHI,nH*xHII)
          ff=metal_freefree(1,T,nH,ne,nH*xHI,nH*xHII)
          cx=metal_charge_exchange_heating(1,T,nH,ne,nH*xHI,nH*xHII)
          if(ff<=0) stop 2
        endif
        if(par%abund_He>0) cx=cx+helium_charge_exchange_heating(T,nH*xHI, &
             nH*xHII,par%abund_He*nH*xHeI,par%abund_He*nH*xHeII)
        expected=base
        if(par%metal_freefree) expected=expected-ff
        if(par%charge_exchange_heating) expected=expected+cx
        if(abs(net-expected)>1e-12_wp*max(abs(base),abs(cx),abs(ff),1e-35_wp)) stop 3
        if(cx==0) stop 4
        count=count+1
        write(*,'(a,1x,i1,7(1x,es24.16))') 'THERMAL',k,T,radiation(j),ne,net,base,cx,ff
      enddo
    enddo
  enddo
  write(*,'(a,i0)') 'PASS: optional thermal consumers cases=',count
  call MPI_FINALIZE(ierr)
end program check_optional_thermal_terms
! thermal_mod::net_rate against the sum of its consumers, assembled here
! term by term, for metal-only, helium-only and combined gas, with the
! switches par%charge_exchange_heating, par%metal_freefree and
! par%neutral_hydrogen_excitation as the input sets them.
