program check_outer_shell_photoionization
  use define
  use species_mod, only: species_setup,n_elements,elem_nstage,elem_eth, &
                         species_cross_section,outer_shell_photoionization, species_packet_cross_sections, &
                         species_cached_opacity, species_rate_path_setup, species_rate_path_reset, &
                         species_rate_path_score, elem_abund, flat_ie, flat_it, stage_frac, metal_rate_path
  use ion_packet_mod, only: build_ion_packet_physics
  use mpi
  implicit none
  character(len=PATH_LEN) :: directory,output,input,mode
  real(kind=wp),parameter :: energy(15)=[10._wp,13.6_wp,20._wp,30._wp, &
       40._wp,50._wp,60._wp,70._wp,80._wp,90._wp,100._wp,200._wp, &
       500._wp,1000._wp,10000._wp]
  real(kind=wp),parameter :: threshold_factor(3)=[0.99_wp,1._wp,1.01_wp]
  real(kind=wp) :: E,outer,outer_heat,total,total_heat
  type(photon_type) :: photon
  real(kind=wp) :: opacity
  integer :: ierr,ie,it,j,u,k,ii,limit
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%p_rank,ierr)
  mpar%hostcomm=MPI_COMM_WORLD
  call get_command_argument(1,directory)
  call get_command_argument(2,output)
  call get_command_argument(3,input)
  if(len_trim(input)>0) then
    block
      use setup_mod, only: read_parameters_with_implicit_defaults
      call read_parameters_with_implicit_defaults(trim(input))
    end block
  endif
  call get_command_argument(4,mode)
  limit=18
  if(trim(mode)=="below") limit=3
  par%use_metals=.true.
  par%atomic_dir=trim(directory)
  par%cooling_model='low_density'
  par%abund_C=1e-4_wp;par%abund_N=1e-4_wp;par%abund_O=1e-4_wp
  par%abund_Ne=1e-4_wp;par%abund_S=1e-4_wp;par%abund_Ar=1e-4_wp
  par%abund_Mg=1e-4_wp;par%abund_Fe=1e-4_wp;par%abund_Si=1e-4_wp
  par%abund_Cl=1e-4_wp;par%abund_Ca=1e-4_wp
  call species_setup(1)
  if(n_elements/=11) stop 1
  call species_rate_path_setup(1)
  do ii=1,n_elements
    stage_frac(ii)%g=0._wp
  enddo
  open(newunit=u,file=trim(output),status='new')
  do ie=1,n_elements
    do it=1,elem_nstage(ie)-1
      do j=1,limit
        if(j<=15) then
          E=energy(j)
        else
          E=elem_eth(ie,it)*threshold_factor(j-15)
        endif
        call outer_shell_photoionization(ie,it,E,outer,outer_heat)
        call species_cross_section(ie,it,E,total,total_heat)
        if(outer<0.or.outer_heat<0.or.total<outer) stop 2
        if(total/=outer.or.total_heat/=outer_heat) error stop 'Selected outer component differs'
        photon%energy_eV=E
        call build_ion_packet_physics(photon)
        call species_rate_path_reset()
        call species_rate_path_score(E,photon%ionphys%sigma_metal,photon%ionphys%heat_metal, &
                                     photon%ionphys%n_metal,1,.75_wp)
        k=0
        do ii=1,photon%ionphys%n_metal
          if(flat_ie(ii)==ie.and.flat_it(ii)==it) k=ii
        enddo
        if(k>0) then
          if(photon%ionphys%sigma_metal(k)/=total.or.photon%ionphys%heat_metal(k)/=total_heat) &
             error stop 'Packet absorption or heating differs from direct coefficient'
        else
          if(total/=0.or.total_heat/=0) error stop 'Packet omitted active transition'
        endif
        if(abs(metal_rate_path(ie)%raw(1,it,1)-.75_wp*total/(E*ev2erg)) &
           >1e-12_wp*max(total/(E*ev2erg),1e-300_wp)) error stop 'Ionization estimator differs'
        if(abs(metal_rate_path(ie)%raw(2,it,1)-.75_wp*total_heat) &
           >1e-12_wp*max(total_heat,1e-300_wp)) error stop 'Heating estimator differs'
        stage_frac(ie)%g(it,1)=.3_wp
        opacity=species_cached_opacity(photon%ionphys%sigma_metal,photon%ionphys%n_metal,1)
        if(abs(opacity-.3_wp*elem_abund(ie)*total)>1e-12_wp*max(elem_abund(ie)*total,1e-300_wp)) &
           error stop 'Cached metal opacity differs'
        stage_frac(ie)%g(it,1)=0._wp
        write(u,'(2(i0,1x),5(es26.16e3,1x))') ie,it,E,outer,outer_heat,total,total_heat
      enddo
    enddo
  enddo
  close(u)
  call MPI_FINALIZE(ierr)
end program check_outer_shell_photoionization
