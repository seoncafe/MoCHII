program check_ionization_thresholds
! Loaded ionization thresholds against the compiled potentials, absorption at
! and below each threshold, the photoelectron energy at threshold, the
! threshold ordering of the packet cross sections, and the terminal-stage
! refusal of a band that reaches above the last carried stage.
  use define
  use species_mod
  use metal_ionization_potential_mod
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: ierr, ie, it, z, stage, count, nt
  integer, parameter :: numbers(11) = [6,7,8,10,16,18,12,26,14,17,20]
  real(kind=wp) :: eth, sig, heat, energy, direct, sigma(MAX_METAL_TRANSITIONS)
  real(kind=wp) :: reject_energy
  ! A two-bin band, 13.598-50-200 eV, equal luminosity in each bin.
  real(kind=wp), parameter :: eedge(3) = [13.598_wp, 50.0_wp, 200.0_wp]
  real(kind=wp), parameter :: lum(2) = [1.0_wp, 1.0_wp]
  character(len=PATH_LEN) :: directory, output
  character(len=32) :: energy_text
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  mpar%hostcomm = MPI_COMM_WORLD
  mpar%h_rank = 0
  call get_command_argument(1, directory)
  call get_command_argument(2, output)
  call get_command_argument(3, energy_text)
  par%atomic_dir = trim(directory)
  par%cooling_model = 'low_density'
  par%abund_C=1e-4_wp; par%abund_N=1e-4_wp; par%abund_O=1e-4_wp
  par%abund_Ne=1e-4_wp; par%abund_S=1e-4_wp; par%abund_Ar=1e-4_wp
  par%abund_Mg=1e-4_wp; par%abund_Fe=1e-4_wp; par%abund_Si=1e-4_wp
  par%abund_Cl=1e-4_wp; par%abund_Ca=1e-4_wp
  call species_setup(1)
  if (n_elements /= 11) error stop 'incomplete element registry'
  if (len_trim(energy_text) > 0) then
     read(energy_text, *) reject_energy
     call species_terminal_photoionization_check(reject_energy, eedge, lum)
     error stop 'band above a terminal stage was accepted'
  end if
  ! The production band reaches 100 eV; every terminal stage lies above it.
  call species_terminal_photoionization_check(100.0_wp, eedge, lum)
  open(unit=71, file=trim(output), status='replace')
  do ie=1,11
     z=numbers(ie)
     do stage=1,7
        write(71,'(2i5,es27.17e3)') z,stage,metal_ionization_potential(z,stage)
     end do
  end do
  close(71)
  count=0
  do ie=1,n_elements
     do it=1,elem_nstage(ie)-1
        count=count+1
        eth=elem_eth(ie,it)
        if (abs(eth/metal_ionization_potential(numbers(ie),it)-1.0_wp) > 1e-12_wp) &
             error stop 'loaded ETH does not match physical potential'
        call species_cross_section(ie,it,eth*(1.0_wp-1e-9_wp),sig,heat)
        if (sig /= 0 .or. heat /= 0) error stop 'subthreshold absorption'
        energy=eth*(1.0_wp+1e-9_wp)
        call species_cross_section(ie,it,energy,sig,heat)
        if (.not. ieee_is_finite(sig) .or. sig <= 0 .or. heat <= 0) &
             error stop 'invalid absorption immediately above threshold'
        ! Just above the outer edge only the outer electron can be removed:
        ! the photoelectron carries E - IP.
        if (abs(heat/sig-(1.0_wp-eth/energy)) > 1e-14_wp) &
             error stop 'threshold photoelectron energy budget'
        direct=0
        do z=1,n_elements
           do stage=1,elem_nstage(z)-1
              call species_cross_section(z,stage,energy,sig,heat)
              direct=direct+sig
           end do
        end do
        call species_packet_cross_sections(energy,sigma,nt)
        if (abs(sum(sigma)/direct-1.0_wp) > 1e-13_wp) error stop 'packet threshold ordering'
     end do
  end do
  write(*,'(a,i0,a)') 'PASS: ', count, ' loaded transitions, 77 compiled potentials'
  call MPI_FINALIZE(ierr)
end program check_ionization_thresholds
