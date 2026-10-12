program check_stage_registry
! Synthetic seven-stage Fe tests storage and consumers, not Fe atomic physics.
  use define
  use species_mod
  use gas_state_mod, only : gas_nH, gas_ne, gas_Te, gas_xHI, gas_xHII, &
                           gas_xHeI, gas_xHeII, gas_xHeIII, gas_nleaf
  use octree_mod, only : amr_grid
  use previous_state_mod, only : load_previous_state
  use iofile_mod
  use mpi
  implicit none
  real(kind=wp) :: frac(MAX_ION_STAGES), gamma(MAX_ION_STAGES-1,1)
  real(kind=wp) :: nem, direct_ne, cxi, cxr, sig, heat
  real(kind=wp) :: sin(MAX_ION_STAGES,MAX_METAL_ELEMENTS), sout(MAX_ION_STAGES,MAX_METAL_ELEMENTS)
  real(kind=wp), allocatable :: readback(:,:)
  type(io_file_type) :: file
  type(ion_packet_physics_type) :: packet
  character(len=PATH_LEN) :: directory, output
  character(len=16) :: mode
  integer :: ierr, i, status

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  call get_command_argument(1, directory)
  call get_command_argument(2, output)
  call get_command_argument(3, mode)
  mpar%hostcomm = MPI_COMM_WORLD
  mpar%h_rank = 0
  par%atomic_dir = trim(directory)
  par%abund_C = 0.0_wp; par%abund_N = 0.0_wp; par%abund_O = 0.0_wp
  par%abund_Ne = 0.0_wp; par%abund_S = 0.0_wp; par%abund_Ar = 0.0_wp
  par%abund_Mg = 0.0_wp; par%abund_Si = 0.0_wp; par%abund_Cl = 0.0_wp
  par%abund_Ca = 0.0_wp; par%abund_Fe = 1.0e-4_wp
  call species_setup(1)
  call require(n_elements == 1 .and. elem_nstage(1) == MAX_ION_STAGES, 'seven-stage registry')
  gamma = 2.0e-10_wp
  call species_set_photoionization_rates(1, gamma)
  call species_fractions(1, 1, 1.0e4_wp, 100.0_wp, 0.0_wp, 100.0_wp, frac)
  call require(abs(sum(frac)-1.0_wp) < 1.0e-14_wp .and. frac(7) > 0.0_wp, 'terminal population')
  direct_ne = 0.0_wp
  do i = 1, MAX_ION_STAGES
     direct_ne = direct_ne + 100.0_wp*elem_abund(1)*frac(i)*real(i-1,wp)
  end do
  call species_ne_prepare(1, 1.0e4_wp)
  call species_metal_sums(100.0_wp, 100.0_wp, 0.0_wp, 100.0_wp, nem, cxi, cxr, sin, sout)
  call require(abs(nem-direct_ne) < 1.0e-14_wp, 'direct and cached charge closure')
  call require(cxi == 0.0_wp .and. cxr == 0.0_wp, 'synthetic charge exchange absent')
  sin = -1.0_wp
  call species_metal_sums(100.0_wp,100.0_wp,0.0_wp,100.0_wp,nem,cxi,cxr,cx_in_stage=sin)
  call require(all(sin == 0.0_wp), 'ionizing reaction output alone')
  sout = -1.0_wp
  call species_metal_sums(100.0_wp,100.0_wp,0.0_wp,100.0_wp,nem,cxi,cxr,cx_out_stage=sout)
  call require(all(sout == 0.0_wp), 'recombining reaction output alone')

  allocate(gas_nH(1), gas_ne(1), gas_Te(1), gas_xHI(1), gas_xHII(1))
  gas_nH = 100.0_wp; gas_ne = 100.0_wp; gas_Te = 1.0e4_wp
  gas_xHI = 0.0_wp; gas_xHII = 1.0_wp
  amr_grid%nleaf = 1
  call species_stage_cache_refresh()
  packet%energy_eV = 90.0_wp
  call species_packet_cross_sections(packet%energy_eV, packet%sigma_metal, packet%n_metal, packet%heat_metal)
  call require(packet%n_metal == 6, 'all six packet transitions')
  call species_cross_section(1, 6, packet%energy_eV, sig, heat)
  call require(sig > 0.0_wp .and. heat > 0.0_wp, 'sixth transition cross section')

  status = 0
  call io_open_new(file, trim(output), status)
  call species_write(file, status)
  call io_close(file, status)
  call require(status == 0, 'write seven-stage state')
  call io_open_old(file, trim(output), status)
  ! The first appended image is the stage array (no empty primary image).
  allocate(readback(7,1))
  call io_read_image(file, readback, status)
  call require(status == 0, 'read seven-stage state')
  call require(size(readback,1) == 7 .and. size(readback,2) == 1, 'stage output shape')
  call require(maxval(abs(readback(:,1)-frac)) < 1.0e-14_wp, 'stage output values')
  deallocate(readback)
  call io_move_to_next_section(file, status)
  allocate(readback(6,1))
  call io_read_image(file, readback, status)
  call require(status == 0 .and. size(readback,1) == 6, 'transition output shape')
  call require(maxval(abs(readback-gamma)) < 1.0e-24_wp, 'transition output values')
  call io_close(file, status)
  deallocate(readback)

  ! An old file can contain disabled elements before an active element.
  ! Three Fe transitions from the old run must survive; three new ones
  ! start at zero and force fresh transport and reconvergence.
  gas_nleaf = 1
  allocate(gas_xHeI(1), gas_xHeII(1), gas_xHeIII(1))
  gas_xHeI = 1.0_wp; gas_xHeII = 0.0_wp; gas_xHeIII = 0.0_wp
  amr_grid%car = .true.
  amr_grid%nx = 1; amr_grid%ny = 1; amr_grid%nz = 1
  amr_grid%dcell = 1.0_wp
  call io_open_new(file, trim(output)//'.old.fits', status)
  call io_append_image(file, reshape([0.5_wp,0.5_wp,0.5_wp],[1,3]), status, bitpix=-64)
  call io_put_keyword(file, 'EXTNAME', 'LeafXYZ', 'synthetic center', status)
  call scalar_field('LeafSize',1.0_wp)
  call scalar_field('x_HI',0.0_wp)
  call scalar_field('x_HII',1.0_wp)
  call scalar_field('x_HeI',1.0_wp)
  call scalar_field('x_HeII',0.0_wp)
  call scalar_field('x_HeIII',0.0_wp)
  call scalar_field('n_e',100.0_wp)
  call scalar_field('T_e',1.0e4_wp)
  call io_append_image(file, gamma(1:3,:), status, bitpix=-64)
  call io_put_keyword(file, 'EXTNAME', 'Gamma_c_stages', 'disabled element', status)
  call io_append_image(file, gamma(1:3,:), status, bitpix=-64)
  call io_put_keyword(file, 'EXTNAME', 'Gamma_fe_stages', 'old transition rates', status)
  call io_close(file, status)
  call require(status == 0, 'write old restart fixture')
  par%require_convergence = .false.
  par%gas_niter = 8
  if (trim(mode) == 'reject') par%gas_niter = 0
  call load_previous_state(trim(output)//'.old.fits')
  call require(par%require_convergence, 'expanded restart requires convergence')
  call io_open_new(file, trim(output)//'.restarted.fits', status)
  call species_write(file, status)
  call io_close(file, status)
  call io_open_old(file, trim(output)//'.restarted.fits', status)
  call io_move_to_next_section(file, status)
  allocate(readback(6,1))
  call io_read_image(file, readback, status)
  call require(status == 0, 'read restarted transition rates')
  call require(maxval(abs(readback(1:3,:)-gamma(1:3,:))) < 1.0e-24_wp, 'common restart rates preserved')
  call require(all(readback(4:6,:) == 0.0_wp), 'new restart rates zero')
  call io_close(file, status)
  write(*,'(a)') 'PASS: seven-stage registry, charge closure, packet transitions, and output readback'
  write(*,'(a)') 'PASS: old restart rate mapping and reconvergence requirement'
  call MPI_FINALIZE(ierr)
contains
  subroutine scalar_field(name, value)
    character(len=*), intent(in) :: name
    real(kind=wp), intent(in) :: value
    call io_append_image(file, [value], status, bitpix=-64)
    call io_put_keyword(file, 'EXTNAME', name, 'synthetic restart field', status)
  end subroutine scalar_field
  subroutine require(condition, name)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name
    if (.not. condition) then
       write(*,'(2a)') 'FAIL: ', name
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
  end subroutine require
end program check_stage_registry
