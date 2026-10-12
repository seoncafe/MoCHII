program check_hei_metastable_input
! The He I 2^3S input defaults as the namelist reader resolves them, and,
! with the coupling off, the ion solve independent of the leaf whose
! metastable estimator it would read.
  use define
  use setup_mod, only: read_parameters_with_implicit_defaults
  use ion_balance_mod, only: solve_ion_cell
  use mpi
  implicit none
  character(len=PATH_LEN) :: input
  real(kind=wp) :: x, h1, h2, ne, x0, h10, h20, ne0
  integer :: ierr
  call MPI_INIT(ierr)
  ! Compiled defaults (define.f90).
  if (.not. par%hei_metastable_coupling) stop 1
  if (.not. par%helium_charge_exchange .or. .not. par%metal_freefree) stop 2
  call get_command_argument(1,input)
  call read_parameters_with_implicit_defaults(trim(input))
  write(*,'(a,6(1x,l1),1x,es24.16)') 'INPUT', par%hei_metastable_coupling, &
       par%hei_metastable, par%add_fuv, par%helium_charge_exchange, &
       par%metal_freefree, par%charge_exchange_heating, par%efuv_min
  ! The disabled coupling must not use a leaf-indexed metastable estimator.
  ! Compare the invoked solver with and without its optional leaf argument.
  if (.not. par%hei_metastable_coupling) then
     par%use_metals=.false.
     par%heii_on_the_spot=.false.
     par%hei_on_the_spot=.false.
     call solve_ion_cell(2.e-10_wp,1.e-10_wp,1.e-12_wp,100._wp,1.e4_wp, &
                        .true.,x,h1,h2,ne,1)
     call solve_ion_cell(2.e-10_wp,1.e-10_wp,1.e-12_wp,100._wp,1.e4_wp, &
                        .true.,x0,h10,h20,ne0)
     if (x/=x0 .or. h1/=h10 .or. h2/=h20 .or. ne/=ne0) stop 3
     write(*,'(a)') 'PASS: disabled metastable population feedback'
  end if
  call MPI_FINALIZE(ierr)
end program check_hei_metastable_input
! The He I 2^3S input defaults as the namelist reader resolves them, and,
! with the coupling off, the ion solve independent of the leaf whose
! metastable estimator it would read.
