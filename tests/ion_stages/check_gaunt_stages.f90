program check_gaunt_stages
  use define
  use cooling_mod, only : cooling_setup, gbar_ff
  use gaunt_vh14_mod, only : gauntff_vh14
  use mpi
  implicit none
  real(kind=wp) :: t, value, reference, worst, lnu, du, u, rykt
  ! Reference: the full thermal average Int_0^inf g_ff(u) exp(-u) du of the
  ! van Hoof (2014) table, which is the rule for Z = 5, 6 at 1e2 and 1e3 K
  ! (gamma^2 > 1e3), by the midpoint rule in ln u on [1e-12, 700]. The
  ! earlier reference averaged over 1e-4 <= u <= 20 only, normalized by
  ! exp(-1e-4) - exp(-20); that is not the quantity gbar_ff defines and
  ! differed from the full integral by 1.7e-4 .. 2.2e-4 here. Tolerance:
  ! the accuracy of cooling_mod's 400-point quadrature, measured against an
  ! independent integral (tools/python/make_gbar_ff_table.py) at 8.1e-6 at
  ! worst over all charges and 1e2-1e6 K.
  integer, parameter :: NQ = 400000
  real(kind=wp), parameter :: U_LOW = 1.0e-12_wp, U_HIGH = 700.0_wp
  integer :: ierr, i, z, j, unit
  character(len=PATH_LEN) :: output
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  call get_command_argument(1, output)
  par%atomic_dir = 'data/atomic'
  par%data_dir = 'data'
  par%cooling_model = 'low_density'
  par%gaunt_vh14 = .false.
  call cooling_setup()
  open(newunit=unit, file=trim(output), status='replace')
  do i = 0, 40
     t = 10.0_wp**(2.0_wp+0.1_wp*real(i,wp))
     do z = 1, MAX_ION_STAGES-1
        value = gbar_ff(t,z)
        if (value <= 0.0_wp) stop 1
     end do
     write(unit,'(7es25.16)') t, (gbar_ff(t,z),z=1,MAX_ION_STAGES-1)
  end do
  close(unit)
  worst = 0.0_wp
  do z = 5, 6
     do i = 2, 3
        t = 10.0_wp**i
        rykt = 13.605693_wp*ev2erg/(kboltz_cgs*t)
        reference = 0.0_wp
        do j = 1, NQ
           lnu = log(U_LOW) + (real(j,wp)-0.5_wp)*log(U_HIGH/U_LOW)/real(NQ,wp)
           u = exp(lnu)
           du = u*log(U_HIGH/U_LOW)/real(NQ,wp)
           reference = reference + gauntff_vh14(real(z,wp),t,u/rykt)*exp(-u)*du
        end do
        value = gbar_ff(t,z)
        worst = max(worst, abs(value/reference-1.0_wp))
        write(*,'(a,i0,a,es12.4,a,es12.4)') 'GAUNT Z=',z,' T=',t,' relative quadrature difference=', &
                                        abs(value/reference-1.0_wp)
     end do
  end do
  if (worst > 1.0e-5_wp) stop 1
  write(*,'(a,es12.4)') 'PASS: 246 positive Gaunt nodes; maximum high-charge quadrature difference=', worst
  call MPI_FINALIZE(ierr)
end program check_gaunt_stages
