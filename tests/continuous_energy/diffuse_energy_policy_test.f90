program diffuse_energy_policy_test
  use, intrinsic :: iso_fortran_env, only : real64
  use ion_energy_policy_mod, only : assign_ion_packet_energy
  implicit none
  integer, parameter :: wp = real64
  real(wp) :: packet_energy
  real(wp), parameter :: sampled(5) = [13.731_wp, 24.400_wp, 24.600_wp, &
                                       54.900_wp, 19.913_wp]
  real(wp), parameter :: grouped(5) = [13.900_wp, 23.950_wp, 25.100_wp, &
                                       58.000_wp, 19.500_wp]
  integer :: i

  do i = 1, size(sampled)
     call assign_ion_packet_energy(packet_energy, sampled(i), grouped(i), 'grouped')
     if (packet_energy /= grouped(i)) &
        error stop 'grouped diffuse energy did not retain bin center'

     call assign_ion_packet_energy(packet_energy, sampled(i), grouped(i), 'continuous')
     if (packet_energy /= sampled(i)) &
        error stop 'continuous diffuse energy did not retain sampled eph'
  end do

  write(*,'(a)') 'DIFFUSE ENERGY SURVIVAL GATE: PASS'
end program diffuse_energy_policy_test
