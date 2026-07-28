module ion_energy_policy_mod
!---------------------------------------------------------------------------
! Packet-energy assignment during the grouped-to-continuous migration.
!
! Keeping this small policy in one place prevents a launch path from silently
! retaining a sampled energy in grouped mode or discarding it in continuous
! mode.  The latter mode remains guarded at setup until the complete H/He
! transport/rate vertical slice is activated.
!---------------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
  private

  integer, parameter :: wp = real64
  public :: assign_ion_packet_energy

contains

  subroutine assign_ion_packet_energy(packet_energy, sampled_energy, &
       grouped_energy, mode)
    real(wp), intent(out) :: packet_energy
    real(wp), intent(in) :: sampled_energy, grouped_energy
    character(len=*), intent(in) :: mode

    select case (trim(mode))
    case ('grouped')
       packet_energy = grouped_energy
    case ('continuous')
       packet_energy = sampled_energy
    case default
       error stop 'ion energy policy: mode must be grouped or continuous'
    end select
  end subroutine assign_ion_packet_energy

end module ion_energy_policy_mod
