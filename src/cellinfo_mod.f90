! MoCHII: copied from MoCafe_v2.00/src/cellinfo_mod.f90 (2026-07-11)
module cellinfo_mod
!--- Grid-agnostic cell accessors (MoCafe v2.00) so the dust-emission stages
!--- (J_lambda tally, dust emission, B&W, MRW, all-sky) treat the Cartesian
!--- grid and the AMR octree the same way.  A "cell" is addressed by a linear
!--- id: for 'car' it is (k-1)*nx*ny + (j-1)*nx + i over the regular grid; for
!--- 'amr' it is the leaf index (1..nleaf).  Opacity, code-unit volume, and
!--- center are returned per cell; the photon's current cell id is derived from
!--- its Cartesian indices ('car') or its leaf index ('amr').
  use define
  use octree_mod, only : amr_grid
  implicit none
  public

contains
  !---------------------------------------------------------------
  integer function ncell_total(grid) result(n)
  type(grid_type), intent(in) :: grid
  if (trim(par%grid_type) == 'amr') then
     n = amr_grid%nleaf
  else
     n = grid%nx*grid%ny*grid%nz
  endif
  end function ncell_total

  !---------------------------------------------------------------
  !--- linear cell id of a photon's current cell.
  integer function cell_id_of_photon(photon, grid) result(ic)
  type(photon_type), intent(in) :: photon
  type(grid_type),   intent(in) :: grid
  if (trim(par%grid_type) == 'amr') then
     ic = photon%icell_amr
  else
     ic = (photon%kcell-1)*grid%nx*grid%ny + (photon%jcell-1)*grid%nx + photon%icell
  endif
  end function cell_id_of_photon

  !---------------------------------------------------------------
  real(kind=wp) function cell_rhokap(grid, ic) result(rhk)
  type(grid_type), intent(in) :: grid
  integer,         intent(in) :: ic
  integer :: i, j, k
  if (trim(par%grid_type) == 'amr') then
     rhk = amr_grid%rhokap(ic)
  else
     call car_ijk(grid, ic, i, j, k)
     rhk = grid%rhokap(i,j,k)
  endif
  end function cell_rhokap

  !---------------------------------------------------------------
  real(kind=wp) function cell_volume(grid, ic) result(vol)
  type(grid_type), intent(in) :: grid
  integer,         intent(in) :: ic
  integer :: icell
  if (trim(par%grid_type) == 'amr') then
     icell = amr_grid%icell_of_leaf(ic)
     vol   = (2.0_wp*amr_grid%ch(icell))**3
  else
     vol   = grid%dx*grid%dy*grid%dz
  endif
  end function cell_volume

  !---------------------------------------------------------------
  subroutine cell_center(grid, ic, x, y, z)
  type(grid_type), intent(in)  :: grid
  integer,         intent(in)  :: ic
  real(kind=wp),   intent(out) :: x, y, z
  integer :: i, j, k, icell
  if (trim(par%grid_type) == 'amr') then
     icell = amr_grid%icell_of_leaf(ic)
     x = amr_grid%cx(icell);  y = amr_grid%cy(icell);  z = amr_grid%cz(icell)
  else
     call car_ijk(grid, ic, i, j, k)
     x = (grid%xface(i)+grid%xface(i+1))*0.5_wp
     y = (grid%yface(j)+grid%yface(j+1))*0.5_wp
     z = (grid%zface(k)+grid%zface(k+1))*0.5_wp
  endif
  end subroutine cell_center

  !---------------------------------------------------------------
  !--- sample a uniform random position inside cell ic and set the photon's
  !--- cell indices (car) or leaf index (amr) accordingly.
  subroutine cell_random_position(grid, ic, photon)
  use random, only : rand_number
  type(grid_type),   intent(in)    :: grid
  integer,           intent(in)    :: ic
  type(photon_type), intent(inout) :: photon
  integer :: i, j, k, icell
  real(kind=wp) :: h
  if (trim(par%grid_type) == 'amr') then
     icell = amr_grid%icell_of_leaf(ic)
     h = amr_grid%ch(icell)
     photon%x = amr_grid%cx(icell) + (2.0_wp*rand_number()-1.0_wp)*h
     photon%y = amr_grid%cy(icell) + (2.0_wp*rand_number()-1.0_wp)*h
     photon%z = amr_grid%cz(icell) + (2.0_wp*rand_number()-1.0_wp)*h
     photon%icell_amr = ic
  else
     call car_ijk(grid, ic, i, j, k)
     photon%x = grid%xface(i) + grid%dx*rand_number()
     photon%y = grid%yface(j) + grid%dy*rand_number()
     photon%z = grid%zface(k) + grid%dz*rand_number()
     photon%icell = i;  photon%jcell = j;  photon%kcell = k
  endif
  end subroutine cell_random_position

  !---------------------------------------------------------------
  !--- place a photon at the center of cell ic and set its cell/leaf index
  !--- (used by the emergent-SED ray tracing, which then walks via the
  !--- raytrace_to_edge procedure pointer).
  subroutine cell_center_photon(grid, ic, photon)
  type(grid_type),   intent(in)    :: grid
  integer,           intent(in)    :: ic
  type(photon_type), intent(inout) :: photon
  integer :: i, j, k
  call cell_center(grid, ic, photon%x, photon%y, photon%z)
  if (trim(par%grid_type) == 'amr') then
     photon%icell_amr = ic
  else
     call car_ijk(grid, ic, i, j, k)
     photon%icell = i;  photon%jcell = j;  photon%kcell = k
  endif
  end subroutine cell_center_photon

  !---------------------------------------------------------------
  pure subroutine car_ijk(grid, ic, i, j, k)
  type(grid_type), intent(in)  :: grid
  integer,         intent(in)  :: ic
  integer,         intent(out) :: i, j, k
  k = (ic-1)/(grid%nx*grid%ny) + 1
  j = mod((ic-1)/grid%nx, grid%ny) + 1
  i = mod(ic-1, grid%nx) + 1
  end subroutine car_ijk

end module cellinfo_mod
