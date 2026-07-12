! MoCHII: copied from MoCafe_v2.00/src/grid_mod_car.f90 (2026-07-11)
module grid_mod
  use define
  use utility
  use memory_mod
  !use read_mod
  implicit none
contains

  !---------------------------------------------------------------------------------
  subroutine grid_create(grid)
  use define
  use mpi
  use read_mod
  implicit none
  type(grid_type), intent(inout) :: grid

  ! local variables
  integer :: status, unit
  integer :: i,j,k
  integer :: idx1,idx2,idx3,icen,jcen,kcen
  integer :: nxcen,nycen,nzcen
  real(kind=wp) :: taueq,taupole
  real(kind=wp) :: dens
  real(kind=wp) :: xx,yy,zz,rr,pp,tanp
  real(kind=wp) :: rmin_positive
  character(len=128) :: dens_out

  integer :: ierr

  ! external density, temperature, and velocity field files.
  if (len_trim(par%density_file) > 0) then
     !call get_dimension(trim(par%density_file),par%nx,par%ny,par%nz,par%xmax,par%ymax,par%zmax,status,&
     call get_dimension(trim(par%density_file),par%nx,par%ny,par%nz,status,&
                                               reduce_factor=par%reduce_factor)
     if (status /= 0 .and. mpar%p_rank == 0) then
        write(*,*) 'Density File ',trim(par%density_file),' does not have Dimensional Information!'
        write(*,*) 'You need to specify xmax, ymax, zmax, nx, ny, nz in input parameter file!'
     endif
  endif

  !--- grid's cell faces.
  grid%nx   = par%nx
  grid%ny   = par%ny
  grid%nz   = par%nz
  grid%xmax = par%xmax
  grid%ymax = par%ymax
  grid%zmax = par%zmax
  grid%rmax = par%rmax

  ! grid%i0, grid%j0, grid%k0 = cell indices when the photon met a boundary and refelected into the grid system.
  ! center of grid system is always located at (x,y,z) = (0,0,0).
  if (par%xyz_symmetry .or. par%z_symmetry) then
     if ((par%nz/2)*2 == par%nz) then
        grid%dz   = grid%zmax/grid%nz
        grid%zmin = 0.0_wp
        grid%k0   = 1
     else
        grid%dz   = grid%zmax/(grid%nz-0.5_wp)
        grid%zmin = -grid%dz/2.0_wp
        grid%k0   = 2
     endif
  else
     grid%dz   = 2.0_wp*par%zmax/par%nz
     grid%zmin = -grid%zmax
     grid%k0   = 0
  endif
  if (par%xyz_symmetry) then
     if ((par%nx/2)*2 == par%nx) then
        grid%dx   = grid%xmax/grid%nx
        grid%xmin = 0.0_wp
        grid%i0   = 1
     else
        grid%dx   = grid%xmax/(grid%nx-0.5_wp)
        grid%xmin = -grid%dx/2.0_wp
        grid%i0   = 2
     endif
     if ((par%ny/2)*2 == par%ny) then
        grid%dy   = grid%ymax/grid%ny
        grid%ymin = 0.0_wp
        grid%j0   = 1
     else
        grid%dy   = grid%ymax/(grid%ny-0.5_wp)
        grid%ymin = -grid%dy/2.0_wp
        grid%j0   = 2
     endif
  else
     grid%dx   = 2.0_wp*par%xmax/par%nx
     grid%dy   = 2.0_wp*par%ymax/par%ny
     grid%xmin = -grid%xmax
     grid%ymin = -grid%ymax
     grid%i0   = 0
     grid%j0   = 0
  endif
  grid%xrange = grid%xmax - grid%xmin
  grid%yrange = grid%ymax - grid%ymin
  grid%zrange = grid%zmax - grid%zmin

  !--- center indices.
  if (par%xyz_symmetry) then
     nxcen = 1
     nycen = 1
     nzcen = 1
  else
     nxcen = (grid%nx+1)/2
     nycen = (grid%ny+1)/2
     nzcen = (grid%nz+1)/2
  endif

  ! allocate shared memories for xface, yface, zface
  call create_shared_mem(grid%xface, [grid%nx+1])
  call create_shared_mem(grid%yface, [grid%ny+1])
  call create_shared_mem(grid%zface, [grid%nz+1])
  if (mpar%h_rank == 0) then
     grid%xface(:) = [ ((i-1)*grid%dx + grid%xmin, i=1,grid%nx+1) ]
     grid%yface(:) = [ ((j-1)*grid%dy + grid%ymin, j=1,grid%ny+1) ]
     grid%zface(:) = [ ((k-1)*grid%dz + grid%zmin, k=1,grid%nz+1) ]
  endif

  ! allocate shared memories for rhokap
  call create_shared_mem(grid%rhokap, [grid%nx,grid%ny,grid%nz])
  if (mpar%h_rank == 0) then
     grid%rhokap(:,:,:) = 0.0_wp
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !---- Density Setup
  if (mpar%h_rank == 0) then
     if (len_trim(par%density_file) == 0) then
        !--- Note distance unit = 1 when no physical distance unit is provided. (2017-06-27)
        !--- SED/dust-emission mode needs a physical length scale for absolute
        !--- J_lambda and dust temperatures, so it keeps the user's distance_unit
        !--- (1 code unit = 1 distance_unit); the legacy scattering path stays
        !--- unit-agnostic (forced to 1) to preserve its output normalization.
        if (.not. par%use_sed) then
           par%distance_unit  = ''
           par%distance2cm    = 1.0_wp
        endif
        grid%rhokap(:,:,:) = 1.0_wp
     else
        call read_3D(trim(par%density_file),grid%rhokap,reduce_factor=par%reduce_factor,centering=par%centering)
        do k=1,grid%nz
           grid%rhokap(:,:,k) = grid%rhokap(:,:,k) * par%cext_dust * par%distance2cm
        enddo
     endif
  endif

  if (mpar%h_rank == 0) then
    ! Cylindrical Geometry
    if (trim(par%geometry) == 'cylinder') then
       if (is_finite(par%density_zscale)) then
         do k=1,grid%nz
            zz = (grid%zface(k)+grid%zface(k+1))/2.0_wp
            if (trim(par%density_zprofile) == 'sech') then
               grid%rhokap(:,:,k) = grid%rhokap(:,:,k)/cosh(zz/par%density_zscale)**2
            else
               grid%rhokap(:,:,k) = grid%rhokap(:,:,k)*exp(-abs(zz) / par%density_zscale)
            endif
         enddo
       endif
       tanp = tan(par%spiral_pitch*pi/180.0_wp)
       do j=1,grid%ny
         yy = (grid%yface(j)+grid%yface(j+1))/2.0_wp
         do i=1,grid%nx
           xx = (grid%xface(i)+grid%xface(i+1))/2.0_wp
           rr = sqrt(xx**2 + yy**2)
           if (par%rmax > 0.0 .and. rr > par%rmax) grid%rhokap(i,j,:) = 0.0_wp
           if (is_finite(par%density_rscale))      grid%rhokap(i,j,:) = grid%rhokap(i,j,:)*exp(-rr / par%density_rscale)
           if (par%spiral_m > 0 .and. rr > 0.0_wp) then
              pp = atan2(yy, xx)
              grid%rhokap(i,j,:) = grid%rhokap(i,j,:) * &
                 (1.0_wp + par%spiral_amp*sin(par%spiral_m*(log(rr)/tanp - pp)))
           endif
         enddo
       enddo
    ! Spherical Geometry
    else if (trim(par%geometry) == 'sphere') then
       do k=1,grid%nz
         zz = (grid%zface(k)+grid%zface(k+1))/2.0_wp
         do j=1,grid%ny
           yy = (grid%yface(j)+grid%yface(j+1))/2.0_wp
           do i=1,grid%nx
             xx = (grid%xface(i)+grid%xface(i+1))/2.0_wp
             rr = sqrt(xx**2 + yy**2 + zz**2)
             if (par%rmax > 0.0 .and. rr > par%rmax) grid%rhokap(i,j,k) = 0.0_wp
             if (is_finite(par%density_rscale)) then
                grid%rhokap(i,j,k) = grid%rhokap(i,j,k)*exp(-rr / par%density_rscale)
             else if (is_finite(par%density_powerlaw_index)) then
                if (par%density_powerlaw_index > 0.0_wp .or. rr > 0.0_wp) then
                   grid%rhokap(i,j,k) = grid%rhokap(i,j,k)*rr**(par%density_powerlaw_index)
                endif
                if (par%density_powerlaw_index < 0.0_wp) then
                   if (i==1 .and. j==1 .and. k==1) then
                      icen = -1
                      jcen = -1
                      kcen = -1
                      rmin_positive = sqrt(par%xmax**2 + par%ymax**2 + par%zmax**2)
                   endif
                   if (rr > 0.0_wp .and. rr < rmin_positive) then
                      idx1 = i
                      idx2 = j
                      idx3 = k
                      rmin_positive = rr
                   endif
                   if (rr == 0.0_wp) then
                      icen = i
                      jcen = j
                      kcen = k
                   endif
                endif
             endif
           enddo
         enddo
       enddo
       if (par%density_powerlaw_index < 0.0_wp .and. icen > 0 .and. jcen > 0 .and. kcen > 0) then
          grid%rhokap(icen,jcen,kcen) = grid%rhokap(idx1,idx2,idx3)
       endif
    else if (trim(par%geometry) == 'Plummer' .and. par%density_rscale > 0.0_wp) then
       if (.not. is_finite(par%density_powerlaw_index)) par%density_powerlaw_index = 1.8
       do k=1,grid%nz
         zz = (grid%zface(k)+grid%zface(k+1))/2.0_wp
         do j=1,grid%ny
           yy = (grid%yface(j)+grid%yface(j+1))/2.0_wp
           do i=1,grid%nx
             xx = (grid%xface(i)+grid%xface(i+1))/2.0_wp
             rr = sqrt(xx**2 + yy**2 + zz**2)
             if (par%rmax > 0.0 .and. rr > par%rmax) grid%rhokap(i,j,k) = 0.0_wp
             grid%rhokap(i,j,k) = grid%rhokap(i,j,k)/(1.0_wp + (rr/par%density_rscale)**2)**(par%density_powerlaw_index/2.0_wp)
           enddo
         enddo
       enddo 
    endif
  endif

  ! Normalize density if density is not given in a physical unit (hydrogen number per cm^3).
  ! Note that normalization is done even when the density is not a constant.
  if (mpar%h_rank == 0) then
     if (par%taumax > 0.0_wp) then
        dens = (2.0_wp*par%taumax)/sum(grid%rhokap(nxcen,nycen,:)*grid%dz)
        grid%rhokap(:,:,:) = grid%rhokap(:,:,:) * dens
     else if (par%tauhomo > 0.0_wp) then
        dens = sum(grid%rhokap)/count(grid%rhokap > 0.0)
        do k=1,grid%nz
           grid%rhokap(:,:,k) = grid%rhokap(:,:,k)*(par%tauhomo/grid%zmax)/dens
        enddo
     endif
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  !--- Calculate equatorial and polar optical depths
  if (mpar%p_rank == 0) then
     if (par%xyz_symmetry) then
        taueq        = sum(grid%rhokap(:,nycen,nzcen))*grid%dx
        taupole      = sum(grid%rhokap(nxcen,nycen,:))*grid%dz
        if ((grid%nx/2)*2 /= grid%nx) taueq   = taueq   - grid%rhokap(1,nycen,nzcen)*(grid%dx/2.0_wp)
        if ((grid%nz/2)*2 /= grid%nz) taupole = taupole - grid%rhokap(nxcen,nycen,1)*(grid%dz/2.0_wp)
     else
        taueq        = sum(grid%rhokap(:,nycen,nzcen))*(grid%dx/2.0_wp)
        taupole      = sum(grid%rhokap(nxcen,nycen,:))*(grid%dz/2.0_wp)
     endif
     dens        = sum(grid%rhokap)/count(grid%rhokap > 0.0)
     par%tauhomo = dens*grid%zmax
     par%taumax  = taupole
  endif
  !call MPI_BCAST(par%taumax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  !--- print out important parameters.
  if (mpar%p_rank == 0) then
     write(6,'(a,es16.3)') 'tau_equat (dust) : ',taueq
     write(6,'(a,es16.3)') 'tau_pole  (dust) : ',taupole
     write(6,'(a,es16.3)') 'tau_homo  (dust) : ',par%tauhomo
  endif
  end subroutine grid_create
  !------------------------------------
  subroutine grid_destroy(grid)
  use define
  use mpi
  implicit none
  type(grid_type), intent(inout) :: grid
  integer :: ierr

  call destroy_shared_mem_all()
  end subroutine grid_destroy
  !=================================================
end module grid_mod
