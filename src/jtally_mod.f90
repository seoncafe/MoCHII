! MoCHII: copied from MoCafe_v2.00/src/jtally_mod.f90 (2026-07-11)
module jtally_mod
!--- Mean-intensity tally J_lambda(cell).
!--- Lucy (1999) pathlength estimator: every actual photon path segment
!--- contributes wgt*dl to its (wavelength-bin, cell) slot; the mean
!--- intensity follows as J_lambda = E_packet * Sum(wgt*dl) / (4 pi V dlam).
!---
!--- The forced first scattering biases the free-path distribution of the
!--- nscatt = 0 flight, so that flight is tallied ANALYTICALLY instead:
!--- raytrace_to_edge_car already walks the full ray to the grid edge, and
!--- (with jt_first = .true.) accumulates the exact expectation
!---    wgt * Int exp(-tau_lambda(l)) dl
!---      = wgt * (exp(-s*tau_in) - exp(-s*tau_out)) / (rhokap*s)
!--- per cell (zero-variance direct component).  Flights with nscatt >= 1
!--- sample the standard exponential free path and are tallied along the
!--- actual walked segments in raytrace_to_tau_car (unbiased).
!---
!--- The tally array is private to each MPI rank and MPI_REDUCEd once at the
!--- end of the run.  SED-mode scope: Cartesian grid.
  use define
  implicit none
  public

  logical :: jt_on    = .false.  ! master switch (par%save_jlam, set in jtally_setup)
  logical :: jt_first = .false.  ! analytic first-flight tally inside the raytrace_to_edge routines
  !--- (nlambda, ncell): Sum(Lpacket*wgt*dl).  ncell = nx*ny*nz ('car') or
  !--- nleaf ('amr'); the linear cell id is defined in cellinfo_mod.
  real(kind=wp), pointer :: jt_sum(:,:) => null()
  integer :: jt_ncell = 0
  real(kind=wp) :: jt_eabs = 0.0_wp  ! independent absorbed-energy counter: Sum Lpacket*wgt*(1-albedo)

  !--- MoCHII: ionizing-band J tally (nnu_ion, nleaf) alongside the SED tally.
  !--- Same Lucy pathlength estimator; filled by raytrace_ion_to_edge_amr
  !--- (analytic first flight; no scattered ionizing flights without dust).
  logical :: jt_ion_on = .false.
  real(kind=wp), pointer :: jt_ion(:,:) => null()
  integer :: jt_ion_nleaf = 0

  !--- MoCHII plane-parallel slab: emergent-intensity boundary tally.  A packet
  !--- that escapes a z-face carries its surviving luminosity into (mu-bin, face)
  !--- with mu = |kz|; face 1 = top (+z, kz>0), face 2 = bottom (-z, kz<0).
  !--- slab_Iesc is the rank-local escaping energy [erg/s] per bin (ALLREDUCEd
  !--- at output).  slab_Labs / slab_Lin track the absorbed / entering totals for
  !--- the energy budget.
  logical :: slab_tally_on = .false.
  integer :: slab_nmu      = 0
  real(kind=wp), allocatable :: slab_Iesc(:,:)     ! (nmu, 2)
  real(kind=wp) :: slab_Lin = 0.0_wp
  !--- which z-faces are illuminated (top = +z, bottom = -z).  The
  !--- reflection/transmission split of the escaping luminosity is only
  !--- defined when EXACTLY ONE face is lit (escape at the lit face = reflected,
  !--- at the opposite face = transmitted).
  logical :: slab_lit_top = .false., slab_lit_bot = .false.

contains

  !---------------------------------------------------------------
  ! Allocate the slab boundary tally (nmu uniform mu-bins in (0,1]).
  ! lit = [top, bottom] illuminated flags (drives the refl/tran labelling).
  subroutine slab_tally_setup(nmu, Lin, lit)
  implicit none
  integer,       intent(in) :: nmu
  real(kind=wp), intent(in) :: Lin
  logical,       intent(in) :: lit(2)
  if (allocated(slab_Iesc)) deallocate(slab_Iesc)
  slab_nmu = nmu
  allocate(slab_Iesc(nmu, 2));  slab_Iesc = 0.0_wp
  slab_Lin = Lin
  slab_lit_top = lit(1);  slab_lit_bot = lit(2)
  slab_tally_on = .true.
  end subroutine slab_tally_setup

  !---------------------------------------------------------------
  ! Add an escaping packet: mu = |kz|, face by sign(kz).  ergs = surviving
  ! luminosity (expo * wgt * Lpacket).
  subroutine slab_bnd_add(kz, ergs)
  implicit none
  real(kind=wp), intent(in) :: kz, ergs
  integer :: ib, iface
  real(kind=wp) :: mu
  if (.not. slab_tally_on) return
  mu = abs(kz)
  ib = min(max(int(mu*slab_nmu) + 1, 1), slab_nmu)
  iface = merge(1, 2, kz > 0.0_wp)      ! +z escape = top, -z = bottom
  slab_Iesc(ib, iface) = slab_Iesc(ib, iface) + ergs
  end subroutine slab_bnd_add

  !---------------------------------------------------------------
  subroutine slab_tally_reduce()
  use mpi
  implicit none
  integer :: ierr
  if (.not. slab_tally_on) return
  call MPI_ALLREDUCE(MPI_IN_PLACE, slab_Iesc, size(slab_Iesc), &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine slab_tally_reduce

  !---------------------------------------------------------------
  ! Write the emergent intensity I(mu) at the two z-boundaries plus the
  ! reflected/transmitted/absorbed energy budget (all normalized by the
  ! entering luminosity).  I(mu) [erg/s/cm^2/sr] uses F = 2 pi int I(mu) mu dmu:
  !   I(mu, face) = Iesc(bin, face) / (A_zface * 2 pi * mu * dmu).
  subroutine slab_write_Imu(fname)
  use octree_mod, only : amr_grid
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=wp) :: azface, dmu, mu, Lesc_top, Lesc_bot, Labs, Lrefl, Ltran
  integer :: ib, u, nlit
  if (.not. slab_tally_on) return
  azface   = (amr_grid%xrange*par%distance2cm)*(amr_grid%yrange*par%distance2cm)
  dmu      = 1.0_wp/real(slab_nmu, wp)
  !--- always-valid base quantities: the raw luminosity escaping each z-face.
  Lesc_top = sum(slab_Iesc(:,1))      ! escaping the top    (+z)
  Lesc_bot = sum(slab_Iesc(:,2))      ! escaping the bottom (-z)
  Labs     = slab_Lin - Lesc_top - Lesc_bot
  nlit     = merge(1,0,slab_lit_top) + merge(1,0,slab_lit_bot)

  open(newunit=u, file=trim(fname), status='replace', action='write')
  write(u,'(a)')        '# MoCHII plane-parallel slab: emergent intensity I(mu) at the z-boundaries.'
  write(u,'(a,es14.6)') '# L_in           [erg/s] = ', slab_Lin
  write(u,'(a,es14.6,a,f8.5)') '# L_escape_top   [erg/s] = ', Lesc_top, '   fraction = ', Lesc_top/slab_Lin
  write(u,'(a,es14.6,a,f8.5)') '# L_escape_bottom[erg/s] = ', Lesc_bot, '   fraction = ', Lesc_bot/slab_Lin
  write(u,'(a,es14.6,a,f8.5)') '# L_absorbed     [erg/s] = ', Labs,     '   fraction = ', Labs/slab_Lin

  !--- reflection/transmission split is only meaningful with a single lit face.
  if (nlit == 1) then
     if (slab_lit_top) then
        write(u,'(a)') '# illuminated face = top: escape at top = reflected, at bottom = transmitted.'
        Lrefl = Lesc_top;  Ltran = Lesc_bot
     else
        write(u,'(a)') '# illuminated face = bottom: escape at bottom = reflected, at top = transmitted.'
        Lrefl = Lesc_bot;  Ltran = Lesc_top
     end if
     write(u,'(a,es14.6,a,f8.5)') '# L_reflected    [erg/s] = ', Lrefl, '   fraction = ', Lrefl/slab_Lin
     write(u,'(a,es14.6,a,f8.5)') '# L_transmitted  [erg/s] = ', Ltran, '   fraction = ', Ltran/slab_Lin
  else
     write(u,'(a)') '# both (or neither) faces illuminated: the top and bottom escape each mix'
     write(u,'(a)') '# reflection and transmission and are not separable without an incident-face tag.'
  end if

  write(u,'(a,es14.6)') '# A_zface        [cm^2]  = ', azface
  write(u,'(a)')        '#   mu        I_top[erg/s/cm^2/sr]   I_bot[erg/s/cm^2/sr]'
  do ib = 1, slab_nmu
     mu = (real(ib,wp) - 0.5_wp)*dmu
     write(u,'(f8.5,2es22.8)') mu, &
        slab_Iesc(ib,1)/(azface*twopi*mu*dmu), &
        slab_Iesc(ib,2)/(azface*twopi*mu*dmu)
  end do
  close(u)
  end subroutine slab_write_Imu

  !---------------------------------------------------------------
  subroutine jtally_setup(grid)
  use sed_mod,      only : sed_nlam
  use memory_mod,   only : create_mem
  use cellinfo_mod, only : ncell_total
  implicit none
  type(grid_type), intent(in) :: grid
  real(kind=wp) :: mem_gb

  jt_ncell = ncell_total(grid)
  mem_gb = real(sed_nlam,wp)*jt_ncell*8.0_wp/1024.0_wp**3
  if (mpar%p_rank == 0) then
     write(*,'(a,i0,a,f8.3,a)') 'J_lambda tally: ', jt_ncell, ' cells, ', mem_gb, ' GB per MPI rank'
     if (mem_gb > 4.0_wp) write(*,'(a)') &
        'WARNING: J_lambda tally is large; consider fewer wavelength bins or cells.'
  endif
  call create_mem(jt_sum, [sed_nlam, jt_ncell])
  jt_sum(:,:) = 0.0_wp
  jt_eabs = 0.0_wp
  jt_on   = .true.
  end subroutine jtally_setup

  !---------------------------------------------------------------
  subroutine jtally_reduce()
  use mpi
  implicit none
  integer :: ierr, ic, nchunk, i0, n
  !--- ALLREDUCE (not reduce-to-0) so every rank holds the full tally: the
  !--- dust-emission stage distributes cells across ranks and reads jt_sum
  !--- locally.  Reduce in chunks of cells to keep each MPI count within int32.
  if (.not. jt_on) return
  nchunk = max(1, 100000000/max(size(jt_sum,1),1))   ! ~1e8 elements per call
  i0 = 1
  do while (i0 <= jt_ncell)
     n = min(nchunk, jt_ncell-i0+1)
     call MPI_ALLREDUCE(MPI_IN_PLACE, jt_sum(:,i0:i0+n-1), int(size(jt_sum,1)*n), &
                        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     i0 = i0 + n
  enddo
  call MPI_ALLREDUCE(MPI_IN_PLACE, jt_eabs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end subroutine jtally_reduce

  !---------------------------------------------------------------
  !--- MoCHII: ionizing-band tally setup/reduce.  The conversion to J_nu and
  !--- the rate integrals live in gas_rates_mod (they need the gas state).
  subroutine jtally_ion_setup(nleaf)
  use memory_mod,   only : create_mem
  use ion_band_mod, only : nnu_band
  implicit none
  integer, intent(in) :: nleaf
  real(kind=wp) :: mem_gb
  jt_ion_nleaf = nleaf
  mem_gb = real(nnu_band,wp)*nleaf*8.0_wp/1024.0_wp**3
  if (mpar%p_rank == 0) write(*,'(a,i0,a,f8.3,a)') &
     ' ION: J tally: ', nleaf, ' leaves, ', mem_gb, ' GB per MPI rank'
  call create_mem(jt_ion, [nnu_band, nleaf])
  jt_ion(:,:) = 0.0_wp
  jt_ion_on   = .true.
  end subroutine jtally_ion_setup

  !---------------------------------------------------------------
  !--- MoCHII: resize after octree re-refinement.
  subroutine jtally_ion_resize(nleaf)
  implicit none
  integer, intent(in) :: nleaf
  if (associated(jt_ion)) deallocate(jt_ion)
  jt_ion => null()
  call jtally_ion_setup(nleaf)
  end subroutine jtally_ion_resize

  !---------------------------------------------------------------
  subroutine jtally_ion_reduce()
  use mpi
  implicit none
  integer :: ierr, nchunk, i0, n
  if (.not. jt_ion_on) return
  nchunk = max(1, 100000000/max(size(jt_ion,1),1))
  i0 = 1
  do while (i0 <= jt_ion_nleaf)
     n = min(nchunk, jt_ion_nleaf-i0+1)
     call MPI_ALLREDUCE(MPI_IN_PLACE, jt_ion(:,i0:i0+n-1), int(size(jt_ion,1)*n), &
                        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     i0 = i0 + n
  enddo
  end subroutine jtally_ion_reduce

  !---------------------------------------------------------------
  !--- convert the tally to J_lambda and write '<base>_jlam.<ext>'.
  !--- J_lambda units: [luminosity unit] / dist_cm^2 / um / sr, where
  !--- dist_cm = par%distance2cm (1 when no physical unit is given).
  !--- Also writes the wavelength-integrated J_bol(x,y,z) and, for the
  !--- energy-conservation check, prints the absorbed luminosity from the
  !--- tally (A) against the independent event counter (B).
  subroutine jtally_write(grid)
  use sed_mod,      only : sed_nlam, sed_wave, sed_dwave, sed_sext, sed_albedo, sed_cext_ref
  use cellinfo_mod, only : cell_rhokap, cell_volume, cell_center, car_ijk
  use iofile_mod
  use utility,      only : get_base_name
  implicit none
  type(grid_type), intent(in) :: grid
  type(io_file_type) :: file
  character(len=192) :: filename
  real(kind=wp), allocatable :: jbol(:), jcube(:,:,:,:), jbolc(:,:,:), leafxyz(:,:)
  real(kind=wp) :: fac, eabs_A, eabs_B, rhk, cx, cy, cz
  integer :: status, il, ic, i, j, k
  logical :: is_amr

  if (mpar%p_rank /= 0 .or. .not. associated(jt_sum)) return
  is_amr = trim(par%grid_type) == 'amr' .or. trim(par%grid_type) == 'car'

  !--- energy-conservation check (jt_sum carries Lpacket; A = tally, B = counter).
  eabs_A = 0.0_wp
  do ic = 1, jt_ncell
     rhk = cell_rhokap(grid, ic)
     if (rhk > 0.0_wp) eabs_A = eabs_A + rhk*sum(jt_sum(:,ic)*sed_sext(:)*(1.0_wp - sed_albedo(:)))
  enddo
  eabs_B = jt_eabs
  write(*,'(a)')        '--- J_lambda tally: energy conservation check ---'
  write(*,'(a,es14.6)') 'absorbed L (pathlength tally, A): ', eabs_A
  write(*,'(a,es14.6)') 'absorbed L (event counter,    B): ', eabs_B
  if (eabs_B > 0.0_wp) write(*,'(a,f10.6)') 'ratio A/B                       : ', eabs_A/eabs_B

  !--- convert Sum(Lpacket*wgt*dl) -> J_lambda [erg/s/cm^2/sr/um] per cell.
  do ic = 1, jt_ncell
     fac = 1.0_wp/(fourpi*cell_volume(grid,ic)*par%distance2cm**2)
     do il = 1, sed_nlam
        jt_sum(il,ic) = jt_sum(il,ic)*(fac/sed_dwave(il))
     enddo
  enddo
  allocate(jbol(jt_ncell))
  do ic = 1, jt_ncell
     jbol(ic) = sum(jt_sum(:,ic)*sed_dwave(:))
  enddo

  status = 0
  filename = trim(get_base_name(par%out_file))//'_jlam'//trim(io_file_extension(par%file_format))
  call io_open_new(file, trim(filename), status)

  if (is_amr) then
     !--- AMR: write the J_lambda(nlam,nleaf), J_bol, and leaf x,y,z.
     call io_append_image(file, jt_sum, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','J_lambda','J(lambda,leaf) mean intensity (AMR)',status)
     call io_put_keyword(file,'J_UNIT','luminosity/dist_cm^2/um/sr','J_lambda unit',status)
     call write_jlam_keys(file, sed_cext_ref, eabs_A, eabs_B)
     call io_append_image(file, jbol, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','J_bol','wavelength-integrated J per leaf',status)
     allocate(leafxyz(jt_ncell,3))
     do ic = 1, jt_ncell
        call cell_center(grid, ic, cx, cy, cz)
        leafxyz(ic,1) = cx;  leafxyz(ic,2) = cy;  leafxyz(ic,3) = cz
     enddo
     call io_append_image(file, leafxyz, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','LeafXYZ','leaf center x,y,z (code units)',status)
     deallocate(leafxyz)
  else
     !--- Cartesian: reshape to the (nlam,nx,ny,nz) cube for backward-compatible output.
     allocate(jcube(sed_nlam,grid%nx,grid%ny,grid%nz), jbolc(grid%nx,grid%ny,grid%nz))
     do ic = 1, jt_ncell
        call car_ijk(grid, ic, i, j, k)
        jcube(:,i,j,k) = jt_sum(:,ic);  jbolc(i,j,k) = jbol(ic)
     enddo
     call io_append_image(file, jcube, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','J_lambda','J(lambda,x,y,z) mean intensity',status)
     call io_put_keyword(file,'J_UNIT','luminosity/dist_cm^2/um/sr','J_lambda unit',status)
     call write_jlam_keys(file, sed_cext_ref, eabs_A, eabs_B)
     call io_append_image(file, jbolc, status, bitpix=-64)
     call io_put_keyword(file,'EXTNAME','J_bol','wavelength-integrated J(x,y,z)',status)
     call io_put_keyword(file,'J_UNIT','luminosity/dist_cm^2/sr','J_bol unit',status)
     deallocate(jcube, jbolc)
  endif
  call io_append_image(file, sed_wave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Wavelength','bin centers [um]',status)
  call io_append_image(file, sed_dwave, status, bitpix=-64)
  call io_put_keyword(file,'EXTNAME','Dwavelength','bin widths [um]',status)
  call io_close(file, status)
  write(*,'(2a)') 'J_lambda written to: ', trim(filename)
  deallocate(jbol)
  end subroutine jtally_write

  !---------------------------------------------------------------
  subroutine write_jlam_keys(file, cext_ref, eabs_A, eabs_B)
  use sed_mod,   only : sed_nlam
  use iofile_mod
  implicit none
  type(io_file_type), intent(inout) :: file
  real(kind=wp),      intent(in)    :: cext_ref, eabs_A, eabs_B
  integer :: status
  status = 0
  call io_put_keyword(file,'SED_NLAM', sed_nlam,       'number of wavelength bins',     status)
  call io_put_keyword(file,'SED_LREF', par%lambda_ref, 'reference wavelength [um]',     status)
  call io_put_keyword(file,'SED_CREF', cext_ref,       'C_ext/H at lambda_ref [cm^2/H]',status)
  call io_put_keyword(file,'TOT_LUM',  par%luminosity, 'total luminosity',              status)
  call io_put_keyword(file,'DIST_CM',  par%distance2cm,'distance unit (cm)',            status)
  call io_put_keyword(file,'nphotons', par%no_photons, 'number of photons',             status)
  call io_put_keyword(file,'taumax',   par%taumax,     'tau_max at lambda_ref',         status)
  call io_put_keyword(file,'EABS_A',   eabs_A, 'absorbed L (pathlength tally)',         status)
  call io_put_keyword(file,'EABS_B',   eabs_B, 'absorbed L (event counter)',            status)
  end subroutine write_jlam_keys

end module jtally_mod
