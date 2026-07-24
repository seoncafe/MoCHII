! MoCHII: copied from MoCafe_v2.00/src/raytrace_amr.f90 (2026-07-11)
module raytrace_amr_mod
!---------------------------------------------------------------------------
! AMR octree raytrace (dust-only).
!
!   Positions and path lengths are in CODE UNITS; amr_grid%rhokap(il) is the
!   grey dust opacity per code unit, so tau = rhokap*ds is dimensionless.
!   Octant traversal uses the precomputed neighbor table (amr_next_leaf,
!   O(1) per face crossing).  Face index iface: 1=+x 2=-x 3=+y 4=-y 5=+z 6=-z.
!
! Bound to the raytrace_to_tau / raytrace_to_edge procedure pointers for
! grid_type='amr'; scattering and peeling-off are SHARED with the Cartesian
! dust path (dust scattering is leaf-independent).  Boundary handling honors
! par%xy_periodic / xy_symmetry /
! xyz_symmetry (no frequency shift for grey dust); a space-filling octree has
! no refinement gaps, so amr_next_leaf alone suffices for the transport walk.
!---------------------------------------------------------------------------
  use octree_mod
  use jtally_mod, only : jt_on, jt_first, jt_sum, jt_ion, &
                         slab_tally_on, slab_bnd_add
  implicit none
  private

  public :: raytrace_to_tau_amr
  public :: raytrace_to_edge_amr
  public :: raytrace_ion_to_edge_amr
  public :: raytrace_ion_tau_only_amr
  public :: transport_ion_packet
  public :: slab_walk_report

  real(wp), parameter :: tau_huge = 745.2_wp  ! exp(-tau_huge) ~ 0 in double

  !--- Safety net for the xy-periodic slab walks.  A near-grazing ray (kz -> 0)
  !--- through a nearly transparent medium can wrap the x/y faces indefinitely
  !--- without ever reaching a z-face or accumulating tau.  The source-side
  !--- guards (beam theta strictly < 90 deg; isotropic mu floored away from 0)
  !--- prevent kz = 0 exactly, so any legitimate walk terminates; MAX_WALK_STEPS
  !--- caps the residual pathology far above any physical traversal length (a
  !--- normal slab crosses O(nz/mu) cells).  A capped packet is terminated
  !--- safely and counted in slab_walk_kills (rank-local, reported once).
  integer, parameter :: MAX_WALK_STEPS = 20000000
  integer :: slab_walk_kills = 0

contains

  !=========================================================================
  ! Reduce and report the near-grazing max-step guard hits once.  Prints
  ! nothing when the guard never fired (the normal case).
  !=========================================================================
  subroutine slab_walk_report()
    use define
    use mpi
    implicit none
    integer :: ierr, total
    call MPI_ALLREDUCE(slab_walk_kills, total, 1, MPI_INTEGER, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    if (mpar%p_rank == 0 .and. total > 0) write(*,'(a,i0,a)') &
       ' WARNING: ', total, ' packet(s) hit the slab max-step guard '// &
       '(near-grazing ray in the periodic slab) and were terminated early.'
    slab_walk_kills = 0
  end subroutine slab_walk_report

  !=========================================================================
  subroutine raytrace_to_tau_amr(photon, grid, tau_in)
    use define
    implicit none
    type(photon_type), intent(inout) :: photon
    type(grid_type),   intent(inout) :: grid   ! unused for AMR; interface conformance
    real(wp),          intent(in)    :: tau_in

    integer  :: il, il_new, icell, iface
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: tau, t_exit, d_step, rhokap
    !--- Lucy pathlength J tally for unforced flights (nscatt > 0); the forced
    !--- first flight is tallied analytically in raytrace_to_edge_amr.
    logical  :: do_tally

    do_tally = jt_on .and. photon%nscatt > 0

    x  = photon%x;    y  = photon%y;    z  = photon%z
    kx = photon%kx;   ky = photon%ky;   kz = photon%kz
    il = photon%icell_amr

    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) then
        photon%inside = .false.
        return
      end if
    end if

    tau = 0.0_wp
    do while (photon%inside)
      icell = leaf_cell(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      rhokap = amr_grid%rhokap(il)

      if (tau + t_exit * rhokap >= tau_in) then
        if (rhokap > 0.0_wp) then
          d_step = (tau_in - tau) / rhokap
        else
          d_step = t_exit
        end if
        if (do_tally) jt_sum(photon%il,il) = jt_sum(photon%il,il) + photon%wgt*photon%Lpacket*d_step
        x = x + d_step * kx
        y = y + d_step * ky
        z = z + d_step * kz
        tau = tau_in
        exit
      end if

      if (do_tally) jt_sum(photon%il,il) = jt_sum(photon%il,il) + photon%wgt*photon%Lpacket*t_exit
      tau = tau + t_exit * rhokap
      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz

      il_new = amr_next_leaf(icell, iface, x, y, z)

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx;  il_new = il
          case (4);  ky = -ky;  il_new = il
          case (6);  kz = -kz;  il_new = il
          end select
        end if
      end if

      if (il_new <= 0) then
        photon%inside = .false.
        exit
      end if
      il = il_new
    end do

    photon%x  = x;   photon%y  = y;   photon%z  = z
    photon%kx = kx;  photon%ky = ky;  photon%kz = kz
    photon%icell_amr = il
  end subroutine raytrace_to_tau_amr

  !=========================================================================
  ! Integrate the dust optical depth from photon0 to the box edge.  photon0
  ! is read-only.  Bound to raytrace_to_edge for the peeling-off estimator.
  !=========================================================================
  subroutine raytrace_to_edge_amr(photon0, grid, tau)
    use define
    implicit none
    type(photon_type), intent(in)  :: photon0
    type(grid_type),   intent(in)  :: grid
    real(wp),          intent(out) :: tau

    integer  :: il, il_new, icell, iface
    real(wp) :: x, y, z, kx, ky, kz, t_exit, rhokap
    !--- analytic first-flight J tally (jt_first): running exp(-s*tau).
    real(wp) :: jt_expo, alpha, expo_out

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il = photon0%icell_amr
    jt_expo = 1.0_wp

    tau = 0.0_wp
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) return
    end if

    do
      icell  = leaf_cell(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      rhokap = amr_grid%rhokap(il)
      if (jt_first) then
        alpha = rhokap*photon0%s_ext
        if (alpha*t_exit > 0.0_wp) then
          expo_out = jt_expo*exp(-alpha*t_exit)
          jt_sum(photon0%il,il) = jt_sum(photon0%il,il) + photon0%wgt*photon0%Lpacket*(jt_expo-expo_out)/alpha
          jt_expo = expo_out
        else
          jt_sum(photon0%il,il) = jt_sum(photon0%il,il) + photon0%wgt*photon0%Lpacket*jt_expo*t_exit
        end if
      end if
      tau = tau + t_exit * rhokap
      if (tau >= tau_huge) return

      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new = amr_next_leaf(icell, iface, x, y, z)

      if (il_new <= 0) then
        if (par%xy_periodic) then
          select case (iface)
          case (1); x = x - amr_grid%xrange
          case (2); x = x + amr_grid%xrange
          case (3); y = y - amr_grid%yrange
          case (4); y = y + amr_grid%yrange
          end select
          if (iface <= 4) il_new = amr_find_leaf(x, y, z)
        else if (par%xyz_symmetry) then
          select case (iface)
          case (2);  kx = -kx;  il_new = il
          case (4);  ky = -ky;  il_new = il
          case (6);  kz = -kz;  il_new = il
          end select
        end if
      end if

      if (il_new <= 0) return
      il = il_new
    end do
  end subroutine raytrace_to_edge_amr

  !=========================================================================
  ! MoCHII: ionizing-band analytic edge walk.  Walks the packet
  ! from its birth point to the domain edge, accumulating the exact
  ! expectation of the pathlength tally in each leaf,
  !     jt_ion(inu, il) += wgt * Lpacket * (e^{-tau_in} - e^{-tau_out}) / kap
  ! with kap = kap_ion(inu, il) (per code length; gas + optional dust).
  ! With no ionizing-band scattering, this zero-variance direct
  ! component is the complete tally; the packet needs no interaction
  ! sampling.  For the plane-parallel slab boundary condition (PP1) the walk
  ! wraps the x/y faces xy-periodically (slab_wrap_xy on the octree, index
  ! modulo on the car/DDA path); z faces remain open escapes.
  !=========================================================================
  !=========================================================================
  ! MoCHII: full ionizing-band transport of one packet.
  ! (1) analytic zero-variance direct tally to the edge (always);
  ! (2) when par%ion_dust_scatter: forced first interaction along the ray
  !     (weight w1 = 1 - e^-tau_edge), albedo weight at each interaction
  !     (a_eff = kappa_sca,dust / kappa_ext), Henyey-Greenstein direction
  !     with g(E_b) from the kext table, Russian roulette below w_min;
  !     the scattered flights tally jt_ion along their sampled segments
  !     (the Lucy pathlength estimator, nscatt > 0).
  ! Gas heating/ionization and dust heating derive from the same J tally,
  ! so the energy split is bookkept by construction.
  !=========================================================================
  subroutine transport_ion_packet(photon)
    use define
    use random,          only : rand_number
    use gas_opacity_mod, only : kap_ion, ion_dust_ssca, ion_dust_g
    implicit none
    type(photon_type), intent(inout) :: photon
    real(wp), parameter :: WMIN = 1.0e-4_wp, PSURV = 0.1_wp
    real(wp) :: tau_edge, wgt1, tau, a_eff, gg, cost, sint, phi
    real(wp) :: kx, ky, kz, ux, uy, uz, vx, vy, vz, norm
    integer  :: il

    call raytrace_ion_to_edge_amr(photon, tau_edge)
    if (.not. (par%ion_add_dust .and. par%ion_dust_scatter)) return
    if (tau_edge <= 1.0e-12_wp) return

    !--- forced first interaction (direct J already tallied analytically)
    wgt1 = 1.0_wp - exp(-tau_edge)
    photon%wgt = photon%wgt*wgt1
    tau = -log(1.0_wp - rand_number()*wgt1)
    call raytrace_ion_to_tau_amr(photon, tau)

    do while (photon%inside)
       il = photon%icell_amr
       if (il <= 0) exit
       !--- effective albedo at (bin, leaf): dust scattering / extinction
       a_eff = amr_grid%rhokap(il)*ion_dust_ssca(photon%inu) &
               / max(kap_ion(photon%inu, il), tinest)
       photon%wgt = photon%wgt*a_eff
       !--- peel-off of the scattered contribution (imaging pass only;
       !--- photon%k is still the incident direction here).
       if (associated(ion_peel_scatter_hook)) &
          call ion_peel_scatter_hook(photon)
       if (photon%wgt < WMIN) then
          if (rand_number() > PSURV) exit
          photon%wgt = photon%wgt/PSURV
       end if
       !--- Henyey-Greenstein scattering about the current direction
       gg = ion_dust_g(photon%inu)
       if (abs(gg) > 1.0e-3_wp) then
          cost = (1.0_wp + gg*gg - ((1.0_wp - gg*gg) &
                 /(1.0_wp + gg*(2.0_wp*rand_number() - 1.0_wp)))**2) &
                 /(2.0_wp*gg)
       else
          cost = 2.0_wp*rand_number() - 1.0_wp
       end if
       cost = min(max(cost, -1.0_wp), 1.0_wp)
       sint = sqrt(1.0_wp - cost*cost)
       phi  = twopi*rand_number()
       kx = photon%kx;  ky = photon%ky;  kz = photon%kz
       norm = sqrt(kx*kx + ky*ky + kz*kz)
       if (norm <= 0.0_wp) exit
       kx = kx/norm;  ky = ky/norm;  kz = kz/norm
       !--- orthonormal frame (u, v, k)
       if (abs(kz) < 0.99_wp) then
          ux = -ky;  uy = kx;  uz = 0.0_wp
       else
          ux = 0.0_wp;  uy = -kz;  uz = ky
       end if
       norm = sqrt(ux*ux + uy*uy + uz*uz)
       ux = ux/norm;  uy = uy/norm;  uz = uz/norm
       vx = ky*uz - kz*uy;  vy = kz*ux - kx*uz;  vz = kx*uy - ky*ux
       photon%kx = sint*(cos(phi)*ux + sin(phi)*vx) + cost*kx
       photon%ky = sint*(cos(phi)*uy + sin(phi)*vy) + cost*ky
       photon%kz = sint*(cos(phi)*uz + sin(phi)*vz) + cost*kz
       norm = sqrt(photon%kx**2 + photon%ky**2 + photon%kz**2)
       photon%kx = photon%kx/norm
       photon%ky = photon%ky/norm
       photon%kz = photon%kz/norm
       photon%nscatt = photon%nscatt + 1
       !--- next flight (unforced); segments tally into jt_ion
       tau = -log(max(rand_number(), tinest))
       call raytrace_ion_to_tau_amr(photon, tau)
    end do
  end subroutine transport_ion_packet

  !=========================================================================
  ! Ionizing-band walk to optical depth tau_in (extinction kap_ion);
  ! scattered flights (nscatt > 0) tally jt_ion along their segments.
  !=========================================================================
  !=========================================================================
  ! Robust xy-periodic wrap for the octree walk.  A packet that leaves an
  ! x/y face (amr_next_leaf returned 0) is re-entered at the OPPOSITE face by
  ! setting the crossed coordinate to that face's EXACT boundary value, not by
  ! x -/+ xrange: the latter can land a rounding-epsilon outside [min,max], and
  ! amr_find_leaf then drops the packet (a photon loss that, at nx=1 where every
  ! lateral crossing wraps, over-attenuates the slab ~10x).  z faces escape.
  !=========================================================================
  subroutine slab_wrap_xy(iface, x, y)
    integer,  intent(in)    :: iface
    real(wp), intent(inout) :: x, y
    select case (iface)
    case (1); x = amr_grid%xmin    ! exited +x -> re-enter at -x
    case (2); x = amr_grid%xmax    ! exited -x -> re-enter at +x
    case (3); y = amr_grid%ymin
    case (4); y = amr_grid%ymax
    end select
  end subroutine slab_wrap_xy

  subroutine raytrace_ion_to_tau_amr(photon, tau_in)
    use define
    use gas_opacity_mod, only : kap_ion
    implicit none
    type(photon_type), intent(inout) :: photon
    real(wp),          intent(in)    :: tau_in

    integer  :: il, il_new, icell, iface, inu, nstep
    real(wp) :: x, y, z, kx, ky, kz
    real(wp) :: tau, t_exit, d_step, kap, wl
    logical  :: do_tally

    if (amr_grid%car .and. trim(par%car_walk) /= 'shared') then
       call raytrace_ion_to_tau_car(photon, tau_in)
       return
    end if

    inu = photon%inu
    do_tally = photon%nscatt > 0
    wl = photon%wgt*photon%Lpacket

    x  = photon%x;    y  = photon%y;    z  = photon%z
    kx = photon%kx;   ky = photon%ky;   kz = photon%kz
    il = photon%icell_amr
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) then
        photon%inside = .false.
        return
      end if
    end if

    tau = 0.0_wp
    nstep = 0
    do while (photon%inside)
      nstep = nstep + 1
      if (nstep > MAX_WALK_STEPS) then
        photon%inside = .false.;  slab_walk_kills = slab_walk_kills + 1;  exit
      end if
      icell = leaf_cell(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      kap = kap_ion(inu, il)
      if (tau + t_exit*kap >= tau_in) then
        if (kap > 0.0_wp) then
          d_step = (tau_in - tau)/kap
        else
          d_step = t_exit
        end if
        if (do_tally) jt_ion(inu,il) = jt_ion(inu,il) + wl*d_step
        x = x + d_step*kx;  y = y + d_step*ky;  z = z + d_step*kz
        exit
      end if
      if (do_tally) jt_ion(inu,il) = jt_ion(inu,il) + wl*t_exit
      tau = tau + t_exit*kap
      x = x + t_exit*kx;  y = y + t_exit*ky;  z = z + t_exit*kz
      il_new = amr_next_leaf(icell, iface, x, y, z)
      if (il_new <= 0 .and. par%xy_periodic) then
        call slab_wrap_xy(iface, x, y)
        if (iface <= 4) il_new = amr_find_leaf(x, y, z)
      end if
      if (il_new <= 0) then
        photon%inside = .false.
        !--- a scattered flight escaping a z-face is the diffuse/reflected
        !--- boundary contribution (the direct beam is tallied in the edge walk).
        if (slab_tally_on .and. iface >= 5) call slab_bnd_add(kz, wl)
        exit
      end if
      il = il_new
    end do

    photon%x = x;   photon%y = y;   photon%z = z
    photon%kx = kx; photon%ky = ky; photon%kz = kz
    photon%icell_amr = il
  end subroutine raytrace_ion_to_tau_amr

  subroutine raytrace_ion_to_edge_amr(photon0, tau_edge_out)
    use define
    use gas_opacity_mod, only : kap_ion
    implicit none
    type(photon_type), intent(in) :: photon0
    real(wp), intent(out), optional :: tau_edge_out

    integer  :: il, il_new, icell, iface, inu, nstep
    real(wp) :: x, y, z, kx, ky, kz, t_exit, kap
    real(wp) :: tau, expo, expo_out, wl

    if (amr_grid%car .and. trim(par%car_walk) /= 'shared') then
       call raytrace_ion_to_edge_car(photon0, tau_edge_out)
       return
    end if

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il  = photon0%icell_amr
    inu = photon0%inu
    wl  = photon0%wgt * photon0%Lpacket
    expo = 1.0_wp
    tau  = 0.0_wp
    if (present(tau_edge_out)) tau_edge_out = 0.0_wp

    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) return
    end if

    nstep = 0
    do
      nstep = nstep + 1
      if (nstep > MAX_WALK_STEPS) then
        slab_walk_kills = slab_walk_kills + 1;  exit
      end if
      icell = leaf_cell(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      kap = kap_ion(inu, il)
      if (kap*t_exit > 0.0_wp) then
        expo_out = expo*exp(-kap*t_exit)
        jt_ion(inu,il) = jt_ion(inu,il) + wl*(expo - expo_out)/kap
        expo = expo_out
      else
        jt_ion(inu,il) = jt_ion(inu,il) + wl*expo*t_exit
      end if
      tau = tau + kap*t_exit
      if (tau >= tau_huge) exit

      x = x + t_exit * kx
      y = y + t_exit * ky
      z = z + t_exit * kz
      il_new = amr_next_leaf(icell, iface, x, y, z)
      if (il_new <= 0 .and. par%xy_periodic) then
        call slab_wrap_xy(iface, x, y)
        if (iface <= 4) il_new = amr_find_leaf(x, y, z)
      end if
      if (il_new <= 0) then
        !--- a z-face escape carries the surviving luminosity to the boundary.
        if (slab_tally_on .and. iface >= 5) call slab_bnd_add(kz, expo*wl)
        exit
      end if
      il = il_new
    end do
    if (present(tau_edge_out)) tau_edge_out = tau
  end subroutine raytrace_ion_to_edge_amr

  !=========================================================================
  ! Incremental Amanatides-Woo walks for the Cartesian (car) grid
  ! (par%car_walk = 'dda'): per-ray tMax/tDelta per axis, the inner loop
  ! advances by comparisons and additions only — no cell-geometry reads.
  ! The tally expressions are identical to the shared walks; only the
  ! step computation differs (results agree statistically, not bitwise).
  !=========================================================================
  subroutine dda_init(x, y, z, kx, ky, kz, ix, iy, iz, stp, tmax, tdel, ok)
    use define
    implicit none
    real(wp), intent(in)  :: x, y, z, kx, ky, kz
    integer,  intent(out) :: ix, iy, iz, stp(3)
    real(wp), intent(out) :: tmax(3), tdel(3)
    logical,  intent(out) :: ok
    real(wp) :: d, pos(3), kk(3), org(3)
    integer  :: n(3), a

    d = amr_grid%dcell
    pos = [x, y, z];  kk = [kx, ky, kz]
    org = [amr_grid%xmin, amr_grid%ymin, amr_grid%zmin]
    n   = [amr_grid%nx, amr_grid%ny, amr_grid%nz]
    ok = .true.
    ix = min(max(int((x - org(1))/d), 0), n(1)-1)
    iy = min(max(int((y - org(2))/d), 0), n(2)-1)
    iz = min(max(int((z - org(3))/d), 0), n(3)-1)
    block
      integer :: idx(3)
      idx = [ix, iy, iz]
      do a = 1, 3
         if (kk(a) > 0.0_wp) then
            stp(a)  = 1
            tmax(a) = (org(a) + real(idx(a)+1, wp)*d - pos(a))/kk(a)
            tdel(a) = d/kk(a)
         else if (kk(a) < 0.0_wp) then
            stp(a)  = -1
            tmax(a) = (org(a) + real(idx(a), wp)*d - pos(a))/kk(a)
            tdel(a) = -d/kk(a)
         else
            stp(a)  = 0
            tmax(a) = hugest
            tdel(a) = hugest
         end if
      end do
    end block
  end subroutine dda_init

  subroutine raytrace_ion_to_edge_car(photon0, tau_edge_out)
    use define
    use gas_opacity_mod, only : kap_ion
    implicit none
    type(photon_type), intent(in) :: photon0
    real(wp), intent(out), optional :: tau_edge_out
    integer  :: ix, iy, iz, stp(3), a, il, inu, nx, ny, nz, nstep
    real(wp) :: tmax(3), tdel(3), t_cur, seg, kap, tau
    real(wp) :: wl, expo, expo_out
    logical  :: ok

    call dda_init(photon0%x, photon0%y, photon0%z, photon0%kx, &
                  photon0%ky, photon0%kz, ix, iy, iz, stp, tmax, tdel, ok)
    nx = amr_grid%nx;  ny = amr_grid%ny;  nz = amr_grid%nz
    inu = photon0%inu
    wl  = photon0%wgt*photon0%Lpacket
    expo = 1.0_wp;  tau = 0.0_wp;  t_cur = 0.0_wp
    if (present(tau_edge_out)) tau_edge_out = 0.0_wp

    nstep = 0
    do
      nstep = nstep + 1
      if (nstep > MAX_WALK_STEPS) then
        slab_walk_kills = slab_walk_kills + 1;  exit
      end if
      il  = 1 + ix + nx*(iy + ny*iz)
      a   = minloc(tmax, dim=1)
      seg = tmax(a) - t_cur
      kap = kap_ion(inu, il)
      if (kap*seg > 0.0_wp) then
        expo_out = expo*exp(-kap*seg)
        jt_ion(inu,il) = jt_ion(inu,il) + wl*(expo - expo_out)/kap
        expo = expo_out
      else
        jt_ion(inu,il) = jt_ion(inu,il) + wl*expo*seg
      end if
      tau = tau + kap*seg
      if (tau >= tau_huge) exit
      t_cur = tmax(a)
      tmax(a) = tmax(a) + tdel(a)
      select case (a)
      case (1); ix = ix + stp(1)
                if (ix < 0 .or. ix >= nx) then
                  !--- slab: wrap x index (uniform grid: tmax/tdel unchanged).
                  if (par%xy_periodic) then; ix = modulo(ix, nx); else; exit; end if
                end if
      case (2); iy = iy + stp(2)
                if (iy < 0 .or. iy >= ny) then
                  if (par%xy_periodic) then; iy = modulo(iy, ny); else; exit; end if
                end if
      case (3); iz = iz + stp(3)
                if (iz < 0 .or. iz >= nz) then
                  !--- z-face escape: surviving luminosity to the boundary.
                  if (slab_tally_on) call slab_bnd_add(photon0%kz, expo*wl)
                  exit
                end if
      end select
    end do
    if (present(tau_edge_out)) tau_edge_out = tau
  end subroutine raytrace_ion_to_edge_car

  subroutine raytrace_ion_to_tau_car(photon, tau_in)
    use define
    use gas_opacity_mod, only : kap_ion
    implicit none
    type(photon_type), intent(inout) :: photon
    real(wp),          intent(in)    :: tau_in
    integer  :: ix, iy, iz, stp(3), a, il, inu, nx, ny, nz, nstep
    real(wp) :: tmax(3), tdel(3), t_cur, seg, kap, tau, wl, d_step
    logical  :: ok, do_tally

    call dda_init(photon%x, photon%y, photon%z, photon%kx, photon%ky, &
                  photon%kz, ix, iy, iz, stp, tmax, tdel, ok)
    nx = amr_grid%nx;  ny = amr_grid%ny;  nz = amr_grid%nz
    inu = photon%inu
    do_tally = photon%nscatt > 0
    wl = photon%wgt*photon%Lpacket
    tau = 0.0_wp;  t_cur = 0.0_wp

    nstep = 0
    do
      nstep = nstep + 1
      if (nstep > MAX_WALK_STEPS) then
        photon%inside = .false.;  slab_walk_kills = slab_walk_kills + 1;  exit
      end if
      il  = 1 + ix + nx*(iy + ny*iz)
      a   = minloc(tmax, dim=1)
      seg = tmax(a) - t_cur
      kap = kap_ion(inu, il)
      if (tau + seg*kap >= tau_in) then
        if (kap > 0.0_wp) then
          d_step = (tau_in - tau)/kap
        else
          d_step = seg
        end if
        if (do_tally) jt_ion(inu,il) = jt_ion(inu,il) + wl*d_step
        t_cur = t_cur + d_step
        exit
      end if
      if (do_tally) jt_ion(inu,il) = jt_ion(inu,il) + wl*seg
      tau = tau + seg*kap
      t_cur = tmax(a)
      tmax(a) = tmax(a) + tdel(a)
      select case (a)
      case (1); ix = ix + stp(1)
                if (ix < 0 .or. ix >= nx) then
                  !--- slab: wrap x index (uniform grid: tmax/tdel unchanged).
                  if (par%xy_periodic) then
                    ix = modulo(ix, nx)
                  else
                    photon%inside = .false.;  exit
                  end if
                end if
      case (2); iy = iy + stp(2)
                if (iy < 0 .or. iy >= ny) then
                  if (par%xy_periodic) then
                    iy = modulo(iy, ny)
                  else
                    photon%inside = .false.;  exit
                  end if
                end if
      case (3); iz = iz + stp(3); if (iz < 0 .or. iz >= nz) then
                   photon%inside = .false.
                   !--- scattered flight escaping a z-face: diffuse/reflected
                   !--- boundary contribution (direct beam is in the edge walk).
                   if (slab_tally_on) call slab_bnd_add(photon%kz, wl)
                   exit
                end if
      end select
    end do

    !--- final position from the total path length travelled.
    photon%x = photon%x + t_cur*photon%kx
    photon%y = photon%y + t_cur*photon%ky
    photon%z = photon%z + t_cur*photon%kz
    !--- the accumulated x/y are un-wrapped; fold them back into the box so
    !--- they match the wrapped index cell (icell_amr below).
    if (par%xy_periodic) then
       photon%x = photon%x - floor((photon%x - amr_grid%xmin)/amr_grid%xrange)*amr_grid%xrange
       photon%y = photon%y - floor((photon%y - amr_grid%ymin)/amr_grid%yrange)*amr_grid%yrange
    end if
    if (photon%inside) then
       photon%icell_amr = 1 + ix + nx*(iy + ny*iz)
    end if
  end subroutine raytrace_ion_to_tau_car

  subroutine raytrace_ion_tau_only_car(photon0, tau_out)
    use define
    use gas_opacity_mod, only : kap_ion
    implicit none
    type(photon_type), intent(in)  :: photon0
    real(wp),          intent(out) :: tau_out
    integer  :: ix, iy, iz, stp(3), a, il, inu, nx, ny, nz, nstep
    real(wp) :: tmax(3), tdel(3), t_cur, seg, tau
    logical  :: ok

    call dda_init(photon0%x, photon0%y, photon0%z, photon0%kx, &
                  photon0%ky, photon0%kz, ix, iy, iz, stp, tmax, tdel, ok)
    nx = amr_grid%nx;  ny = amr_grid%ny;  nz = amr_grid%nz
    inu = photon0%inu
    tau = 0.0_wp;  t_cur = 0.0_wp
    nstep = 0
    do
      nstep = nstep + 1
      if (nstep > MAX_WALK_STEPS) then
        slab_walk_kills = slab_walk_kills + 1;  exit
      end if
      il  = 1 + ix + nx*(iy + ny*iz)
      a   = minloc(tmax, dim=1)
      seg = tmax(a) - t_cur
      tau = tau + kap_ion(inu, il)*seg
      if (tau >= tau_huge) exit
      t_cur = tmax(a)
      tmax(a) = tmax(a) + tdel(a)
      select case (a)
      case (1); ix = ix + stp(1)
                if (ix < 0 .or. ix >= nx) then
                  !--- slab: wrap x index (uniform grid: tmax/tdel unchanged).
                  if (par%xy_periodic) then; ix = modulo(ix, nx); else; exit; end if
                end if
      case (2); iy = iy + stp(2)
                if (iy < 0 .or. iy >= ny) then
                  if (par%xy_periodic) then; iy = modulo(iy, ny); else; exit; end if
                end if
      case (3); iz = iz + stp(3); if (iz < 0 .or. iz >= nz) exit
      end select
    end do
    tau_out = tau
  end subroutine raytrace_ion_tau_only_car

  !=========================================================================
  ! MoCHII peel-off: optical depth from the photon position to the box
  ! edge along the photon direction at its band bin — NO tally (the peel
  ! contribution is virtual; tallying it would double count).
  !=========================================================================
  subroutine raytrace_ion_tau_only_amr(photon0, tau_out)
    use define
    use gas_opacity_mod, only : kap_ion
    implicit none
    type(photon_type), intent(in)  :: photon0
    real(wp),          intent(out) :: tau_out

    integer  :: il, il_new, icell, iface, inu, nstep
    real(wp) :: x, y, z, kx, ky, kz, t_exit, tau

    if (amr_grid%car .and. trim(par%car_walk) /= 'shared') then
       call raytrace_ion_tau_only_car(photon0, tau_out)
       return
    end if

    x  = photon0%x;    y  = photon0%y;    z  = photon0%z
    kx = photon0%kx;   ky = photon0%ky;   kz = photon0%kz
    il  = photon0%icell_amr
    inu = photon0%inu
    tau = 0.0_wp
    if (il <= 0) then
      il = amr_find_leaf(x, y, z)
      if (il <= 0) then
        tau_out = 0.0_wp
        return
      end if
    end if

    nstep = 0
    do
      nstep = nstep + 1
      if (nstep > MAX_WALK_STEPS) then
        slab_walk_kills = slab_walk_kills + 1;  exit
      end if
      icell = leaf_cell(il)
      call amr_cell_exit(icell, x, y, z, kx, ky, kz, t_exit, iface)
      tau = tau + kap_ion(inu, il)*t_exit
      if (tau >= tau_huge) exit
      x = x + t_exit*kx
      y = y + t_exit*ky
      z = z + t_exit*kz
      il_new = amr_next_leaf(icell, iface, x, y, z)
      if (il_new <= 0 .and. par%xy_periodic) then
        call slab_wrap_xy(iface, x, y)
        if (iface <= 4) il_new = amr_find_leaf(x, y, z)
      end if
      if (il_new <= 0) exit
      il = il_new
    end do
    tau_out = tau
  end subroutine raytrace_ion_tau_only_amr

end module raytrace_amr_mod
