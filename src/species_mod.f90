module species_mod
!---------------------------------------------------------------------------
! MoCHII: trace-metal species registry + element ionization cascade.
!
! Design rule: adding an ion is a data operation.
! Each registered element is described by data/atomic/element_<el>.txt
! (produced by tools/fitting/make_element_data.py): VFKY96 photoionization
! parameters, Voronov collisional ionization, Badnell RR + DR of the
! recombining ion, and charge exchange with H (first transition).
!
! Trace approximation: metals never feed back on n_e or on the ionizing
! opacity.  After the H/He + thermal state of a cell is known, the cascade
!     n_i (Gamma_i + C_i n_e + k_CXI n_HII) =
!     n_{i+1} (alpha_{i+1} n_e + k_CXR n_HI)
! gives the stage fractions as a product chain (same structure as He).
! Metal Gamma_i per leaf are bin sums over the ionizing-band J tally with
! the element cross sections (computed once per iteration, all leaves).
!
! Line cooling: each (element, stage) with a coefficient file
! data/atomic/cooling_<el>_<i>.txt contributes
! n_e n_stage Lambda(T) to the thermal balance (metal_cooling).
!---------------------------------------------------------------------------
  use define
  use photo_xsec, only : sigma_vfky96
  use cooling_mod, only : cooling_fit_type, cooling_load, cooling_eval, gbar_ff
  use recomb_mod,  only : ci_dere_ratio
  use nlevel_cooling_mod, only : nlevel_cooling_add, nlevel_cooling_apply
  implicit none
  private

  public :: species_setup, species_gamma_compute, species_fractions
  public :: metal_cooling, metal_freefree, species_write, species_resize
  public :: species_opacity_add, species_ne, metal_heating
  public :: metal_cooling_H
  public :: species_ne_prepare, species_ne_cached
  public :: n_elements, elem_name, elem_nstage, elem_abund, elem_eth
  public :: species_shadow_setup, species_shadow_reset, species_shadow_score
  public :: species_shadow_score_exact
  public :: species_shadow_reduce, species_shadow_compare
  public :: species_stage_cache_refresh
  public :: species_packet_cross_sections, species_cached_opacity

  integer, parameter :: MAX_EL = 11, MAX_ST = 6

  type element_type
     character(len=2) :: name = ''
     integer  :: Z = 0
     integer  :: nstage = 0
     real(kind=wp) :: abund = 0.0_wp
     !--- transition data (index 1..nstage-1)
     real(kind=wp) :: eth(MAX_ST)      = 0.0_wp
     !--- photo_form 1: VFKY96 (Eth,E0,s0,ya,P,yw,y0,y1); 2: VY95 outer
     !--- shell (Eth,E0,s0,ya,P,yw,l) for the elements absent from the
     !--- VFKY96 table (Cl).
     integer       :: photo_form(MAX_ST) = 1
     real(kind=wp) :: photo(8,MAX_ST)  = 0.0_wp
     real(kind=wp) :: ci(5,MAX_ST)     = 0.0_wp   ! dE,P,A,X,K (Voronov)
     real(kind=wp) :: rr(6,MAX_ST)     = 0.0_wp   ! A,B,T0,T1,C,T2 (Badnell)
     integer       :: ndr(MAX_ST)      = 0
     real(kind=wp) :: dr_c(9,MAX_ST)   = 0.0_wp
     real(kind=wp) :: dr_e(9,MAX_ST)   = 0.0_wp
     !--- DR2 = Shull & Van Steenberg (1982) form (used where CHIANTI has
     !--- no Badnell DR, e.g. Ar II): A T^-3/2 e^-T0/T (1 + B e^-T1/T).
     logical       :: has_dr2(MAX_ST)  = .false.
     real(kind=wp) :: dr2(4,MAX_ST)    = 0.0_wp
     !--- RR2 = power-law RR alpha = A (T/1e4)^-eta (CHIANTI type 3,
     !--- e.g. Fe II); replaces the Badnell RR for that transition.
     logical       :: has_rr2(MAX_ST)  = .false.
     real(kind=wp) :: rr2(2,MAX_ST)    = 0.0_wp
     integer       :: cxi_form(MAX_ST) = 0
     real(kind=wp) :: cxi(6,MAX_ST)    = 0.0_wp
     integer       :: cxr_form(MAX_ST) = 0
     real(kind=wp) :: cxr(6,MAX_ST)    = 0.0_wp
     !--- cooling fits per stage (loaded when the file exists)
     logical :: has_cool(MAX_ST) = .false.
     type(cooling_fit_type) :: cool(MAX_ST)
     !--- index into the nlevel_cooling_mod (T, n_e) table for each stage
     !--- (par%cooling_model = 'local_ne'); 0 = no table, keep the fit.
     integer :: cool_idx(MAX_ST) = 0
  end type element_type

  integer :: n_elements = 0
  type(element_type), target :: elems(MAX_EL)

  !--- Gamma_i(transition, leaf) per element [s^-1]; filled each iteration.
  type gamma_block
     real(kind=wp), allocatable :: g(:,:)
  end type gamma_block
  type(gamma_block) :: egam(MAX_EL)
  !--- photoheating per particle [erg/s] (transition, leaf); filled with
  !--- the Gamma's, consumed by metal_heating (par%metal_heat).
  type(gamma_block) :: eheat(MAX_EL)
  !--- Direct metal-rate estimator blocks.  coeff(1,:,:) is the grouped
  !--- validation coefficient sigma/E and coeff(2,:,:) its heating
  !--- coefficient. raw carries unnormalized path sums and is also the
  !--- authoritative continuous-energy estimator.
  type metal_shadow_block
     real(kind=wp), allocatable :: coeff(:,:,:)  ! (2,ntransition,nnu)
     real(kind=wp), allocatable :: raw(:,:,:)    ! (2,ntransition,nleaf)
  end type metal_shadow_block
  type(metal_shadow_block) :: metal_shadow(MAX_EL)
  type(gamma_block) :: stage_frac(MAX_EL)  ! (stage,leaf), refreshed per gas state
  integer :: flat_ntransition = 0
  integer :: flat_ie(MAX_METAL_TRANSITIONS) = 0
  integer :: flat_it(MAX_METAL_TRANSITIONS) = 0
  real(kind=wp) :: flat_eth(MAX_METAL_TRANSITIONS) = 0.0_wp

  !--- rate cache for the n_e fixed point at one (leaf, T): all the
  !--- transition coefficients are functions of T (and the leaf Gammas)
  !--- only, so species_ne_prepare evaluates them ONCE and
  !--- species_ne_cached reduces each fixed-point iteration to
  !--- multiply-adds.  Without this, solve_ion_cell with par%metal_ne
  !--- re-evaluated every Badnell/Voronov/CX fit up to 200 times per
  !--- trial temperature (the 185-minute L6 PDR run).
  real(kind=wp) :: cch_gam(MAX_ST, MAX_EL), cch_ci(MAX_ST, MAX_EL)
  real(kind=wp) :: cch_rec(MAX_ST, MAX_EL), cch_cxi(MAX_ST, MAX_EL)
  real(kind=wp) :: cch_cxr(MAX_ST, MAX_EL)

contains

  character(len=2) function elem_name(ie);  integer, intent(in) :: ie
    elem_name = elems(ie)%name;  end function elem_name
  integer function elem_nstage(ie);  integer, intent(in) :: ie
    elem_nstage = elems(ie)%nstage;  end function elem_nstage
  real(kind=wp) function elem_abund(ie);  integer, intent(in) :: ie
    elem_abund = elems(ie)%abund;  end function elem_abund
  !--- photoionization threshold [eV] of transition it (stage it -> it+1),
  !--- it = 1..elem_nstage(ie)-1.  Exposed for the threshold-aligned band.
  real(kind=wp) function elem_eth(ie, it);  integer, intent(in) :: ie, it
    elem_eth = elems(ie)%eth(it);  end function elem_eth

  !=========================================================================
  subroutine species_setup(nleaf)
    use mpi
    implicit none
    integer, intent(in) :: nleaf
    character(len=8)   :: names(11)
    real(kind=wp)      :: abunds(11)
    character(len=256) :: fname
    logical :: exists
    integer :: k, ie, i, ierr, j, ie_tmp, it_tmp
    real(kind=wp) :: eth_tmp

    names  = [character(len=8) :: 'c', 'n', 'o', 'ne', 's', 'ar', 'mg', &
              'fe', 'si', 'cl', 'ca']
    abunds = [par%abund_C, par%abund_N, par%abund_O, par%abund_Ne, &
              par%abund_S, par%abund_Ar, par%abund_Mg, par%abund_Fe, &
              par%abund_Si, par%abund_Cl, par%abund_Ca]

    n_elements = 0
    do k = 1, size(names)
       if (abunds(k) <= 0.0_wp) cycle
       n_elements = n_elements + 1
       ie = n_elements
       call read_element(trim(par%atomic_dir)//'/element_'// &
                         trim(names(k))//'.txt', elems(ie))
       elems(ie)%abund = abunds(k)
       !--- cooling fits for each stage where a file exists.
       do i = 1, elems(ie)%nstage
          write(fname,'(a,a,a,i0,a)') trim(par%atomic_dir)// &
             '/cooling_', trim(names(k)), '_', i, '.txt'
          inquire(file=trim(fname), exist=exists)
          if (exists) then
             call cooling_load(trim(fname), elems(ie)%cool(i))
             elems(ie)%has_cool(i) = .true.
             !--- (T, n_e) suppression table for the local-n_e cooling model.
             if (trim(par%cooling_model) == 'local_ne') &
                elems(ie)%cool_idx(i) = nlevel_cooling_add(trim(names(k)), i)
          end if
       end do
       allocate(egam(ie)%g(elems(ie)%nstage-1, nleaf))
       egam(ie)%g = 0.0_wp
       allocate(eheat(ie)%g(elems(ie)%nstage-1, nleaf))
       eheat(ie)%g = 0.0_wp
       allocate(stage_frac(ie)%g(elems(ie)%nstage, nleaf))
       stage_frac(ie)%g = 0.0_wp
       if (mpar%p_rank == 0) write(*,'(3a,es10.3,a,i2,a,i2,a)') &
          ' SPEC: ', trim(names(k)), ' loaded (abund=', elems(ie)%abund, &
          ', ', elems(ie)%nstage, ' stages, ', &
          count(elems(ie)%has_cool(1:elems(ie)%nstage)), ' cooling fits)'
    end do
    flat_ntransition = 0
    do ie = 1, n_elements
       do i = 1, elems(ie)%nstage-1
          flat_ntransition = flat_ntransition + 1
          if (flat_ntransition > MAX_METAL_TRANSITIONS) then
             if (mpar%p_rank == 0) write(*,'(a,i0)') &
                'ERROR: active metal transitions exceed packet-cache capacity ', &
                MAX_METAL_TRANSITIONS
             call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          end if
          flat_ie(flat_ntransition) = ie
          flat_it(flat_ntransition) = i
          flat_eth(flat_ntransition) = elems(ie)%eth(i)
       end do
    end do
    !--- Packet energies are strongly concentrated near the H edge. Sort the
    !--- compact transition registry once so packet-cache construction can
    !--- stop at the first threshold above the sampled energy.
    do i = 2, flat_ntransition
       ie_tmp = flat_ie(i)
       it_tmp = flat_it(i)
       eth_tmp = flat_eth(i)
       j = i - 1
       do while (j >= 1)
          if (flat_eth(j) <= eth_tmp) exit
          flat_ie(j+1) = flat_ie(j)
          flat_it(j+1) = flat_it(j)
          flat_eth(j+1) = flat_eth(j)
          j = j - 1
       end do
       flat_ie(j+1) = ie_tmp
       flat_it(j+1) = it_tmp
       flat_eth(j+1) = eth_tmp
    end do
  end subroutine species_setup

  !=========================================================================
  !--- resize the Gamma blocks after octree re-refinement.
  subroutine species_resize(nleaf)
    implicit none
    integer, intent(in) :: nleaf
    integer :: ie
    do ie = 1, n_elements
       if (allocated(egam(ie)%g)) deallocate(egam(ie)%g)
       allocate(egam(ie)%g(elems(ie)%nstage-1, nleaf))
       egam(ie)%g = 0.0_wp
       if (allocated(eheat(ie)%g)) deallocate(eheat(ie)%g)
       allocate(eheat(ie)%g(elems(ie)%nstage-1, nleaf))
       eheat(ie)%g = 0.0_wp
       if (allocated(stage_frac(ie)%g)) deallocate(stage_frac(ie)%g)
       allocate(stage_frac(ie)%g(elems(ie)%nstage, nleaf))
       stage_frac(ie)%g = 0.0_wp
    end do
  end subroutine species_resize

  subroutine species_stage_cache_refresh()
    use gas_state_mod, only : gas_nH, gas_xHI, gas_ne, gas_Te
    implicit none
    real(kind=wp) :: frac(MAX_ST), nHI, nHII
    integer :: ie, il
    do ie = 1, n_elements
       do il = 1, size(stage_frac(ie)%g,2)
          nHI = gas_nH(il)*gas_xHI(il)
          nHII = gas_nH(il)*(1.0_wp-gas_xHI(il))
          call species_fractions(ie, il, gas_Te(il), gas_ne(il), &
                                 nHI, nHII, frac)
          stage_frac(ie)%g(:,il) = frac(1:elems(ie)%nstage)
       end do
    end do
  end subroutine species_stage_cache_refresh

  subroutine species_packet_cross_sections(energy, sigma, ntransition)
    real(kind=wp), intent(in)  :: energy
    real(kind=wp), intent(out) :: sigma(MAX_METAL_TRANSITIONS)
    integer,       intent(out) :: ntransition
    integer :: k
    sigma = 0.0_wp
    ntransition = 0
    do k = 1, flat_ntransition
       if (energy < flat_eth(k)) exit
       sigma(k) = species_sigma(flat_ie(k), flat_it(k), energy)
       ntransition = k
    end do
  end subroutine species_packet_cross_sections

  real(kind=wp) function species_cached_opacity(sigma, ntransition, il) result(opacity)
    real(kind=wp), intent(in) :: sigma(MAX_METAL_TRANSITIONS)
    integer,       intent(in) :: ntransition, il
    integer :: k, ie, it
    opacity = 0.0_wp
    do k = 1, ntransition
       ie = flat_ie(k)
       it = flat_it(k)
       opacity = opacity + elems(ie)%abund*stage_frac(ie)%g(it,il)*sigma(k)
    end do
  end function species_cached_opacity

  !=========================================================================
  ! Grouped metal-rate shadow estimator.  These routines intentionally live
  ! in species_mod so they can use the private registry and compare directly
  ! with egam/eheat without exposing implementation storage.
  !=========================================================================
  subroutine species_shadow_setup(nleaf)
    use ion_band_mod, only : ion_e, nnu_band
    implicit none
    integer, intent(in) :: nleaf
    real(kind=wp) :: energy, sig
    integer :: ie, it, inu, nt

    do ie = 1, n_elements
       if (allocated(metal_shadow(ie)%coeff)) deallocate(metal_shadow(ie)%coeff)
       if (allocated(metal_shadow(ie)%raw)) deallocate(metal_shadow(ie)%raw)
       nt = elems(ie)%nstage - 1
       allocate(metal_shadow(ie)%raw(2,nt,nleaf))
       metal_shadow(ie)%raw = 0.0_wp
       if (par%ion_shadow_rates) then
          allocate(metal_shadow(ie)%coeff(2,nt,nnu_band))
          do it = 1, nt
             do inu = 1, nnu_band
                energy = ion_e(inu)
                sig = species_sigma(ie, it, energy)
                metal_shadow(ie)%coeff(1,it,inu) = sig/(energy*ev2erg)
                metal_shadow(ie)%coeff(2,it,inu) = sig*max( &
                   1.0_wp - elems(ie)%eth(it)/energy, 0.0_wp)
             end do
          end do
       end if
    end do
  end subroutine species_shadow_setup

  subroutine species_shadow_reset()
    integer :: ie
    do ie = 1, n_elements
       if (allocated(metal_shadow(ie)%raw)) metal_shadow(ie)%raw = 0.0_wp
    end do
  end subroutine species_shadow_reset

  subroutine species_shadow_score(inu, il, path_lum)
    integer,       intent(in) :: inu, il
    real(kind=wp), intent(in) :: path_lum
    integer :: ie, it
    do ie = 1, n_elements
       do it = 1, elems(ie)%nstage-1
          metal_shadow(ie)%raw(:,it,il) = metal_shadow(ie)%raw(:,it,il) &
             + path_lum*metal_shadow(ie)%coeff(:,it,inu)
       end do
    end do
  end subroutine species_shadow_score

  subroutine species_shadow_score_exact(energy, sigma, ntransition, il, path_lum)
    real(kind=wp), intent(in) :: energy
    real(kind=wp), intent(in) :: sigma(MAX_METAL_TRANSITIONS)
    integer,       intent(in) :: ntransition, il
    real(kind=wp), intent(in) :: path_lum
    integer :: k, ie, it
    do k = 1, ntransition
       ie = flat_ie(k)
       it = flat_it(k)
       metal_shadow(ie)%raw(1,it,il) = metal_shadow(ie)%raw(1,it,il) &
          + path_lum*sigma(k)/(energy*ev2erg)
       metal_shadow(ie)%raw(2,it,il) = metal_shadow(ie)%raw(2,it,il) &
          + path_lum*sigma(k)*max(1.0_wp-elems(ie)%eth(it)/energy, 0.0_wp)
    end do
  end subroutine species_shadow_score_exact

  subroutine species_shadow_reduce()
    use mpi
    implicit none
    integer :: ie, ierr
    do ie = 1, n_elements
       call MPI_ALLREDUCE(MPI_IN_PLACE, metal_shadow(ie)%raw, &
                          size(metal_shadow(ie)%raw), MPI_DOUBLE_PRECISION, &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
    end do
  end subroutine species_shadow_reduce

  subroutine species_shadow_compare(rtol)
    use octree_mod, only : leaf_half
    use mpi
    implicit none
    real(kind=wp), intent(in) :: rtol
    real(kind=wp) :: max_abs(2), scale(2), rel(2), reference(2)
    real(kind=wp) :: fac, shadow
    integer :: ie, it, il, i, ierr
    logical :: passed

    max_abs = 0.0_wp
    scale = 0.0_wp
    do ie = 1, n_elements
       do it = 1, elems(ie)%nstage-1
          do il = 1, size(egam(ie)%g,2)
             fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
             reference = [egam(ie)%g(it,il), eheat(ie)%g(it,il)]
             do i = 1, 2
                shadow = metal_shadow(ie)%raw(i,it,il)*fac
                max_abs(i) = max(max_abs(i), abs(shadow-reference(i)))
                scale(i) = max(scale(i), abs(reference(i)))
             end do
          end do
       end do
    end do
    do i = 1, 2
       if (scale(i) > 0.0_wp) then
          rel(i) = max_abs(i)/scale(i)
       else
          rel(i) = max_abs(i)
       end if
    end do
    passed = maxval(rel) <= rtol
    if (mpar%p_rank == 0) then
       write(*,'(a,2es11.3)') ' ION: metal shadow relative errors: ', rel
       write(*,'(a,a)') ' ION: grouped metal direct-rate shadow: ', &
                        merge('PASS', 'FAIL', passed)
    end if
    if (.not. passed) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end subroutine species_shadow_compare

  !=========================================================================
  subroutine read_element(fname, el)
    use mpi
        use utility, only : fatal_error
implicit none
    character(len=*),   intent(in)  :: fname
    type(element_type), intent(out) :: el
    character(len=512) :: line, key
    integer :: unit, ios, it, ierr, n

    open(newunit=unit, file=fname, status='old', iostat=ios)
    if (ios /= 0) then
       call fatal_error('cannot open '//trim(fname))
    end if
    it = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) key
       select case (trim(key))
       case ('ELEMENT')
          read(line,*) key, el%name, el%Z, el%nstage
       case ('TRANSITION')
          read(line,*) key, it
       case ('ETH')
          read(line,*) key, el%eth(it)
       case ('PHOTO')
          read(line,*) key, el%photo(1:8,it)
          el%photo_form(it) = 1
       case ('PHOTO2')
          read(line,*) key, el%photo(1:7,it)
          el%photo_form(it) = 2
       case ('CI')
          read(line,*) key, el%ci(1:5,it)
       case ('RR')
          read(line,*) key, el%rr(1:6,it)
       case ('DR')
          read(line,*) key, n
          el%ndr(it) = n
          if (n > 0) then
             backspace(unit)
             read(unit,*) key, n, el%dr_c(1:n,it), el%dr_e(1:n,it)
          end if
       case ('DR2')
          read(line,*) key, el%dr2(1:4,it)
          el%has_dr2(it) = .true.
       case ('RR2')
          read(line,*) key, el%rr2(1:2,it)
          el%has_rr2(it) = .true.
       case ('CXI')
          read(line,*) key, el%cxi_form(it)
          if (el%cxi_form(it) > 0) then
             backspace(unit)
             read(unit,*) key, el%cxi_form(it), el%cxi(1:6,it)
          end if
       case ('CXR')
          read(line,*) key, el%cxr_form(it)
          if (el%cxr_form(it) > 0) then
             backspace(unit)
             read(unit,*) key, el%cxr_form(it), el%cxr(1:6,it)
          end if
       end select
    end do
    close(unit)
  end subroutine read_element

  !=========================================================================
  ! Badnell RR + DR recombination of the recombining ion of transition it.
  !=========================================================================
  real(kind=wp) function alpha_rec(el, it, T) result(a)
    type(element_type), intent(in) :: el
    integer,            intent(in) :: it
    real(kind=wp),      intent(in) :: T
    real(kind=wp) :: b, s0, s1
    integer :: k
    if (el%has_rr2(it)) then
       a = el%rr2(1,it)*(T/1.0e4_wp)**(-el%rr2(2,it))
    else
       s0 = sqrt(T/el%rr(3,it));  s1 = sqrt(T/el%rr(4,it))
       b  = el%rr(2,it) + el%rr(5,it)*exp(-el%rr(6,it)/T)
       a  = el%rr(1,it) / ( s0 * (1.0_wp+s0)**(1.0_wp-b) * (1.0_wp+s1)**(1.0_wp+b) )
    end if
    do k = 1, el%ndr(it)
       a = a + T**(-1.5_wp)*el%dr_c(k,it)*exp(-el%dr_e(k,it)/T)
    end do
    if (el%has_dr2(it)) &
       a = a + el%dr2(1,it)*T**(-1.5_wp)*exp(-el%dr2(3,it)/T) &
               *(1.0_wp + el%dr2(2,it)*exp(-el%dr2(4,it)/T))
  end function alpha_rec

  !=========================================================================
  ! Charge-exchange rate [cm^3 s^-1]; forms per make_element_data.py.
  !=========================================================================
  real(kind=wp) function cx_rate(form, c, T) result(k)
    integer,       intent(in) :: form
    real(kind=wp), intent(in) :: c(6), T
    real(kind=wp) :: t4, lnT
    t4 = T/1.0e4_wp
    select case (form)
    case (1)
       k = c(1)*(T/300.0_wp)**c(2)
       if (c(3) /= 0.0_wp) k = k*exp(-c(3)/t4)
    case (2)
       k = c(1)*t4**c(2) + c(3)*t4**c(4)
    case (3)
       k = (c(1)*t4**c(2) + c(3)*t4**c(4))*exp(-c(5)/T)
    case (4)
       lnT = log(T)
       k = exp(c(1) + lnT*(c(2) + lnT*(c(3) + lnT*(c(4) + lnT*c(5)))))
    case (5)
       lnT = log(T)
       k = exp(c(1) + lnT*(c(2) + lnT*(c(3) + lnT*(c(4) + lnT*c(5)))) &
               - c(6)/t4)
    case (6)
       !--- Kingdon & Ferland (1996), evaluated at T clamped into the fit
       !--- validity range (~5e3-5e4 K; MOCASSIN zeroes outside instead).
       t4 = min(max(t4, 0.5_wp), 5.0_wp)
       k = c(1)*t4**c(2)*(1.0_wp + c(3)*exp(c(4)*t4))
       k = max(k, 0.0_wp)
    case (7)
       !--- KF96 form with an activation barrier (endothermic reverse).
       t4 = min(max(t4, 0.5_wp), 5.0_wp)
       k = c(1)*t4**c(2)*(1.0_wp + c(3)*exp(c(4)*t4))*exp(-c(5)/t4)
       k = max(k, 0.0_wp)
    case default
       k = 0.0_wp
    end select
  end function cx_rate

  !=========================================================================
  ! Metal photoionization absorption added to kap_ion (par%ion_metal_abs):
  ! kap(inu,il) += n_H abund Sum_i frac_i sigma_VFKY96,i(E_inu) per code
  ! length.  Negligible next to H/He above 13.6 eV; with par%add_fuv it is
  ! the only GAS opacity in the FUV bins (Mg I 7.65, C I 11.26, S I 10.36,
  ! Fe I 7.90 eV thresholds).  Stage fractions come from the current state
  ! (same product chain as the cooling and the output), so the opacity
  ! feedback iterates them exactly like x_HI.  Caller: gas_opacity_fill on
  ! h_rank 0 (kap lives in node-shared memory).
  !=========================================================================
  subroutine species_opacity_add(kap, nnu, nleaf)
    use ion_band_mod,  only : ion_e
    use gas_state_mod, only : gas_nH
    implicit none
    integer,       intent(in)    :: nnu, nleaf
    real(kind=wp), intent(inout) :: kap(nnu, nleaf)
    real(kind=wp) :: sig(nnu, MAX_ST), add
    integer :: ie, it, il, inu

    do ie = 1, n_elements
       do it = 1, elems(ie)%nstage-1
          do inu = 1, nnu
             sig(inu,it) = species_sigma(ie, it, ion_e(inu))
          end do
       end do
       do il = 1, nleaf
          if (gas_nH(il) <= 0.0_wp) cycle
          do inu = 1, nnu
             add = 0.0_wp
             do it = 1, elems(ie)%nstage-1
                add = add + stage_frac(ie)%g(it,il)*sig(inu,it)
             end do
             kap(inu,il) = kap(inu,il) &
                + gas_nH(il)*elems(ie)%abund*add*par%distance2cm
          end do
       end do
    end do
  end subroutine species_opacity_add

  !=========================================================================
  ! Metal photoionization rate integrals from the ionizing-band J tally
  ! (all elements, all leaves); call after jtally_ion_reduce each iteration.
  !=========================================================================
  subroutine species_gamma_compute()
    use octree_mod, only : amr_grid, leaf_half
    use jtally_mod,   only : jt_ion
    use ion_band_mod, only : ion_e, nnu_band
    implicit none
    integer :: ie, it, il, inu, ic, nleaf
    real(kind=wp) :: sig(nnu_band), vol, fac, fJ, fH

    nleaf = amr_grid%nleaf
    if (trim(par%ion_energy_mode) == 'continuous') then
       do ie = 1, n_elements
          if (.not. allocated(metal_shadow(ie)%raw)) &
             error stop 'continuous metal path estimators are not allocated'
          do it = 1, elems(ie)%nstage-1
             do il = 1, nleaf
                vol = (2.0_wp*leaf_half(il))**3
                fac = 1.0_wp/(vol*par%distance2cm**2)
                egam(ie)%g(it,il) = metal_shadow(ie)%raw(1,it,il)*fac
                eheat(ie)%g(it,il) = metal_shadow(ie)%raw(2,it,il)*fac
             end do
          end do
       end do
       return
    end if

    do ie = 1, n_elements
       do it = 1, elems(ie)%nstage-1
          do inu = 1, nnu_band
             sig(inu) = species_sigma(ie, it, ion_e(inu))
          end do
          do il = 1, nleaf
             vol = (2.0_wp*leaf_half(il))**3
             fac = 1.0_wp/(vol*par%distance2cm**2)
             fJ  = 0.0_wp
             fH  = 0.0_wp
             do inu = 1, nnu_band
                fJ = fJ + jt_ion(inu,il)*sig(inu)/(ion_e(inu)*ev2erg)
                fH = fH + jt_ion(inu,il)*sig(inu) &
                     *max(1.0_wp - el8(ie,it,1)/ion_e(inu), 0.0_wp)
             end do
             egam(ie)%g(it,il)  = fJ*fac
             eheat(ie)%g(it,il) = fH*fac
          end do
       end do
    end do
  end subroutine species_gamma_compute

  !=========================================================================
  ! Fill the (leaf, T) rate cache for the n_e fixed point.
  !=========================================================================
  subroutine species_ne_prepare(il, T)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T
    integer :: ie, it
    do ie = 1, n_elements
       do it = 1, elems(ie)%nstage-1
          cch_gam(it,ie) = egam(ie)%g(it,il)
          cch_ci(it,ie)  = voronov_ci(elems(ie)%ci(1:5,it), T) &
                           * ci_dere_ratio(elems(ie)%Z, it)
          cch_rec(it,ie) = alpha_rec(elems(ie), it, T)
          cch_cxi(it,ie) = 0.0_wp
          cch_cxr(it,ie) = 0.0_wp
          if (elems(ie)%cxi_form(it) > 0) cch_cxi(it,ie) = &
             cx_rate(elems(ie)%cxi_form(it), elems(ie)%cxi(1:6,it), T)
          if (elems(ie)%cxr_form(it) > 0) cch_cxr(it,ie) = &
             cx_rate(elems(ie)%cxr_form(it), elems(ie)%cxr(1:6,it), T)
       end do
    end do
  end subroutine species_ne_prepare

  !=========================================================================
  ! Metal electrons from the cached rates (same product chain and the
  ! same arithmetic order as species_fractions, so the results are
  ! identical to the uncached path).
  !=========================================================================
  real(kind=wp) function species_ne_cached(nH, ne, nHI, nHII) result(nem)
    implicit none
    real(kind=wp), intent(in) :: nH, ne, nHI, nHII
    real(kind=wp) :: r(MAX_ST), rion, rrec, s, prod, fr
    integer :: ie, it, i, ns
    nem = 0.0_wp
    do ie = 1, n_elements
       ns = elems(ie)%nstage
       do it = 1, ns-1
          rion = cch_gam(it,ie) + ne*cch_ci(it,ie)
          if (cch_cxi(it,ie) > 0.0_wp) rion = rion + nHII*cch_cxi(it,ie)
          rrec = ne*cch_rec(it,ie)
          if (cch_cxr(it,ie) > 0.0_wp) rrec = rrec + nHI*cch_cxr(it,ie)
          r(it) = rion/max(rrec, tinest)
       end do
       s = 1.0_wp;  prod = 1.0_wp
       do i = 1, ns-1
          prod = prod*r(i)
          s = s + prod
       end do
       fr = 1.0_wp/s
       do i = 2, ns
          fr = fr*r(i-1)
          nem = nem + elems(ie)%abund*nH*fr*real(i-1, wp)
       end do
    end do
  end function species_ne_cached

  !=========================================================================
  ! Electrons contributed by the metal cascade [cm^-3] at the current
  ! (T, n_e) of one leaf — the n_e closure extension of par%metal_ne
  ! (beyond the I-front these are the only electrons).
  !=========================================================================
  real(kind=wp) function species_ne(il, T, ne, nHI, nHII) result(nem)
    use gas_state_mod, only : gas_nH
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, ne, nHI, nHII
    real(kind=wp) :: frac(MAX_ST)
    integer :: ie, i
    nem = 0.0_wp
    do ie = 1, n_elements
       call species_fractions(ie, il, T, ne, nHI, nHII, frac)
       do i = 2, elems(ie)%nstage
          nem = nem + elems(ie)%abund*gas_nH(il)*frac(i)*real(i-1, wp)
       end do
    end do
  end function species_ne

  !=========================================================================
  ! H-impact fine-structure cooling [erg s^-1 cm^-3]: [C II] 158 um and
  ! [O I] 63 um excited by NEUTRAL HYDROGEN collisions — the dominant
  ! PDR-zone coolant (n_HI/n_e ~ 10^3-10^4 there; the n-level cooling table
  ! is electron-impact only, so the H channel is added separately).  The
  ! two-level population saturates above the H critical density, so the
  ! cooling is the low-density excitation power divided by the two-level
  ! saturation factor.  The steady state n_u (A_ul + n_HI q_ul^H)
  ! = n_l n_HI q_lu^H with n_ion = n_l + n_u gives
  !   Lambda = n_ion n_HI q_lu^H dE / (1 + n_HI (q_ul^H + q_lu^H)/A_ul)
  ! and hence n_crit,H = A_ul / (q_ul^H + q_lu^H).  BOTH collision
  ! directions drain the lower level, so q_lu^H belongs in the denominator:
  ! omitting it leaves the low-density limit right but overestimates the
  ! dense limit by 1 + (g_u/g_l) e^{-dE/kT}, i.e. 1.8 at 100 K for [C II]
  ! rising to 3 when kT >> dE, and up to 1.6 for [O I].
  ! With q_lu^H = (g_u/g_l) q_ul^H e^{-dE/kT},
  ! q_ul^H([C II]) = 7.6e-10 (T/100)^0.14 (Barinovs et al. 2005),
  ! q_ul^H([O I]) = 9.2e-11 (T/100)^0.67 cm^3/s, A([C II] 158um) = 2.321e-6
  ! and A([O I] 63um) = 8.865e-5 s^-1.  n_crit,H is ~1.7e3 for [C II] at
  ! 100 K (so the factor bites in dense PDR gas, nH ~ 5e3) and ~9e5 for
  ! [O I] (rarely saturated); at low n_HI it reduces to the low-density
  ! limit.  The electron channel is left to the density-suppressed n-level
  ! table (its own saturation), so only n_HI enters the factor here.  That
  ! split is valid while n_e << n_HI, which is the PDR regime this channel
  ! matters in; where both colliders are comparable the two saturation
  ! denominators should be merged into one.  Called from the par%use_metals
  ! cooling path: it is the coolant that keeps the neutral zone off the
  ! ~10^4 K runaway under grain photoelectric heating, but it does not
  ! depend on that heating being enabled.
  !=========================================================================
  real(kind=wp) function metal_cooling_H(il, T, nH, ne, nHI, nHII) result(cool)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne, nHI, nHII
    real(kind=wp), parameter :: kb = 1.380649e-16_wp
    real(kind=wp), parameter :: A_CII = 2.321e-6_wp, A_OI = 8.865e-5_wp
    real(kind=wp) :: frac(MAX_ST), nion, qlu, qul
    integer :: ie
    cool = 0.0_wp
    do ie = 1, n_elements
       select case (trim(elems(ie)%name))
       case ('c')
          if (elems(ie)%nstage < 2) cycle
          call species_fractions(ie, il, T, ne, nHI, nHII, frac)
          nion = elems(ie)%abund*nH*frac(2)                 ! C II
          qul  = 7.6e-10_wp*(T/100.0_wp)**0.14_wp           ! H de-excitation
          qlu  = 2.0_wp*qul*exp(-91.25_wp/T)                ! g_u/g_l = 4/2
          cool = cool + nion*nHI*qlu*(91.25_wp*kb) &
                 / (1.0_wp + nHI*(qul + qlu)/A_CII)
       case ('o')
          call species_fractions(ie, il, T, ne, nHI, nHII, frac)
          nion = elems(ie)%abund*nH*frac(1)                 ! O I
          qul  = 9.2e-11_wp*(T/100.0_wp)**0.67_wp           ! H de-excitation
          qlu  = 0.6_wp*qul*exp(-227.7_wp/T)                ! g_u/g_l = 3/5
          cool = cool + nion*nHI*qlu*(227.7_wp*kb) &
                 / (1.0_wp + nHI*(qul + qlu)/A_OI)
       end select
    end do
  end function metal_cooling_H

  !=========================================================================
  ! Metal photoheating [erg s^-1 cm^-3] (par%metal_heat): the heating
  ! integrals accumulated with the Gamma's, weighted by the cascade
  ! stage fractions.
  !=========================================================================
  real(kind=wp) function metal_heating(il, T, nH, ne, nHI, nHII) result(heat)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne, nHI, nHII
    real(kind=wp) :: frac(MAX_ST)
    integer :: ie, i
    heat = 0.0_wp
    do ie = 1, n_elements
       call species_fractions(ie, il, T, ne, nHI, nHII, frac)
       do i = 1, elems(ie)%nstage-1
          heat = heat + elems(ie)%abund*nH*frac(i)*eheat(ie)%g(i,il)
       end do
    end do
  end function metal_heating

  real(kind=wp) function el8(ie, it, k)
    integer, intent(in) :: ie, it, k
    el8 = elems(ie)%photo(k, it)
  end function el8

  !=========================================================================
  ! Photoionization cross section of (element ie, stage it) at E [eV]:
  ! VFKY96 or the VY95 outer-shell form by photo_form.
  !=========================================================================
  real(kind=wp) function species_sigma(ie, it, E) result(sig)
    use photo_xsec, only : sigma_vy95
    integer,       intent(in) :: ie, it
    real(kind=wp), intent(in) :: E
    if (elems(ie)%photo_form(it) == 2) then
       sig = sigma_vy95(E, el8(ie,it,1), el8(ie,it,2), el8(ie,it,3), &
                        el8(ie,it,4), el8(ie,it,5), el8(ie,it,6), &
                        el8(ie,it,7))
    else
       sig = sigma_vfky96(E, el8(ie,it,1), el8(ie,it,2), el8(ie,it,3), &
                          el8(ie,it,4), el8(ie,it,5), el8(ie,it,6), &
                          el8(ie,it,7), el8(ie,it,8))
    end if
  end function species_sigma

  !=========================================================================
  ! Stage fractions of element ie in leaf il at (T, n_e, x_HI): the product
  ! chain n_{i+1}/n_i = R_ion(i) / R_rec(i).
  !=========================================================================
  subroutine species_fractions(ie, il, T, ne, nHI, nHII, frac)
    implicit none
    integer,       intent(in)  :: ie, il
    real(kind=wp), intent(in)  :: T, ne, nHI, nHII
    real(kind=wp), intent(out) :: frac(MAX_ST)
    real(kind=wp) :: r(MAX_ST), rion, rrec, s, prod
    integer :: it, i, ns

    ns = elems(ie)%nstage
    do it = 1, ns-1
       rion = egam(ie)%g(it,il) &
            + ne*voronov_ci(elems(ie)%ci(1:5,it), T) &
               *ci_dere_ratio(elems(ie)%Z, it)
       if (elems(ie)%cxi_form(it) > 0) &
          rion = rion + nHII*cx_rate(elems(ie)%cxi_form(it), elems(ie)%cxi(1:6,it), T)
       rrec = ne*alpha_rec(elems(ie), it, T)
       if (elems(ie)%cxr_form(it) > 0) &
          rrec = rrec + nHI*cx_rate(elems(ie)%cxr_form(it), elems(ie)%cxr(1:6,it), T)
       r(it) = rion/max(rrec, tinest)
    end do
    !--- fractions from the ratio chain, normalized.
    s = 1.0_wp;  prod = 1.0_wp
    do i = 1, ns-1
       prod = prod*r(i)
       s = s + prod
    end do
    frac(1) = 1.0_wp/s
    do i = 2, ns
       frac(i) = frac(i-1)*r(i-1)
    end do
  end subroutine species_fractions

  !=========================================================================
  pure real(kind=wp) function voronov_ci(c, T) result(k)
    real(kind=wp), intent(in) :: c(5), T
    real(kind=wp) :: U
    U = c(1)*ev2erg/(kboltz_cgs*T)
    k = c(3)*(1.0_wp + c(2)*sqrt(U))*U**c(5)*exp(-U)/(c(4) + U)
  end function voronov_ci

  !=========================================================================
  ! Metal line cooling per unit volume [erg cm^-3 s^-1] at trial T.
  !=========================================================================
  real(kind=wp) function metal_cooling(il, T, nH, ne, nHI, nHII) result(cool)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne, nHI, nHII
    real(kind=wp) :: frac(MAX_ST), nel
    integer :: ie, i

    cool = 0.0_wp
    do ie = 1, n_elements
       call species_fractions(ie, il, T, ne, nHI, nHII, frac)
       nel = elems(ie)%abund*nH
       do i = 1, elems(ie)%nstage
          if (elems(ie)%has_cool(i)) &
             cool = cool + ne*nel*frac(i)*nlevel_cooling_apply(elems(ie)%cool_idx(i), &
                    cooling_eval(elems(ie)%cool(i), T), T, ne)
       end do
    end do
  end function metal_cooling

  !=========================================================================
  ! Metal free-free (thermal bremsstrahlung) cooling per unit volume
  ! [erg cm^-3 s^-1] at trial T (par%metal_freefree): the trace completion
  ! of the H II + He II/III free-free term.  A metal ion in stage j
  ! (j = 1 neutral, j = 2 = X+, ...) has net charge Z_ion = j - 1 and
  ! contributes n_ion Z_ion^2 gbar_ff(T, Z_ion); only j >= 2 radiate.
  ! Same prefactor 1.42554e-27 sqrt(T) n_e as the H/He free-free.
  !=========================================================================
  real(kind=wp) function metal_freefree(il, T, nH, ne, nHI, nHII) result(cool_ff)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne, nHI, nHII
    real(kind=wp) :: frac(MAX_ST), s
    integer :: ie, j, zion

    s = 0.0_wp
    do ie = 1, n_elements
       call species_fractions(ie, il, T, ne, nHI, nHII, frac)
       do j = 2, elems(ie)%nstage
          zion = j - 1
          s = s + elems(ie)%abund*frac(j)*real(zion*zion,wp)*gbar_ff(T, zion)
       end do
    end do
    cool_ff = 1.42554e-27_wp*sqrt(T)*ne*nH*s
  end function metal_freefree

  !=========================================================================
  ! Converged stage fractions of every element -> output file blocks.
  !=========================================================================
  subroutine species_write(file, status)
    use iofile_mod
    use octree_mod,    only : amr_grid
    use gas_state_mod, only : gas_nH, gas_xHI, gas_ne, gas_Te
    implicit none
    type(io_file_type), intent(inout) :: file
    integer, intent(inout) :: status
    real(kind=wp), allocatable :: fr(:,:)
    real(kind=wp) :: frac(MAX_ST), nHI, nHII
    character(len=64) :: extname
    integer :: ie, il, i, nleaf

    if (status /= 0) return
    nleaf = amr_grid%nleaf
    do ie = 1, n_elements
       allocate(fr(elems(ie)%nstage, nleaf))
       do il = 1, nleaf
          nHI  = gas_nH(il)*gas_xHI(il)
          nHII = gas_nH(il)*(1.0_wp - gas_xHI(il))
          if (gas_nH(il) > 0.0_wp) then
             call species_fractions(ie, il, gas_Te(il), gas_ne(il), &
                                    nHI, nHII, frac)
          else
             frac = 0.0_wp;  frac(1) = 1.0_wp
          end if
          fr(:,il) = frac(1:elems(ie)%nstage)
       end do
       call io_append_image(file, fr, status, bitpix=-64)
       write(extname,'(3a)') 'x_', trim(elems(ie)%name), '_stages'
       call io_put_keyword(file,'EXTNAME',trim(extname), &
            'ion fractions (stage, leaf)',status)
       call io_put_keyword(file,'ABUND',elems(ie)%abund,'n(X)/n(H)',status)
       deallocate(fr)

       call io_append_image(file, egam(ie)%g, status, bitpix=-64)
       write(extname,'(3a)') 'Gamma_', trim(elems(ie)%name), '_stages'
       call io_put_keyword(file,'EXTNAME',trim(extname), &
            'photoionization rates (transition, leaf) [s^-1]',status)
       call io_append_image(file, eheat(ie)%g, status, bitpix=-64)
       write(extname,'(3a)') 'Heat_', trim(elems(ie)%name), '_stages'
       call io_put_keyword(file,'EXTNAME',trim(extname), &
            'photoheating rates (transition, leaf) [erg/s per ion]',status)
    end do
  end subroutine species_write

end module species_mod
