module species_mod
!---------------------------------------------------------------------------
! MoCHII: trace-metal species registry + element ionization cascade.
!
! Design rule: adding an ion is a data operation.
! Each registered element is described by data/atomic/element_<el>.txt
! (produced by tools/fitting/make_element_data.py): VFKY96 photoionization
! and VY95 subshell parameters, Voronov collisional ionization, RR + DR of
! the recombining ion (Badnell, other sources where Badnell has none), and
! charge exchange with H for the transitions that have published rates.
!
! Trace approximation, and its three exceptions.  The cascade is evaluated on
! a known H/He and thermal state rather than solved jointly with it, but the
! metals are NOT excluded from the ionizing opacity, from n_e, or from
! hydrogen's own balance: par%ion_metal_abs is on by default and removing it
! moves the HII40 He front 0.615 pc; par%metal_ne feeds the stage charges into
! the n_e closure as a fixed point (charge_neutrality_ne); and with
! par%charge_exchange = 'two_way' the k_CXI and k_CXR terms below enter
! hydrogen's balance as well as the cascade's, so the charge each reaction
! moves is conserved in the reacting pair (species_metal_sums,
! ion_balance_mod::solve_ion_cell; md/METAL_COUPLING_PLAN.md, Stage M1).
! An earlier version of this comment claimed metals never feed back on n_e or
! on the opacity; it is withdrawn.
! After the H/He + thermal state of a cell is known, the cascade
!     n_i (Gamma_i + C_i n_e + k_CXI n_HII) =
!     n_{i+1} (alpha_{i+1} n_e + k_CXR n_HI)
! gives the stage fractions (ion_stage_equilibrium_mod::adjacent_ion_fractions).
! Metal Gamma_i per leaf are bin sums over the ionizing-band J tally with
! the element cross sections (computed once per iteration, all leaves).
!
! Line cooling: each (element, stage) with an atom in the selected excitation
! model, data/atomic/excitation/<par%metal_excitation_model>/nlevel_<el>_<i>.txt
! (with its fit cooling_<el>_<i>.txt), contributes n_e n_stage Lambda(T, n_e)
! to the thermal balance (metal_cooling).
!---------------------------------------------------------------------------
  use define
  use mpi
  use memory_mod, only : create_shared_mem, destroy_shared_mem, get_window, &
                        shared_host_rank, shared_host_size, shared_host_comm
  use metal_ionization_potential_mod, only: metal_ionization_potential
  use ion_stage_equilibrium_mod, only : adjacent_ion_fractions
  use utility, only : fatal_error
  use photo_xsec, only : sigma_vfky96
  use cooling_mod, only : cooling_fit_type, cooling_load, cooling_eval, gbar_ff
  use free_free_mod, only : FF_TOTAL
  use recomb_mod,  only : ci_dere_ratio
  use nlevel_cooling_mod, only : nlevel_cooling_add, nlevel_cooling_apply, &
                                nlevel_cooling_lam_se
  use nlevel_mod, only : nlevel_atom_type, nlevel_load, nlevel_has_hcoll, &
                         nlevel_hcoll_cooling, metal_excitation_data_file
  implicit none
  private

  public :: species_setup, species_gamma_compute, species_fractions
  public :: metal_cooling, metal_freefree, species_write, species_resize
  public :: species_ne, metal_heating, metal_ionization_cooling
  public :: metal_cooling_H
  public :: species_ne_prepare, species_ne_cached, species_metal_sums
  public :: species_cx_coefficients, species_cx_kinetic_energy
  public :: metal_charge_exchange_heating
  public :: species_terminal_photoionization_check
  public :: n_elements, elem_name, elem_nstage, elem_abund, elem_eth
  public :: elem_atomic_number
  public :: elem_cx_form
  public :: species_rate_path_setup, species_rate_path_reset
  public :: species_rate_path_score, species_rate_path_reduce
  public :: species_stage_cache_refresh, species_set_photoionization_rates
  public :: species_photoionization_rates
  public :: species_packet_cross_sections, species_cached_opacity
  public :: species_cross_section, species_outer_shell_cross_section
  public :: species_transition_count, species_transition
  public :: species_stage_fraction, species_recombination_coefficient
  public :: species_radiative_recombination, species_outer_shell_break
  public :: species_photo_table_end, species_photo_table_nodes
  public :: species_photo_table_node_energy

  integer, parameter :: MAX_EL = MAX_METAL_ELEMENTS, MAX_ST = MAX_ION_STAGES
  !--- subshells of a ground-state ion in phfit2 (1s 2s 2p 3s 3p 3d 4s)
  integer, parameter :: MAX_SHELL = 7
  !--- points of a tabulated charge-exchange rate coefficient (CXRTAB)
  integer, parameter :: MAX_CX_TAB = 40
  ! Rate grids resolve the curvature of published KF96 fits to 1e-5.
  ! Heat-defect tables retain their separate 40-point capacity.
  integer, parameter :: MAX_CXR_TAB = 1024
  !--- charge exchange resolved by product term (form 9, CXSTATE rows):
  !--- product terms and tabulated temperatures
  integer, parameter :: MAX_CXS_TERM = 4, MAX_CXS_T = 40
  !--- nodes of a tabulated ground-term photoionization cross section
  !--- (PHTABLE; Ne III carries 25013)
  integer, parameter :: MAX_PHTABLE = 200000

  !--- PHTABLE: photoionization cross section of the ground term tabulated
  !--- against the photoelectron energy eps = E - E_th [eV], linear in eps;
  !--- it replaces the VFKY96 fit below E_th + eps(n), where it ends on the
  !--- fit's value (the reader checks both).  Source, choices and validity:
  !--- tools/fitting/ground_term_photoionization_table.py.
  type photo_table_type
     real(kind=wp), allocatable :: eps(:)     ! [eV], eps(1) = 0
     real(kind=wp), allocatable :: sig(:)     ! [cm^2]
  end type photo_table_type
  !--- the (log T, log n_e) grid of their term populations
  integer, parameter :: CXS_NT = 201, CXS_NNE = 161
  real(kind=wp), parameter :: CXS_T_LO = 2.0_wp, CXS_DT = 0.02_wp
  real(kind=wp), parameter :: CXS_NE_LO = -6.0_wp, CXS_DNE = 0.1_wp

  type element_type
     character(len=2) :: name = ''
     integer  :: Z = 0
     integer  :: nstage = 0
     real(kind=wp) :: abund = 0.0_wp
     !--- transition data (index 1..nstage-1)
     real(kind=wp) :: eth(MAX_ST)      = 0.0_wp
     !--- Photoionization cross section, Verner's phfit2 prescription in
     !--- full: below the first inner-shell edge einn the VFKY96 outer-shell
     !--- fit photo(1:8) = (Eth,E0,s0,ya,P,yw,y0,y1) (has_ph2); at E >= einn,
     !--- and at every E when einn = 0 (elements absent from the VFKY96
     !--- table: Cl), the sum over the occupied subshells of the Verner &
     !--- Yakovlev (1995) fits shell(1:7,k) = (l,Eth,E0,s0,ya,P,yw).
     logical       :: has_ph2(MAX_ST)  = .false.
     real(kind=wp) :: photo(8,MAX_ST)  = 0.0_wp
     real(kind=wp) :: einn(MAX_ST)     = 0.0_wp
     integer       :: nshell(MAX_ST)   = 0
     real(kind=wp) :: shell(7,MAX_SHELL,MAX_ST) = 0.0_wp
     !--- PHTABLE (photo_table_type above): nphtab > 0 replaces the VFKY96
     !--- fit below eth + phtab%eps(nphtab)
     integer       :: nphtab(MAX_ST)   = 0
     type(photo_table_type) :: phtab(MAX_ST)
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
     !--- cxr_form 8: tabulated recombination rate coefficient (CXRTAB row),
     !--- interpolated log-log in T and clamped to the table's range.
     integer       :: ncxr_tab(MAX_ST) = 0
     real(kind=wp) :: cxr_tab_T(MAX_CXR_TAB,MAX_ST) = 0.0_wp
     real(kind=wp) :: cxr_tab_k(MAX_CXR_TAB,MAX_ST) = 0.0_wp
     !--- cxi_form = cxr_form = 9: charge exchange resolved by the product
     !--- term of the ion X^(it+1) (CXSTATE rows; Zhao et al. 2005 for S).
     !--- cxs_k(:,j) is the forward rate X^it + H+ -> X^(it+1)(term j) + H
     !--- on the grid cxs_T; cxs_g, cxs_E the term's weight and energy [eV];
     !--- cxs_lev the first and last level of the term in the ion's n-level
     !--- atom; cxs_greact the weight of the X^it reactant; cxs_dip =
     !--- IP(H) - IP(X^it) [eV], the energy released into the ground term.
     integer       :: ncxs_term(MAX_ST) = 0, ncxs_T(MAX_ST) = 0
     real(kind=wp) :: cxs_greact(MAX_ST) = 0.0_wp, cxs_dip(MAX_ST) = 0.0_wp
     real(kind=wp) :: cxs_T(MAX_CXS_T,MAX_ST) = 0.0_wp
     real(kind=wp) :: cxs_k(MAX_CXS_T,MAX_CXS_TERM,MAX_ST) = 0.0_wp
     real(kind=wp) :: cxs_g(MAX_CXS_TERM,MAX_ST) = 0.0_wp
     real(kind=wp) :: cxs_E(MAX_CXS_TERM,MAX_ST) = 0.0_wp
     integer       :: cxs_lev(2,MAX_CXS_TERM,MAX_ST) = 0
     !--- CXHEAT row: kinetic energy [eV] one event gives the gas,
     !--- cxheat(1) for X^it + H+ -> X^(it+1) + H0, cxheat(2) for
     !--- X^(it+1) + H0 -> X^it + H+, and cxheat(3) > 0 when the ionization
     !--- rate includes radiative charge transfer (the statistical-weight
     !--- ratio of its non-radiative part; metal_charge_exchange_heating).
     logical       :: has_cxheat(MAX_ST) = .false.
     real(kind=wp) :: cxheat(3,MAX_ST) = 0.0_wp
     !--- CXHEATTAB: dE_ion(T) of an ionization whose product channels
     !--- branch with temperature (Mg + H+, Allan et al. 1988); nonzero
     !--- ncxheat_tab replaces cxheat(1) by linear interpolation in T
     integer       :: ncxheat_tab(MAX_ST) = 0
     real(kind=wp) :: cxheat_tab_T(MAX_CX_TAB,MAX_ST) = 0.0_wp
     real(kind=wp) :: cxheat_tab_e(MAX_CX_TAB,MAX_ST) = 0.0_wp
     !--- cooling fits per stage (loaded when the file exists)
     logical :: has_cool(MAX_ST) = .false.
     type(cooling_fit_type) :: cool(MAX_ST)
     !--- index into the nlevel_cooling_mod (T, n_e) table for each stage
     !--- (par%cooling_model = 'local_ne'); 0 = no table, keep the fit.
     integer :: cool_idx(MAX_ST) = 0
     !--- the stage's n-level atom carries atomic-hydrogen impact rates
     !--- (hatom); used by metal_cooling_H.
     logical :: has_hcoll(MAX_ST) = .false.
  end type element_type

  integer :: n_elements = 0
  type(element_type), target :: elems(MAX_EL)
  !--- n-level atoms of the stages with atomic-hydrogen impact rates (C I,
  !--- C II, O I); allocated only for those stages.
  type hcoll_atom_box
     type(nlevel_atom_type), allocatable :: a
  end type hcoll_atom_box
  type(hcoll_atom_box) :: hatom(MAX_EL, MAX_ST)

  !--- Gamma_i(transition, leaf) per element [s^-1]; filled each iteration.
  !--- egam, eheat and stage_frac hold the same values on every rank, so
  !--- they live in node-shared windows (one copy per node, 2026-10-12;
  !--- they were one copy per rank, 3.5 GB each at 2.1M leaves with the 77
  !--- stages).  Each node's ranks fill disjoint leaves (h_rank + 1, step
  !--- h_size) and species_shared_sync publishes them; the transport path
  !--- estimators metal_rate_path below are accumulated separately by every
  !--- rank and stay private.
  type gamma_block
     real(kind=wp), pointer, contiguous :: g(:,:) => null()
     integer :: win = MPI_WIN_NULL
  end type gamma_block
  type(gamma_block) :: egam(MAX_EL)
  !--- photoheating per particle [erg/s] (transition, leaf); filled with
  !--- the Gamma's, consumed by metal_heating (par%metal_heat).
  type(gamma_block) :: eheat(MAX_EL)
  !--- Metal photoionization/photoheating path estimators: unnormalized
  !--- Sum(path_lum * coeff(E)) accumulated at each packet's exact energy,
  !--- raw(1,:,:) for sigma/(h nu) and raw(2,:,:) for the heating coefficient
  !--- sum_s sigma_s (1 - E_th,s/E) (species_cross_section).
  !--- species_gamma_compute converts them into egam/eheat.
  type metal_rate_path_block
     real(kind=wp), allocatable :: raw(:,:,:)    ! (2,ntransition,nleaf)
  end type metal_rate_path_block
  type(metal_rate_path_block) :: metal_rate_path(MAX_EL)
  type(gamma_block) :: stage_frac(MAX_EL)  ! (stage,leaf), refreshed per gas state
  integer :: flat_ntransition = 0
  integer :: flat_ie(MAX_METAL_TRANSITIONS) = 0
  integer :: flat_it(MAX_METAL_TRANSITIONS) = 0
  real(kind=wp) :: flat_eth(MAX_METAL_TRANSITIONS) = 0.0_wp

  !--- rate cache for the n_e fixed point at one (leaf, T): all the
  !--- transition coefficients are functions of T (and the leaf Gammas)
  !--- only, so species_ne_prepare evaluates them ONCE and
  !--- species_metal_sums reduces each fixed-point iteration to
  !--- multiply-adds.  Without this, solve_ion_cell with par%metal_ne
  !--- re-evaluated every Badnell/Voronov/CX fit up to 200 times per
  !--- trial temperature (the 185-minute L6 PDR run).
  real(kind=wp) :: cch_gam(MAX_ST, MAX_EL), cch_ci(MAX_ST, MAX_EL)
  real(kind=wp) :: cch_rec(MAX_ST, MAX_EL), cch_cxi(MAX_ST, MAX_EL)
  real(kind=wp) :: cch_cxr(MAX_ST, MAX_EL)
  !--- form 9: the reverse rate of each product term at the cached T (its
  !--- population, which depends on n_e, is applied in species_metal_sums)
  real(kind=wp) :: cch_cxs(MAX_CXS_TERM, MAX_ST, MAX_EL), cch_T = 0.0_wp

  !--- form 9: term populations of the product ion on the (log T, log n_e)
  !--- grid, from its n-level atom at setup
  type cxs_pop_type
     real(kind=wp), allocatable :: f(:,:,:)      ! ln pop (CXS_NT, CXS_NNE, term)
  end type cxs_pop_type
  type(cxs_pop_type) :: cxs_pop(MAX_ST, MAX_EL)

contains

  character(len=2) function elem_name(ie);  integer, intent(in) :: ie
    elem_name = elems(ie)%name;  end function elem_name
  integer function elem_atomic_number(ie);  integer, intent(in) :: ie
    elem_atomic_number = elems(ie)%Z;  end function elem_atomic_number
  integer function elem_nstage(ie);  integer, intent(in) :: ie
    elem_nstage = elems(ie)%nstage;  end function elem_nstage
  real(kind=wp) function elem_abund(ie);  integer, intent(in) :: ie
    elem_abund = elems(ie)%abund;  end function elem_abund
  !--- photoionization threshold [eV] of transition it (stage it -> it+1),
  !--- it = 1..elem_nstage(ie)-1.  Exposed for the threshold-aligned band.
  real(kind=wp) function elem_eth(ie, it);  integer, intent(in) :: ie, it
    elem_eth = elems(ie)%eth(it);  end function elem_eth
  !--- the charge-exchange form of transition it (9: resolved by product term)
  integer function elem_cx_form(ie, it);  integer, intent(in) :: ie, it
    elem_cx_form = elems(ie)%cxr_form(it);  end function elem_cx_form

  ! A transported band that reaches above the ionization potential of an
  ! element's highest carried stage would photoionize that stage into one
  ! the registry does not have; the stage chain cannot represent it, so the
  ! run is refused, naming the element, its terminal ionization potential and
  ! the share of the ionizing luminosity (above 13.598 eV) that lies above it
  ! (band bins; a bin straddling the potential split linearly in energy).
  ! With stages I-VII every terminal potential is above 108 eV (Ar VII 124.4,
  ! Cl VII 114.2, Fe VII 125.0 eV; C VII is bare and has none).
  subroutine species_terminal_photoionization_check(upper_energy, eedge, lum)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(kind=wp), intent(in) :: upper_energy
    real(kind=wp), intent(in) :: eedge(:), lum(:)
    real(kind=wp) :: terminal_energy, l_above, l_ionizing, share
    character(len=512) :: message
    integer :: ie, k
    if (.not. ieee_is_finite(upper_energy) .or. upper_energy <= 0.0_wp) &
         call fatal_error('invalid transported upper photon energy')
    l_ionizing = band_luminosity_above(eth_HI)
    do ie = 1, n_elements
       if (elems(ie)%abund <= 0.0_wp) cycle
       ! A fully stripped ion has no electron and cannot be photoionized.
       if (elems(ie)%nstage == elems(ie)%Z + 1) cycle
       terminal_energy = metal_ionization_potential(elems(ie)%Z, elems(ie)%nstage)
       if (terminal_energy <= 0.0_wp) &
            call fatal_error('unknown terminal ionization potential for '//trim(elems(ie)%name))
       if (upper_energy <= terminal_energy) cycle
       l_above = band_luminosity_above(terminal_energy)
       share = 0.0_wp
       if (l_ionizing > 0.0_wp) share = l_above/l_ionizing
       write(message,'(3a,i0,a,f8.3,a,f8.3,a,es10.3,a)') &
          'the band reaches above the terminal stage of ', trim(elems(ie)%name), &
          ' (stage ', elems(ie)%nstage, ', ionization potential ', terminal_energy, &
          ' eV < eion_max ', upper_energy, ' eV; ', share, &
          ' of the ionizing luminosity lies above it); lower par%eion_max'
       call fatal_error(trim(message))
    end do
  contains
    real(kind=wp) function band_luminosity_above(energy) result(l)
      real(kind=wp), intent(in) :: energy
      l = 0.0_wp
      do k = 1, size(lum)
         if (eedge(k) >= energy) then
            l = l + lum(k)
         else if (eedge(k+1) > energy) then
            l = l + lum(k)*(eedge(k+1) - energy)/(eedge(k+1) - eedge(k))
         end if
      end do
    end function band_luminosity_above
  end subroutine species_terminal_photoionization_check

  !=========================================================================
  subroutine species_setup(nleaf)
    use mpi
    implicit none
    integer, intent(in) :: nleaf
    character(len=8)   :: names(11)
    integer, parameter :: atomic_numbers(MAX_EL) = [6,7,8,10,16,18,12,26,14,17,20]
    real(kind=wp)      :: abunds(11)
    character(len=PATH_LEN) :: fname
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
       if (trim(elems(ie)%name) /= trim(names(k)) .or. elems(ie)%Z /= atomic_numbers(k)) &
          call fatal_error('element name or atomic number disagrees with registered species '//trim(names(k)))
       elems(ie)%abund = abunds(k)
       call cx_state_setup(ie, trim(names(k)))
       !--- cooling fits for each stage where a file exists.
       do i = 1, elems(ie)%nstage
          fname=metal_excitation_data_file('cooling',trim(names(k)),i)
          inquire(file=trim(fname), exist=exists)
          if (exists) then
             call cooling_load(trim(fname), elems(ie)%cool(i))
             elems(ie)%has_cool(i) = .true.
          else
             if (mpar%p_rank == 0) write(*,'(3a,i0,a)') &
                ' SPEC: ', trim(names(k)), ' stage ', i, ' has no line-cooling coefficient file'
          end if
          ! A supplied level atom computes absolute line cooling even without
          ! an exponential fit. A missing fit must not discard its physics.
          ! Existing fitted cooling keeps its previous combination rule.
          if (trim(par%cooling_model) == 'local_ne' .or. .not. exists) &
             elems(ie)%cool_idx(i) = nlevel_cooling_add(trim(names(k)), i)
          if (.not. exists .and. elems(ie)%cool_idx(i) > 0 .and. mpar%p_rank == 0) &
             write(*,'(3a,i0,a)') ' SPEC: ', trim(names(k)), ' stage ', i, &
                  ' uses absolute n-level cooling without an exponential fit'
          !--- atomic-hydrogen impact: the stage's n-level atom, kept when
          !--- nlevel_load found hydrogen rates beside it.
          block
            type(nlevel_atom_type), allocatable :: atom
            logical :: ok
            allocate(atom)
            fname=metal_excitation_data_file('nlevel',trim(names(k)),i)
            call nlevel_load(trim(fname), atom, ok)
            if (ok .and. nlevel_has_hcoll(atom) .and. par%neutral_hydrogen_excitation) then
               call move_alloc(atom, hatom(ie,i)%a)
               elems(ie)%has_hcoll(i) = .true.
            end if
          end block
       end do
       call species_shared_create(ie, nleaf)
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
       call destroy_shared_mem(egam(ie)%g)
       call destroy_shared_mem(eheat(ie)%g)
       call destroy_shared_mem(stage_frac(ie)%g)
       call species_shared_create(ie, nleaf)
    end do
  end subroutine species_resize

  !--- the node-shared windows of element ie (zeroed by create_shared_mem)
  subroutine species_shared_create(ie, nleaf)
    integer, intent(in) :: ie, nleaf
    call create_shared_mem(egam(ie)%g, [elems(ie)%nstage-1, nleaf])
    call get_window(egam(ie)%g, egam(ie)%win)
    call create_shared_mem(eheat(ie)%g, [elems(ie)%nstage-1, nleaf])
    call get_window(eheat(ie)%g, eheat(ie)%win)
    call create_shared_mem(stage_frac(ie)%g, [elems(ie)%nstage, nleaf])
    call get_window(stage_frac(ie)%g, stage_frac(ie)%win)
  end subroutine species_shared_create

  !--- publish direct stores to the shared windows (collective on the node,
  !--- the gas_state_sync pattern)
  subroutine species_shared_sync()
    integer :: ie, ierr
    do ie = 1, n_elements
       if (egam(ie)%win /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, egam(ie)%win, ierr)
       if (eheat(ie)%win /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, eheat(ie)%win, ierr)
       if (stage_frac(ie)%win /= MPI_WIN_NULL) &
          call MPI_WIN_FENCE(0, stage_frac(ie)%win, ierr)
    end do
    call MPI_BARRIER(shared_host_comm(), ierr)
  end subroutine species_shared_sync

  !--- Metal photoionization rates of element ie on every leaf, g(transition,
  !--- leaf) [s^-1], as a previous run wrote them (previous_state_mod) or
  !--- carried onto a re-refined grid (amr_refine_mod); the
  !--- next stage-cache refresh builds the fractions the opacity uses from
  !--- them.  Overwritten by the first transport pass.  Collective on the
  !--- node; every rank passes the same g, the node's rank 0 stores it.
  subroutine species_set_photoionization_rates(ie, g)
    integer,       intent(in) :: ie
    real(kind=wp), intent(in) :: g(:,:)
    integer :: ierr
    if (shared_host_rank() == 0) egam(ie)%g(:,:) = g
    if (egam(ie)%win /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, egam(ie)%win, ierr)
    call MPI_BARRIER(shared_host_comm(), ierr)
  end subroutine species_set_photoionization_rates

  !--- the same rates as they stand (amr_refine_mod carries them onto a
  !--- re-refined grid)
  subroutine species_photoionization_rates(ie, g)
    integer,       intent(in)  :: ie
    real(kind=wp), intent(out) :: g(:,:)
    g = egam(ie)%g(:,:)
  end subroutine species_photoionization_rates

  !--- Stage fractions read by the packet opacity.  A vacuum leaf (n_H = 0)
  !--- carries the same convention as the written state (all in the first
  !--- stage): there species_fractions has no recombination (n_e = 0) against
  !--- a nonzero photoionization rate, its ratio chain overflows to Inf and
  !--- the normalized fractions come out NaN.  Until 2026-09-25 this cache
  !--- took them anyway; the opacity n_H(...+ metals) was then 0 x NaN in
  !--- every vacuum leaf, which the ray walk passed as unattenuated -- right
  !--- for vacuum, so the tallies were unaffected -- but which turned the
  !--- accumulated optical depth NaN and made the energy ledger book every
  !--- packet through a hollow as escaped (97% at HII20).
  subroutine species_stage_cache_refresh()
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHII, gas_ne, gas_Te
    implicit none
    real(kind=wp) :: frac(MAX_ST), nHI, nHII
    integer :: ie, il
    do ie = 1, n_elements
       do il = shared_host_rank() + 1, size(stage_frac(ie)%g,2), shared_host_size()
          if (gas_nH(il) > 0.0_wp) then
             nHI = gas_nH(il)*gas_xHI(il)
             nHII = gas_nH(il)*gas_xHII(il)
             call species_fractions(ie, il, gas_Te(il), gas_ne(il), &
                                    nHI, nHII, frac)
          else
             frac = 0.0_wp;  frac(1) = 1.0_wp
          end if
          stage_frac(ie)%g(:,il) = frac(1:elems(ie)%nstage)
       end do
    end do
    call species_shared_sync()
  end subroutine species_stage_cache_refresh

  !--- cross sections (and, optionally, photoheating coefficients) of every
  !--- active transition at one packet energy, ordered by threshold.
  subroutine species_packet_cross_sections(energy, sigma, ntransition, heat)
    real(kind=wp), intent(in)  :: energy
    real(kind=wp), intent(out) :: sigma(MAX_METAL_TRANSITIONS)
    integer,       intent(out) :: ntransition
    real(kind=wp), intent(out), optional :: heat(MAX_METAL_TRANSITIONS)
    real(kind=wp) :: h
    integer :: k
    sigma = 0.0_wp
    if (present(heat)) heat = 0.0_wp
    ntransition = 0
    do k = 1, flat_ntransition
       if (energy < flat_eth(k)) exit
       call species_cross_section(flat_ie(k), flat_it(k), energy, sigma(k), h)
       if (present(heat)) heat(k) = h
       ntransition = k
    end do
  end subroutine species_packet_cross_sections

  !=========================================================================
  ! Metal photoionization absorption of a packet (par%ion_metal_abs), per unit
  ! n_H and per code length: Sum_i abund frac_i sigma_i(E) (species_cross_
  ! section: the VFKY96 outer shell, or the VY95 subshell sum above the first
  ! inner-shell edge), with the cross sections already evaluated at the packet's exact energy by
  ! species_packet_cross_sections.  Small next to H/He above 13.6 eV, but NOT
  ! negligible for the He front: removing it moves r(He^0 = 0.5) inward by
  ! 0.615 pc at HII40 by softening the He-ionizing field (see the
  ! par%ion_metal_abs comment in define.f90 for the mechanism and the measured
  ! derivatives).  With par%add_fuv it is also the only GAS opacity below the
  ! Lyman limit (Mg I 7.65, C I 11.26, S I 10.36, Fe I 7.90 eV thresholds).
  ! Stage fractions come from the current state (same product chain as the
  ! cooling and the output), so the opacity feedback iterates them exactly like
  ! x_HI.
  !=========================================================================
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
  ! Metal photoionization/photoheating path estimators.  These routines
  ! intentionally live in species_mod so they can use the private registry and
  ! fill egam/eheat without exposing implementation storage.
  !=========================================================================
  subroutine species_rate_path_setup(nleaf)
    implicit none
    integer, intent(in) :: nleaf
    integer :: ie, nt

    do ie = 1, n_elements
       if (allocated(metal_rate_path(ie)%raw)) deallocate(metal_rate_path(ie)%raw)
       nt = elems(ie)%nstage - 1
       allocate(metal_rate_path(ie)%raw(2,nt,nleaf))
       metal_rate_path(ie)%raw = 0.0_wp
    end do
  end subroutine species_rate_path_setup

  subroutine species_rate_path_reset()
    integer :: ie
    do ie = 1, n_elements
       if (allocated(metal_rate_path(ie)%raw)) metal_rate_path(ie)%raw = 0.0_wp
    end do
  end subroutine species_rate_path_reset


  subroutine species_rate_path_score(energy, sigma, heat, ntransition, il, path_lum)
    real(kind=wp), intent(in) :: energy
    real(kind=wp), intent(in) :: sigma(MAX_METAL_TRANSITIONS)
    real(kind=wp), intent(in) :: heat(MAX_METAL_TRANSITIONS)
    integer,       intent(in) :: ntransition, il
    real(kind=wp), intent(in) :: path_lum
    integer :: k, ie, it
    do k = 1, ntransition
       ie = flat_ie(k)
       it = flat_it(k)
       metal_rate_path(ie)%raw(1,it,il) = metal_rate_path(ie)%raw(1,it,il) &
          + path_lum*sigma(k)/(energy*ev2erg)
       metal_rate_path(ie)%raw(2,it,il) = metal_rate_path(ie)%raw(2,it,il) &
          + path_lum*heat(k)
    end do
  end subroutine species_rate_path_score

  subroutine species_rate_path_reduce()
    use mpi
    implicit none
    integer :: ie, ierr
    do ie = 1, n_elements
       call MPI_ALLREDUCE(MPI_IN_PLACE, metal_rate_path(ie)%raw, &
                          size(metal_rate_path(ie)%raw), MPI_DOUBLE_PRECISION, &
                          MPI_SUM, MPI_COMM_WORLD, ierr)
    end do
  end subroutine species_rate_path_reduce


  !=========================================================================
  subroutine read_element(fname, el)
    use mpi
    use utility, only : fatal_error
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    character(len=*),   intent(in)  :: fname
    type(element_type), intent(out) :: el
    character(len=512) :: line, key
    character(len=8)   :: label
    integer :: unit, ios, it, ierr, n, k, ncxs_seen, next_transition, record_column
    logical :: declared, required(7,MAX_ST)

    open(newunit=unit, file=fname, status='old', iostat=ios)
    if (ios /= 0) then
       call fatal_error('cannot open '//trim(fname))
    end if
    it = 0
    next_transition = 1
    declared = .false.
    required = .false.
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*) key
       select case (trim(key))
       case ('ETH', 'PHOTO', 'PHINNER', 'PHSHELL', 'PHTABLE', 'CI', 'RR', 'DR', &
             'DR2', 'RR2', 'CXI', 'CXR', 'CXRTAB', 'CXSTATE', &
             'CXSTATET', 'CXSTATEK', 'CXHEAT', 'CXHEATTAB')
          if (it < 1 .or. it >= el%nstage) &
             call fatal_error('atomic rate record without a valid TRANSITION in '//trim(fname))
       end select
       record_column = 0
       select case (trim(key))
       case ('ETH'); record_column = 1
       case ('CI'); record_column = 2
       case ('RR', 'RR2'); record_column = 3
       case ('DR', 'DR2'); record_column = 4
       case ('PHINNER'); record_column = 5
       case ('CXI'); record_column = 6
       case ('CXR'); record_column = 7
       end select
       if (record_column > 0) then
          if (required(record_column,it)) call fatal_error('duplicate atomic rate record in '//trim(fname))
          required(record_column,it) = .true.
       end if
       select case (trim(key))
       case ('ELEMENT')
          if (declared) call fatal_error('duplicate ELEMENT record in '//trim(fname))
          declared = .true.
          read(line,*) key, el%name, el%Z, el%nstage
          if (el%Z < 1 .or. el%nstage < 1 .or. &
              el%nstage > min(MAX_ST, el%Z + 1)) &
             call fatal_error('ELEMENT stage count out of range in '//trim(fname))
       case ('TRANSITION')
          read(line,*) key, it
          if (it < 1 .or. it >= el%nstage) &
             call fatal_error('TRANSITION index out of range in '//trim(fname))
          if (it /= next_transition) call fatal_error('TRANSITION records must be consecutive in '//trim(fname))
          next_transition = next_transition + 1
       case ('ETH')
          read(line,*) key, el%eth(it)
       case ('PHOTO')
          read(line,*) key, el%photo(1:8,it)
          el%has_ph2(it) = .true.
       case ('PHINNER')
          read(line,*) key, el%einn(it)
       case ('PHSHELL')
          if (el%nshell(it) >= MAX_SHELL) &
             call fatal_error('more than 7 PHSHELL rows in '//trim(fname))
          el%nshell(it) = el%nshell(it) + 1
          read(line,*) key, el%shell(1:7,el%nshell(it),it)
       case ('PHTABLE')
          !--- n rows "eps[eV] sigma[Mb]" follow, one per line
          if (el%nphtab(it) > 0) call fatal_error('duplicate PHTABLE in '//trim(fname))
          read(line,*,iostat=ios) key, n
          if (ios /= 0 .or. n < 2 .or. n > MAX_PHTABLE) &
             call fatal_error('PHTABLE size out of range in '//trim(fname))
          allocate(el%phtab(it)%eps(n), el%phtab(it)%sig(n))
          do k = 1, n
             read(unit,*,iostat=ios) el%phtab(it)%eps(k), el%phtab(it)%sig(k)
             if (ios /= 0) call fatal_error('PHTABLE shorter than its count or '// &
                'with a malformed row in '//trim(fname))
          end do
          if (.not. all(ieee_is_finite(el%phtab(it)%eps)) .or. &
              .not. all(ieee_is_finite(el%phtab(it)%sig))) &
             call fatal_error('nonfinite PHTABLE values in '//trim(fname))
          if (el%phtab(it)%eps(1) /= 0.0_wp .or. &
              any(el%phtab(it)%eps(2:n) <= el%phtab(it)%eps(1:n-1))) &
             call fatal_error('PHTABLE photoelectron energies must start at 0 '// &
                'and increase in '//trim(fname))
          if (any(el%phtab(it)%sig < 0.0_wp)) &
             call fatal_error('negative PHTABLE cross section in '//trim(fname))
          el%phtab(it)%sig = el%phtab(it)%sig*1.0e-18_wp   ! Mb -> cm^2
          el%nphtab(it) = n
       case ('CI')
          read(line,*) key, el%ci(1:5,it)
       case ('RR')
          read(line,*) key, el%rr(1:6,it)
       case ('DR')
          read(line,*) key, n
          if (n < 0 .or. n > size(el%dr_c, 1)) &
             call fatal_error('DR term count out of range in '//trim(fname))
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
       case ('CXRTAB')
          read(line,*) key, n
          if (n < 2 .or. n > MAX_CXR_TAB) &
             call fatal_error('CXRTAB size out of range in '//trim(fname))
          el%ncxr_tab(it) = n
          backspace(unit)
          read(unit,*) key, n, el%cxr_tab_T(1:n,it), el%cxr_tab_k(1:n,it)
          if (.not. all(ieee_is_finite(el%cxr_tab_T(1:n,it))) .or. &
              .not. all(ieee_is_finite(el%cxr_tab_k(1:n,it)))) &
               call fatal_error('nonfinite CXRTAB values in '//trim(fname))
          if (any(el%cxr_tab_T(1:n,it) <= 0.0_wp) .or. &
              any(el%cxr_tab_k(1:n,it) <= 0.0_wp)) &
               call fatal_error('nonpositive CXRTAB values in '//trim(fname))
          if (any(el%cxr_tab_T(2:n,it) <= el%cxr_tab_T(1:n-1,it))) &
               call fatal_error('nonincreasing CXRTAB temperature in '//trim(fname))
       case ('CXSTATE')
          read(line,*) key, el%ncxs_term(it), el%ncxs_T(it), &
                       el%cxs_greact(it), el%cxs_dip(it)
          if (el%ncxs_term(it) < 1 .or. el%ncxs_term(it) > MAX_CXS_TERM .or. &
              el%ncxs_T(it) < 2 .or. el%ncxs_T(it) > MAX_CXS_T) &
             call fatal_error('CXSTATE size out of range in '//trim(fname))
          ncxs_seen = 0
       case ('CXSTATET')
          backspace(unit)
          read(unit,*) key, el%cxs_T(1:el%ncxs_T(it),it)
       case ('CXSTATEK')
          ncxs_seen = ncxs_seen + 1
          if (ncxs_seen > el%ncxs_term(it)) &
             call fatal_error('more CXSTATEK rows than CXSTATE terms in '// &
                              trim(fname))
          backspace(unit)
          read(unit,*) key, label, el%cxs_lev(1:2,ncxs_seen,it), &
                       el%cxs_g(ncxs_seen,it), el%cxs_E(ncxs_seen,it), &
                       el%cxs_k(1:el%ncxs_T(it),ncxs_seen,it)
       case ('CXHEAT')
          read(line,*) key, el%cxheat(1:3,it)
          el%has_cxheat(it) = .true.
       case ('CXHEATTAB')
          read(line,*) key, n
          if (n < 2 .or. n > MAX_CX_TAB) &
             call fatal_error('CXHEATTAB size out of range in '//trim(fname))
          !--- the row can outgrow the line buffer (n pairs of reals), so it
          !--- is read again from the unit, as CXRTAB is
          backspace(unit)
          read(unit,*) key, n, (el%cxheat_tab_T(k,it), el%cxheat_tab_e(k,it), k = 1, n)
          el%ncxheat_tab(it) = n
       case default
          call fatal_error('unknown atomic record '//trim(key)//' in '//trim(fname))
       end select
    end do
    close(unit)
    if (.not. declared) call fatal_error('missing ELEMENT record in '//trim(fname))
    if (next_transition /= el%nstage .or. .not. all(required(:,1:el%nstage-1))) &
       call fatal_error('incomplete atomic transition records in '//trim(fname))
    if (.not. all(ieee_is_finite(el%eth)) .or. .not. all(ieee_is_finite(el%photo)) .or. &
        .not. all(ieee_is_finite(el%einn)) .or. .not. all(ieee_is_finite(el%shell)) .or. &
        .not. all(ieee_is_finite(el%ci)) .or. .not. all(ieee_is_finite(el%rr)) .or. &
        .not. all(ieee_is_finite(el%rr2)) .or. .not. all(ieee_is_finite(el%dr_c)) .or. &
        .not. all(ieee_is_finite(el%dr_e)) .or. .not. all(ieee_is_finite(el%dr2)) .or. &
        .not. all(ieee_is_finite(el%cxi)) .or. .not. all(ieee_is_finite(el%cxr)) .or. &
        .not. all(ieee_is_finite(el%cxheat))) &
       call fatal_error('nonfinite atomic rate coefficients in '//trim(fname))
    !--- every transition needs its subshell fits; a file written before
    !--- they were added (PHOTO2 or PHOTO alone) would silently leave the
    !--- inner shells out, so refuse it and point at the generator.
    do it = 1, el%nstage - 1
       if (el%eth(it) <= 0.0_wp .or. el%ci(1,it) <= 0.0_wp .or. &
           el%ci(3,it) < 0.0_wp .or. el%ci(4,it) < 0.0_wp .or. el%einn(it) < 0.0_wp) &
          call fatal_error('invalid atomic threshold or CI coefficients in '//trim(fname))
       if (el%has_rr2(it)) then
          if (el%rr2(1,it) <= 0.0_wp) call fatal_error('nonpositive RR2 coefficient in '//trim(fname))
       else
          if (el%rr(1,it) <= 0.0_wp .or. el%rr(3,it) <= 0.0_wp .or. el%rr(4,it) <= 0.0_wp) &
             call fatal_error('nonpositive RR coefficient or temperature in '//trim(fname))
       end if
       if (el%nshell(it) == 0) call fatal_error(trim(fname)// &
          ' has no PHSHELL rows: regenerate it with '// &
          'tools/fitting/make_element_data.py')
       if (el%einn(it) > 0.0_wp .and. .not. el%has_ph2(it)) &
          call fatal_error(trim(fname)//': PHINNER > 0 without a PHOTO row')
       !--- a PHTABLE replaces the PHOTO fit below its last node, so it needs
       !--- the fit, must end below the first inner edge, and must end ON the
       !--- fit (continuity at the join; the generator writes the fit's value)
       if (el%nphtab(it) > 0) then
          if (.not. el%has_ph2(it)) &
             call fatal_error(trim(fname)//': PHTABLE without a PHOTO row')
          n = el%nphtab(it)
          if (el%eth(it) + el%phtab(it)%eps(n) >= el%einn(it)) &
             call fatal_error(trim(fname)//': PHTABLE reaches the first inner-shell edge')
          if (abs(el%phtab(it)%sig(n)/max(sigma_vfky96(el%eth(it) + el%phtab(it)%eps(n), &
                  el%photo(1,it), el%photo(2,it), el%photo(3,it), el%photo(4,it), &
                  el%photo(5,it), el%photo(6,it), el%photo(7,it), el%photo(8,it)), &
                  tinest) - 1.0_wp) > 1.0e-8_wp) &
             call fatal_error(trim(fname)//': PHTABLE does not end on the PHOTO fit')
       end if
       if (el%cxr_form(it) == 8 .and. el%ncxr_tab(it) < 2) &
          call fatal_error(trim(fname)//': CXR form 8 without a CXRTAB row')
       if ((el%cxi_form(it) == 9 .or. el%cxr_form(it) == 9) .and. &
           (el%cxi_form(it) /= el%cxr_form(it) .or. el%ncxs_term(it) < 1 .or. &
            any(el%cxs_k(1:el%ncxs_T(it),1:el%ncxs_term(it),it) <= 0.0_wp))) &
          call fatal_error(trim(fname)//': CX form 9 needs both directions '// &
                           'and complete CXSTATE rows')
       !--- charge-exchange heating needs the energy of every carried
       !--- reaction; a file written before the CXHEAT rows existed would
       !--- leave it silently at zero.
       if (par%charge_exchange_heating .and. .not. el%has_cxheat(it)) &
          call fatal_error(trim(fname)//' has no CXHEAT row for a '// &
             'transition (par%charge_exchange_heating): regenerate it '// &
             'with tools/fitting/make_element_data.py')
    end do
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
  !--- charge-exchange rate coefficients of transition it of element ie at T:
  !--- ionization X^it + H+ -> X^(it+1) + H0 and recombination
  !--- X^(it+1) + H0 -> X^it + H+.  The recombination may be tabulated.
  real(kind=wp) function cxi_rate(ie, it, T) result(k)
    integer,       intent(in) :: ie, it
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: kf(MAX_CXS_TERM)
    k = 0.0_wp
    if (it >= 5 .and. .not. par%langevin_charge_exchange) return
    if (elems(ie)%cxi_form(it) == 9) then
       call cxs_forward_terms(ie, it, T, kf)
       k = sum(kf(1:elems(ie)%ncxs_term(it)))
       return
    end if
    if (elems(ie)%Z == 7 .and. it == 1) then
       ! KF96 Table 3: N + H+ base fit is valid from 100 to 50000 K.
       k = cx_rate(elems(ie)%cxi_form(it), elems(ie)%cxi(1:6,it), T, 100.0_wp)
    else
       k = cx_rate(elems(ie)%cxi_form(it), elems(ie)%cxi(1:6,it), T)
    end if
  end function cxi_rate

  !--- ne: the electron density that sets the product-term populations of a
  !--- form-9 reverse rate; left out, their n_e -> 0 limit (the ground term)
  real(kind=wp) function cxr_rate(ie, it, T, ne) result(k)
    integer,       intent(in) :: ie, it
    real(kind=wp), intent(in) :: T
    real(kind=wp), intent(in), optional :: ne
    real(kind=wp) :: lt, w, kr(MAX_CXS_TERM), f(MAX_CXS_TERM)
    integer :: n, j
    k = 0.0_wp
    if (it >= 5 .and. .not. par%langevin_charge_exchange) return
    if (elems(ie)%cxr_form(it) == 9) then
       call cxs_reverse_terms(ie, it, T, kr)
       if (present(ne)) then
          call cxs_populations(ie, it, T, ne, f)
       else
          call cxs_populations(ie, it, T, 0.0_wp, f)
       end if
       n = elems(ie)%ncxs_term(it)
       k = sum(f(1:n)*kr(1:n))
       return
    end if
    if (elems(ie)%cxr_form(it) /= 8) then
       if (elems(ie)%Z == 7 .and. it == 1) then
          ! KF96 Table 1: N+ + H has the same 100--50000 K base-fit interval.
          k = cx_rate(elems(ie)%cxr_form(it), elems(ie)%cxr(1:6,it), T, 100.0_wp)
       else
          k = cx_rate(elems(ie)%cxr_form(it), elems(ie)%cxr(1:6,it), T)
       end if
       return
    end if
    n = elems(ie)%ncxr_tab(it)
    associate (tt => elems(ie)%cxr_tab_T(1:n,it), kk => elems(ie)%cxr_tab_k(1:n,it))
      lt = log(min(max(T, tt(1)), tt(n)))
      j = 1
      do while (j < n-1 .and. lt > log(tt(j+1)))
         j = j + 1
      end do
      w = (lt - log(tt(j)))/(log(tt(j+1)) - log(tt(j)))
      k = exp(log(kk(j)) + w*(log(kk(j+1)) - log(kk(j))))
    end associate
  end function cxr_rate

  !=========================================================================
  ! Charge exchange resolved by product term (form 9).  Forward
  !   X^it + H+ -> X^(it+1)(j) + H,  k_j(T) from the table (log-log in T,
  !   clamped to it);
  ! reverse from term j of the ion by detailed balance of that channel,
  !   k_rj = k_j g_react / (g_j g_H) exp[-(dIP - E_j)/kT],  g_H = 2, g_H+ = 1,
  ! dIP - E_j being the energy the reverse reaction from j needs.  The rate
  ! of the ion is sum_j pop_j(T, n_e) k_rj with its ACTUAL term populations
  ! (cxs_populations).  A reverse formed from the summed forward rate by one
  ! detailed balance -- Huang et al. (2023, sec. 2.4.2) for S+ + H -- is
  ! this sum with the populations in LTE; at nebular densities S+ sits in
  ! 4S and that overstates the reverse by ~1e4 (2.3e-17 against 3.7e-13 at
  ! 1e4 K and n_e -> 0).
  !=========================================================================
  subroutine cxs_forward_terms(ie, it, T, kf)
    integer,       intent(in)  :: ie, it
    real(kind=wp), intent(in)  :: T
    real(kind=wp), intent(out) :: kf(MAX_CXS_TERM)
    real(kind=wp) :: lt, w
    integer :: n, j, i
    kf = 0.0_wp
    n = elems(ie)%ncxs_T(it)
    associate (tt => elems(ie)%cxs_T(1:n,it))
      lt = log(min(max(T, tt(1)), tt(n)))
      i = 1
      do while (i < n-1 .and. lt > log(tt(i+1)))
         i = i + 1
      end do
      w = (lt - log(tt(i)))/(log(tt(i+1)) - log(tt(i)))
    end associate
    do j = 1, elems(ie)%ncxs_term(it)
       kf(j) = exp(log(elems(ie)%cxs_k(i,j,it)) &
                   + w*(log(elems(ie)%cxs_k(i+1,j,it)) - log(elems(ie)%cxs_k(i,j,it))))
    end do
  end subroutine cxs_forward_terms

  subroutine cxs_reverse_terms(ie, it, T, kr)
    integer,       intent(in)  :: ie, it
    real(kind=wp), intent(in)  :: T
    real(kind=wp), intent(out) :: kr(MAX_CXS_TERM)
    real(kind=wp) :: kf(MAX_CXS_TERM), x
    integer :: j
    call cxs_forward_terms(ie, it, T, kf)
    kr = 0.0_wp
    do j = 1, elems(ie)%ncxs_term(it)
       x = (elems(ie)%cxs_dip(it) - elems(ie)%cxs_E(j,it))*ev2erg/(kboltz_cgs*T)
       kr(j) = kf(j)*elems(ie)%cxs_greact(it)/(2.0_wp*elems(ie)%cxs_g(j,it)) &
               *exp(-min(x, 600.0_wp))
    end do
  end subroutine cxs_reverse_terms

  !--- populations of the represented product-ion terms, bilinear in
  !--- (log T, log n_e) of their logarithm on the setup table. The selected
  !--- atom is normalized over all its levels; represented terms need not
  !--- sum to one when other levels exist. Do not renormalize this subset:
  !--- those other levels do not have rates in the resolved reaction table.
  !--- An excited
  !--- term goes as exp(-E/kT), which a linear interpolation of the
  !--- population itself misses by ~1% at 3000 K on a 0.02-dex grid -- T and
  !--- n_e clamped to the table (n_e below 1e-6 cm^-3 is the low-density
  !--- column, the ground term)
  subroutine cxs_populations(ie, it, T, ne, f)
    integer,       intent(in)  :: ie, it
    real(kind=wp), intent(in)  :: T, ne
    real(kind=wp), intent(out) :: f(MAX_CXS_TERM)
    real(kind=wp) :: x, y, wx, wy
    integer :: i, j, n
    n = elems(ie)%ncxs_term(it)
    x = (log10(max(T, 1.0_wp)) - CXS_T_LO)/CXS_DT
    x = min(max(x, 0.0_wp), real(CXS_NT-1, wp) - 1.0e-9_wp)
    y = (log10(max(ne, 1.0e-30_wp)) - CXS_NE_LO)/CXS_DNE
    y = min(max(y, 0.0_wp), real(CXS_NNE-1, wp) - 1.0e-9_wp)
    i = int(x) + 1;  wx = x - real(i-1, wp)
    j = int(y) + 1;  wy = y - real(j-1, wp)
    f = 0.0_wp
    associate (p => cxs_pop(it,ie)%f)
      f(1:n) = exp((1.0_wp-wx)*(1.0_wp-wy)*p(i,j,1:n) + wx*(1.0_wp-wy)*p(i+1,j,1:n) &
                   + (1.0_wp-wx)*wy*p(i,j+1,1:n) + wx*wy*p(i+1,j+1,1:n))
    end associate
  end subroutine cxs_populations

  !--- build the population table of every form-9 transition of element ie
  !--- from the n-level atom of the product ion, after checking that the
  !--- term energies of the element file are that atom's
  subroutine cx_state_setup(ie, name)
    use utility,    only : fatal_error
    use nlevel_mod, only : nlevel_cooling_prepT, nlevel_populations
    integer,          intent(in) :: ie
    character(len=*), intent(in) :: name
    type(nlevel_atom_type), allocatable :: atom
    real(kind=wp), allocatable :: Rcoll(:,:), Rrad(:,:), pop(:)
    real(kind=wp) :: sg, se, Eterm
    character(len=PATH_LEN) :: fname
    logical :: ok
    integer :: it, i, j, jt, lo, hi, n
    real(kind=wp), parameter :: HC_EV_CM = 1.239841984e-4_wp
    do it = 1, elems(ie)%nstage - 1
       if (elems(ie)%cxr_form(it) /= 9) cycle
       n = elems(ie)%ncxs_term(it)
       fname=metal_excitation_data_file('nlevel',trim(name),it+1)
       allocate(atom)
       call nlevel_load(trim(fname), atom, ok)
       if (.not. ok) call fatal_error('charge exchange form 9 needs '// &
          trim(fname)//' (the product ion''s levels)')
       do jt = 1, n
          lo = elems(ie)%cxs_lev(1,jt,it);  hi = elems(ie)%cxs_lev(2,jt,it)
          if (lo < 1 .or. hi > atom%nlev .or. hi < lo) call fatal_error( &
             'CXSTATEK levels outside '//trim(fname))
          sg = sum(atom%g(lo:hi))
          se = sum(atom%g(lo:hi)*atom%e_cm1(lo:hi))
          Eterm = se/sg*HC_EV_CM
          if (abs(Eterm - elems(ie)%cxs_E(jt,it)) > 2.0e-3_wp .or. &
              abs(sg - elems(ie)%cxs_g(jt,it)) > 1.0e-6_wp) call fatal_error( &
             'CXSTATEK term energy or weight disagrees with '//trim(fname))
       end do
       allocate(cxs_pop(it,ie)%f(CXS_NT, CXS_NNE, n))
       allocate(Rcoll(atom%nlev,atom%nlev), Rrad(atom%nlev,atom%nlev), pop(atom%nlev))
       do i = 1, CXS_NT
          call nlevel_cooling_prepT(atom, 10.0_wp**(CXS_T_LO + (i-1)*CXS_DT), &
                                    Rcoll, Rrad)
          do j = 1, CXS_NNE
             call nlevel_populations(atom, Rcoll, Rrad, &
                  10.0_wp**(CXS_NE_LO + (j-1)*CXS_DNE), pop)
             do jt = 1, n
                lo = elems(ie)%cxs_lev(1,jt,it);  hi = elems(ie)%cxs_lev(2,jt,it)
                cxs_pop(it,ie)%f(i,j,jt) = log(max(sum(pop(lo:hi)), 1.0e-300_wp))
             end do
          end do
       end do
       deallocate(Rcoll, Rrad, pop, atom)
       if (mpar%p_rank == 0) write(*,'(3a,i0,a,es10.3,a,es10.3)') &
          ' SPEC: ', trim(name), ' transition ', it, &
          ' charge exchange by product term; reverse at 1e4 K, n_e -> 0: ', &
          cxr_rate(ie, it, 1.0e4_wp), ', n_e = 1e2: ', cxr_rate(ie, it, 1.0e4_wp, 1.0e2_wp)
    end do
  end subroutine cx_state_setup

  real(kind=wp) function cx_rate(form, c, T, minimum_temperature) result(k)
    integer,       intent(in) :: form
    real(kind=wp), intent(in) :: c(6), T
    real(kind=wp), intent(in), optional :: minimum_temperature
    real(kind=wp) :: t4, lnT, fitted_t4, lower_t4
    t4 = T/1.0e4_wp
    lower_t4 = 0.5_wp
    if (present(minimum_temperature)) lower_t4 = minimum_temperature/1.0e4_wp
    fitted_t4 = min(max(t4, lower_t4), 5.0_wp)
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
       ! KF96 Mg + H+ uses 5000--50000 K; N uses 100--50000 K.
       ! Endpoint continuation is an approximation outside the fit interval.
       ! Fe's temperature-independent coefficient is unaffected by clamping.
       ! Higher boundaries use their individually resolved CXRTAB intervals.
       k = c(1)*fitted_t4**c(2)*(1.0_wp + c(3)*exp(c(4)*fitted_t4))
       k = max(k, 0.0_wp)
    case (7)
       !--- KF96 form with an activation barrier (the endothermic direction,
       !--- KF96 Table 3).  Only the fitted form-6 factor is clamped into its
       !--- validity range; the Boltzmann factor exp(-dE/kT) is exact physics
       !--- and takes the actual T.  Clamping it too (as until 2026-09-24)
       !--- froze exp(-dE/kT) at its 5000 K value below 5000 K, overstating
       !--- the endothermic rate -- 26x for N + H+ at 2000 K -- and breaking
       !--- detailed balance against the form-6 exothermic partner.
       ! Mg's inferred reverse retains the forward fit's interval; its
       ! product-term assumptions and the published reverse range are unverified.
       k = c(1)*fitted_t4**c(2) &
           *(1.0_wp + c(3)*exp(c(4)*fitted_t4)) &
           *exp(-c(5)/t4)
       k = max(k, 0.0_wp)
    case default
       k = 0.0_wp
    end select
  end function cx_rate

  !=========================================================================
  ! Metal photoionization rate and photoheating integrals per leaf, from the
  ! reduced transport-path estimators (all elements, all stages); call after
  ! ion_score_reduce each iteration.
  !=========================================================================
  subroutine species_gamma_compute()
    use octree_mod, only : amr_grid, leaf_half
    implicit none
    integer :: ie, it, il, nleaf
    real(kind=wp) :: vol, fac

    nleaf = amr_grid%nleaf
    do ie = 1, n_elements
       if (.not. allocated(metal_rate_path(ie)%raw)) &
          error stop 'metal path estimators are not allocated'
       do il = shared_host_rank() + 1, nleaf, shared_host_size()
          vol = (2.0_wp*leaf_half(il))**3
          fac = 1.0_wp/(vol*par%distance2cm**2)
          do it = 1, elems(ie)%nstage-1
             egam(ie)%g(it,il) = metal_rate_path(ie)%raw(1,it,il)*fac
             eheat(ie)%g(it,il) = metal_rate_path(ie)%raw(2,it,il)*fac
          end do
       end do
    end do
    call species_shared_sync()
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
          cch_ci(it,ie) = metal_electron_impact_ionization(elems(ie), it, T)
          cch_rec(it,ie) = alpha_rec(elems(ie), it, T)
          cch_cxi(it,ie) = 0.0_wp
          cch_cxr(it,ie) = 0.0_wp
          if (elems(ie)%cxi_form(it) > 0) cch_cxi(it,ie) = &
             cxi_rate(ie, it, T)
          if (elems(ie)%cxr_form(it) > 0 .and. elems(ie)%cxr_form(it) /= 9) &
             cch_cxr(it,ie) = cxr_rate(ie, it, T)
          if (elems(ie)%cxr_form(it) == 9) &
             call cxs_reverse_terms(ie, it, T, cch_cxs(:,it,ie))
       end do
    end do
    cch_T = T
  end subroutine species_ne_prepare

  !=========================================================================
  ! The three sums over the metal cascade that hydrogen's balance and the
  ! n_e closure need, from the same normalized adjacent-stage equilibrium
  ! calculation as species_fractions:
  !
  !   nem      = sum_el n_H A_el sum_i (i-1) f_i          [cm^-3]
  !              electrons released by the cascade (par%metal_ne);
  !   r_cx_in  = sum_el sum_it k_CXI(el,it) n(X_el^it)    [s^-1]
  !              rate at which one H^+ is neutralized by charge exchange,
  !              X^it + H^+ -> X^(it+1) + H^0;
  !   r_cx_out = sum_el sum_it k_CXR(el,it) n(X_el^(it+1)) [s^-1]
  !              rate at which one H^0 is ionized by the reverse reaction,
  !              X^(it+1) + H^0 -> X^it + H^+.
  !
  ! The two charge-exchange sums are hydrogen's side of the reactions whose
  ! metal side the chain above already carries (the k_CXI n_HII and
  ! k_CXR n_HI terms in rion/rrec).  Both sides read the same cached
  ! coefficient here, so the pair conserves charge by construction rather
  ! than by a tolerance (par%charge_exchange = 'two_way',
  ! ion_balance_mod::solve_ion_cell).
  !
  ! cx_in_stage / cx_out_stage optionally return the same two sums resolved
  ! by reaction, cx_in_stage(it,ie) = k_CXI(it,ie) n(X_el^it) and
  ! cx_out_stage(it,ie) = k_CXR(it,ie) n(X_el^(it+1)) [s^-1]: which reaction
  ! moves the charge, and the form the conservation gate compares against the
  ! cascade's own side (tests/charge_exchange).
  !=========================================================================
  subroutine species_metal_sums(nH, ne, nHI, nHII, nem, r_cx_in, r_cx_out, &
                                cx_in_stage, cx_out_stage)
    implicit none
    real(kind=wp), intent(in)  :: nH, ne, nHI, nHII
    real(kind=wp), intent(out) :: nem, r_cx_in, r_cx_out
    real(kind=wp), intent(out), optional :: cx_in_stage(MAX_ST, MAX_EL)
    real(kind=wp), intent(out), optional :: cx_out_stage(MAX_ST, MAX_EL)
    real(kind=wp) :: ionization(MAX_ST-1), recombination(MAX_ST-1)
    real(kind=wp) :: frac(MAX_ST), fr, nX, kcr(MAX_ST)
    real(kind=wp) :: fpop(MAX_CXS_TERM)
    integer :: ie, it, i, ns, nt, status
    nem = 0.0_wp
    r_cx_in  = 0.0_wp
    r_cx_out = 0.0_wp
    if (present(cx_in_stage)) cx_in_stage = 0.0_wp
    if (present(cx_out_stage)) cx_out_stage = 0.0_wp
    do ie = 1, n_elements
       ns = elems(ie)%nstage
       !--- reverse charge-exchange rates; a form-9 transition's depends on
       !--- n_e through the product-term populations, applied here
       kcr(1:ns-1) = cch_cxr(1:ns-1,ie)
       do it = 1, ns-1
          if (elems(ie)%cxr_form(it) == 9) then
             nt = elems(ie)%ncxs_term(it)
             call cxs_populations(ie, it, cch_T, ne, fpop)
             kcr(it) = sum(fpop(1:nt)*cch_cxs(1:nt,it,ie))
          end if
       end do
       do it = 1, ns-1
          ionization(it) = cch_gam(it,ie) + ne*cch_ci(it,ie) &
                           + nHII*cch_cxi(it,ie)
          recombination(it) = ne*cch_rec(it,ie) + nHI*kcr(it)
       end do
       if (ne == 0.0_wp .and. nHII == 0.0_wp .and. &
           all(cch_gam(1:ns-1,ie) == 0.0_wp)) then
          ! The same electron-free neutral limit as species_fractions.
          frac = 0.0_wp; frac(1) = 1.0_wp
       else
          call adjacent_ion_fractions(ionization(1:ns-1), recombination(1:ns-1), frac(1:ns), status)
          if (status /= 0) call ion_stage_balance_error(ie, status)
       end if
       nX = elems(ie)%abund*nH
       fr = frac(1)
       !--- fr is now the stage-1 fraction; the CXI reaction of transition 1
       !--- ionizes that stage.
       if (ns > 1) then
          r_cx_in = r_cx_in + cch_cxi(1,ie)*nX*fr
          if (present(cx_in_stage)) cx_in_stage(1,ie) = cch_cxi(1,ie)*nX*fr
       end if
       do i = 2, ns
          fr = frac(i)                             ! fraction of stage i
          nem = nem + elems(ie)%abund*nH*fr*real(i-1, wp)
          !--- stage i is the recombining partner of transition i-1 and the
          !--- ionized partner of transition i.
          r_cx_out = r_cx_out + kcr(i-1)*nX*fr
          if (i < ns) r_cx_in = r_cx_in + cch_cxi(i,ie)*nX*fr
          if (present(cx_out_stage)) cx_out_stage(i-1,ie) = kcr(i-1)*nX*fr
          if (present(cx_in_stage) .and. i < ns) cx_in_stage(i,ie) = cch_cxi(i,ie)*nX*fr
       end do
    end do
  end subroutine species_metal_sums

  !=========================================================================
  ! Charge-exchange rate coefficients [cm^3 s^-1] of transition it of element
  ! ie at temperature T, evaluated from the registry data:
  !   k_cxi  for  X^it + H^+ -> X^(it+1) + H^0,
  !   k_cxr  for  X^(it+1) + H^0 -> X^it + H^+,
  ! zero where the registry carries no determination for that direction.
  !=========================================================================
  subroutine species_cx_coefficients(ie, it, T, k_cxi, k_cxr, ne)
    implicit none
    integer,       intent(in)  :: ie, it
    real(kind=wp), intent(in)  :: T
    real(kind=wp), intent(out) :: k_cxi, k_cxr
    !--- the electron density a form-9 reverse rate depends on (its product-
    !--- term populations); left out, the n_e -> 0 limit
    real(kind=wp), intent(in), optional :: ne
    k_cxi = 0.0_wp
    k_cxr = 0.0_wp
    if (ie < 1 .or. ie > n_elements) return
    if (it < 1 .or. it > elems(ie)%nstage-1) return
    if (elems(ie)%cxi_form(it) > 0) &
       k_cxi = cxi_rate(ie, it, T)
    if (elems(ie)%cxr_form(it) > 0) &
       k_cxr = cxr_rate(ie, it, T, ne)
  end subroutine species_cx_coefficients

  !=========================================================================
  ! Kinetic energy [eV] one event of each charge-exchange direction of
  ! transition it of element ie gives the gas (the CXHEAT row):
  !   de_ion  for  X^it + H^+ -> X^(it+1) + H^0,
  !   de_rec  for  X^(it+1) + H^0 -> X^it + H^+,
  ! and g_split > 0 when the ionization rate includes radiative charge
  ! transfer (metal_charge_exchange_heating).  Zero where no rate is carried.
  ! A transition resolved by product term (form 9) returns the rate-weighted
  ! mean over the terms at (T, ne) -- the forward into term j releases
  ! dIP - E_j, the reverse from j takes it back -- which is what the heating
  ! uses; without T the CXHEAT row (its 1e4 K, n_e -> 0 values).
  !=========================================================================
  subroutine species_cx_kinetic_energy(ie, it, de_ion, de_rec, g_split, T, ne)
    implicit none
    integer,       intent(in)  :: ie, it
    real(kind=wp), intent(out) :: de_ion, de_rec, g_split
    real(kind=wp), intent(in), optional :: T   ! for a CXHEATTAB defect
    real(kind=wp), intent(in), optional :: ne  ! for a form-9 reverse
    real(kind=wp) :: kf(MAX_CXS_TERM), kr(MAX_CXS_TERM), f(MAX_CXS_TERM)
    real(kind=wp) :: de(MAX_CXS_TERM)
    integer :: nt
    de_ion = 0.0_wp;  de_rec = 0.0_wp;  g_split = 0.0_wp
    if (ie < 1 .or. ie > n_elements) return
    if (it < 1 .or. it > elems(ie)%nstage-1) return
    if (it >= 5 .and. .not. par%langevin_charge_exchange) return
    if (elems(ie)%cxi_form(it) == 9 .and. present(T)) then
       nt = elems(ie)%ncxs_term(it)
       call cxs_forward_terms(ie, it, T, kf)
       call cxs_reverse_terms(ie, it, T, kr)
       if (present(ne)) then
          call cxs_populations(ie, it, T, ne, f)
       else
          call cxs_populations(ie, it, T, 0.0_wp, f)
       end if
       de(1:nt) = elems(ie)%cxs_dip(it) - elems(ie)%cxs_E(1:nt,it)
       de_ion = sum(kf(1:nt)*de(1:nt))/sum(kf(1:nt))
       de_rec = -sum(f(1:nt)*kr(1:nt)*de(1:nt))/max(sum(f(1:nt)*kr(1:nt)), tinest)
       return
    end if
    de_ion  = elems(ie)%cxheat(1,it)
    if (present(T)) de_ion = cx_ionization_defect(ie, it, T)
    de_rec  = elems(ie)%cxheat(2,it)
    g_split = elems(ie)%cxheat(3,it)
  end subroutine species_cx_kinetic_energy

  !=========================================================================
  ! Metal electrons alone, for the n_e closure of charge_neutrality_ne.
  !=========================================================================
  real(kind=wp) function species_ne_cached(nH, ne, nHI, nHII) result(nem)
    implicit none
    real(kind=wp), intent(in) :: nH, ne, nHI, nHII
    real(kind=wp) :: r_cx_in, r_cx_out
    call species_metal_sums(nH, ne, nHI, nHII, nem, r_cx_in, r_cx_out)
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
  ! Atomic-hydrogen impact line cooling [erg s^-1 cm^-3]: the ground-term
  ! fine-structure lines of C I (609/370 um), C II (158 um) and O I (63/145
  ! um) excited by NEUTRAL HYDROGEN -- the dominant coolant of neutral gas,
  ! where n_HI/n_e ~ 10^3-10^4.  metal_cooling carries electron impact from
  ! the (T, n_e) table; this adds, for every stage whose n-level atom has
  ! hydrogen rates (data/atomic/hcoll_<el>_<stage>.txt), n_ion times the
  ! increment of the radiated power when n_HI q^H joins n_e q^e in one
  ! statistical equilibrium (nlevel_mod::nlevel_hcoll_cooling).  Both
  ! colliders therefore saturate the populations together, all ground-term
  ! transitions are included, and the level energies, A-values and
  ! populations are those of the atom that writes the line list.
  ! Rates: Abrahamsson, Krems & Dalgarno (2007) for C I and O I, Barinovs et
  ! al. (2005) for C II (make_hcoll_fine_structure.py).  Called from the
  ! par%use_metals cooling path whatever heats the gas.
  !=========================================================================
  real(kind=wp) function metal_cooling_H(il, T, nH, ne, nHI, nHII) result(cool)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne, nHI, nHII
    real(kind=wp) :: frac(MAX_ST), nion
    integer :: ie, i
    cool = 0.0_wp
    if (.not. par%neutral_hydrogen_excitation) return
    if (nHI <= 0.0_wp) return
    do ie = 1, n_elements
       if (.not. any(elems(ie)%has_hcoll(1:elems(ie)%nstage))) cycle
       call species_fractions(ie, il, T, ne, nHI, nHII, frac)
       do i = 1, elems(ie)%nstage
          if (.not. elems(ie)%has_hcoll(i)) cycle
          nion = elems(ie)%abund*nH*frac(i)
          if (nion <= 0.0_wp) cycle
          cool = cool + nion*nlevel_hcoll_cooling(hatom(ie,i)%a, T, ne, nHI)
       end do
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

  !=========================================================================
  ! Charge-exchange heating of the metals [erg s^-1 cm^-3]
  ! (par%charge_exchange_heating): every carried reaction with hydrogen,
  !   X^it + H^+ -> X^(it+1) + H^0      rate k_CXI n(X^it) n_HII,
  !   X^(it+1) + H^0 -> X^it + H^+      rate k_CXR n(X^(it+1)) n_HI,
  ! times the kinetic energy one event gives the gas,
  !   dE = [IP of the stage the electron joins - IP of the stage it
  !         leaves] - E_rad,
  ! i.e. the binding energy the electron gains (O2+ + H: 35.121 - 13.598
  ! = 21.52 eV less the O+ excitation that radiates),
  ! E_rad being the excitation of the product state that leaves as
  ! radiation.  The ionization energy the reaction moves between the two
  ! species is not thermal energy in this code (photoionization deposits
  ! h nu - IP, recombination removes the electron's kinetic energy), so this
  ! is the whole change of the gas energy; an endothermic reaction has
  ! dE < 0.  The dE of every reaction, its product states and their sources
  ! are in the CXHEAT rows written by tools/fitting/make_element_data.py.
  ! E_rad counts a metastable product (S+ 2D, 2P) as radiated, which holds
  ! below its critical density (~1e3 cm^-3 for S+ 2D); above it the
  ! excitation is partly returned to the gas by collisional de-excitation.
  ! Where the ionization rate includes radiative charge transfer (carbon,
  ! cxheat(3) > 0) only its non-radiative part heats, and that part is the
  ! detailed-balance reverse of the recombination rate,
  !   k_nr = min[k_CXI, g k_CXR exp(dE_ion/kT)],
  ! so the pair deposits nothing when the two reactions balance.  Stage
  ! fractions and rate coefficients are those of the cascade.
  !=========================================================================
  real(kind=wp) function metal_charge_exchange_heating(il, T, nH, ne, nHI, nHII) &
                                                       result(heat)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne, nHI, nHII
    real(kind=wp) :: frac(MAX_ST), nel, k_ion, k_rec
    real(kind=wp) :: kf(MAX_CXS_TERM), kr(MAX_CXS_TERM), fpop(MAX_CXS_TERM)
    real(kind=wp) :: de(MAX_CXS_TERM)
    integer :: ie, it, nt
    heat = 0.0_wp
    do ie = 1, n_elements
       if (all(elems(ie)%cxi_form(1:elems(ie)%nstage-1) <= 0) .and. &
           all(elems(ie)%cxr_form(1:elems(ie)%nstage-1) <= 0)) cycle
       call species_fractions(ie, il, T, ne, nHI, nHII, frac)
       nel = elems(ie)%abund*nH
       do it = 1, elems(ie)%nstage-1
          if (elems(ie)%cxi_form(it) == 9) then
             !--- by product term: into term j the reaction releases
             !--- dIP - E_j (the term's excitation radiated, as for the
             !--- CXHEAT products), and the reverse from j takes it back
             call cxs_forward_terms(ie, it, T, kf)
             call cxs_reverse_terms(ie, it, T, kr)
             call cxs_populations(ie, it, T, ne, fpop)
             nt = elems(ie)%ncxs_term(it)
             de(1:nt) = elems(ie)%cxs_dip(it) - elems(ie)%cxs_E(1:nt,it)
             heat = heat + nel*( frac(it)*nHII*sum(kf(1:nt)*de(1:nt)) &
                               - frac(it+1)*nHI*sum(fpop(1:nt)*kr(1:nt)*de(1:nt)) )
             cycle
          end if
          k_ion = 0.0_wp;  k_rec = 0.0_wp
          if (elems(ie)%cxi_form(it) > 0) k_ion = cxi_rate(ie, it, T)
          if (elems(ie)%cxr_form(it) > 0) k_rec = cxr_rate(ie, it, T)
          if (elems(ie)%cxheat(3,it) > 0.0_wp) &
             k_ion = min(k_ion, elems(ie)%cxheat(3,it)*k_rec &
                     *exp(min(elems(ie)%cxheat(1,it)*ev2erg/(kboltz_cgs*T), 600.0_wp)))
          heat = heat + nel*( frac(it)*nHII*k_ion*cx_ionization_defect(ie, it, T) &
                            + frac(it+1)*nHI*k_rec*elems(ie)%cxheat(2,it) )
       end do
    end do
    heat = heat*ev2erg
  end function metal_charge_exchange_heating

  !--- kinetic energy [eV] of one X^it + H+ -> X^(it+1) + H event at T: the
  !--- CXHEAT value, or the CXHEATTAB interpolation (linear in T, end values
  !--- held) where the product channels branch with temperature
  real(kind=wp) function cx_ionization_defect(ie, it, T) result(de)
    integer,       intent(in) :: ie, it
    real(kind=wp), intent(in) :: T
    integer :: n, k
    real(kind=wp) :: w
    n = elems(ie)%ncxheat_tab(it)
    if (n == 0) then
       de = elems(ie)%cxheat(1,it)
       return
    end if
    if (T <= elems(ie)%cxheat_tab_T(1,it)) then
       de = elems(ie)%cxheat_tab_e(1,it);  return
    end if
    if (T >= elems(ie)%cxheat_tab_T(n,it)) then
       de = elems(ie)%cxheat_tab_e(n,it);  return
    end if
    k = 1
    do while (elems(ie)%cxheat_tab_T(k+1,it) < T)
       k = k + 1
    end do
    w  = (T - elems(ie)%cxheat_tab_T(k,it)) &
         /(elems(ie)%cxheat_tab_T(k+1,it) - elems(ie)%cxheat_tab_T(k,it))
    de = (1.0_wp - w)*elems(ie)%cxheat_tab_e(k,it) + w*elems(ie)%cxheat_tab_e(k+1,it)
  end function cx_ionization_defect

  !=========================================================================
  ! The loss side of the metal ionization energy budget [erg s^-1 cm^-3]
  ! (par%metal_heat, the counterpart of metal_heating):
  !  - recombination removes the kinetic energy of the captured electron,
  !    beta = k T [ (3/2) alpha + d alpha/d ln T ], exact for a Maxwellian and
  !    any capture cross section (radiative and dielectronic), from the same
  !    alpha_rec the cascade uses (see cooling_mod's recombination cooling);
  !  - collisional ionization removes the ionization energy E_th of the stage.
  ! Each transition it (stage it -> it+1) contributes
  !   n_e [ n_(it+1) beta_it + n_it C_it E_th,it ].
  ! The kinetic energy the charge-exchange reactions exchange with the gas is
  ! metal_charge_exchange_heating (par%charge_exchange_heating).
  !=========================================================================
  real(kind=wp) function metal_ionization_cooling(il, T, nH, ne, nHI, nHII) &
                                                  result(cool)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne, nHI, nHII
    real(kind=wp), parameter :: H = 1.0e-4_wp       ! ln T step of d alpha/d ln T
    real(kind=wp) :: frac(MAX_ST), nel, a0, beta, cci
    integer :: ie, it
    cool = 0.0_wp
    do ie = 1, n_elements
       call species_fractions(ie, il, T, ne, nHI, nHII, frac)
       nel = elems(ie)%abund*nH
       do it = 1, elems(ie)%nstage-1
          a0   = alpha_rec(elems(ie), it, T)
          beta = kboltz_cgs*T*(1.5_wp*a0 + (alpha_rec(elems(ie), it, T*exp(H)) &
                 - alpha_rec(elems(ie), it, T*exp(-H)))/(2.0_wp*H))
          cci = metal_electron_impact_ionization(elems(ie), it, T)
          cool = cool + ne*nel*( frac(it+1)*beta &
                 + frac(it)*cci*elems(ie)%eth(it)*ev2erg )
       end do
    end do
  end function metal_ionization_cooling

  real(kind=wp) function metal_electron_impact_ionization(element, it, T) result(rate)
    type(element_type), intent(in) :: element
    integer, intent(in) :: it
    real(kind=wp), intent(in) :: T
    ! Use the same coefficient in direct balance, cached balance, and the
    ! ionization-energy cooling term. Disabled means no rate calculation.
    rate = 0.0_wp
    if (.not. par%metal_collisional_ionization) return
    rate = voronov_ci(element%ci(1:5,it), T)*ci_dere_ratio(element%Z, it)
  end function metal_electron_impact_ionization

  real(kind=wp) function el8(ie, it, k)
    integer, intent(in) :: ie, it, k
    el8 = elems(ie)%photo(k, it)
  end function el8

  !=========================================================================
  ! Photoionization cross section sig [cm^2] of (element ie, stage it) at
  ! E [eV], and the photoheating coefficient heat [cm^2] = sum over the
  ! absorbing shells of sigma_s (1 - E_th,s/E), phfit2's prescription:
  !   E <  einn : the VFKY96 outer-shell fit, heat = sig (1 - E_th/E);
  !   E >= einn : sum of the VY95 fits of every occupied subshell with
  !               E >= E_th,s (einn = 0: at every E);
  ! except where the element file carries a PHTABLE (Ne III): below its
  ! last node E_join the tabulated R-matrix ground-term cross section,
  ! linear in E between its nodes, at the exact E; heat = sig (1 - E_th/E),
  ! the convention of the L-shell fit it replaces (the table is a total over
  ! the Ne IV exit terms; an exit into an excited term counts its excitation,
  ! >= 5.1 eV above E_th, as heat -- written at the generator,
  ! tools/fitting/ground_term_photoionization_table.py).
  ! The photoelectron of an inner-shell ionization carries E - E_th,s; the
  ! excitation E_th,s - E_th left in the ion is assumed to leave the gas as
  ! fluorescence (not transported) and the ion to move up ONE stage.
  ! KNOWN DEFECT (flagged 2026-09-24): that holds only where the hole lies
  ! below the double-ionization threshold IP_i + IP_(i+1).  Inside the
  ! 13.6-100 eV band it does for Mg II 2p (65.7 < 95.2 eV), Ca II 3p
  ! (40.9 < 62.8) and Cl I-IV 3s, but NOT for Mg I 2p (54.9 > 22.7 eV),
  ! Ca I 3p (34.4 > 18.0), Fe I 3p (66.0 > 24.1), Fe II 3p (76.2 > 46.9) and
  ! Fe III 3p (87.1 > 85.6): those holes autoionize (Auger yields 0.995-0.9999,
  ! Kaastra & Mewe 1993 as tabulated in Cloudy's mewe_nelectron.dat), so the
  ! ion should move up TWO stages and the Auger electron, ~E_th,s - IP_i -
  ! IP_(i+1), should heat the gas.  Carrying that needs the stage matrix
  ! (md/METAL_COUPLING_PLAN.md M4): the nearest-neighbor ratio chain in
  ! species_fractions cannot represent a two-stage jump.
  !=========================================================================
  subroutine species_cross_section(ie, it, E, sig, heat)
    use photo_xsec, only : sigma_vy95
    integer,       intent(in)  :: ie, it
    real(kind=wp), intent(in)  :: E
    real(kind=wp), intent(out) :: sig, heat
    real(kind=wp) :: s
    integer :: k
    sig = 0.0_wp;  heat = 0.0_wp
    if (E < elems(ie)%eth(it)) return
    if (photo_table_covers(ie, it, E)) then
       sig  = photo_table_sigma(ie, it, E)
       heat = sig*max(1.0_wp - elems(ie)%eth(it)/E, 0.0_wp)
    else if (E < elems(ie)%einn(it)) then
       sig = sigma_vfky96(E, el8(ie,it,1), el8(ie,it,2), el8(ie,it,3), &
                          el8(ie,it,4), el8(ie,it,5), el8(ie,it,6), &
                          el8(ie,it,7), el8(ie,it,8))
       heat = sig*max(1.0_wp - elems(ie)%eth(it)/E, 0.0_wp)
    else
       do k = 1, elems(ie)%nshell(it)
          associate (sh => elems(ie)%shell(:,k,it))
            if (E < sh(2)) cycle
            s = sigma_vy95(E, sh(2), sh(3), sh(4), sh(5), sh(6), sh(7), sh(1))
            sig  = sig + s
            heat = heat + s*max(1.0_wp - sh(2)/E, 0.0_wp)
          end associate
       end do
    end if
  end subroutine species_cross_section

  !=========================================================================
  ! Ground-state photoionization cross section [cm^2] of (element ie, stage
  ! it) at E [eV] that the Milne relation turns into the ground-state
  ! recombination X^(it+1) + e -> X^it(ground term) + h nu (the metal ground
  ! continua of the diffuse field, milne_recomb_spectrum_mod).
  !
  ! WHICH PART OF THE phfit2 SUM.  Only the OUTER shell: removing an inner
  ! electron leaves the ion with a hole, whose inverse is a capture into an
  ! excited configuration of X^it, not into its ground term.  The outer
  ! shell is what species_cross_section uses for it: below the first
  ! inner-shell edge einn the VFKY96 outer-shell fit (phfit2's PH2), at and
  ! above einn the Verner & Yakovlev (1995) fit of the subshell whose
  ! threshold is the ionization potential (the PHSHELL row the element file
  ! writes with ETH; einn = 0, Cl: that subshell at every E).  So the
  ! emitted photons are absorbed with the outer-shell part of the cross
  ! section the transport uses.
  ! VALIDITY.  The outer-shell fit is the photoionization of the ground term
  ! of X^it summed over the final terms of X^(it+1).  Its inverse is the
  ! capture from the ground term of X^(it+1) only below the ion's first
  ! excited term E_1 (data/atomic/ground_term_weights.txt): above
  ! E_th + E_1 the fit also carries photoionization into excited terms of
  ! the ion, whose captures (from ions that are not in those terms) it
  ! counts as ground-term captures, an overestimate bounded by the share of
  ! the Milne integral above E_th + E_1 (printed at setup; e.g. O+ E_1 =
  ! 3.32 eV, 2% at 1e4 K).  Where the ion has terms within kT of its ground
  ! (Fe+ a4F at 0.23 eV) the bound is O(1) and the channel is an upper limit.
  ! Resonances are smoothed by the fits (VFKY96, VY95).  Where the element
  ! file carries a PHTABLE (Ne III), the tabulated R-matrix ground-term cross
  ! section below its last node, the same function species_cross_section
  ! transports with; it resolves the resonances, so a quadrature over it
  ! must step on its nodes (species_photo_table_nodes).
  !=========================================================================
  real(kind=wp) function species_outer_shell_cross_section(ie, it, E) result(sig)
    use photo_xsec, only : sigma_vy95
    integer,       intent(in) :: ie, it
    real(kind=wp), intent(in) :: E
    integer :: k
    sig = 0.0_wp
    if (E < elems(ie)%eth(it)) return
    if (photo_table_covers(ie, it, E)) then
       sig = photo_table_sigma(ie, it, E)
       return
    end if
    if (E < elems(ie)%einn(it)) then
       sig = sigma_vfky96(E, el8(ie,it,1), el8(ie,it,2), el8(ie,it,3), &
                          el8(ie,it,4), el8(ie,it,5), el8(ie,it,6), &
                          el8(ie,it,7), el8(ie,it,8))
       return
    end if
    k = outer_shell_index(ie, it)
    associate (sh => elems(ie)%shell(:,k,it))
      sig = sigma_vy95(E, sh(2), sh(3), sh(4), sh(5), sh(6), sh(7), sh(1))
    end associate
  end function species_outer_shell_cross_section

  !--- energy [eV] where species_outer_shell_cross_section changes from the
  !--- VFKY96 fit to the VY95 outer subshell (einn), 0 when it does not; a
  !--- quadrature across it must break there (the two fits do not join)
  real(kind=wp) function species_outer_shell_break(ie, it) result(e)
    integer, intent(in) :: ie, it
    e = 0.0_wp
    if (elems(ie)%einn(it) > elems(ie)%eth(it) .and. &
        elems(ie)%einn(it) < 1.0e29_wp) e = elems(ie)%einn(it)
  end function species_outer_shell_break

  !--- E [eV] lies on the PHTABLE of (ie, it): E_th <= E < E_join
  logical function photo_table_covers(ie, it, E) result(yes)
    integer,       intent(in) :: ie, it
    real(kind=wp), intent(in) :: E
    integer :: n
    yes = .false.
    n = elems(ie)%nphtab(it)
    if (n == 0) return
    yes = E - elems(ie)%eth(it) < elems(ie)%phtab(it)%eps(n)
  end function photo_table_covers

  !--- PHTABLE cross section [cm^2] at E (E_th <= E < E_join), linear in
  !--- the photoelectron energy between the two nodes that bracket it
  real(kind=wp) function photo_table_sigma(ie, it, E) result(sig)
    integer,       intent(in) :: ie, it
    real(kind=wp), intent(in) :: E
    real(kind=wp) :: x, w
    integer :: lo, hi, mid
    associate (tb => elems(ie)%phtab(it))
      x  = max(E - elems(ie)%eth(it), 0.0_wp)
      lo = 1;  hi = elems(ie)%nphtab(it)
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (tb%eps(mid) <= x) then
            lo = mid
         else
            hi = mid
         end if
      end do
      w   = (x - tb%eps(lo))/(tb%eps(hi) - tb%eps(lo))
      sig = tb%sig(lo) + w*(tb%sig(hi) - tb%sig(lo))
    end associate
  end function photo_table_sigma

  !--- E_join [eV], the photon energy of the last PHTABLE node of (ie, it);
  !--- 0 when the transition has no table
  real(kind=wp) function species_photo_table_end(ie, it) result(e)
    integer, intent(in) :: ie, it
    e = 0.0_wp
    if (elems(ie)%nphtab(it) > 0) &
       e = elems(ie)%eth(it) + elems(ie)%phtab(it)%eps(elems(ie)%nphtab(it))
  end function species_photo_table_end

  !--- indices j1..j2 of the PHTABLE nodes of (ie, it) with photon energies
  !--- strictly inside (ea, eb) [eV]; j1 > j2 when there is none
  subroutine species_photo_table_nodes(ie, it, ea, eb, j1, j2)
    integer,       intent(in)  :: ie, it
    real(kind=wp), intent(in)  :: ea, eb
    integer,       intent(out) :: j1, j2
    integer :: n
    j1 = 1;  j2 = 0
    n = elems(ie)%nphtab(it)
    if (n == 0) return
    j1 = first_node_above(ea - elems(ie)%eth(it))
    j2 = first_node_above(eb - elems(ie)%eth(it)) - 1
    !--- the same comparison in photon energy, so a node on an end is out
    if (j1 <= n) then
       if (elems(ie)%eth(it) + elems(ie)%phtab(it)%eps(j1) <= ea) j1 = j1 + 1
    end if
    if (j2 >= 1) then
       if (elems(ie)%eth(it) + elems(ie)%phtab(it)%eps(j2) >= eb) j2 = j2 - 1
    end if
  contains
    !--- the first node with eps > x (n + 1 when none)
    integer function first_node_above(x) result(j)
      real(kind=wp), intent(in) :: x
      integer :: lo, hi, mid
      associate (ep => elems(ie)%phtab(it)%eps)
        if (ep(1) > x) then
           j = 1;  return
        end if
        if (ep(n) <= x) then
           j = n + 1;  return
        end if
        lo = 1;  hi = n
        do while (hi - lo > 1)
           mid = (lo + hi)/2
           if (ep(mid) <= x) then
              lo = mid
           else
              hi = mid
           end if
        end do
        j = hi
      end associate
    end function first_node_above
  end subroutine species_photo_table_nodes

  !--- photon energy [eV] of PHTABLE node j of (ie, it)
  real(kind=wp) function species_photo_table_node_energy(ie, it, j) result(e)
    integer, intent(in) :: ie, it, j
    e = elems(ie)%eth(it) + elems(ie)%phtab(it)%eps(j)
  end function species_photo_table_node_energy

  !--- the PHSHELL row of (ie, it) whose threshold is the ionization
  !--- potential ETH (make_element_data.py writes the outer edge with ETH)
  integer function outer_shell_index(ie, it) result(kout)
    integer, intent(in) :: ie, it
    integer :: k
    kout = 0
    do k = 1, elems(ie)%nshell(it)
       if (abs(elems(ie)%shell(2,k,it) - elems(ie)%eth(it)) <= &
           1.0e-9_wp*elems(ie)%eth(it)) kout = k
    end do
    if (kout == 0) call fatal_error('no outer subshell at the ionization '// &
       'potential for '//trim(elems(ie)%name)//' (element file PHSHELL rows)')
  end function outer_shell_index

  !--- the metal transitions in threshold order (the packet-cache order)
  integer function species_transition_count() result(n)
    n = flat_ntransition
  end function species_transition_count

  subroutine species_transition(k, ie, it)
    integer, intent(in)  :: k
    integer, intent(out) :: ie, it
    ie = flat_ie(k);  it = flat_it(k)
  end subroutine species_transition

  !--- stage fraction of stage istage of element ie in leaf il, as cached for
  !--- the packet opacity from the current gas state
  !--- (species_stage_cache_refresh)
  real(kind=wp) function species_stage_fraction(ie, istage, il) result(f)
    integer, intent(in) :: ie, istage, il
    f = stage_frac(ie)%g(istage, il)
  end function species_stage_fraction

  !--- total (RR + DR) recombination coefficient of transition it [cm^3/s]
  real(kind=wp) function species_recombination_coefficient(ie, it, T) result(a)
    integer,       intent(in) :: ie, it
    real(kind=wp), intent(in) :: T
    a = alpha_rec(elems(ie), it, T)
  end function species_recombination_coefficient

  !--- its radiative part alone [cm^3/s] (the RR or RR2 row of alpha_rec)
  real(kind=wp) function species_radiative_recombination(ie, it, T) result(a)
    integer,       intent(in) :: ie, it
    real(kind=wp), intent(in) :: T
    type(element_type) :: el
    el = elems(ie)
    el%ndr(it) = 0
    el%has_dr2(it) = .false.
    a = alpha_rec(el, it, T)
  end function species_radiative_recombination

  !=========================================================================
  ! Stage fractions of element ie in leaf il at (T, n_e, x_HI): stationary
  ! adjacent-stage balance, normalized without cumulative rate products.
  !=========================================================================
  subroutine species_fractions(ie, il, T, ne, nHI, nHII, frac)
    implicit none
    integer,       intent(in)  :: ie, il
    real(kind=wp), intent(in)  :: T, ne, nHI, nHII
    real(kind=wp), intent(out) :: frac(MAX_ST)
    real(kind=wp) :: ionization(MAX_ST-1), recombination(MAX_ST-1)
    integer :: it, ns, status

    frac = 0.0_wp
    ns = elems(ie)%nstage
    if (ne == 0.0_wp .and. nHII == 0.0_wp .and. &
        all(egam(ie)%g(:,il) == 0.0_wp)) then
       ! In unilluminated gas with no electrons or protons, charge
       ! neutrality requires neutral metals: negative ions are not included
       ! in this model. No upward reaction can populate a higher stage.
       ! This constraint fixes the neutral initial state even when the
       ! stage equations alone have disconnected zero-rate equilibria.
       frac(1) = 1.0_wp
       return
    end if
    do it = 1, ns-1
       ionization(it) = egam(ie)%g(it,il) &
            + ne*metal_electron_impact_ionization(elems(ie), it, T)
       if (elems(ie)%cxi_form(it) > 0) &
          ionization(it) = ionization(it) + nHII*cxi_rate(ie, it, T)
       recombination(it) = ne*alpha_rec(elems(ie), it, T)
       if (elems(ie)%cxr_form(it) > 0) &
          recombination(it) = recombination(it) + nHI*cxr_rate(ie, it, T, ne)
    end do
    call adjacent_ion_fractions(ionization(1:ns-1), recombination(1:ns-1), frac(1:ns), status)
    if (status /= 0) call ion_stage_balance_error(ie, status)
  end subroutine species_fractions

  subroutine ion_stage_balance_error(ie, status)
    integer, intent(in) :: ie, status
    if (status == 2) then
       call fatal_error('nonunique stationary ion populations for '//trim(elems(ie)%name)// &
                        ': disconnected zero-rate stages require initial populations')
    else
       call fatal_error('invalid ionization or recombination rates for '//trim(elems(ie)%name))
    end if
  end subroutine ion_stage_balance_error

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
    real(kind=wp) :: frac(MAX_ST), nel, line_ne, coefficient
    integer :: ie, i

    cool = 0.0_wp
    ! No electron-impact line radiation is produced in vacuum or without
    ! electrons. Do not form an indeterminate stationary ion balance there.
    if (nH == 0.0_wp .or. ne == 0.0_wp) return
    do ie = 1, n_elements
       call species_fractions(ie, il, T, ne, nHI, nHII, frac)
       nel = elems(ie)%abund*nH
       do i = 1, elems(ie)%nstage
          coefficient = 0.0_wp
          if (elems(ie)%has_cool(i)) then
             coefficient = nlevel_cooling_apply(elems(ie)%cool_idx(i), &
                    cooling_eval(elems(ie)%cool(i), T), T, ne)
          else if (elems(ie)%cool_idx(i) > 0) then
             line_ne = ne
             ! The table clamps below its 1e-6 cm^-3 density column. For the
             ! requested low-density model use that numerical reference,
             ! not the cell density. Source accuracy is still conditional.
             if (trim(par%cooling_model) == 'low_density') line_ne = tiny(1.0_wp)
             coefficient = nlevel_cooling_lam_se(elems(ie)%cool_idx(i), T, line_ne)
          end if
          cool = cool + ne*nel*frac(i)*coefficient
       end do
    end do
  end function metal_cooling

  !=========================================================================
  ! Metal free-free (thermal bremsstrahlung) cooling per unit volume
  ! [erg cm^-3 s^-1] at trial T (par%metal_freefree): the trace completion
  ! of the H II + He II/III free-free term.  A metal ion in stage j
  ! (j = 1 neutral, j = 2 = X+, ...) has net charge Z_ion = j - 1 and
  ! contributes n_ion Z_ion^2 gbar_ff(T, Z_ion); only j >= 2 radiate.
  ! Same prefactor FF_TOTAL sqrt(T) n_e as the H/He free-free (free_free_mod).
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
    cool_ff = FF_TOTAL*sqrt(T)*ne*nH*s
  end function metal_freefree

  !=========================================================================
  ! Converged stage fractions of every element -> output file blocks.
  !=========================================================================
  subroutine species_write(file, status)
    use iofile_mod
    use octree_mod,    only : amr_grid
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHII, gas_ne, gas_Te
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
          nHII = gas_nH(il)*gas_xHII(il)
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
