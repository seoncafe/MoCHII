module species_mod
!---------------------------------------------------------------------------
! MoCHII: trace-metal species registry + element ionization cascade (G2b).
!
! Design rule (docs/PLAN.md section 8): adding an ion is a data operation.
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
! Tier-1 line cooling: each (element, stage) with a coefficient file
! data/atomic/cooling_tier1_<el>_<i>.txt contributes
! n_e n_stage Lambda(T) to the thermal balance (metal_cooling).
!---------------------------------------------------------------------------
  use define
  use photo_xsec, only : sigma_vfky96
  use cooling_mod, only : tier1_fit_type, tier1_load, tier1_eval
  implicit none
  private

  public :: species_setup, species_gamma_compute, species_fractions
  public :: metal_cooling, species_write, species_resize
  public :: species_opacity_add, species_ne, metal_heating
  public :: metal_cooling_H
  public :: n_elements, elem_name, elem_nstage, elem_abund

  integer, parameter :: MAX_EL = 8, MAX_ST = 6

  type element_type
     character(len=2) :: name = ''
     integer  :: Z = 0
     integer  :: nstage = 0
     real(kind=wp) :: abund = 0.0_wp
     !--- transition data (index 1..nstage-1)
     real(kind=wp) :: eth(MAX_ST)      = 0.0_wp
     real(kind=wp) :: photo(8,MAX_ST)  = 0.0_wp   ! Eth,E0,s0,ya,P,yw,y0,y1
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
     !--- Tier-1 cooling fits per stage (loaded when the file exists)
     logical :: has_cool(MAX_ST) = .false.
     type(tier1_fit_type) :: cool(MAX_ST)
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

contains

  character(len=2) function elem_name(ie);  integer, intent(in) :: ie
    elem_name = elems(ie)%name;  end function elem_name
  integer function elem_nstage(ie);  integer, intent(in) :: ie
    elem_nstage = elems(ie)%nstage;  end function elem_nstage
  real(kind=wp) function elem_abund(ie);  integer, intent(in) :: ie
    elem_abund = elems(ie)%abund;  end function elem_abund

  !=========================================================================
  subroutine species_setup(nleaf)
    use mpi
    implicit none
    integer, intent(in) :: nleaf
    character(len=8)   :: names(8)
    real(kind=wp)      :: abunds(8)
    character(len=256) :: fname
    logical :: exists
    integer :: k, ie, i, ierr

    names  = [character(len=8) :: 'c', 'n', 'o', 'ne', 's', 'ar', 'mg', 'fe']
    abunds = [par%abund_C, par%abund_N, par%abund_O, par%abund_Ne, &
              par%abund_S, par%abund_Ar, par%abund_Mg, par%abund_Fe]

    n_elements = 0
    do k = 1, size(names)
       if (abunds(k) <= 0.0_wp) cycle
       n_elements = n_elements + 1
       ie = n_elements
       call read_element(trim(par%atomic_dir)//'/element_'// &
                         trim(names(k))//'.txt', elems(ie))
       elems(ie)%abund = abunds(k)
       !--- Tier-1 cooling fits for each stage where a file exists.
       do i = 1, elems(ie)%nstage
          write(fname,'(a,a,a,i0,a)') trim(par%atomic_dir)// &
             '/cooling_tier1_', trim(names(k)), '_', i, '.txt'
          inquire(file=trim(fname), exist=exists)
          if (exists) then
             call tier1_load(trim(fname), elems(ie)%cool(i))
             elems(ie)%has_cool(i) = .true.
          end if
       end do
       allocate(egam(ie)%g(elems(ie)%nstage-1, nleaf))
       egam(ie)%g = 0.0_wp
       allocate(eheat(ie)%g(elems(ie)%nstage-1, nleaf))
       eheat(ie)%g = 0.0_wp
       if (mpar%p_rank == 0) write(*,'(3a,es10.3,a,i2,a,i2,a)') &
          ' SPEC: ', trim(names(k)), ' loaded (abund=', elems(ie)%abund, &
          ', ', elems(ie)%nstage, ' stages, ', &
          count(elems(ie)%has_cool(1:elems(ie)%nstage)), ' cooling fits)'
    end do
  end subroutine species_setup

  !=========================================================================
  !--- G4: resize the Gamma blocks after octree re-refinement.
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
    end do
  end subroutine species_resize

  !=========================================================================
  subroutine read_element(fname, el)
    use mpi
    implicit none
    character(len=*),   intent(in)  :: fname
    type(element_type), intent(out) :: el
    character(len=512) :: line, key
    integer :: unit, ios, it, ierr, n

    open(newunit=unit, file=fname, status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(2a)') 'ERROR: cannot open ', trim(fname)
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
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
    use gas_state_mod, only : gas_nH, gas_xHI, gas_ne, gas_Te
    implicit none
    integer,       intent(in)    :: nnu, nleaf
    real(kind=wp), intent(inout) :: kap(nnu, nleaf)
    real(kind=wp) :: sig(nnu, MAX_ST), frac(MAX_ST), nHI, nHII, add
    integer :: ie, it, il, inu

    do ie = 1, n_elements
       do it = 1, elems(ie)%nstage-1
          do inu = 1, nnu
             sig(inu,it) = sigma_vfky96(ion_e(inu), el8(ie,it,1), &
                el8(ie,it,2), el8(ie,it,3), el8(ie,it,4), el8(ie,it,5), &
                el8(ie,it,6), el8(ie,it,7), el8(ie,it,8))
          end do
       end do
       do il = 1, nleaf
          if (gas_nH(il) <= 0.0_wp) cycle
          nHI  = gas_nH(il)*gas_xHI(il)
          nHII = gas_nH(il)*(1.0_wp - gas_xHI(il))
          call species_fractions(ie, il, gas_Te(il), gas_ne(il), &
                                 nHI, nHII, frac)
          do inu = 1, nnu
             add = 0.0_wp
             do it = 1, elems(ie)%nstage-1
                add = add + frac(it)*sig(inu,it)
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
    use octree_mod,   only : amr_grid
    use jtally_mod,   only : jt_ion
    use ion_band_mod, only : ion_e
    implicit none
    integer :: ie, it, il, inu, ic, nleaf
    real(kind=wp) :: sig(par%nnu_ion), vol, fac, fJ, fH

    nleaf = amr_grid%nleaf
    do ie = 1, n_elements
       do it = 1, elems(ie)%nstage-1
          do inu = 1, par%nnu_ion
             sig(inu) = sigma_vfky96(ion_e(inu), el8(ie,it,1), el8(ie,it,2), &
                        el8(ie,it,3), el8(ie,it,4), el8(ie,it,5), &
                        el8(ie,it,6), el8(ie,it,7), el8(ie,it,8))
          end do
          do il = 1, nleaf
             ic  = amr_grid%icell_of_leaf(il)
             vol = (2.0_wp*amr_grid%ch(ic))**3
             fac = 1.0_wp/(vol*par%distance2cm**2)
             fJ  = 0.0_wp
             fH  = 0.0_wp
             do inu = 1, par%nnu_ion
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
  ! PDR-zone coolant (n_HI/n_e ~ 10^3-10^4 there; the Tier-1 fits are
  ! electron-impact only).  Two-level, low-density limit:
  !   Lambda = n_ion n_HI (g_u/g_l) q_ul^H e^{-dE/kT} dE
  ! with q_ul^H([C II]) = 7.6e-10 (T/100)^0.14 (Barinovs et al. 2005)
  ! and  q_ul^H([O I])  = 9.2e-11 (T/100)^0.67 cm^3/s.  Part of the
  ! par%grain_pe PDR thermal package (without it the photoelectric
  ! heating has no coolant below the Ly-alpha regime and the PDR zone
  ! runs away to ~10^4 K).
  !=========================================================================
  real(kind=wp) function metal_cooling_H(il, T, nH, ne, nHI, nHII) result(cool)
    implicit none
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne, nHI, nHII
    real(kind=wp), parameter :: kb = 1.380649e-16_wp
    real(kind=wp) :: frac(MAX_ST), nion, qlu
    integer :: ie
    cool = 0.0_wp
    do ie = 1, n_elements
       select case (trim(elems(ie)%name))
       case ('c')
          if (elems(ie)%nstage < 2) cycle
          call species_fractions(ie, il, T, ne, nHI, nHII, frac)
          nion = elems(ie)%abund*nH*frac(2)                 ! C II
          qlu  = 2.0_wp*7.6e-10_wp*(T/100.0_wp)**0.14_wp &  ! g_u/g_l = 4/2
                 *exp(-91.25_wp/T)
          cool = cool + nion*nHI*qlu*(91.25_wp*kb)
       case ('o')
          call species_fractions(ie, il, T, ne, nHI, nHII, frac)
          nion = elems(ie)%abund*nH*frac(1)                 ! O I
          qlu  = 0.6_wp*9.2e-11_wp*(T/100.0_wp)**0.67_wp &  ! g_u/g_l = 3/5
                 *exp(-227.7_wp/T)
          cool = cool + nion*nHI*qlu*(227.7_wp*kb)
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
            + ne*voronov_ci(elems(ie)%ci(1:5,it), T)
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
             cool = cool + ne*nel*frac(i)*tier1_eval(elems(ie)%cool(i), T)
       end do
    end do
  end function metal_cooling

  !=========================================================================
  ! Converged stage fractions of every element -> output file blocks.
  !=========================================================================
  subroutine species_write(file)
    use iofile_mod
    use octree_mod,    only : amr_grid
    use gas_state_mod, only : gas_nH, gas_xHI, gas_ne, gas_Te
    implicit none
    type(io_file_type), intent(inout) :: file
    real(kind=wp), allocatable :: fr(:,:)
    real(kind=wp) :: frac(MAX_ST), nHI, nHII
    character(len=64) :: extname
    integer :: ie, il, i, status, nleaf

    status = 0
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
    end do
  end subroutine species_write

end module species_mod
