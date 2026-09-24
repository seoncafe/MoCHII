program check_photoelectron_partition
!---------------------------------------------------------------------------
! Gate for photoelectron_degradation_mod (par%use_sec_ion):
!  (a) at the top node (E0 >= 1 keV) and SvS85's own composition
!      (n_HeI/n_HI = 0.1) the secondary-ionization budget is exactly Shull &
!      van Steenberg (1985) Table 2: n_HI c_HI = phi_HI(x), n_HeI c_HeI =
!      phi_HeI(x);
!  (b) energy is conserved at every node and composition,
!      n_HI c_HI + n_HeI c_HeI = f_ion, and f_heat + f_ion <= 1;
!  (c) the heat share is Dalgarno, Yan & Liu (1999) eq. (14) with the (1-x)
!      closure: eta0 of their Table 7 at x -> 0 and exactly 1 at x = 1;
!  (d) the energy weights sum to one above 30 eV, reproduce each node and
!      are linear in log E0 between nodes;
!  (e) record: heat share and ionizations per photoelectron at E0 = 100 eV,
!      old SvS85-only form against this one.
!---------------------------------------------------------------------------
  use define, only : wp
  use photoelectron_degradation_mod
  implicit none
  real(kind=wp) :: fh(N_EDEP), cH(N_EDEP), cHe(N_EDEP), w(N_EDEP), x, worst, r
  real(kind=wp), parameter :: XS(6) = [1.0e-4_wp, 1.0e-3_wp, 1.0e-2_wp, &
                                       0.1_wp, 0.5_wp, 0.9_wp]
  real(kind=wp), parameter :: RATIO(4) = [0.0_wp, 0.1_wp, 1.0_wp, 10.0_wp]
  real(kind=wp), parameter :: ETA0(N_EDEP) = [0.259_wp, 0.186_wp, 0.148_wp, &
                                              0.132_wp, 0.122_wp, 0.117_wp]
  real(kind=wp), parameter :: ENODE(N_EDEP) = [30.0_wp, 50.0_wp, 100.0_wp, &
                                               200.0_wp, 500.0_wp, 1000.0_wp]
  real(kind=wp) :: fion(N_EDEP), nHI, nHeI
  integer :: i, j, k
  logical :: ok

  ok = .true.
  !--- (a)
  worst = 0.0_wp
  do i = 1, size(XS)
     x = XS(i)
     call photoelectron_partition(x, 1.0_wp, 0.1_wp, fh, cH, cHe)
     worst = max(worst, abs(cH(N_EDEP)*1.0_wp/svs85_fion_HI(x) - 1.0_wp), &
                        abs(cHe(N_EDEP)*0.1_wp/svs85_fion_HeI(x) - 1.0_wp))
  end do
  write(*,'(a,es9.2)') ' (a) 1 keV, SvS85 composition = Table 2: worst ', worst
  if (worst > 1.0e-13_wp) ok = .false.

  !--- (b)
  worst = 0.0_wp
  do i = 1, size(XS)
     do j = 1, size(RATIO)
        nHI = 1.0_wp;  nHeI = RATIO(j)
        call photoelectron_partition(XS(i), nHI, nHeI, fh, cH, cHe)
        fion = nHI*cH + nHeI*cHe
        !--- f_ion from the budget: recompute at the SvS85 composition
        block
          real(kind=wp) :: fh0(N_EDEP), cH0(N_EDEP), cHe0(N_EDEP)
          call photoelectron_partition(XS(i), 1.0_wp, 0.1_wp, fh0, cH0, cHe0)
          worst = max(worst, maxval(abs(fion - (cH0 + 0.1_wp*cHe0)) &
                                    /max(cH0 + 0.1_wp*cHe0, 1.0e-300_wp)))
        end block
        if (any(fh + fion > 1.0_wp + 1.0e-14_wp)) ok = .false.
     end do
  end do
  write(*,'(a,es9.2)') ' (b) n_HI c_HI + n_HeI c_HeI = f_ion at any He/H: worst ', worst
  if (worst > 1.0e-13_wp) ok = .false.

  !--- (c)
  call photoelectron_partition(0.0_wp, 1.0_wp, 0.1_wp, fh, cH, cHe)
  worst = maxval(abs(fh/ETA0 - 1.0_wp))
  call photoelectron_partition(1.0_wp, 1.0_wp, 0.1_wp, fh, cH, cHe)
  worst = max(worst, maxval(abs(fh - 1.0_wp)))
  write(*,'(a,es9.2)') ' (c) f_heat = eta0 at x=0 (Table 7), 1 at x=1: worst ', worst
  if (worst > 1.0e-13_wp) ok = .false.

  !--- (d)
  worst = 0.0_wp
  do k = 1, N_EDEP
     call photoelectron_energy_weights(ENODE(k)*(1.0_wp + 1.0e-12_wp), w)
     worst = max(worst, abs(w(k) - 1.0_wp), abs(sum(w) - 1.0_wp))
  end do
  call photoelectron_energy_weights(sqrt(50.0_wp*100.0_wp), w)
  worst = max(worst, abs(w(2) - 0.5_wp), abs(w(3) - 0.5_wp))
  call photoelectron_energy_weights(5000.0_wp, w)
  worst = max(worst, abs(w(N_EDEP) - 1.0_wp))
  call photoelectron_energy_weights(29.9_wp, w)
  worst = max(worst, sum(w))
  write(*,'(a,es9.2)') ' (d) energy weights: worst ', worst
  if (worst > 1.0e-10_wp) ok = .false.

  !--- (e) record at E0 = 100 eV (node 3), SvS85 composition, H I absorber
  write(*,'(a)') ' (e) E0 = 100 eV, n_HeI/n_HI = 0.1: f_heat, N_HI, N_HeI per photoelectron'
  write(*,'(a)') '       x      SvS85-only (old)            this module'
  do i = 1, size(XS)
     x = XS(i)
     call photoelectron_partition(x, 1.0_wp, 0.1_wp, fh, cH, cHe)
     r = 0.9971_wp*(1.0_wp - (1.0_wp - x**0.2663_wp)**1.3163_wp)
     write(*,'(es9.1,3f8.3,4x,3f8.3)') x, r, svs85_fion_HI(x)*100.0_wp/13.598_wp, &
        svs85_fion_HeI(x)*100.0_wp/24.587_wp, fh(3), cH(3)*100.0_wp/13.598_wp, &
        0.1_wp*cHe(3)*100.0_wp/24.587_wp
  end do

  if (ok) then
     write(*,'(a)') 'PASS'
  else
     write(*,'(a)') 'FAIL'
     stop 1
  end if
end program check_photoelectron_partition
