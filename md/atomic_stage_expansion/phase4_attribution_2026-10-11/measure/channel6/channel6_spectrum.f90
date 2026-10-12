program channel6_spectrum
!---------------------------------------------------------------------------
! Diagnostic (Phase 4 attribution): the channel-6 luminosity of a written
! state, by metal transition, radial shell and energy window.  Leaf weight
! n_e V n(X^(i+1)) as diffuse_mod::metal_ground_weights; the band total uses
! milne_metal_band_coefficient at the leaf T (what diffuse_build sums), the
! windows milne_metal_energy_integral at the center of 0.01 dex T bins.
! Input: a stream file from dump_state.py.  Usage: channel6_spectrum.x <file>
!---------------------------------------------------------------------------
  use define
  use species_mod, only : species_setup, n_elements, elem_name, elem_abund, &
                          elem_eth, species_transition
  use milne_recomb_spectrum_mod, only : milne_setup, milne_metal_count, &
                          milne_metal_band_coefficient, milne_metal_energy_integral
  use mpi
  implicit none
  integer, parameter :: NSH = 4, NW = 5, NTB = 130
  real(kind=wp), parameter :: WE(NW+1) = [13.598_wp, 23.338_wp, 24.587_wp, 34.856_wp, 54.418_wp, 150.0_wp]
  real(kind=wp), parameter :: LT0 = 3.4_wp, DLT = 0.01_wp
  character(len=2), parameter :: ELN(5) = ['c ','n ','o ','ne','s ']
  integer :: ierr, nleaf, il, k, ie, it, j, ish, itb, iw, iel, kk
  real(kind=wp), allocatable :: T(:), ne(:), nH(:), V(:), sh(:), x(:,:,:)
  real(kind=wp), allocatable :: wk(:,:,:), lband(:,:), lwin(:,:,:)
  real(kind=wp) :: w, Tc
  character(len=512) :: fn
  character(len=8) :: ion

  call MPI_INIT(ierr);  mpar%p_rank = 0
  par%atomic_dir = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib/data_snapshot/atomic'
  call get_command_argument(2, fn)
  if (len_trim(fn) > 0) par%atomic_dir = fn
  par%use_metals = .true.
  par%cooling_model = 'low_density'
  par%eion_min = 13.598_wp;  par%eion_max = 150.0_wp
  par%abund_C = 2.2e-4_wp;  par%abund_N = 4.0e-5_wp;  par%abund_O = 3.3e-4_wp
  par%abund_Ne = 5.0e-5_wp;  par%abund_S = 9.0e-6_wp
  call species_setup(1)
  call milne_setup()

  call get_command_argument(1, fn)
  open(10, file=trim(fn), access='stream', form='unformatted', status='old')
  read(10) nleaf
  allocate(T(nleaf), ne(nleaf), nH(nleaf), V(nleaf), sh(nleaf), x(nleaf, 7, 5))
  read(10) T, ne, nH, V, sh
  do iel = 1, 5
     do j = 1, 7
        read(10) x(:, j, iel)
     end do
  end do
  close(10)

  allocate(wk(NTB, NSH, milne_metal_count()), lband(NSH, milne_metal_count()), &
           lwin(NW, NSH, milne_metal_count()))
  wk = 0;  lband = 0;  lwin = 0
  do k = 1, milne_metal_count()
     call species_transition(k, ie, it)
     iel = 0
     do j = 1, 5
        if (trim(ELN(j)) == trim(elem_name(ie))) iel = j
     end do
     if (iel == 0) cycle
     do il = 1, nleaf
        if (nH(il) <= 0.0_wp) cycle
        ish = nint(sh(il))
        w = ne(il)*V(il)*nH(il)*elem_abund(ie)*x(il, it+1, iel)
        lband(ish, k) = lband(ish, k) + w*milne_metal_band_coefficient(k, T(il))
        itb = min(max(int((log10(T(il)) - LT0)/DLT) + 1, 1), NTB)
        wk(itb, ish, k) = wk(itb, ish, k) + w
     end do
     do itb = 1, NTB
        Tc = 10.0_wp**(LT0 + (real(itb,wp) - 0.5_wp)*DLT)
        if (all(wk(itb, :, k) <= 0.0_wp)) cycle
        do iw = 1, NW
           lwin(iw, :, k) = lwin(iw, :, k) + wk(itb, :, k) &
              *milne_metal_energy_integral(k, Tc, WE(iw), WE(iw+1))
        end do
     end do
  end do

  write(*,'(a,es12.5)') 'channel 6 total (band coefficient at leaf T) [erg/s]: ', sum(lband)
  write(*,'(a,es12.5)') 'channel 6 total (windows, 0.01 dex T bins)   [erg/s]: ', sum(lwin)
  write(*,'(a)') 'windows [eV]: 13.598-23.338 | 23.338-24.587 | 24.587-34.856 | 34.856-54.418 | 54.418-150'
  do ish = 1, NSH
     write(*,'(/,a,i0,a,es11.4,a,5es11.3)') 'shell ', ish, ': L6 ', sum(lband(ish,:)), &
        '  windows', (sum(lwin(iw, ish, :)), iw = 1, NW)
     do k = 1, milne_metal_count()
        if (lband(ish, k) < 1.0e-3_wp*sum(lband(ish, :))) cycle
        call species_transition(k, ie, it)
        write(ion, '(a,i0)') trim(elem_name(ie)), it
        write(*,'(3x,a6,f8.3,a,es11.4,a,5es11.3)') ion, elem_eth(ie, it), ' eV  L ', &
           lband(ish, k), '  windows', (lwin(iw, ish, k), iw = 1, NW)
     end do
  end do
  call MPI_FINALIZE(ierr)
end program channel6_spectrum
