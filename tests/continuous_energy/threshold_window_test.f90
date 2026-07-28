program threshold_window_test
!---------------------------------------------------------------------------
! Phase-0 gate for the physically important interval between the C II
! (24.383 eV) and He I (24.587 eV) photoionization thresholds.
!
! This test deliberately uses the production VFKY96 and He I cross-section
! routines.  The C II VFKY96 parameters are read from the production atomic
! registry file rather than duplicated here.
!---------------------------------------------------------------------------
  use define,      only : wp
  use photo_xsec,  only : sigma_vfky96, sigma_HeI
  implicit none

  real(kind=wp), parameter :: ETEST(4) = [24.38_wp, 24.40_wp, 24.58_wp, 24.60_wp]
  logical, parameter :: EXPECT_CII(4) = [.false., .true., .true., .true.]
  logical, parameter :: EXPECT_HEI(4) = [.false., .false., .false., .true.]
  character(len=512) :: atomic_file
  real(kind=wp) :: photo(8), eth_cii, sig_cii, sig_hei
  integer :: i, nfail

  call get_command_argument(1, atomic_file)
  if (len_trim(atomic_file) == 0) then
     write(*,'(a)') 'usage: threshold_window_test <element_c.txt>'
     error stop 2
  end if
  call read_cii_photo(trim(atomic_file), eth_cii, photo)

  nfail = 0
  write(*,'(a)') '=== C II / He I threshold-window gate ==='
  write(*,'(a,f10.6,a)') 'C II threshold from registry: ', eth_cii, ' eV'
  write(*,'(a)') ' E [eV]       sigma(C II)       sigma(He I)    expected'
  do i = 1, size(ETEST)
     sig_cii = cii_sigma(ETEST(i), photo)
     sig_hei = sigma_HeI(ETEST(i))
     write(*,'(f8.3,2es18.8,4x,l1,1x,l1)') ETEST(i), sig_cii, sig_hei, &
          EXPECT_CII(i), EXPECT_HEI(i)
     if ((sig_cii > 0.0_wp) .neqv. EXPECT_CII(i)) nfail = nfail + 1
     if ((sig_hei > 0.0_wp) .neqv. EXPECT_HEI(i)) nfail = nfail + 1
  end do

  ! Exact-edge convention: the production fits are active at E == E_th.
  if (cii_sigma(eth_cii, photo) <= 0.0_wp) nfail = nfail + 1
  if (sigma_HeI(24.587_wp) <= 0.0_wp) nfail = nfail + 1

  ! The carbon window itself: a 24.40-eV packet must contribute C II
  ! opacity/ionization while contributing exactly zero He I opacity.
  sig_cii = cii_sigma(24.40_wp, photo)
  sig_hei = sigma_HeI(24.40_wp)
  if (sig_cii <= 0.0_wp .or. sig_hei /= 0.0_wp) nfail = nfail + 1

  ! Guard against reading the wrong transition or a corrupted registry.
  if (abs(eth_cii - 24.383_wp) > 1.0e-10_wp) nfail = nfail + 1
  if (abs(sig_cii/1.0e-18_wp - 4.934570_wp) > 5.0e-3_wp) nfail = nfail + 1

  if (nfail /= 0) then
     write(*,'(a,i0)') 'THRESHOLD-WINDOW GATE: FAIL, checks failed = ', nfail
     error stop 1
  end if
  write(*,'(a)') 'THRESHOLD-WINDOW GATE: PASS'

contains

  real(kind=wp) function cii_sigma(E, p) result(sig)
    real(kind=wp), intent(in) :: E, p(8)
    sig = sigma_vfky96(E, p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))
  end function cii_sigma

  subroutine read_cii_photo(fname, eth, p)
    character(len=*), intent(in) :: fname
    real(kind=wp), intent(out) :: eth, p(8)
    character(len=512) :: line
    character(len=32) :: key
    integer :: unit, ios, transition
    logical :: got_eth, got_photo

    eth = 0.0_wp
    p = 0.0_wp
    transition = 0
    got_eth = .false.
    got_photo = .false.
    open(newunit=unit, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,'(2a)') 'cannot open ', trim(fname)
       error stop 2
    end if
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       line = adjustl(line)
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       read(line,*,iostat=ios) key
       if (ios /= 0) cycle
       select case (trim(key))
       case ('TRANSITION')
          read(line,*) key, transition
       case ('ETH')
          if (transition == 2) then
             read(line,*) key, eth
             got_eth = .true.
          end if
       case ('PHOTO')
          if (transition == 2) then
             read(line,*) key, p
             got_photo = .true.
          end if
       end select
    end do
    close(unit)
    if (.not. got_eth .or. .not. got_photo) then
       write(*,'(2a)') 'C II transition data not found in ', trim(fname)
       error stop 2
    end if
  end subroutine read_cii_photo

end program threshold_window_test
