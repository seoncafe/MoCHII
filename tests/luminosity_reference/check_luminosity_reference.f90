program check_luminosity_reference
  use mpi
  use define, only : par, mpar, wp, int64
  use ion_band_mod, only : ion_setup, ion_Ltot, ion_lum, ion_L_ionizing, &
       ion_L_fuv, source_luminosity, band_packet_luminosity, band_packet_count, &
       LAUNCH_BAND_WHOLE, LAUNCH_BAND_IONIZING, LAUNCH_BAND_FUV
  use octree_mod, only : amr_grid
  implicit none
  character(len=512) :: input_file
  integer :: ierr, unit, ios, is, j
  namelist /parameters/ par

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpar%p_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpar%nproc, ierr)
  call get_command_argument(1, input_file)
  if (len_trim(input_file) == 0) error stop 'input file required'
  open(newunit=unit, file=trim(input_file), status='old', iostat=ios)
  if (ios /= 0) error stop 'cannot open input file'
  read(unit, nml=parameters, iostat=ios)
  if (ios /= 0) error stop 'cannot read parameters namelist'
  close(unit)
  amr_grid%xrange = 1.0_wp
  amr_grid%yrange = 1.0_wp
  amr_grid%zrange = 1.0_wp
  call ion_setup()
  if (mpar%p_rank == 0) then
     write(*,'(a,es26.17)') 'RESULT Ltot ', ion_Ltot
     write(*,'(a,es26.17)') 'RESULT binned ', sum(ion_lum)
     write(*,'(a,es26.17)') 'RESULT Lion ', ion_L_ionizing
     write(*,'(a,es26.17)') 'RESULT Lfuv ', ion_L_fuv
     write(*,'(a,es26.17)') 'RESULT parameter_luminosity ', par%luminosity
     write(*,'(a,3(es26.17,1x))') 'RESULT packet_luminosities ', &
          band_packet_luminosity(LAUNCH_BAND_WHOLE), &
          band_packet_luminosity(LAUNCH_BAND_IONIZING), &
          band_packet_luminosity(LAUNCH_BAND_FUV)
     write(*,'(a,3(i0,1x))') 'RESULT packet_counts ', &
          band_packet_count(LAUNCH_BAND_WHOLE), &
          band_packet_count(LAUNCH_BAND_IONIZING), &
          band_packet_count(LAUNCH_BAND_FUV)
     if (allocated(source_luminosity)) then
        write(*,'(a,i0)') 'RESULT source_count ', size(source_luminosity,2)
        do is = 1, size(source_luminosity,2)
           write(*,'(a,i0,1x,7(es26.17,1x))') 'RESULT source ', is, &
                (source_luminosity(j,is),j=1,7)
        end do
     else
        write(*,'(a,i0)') 'RESULT source_count ', 0
     end if
  end if
  call MPI_FINALIZE(ierr)
end program check_luminosity_reference
