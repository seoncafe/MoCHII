program fortran_spline_values
  use interpolation_mod, only: spline, splint
  implicit none
  real :: nodes(11), strength(11), second(11), temperature, value
  integer :: i, stat
  character(len=512) :: node_path, query_path, result_path

  call get_command_argument(1, node_path)
  call get_command_argument(2, query_path)
  call get_command_argument(3, result_path)
  open(10,file=trim(node_path),status='old',action='read')
  do i=1,11
     read(10,*) nodes(i), strength(i)
     nodes(i)=log10(nodes(i))
  end do
  close(10)
  call spline(nodes,strength,1.e30,1.e30,second)
  open(11,file=trim(query_path),status='old',action='read')
  open(12,file=trim(result_path),status='replace',action='write')
  do
     read(11,*,iostat=stat) temperature
     if (stat < 0) exit
     if (stat /= 0) error stop 'invalid temperature input'
     call splint(nodes,strength,second,log10(temperature),value)
     write(12,'(ES24.16,1X,ES24.16)') temperature,value
  end do
  close(11)
  close(12)
end program fortran_spline_values
