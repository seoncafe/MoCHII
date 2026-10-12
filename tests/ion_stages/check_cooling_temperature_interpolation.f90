program check_cooling_temperature_interpolation
  use define,only: wp
  use nlevel_cooling_mod,only: reg,col_cooling
  use, intrinsic :: ieee_arithmetic,only: ieee_is_finite
  implicit none
  character(len=1024) :: input,output
  real(kind=wp),allocatable :: temperature(:)
  real(kind=wp) :: T,w1,wt,value
  integer :: unit,results,ncase,n,nquery,i,k,j,it
  call get_command_argument(1,input)
  call get_command_argument(2,output)
  open(newunit=unit,file=trim(input),status='old')
  open(newunit=results,file=trim(output),status='new')
  read(unit,*) ncase
  do k=1,ncase
    read(unit,*) n,nquery
    allocate(temperature(n),reg(1)%logT(n),reg(1)%lambda(n,1))
    do i=1,n
      read(unit,*) temperature(i),reg(1)%lambda(i,1)
    enddo
    reg(1)%logT=log10(temperature)
    do j=1,nquery
      read(unit,*) T
      it=1
      do while(it<n-1)
        if(T<temperature(it+1)) exit
        it=it+1
      enddo
      w1=(1/T-1/temperature(it))/(1/temperature(it+1)-1/temperature(it))
      wt=(log10(T)-reg(1)%logT(it))/(reg(1)%logT(it+1)-reg(1)%logT(it))
      value=col_cooling(1,it,1,w1,wt)
      if(.not.ieee_is_finite(value).or.value<0) stop 1
      write(results,'(i4,2(1x,es26.16e3))') k,T,value
    enddo
    deallocate(temperature,reg(1)%logT,reg(1)%lambda)
  enddo
  close(unit);close(results)
  print *, 'PASS: cooling interpolation finite and nonnegative'
end program check_cooling_temperature_interpolation
