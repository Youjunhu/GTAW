program read
  use precision,only:p_
  implicit none
  integer,parameter:: max_num=3000
  real(p_):: tmp_density(max_num)
  integer:: j,ndata
  open(11,file='ne_g038300.03900')
  do j=1,max_num
     read(11,*,end=111) tmp_density(j) !this should be the electron number density
     !number_density_old(j)=0.2d0 !set analytical form for the electron number density profile
  enddo
111 close(11)

  ndata=j-1
  write(*,*) 'number of data of the density radial profile=',ndata

  open(11,file='ne_g048916.04500')
  do j=1,ndata
     write(11,*) tmp_density(j)/tmp_density(1)*4.0d19
  enddo
  close(11)
end program read

