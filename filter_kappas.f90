subroutine filter_theta(mpoloidal,nflux,kappas_in,kappas_out) 
!filter out some Fourier components to generate up-down asymmetric kappas. This routine can actually filter any two dimesion quantities, such as sigma_mu0(mpoloidal,nflux)
  use precision,only:p_
  use constants,only: zero,one,two,pi,twopi
  implicit none
  integer,intent(in)::mpoloidal,nflux
  real(p_),intent(in):: kappas_in(mpoloidal,nflux)
  real(p_),intent(out):: kappas_out(mpoloidal,nflux)

  integer,parameter:: n=50 !number of terms included in the Fourier series
  real(p_):: a(0:n),b(0:n)  !fourier coefficients
  real(p_):: dth,theta(mpoloidal),data(mpoloidal),output(mpoloidal) 
  real(p_):: sum1,sum2
  integer:: i,j,k


  dth=twopi/(mpoloidal-1)
  do i=1,mpoloidal
     theta(i)=zero+dth*(i-1)
  enddo


  do j=1,nflux

     do i=1,mpoloidal
        data(i)=kappas_in(i,j)
     enddo

     !calculate the Fourier coefficients
     do k=0,n
        sum1=zero
        sum2=zero
        do i=1,mpoloidal-1 !the integration in the Fourier coefficients
           sum1=sum1+data(i)*cos(k*theta(i))*dth
           sum2=sum2+data(i)*sin(k*theta(i))*dth
        enddo
        a(k)=sum1/pi
        b(k)=sum2/pi
     enddo

!!$     do k=0,6
!!$        write(*,*) k,a(k),b(k)
!!$     enddo

     b(1)=b(1)*(-120.0) !artificially change the coefficient, alter this to get different up-down asymmetricity

     do i=1,mpoloidal !reconstruct the original function
        output(i)=a(0)/two
        do k=1,n
           output(i)=output(i)+a(k)*cos(k*theta(i))+b(k)*sin(k*theta(i))
        enddo
     enddo

     do i=1,mpoloidal !store the filtered datas
        kappas_out(i,j)=output(i)
     enddo

  enddo

end subroutine filter_theta
