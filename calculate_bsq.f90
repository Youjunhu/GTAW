subroutine calculate_bsq(mpoloidal,nflux,r_new,z_new,psi_gradient_fcg,fpsi,bsq,&
     & b0_th,b_toroidal,bp_boundary_av)
  !calculate the magnetic field
  use precision,only:p_
  use constants,only:twopi,two
  implicit none
  integer,intent(in)::mpoloidal,nflux
  real(p_),intent(in)::r_new(mpoloidal,nflux),z_new(mpoloidal,nflux),psi_gradient_fcg(mpoloidal,nflux),fpsi(nflux)
  real(p_),intent(out)::bsq(mpoloidal,nflux),b0_th(mpoloidal,nflux),b_toroidal(mpoloidal,nflux)
  real(p_)::b0(mpoloidal,nflux),dth,bp_boundary_av
  integer:: i,j,k1,k2
  real(p_):: sum

  do i=1,mpoloidal
     do j=1,nflux
        bsq(i,j)=(fpsi(j)**2+psi_gradient_fcg(i,j)**2)/r_new(i,j)**2
     enddo
  enddo

  b0=sqrt(bsq)

  dth=twopi/(mpoloidal-1)
  do i=1,mpoloidal
     do j=1,nflux
        k1=i+1
        k2=i-1
        if (i .eq. mpoloidal) k1=2 !deal with boundary points
        if (i .eq. 1) k2=mpoloidal-1 !deal with boundary points
        b0_th(i,j)=(b0(k1,j)-b0(k2,j))/(two*dth)
     enddo
  enddo


  do i=1,mpoloidal
     do j=1,nflux
        b_toroidal(i,j)=fpsi(j)/r_new(i,j)
     enddo
  enddo

!---test, output averaged poloidal magentic field
  j=nflux
  sum=0.
  do i=1,mpoloidal
     sum=sum+psi_gradient_fcg(i,j)/r_new(i,j)
  enddo
 write(*,*) 'at surface nflux=',j,'averaged poloidal magentic field (Tesla)=',sum/mpoloidal

  j=nflux-2
  sum=0.
  do i=1,mpoloidal
     sum=sum+psi_gradient_fcg(i,j)/r_new(i,j)
  enddo
 write(*,*) 'at surface nflux=',j,'averaged poloidal magentic field (Tesla)=',sum/mpoloidal
 bp_boundary_av=sum/mpoloidal

 j=nflux-3
  sum=0.
  do i=1,mpoloidal
     sum=sum+psi_gradient_fcg(i,j)/r_new(i,j)
  enddo
 write(*,*) 'at surface nflux=',j,'averaged poloidal magentic field (Tesla)=',sum/mpoloidal
!--test end---


end subroutine calculate_bsq


subroutine surface_average(mpoloidal,nflux,bp,dl,weight,averaged) !not used presently
  use precision,only:p_
!  use constants,only:twopi
  implicit none
  integer,intent(in)::mpoloidal,nflux
  real(p_),intent(in):: bp(mpoloidal,nflux),dl(mpoloidal,nflux),weight(mpoloidal,nflux)
  real(p_),intent(out):: averaged(nflux)
  real(p_):: sum1,sum2
  integer:: i,j

!  do j=2,nflux-1 !calculate the surface averaged quantity: <B^2>
  do j=1,nflux !calculate the surface averaged quantity: <B^2>
     sum1=0.
     sum2=0.
     do i=1,mpoloidal-1
        sum1=sum1+weight(i,j)*dl(i,j)/bp(i,j)
        sum2=sum2+dl(i,j)/bp(i,j)
     enddo
     averaged(j)=sum1/sum2
  enddo

end subroutine surface_average
