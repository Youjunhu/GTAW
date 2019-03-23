subroutine calculate_jacobian(mpoloidal,nflux,r_new,z_new,ra,dpsidra,poloidal_angle_type,psi_gradient_fcg,&
     & circumference,r_grad_psi_s,jacobian,jacobian_th)
  !refer to my notes (tokamake_equilibrium.tm) for the details about the form of the Jacobian
  use precision,only:p_
  use constants,only: zero,one,two,twopi
  implicit none
  integer,intent(in)::mpoloidal,nflux
  real(p_),intent(in)::r_new(mpoloidal,nflux),z_new(mpoloidal,nflux),psi_gradient_fcg(mpoloidal,nflux)!,sign_gradient_psi
  real(p_),intent(in)::circumference(nflux),r_grad_psi_s(nflux),ra(nflux),dpsidra(nflux)
  character(100),intent(in):: poloidal_angle_type
  real(p_),intent(out):: jacobian(mpoloidal,nflux), jacobian_th(mpoloidal,nflux)
!  real(p_):: s(nflux),hfunc(mpoloidal,nflux) !hfunc is a function in the expression of jacobian
  integer:: i,j,k1,k2
  real(p_)::dtheta,sign_jacobian

  sign_jacobian=1._p_
  if(ra(nflux)>ra(1)) sign_jacobian=-1._p_

  if(trim(poloidal_angle_type).eq.'equal-arc') then   !equal arc length theta coordinate, and the radial coordinate here is chosen to be equal to A_phi*R (GS poloidal flux function)
     do i=1,mpoloidal
        do j=1,nflux         !     do j=2,nflux
           jacobian(i,j)=sign_jacobian*r_new(i,j)/psi_gradient_fcg(i,j)*circumference(j)/twopi*abs(dpsidra(j))
        enddo
     enddo

  else if(trim(poloidal_angle_type).eq.'straight-line') then 
     do i=1,mpoloidal
        do j=1,nflux
           jacobian(i,j)=sign_jacobian*r_new(i,j)**2*r_grad_psi_s(j)*dpsidra(j)/twopi
        enddo
     enddo
  else
     stop "please choose poloidal_angle_type between 'equal-arc' and 'straight-line'"
  endif

  !calculate jacobian_th:
  dtheta=twopi/(mpoloidal-1)
  ! do j=2,nflux
  do j=1,nflux
     do i=1,mpoloidal
        k1=i+1
        k2=i-1
        if (i .eq. mpoloidal) k1=2 !deal with boundary points
        if (i .eq. 1) k2=mpoloidal-1 !deal with boundary points
        jacobian_th(i,j)=(jacobian(k1,j)-jacobian(k2,j))/(two*dtheta)
     enddo
  enddo
end subroutine calculate_jacobian


subroutine calculate_jacobian_directly(m,n,r,z,dtheta,dpsi,jacobian,psi_dot_theta,rpsi,rth,zpsi,zth,jacobian_th)
  !wrapper of subroutine partial_derivative
  !to output selected quantities, such as jacobian, jacobian_th, etc.
  use precision,only:p_
  use constants,only:zero,one,two
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  real(p_),intent(in):: dtheta, dpsi
  real(p_),intent(out):: jacobian(m,n),psi_dot_theta(m,n),jacobian_th(m,n),rth(m,n),zth(m,n)
  real(p_):: rpsi(m,n),zpsi(m,n)
  integer:: i,j,k1,k2

  call partial_derivative(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacobian) 

  do i=1,m
!     do j=2,n
     do j=1,n
        !here psi_dot_theta=(scalar product between grad_psi and grad_theta)/grad_psi^2
        psi_dot_theta(i,j)=-(zth(i,j)*zpsi(i,j)+rth(i,j)*rpsi(i,j))/(zth(i,j)**2+rth(i,j)**2)
     enddo
  enddo

  do j=1,n
     do i=1,m
        k1=i+1
        k2=i-1
        if (i .eq. m) k1=2 !deal with boundary points
        if (i .eq. 1) k2=m-1 !deal with boundary points
        jacobian_th(i,j)=(jacobian(k1,j)-jacobian(k2,j))/(two*dtheta)
     enddo
  enddo

end subroutine calculate_jacobian_directly


subroutine partial_derivative(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacob)

  !calculate the partial derivative of R and Z with respect to theta and psi
  !jacob is also calculated in this subroutine
  use precision,only:p_
  use constants,only:zero,one,two,twopi,one_half
  implicit none

  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: rpsi(m,n),rth(m,n),zpsi(m,n),zth(m,n),jacob(m,n)
  real(p_):: tmp0
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        rpsi(i,j)=(r(i,j+1)-r(i,j-1))/(two*dpsi)
        zpsi(i,j)=(z(i,j+1)-z(i,j-1))/(two*dpsi)
     enddo

     !use linear interpolation to get the value  j=n
     tmp0=(r(i,n)-r(i,n-1))/dpsi
     rpsi(i,n)=two*tmp0-rpsi(i,n-1)
     tmp0=(z(i,n)-z(i,n-1))/dpsi
     zpsi(i,n)=two*tmp0-zpsi(i,n-1)

     !use linear interpolation to get the value j=1
     tmp0=(r(i,2)-r(i,1))/dpsi
     rpsi(i,1)=two*tmp0-rpsi(i,2)

     tmp0=(z(i,2)-z(i,1))/dpsi
     zpsi(i,1)=two*tmp0-zpsi(i,2)

  enddo

  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        rth(i,j)= (r(i+1,j)-r(i-1,j))/(two*dtheta)
        zth(i,j)=(z(i+1,j)-z(i-1,j))/(two*dtheta)
     enddo

     !use peroidic property of r and z to calculate the partial derivative for boundary points at theta=0 and 2pi
     rth(1,j)=(r(2,j)-r(m-1,j))/(two*dtheta)
     zth(1,j)=(z(2,j)-z(m-1,j))/(two*dtheta)
     rth(m,j)=rth(1,j)
     zth(m,j)=zth(1,j)
  enddo

  !calculate the Jacobian:
  do i=1,m
!          do j=2,n !the jacobian at the magnetic axis is zero
     do j=1,n
        jacob(i,j)=r(i,j)*(rth(i,j)*zpsi(i,j)-rpsi(i,j)*zth(i,j))   !Jacobain of coordinate system (psi,theta,fai)
     enddo
     !jacob(i,n)=two*jacob(i,n-1)-jacob(i,n-2)
     !use linear interpolation to get the value of jacobian at the magnetic surface near the magnetic axis
  enddo

  !write(*,*) jacob(1:m,5)
end subroutine partial_derivative



subroutine safety_factor2(m,n,r,z,jacob,g,dpsidra)
  use precision,only:p_
  use constants,only:one,two,zero,twopi
  use radial_module,only: pfn,tfn
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n),jacob(m,n)
  real(p_),intent(in):: g(n),dpsidra(n)
  real(p_)::q(n),sum0,dth
  integer:: i,j

  dth=twopi/(m-1)

!  do j=3,n-2
  do j=1,n
     sum0=zero
     do i=1,m-1
        sum0=sum0+jacob(i,j)/r(i,j)**2*dth
     enddo
     q(j)=-g(j)/dpsidra(j)*sum0/twopi
  enddo

!!$   q(2)=two*q(3)-q(4) !use linear expolation to obtain q
!!$   q(1)=two*q(2)-q(3) !use linear expolation to obtain q 
!!$
!!$   q(n-1)=two*q(n-2)-q(n-3)
!!$   q(n)=two*q(n-1)-q(n-2)
  open(111,file='q.txt')
  do j=1,n
     write(111,*) sqrt(pfn(j)),sqrt(tfn(j)),q(j),g(j)
  enddo
  close(11)
end subroutine safety_factor2




!!$subroutine calculate_psi_dot_theta(mpoloidal,nflux,rpsi,rth,zpsi,zth,circumference,psi_dot_theta)
!!$!calculate psi_dot_theta, which is defined by (scalar product between grad_psi and grad_theta)/grad_psi^2
!!$  use precision,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  integer,intent(in)::mpoloidal,nflux
!!$  real(p_),intent(in)::rpsi(mpoloidal,nflux),rth(mpoloidal,nflux),zpsi(mpoloidal,nflux),zth(mpoloidal,nflux)
!!$  real(p_),intent(in)::circumference(nflux)
!!$  real(p_),intent(out):: psi_dot_theta(mpoloidal,nflux) 
!!$  integer:: i,j
!!$  
!!$  do i=1,mpoloidal
!!$     do j=1,nflux
!!$        psi_dot_theta(i,j)=-(zth(i,j)*zpsi(i,j)+rth(i,j)*rpsi(i,j))/circumference(j)**2*twopi**2 
!!$        !the above formula is for equal arc length Jacobian, refer to my notes for the formual
!!$     enddo
!!$  enddo
!!$
!!$end subroutine calculate_psi_dot_theta

subroutine calculate_psi_dot_theta(mpoloidal,nflux,rpsi,rth,zpsi,zth,jacobian,r_new,psi_dot_theta,psi_dot_theta_psi)
!calculate psi_dot_theta, which is defined by (scalar product between grad_psi and grad_theta)
  use precision,only:p_
  use constants,only:twopi
  implicit none
  integer,intent(in)::mpoloidal,nflux
  real(p_),intent(in)::rpsi(mpoloidal,nflux),rth(mpoloidal,nflux),zpsi(mpoloidal,nflux),zth(mpoloidal,nflux)
  real(p_),intent(in)::jacobian(mpoloidal,nflux),r_new(mpoloidal,nflux)
  real(p_),intent(out):: psi_dot_theta(mpoloidal,nflux) ,psi_dot_theta_psi(mpoloidal,nflux) 
  integer:: i,j
  
  do i=1,mpoloidal
     do j=1,nflux
        psi_dot_theta(i,j)=-(zth(i,j)*zpsi(i,j)+rth(i,j)*rpsi(i,j))*r_new(i,j)**2/jacobian(i,j)**2
        psi_dot_theta_psi(i,j)=-(zth(i,j)*zpsi(i,j)+rth(i,j)*rpsi(i,j))/(zth(i,j)**2+rth(i,j)**2)
     enddo
  enddo

end subroutine calculate_psi_dot_theta
