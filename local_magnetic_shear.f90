subroutine local_magnetic_shear0(mpoloidal,nflux,r_new,fpsi,ffprime,psi_x_fcg,psi_z_fcg,&
     & psi_xx_fcg,psi_zz_fcg,psi_xz_fcg,local_shear)
  !S*grad_psi**2
  use precision,only:p_
  use constants,only: one,two,twopi,four,pi
  implicit none

  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: r_new(mpoloidal,nflux),fpsi(nflux),ffprime(nflux)
  real(p_),intent(in):: psi_x_fcg(mpoloidal,nflux),psi_z_fcg(mpoloidal,nflux),&
       & psi_xx_fcg(mpoloidal,nflux),psi_zz_fcg(mpoloidal,nflux),psi_xz_fcg(mpoloidal,nflux)
  real(p_),intent(out):: local_shear(mpoloidal,nflux) !,local_shear3(mpoloidal,nflux)
  real(p_):: fprime(nflux)
  integer:: i,j
  real(p_):: psi_grad_sq

  fprime=ffprime/fpsi

  do i=1,mpoloidal
     do j=1,nflux
        psi_grad_sq= psi_x_fcg(i,j)**2+psi_z_fcg(i,j)**2

        local_shear(i,j)=one/r_new(i,j)*(fpsi(j)*psi_zz_fcg(i,j)+ &
             & fprime(j)*psi_z_fcg(i,j)**2+fpsi(j)*psi_xx_fcg(i,j)+&
             & (fprime(j)*psi_x_fcg(i,j)**2*r_new(i,j)-fpsi(j)*psi_x_fcg(i,j))/r_new(i,j))&
             & -one/psi_grad_sq*fpsi(j)/r_new(i,j)*(four*psi_x_fcg(i,j)*psi_xz_fcg(i,j)*psi_z_fcg(i,j)&
             & +two*psi_z_fcg(i,j)**2*psi_zz_fcg(i,j)+two*psi_x_fcg(i,j)**2*psi_xx_fcg(i,j))

        local_shear(i,j)=local_shear(i,j)/r_new(i,j) 

     enddo
  enddo

end subroutine local_magnetic_shear0


subroutine local_magnetic_shear(mpoloidal,nflux,r_new,fpsi,ffprime,psi_x_fcg,psi_z_fcg,&
     & psi_xx_fcg,psi_zz_fcg,psi_xz_fcg,local_shear)
  !this subroutine calculate the local magnetic shear, refer to my notes for the calculating formula in cylindrical coordinates
  use precision,only:p_
  use constants,only: one,two,twopi,four,pi
  implicit none

  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: r_new(mpoloidal,nflux),fpsi(nflux),ffprime(nflux)
  real(p_),intent(in):: psi_x_fcg(mpoloidal,nflux),psi_z_fcg(mpoloidal,nflux),&
       & psi_xx_fcg(mpoloidal,nflux),psi_zz_fcg(mpoloidal,nflux),psi_xz_fcg(mpoloidal,nflux)
  real(p_),intent(out):: local_shear(mpoloidal,nflux) !,local_shear3(mpoloidal,nflux)
  real(p_):: fprime(nflux)
  integer:: i,j
  real(p_):: psi_grad_sq

  fprime=ffprime/fpsi

  do i=1,mpoloidal
     do j=1,nflux
    ! do j=2,nflux
        psi_grad_sq= psi_x_fcg(i,j)**2+psi_z_fcg(i,j)**2

        local_shear(i,j)=one/psi_grad_sq/r_new(i,j)*(fpsi(j)*psi_zz_fcg(i,j)+ &
             & fprime(j)*psi_z_fcg(i,j)**2+fpsi(j)*psi_xx_fcg(i,j)+&
             & (fprime(j)*psi_x_fcg(i,j)**2*r_new(i,j)-fpsi(j)*psi_x_fcg(i,j))/r_new(i,j))&
             & -one/psi_grad_sq**2*fpsi(j)/r_new(i,j)*(four*psi_x_fcg(i,j)*psi_xz_fcg(i,j)*psi_z_fcg(i,j)&
             & +two*psi_z_fcg(i,j)**2*psi_zz_fcg(i,j)+two*psi_x_fcg(i,j)**2*psi_xx_fcg(i,j))

        local_shear(i,j)=local_shear(i,j)/r_new(i,j) 
        
!!$        local_shear3(i,j)=(r_new(i,j)*fprime(j)*psi_grad_sq**2-&
!!$             & fpsi(j)*(four*r_new(i,j)*psi_z_fcg(i,j)*psi_x_fcg(i,j)*psi_xz_fcg(i,j)+&
!!$             & psi_z_fcg(i,j)**2*(r_new(i,j)*psi_zz_fcg(i,j)+psi_x_fcg(i,j)-r_new(i,j)*psi_xx_fcg(i,j))+&
!!$             & psi_x_fcg(i,j)**2*(-r_new(i,j)*psi_zz_fcg(i,j)+psi_x_fcg(i,j)+r_new(i,j)*psi_xx_fcg(i,j))))&
!!$             & /(r_new(i,j)**3*psi_grad_sq**2)
     enddo
     !local_shear(i,1)=two*local_shear(i,2)-local_shear(i,3) !linear interpolation to the first flux surface
  enddo

end subroutine local_magnetic_shear


subroutine local_magnetic_shear2(mpoloidal,nflux,dpsi,fpsi,jacobian,r_new,psi_dot_theta_psi,dpsidra,local_shear)
  !this subroutine calculate the local magnetic shear in the flux coordinates, refer to my notes (mhd.tm) for the calculating formula
  use precision,only:p_
  use constants,only: one,two,twopi,four,pi
  implicit none

  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: r_new(mpoloidal,nflux),jacobian(mpoloidal,nflux)
  real(p_),intent(in):: psi_dot_theta_psi(mpoloidal,nflux),dpsidra(nflux)
  real(p_),intent(in):: fpsi(nflux),dpsi
  real(p_),intent(out):: local_shear(mpoloidal,nflux)

  real(p_):: w1(mpoloidal,nflux),w2(mpoloidal,nflux),w1_psi(mpoloidal,nflux),w2_th(mpoloidal,nflux)
  integer:: i,j,k1,k2
  real(p_):: dth
  !real(p_)::sum1,sum2,sum3,sum4,sum5

  dth=twopi/(mpoloidal-1)

  do j=1,nflux
     do i=1,mpoloidal
        w1(i,j)=fpsi(j)*jacobian(i,j)/(r_new(i,j)**2*dpsidra(j))
        w2(i,j)=w1(i,j)*psi_dot_theta_psi(i,j) !here psi_dot_theta=(scalar product between grad_psi and grad_theta)/grad_psi^2
     enddo
  enddo

  do i=1,mpoloidal
!     do j=3,nflux-1
     do j=2,nflux-1
        w1_psi(i,j)=(w1(i,j+1)-w1(i,j-1))/(two*dpsi) !radial derivative
     enddo
     w1_psi(i,2)=two*(w1(i,3)-w1(i,2))/dpsi-w1_psi(i,3) !radial derivative

     w1_psi(i,nflux)=two*w1_psi(i,nflux-1)-w1_psi(i,nflux-2) !linear extrapolate
     w1_psi(i,1)=two*w1_psi(i,2)-w1_psi(i,3) !linear exptrapolate
  enddo

  do j=1,nflux
     do i=1,mpoloidal
        k1=i+1
        k2=i-1
        if (i .eq. mpoloidal) k1=2 !deal with boundary points
        if (i .eq. 1) k2=mpoloidal-1 !deal with boundary points
        w2_th(i,j)=(w2(k1,j)-w2(k2,j))/(two*dth)
     enddo
  enddo

  do i=1,mpoloidal
!     do j=2,nflux
     do j=1,nflux
        local_shear(i,j)=(w1_psi(i,j)+w2_th(i,j))/jacobian(i,j)
     enddo
     !local_shear(i,1)=two*local_shear(i,2)-local_shear(i,3) !linear extrapolate
  enddo

  !write(*,*) 'local_shear',local_shear

!!$  do i=1,mpoloidal
!!$     local_shear(i,1)=0.
!!$     local_shear(i,nflux)=0.
!!$  enddo
  !write(*,*) 'local_shear=', w1,w2
!!$  open(111,file='local_magnetic_shear.txt')
!!$  do j=1,nflux
!!$     do i=1,mpoloidal
!!$        if (isnan(local_shear(i,j))) then
!!$           write(*,*) 'local_shear',' i=',i,'j=',j
!!$           stop 'a NaN is found'
!!$        endif
!!$        write(111,*) i, 0.+dth*(i-1), local_shear(i,j)
!!$     enddo
!!$     write(111,*)
!!$     write(111,*)
!!$  enddo
!!$  close(111)

!!$  open(12345,file='global_magnetic_shear2.txt')
!!$  do j=2,nflux
!!$     sum1=0.
!!$     sum2=0.
!!$     sum3=0.
!!$     sum4=0.
!!$     sum5=0.
!!$     do i=1,mpoloidal-1
!!$!        sum1=sum1+local_shear(i,j)*dth
!!$        sum2=sum2+w1_psi(i,j)*dth
!!$ !       sum3=sum3+w2_th(i,j)/jacobian(i,j)*dth
!!$        sum4=sum4+local_shear(i,j)*jacobian(i,j)*dth
!!$        sum5=sum5+jacobian(i,j)*dth
!!$     enddo
!!$     write(12345,*) j,-sum4/twopi !the results will be compared with qprime
!!$  enddo
!!$  close(12345)
!!$
!!$do j=1,nflux
!!$write(*,*) w1(10,j),w2(10,j)
!!$enddo

end subroutine local_magnetic_shear2

subroutine benchmark_local_gloabal_shear(mpoloidal,nflux,jacobian,qprime,dpsidra,local_shear,filename)
  use precision,only:p_
  use constants,only: two,twopi
  implicit none
   integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in)::jacobian(mpoloidal,nflux), local_shear(mpoloidal,nflux),qprime(nflux),dpsidra(nflux)
  character(len=*),intent(in):: filename
  real(p_):: sum,dth
  integer:: i,j

  dth=twopi/(mpoloidal-1)
  open(12345,file=filename)
!  do j=2,nflux
  do j=1,nflux
     sum=0.
     do i=1,mpoloidal-1
        !sum=sum+local_shear(i,j)*jacobian(i,j)*dth
         sum=sum+(local_shear(i,j)*jacobian(i,j)+local_shear(i+1,j)*jacobian(i+1,j))/two*dth
     enddo
     write(12345,*) j,-sum/twopi,qprime(j)*dpsidra(j) !the latter is dq/dra. the two quantities should agree with each other
  enddo
  close(12345)
end subroutine benchmark_local_gloabal_shear









