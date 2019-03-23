subroutine calculate_weight_functions(mpoloidal,nflux,r_new,circumference,jacobian, &
     & psi_gradient_fcg,psi_x_fcg,psi_z_fcg,psi_xx_fcg,psi_zz_fcg,psi_xz_fcg,psi_dot_theta, &
     & bsq,kappas,kappa_psi,sigma_mu0,local_shear, &
     & fpsi,ffprime,qpsi,qprime,dra,dpsidra, &
     & w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,w23,w24)
  !this subroutine calculates the weight functions appearing in the Fourier integrals
  use precision,only:p_
  use constants,only: one,two,twopi,four,pi
  use flux_grids,only: delta_q !as output
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: r_new(mpoloidal,nflux),circumference(nflux),jacobian(mpoloidal,nflux)
  real(p_),intent(in):: psi_gradient_fcg(mpoloidal,nflux), psi_dot_theta(mpoloidal,nflux)
  real(p_),intent(in):: psi_x_fcg(mpoloidal,nflux),psi_z_fcg(mpoloidal,nflux),psi_xx_fcg(mpoloidal,nflux)
  real(p_),intent(in):: psi_zz_fcg(mpoloidal,nflux),psi_xz_fcg(mpoloidal,nflux)
  real(p_),intent(in):: bsq(mpoloidal,nflux),kappas(mpoloidal,nflux),kappa_psi(mpoloidal,nflux)
  real(p_),intent(in):: sigma_mu0(mpoloidal,nflux),local_shear(mpoloidal,nflux)
  real(p_),intent(in):: fpsi(nflux),ffprime(nflux),qpsi(nflux),qprime(nflux),dpsidra(nflux)
  real(p_),intent(in):: dra !interval of GS psi between two neighbour magnetic surfaces
  real(p_),intent(out)::w1(mpoloidal,nflux),w2(mpoloidal,nflux),w3(mpoloidal,nflux),w4(mpoloidal,nflux)
  real(p_),intent(out)::w5(mpoloidal,nflux),w6(mpoloidal,nflux),w7(mpoloidal,nflux)
  real(p_),intent(out)::w8(mpoloidal,nflux),w9(mpoloidal,nflux),w10(mpoloidal,nflux),w11(mpoloidal,nflux)
  real(p_),intent(out)::w12(mpoloidal,nflux),w13(mpoloidal,nflux),w14(mpoloidal,nflux),w15(mpoloidal,nflux)
  real(p_),intent(out)::w16(mpoloidal,nflux),w17(mpoloidal,nflux),w18(mpoloidal,nflux),w19(mpoloidal,nflux)
  real(p_),intent(out)::w20(mpoloidal,nflux),w21(mpoloidal,nflux),w22(mpoloidal,nflux)
  real(p_),intent(out)::w23(mpoloidal,nflux),w24(mpoloidal,nflux)

  real(p_):: w2_tmp(mpoloidal,nflux),w4_tmp(mpoloidal,nflux),w11_tmp(mpoloidal,nflux),w21_tmp(mpoloidal,nflux)
  real(p_):: w15_tmp(mpoloidal,nflux)
  integer:: i,j,k1,k2
  real(p_):: dth
  real(p_):: a(mpoloidal,nflux),a_psi(mpoloidal,nflux)
  real(p_):: delta_q_psi(mpoloidal,nflux),delta_q_th(mpoloidal,nflux)
  !  real(p_):: delta_q_psi2(mpoloidal,nflux),delta_q_th2(mpoloidal,nflux)

  dth=twopi/(mpoloidal-1) !interval of poloidal grids


  w1=one/(jacobian**2*bsq)
  w2_tmp=one/(jacobian*bsq)
  w3=psi_gradient_fcg**2/(jacobian**2*bsq)
  w4_tmp= psi_gradient_fcg**2/(jacobian*bsq)
  w5=psi_gradient_fcg**2/bsq
  w6=kappas
  w7=one/bsq
  w8=one/(bsq*jacobian)
  w9=sigma_mu0/jacobian
  w10=local_shear/(bsq*jacobian)
  w11_tmp=local_shear/bsq
  w12=kappa_psi/psi_gradient_fcg**2

  w13=psi_gradient_fcg**2
  w15=kappa_psi

 w16= psi_x_fcg/r_new+(psi_xx_fcg+psi_zz_fcg) &
       & -(two*psi_x_fcg**2*psi_xx_fcg+four*psi_x_fcg*psi_z_fcg*psi_xz_fcg+two*psi_z_fcg**2*psi_zz_fcg) &
       & /psi_gradient_fcg**2

        w17=(local_shear-bsq*sigma_mu0)*local_shear/bsq
        w18=(local_shear-bsq*sigma_mu0)*psi_gradient_fcg**2/(bsq*jacobian)
        w19=psi_dot_theta !here psi_dot_theta=(scalar product between grad_psi and grad_theta)
 
        w21_tmp=one/(jacobian*psi_gradient_fcg**2)
        w22=one/(jacobian**2) !for general Jacobian
        w23=psi_gradient_fcg**2/(jacobian*bsq)
        w24=psi_gradient_fcg**2*kappas

  do j=1,nflux
     do i=1,mpoloidal
        k1=i+1
        k2=i-1
        if (i .eq. mpoloidal) k1=2 !deal with boundary points when calculating derivatives over poloidal angle
        if (i .eq. 1) k2=mpoloidal-1 !deal with boundary points
        w2(i,j)=(w2_tmp(k1,j)-w2_tmp(k2,j))/(two*dth)/jacobian(i,j)
        w4(i,j)=(w4_tmp(k1,j)-w4_tmp(k2,j))/(two*dth)/jacobian(i,j)
        w11(i,j)=(w11_tmp(k1,j)-w11_tmp(k2,j))/(two*dth)/jacobian(i,j)
        w21(i,j)=psi_gradient_fcg(i,j)**2/jacobian(i,j)*(w21_tmp(k1,j)-w21_tmp(k2,j))/(two*dth)
     enddo
  enddo

  w14=w2*psi_gradient_fcg**2

  call toroidal_angle_factor(mpoloidal,nflux,dra,r_new,jacobian,qpsi,fpsi,ffprime,qprime,dpsidra,delta_q,delta_q_th,delta_q_psi)

  do i=1,mpoloidal
     do j=1,nflux
        w20(i,j)=-delta_q_psi(i,j)*psi_gradient_fcg(i,j)**2/dpsidra(j)**2-delta_q_th(i,j)*psi_dot_theta(i,j)
     enddo
  enddo

end subroutine calculate_weight_functions


subroutine toroidal_angle_factor(m,n,dra,r,jacobian,qpsi,fpsi,ffprime,qprime,dpsidra,delta_q,delta_q_theta,delta_q_psi)
  !this subroutine calculates the toroidal angle shift factor delta*q and its derivatives, where delta is introduced in the definition of the generalized toroidal angle: zeta=fai-q*delta
  !delta_q_theta is the derivative of (delta*q) with respect to theta, !delta_q_psi is the derivative of (delta*q) with respect to psi
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: dra
  real(p_),intent(in):: r(m,n),jacobian(m,n)
  real(p_),intent(in):: qpsi(n),fpsi(n),ffprime(n),qprime(n),dpsidra(n)
  real(p_),intent(out):: delta_q(m,n),delta_q_theta(m,n),delta_q_psi(m,n) !delta_q_theta is the derivative of (delta*q) with respect to theta, !delta_q_psi is the derivative of (delta*q) with respect to psi
  real(p_):: delta_q_part1(m,n),delta_q_part2(m,n)
  !  real(p_):: dpsidradra(n)
  real(p_):: tmp(n)
  integer:: i,j
  real(p_):: dth,sum
  real(p_):: theta(m),a(m,n),a_psi(m,n)

  dth=twopi/(m-1) !polodal (theta) grid interval, uniform grid is assumed
  do i=1,m
     theta(i)=0._p_+dth*(i-1)
  enddo

  do i=1,m
     do j=1,n
        delta_q_theta(i,j)=-fpsi(j)*jacobian(i,j)/(dpsidra(j)*r(i,j)**2)-qpsi(j) !delta_q_theta is the derivative of (delta*q) with respect to theta
     enddo
  enddo

  do j=1,n
     a(1,j)=0.
     do i=2,m
        a(i,j)=a(i-1,j)+(jacobian(i-1,j)/r(i-1,j)**2+jacobian(i,j)/r(i,j)**2)/two*dth !a(i,j) is the variable upper-limit poloidal integration of jacobian/R**2
     enddo
  enddo

  do i=1,m
     do j=2,n-1
        a_psi(i,j)=(a(i,j+1)-a(i,j-1))/(two*dra)
     enddo
     a_psi(i,1)=two*a_psi(i,2)-a_psi(i,3)
     a_psi(i,n)=two*a_psi(i,n-1)-a_psi(i,n-2)
  enddo

!!$  do j=2,n-1
!!$     dpsidradra(j)=(dpsidra(j+1)-dpsidra(j-1))/(two*dra)
!!$  enddo
!!$     dpsidradra(1)=two*dpsidradra(2)-dpsidradra(3)
!!$     dpsidradra(n)=two*dpsidradra(n-1)-dpsidradra(n-2)

 do j=2,n-1
     tmp(j)=(fpsi(j+1)/dpsidra(j+1)-fpsi(j-1)/dpsidra(j-1))/(two*dra) !d(fpsi/dpsidra)/dra
  enddo
   tmp(1)=two*tmp(2)-tmp(3)
   tmp(n)=two*tmp(n-1)-tmp(n-2)

  do i=1,m
     do j=1,n !delta_q_psi is the derivative of (delta*q) with respect to psi
!!$        delta_q_psi(i,j)=-(ffprime(j)/fpsi(j)/dpsidra(j)+fpsi(j)*(-one/dpsidra(j)**2)*dpsidradra(j))*a(i,j) &
!!$             & -fpsi(j)/dpsidra(j)*a_psi(i,j)-qprime(j)*dpsidra(j)*theta(i)
 delta_q_psi(i,j)=-tmp(j)*a(i,j) &
             & -fpsi(j)/dpsidra(j)*a_psi(i,j)-qprime(j)*dpsidra(j)*theta(i)
     enddo
  enddo

  do j=1,n
     do i=1,m
delta_q_part1(i,j)=-fpsi(j)*a(i,j)/dpsidra(j)
delta_q_part2(i,j)=qpsi(j)*theta(i) 
     enddo
  enddo

   delta_q=delta_q_part1-delta_q_part2 !delta_q is delta*q

  call plot_poloidal(m,n,delta_q_part1,delta_q_part2,'delta_q.txt')

end subroutine toroidal_angle_factor

