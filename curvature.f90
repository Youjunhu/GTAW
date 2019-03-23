subroutine curvature(mpoloidal,nflux,r_new,z_new,fpsi,dpsidra,jacobian,jacobian_th,rth,zth,bsq,b0_th,circumference,kappas,kappa_psi)
!calculate the norm curvature and geodesic curvature of the equilibrium magentic field
  use precision,only:p_
  use constants,only: one,two,twopi
  implicit none
  integer,intent(in)::mpoloidal,nflux
  real(p_),intent(in)::r_new(mpoloidal,nflux),z_new(mpoloidal,nflux)
  real(p_),intent(in):: jacobian(mpoloidal,nflux),jacobian_th(mpoloidal,nflux)
  real(p_),intent(in):: rth(mpoloidal,nflux),zth(mpoloidal,nflux)
  real(p_),intent(in):: bsq(mpoloidal,nflux),b0_th(mpoloidal,nflux)
  real(p_),intent(in):: fpsi(nflux),circumference(nflux),dpsidra(nflux)
  real(p_),intent(out):: kappas(mpoloidal,nflux),kappa_psi(mpoloidal,nflux) !kappas is the geodesic curvature
  real(p_):: kappas2(mpoloidal,nflux)
  integer:: i,j,k1,k2
  real(p_):: dth,theta_derivative_rth_over_zth(mpoloidal,nflux),rthth(mpoloidal,nflux), zthth(mpoloidal,nflux)

  do j=1,nflux
     do i=1,mpoloidal
        !for the equal arc Jacobian
        kappas2(i,j)=dpsidra(j)**3*fpsi(j)/bsq(i,j)**2/jacobian(i,j)**2 &
             & *(circumference(j)/twopi)**2*(-one/jacobian(i,j)**2*jacobian_th(i,j)) &
             & -dpsidra(j)*fpsi(j)**3/r_new(i,j)**3/bsq(i,j)**2/jacobian(i,j)*rth(i,j)
        !refer to /home/yj/theory/mhd.tm for the derivation of the above formula
     enddo
  enddo

  do i=1,mpoloidal
     do j=1,nflux
        !a better method of calculating geodesic curvature (using W. Deng's formula, refer to my notes mhd.pdf)
        kappas(i,j)=dpsidra(j)*fpsi(j)/(jacobian(i,j)*sqrt(bsq(i,j))**3)*b0_th(i,j)
     enddo
  enddo

  call plot_poloidal(mpoloidal,nflux,kappas,kappas2,'geodesic.txt')


  dth=twopi/(mpoloidal-1)
!!$  do j=1,nflux
!!$     do i=1,mpoloidal
!!$        k1=i+1
!!$        k2=i-1
!!$        if (i .eq. mpoloidal) k1=2 !deal with boundary points
!!$        if (i .eq. 1) k2=mpoloidal-1 !deal with boundary points
!!$        theta_derivative_rth_over_zth(i,j)=(rth(k1,j)/zth(k1,j)-rth(k2,j)/zth(k2,j))/(two*dth)
!!$     enddo
!!$  enddo

  do j=1,nflux
     do i=1,mpoloidal
        k1=i+1
        k2=i-1
        if (i .eq. mpoloidal) k1=2 !deal with boundary points
        if (i .eq. 1) k2=mpoloidal-1 !deal with boundary points
        rthth(i,j)=(rth(k1,j)-rth(k2,j))/(two*dth)
        zthth(i,j)=(zth(k1,j)-zth(k2,j))/(two*dth)
     enddo
  enddo

  !write(*,*) 'theta_derivative_rth_over_zth',theta_derivative_rth_over_zth

  do j=1,nflux
     do i=1,mpoloidal
        ! kappa_psi(i,j)= -r_new(i,j)/bsq(i,j)/jacobian(i,j)**3*zth(i,j)**2*theta_derivative_rth_over_zth(i,j) &
        !     & +fpsi(j)**2/bsq(i,j)/r_new(i,j)**2/jacobian(i,j)*zth(i,j)
        !the term containing theta_derivative_rth_over_zth needs checking.2013-8-11,
        !it turns out the above formula is not suitable for numerical calculation because it involves a 0/0 term
        ! the suitable formula is as follows, (refer to my note mhd.tm)
        kappa_psi(i,j)= -dpsidra(j)**3*r_new(i,j)/bsq(i,j)/jacobian(i,j)**3*(rthth(i,j)*zth(i,j)-rth(i,j)*zthth(i,j)) &
             & +dpsidra(j)*fpsi(j)**2/(bsq(i,j)*r_new(i,j)**2*jacobian(i,j))*zth(i,j)

     enddo
  enddo

end subroutine curvature
