

!!$subroutine fast_ions_pressure(nflux,mpoloidal,psival,b,b_axis)
!!$!use spherical coordinates in velocity space, for axisymmetrical (about magnetic field) velcoity distriubtion
!!$!for mega code
!!$  use precision,only:p_
!!$  use constants,only:two,twopi,pi
!!$  use fast_ions,only: mass,kev,me
!!$  implicit none
!!$  integer,intent(in):: nflux,mpoloidal
!!$  real(p_),intent(in):: psival(nflux),b_axis,b(mpoloidal,nflux) !psival is the normalized poloidal magnetic flux
!!$  !real(p_),intent(out):: ppara
!!$  real(p_)::ppara0,ppara(mpoloidal,nflux),pperp(mpoloidal,nflux)
!!$  real(p_),parameter::  psi_scale=0.4_p_, energy_birth=50._p_ !in unti of keV
!!$  real(p_),parameter:: te=2._p_   !electron temperature in unit of kev
!!$  real(p_)::vcrit,vbirth,vmax,deltav
!!$  real(p_),parameter:: clambda0=0.5_p_,dclambda=0.2_p_
!!$  integer,parameter::m=100,n=100
!!$  real(p_):: theta(m),v(n),dv,dtheta,dvol
!!$  real(p_):: clambda,kenergy,mu
!!$  integer:: i,j,ii,jj
!!$
!!$  vbirth=sqrt(two*energy_birth*kev/mass)
!!$  vcrit=sqrt(2.*te*kev/me)*(3*sqrt(pi)/(4.*1836))**0.3333 !refer to the formula in my notes
!!$  vmax=vbirth*1.2_p_
!!$  deltav=0.05*vbirth
!!$
!!$  write(*,*) 'vbirth=',vbirth,'vcrit=',vcrit
!!$  dtheta=pi/(m-1)
!!$  do i=1,m
!!$     theta(i)=0._p_+dtheta*(i-1)
!!$  enddo
!!$
!!$  dv=vmax/(n-1)
!!$  do j=1,n
!!$     v(j)=0._p_+dv*(i-1)
!!$  enddo
!!$
!!$  !---------
!!$  ppara0=0._p_
!!$  do i=1,m
!!$     do j=1,n
!!$        kenergy=mass*v(j)**2/two
!!$        mu=mass*(v(j)*sin(theta(i)))**2/(two*b_axis)
!!$        clambda=mu/kenergy*b_axis
!!$        dvol=twopi*v(j)**2*sin(theta(i))*dtheta*dv !volume element in velocity space
!!$        ppara0=ppara0+mass*(v(j)*cos(theta(i)))**2*f(0._p_,v(j),clambda)*dvol
!!$     enddo
!!$  enddo
!!$  write(*,*) 'ppara0=',ppara0 !this will be used as a normalzing factor
!!$  !---------------
!!$
!!$  do ii=1,mpoloidal
!!$     do jj=1,nflux
!!$        ppara(ii,jj)=0._p_
!!$        pperp(ii,jj)=0._p_
!!$        do i=1,m
!!$           do j=1,n
!!$              kenergy=mass*v(j)**2/two
!!$              mu=mass*(v(j)*sin(theta(i)))**2/(two*b(ii,jj))
!!$              clambda=mu/kenergy*b_axis
!!$              dvol=twopi*v(j)**2*sin(theta(i))*dv*dtheta !volume element in velocity space
!!$              ppara(ii,jj)=ppara(ii,jj)+mass*(v(j)*cos(theta(i)))**2*f(psival(jj),v(j),clambda)*dvol
!!$              pperp(ii,jj)=pperp(ii,jj)+mass*(v(j)*sin(theta(i)))**2*f(psival(jj),v(j),clambda)*dvol
!!$           enddo
!!$        enddo
!!$     enddo
!!$  enddo
!!$
!!$  ppara=ppara/ppara0
!!$  pperp=pperp/ppara0
!!$
!!$  call plot_poloidal(mpoloidal,nflux,ppara,pperp,'ppara_fast_ions')
!!$contains
!!$  function f(psi,v,clambda) result (z)
!!$    implicit none
!!$    real(p_):: psi,v,clambda,z
!!$
!!$    z=exp(-psi/psi_scale)/(v**3 + vcrit**3) &
!!$         & *0.50d0*erfc((v-vbirth)/deltav)*exp(-(clambda-clambda0)**2/dclambda**2)
!!$  end function f
!!$end subroutine 
