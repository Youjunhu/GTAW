subroutine interpolate_flux_function(nx,pfn_nx,psival_nx,fpsi_nx,qpsi_nx,press_nx,pprime_nx,ffprime_nx,tfn_nx, &
     & nflux,pfn,fpsi,qpsi,pressure,pprime,ffprime,qprime,tfn)
  use precision,only:p_
  use constants,only:zero,one,two
  implicit none
  integer,intent(in):: nx,nflux
  real(p_),intent(in):: pfn_nx(nx),psival_nx(nx),fpsi_nx(nx),qpsi_nx(nx),press_nx(nx),pprime_nx(nx), &
       & ffprime_nx(nx),tfn_nx(nx),pfn(nflux)
  real(p_),intent(out):: fpsi(nflux),qpsi(nflux),pressure(nflux),pprime(nflux), &
       & ffprime(nflux),qprime(nflux),tfn(nflux)
  real(p_):: qprime_nx(nx)
  real(p_):: tmp_y2(nx),y_tmp,q_y2(nx)
  real(p_):: dpsi_nx,q95
  integer:: j

  !interpolate the array qpsi, fpsi, press, and their radial derivative to the new psival grid points
  call spline(pfn_nx,fpsi_nx,nx,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do j=1,nflux
     call splint(pfn_nx,fpsi_nx,tmp_y2,nx,pfn(j),y_tmp) 
     fpsi(j)=y_tmp !to get fpsi array corresponding to the new radial grid points
  enddo

  call spline(pfn_nx,qpsi_nx,nx,2.d30,2.d30,q_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do j=1,nflux
     call splint(pfn_nx,qpsi_nx,q_y2,nx,pfn(j),y_tmp) 
     qpsi(j)=y_tmp !to get qpsi array corresponding to the new radial grid points
  enddo

  call splint(pfn_nx,qpsi_nx,tmp_y2,nx,0.95_p_,y_tmp) !get the value of q at 95% poloidal magnetic flux
  q95=y_tmp !the value of the safety factor at the 95% poloidal magnetic flux
  write(*,*) 'q(0)=',qpsi_nx(1),'q95=', q95

  call spline(pfn_nx,press_nx,nx,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do j=1,nflux
     call splint(pfn_nx,press_nx,tmp_y2,nx,pfn(j),y_tmp) 
     pressure(j)=y_tmp !to get pressure array corresponding to the new radial grid points
     !pressure(j)=1.*10**9. !test high/zero beta limit
  enddo

  call spline(pfn_nx,pprime_nx,nx,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do j=1,nflux
     call splint(pfn_nx,pprime_nx,tmp_y2,nx,pfn(j),y_tmp) 
     pprime(j)=y_tmp !to get pprime array corresponding to the new radial grid points
     !pprime(j)=0. !test zero beta limit
  enddo

  call spline(pfn_nx,ffprime_nx,nx,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do j=1,nflux
     call splint(pfn_nx,ffprime_nx,tmp_y2,nx,pfn(j),y_tmp) 
     ffprime(j)=y_tmp !to get ffprime array corresponding to the new radial grid points
  enddo

  dpsi_nx=(psival_nx(nx)-psival_nx(1))/(nx-1)
  do j=2,nx-1 !radial derivative of safety factor
     qprime_nx(j)=(qpsi_nx(j+1)-qpsi_nx(j-1))/(two*dpsi_nx)
  enddo
  qprime_nx(nx)=two*qprime_nx(nx-1)-qprime_nx(nx-2) !linear extrapolate
  qprime_nx(1)=two*qprime_nx(2)-qprime_nx(3) !linear extrapolate

  call spline(pfn_nx,qprime_nx,nx,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do j=1,nflux
     call splint(pfn_nx,qprime_nx,tmp_y2,nx,pfn(j),y_tmp) 
     qprime(j)=y_tmp !to get qprime array corresponding to the new radial grid points
  enddo

!!$  do j=1,nflux
!!$     qprime(j)=(q(psival_new(j)+dpsi_nx*0.2_p_)-q(psival_new(j)-dpsi_nx*0.2_p_))/(0.4_p_*dpsi_nx)
!!$  enddo
!!$
!!$  do j=1,nflux
!!$     write(*,*) j,psival_new(j),qprime(j)
!!$  enddo
!!$  stop

  call spline(pfn_nx,tfn_nx,nx,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do j=1,nflux
     call splint(pfn_nx,tfn_nx,tmp_y2,nx,pfn(j),tfn(j)) 
  enddo

!!$contains
!!$
!!$  function q(psival) 
!!$    implicit none
!!$    real(p_):: q,psival
!!$    call splint(psival_nx,qpsi_nx,q_y2,nx,psival,y_tmp) 
!!$    q=y_tmp 
!!$  end function q

end subroutine interpolate_flux_function



 
