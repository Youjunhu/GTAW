subroutine plasma_volume(nflux,mpoloidal,jacobian,dpsi,dtheta,dvdpsi,psival_new,vol)
  !this is to calculate the volume within the LCFS
  use precision,only: p_
  use constants,only: two,twopi
  implicit none
  integer,intent(in)::nflux,mpoloidal
  real(p_),intent(in):: jacobian(mpoloidal,nflux),dpsi,dtheta,dvdpsi(nflux),psival_new(nflux)
  real(p_),intent(out):: vol
  real(p_)::dv,vol2
  integer:: i,j

  vol=0._p_
  do i=1,mpoloidal
     do j=1,nflux 
        dv=-jacobian(i,j)*dpsi*dtheta
        vol=vol+dv*twopi
     enddo
  enddo

  vol2=0._p_
  do j=1,nflux-1
     vol2=vol2+(dvdpsi(j)+dvdpsi(j+1))/two*(psival_new(j+1)-psival_new(j))
  enddo
  write(*,*) 'Plasma volume within LCFS is (m^3)',vol, 'vol2=', vol2

end subroutine plasma_volume



subroutine calculate_volume_averaged_beta(nflux,mpoloidal,pressure_normalized,jacobian,dpsi,dtheta,averaged_pressure)
  use precision,only: p_
  use constants,only: twopi
  implicit none
  integer,intent(in)::nflux,mpoloidal
  real(p_),intent(in):: pressure_normalized(nflux),jacobian(mpoloidal,nflux),dpsi,dtheta
  real(p_),intent(out):: averaged_pressure
  real(p_)::dv,vol
  integer:: i,j

  vol=0._p_
  do i=1,mpoloidal
     do j=1,nflux 
        dv=jacobian(i,j)*dpsi*dtheta
        vol=vol+dv*twopi
     enddo
  enddo

  averaged_pressure=0._p_
  do i=1,mpoloidal
     do j=1,nflux 
        dv=jacobian(i,j)*dpsi*dtheta
        averaged_pressure=averaged_pressure+pressure_normalized(j)*dv*twopi
     enddo
  enddo

  averaged_pressure=averaged_pressure/vol
  write(*,*) 'Volume averaged beta',averaged_pressure

end subroutine calculate_volume_averaged_beta



subroutine calculate_surface_averaged_beta(nflux,mpoloidal,pressure_normalized,jacobian,r,dpsi,dtheta)
  use precision,only: p_
  use constants,only: twopi
  implicit none
  integer,intent(in)::nflux,mpoloidal
  real(p_),intent(in):: pressure_normalized(nflux),jacobian(mpoloidal,nflux),dpsi,dtheta,r(mpoloidal,nflux)
  real(p_)::averaged_pressure,ds,area
  integer:: i,j

  area=0._p_
  do i=1,mpoloidal
     do j=1,nflux 
        ds=jacobian(i,j)*dpsi*dtheta/r(i,j)
        area=area+ds
     enddo
  enddo

  averaged_pressure=0._p_
  do i=1,mpoloidal
     do j=1,nflux 
        ds=jacobian(i,j)*dpsi*dtheta/r(i,j)
        averaged_pressure=averaged_pressure+pressure_normalized(j)*ds
     enddo
  enddo

  write(*,*) 'Surface averaged beta',averaged_pressure/area

end subroutine calculate_surface_averaged_beta



subroutine toroidal_current(mpoloidal,nflux,dpsi,dtheta,jacobian,ffprime,pprime,r,current)
  !calculate the total toroidal current within LCFS using the data of pressure and toroidal field function g
  use precision,only:p_
  use constants,only: mu0
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: dpsi,dtheta
  real(p_),intent(in):: ffprime(nflux),pprime(nflux),jacobian(mpoloidal,nflux),r(mpoloidal,nflux)
  real(p_),intent(out):: current
  real(p_):: it,ds,sum1,sum2
  integer:: i,j

  it=0._p_
  sum1=0._p_
  sum2=0._p_
  do i=1,mpoloidal
     do j=1,nflux 
        !ds=dpsi*dtheta*abs(jacobian(i,j))/r(i,j) !dpsi*dtheta*abs(jacobian(i,j))/r(i,j) is the surface element
         ds=abs(dpsi*dtheta*jacobian(i,j))/r(i,j) !abs(dpsi*dtheta*jacobian(i,j))/r(i,j) is the surface element
        it=it+(mu0*r(i,j)*pprime(j)+ffprime(j)/r(i,j))*ds
        sum1=sum1+mu0*r(i,j)*pprime(j)*ds
        sum2=sum2+ffprime(j)/r(i,j)*ds
     enddo
  enddo
  current=it/mu0
  write(*,*) 'toroidal current within LCFS is (kA) (calculated by using the numerical pressure and toroidal field function)', &
       & current/1000._p_, 'pressure part=',sum1/mu0/1000._p_, 'toroidal field part=',sum2/mu0/1000._p_
!  write(*,*) 'jacobian at (10,10)', jacobian(10,10)
end subroutine toroidal_current



subroutine toroidal_current2(mpoloidal,nflux,dpsi,dtheta,jacobian,psi_r,psi_rr,psi_zz,r,current)
  !calculate the total toroidal current within LCFS, using the magnetic field, instead of the pressure
  use precision,only:p_
  use constants,only: mu0,one
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: dpsi,dtheta
  real(p_),intent(in):: psi_r(mpoloidal,nflux),psi_rr(mpoloidal,nflux),psi_zz(mpoloidal,nflux)
  real(p_),intent(in)::jacobian(mpoloidal,nflux),r(mpoloidal,nflux)
  real(p_),intent(out):: current
  real(p_):: it,ds,laplace_star
  integer:: i,j

  it=0._p_
  do i=1,mpoloidal
     do j=1,nflux 
        !ds=dpsi*dtheta*abs(jacobian(i,j))/r(i,j) !dpsi*dtheta*abs(jacobian(i,j))/r(i,j) is the surface element
        ds=abs(dpsi*dtheta*jacobian(i,j))/r(i,j) !abs(dpsi*dtheta*jacobian(i,j))/r(i,j) is the surface element
        laplace_star=psi_zz(i,j)+r(i,j)*(-one/r(i,j)**2*psi_r(i,j)+one/r(i,j)*psi_rr(i,j))
        it=it+laplace_star/(-r(i,j))*ds
     enddo
  enddo
  current= it/mu0
  write(*,*) 'toroidal current within LCFS is (kA) (calculated by using the numerical poloidal flux)',current/1000._p_
end subroutine toroidal_current2


subroutine force_analysis(mpoloidal,nflux,psi_r,psi_rr,psi_zz,r,z,ffprime,pprime)
  use precision,only:p_
  use constants,only: mu0,one
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: psi_r(mpoloidal,nflux),psi_rr(mpoloidal,nflux),psi_zz(mpoloidal,nflux)
  real(p_),intent(in)::r(mpoloidal,nflux),z(mpoloidal,nflux)
  real(p_),intent(in)::ffprime(nflux),pprime(nflux)
  real(p_):: laplace_star,force1,force2,force3
  integer:: i,j


 open(33,file='force_analysis.txt')
  do j=1,nflux 
     do i=1,mpoloidal
        laplace_star=psi_zz(i,j)+r(i,j)*(-one/r(i,j)**2*psi_r(i,j)+one/r(i,j)*psi_rr(i,j))
        force1=laplace_star
        force2=mu0*r(i,j)**2*pprime(j)
        force3=ffprime(j)
        write(33,*) r(i,j),z(i,j),force1,force2,force3
     enddo
     write(33,*)
  enddo
  close(33)
end subroutine force_analysis



subroutine volume_average(nflux,mpoloidal,pressure,jacobian,dpsi,dtheta)
!this is a general subroutine that calculates the volume averaged equilibrium quantity
  use precision,only: p_
  use constants,only: twopi
  implicit none
  integer,intent(in)::nflux,mpoloidal
  real(p_),intent(in):: pressure(nflux),jacobian(mpoloidal,nflux),dpsi,dtheta
  real(p_):: averaged_pressure
  real(p_)::dv,vol
  integer:: i,j

  vol=0._p_
  do i=1,mpoloidal
     do j=1,nflux 
        dv=jacobian(i,j)*dpsi*dtheta
        vol=vol+dv*twopi
     enddo
  enddo

  averaged_pressure=0._p_
  do i=1,mpoloidal
     do j=1,nflux 
        dv=jacobian(i,j)*dpsi*dtheta
        averaged_pressure=averaged_pressure+pressure(j)*dv*twopi
     enddo
  enddo

  averaged_pressure=averaged_pressure/vol

end subroutine volume_average


subroutine check_directions(current,psi_axis,psi_lcfs,fpsi_axis,fpsi_lcfs)
!check the directions of the poloidal magnetic field and poloidal electric current
  use precision,only: p_
  implicit none

  real(p_),intent(in):: current, psi_axis,psi_lcfs,fpsi_axis,fpsi_lcfs

  if(current<0.) then
     write(*,*) 'poloidal magentic field is anticlockwise viewed in grad(phi) direction'
  else
     write(*,*) 'poloidal magentic field is clockwise viewed in grad(phi) direction'
  endif

  if(psi_axis<psi_lcfs) then
     write(*,*) 'poloidal magentic field is anticlockwise viewed in grad(phi) direction'
  else
     write(*,*) 'poloidal magentic field is clockwise viewed in grad(phi) direction'
  endif

  if(fpsi_axis<fpsi_lcfs) then
     write(*,*) 'poloidal current is anticlockwise viewed in grad(phi) direction'
  else
     write(*,*) 'poloidal current is clockwise viewed in grad(phi) direction'
  endif

end subroutine check_directions




subroutine vacuum_toroidal_magnetic_field(r_axis)
  use precision,only: p_
  use constants,only:twopi,mu0
  implicit none
  real(p_),intent(in):: r_axis
  integer::num_tf_groups,turns
  real(p_):: current_per_turn
  namelist/tf_coils_parameters/num_tf_groups,turns,current_per_turn

  open(11,file='gtaw.in')
  read(11,tf_coils_parameters)
  close(11)
  write(*,*) 'based on the parameters of TF coils: ', "number_tf_coils, turns, current_per_turn (kA)=",&
       & num_tf_groups,turns,current_per_turn
  write(*,*) 'the estimated vacuum toroidal magnetic field (Telsa) at magentic axis is ', &
       & mu0*num_tf_groups*turns*current_per_turn/(twopi*r_axis)

end subroutine vacuum_toroidal_magnetic_field
