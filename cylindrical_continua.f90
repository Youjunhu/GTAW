  subroutine cylindrical_continua(mpoloidal,nflux,psival_normalized,toroidal_flux,qpsi,rho,pressure,jacobian,bsq,dpsidra) 
  !obtain the corresponding cylindrical continuous spectrum, which will be compared with continua in toroidal geometry
  use precision,only:p_
  use constants,only:one,twopi,mu0,gamma !gamma is polytrope index
  use toroidal_harmonics,only:nh !toroidal mode number
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: qpsi(nflux),rho(nflux),pressure(nflux),psival_normalized(nflux),toroidal_flux(nflux)
  real(p_),intent(in):: jacobian(mpoloidal,nflux),bsq(mpoloidal,nflux),dpsidra(nflux)
  integer:: i,j,mh
  real(p_):: sum,jacobian_av(nflux),bsq_av(nflux)
  !--calculate the Alfven continua in cylindrical geometry limit
  !first calculate the average value of Jacobian on every magnetic surface


  do j=1,nflux
     sum=0
     do i=1,mpoloidal
        sum=sum+jacobian(i,j)
     enddo
     jacobian_av(j)=sum/mpoloidal
  enddo
  !then use the formual given in my note to calculate the cylindrical alfven continua for different poloidal mode numbers
  open(88,file='cylindrical_alfven.txt')
!  do j=2,nflux-1
  do j=1,nflux
     write(88,*) sqrt(psival_normalized(j)), sqrt(toroidal_flux(j)),&
          & (sqrt((dpsidra(j)/jacobian_av(j))**2*(mh-nh*qpsi(j))**2/(mu0*rho(j)))/twopi/1000._p_,mh=-10,10)
  enddo
  close(88)
!!$do j=1,nflux
!!$write(*,*) j, dpsidra(j), 'cylindrical'
!!$enddo
!the following calculates the sound continua (the formula is given in my notes)
!calculate the the average value of bsq on a flux surface
  do j=1,nflux
     sum=0.
     do i=1,mpoloidal
        sum=sum+bsq(i,j)
     enddo
     bsq_av(j)=sum/mpoloidal
  enddo

  open(88,file='cylindrical_sound.txt')
!  do j=2,nflux-1
  do j=1,nflux
     write(88,*) sqrt(psival_normalized(j)), sqrt(toroidal_flux(j)),&
          & (sqrt(gamma*pressure(j)/rho(j)/(gamma*pressure(j)/rho(j)+bsq_av(j)/mu0/rho(j)))* &
          & sqrt((dpsidra(j)/jacobian_av(j))**2*(mh-nh*qpsi(j))**2/(mu0*rho(j)))/twopi/1000._p_, mh=-1,20)
  enddo
  close(88)
!write(*,*) 'sound speed at magentic axis (10^6 m/s):', sqrt(gamma*pressure(1)/rho(1))/10**6
end subroutine cylindrical_continua
