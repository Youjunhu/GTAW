subroutine artifial_mode(r_new,z_new,ra,nflux,mpoloidal) 
  !this is to construct artifial mode structure by using analytical expression, for test purpose
  use precision,only:p_
  use constants,only:zero,one,two,twopi,pi
  implicit none
  integer,intent(in):: nflux,mpoloidal
  real(p_),intent(in):: ra(nflux),r_new(nflux,mpoloidal),z_new(nflux,mpoloidal)
  real(p_):: theta(mpoloidal)
  !integer,parameter:: m=7
  integer,parameter:: m=3
  integer:: i,j,k
  real(p_),parameter:: r0=0.3_p_,scale=0.1_p_
  real(p_):: phase(nflux)
  real(p_),parameter:: omega=twopi/one
  real(p_):: time
  integer,parameter:: nt=3 !number of time step
  character(100):: filename='zxxx.txt'

  do i=1,mpoloidal
     theta(i)=zero+twopi/(mpoloidal-1)*(i-1)
  enddo

  do j=1,nflux
     !phase(j)=(ra(j)-r0)*twopi*4
     !phase(j)=-(ra(j)-r0)*twopi*4
     phase(j)=zero
  enddo

  !open(11,file='artifial_mode')

  do k=1,nt !for animation
     time=zero+one/(nt-1)*(k-1)
     write(filename(2:4),'(i3.3)') k
     open(11,file=filename)
     do j=1,nflux
        do i=1,mpoloidal
           !write(11,*) r_new(i,j),z_new(i,j), cos(phase(j)-m*theta(i)-omega*time)*exp(-(r0-ra(j))**2/scale**2)
           if (theta(i).le.pi) then
              write(11,*) theta(i), cos(phase(j)-m*theta(i)-omega*time)*exp(-(r0-ra(j))**2/scale**2)+&
                   &  cos(phase(j)-(m+1)*theta(i)-omega*time)*exp(-(r0-ra(j))**2/scale**2), &
                   & abs(2.*cos(((m+1)*theta(i)-m*theta(i))/two))*exp(-(r0-ra(j))**2/scale**2) !amplitude
           else
              write(11,*) theta(i)-twopi, cos(phase(j)-m*theta(i)-omega*time)*exp(-(r0-ra(j))**2/scale**2)+&
                   &  cos(phase(j)-(m+1)*theta(i)-omega*time)*exp(-(r0-ra(j))**2/scale**2), &
                   & abs(2.*cos(((m+1)*theta(i)-m*theta(i))/two))*exp(-(r0-ra(j))**2/scale**2) !amplitude
           endif
        enddo
        write(11,*) 
        write(11,*)
     enddo
     close(11)
  enddo
end subroutine artifial_mode
