subroutine contour(x_lcfs,z_lcfs,np_lcfs,x_axis,z_axis,psival,x_contour,z_contour)
!given a value of the poloidal flux, psival, this subroutine find the magnetic surface corresponding to this poloidal flux
!these codes are drawn from subroutine calculate_contours() in the file calculate_contours.f90
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  implicit none

  integer,intent(in):: np_lcfs
  real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs),x_axis,z_axis,psival
  real(p_),intent(out):: x_contour(np_lcfs),z_contour(np_lcfs)

  real(p_),parameter:: xacc=1.0d-6 !tolerance used in bi-section root-finder
  real(p_):: x1,x2,z1,z2

  real(p_):: slope(np_lcfs),slope2(np_lcfs)
  real(p_):: rtbis !function name of the root finder using the bisection method
  real(p_):: zfunc,xfunc !equation of the straight line (in poloidal plane) that passing throught the magnetic axis point and one point on LCFS
  real(p_):: one_dim_psi_func,one_dim_psi_func2 !one dimension function [psi(x,z(x)) and psi(x(z),z)]on the straight line mentioned in the above.
  external:: one_dim_psi_func,one_dim_psi_func2 !this two function will be passed to a root-finding subroutine
  integer:: i


  !do i=1,np_lcfs-1
  do i=1,np_lcfs
     slope(i)= (z_lcfs(i)-z_axis)/(x_lcfs(i)-x_axis) !the slope for function Z=Z(X)
     slope2(i)=(x_lcfs(i)-x_axis)/(z_lcfs(i)-z_axis) !the slope for function X=X(Z)
     !write(*,*) i,slope(i),slope2(i)
  enddo

!!$  write(*,*) maxval(slope),minval(slope)
!!$  write(*,*) maxloc(slope),minloc(slope)
!!$  write(*,*) x_axis,x_lcfs( maxloc(slope)),x_lcfs(minloc(slope))
  do i=1,np_lcfs-1  !exclude i=np_lcfs because it is identical to i=1
     if(abs(slope(i)).le.1.0_p_) then !use Z=Z(X) function, the reason that I switch between using function X=X(Z) and Z=Z(X) is to aviod large slope.
        x1=x_axis
        x2=x_lcfs(i) !+0.01 !shift left a little to gurrantee that the range is enough for a root to lie in
        x_contour(i)=rtbis(one_dim_psi_func,x1,x2,xacc,x_axis,z_axis,slope(i),psival)
        z_contour(i)=zfunc(x_axis,z_axis,slope(i),x_contour(i))
     else !switch to using X=X(Z) function
        z1=z_axis
        z2=z_lcfs(i)
        z_contour(i)=rtbis(one_dim_psi_func2,z1,z2,xacc,x_axis,z_axis,slope2(i),psival)
        x_contour(i)=xfunc(x_axis,z_axis,slope2(i),z_contour(i)) 
     endif
  enddo

  x_contour(np_lcfs)=x_contour(1) !i=1 and i=np_lcfs respectively corresponds to theta=0 and theta=2pi, so they are equal
  z_contour(np_lcfs)=z_contour(1) !i=1 and i=np_lcfs respectively corresponds to theta=0 and theta=2pi, so they are equal


end subroutine contour

