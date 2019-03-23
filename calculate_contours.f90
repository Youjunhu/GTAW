subroutine calculate_contours_wrapper(recalculate_flux_coordinates,poloidal_angle_type,mpoloidal,nflux,psival,fpsi,&
     & r_axis,z_axis,x_lcfs,z_lcfs,np_lcfs,r_new,z_new,av_one_over_rsq,dvdpsi,safety_factor,circumference,&
     & r_grad_psi_s,bp_sq_av,b_av,one_over_b_av)
  !wrapper of "calculate_contours" subroutine
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  implicit none
  logical,intent(in):: recalculate_flux_coordinates
  character(100),intent(in):: poloidal_angle_type
  integer,intent(in):: mpoloidal,nflux,np_lcfs
  real(p_),intent(in):: psival(nflux),fpsi(nflux)
  real(p_),intent(in):: r_axis,z_axis,x_lcfs(np_lcfs),z_lcfs(np_lcfs) 
  real(p_),intent(out):: r_new(mpoloidal,nflux),z_new(mpoloidal,nflux) !points on every magnetic surface are evenly distributed about poloidal angle
  real(p_),intent(out):: av_one_over_rsq(nflux),dvdpsi(nflux),bp_sq_av(nflux),b_av(nflux),one_over_b_av(nflux)
  real(p_),intent(out):: safety_factor(nflux),circumference(nflux),r_grad_psi_s(nflux) !circumference of every contour

  integer:: i,j

  if(recalculate_flux_coordinates.eqv..false.) then   !read the saved data 
     open(52,file='flux_coordinates_data.txt')
     do j=1,nflux
        read(52,*) (r_new(i,j),i=1,mpoloidal),(z_new(i,j),i=1,mpoloidal)
        read(52,*) av_one_over_rsq(j),dvdpsi(j),bp_sq_av(j),b_av(j),one_over_b_av(j)
        read(52,*) safety_factor(j),circumference(j) !circumference of every contour
     enddo
     close(52)
  else !claculate the flux coordinates
     call calculate_contours(poloidal_angle_type,mpoloidal,nflux,psival,fpsi,r_axis,z_axis,x_lcfs,z_lcfs,np_lcfs,&
          & r_new,z_new,av_one_over_rsq,dvdpsi,safety_factor,circumference,r_grad_psi_s,bp_sq_av,b_av,one_over_b_av)
  endif
end subroutine calculate_contours_wrapper


subroutine calculate_contours(poloidal_angle_type,mpoloidal,nflux,psival,fpsi,x_axis,z_axis,x_lcfs,z_lcfs,np_lcfs,&
     & r_new,z_new,av_one_over_rsq,dvdpsi,safety_factor,circumference,r_grad_psi_s,bp_sq_av,b_av,one_over_b_av)
  !This routine searches for a series of magnetic surfaces (within LCFS) corresponding to the values of psi specified in psival array
  !this is a key step in constructing a magnetic surface coordinate system
  !Then the code calculates the theta coordinates of points on every magnetic surface
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  implicit none
  integer,intent(in):: mpoloidal,nflux,np_lcfs
  real(p_),intent(in):: psival(nflux),fpsi(nflux)
  real(p_),intent(in):: x_axis,z_axis
  real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs) 
  character(100),intent(in):: poloidal_angle_type
  real(p_),intent(out):: r_new(mpoloidal,nflux),z_new(mpoloidal,nflux)
  real(p_),intent(out):: av_one_over_rsq(nflux),dvdpsi(nflux),bp_sq_av(nflux),b_av(nflux),one_over_b_av(nflux)
  real(p_),intent(out):: safety_factor(nflux),circumference(nflux),r_grad_psi_s(nflux) !circumference of every contour
  !  real(p_):: grad_psi(mpoloidal,nflux)
  real(p_):: grad_psi(np_lcfs,nflux)
  real(p_):: b_total,one_over_b
  real(p_):: theta_uniform(mpoloidal)
  real(p_):: kernel(np_lcfs,nflux)

  real(p_),parameter:: xacc=1.0d-6 !tolerance used in bi-section root-finder
  real(p_):: x1,x2,z1,z2
  real(p_):: x_contour(np_lcfs,nflux),z_contour(np_lcfs,nflux) !magnetic surface found
  real(p_):: theta(np_lcfs,nflux) !values of theta on the grids
  real(p_):: theta_old(np_lcfs),y2(np_lcfs),r_old(np_lcfs),z_old(np_lcfs),y_tmp !used in interpolation for theta on one flux surface
  real(p_):: y_new(np_lcfs) !used in linear interpolation
  real(p_):: slope(np_lcfs),slope2(np_lcfs)
  real(p_):: dl(np_lcfs-1,nflux) !arc lengh between neighbour points on every magnetic surface
  real(p_):: x_tmp,sum_tmp(np_lcfs),sum,sum1,sum2,tmp
  real(p_):: rtbis !function name of the root finder using the bisection method
  real(p_):: zfunc,xfunc !equation of the straight line (in poloidal plane) that passing throught the magnetic axis point and one point on LCFS
  real(p_):: one_dim_psi_func,one_dim_psi_func2 !one dimension function [psi(x,z(x)) and psi(x(z),z)]on the straight line mentioned in the above.
  external:: one_dim_psi_func,one_dim_psi_func2 !this two function will be passed to a root-finding subroutine
  real(p_):: psi_gradient_func !2d interpolating function psi_gradient(x,z)
  real(p_):: bp(np_lcfs,nflux) !poloidal magnetic field
  real(p_):: direction !flag of the direction of theta(i,j) with i increasing, direction>0 for count clockwise, direction<0 for clowckwise
  integer:: i,j



  !Since the points on the LCFS (last closed flux surfac) are given in the G-eqdsk-file,
  !I use these points and the magnetic axis to define a series of straight lines on poloidal plane,
  !which starts from the magnetic axis and ends at one of the points on the LCFS.
  !the number of points on the LCFS is specified by np_lcfs, which is usually a large number (e.g. 349)

!!$  write(*,*) 'starting point of lcfs', x_lcfs(1),z_lcfs(1)
!!$  write(*,*) 'ending point of lcfs', x_lcfs(np_lcfs),z_lcfs(np_lcfs)
!!$  write(*,*) 'increasing', x_lcfs(10),z_lcfs(10)
!!$  stop
  !  write(*,*) 'discrete psi value= ', psival
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
     !do j=1,nflux-1 !exclude LCFS since it is explicitly recorded in g-file
     !$omp parallel do     
     do j=1,nflux
        if(abs(slope(i)).le.1.0_p_) then !use Z=Z(X) function, the reason that I switch between using function X=X(Z) and Z=Z(X) is to aviod large slope.
           x1=x_axis
           x2=x_lcfs(i) !+0.01 !shift left a little to gurrantee that the range is enough for a root to lie in
           x_contour(i,j)=rtbis(one_dim_psi_func,x1,x2,xacc,x_axis,z_axis,slope(i),psival(j))
           z_contour(i,j)=zfunc(x_axis,z_axis,slope(i),x_contour(i,j))
        else !switch to using X=X(Z) function
           z1=z_axis
           z2=z_lcfs(i)
           z_contour(i,j)=rtbis(one_dim_psi_func2,z1,z2,xacc,x_axis,z_axis,slope2(i),psival(j))
           x_contour(i,j)=xfunc(x_axis,z_axis,slope2(i),z_contour(i,j)) 
        endif
     enddo
     !$omp end parallel do
     ! write(*,*) 'finish in one poloidal direction', ' i=',i
  enddo

  do j=1,nflux
     x_contour(np_lcfs,j)=x_contour(1,j) !i=1 and i=np_lcfs respectively corresponds to theta=0 and theta=2pi, so they are equal
     z_contour(np_lcfs,j)=z_contour(1,j) !i=1 and i=np_lcfs respectively corresponds to theta=0 and theta=2pi, so they are equal
  enddo

!!$  do i=1,np_lcfs !include the contour corresponding to the LCFS
!!$     x_contour(i,nflux)=x_lcfs(i)
!!$     z_contour(i,nflux)=z_lcfs(i)
!!$  enddo

  write(*,*) '>>>>Finish determining magnetic surfaces.'
  call wrt_r_z(np_lcfs,nflux,x_contour,z_contour,'surface1.txt','radial1.txt')

  call arc_length(x_contour,z_contour,nflux,np_lcfs,dl)
  call plot_poloidal(np_lcfs-1,nflux,dl,dl,'dl.txt')



  do i=1,np_lcfs
     do j=1,nflux
        grad_psi(i,j)= psi_gradient_func(x_contour(i,j),z_contour(i,j))
        bp(i,j)=grad_psi(i,j)/x_contour(i,j) !polodal magnetic field
     enddo
  enddo


  do j=1,nflux !calculate the surface averaged quantity: <1/R^2>
     sum1=0.
     sum2=0.
     do i=1,np_lcfs-1
        sum1=sum1+one/((x_contour(i,j)+x_contour(i+1,j))/two)**2/((bp(i,j)+bp(i+1,j))/two)*dl(i,j)
        sum2=sum2+dl(i,j)/((bp(i,j)+bp(i+1,j))/two)
     enddo
     av_one_over_rsq(j)=sum1/sum2
  enddo

  do j=1,nflux !calculate the surface averaged quantity: <Bp^2>
     sum1=0.
     sum2=0.
     do i=1,np_lcfs-1
        sum1=sum1+(bp(i,j)+bp(i+1,j))/two*dl(i,j)
        sum2=sum2+dl(i,j)/((bp(i,j)+bp(i+1,j))/two)
     enddo
     bp_sq_av(j)=sum1/sum2
  enddo


  do j=1,nflux !calculate the surface averaged quantity: <B>
     sum1=0.
     sum2=0.
     do i=1,np_lcfs-1
        b_total=sqrt(((bp(i,j)+bp(i+1,j))/two)**2+((fpsi(j)/x_contour(i,j)+fpsi(j)/x_contour(i+1,j))/two)**2)
        sum1=sum1+b_total*dl(i,j)/((bp(i,j)+bp(i+1,j))/two)
        sum2=sum2+dl(i,j)/((bp(i,j)+bp(i+1,j))/two)
     enddo
     b_av(j)=sum1/sum2
  enddo


  do j=1,nflux !calculate the surface averaged quantity: <1/B>
     sum1=0.
     sum2=0.
     do i=1,np_lcfs-1
        one_over_b=one/sqrt(((bp(i,j)+bp(i+1,j))/two)**2+((fpsi(j)/x_contour(i,j)+fpsi(j)/x_contour(i+1,j))/two)**2)
        sum1=sum1+one_over_b/((bp(i,j)+bp(i+1,j))/two)*dl(i,j)
        sum2=sum2+dl(i,j)/((bp(i,j)+bp(i+1,j))/two)
     enddo
     one_over_b_av(j)=sum1/sum2
  enddo

  !check wheter the direction of the sequecne (r(i),z(i)) with i increasing is clockwise or anticlockwise when viewed along grad_phi direction
  !This is achieved by using the determination of the direction matrix (a well known method in graphic theory).
  !Because the contours of Psi considered here are always convex polygons (instead of concave polygons), we can select any vertex on the curve to calculate the direction matrix. (refer to wikipedia about the direction matrix)
  direction=(x_lcfs(2)-x_lcfs(1))*(z_lcfs(3)-z_lcfs(1))-(x_lcfs(3)-x_lcfs(1))*(z_lcfs(2)-z_lcfs(1))
  if(direction .lt. 0.) then
     write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) with i increasing is clockwise'
  else if (direction .gt. 0.) then
     write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) with i increasing is anticlockwise'
  else
     stop 'the three vertex (points) used in calculating the direction matrix is collinear'
  endif

  do j=1,nflux
     sum=0.
     do i=1,np_lcfs-1
        sum=sum+dl(i,j)
     enddo
     circumference(j)=sum
  enddo
  open(123,file='circumference.txt')
  do j=1,nflux
     write(123,*) j, circumference(j)
  enddo
  close(123)
  write(*,*) 'the radius of the first flux surface (adjacent to the magnetic axis) is (m) ',circumference(1)/twopi

  do j=1,nflux
     sum=0.
     do i=1,np_lcfs-1
        tmp=(x_contour(i,j)*grad_psi(i,j)+x_contour(i+1,j)*grad_psi(i+1,j))/two
        sum=sum+dl(i,j)/tmp
     enddo
     r_grad_psi_s(j)=sum !r_grad_psi_s is used in constructing straight-line poloidal angle
  enddo

  !calculate the theta coordinates of (i,j) points
  do j=1,nflux !I want the positive direction of theta to be in the anticlockwise direction when viewed in the grad_phi direction
     if (direction .gt.0.) then !indicates count-clockwise, then the theta value of the starting point is set to be zero
        theta(1,j)=0.0 
     else if(direction .lt.0) then !indicates clockwise, then the theta value of the starting point is set to be twopi.   
        theta(1,j)=twopi 
     endif

     if(trim(poloidal_angle_type).eq.'equal-arc') then
        do i=2,np_lcfs
           theta(i,j)=theta(i-1,j)+sign(one,direction)*dl(i-1,j)/circumference(j)*twopi !equal arc length theta coordinates
        enddo
     else if(trim(poloidal_angle_type).eq.'straight-line') then 
        do i=2,np_lcfs
           tmp=(x_contour(i-1,j)*grad_psi(i-1,j)+x_contour(i,j)*grad_psi(i,j))/two
           theta(i,j)=theta(i-1,j)+ sign(one,direction)*dl(i-1,j)/tmp*twopi/r_grad_psi_s(j)
        enddo
     else
        stop "please choose poloidal_angle_type between 'equal-arc' and 'straight-line'"
     endif
     !if the sequence (x_contour(i,j),z_contour(i,j)) with i increasing is anticlockwise, then theta is increasing when i increases
     !if the sequence (x_contour(i,j),z_contour(i,j)) with i increasing is clockwise, then theta is decreasing when i increases
     !the above treatment ensures that the positive direction of poloidal angle is always in the anticlockwise direction
  enddo

  call plot_poloidal(np_lcfs,nflux,theta,theta,'theta.txt')

  !The theta array obtained above is usually not uniform. Next, interpolate R and Z to uniform theta grids on every magnetic surfac.
  do i=1,mpoloidal !First, constructe a uniform theta grids ranging from 0 to twopi
     theta_uniform(i)=0.0_p_+twopi/(mpoloidal-1)*(i-1)
  enddo

  do j=1,nflux !then, for every magnetic surface, interpolate to get value of R and Z on uniform theta grid points
     if (direction .gt.0.) then !the direction is anti-clockwise, which is what I want.
        do i=1,np_lcfs
           theta_old(i)=theta(i,j)
           r_old(i)=x_contour(i,j)
           z_old(i)=z_contour(i,j)
        enddo
     else if(direction .lt.0) then 
        !the spline interpolating works only for sequence x1<x2<x3...<xn.  In the clockwise case, I need to revert the array obtained in the above because their order is twopi=x1>x2>x3..>xn=0
        do i=1,np_lcfs
           theta_old(i)=theta(np_lcfs-i+1,j)
           r_old(i)=x_contour(np_lcfs-i+1,j)
           z_old(i)=z_contour(np_lcfs-i+1,j)
        enddo
     endif

     call spline(theta_old,r_old,np_lcfs,2.d30,2.d30,y2) !prepare the second order derivative needed in the cubic spline interpolation
     do i=2,mpoloidal-1
        call splint(theta_old,r_old,y2,np_lcfs,theta_uniform(i),y_tmp) !to get R corresponding uniform theta grids
        r_new(i,j)=y_tmp
     enddo

     call spline(theta_old,z_old,np_lcfs,2.d30,2.d30,y2) !prepare the second order derivative needed in the cubic spline interpolation
     do i=2,mpoloidal-1
        call splint(theta_old,z_old,y2,np_lcfs,theta_uniform(i),y_tmp) !to get Z corresponding uniform theta grids
        z_new(i,j)=y_tmp
     enddo

     !the following is to make sure that the points theta=0 and theta=twopi have identical location (since these two points are the end points for the above interpolation, the results obtained above may be inaccurate, I find this from the graphic of z-theta, so I added the following steps)
     r_new(1,j)=x_contour(1,j)
     z_new(1,j)=z_contour(1,j)
     r_new(mpoloidal,j)=x_contour(1,j)
     z_new(mpoloidal,j)=z_contour(1,j)
  enddo

  write(*,*) '>>>Finish constructing the polodial angle coordinate.'
  call wrt_r_z(mpoloidal,nflux,r_new,z_new,'surface2.txt','radial2.txt')


  !calculate dvdpsi by using Eq. (5.39) in S. Jardin's book. The result can be benchmarked with the result obtained by using the equation that involves the safety factor (this benchmark is done in the main subroutine)
  do j=1,nflux
     do i=1,np_lcfs
        kernel(i,j)=x_contour(i,j)/psi_gradient_func(x_contour(i,j),z_contour(i,j))
     enddo
  enddo

  do j=1,nflux
     sum=0.0
     do i=1,np_lcfs-1
        sum=sum+(kernel(i,j)+kernel(i+1,j))/two*dl(i,j)
     enddo
     dvdpsi(j)=sum*twopi
  enddo
  !--------------end------------------

  !calculate safety factor, which will be compared with the value given in G-file (the comparison is done in the main program)
  do j=2,nflux
     !  do j=1,nflux 
     sum=0.
     do i=1,np_lcfs-1
        sum=sum+dl(i,j)/(x_contour(i,j)*psi_gradient_func(x_contour(i,j),z_contour(i,j)) &
             & +x_contour(i+1,j)*psi_gradient_func(x_contour(i+1,j),z_contour(i+1,j)))*two !use trappeziod formula
     enddo
     safety_factor(j)=-one/twopi*sum*fpsi(j) !the term sign(J) will be included in the main program
     !     write(*,*) j,sum,fpsi(j),safety_factor(j)
  enddo
  safety_factor(1)=two*safety_factor(2)-safety_factor(3) !use linear expolation to obtain q at the innermost flux surface
  !--------------end-----------


  !save the data for later use so that we do not need to recalculate this data again when recalculate_flux_coordinates=.false.
  open(52,file='flux_coordinates_data.txt')
  do j=1,nflux
     write(52,*) (r_new(i,j),i=1,mpoloidal),(z_new(i,j),i=1,mpoloidal)
     write(52,*) av_one_over_rsq(j),dvdpsi(j),bp_sq_av(j),b_av(j),one_over_b_av(j)
     write(52,*) safety_factor(j),circumference(j) !circumference of every contour
  enddo
  close(52)

 
end subroutine calculate_contours


function one_dim_psi_func(x_axis,z_axis,slope,psival,x) 
  !poloidal flux as a function of x on a straight line with slope "slope" in poloidal plane
  use precision,only:p_
  implicit none
  real(p_):: one_dim_psi_func,x,x_axis,z_axis,slope,psival
  real(p_):: zfunc,psi_func
  one_dim_psi_func=psi_func(x,zfunc(x_axis,z_axis,slope,x))-psival
end function one_dim_psi_func


function zfunc(x_axis,z_axis,slope,x) !straight line Z=Z(x) with slope "slope" in poloidal plane starting from the location of magnetic axis
  use precision,only:p_
  implicit none
  real(p_):: zfunc,x,x_axis,z_axis,slope
  zfunc=z_axis+slope*(x-x_axis)
end function zfunc


function one_dim_psi_func2(x_axis,z_axis,slope,psival,z) result(fun_val)
  !poloidal flux as a function of z on a straight line with slope "slope" in poloidal plane
  use precision,only:p_
  implicit none
  real(p_):: fun_val,z
  real(p_):: x_axis,z_axis,slope,psival
  real(p_):: xfunc,psi_func
  fun_val=psi_func(xfunc(x_axis,z_axis,slope,z),z)-psival
end function one_dim_psi_func2


function xfunc(x_axis,z_axis,slope,z) !straight line X=X(Z) with slope "slope" in poloidal plane starting from the location of magnetic axis
  use precision,only:p_
  implicit none
  real(p_):: xfunc,z
  real(p_)::x_axis,z_axis,slope
  xfunc=x_axis+slope*(z-z_axis)
end function xfunc




