subroutine calculate_phase(ypath,phase,phase_total,phase_midplane,xir_midplane)
  !calculate the phase of the radial displacement across the minor radius
  use precision,only: p_
  use constants,only:ii,twopi
  use flux_grids,only: nflux,mpoloidal,delta_q,midplane_theta
  use radial_module,only: starting_surface_number,ending_surface_number !radial grid points
  use poloidal_harmonics,only: mh_low,mhtot
  use toroidal_harmonics,only: nh
  implicit none
  complex(p_),intent(in):: ypath(2*mhtot+1,nflux)
  real(p_),intent(out):: phase(mhtot,nflux),phase_total(mpoloidal,nflux),phase_midplane(nflux)
  complex(p_),intent(out):: xir_midplane(mhtot,nflux)
  complex(p_):: sum0(mpoloidal),sum1(mpoloidal),sum3,sum4,tmp
  real(p_):: theta(mpoloidal),tmp_yj
  integer:: i,j,k,mhk
  real(p_):: delta_q_midplane(nflux)
  real(p_):: y(mpoloidal),tmp_y2(mpoloidal)
  real(p_)::  y_tmp

  do j=starting_surface_number,ending_surface_number
     do i=1,mhtot !for radial plasma displacement
        tmp=ypath(i+mhtot,j)/ypath(i+mhtot,starting_surface_number)
        phase(i,j)=acos(real(tmp)/abs(tmp))
        if(imag(tmp).lt.0.) phase(i,j)=-phase(i,j)
     enddo
  enddo

  do i=1,mpoloidal
     theta(i)=0.+(i-1)*twopi/(mpoloidal-1)
  enddo

  !first calculate the values of the radial displacement at the inner most flux surface for all the poloidal directions
  !the values will be used as normalizing factors to let the phase of the radial displacement on the innermost flux surface are all zeros.
  do i=1,mpoloidal
     sum0(i)=0.
     do k=1,mhtot
        mhk=mh_low+(k-1)      
        sum0(i)=sum0(i)+ypath(k+mhtot,starting_surface_number)*exp(ii*mhk*theta(i)) &
             & *exp(ii*nh*delta_q(i,starting_surface_number)) !include the toroidal angle shift factor
     enddo
  enddo

  do j=starting_surface_number,ending_surface_number
     do i=1,mpoloidal
        sum1(i)=0
        do k=1,mhtot
           mhk=mh_low+(k-1) 
           sum1(i)=sum1(i)+ypath(k+mhtot,j)*exp(ii*mhk*theta(i))*exp(ii*nh*delta_q(i,j))
        enddo
        sum1(i)=sum1(i)/sum0(i)
        phase_total(i,j)=acos(real(sum1(i))/abs(sum1(i)))
        if(imag(sum1(i)).lt.0.) phase_total(i,j)=-phase_total(i,j)
     enddo
  enddo


  do j=starting_surface_number,ending_surface_number
     do i=1,mpoloidal
        y(i)=delta_q(i,j)
     enddo
     call spline(theta,y,mpoloidal,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
     call splint(theta,y,tmp_y2,mpoloidal,midplane_theta(j),y_tmp)
     delta_q_midplane(j)=y_tmp

  enddo


  sum3=0.
  j=starting_surface_number+1 ! the reference point to define the phase change
  do k=1,mhtot
     mhk=mh_low+(k-1)      
     sum3=sum3+ypath(k+mhtot,j)*exp(ii*mhk*midplane_theta(j)) &
          & *exp(ii*nh*delta_q_midplane(j)) !include the toroidal angle shift factor
  enddo

  tmp_yj=acos(real(sum3)/abs(sum3))
  if(imag(sum3).lt.0.) tmp_yj=-tmp_yj
  write(*,*) 'phase of the point near the end-point of the computational region is (rad) ', tmp_yj

  do j=starting_surface_number,ending_surface_number
     sum4=0
     do k=1,mhtot
        mhk=mh_low+(k-1) 
        sum4=sum4+ypath(k+mhtot,j)*exp(ii*mhk*midplane_theta(j))*exp(ii*nh*delta_q_midplane(j))
     enddo
     sum4=sum4/sum3
     phase_midplane(j)=acos(real(sum4)/abs(sum4))
     if(imag(sum4).lt.0.) phase_midplane(j)=-phase_midplane(j)
  enddo


  do j=starting_surface_number,ending_surface_number
     do i=1,mhtot
        mhk=mh_low+(i-1)
        xir_midplane(i,j)=ypath(i+mhtot,j)*exp(ii*mhk*midplane_theta(j))*exp(ii*nh*delta_q_midplane(j)) !measured on the low-field side of the midplane,
     enddo
  enddo


end subroutine calculate_phase
