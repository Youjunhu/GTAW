subroutine major_radius_on_lfs_midplane(mpoloidal,rs,zs,r_axis,z_axis,Rout)
  !given a flux surface, this subroutine determines the major radius of the point on the low field side of the middle plane
  use precision,only:p_
  implicit none
  integer,intent(in):: mpoloidal
  real(p_),intent(in):: rs(mpoloidal), zs(mpoloidal),r_axis,z_axis
  real(p_),intent(out):: Rout !the major radius of the point on the low-field-side of the mid-plane
  integer:: i,j,k1,k2,n
  real(p_):: r_select(mpoloidal),z_select(mpoloidal)
  !real(p_),dimension(:),allocatable:: x,z,tmp_y2
  real(p_):: tmp

 !select the low-field side (i.e., the part with r larger than r_axis) of a flux surface
  n=1
  do i=1,mpoloidal-1 
     if(rs(i).gt.r_axis) then
        r_select(n)=rs(i)
        z_select(n)=zs(i)
        n=n+1
     endif
  enddo
!  write(*,*) 'n-1= ', n-1
  !order the array according to the value of z_select
  do i=1,n-1
     do j=i+1,n-1
        if(z_select(j).le.z_select(i)) then
           !exchange the z value 
           tmp=z_select(i)
           z_select(i)=z_select(j)
           z_select(j)=tmp
           !also exchange the r value (I forgot this step in the older version, which causes a serious mistake)
           tmp=r_select(i)
           r_select(i)=r_select(j)
           r_select(j)=tmp
        endif
     enddo
  end do

!!$  allocate(x(n-1))
!!$  allocate(z(n-1))
!!$  allocate(tmp_y2(n-1))
!!$
!!$  do i=1,n-1
!!$     x(i)=r_select(i)
!!$     z(i)=z_select(i)
!!$!     write(*,*) x(i),z(i)
!!$  enddo
!!$
!!$  !the above steps select out the low-field part of a flux surface, so that z-x is a single-valued function
!!$  call spline(z,x,n-1,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
!!$  call splint(z,x,tmp_y2,n-1,z_axis,Rout) !Rout is the output

!!$  deallocate(x)
!!$  deallocate(z)
!!$  deallocate(tmp_y2)

  call linear_1d_interpolation(n-1,z_select,r_select,z_axis,Rout)  

end subroutine 



