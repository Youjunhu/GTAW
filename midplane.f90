subroutine midplane() !calculate the values of poloidal angle on the low-field side of the midplane
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  use flux_grids,only:mpoloidal,nflux,r_new,z_new !input
  use radial_module,only:r_axis,z_axis,pfn,tfn !input
  use flux_grids,only: x0,midplane_theta !as output
  implicit none

  integer:: i,j,k
  real(p_):: ds1,ds2,theta(mpoloidal) !,rs(mpoloidal),zs(mpoloidal)

  do j=1,nflux !determine the major radius of the point on the low-field side of every flux surface
!!$     do i=1,mpoloidal
!!$        rs(i)=r_new(i,j)
!!$        zs(i)=z_new(i,j)
!!$     enddo
     call major_radius_on_lfs_midplane(mpoloidal,r_new(1,j),z_new(1,j),r_axis,z_axis,x0(j))
     !call major_radius_on_lfs_midplane(mpoloidal,rs,zs,r_axis,z_axis,x0(j))
     !write(*,*) 'j=',j,'x0(j)=',x0(j)
  enddo

  do i=1,mpoloidal 
     theta(i)=0.0_p_+twopi/(mpoloidal-1)*(i-1)
  enddo

  do j=1,nflux
     k=-1 ! a flag integer
     do i=1,mpoloidal-1
        if((r_new(i,j).gt.r_axis) .and. ((z_new(i,j)-z_axis)*(z_new(i+1,j)-z_axis).le.zero)) then !check whether it is the low-field side of the midplane
           k=i
           exit
        endif
     enddo
     if(k.eq.-1) stop 'error in finding points on the midplane'
     ds1=sqrt((x0(j)-r_new(k,j))**2+(z_axis-z_new(k,j))**2)
     ds2=sqrt((x0(j)-r_new(k+1,j))**2+(z_axis-z_new(k+1,j))**2) 
     midplane_theta(j)=(theta(k)*ds2+ theta(k+1)*ds1)/(ds1+ds2)
     !write(*,*) 'k=', k, 'ds1=',ds1,'ds2=',ds2
     !write(*,*) 'k=', k, "x0(j)-r_new(k,j)", x0(j)-r_new(k,j), "z_axis-z_new(k,j)", z_axis-z_new(k,j)

     !write(*,*) midplane_theta(j)
!!$    midplane_theta(j)=zero

  enddo

  open(11,file='midplane.txt')
  do j=1,nflux
     write(11,*) x0(j),z_axis,midplane_theta(j), sqrt(pfn(j)),sqrt(tfn(j))
  enddo
  close(11)

end subroutine midplane
