subroutine arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,x_axis,z_axis)
  !replacing one point near the low-field side of the midplane by a point that is exactly on the low-field side of the midplane
  !then arrange the arrays x_lcfs_new and z_lcfs_new so that (x_lcfs_new(1),z_lcfs_new(1)) is exactly on the low-field-side of the midplane
  !midplane is defined as the z=z_axis plane, where z_axis the z coordinate of the magnetic axis
  use precision,only:p_
  use radial_module,only: r_major,r_minor,eps !as output
  implicit none
  integer,intent(in):: np_lcfs
  real(p_),intent(inout):: x_lcfs(np_lcfs),z_lcfs(np_lcfs)
  real(p_),intent(in):: x_axis,z_axis
  real(p_):: x_lcfs_new(np_lcfs),z_lcfs_new(np_lcfs) 
  real(p_):: Rout

  integer:: kk(1),k,i


  open(321,file='lcfs.txt')
  do i=1,np_lcfs
     write(321,*) x_lcfs(i), z_lcfs(i)
  enddo
  close(321)


  !write(*,*) ' Z values of the uppermost point of LCFS: ', maxval(z_lcfs)
  !write(*,*) ' Z values of the lowest point of LCFS: ' ,minval(z_lcfs)

  !set the starting point of LCFS to be at the low-field-side of the midplane
  !maxval(x_lcfs)
  kk=maxloc(x_lcfs) !return the index of the array for which R is the largest,in order to determine the low-field side of the midplane
  k=kk(1)
  !k=10
  !write(*,*) 'index of the point on the lcfs that have the largest R, k=',k
  !if((z_lcfs(k+1)-z_axis)*(z_lcfs(k-1)-z_axis)>0) stop 'error in selecting the point on LCFS'

  call major_radius_on_lfs_midplane(np_lcfs,x_lcfs,z_lcfs,x_axis,z_axis,Rout)

  r_major=(maxval(x_lcfs)+minval(x_lcfs))/2._p_ !by definition
  r_minor=(maxval(x_lcfs)-minval(x_lcfs))/2._p_ !by definition

  ! write(*,*) 'inverse aspect ratio of LCFS is (definied at low-field side) ', (Rout-x_axis)/x_axis
  write(*,*) 'r_axis=',x_axis, 'r_major=', r_major, 'r_minor=',r_minor
  eps= r_minor/r_major !standard definition
  write(*,*) 'inverse aspect ratio of LCFS (i.e., r_minor/r_major) is ', eps
  write(*,*) 'ellipticity (elongation) of LCFS is ', (maxval(z_lcfs)-minval(z_lcfs))/2._p_/r_minor
  write(*,*) 'upper triangularity of LCFS is ', (r_major-x_lcfs(maxloc(z_lcfs)))/r_minor, &
             & 'lower triangularity of LCFS is ', (r_major-x_lcfs(minloc(z_lcfs)))/r_minor
  !replace one point of LCFS with the new point
  x_lcfs(k)=Rout
  z_lcfs(k)=z_axis

  !arrange the arrays x_lcfs_new and z_lcfs_new so that it starts from the low-field-side of the midplane
  do i=1,np_lcfs
     if(k+i-1.le.np_lcfs) then
        x_lcfs_new(i)=x_lcfs(k+i-1)
        z_lcfs_new(i)=z_lcfs(k+i-1)
     else
        x_lcfs_new(i)=x_lcfs(k+i-np_lcfs)
        z_lcfs_new(i)=z_lcfs(k+i-np_lcfs)
     endif
  enddo

  !use x_lcfs and z_lcfs to store the new data
  x_lcfs=x_lcfs_new
  z_lcfs=z_lcfs_new


  open(123,file='lcfs2.txt')
  do i=1,np_lcfs
     write(123,*) x_lcfs(i),z_lcfs(i)
  enddo
  close(123)


end subroutine arrange_lcfs
