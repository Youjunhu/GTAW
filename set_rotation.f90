subroutine set_rotation(rotation_file_name,prof_grids,nflux,psival_new,psi_axis,psi_lcfs,rotation_new)
  !set rotation profile by readin a file and interpolate the data to desired radial location, not tested yet?
  use precision,only:p_
  use constants,only: one,two
  use radial_module,only:poloidal_flux_normalized,toroidal_flux !used only for plotting
  use radial_module,only: r_axis
  implicit none
  integer,intent(in)::nflux
  real(p_),intent(in):: psival_new(nflux),psi_axis,psi_lcfs
  character(*),intent(in):: rotation_file_name,prof_grids
  real(p_),intent(out):: rotation_new(nflux)
  integer,parameter:: max_num=3000
  real(p_):: tmp_rotation(max_num)
  integer::ndata
  real(p_),dimension(:),allocatable::  rotation_old,psival_old,tmp_y2
  real(p_):: y_tmp, tmp_y2b(nflux),tf_sqrt
  integer:: j

  open(11,file=rotation_file_name)
  do j=1,max_num
     read(11,*,end=111) tmp_rotation(j) 
  enddo
111 close(11)

  ndata=j-1
  write(*,*) 'number of data of the density radial profile=',ndata

  allocate(rotation_old(ndata))
  allocate(psival_old(ndata))
  allocate(tmp_y2(ndata))

  do j=1,ndata
     rotation_old(j)=tmp_rotation(j)/r_axis
     !psival_old(j)=psi_axis+(j-1)*(psi_lcfs-psi_axis)/(ndata-1)
  enddo

  if (trim(prof_grids).eq.'tfn_sqrt') then   !for DIIID
     call spline(sqrt(toroidal_flux),psival_new,nflux,2.d30,2.d30,tmp_y2b) !prepare the second order derivative needed in the cubic spline interpolation
     do j=1,ndata
        tf_sqrt=0.+one/(ndata-1)*(j-1) !uniform sqrt(toroidal_flux) grids
        call splint(sqrt(toroidal_flux),psival_new,tmp_y2b,nflux,tf_sqrt,psival_old(j))
     enddo
  else if(trim(prof_grids).eq.'pfn') then   
     psival_old(j)=psi_axis+(j-1)*(psi_lcfs-psi_axis)/(ndata-1)
  else 
     stop 'please specify the type of the radial grids used in the profile file'
  endif

  call spline(psival_old,rotation_old,ndata,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do j=1,nflux
     call splint(psival_old,rotation_old,tmp_y2,ndata,psival_new(j),rotation_new(j)) !interpolate to the new radial grid points
  enddo

  open(12,file='rotation_profile.txt')
  do j=1,nflux
     write(12,*) sqrt(poloidal_flux_normalized(j)),sqrt(toroidal_flux(j)),rotation_new(j)
  enddo
  close(12)
write(*,*) 'r_axis=',r_axis, 'in set rotation subroutine'
end subroutine set_rotation
