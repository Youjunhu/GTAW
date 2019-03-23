subroutine provide_equ_to_GEM_code(nflux,mpoloidal,r_new,z_new,ra,qpsi,fpsi)
  !---write data to be used by equil.f of GEM code----
  use precision,only:p_
  implicit none
  integer,intent(in):: nflux,mpoloidal
  real(p_),intent(in):: r_new(mpoloidal,nflux),z_new(mpoloidal,nflux)
  real(p_),intent(in):: ra(nflux),qpsi(nflux),fpsi(nflux)

  real(p_)::ne(nflux), te(nflux),ti(nflux),n_ep(nflux),t_ep(nflux)

  integer::i,j

  open(52,file='gtcAE_rdata.txt')
  do j = 1,nflux
     write(52,*) (r_new(i,j),i = 1,mpoloidal)
  end do
  close(52)

  open(52,file='gtcAE_zdata.txt')
  do j = 1,nflux
     write(52,*) (z_new(i,j),i = 1,mpoloidal)
  end do
  close(52)

  open(51,file='gtcAE_rho.txt')
  do j = 1,nflux
     write(51,*) ra(j)
  end do
  close(51)

  open(51,file='gtcAE_bdata_tor.txt')
  do j = 1,nflux
     write(51,*) (fpsi(j)/r_new(i,j),i = 1,mpoloidal)
  end do
  close(51)

  call set_profile_for_gem()
  open(51,file='gtcAE_1dProfiles.txt')
  do i = 1,nflux
     write(51,*) qpsi(i),ne(i),n_ep(i),te(i),t_ep(i),ti(i)
  end do
  close(51)

end subroutine provide_equ_to_GEM_code




subroutine set_profile_for_gem(density_file,unit_of_ne,prof_impurity_density_available,&
     & prof_grid_type,prof_grids_data_available,nflux,pfn,psi_axis,psi_lcfs,rho)
  !set radial mass density profile by readin a file of the radial density profile and interpolate the data to desired radial location
  use precision,only:p_
  use constants,only: one,two, twelve
  use radial_module,only:tfn !used only for plotting
  implicit none
  integer,intent(in)::nflux
  real(p_),intent(in):: pfn(nflux),psi_axis,psi_lcfs,unit_of_ne
  character(*),intent(in):: density_file,prof_grid_type
  logical,intent(in)::prof_impurity_density_available,prof_grids_data_available
  real(p_),intent(out):: rho(nflux) !mass density profile array
  real(p_):: number_density_new(nflux),number_density2_new(nflux) !number density of ions
  integer,parameter:: max_num=3000
  real(p_):: radialc(max_num),tmp_density(max_num),tmp_density2(max_num)
  integer::ndata
  real(p_),dimension(:),allocatable::  number_density_old,number_density2_old,pfn_old,tf_sqrt,tmp_y2
  !real(p_):: number_density_old(ndata),tmp_y2(ndata),pfn_old(ndata)
  real(p_):: y_tmp, tmp_y2b(nflux)
  integer:: j


  open(11,file=density_file)
  if(prof_grids_data_available.eqv..true.) then
     do j=1,max_num
        read(11,*,end=111) radialc(j), tmp_density(j),tmp_density2(j)
     enddo
  else
     do j=1,max_num
        read(11,*,end=111)             tmp_density(j),tmp_density2(j)
     enddo
  endif


111 close(11)

  ndata=j-1
  !  write(*,*) 'number of data of the density radial profile=',ndata
  if(ndata.le.1) stop 'please provide the profile of the electron number density'

  allocate(number_density_old(ndata))
  allocate(number_density2_old(ndata))
  allocate(tf_sqrt(ndata))
  allocate(pfn_old(ndata))
  allocate(tmp_y2(ndata))


  do j=1,ndata
     number_density_old(j)=tmp_density(j)*unit_of_ne
     number_density2_old(j)=tmp_density2(j)*unit_of_ne
     !rotate_old(j)=tmp_rotate(j)
     !pfn_old(j)=psi_axis+(j-1)*(psi_lcfs-psi_axis)/(ndata-1)
  enddo


  if (trim(prof_grid_type).eq.'tfn_sqrt') then   !--for ni defined on uniform sqrt(toroidal_flux) grids----for DIIID gaprofile data-
     if (prof_grids_data_available.eqv..true.) then
        do j=1,ndata
           tf_sqrt(j)=radialc(j)
        enddo
     else
        do j=1,ndata 
           tf_sqrt=0.+one/(ndata-1)*(j-1) !uniform sqrt(toroidal_flux) grids
        enddo
     endif

     call spline(sqrt(tfn),pfn,nflux,2.d30,2.d30,tmp_y2b) !prepare the second order derivative needed in the cubic spline interpolation
     do j=1,ndata !interpolating to get the corresponding pfn
        call splint(sqrt(tfn),pfn,tmp_y2b,nflux,tf_sqrt(j),pfn_old(j))
     enddo
  else if(trim(prof_grid_type).eq.'pfn') then

     if (prof_grids_data_available.eqv..true.) then
        do j=1,ndata
           pfn_old(j)=radialc(j)
        enddo
     else
        do j=1,ndata
           pfn_old(j)=0.+one/(ndata-1)*(j-1) !use uniform radial grids
        enddo
     endif
  else 
     stop 'please specify the type of the radial grids used in the profile file'
  endif

  call spline(pfn_old,number_density_old,ndata,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do j=1,nflux
     call splint(pfn_old,number_density_old,tmp_y2,ndata,pfn(j),number_density_new(j)) !interpolate to the new radial grid points pfn()
  enddo

end subroutine set_profile_for_gem
