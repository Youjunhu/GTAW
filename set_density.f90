subroutine set_density(nflux,pfn,psi_axis,psi_lcfs,rho)
  !set number density profile of electrons, ne. then we use charge neutality to estimate the ion number density and then mass density of the plasma
  !If there is only carbon impurity ions, we can use rhom=ne*mass_ion since the ratio of mass to charge of carbon is idential to the deturium.
  !if there are other impurity ions that have differnt ratio of mass to charge, then we need to know their compostion to know thier contribution to the mass density
  !set radial mass density profile by read a file and interpolate the data to desired radial location
  use precision,only:p_
  implicit none
  integer,intent(in)::nflux
  real(p_),intent(in):: pfn(nflux),psi_axis,psi_lcfs
  real(p_),intent(out):: rho(nflux) !mass density profile array
  real(p_),dimension(:),allocatable::  r,z,ne

  logical:: use_cylindrical_data
  character(100):: cylindrical_ne_file,density_file,prof_grid_type
  logical::prof_impurity_density_available,prof_grids_data_available
  integer,parameter:: max_num=300
  real(p_)::r_tmp(max_num),z_tmp(max_num),ne_tmp(max_num)
  integer:: ndata,j
  real(p_):: mass_number_ion,mi,unit_of_ne
  namelist/electron_density_profile_namelist/mass_number_ion,use_cylindrical_data,cylindrical_ne_file,&
       & density_file,unit_of_ne,prof_impurity_density_available,prof_grid_type,prof_grids_data_available

  open(11,file='gtaw.in')
  read(11,electron_density_profile_namelist)
  close(11)
  write(*,electron_density_profile_namelist)
  mi=mass_number_ion*1.6726d-27    !mass of ions (kg)

  if (use_cylindrical_data .eqv. .true.) then  !Thomoson scattering measurement is on the points along the line (r(i),z(i),i=1,ndata)
     open(11,file=cylindrical_ne_file)
     do j=1,max_num
    !    read(11,*,end=111) r_tmp(j),z_tmp(j),ne_tmp(j)
    read(11,*,end=111) r_tmp(j),ne_tmp(j)
     enddo
111  close(11)
     ndata=j-1 !ndata is the number of data points where the the TS measured ne is available and reliable
     z_tmp=0._p_ !assume the measurement is made on the midplane
     r_tmp=r_tmp/100. !change from cm to meter

     allocate(r(ndata))
     allocate(z(ndata))
     allocate(ne(ndata))

     do j=1,ndata
        r(j)=r_tmp(j)
        z(j)=z_tmp(j)
        ne(j)=ne_tmp(j)*unit_of_ne
     enddo

!!$    open(11,file='ts_com.txt')
!!$     do i=1,ndata
!!$        write(11,*) r(i),z(i),ne(i)
!!$     enddo
!!$     close(11)

     call from_cylindrical_to_flux_coordinates(mi,ndata,r,z,ne,nflux,pfn,psi_axis,psi_lcfs,rho)
  else
     call set_density_profile(mi,density_file,unit_of_ne,prof_impurity_density_available,&
          & prof_grid_type,prof_grids_data_available,nflux,pfn,psi_axis,psi_lcfs,rho)
  endif
end subroutine set_density


subroutine set_density_profile(mi,density_file,unit_of_ne,prof_impurity_density_available,&
     & prof_grid_type,prof_grids_data_available,nflux,pfn,psi_axis,psi_lcfs,rho)
  !set radial mass density profile by readin a file of the radial density profile and interpolate the data to desired radial location
  use precision,only:p_
  use constants,only: one,two, twelve
  use radial_module,only:tfn !used only for plotting
  implicit none
  integer,intent(in)::nflux
  real(p_),intent(in):: mi,pfn(nflux),psi_axis,psi_lcfs,unit_of_ne
  character(*),intent(in):: density_file,prof_grid_type
  logical,intent(in)::prof_impurity_density_available,prof_grids_data_available
  real(p_),intent(out):: rho(nflux) !mass density profile array
  real(p_):: number_density_new(nflux),number_density2_new(nflux) !number density of ions
  integer,parameter:: max_num=3000
  real(p_):: radialc(max_num),tmp_density(max_num),tmp_density2(max_num)
  integer::ndata
  real(p_),dimension(:),allocatable::  pfn_old,tfn_sqrt,pfn_sqrt
  real(p_),dimension(:),allocatable::  number_density_old,number_density2_old,tmp_y2
  !real(p_):: number_density_old(ndata),tmp_y2(ndata),pfn_old(ndata)
  real(p_):: y_tmp, tmp_y2b(nflux)
  !  real(p_),parameter:: mi=two*1.6726d-7    !10^20 times mass of deuterium ions (kg)
  !real(p_),parameter:: mi=two*1.6726d-27    !mass of deuterium ions (kg)
  real(p_),parameter:: mi2=twelve*1.6726d-27    !mass of impurity ions, (here carbon-12 ions) (kg), 
  integer:: j


  open(11,file=density_file)

  if(prof_impurity_density_available.eqv..true.) then
     if(prof_grids_data_available.eqv..true.) then
        do j=1,max_num
           read(11,*,end=111) radialc(j), tmp_density(j),tmp_density2(j)
        enddo
     else
        do j=1,max_num
           read(11,*,end=111)             tmp_density(j),tmp_density2(j)
        enddo
     endif

  else
     if(prof_grids_data_available.eqv..true.) then     
        do j=1,max_num
           read(11,*,end=111) radialc(j), tmp_density(j) 
        enddo
     else
        do j=1,max_num
           read(11,*,end=111)             tmp_density(j) 
        enddo
     endif
  endif
111 close(11)

  ndata=j-1
  !  write(*,*) 'number of data of the density radial profile=',ndata
  if(ndata.le.1) stop 'please provide the profile of the electron number density'

  allocate(number_density_old(ndata))
  allocate(number_density2_old(ndata))
  allocate(tfn_sqrt(ndata))
  allocate(pfn_sqrt(ndata))
  allocate(pfn_old(ndata))
  allocate(tmp_y2(ndata))


  do j=1,ndata
     number_density_old(j)=tmp_density(j)*unit_of_ne
     number_density2_old(j)=tmp_density2(j)*unit_of_ne
     !rotate_old(j)=tmp_rotate(j)
     !pfn_old(j)=psi_axis+(j-1)*(psi_lcfs-psi_axis)/(ndata-1)
  enddo


  if (prof_grids_data_available.eqv..false.) then
     do j=1,ndata 
        radialc(j)=0.+one/(ndata-1)*(j-1) !uniform radial grids is assumed if it is not provided in profile files
     enddo
  endif


  if (trim(prof_grid_type).eq.'tfn_sqrt') then   !--for ni defined on uniform sqrt(toroidal_flux) grids----for DIIID gaprofile data-
     tfn_sqrt=radialc

     call spline(sqrt(tfn),pfn,nflux,2.d30,2.d30,tmp_y2b) !prepare the second order derivative needed in the cubic spline interpolation
     do j=1,ndata !interpolating to get the corresponding pfn
        call splint(sqrt(tfn),pfn,tmp_y2b,nflux,tfn_sqrt(j),pfn_old(j))
     enddo

  else if(trim(prof_grid_type).eq.'pfn_sqrt') then   
     pfn_old=radialc**2
  else if(trim(prof_grid_type).eq.'pfn') then   
     pfn_old=radialc

  else 
     stop 'please specify the type of the radial grids used in the profile file'
  endif


  call spline(pfn_old,number_density_old,ndata,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation

  do j=1,nflux
     call splint(pfn_old,number_density_old,tmp_y2,ndata,pfn(j),number_density_new(j)) !interpolate to the new radial grid points pfn()
  enddo

  !stop 'stop in set_density.f90'

  if(prof_impurity_density_available.eqv..true.) then
     call spline(pfn_old,number_density2_old,ndata,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
     do j=1,nflux
        call splint(pfn_old,number_density2_old,tmp_y2,ndata,pfn(j),number_density2_new(j)) !interpolate to the new radial grid points
     enddo
  else
     number_density2_new=0._p_
  endif

  open(12,file='electron_density_compare.txt')
  do j=1,ndata
     write(12,*) pfn_old(j),number_density_old(j)
  enddo
  write(12,*)
  write(12,*)
  do j=1,nflux
     write(12,*) pfn(j),number_density_new(j)
  enddo
  close(12)

  open(12,file='ne_profile.txt')
  do j=1,nflux
     write(12,*) sqrt(pfn(j)),sqrt(tfn(j)),number_density_new(j)
  enddo
  close(12)

  !in the past, the number density of ions is infered from the charge neutral condition: ni=ne, and I did not take into acount of the contribution of impurity ions to the mass
  !later (2015-6-1), when I deal with DIIID data, ion number desnity of ions are read in from gaprofile and the contribution of impurity ions to the mass is taken into account
  !rho=number_density_new*mi !mass density in SI units: kg*m^-3
  rho=number_density_new*mi+number_density2_new*mi2 !mass density in SI units: kg*m^-3
  write(*,*) 'Electron density at the magnetic axis is (10^19m^-3)', number_density_new(1)/(1.d19)
end subroutine set_density_profile


subroutine from_cylindrical_to_flux_coordinates(mi,ndata,r,z,ne,nflux,pfn,psi_axis,psi_lcfs,rho) !need checking
  use precision,only:p_
  use constants,only: one,two
  use radial_module,only:tfn 
  implicit none
  integer,intent(in)::ndata
  real(p_),intent(in):: mi
  real(p_),intent(in):: r(ndata),z(ndata),ne(ndata)
  integer,intent(in)::nflux
  real(p_),intent(in):: pfn(nflux),psi_axis,psi_lcfs
  real(p_),intent(out):: rho(nflux) !mass density profile array
  real(p_):: psi_uniform(ndata),number_density_uniform_psi(ndata)
  real(p_):: number_density_new(nflux) !number density of electrons in unit of 10^20m^-3
  real(p_):: number_density_old(ndata),pfn_old(ndata),tmp_y2(ndata)
!  real(p_),parameter:: mi=two*1.6726d-7    !10^20 times mass of deuterium ions (kg)
  real(p_):: psi_func !2D polodial magnetic flux function on (R,Z) plane
  integer:: j,k
  real(p_):: tmp

  do j=1,ndata
     pfn_old(j)=psi_func(r(j),z(j))
     number_density_old(j)=ne(j)/(1.d20) !normalized by 10^20m^-3
  enddo

  pfn_old=(pfn_old-psi_axis)/(psi_lcfs-psi_axis) !normalized poloidal magnetic flux

  do k=1,ndata !sort the pfn_old so that the value of its element is increasing with index i
     do j=1,ndata-k
        if(pfn_old(j)>pfn_old(j+1)) then
           tmp=pfn_old(j+1)
           pfn_old(j+1)=pfn_old(j)
           pfn_old(j)=tmp
           tmp=number_density_old(j+1)
           number_density_old(j+1)=number_density_old(j)
           number_density_old(j)=tmp
        endif
     enddo
  enddo


  open(11,file='ts_p.txt')
  do j=1,ndata
     write(11,*) pfn_old(j),number_density_old(j)
  enddo
  close(11)

  call spline(pfn_old,number_density_old,ndata,2.d30,10._p_,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation

  do j=1,ndata
     psi_uniform(j)=0._p_+one/(ndata-1)*(j-1)
     call splint(pfn_old,number_density_old,tmp_y2,ndata,psi_uniform(j),number_density_uniform_psi(j)) !to get fpsi array corresponding to the new radial grid points
  enddo

  open(11,file='ne_data_for_mega') !this is the ne (in unit of 10^20m^-3) corresponding to uniform poloidal magnetic flux, this file provide a input file for  MEGA code
  do j=1,ndata
     write(11,*) number_density_uniform_psi(j)
  enddo
  close(11)

  do j=1,nflux
     call splint(pfn_old,number_density_old,tmp_y2,ndata,pfn(j),number_density_new(j)) 
  enddo

!the above spline does not work for the ts density data, so I switch to using the linear interpolation
!!$  do j=1,ndata
!!$     psi_uniform(j)=0._p_+one/(ndata-1)*(j-1)
!!$     call linear_1d_interpolation(ndata,pfn_old,number_density_old,psi_uniform(j),number_density_uniform_psi(j))
!!$  enddo
!!$
!!$  open(11,file='ne_data_for_mega') !this is the ne (in unit of 10^20m^-3) corresponding to uniform poloidal magnetic flux, this file provide a input file for  MEGA code
!!$  do j=1,ndata
!!$     write(11,*) number_density_uniform_psi(j)
!!$  enddo
!!$  close(11)
!!$
!!$  do j=1,nflux
!!$     call linear_1d_interpolation(ndata,pfn_old,number_density_old,pfn(j),number_density_new(j))
!!$  enddo


  open(12,file='ne_profile_ts.txt')
  do j=1,nflux
     write(12,*) sqrt(pfn(j)),sqrt(tfn(j)),number_density_new(j)
  enddo
  close(12)


  !number density of ions is infered from the charge neutral condition: ni=ne
  rho=number_density_new*(1.d20)*mi !mass density in SI units: kg*m^-3


end subroutine from_cylindrical_to_flux_coordinates
