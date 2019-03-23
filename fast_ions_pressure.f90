subroutine fast_ions_pressure(pfn,fpsi,nflux,mpoloidal,r_mc,z_mc,b_axis,psi_axis,psi_lcfs,va0)
  !calculate the parallel and perpendicular pressure of fast ions (for mega code)
  !Use spherical coordinates in velocity space, assume axisymmetrical velcoity distriubtion about the local magnetic field.
  use precision,only:p_
  use constants,only:one,two,twopi,pi,mu0
  use fast_ions,only: md,kev,me,charge
  implicit none
  integer,intent(in):: nflux,mpoloidal
  real(p_),intent(in):: pfn(nflux),fpsi(nflux)
  real(p_),intent(in):: r_mc(mpoloidal,nflux),z_mc(mpoloidal,nflux)
  real(p_),intent(in):: b_axis,psi_axis,psi_lcfs,va0 !va0 is Alfven speed at the magnetic axis

  real(p_):: b_mc(mpoloidal,nflux), ppara_mc(mpoloidal,nflux),pperp_mc(mpoloidal,nflux)
  integer,parameter:: nx=128,nz=128
  real(p_):: r(nx),z(nz)
  real(p_)::pfn_xz(nx,nz),bp_xz(nx,nz),fpsi_xz(nx,nz),b(nx,nz) !pfn_xz is the normalized poloidal magnetic flux
  real(p_):: psi_func, psi_gradient_func !function names
  real(p_),parameter:: left=1.42635d+00,right=2.29851d+00,bottom=-6.97372d-01,top=6.38665d-01 ! computational box on poloidal plane used in MEGA simulation (these values are the output of 'to-MEGA.f' code)
  real(p_)::ppara(nx,nz),pperp(nx,nz),ppara_mega(nx,nz),pperp_mega(nx,nz)
  real(p_),parameter:: psi_scale=0.3_p_, energy_birth=60._p_ !in unti of keV
  real(p_):: stored_eng_desired_in_mega_unit=1.74378d+02
  real(p_):: stored_eng_desired
  real(p_),parameter::clambda0=0.88_p_,dclambda=0.1_p_ 
  real(p_),parameter:: te=2._p_   !electron temperature in unit of kev, used in calculating the critical velocity of fast ions
  real(p_)::  vcrit !critical velocity
  real(p_):: vbirth,vmin,vmax,deltav
  integer,parameter::m=100,n=100 !number of grids in velocity space (spherical coordinates (v, theta))
  real(p_):: theta(m),v(n),dv,dtheta,vac
  integer:: i,j,ii,jj

  real(p_):: dr,dz !rectangle computational box used by MEGA code
  real(p_):: stored_eng,df_norm_factor,density0
  real(p_):: tmp_y2(nflux) !used in spline interpolation
  real(8):: Omega_h ,mega_energy_unit1,mega_energy_unit2,mega_pressure_unit


  omega_h=b_axis*charge/md !check this with that given in file mega_job.oxxxxx, in unit of rad/second
  vcrit=sqrt(two*te*kev/me)*(3*sqrt(pi)/(4.*(md/me)))**0.3333333333_p_ !critical velocity, refer to the formula in my notes
  vbirth=sqrt(two*energy_birth*kev/md)

  deltav=0.05_p_*va0 !according to mega code
  vmin=vbirth*1.0d-1  !according to mega code
  vmax=vbirth + 3.0d0*deltav !according to mega code

  write(*,*) 'fast ions parameters: '
  write(*,*) 'vbirth (10^6m/s)=',vbirth/1.0d6, 'vcrit (10^6m/s)=',vcrit/1.0d6
  write(*,*) 'vbirth/va0=',vbirth/va0,'vcrit/va0=',vcrit/va0
  write(*,*) "cutoff width near vbirth (10^6m/s)",deltav/1.0d6
  mega_energy_unit1=(va0/omega_h)**3*b_axis**2/mu0
  mega_energy_unit2=b_axis**2/(va0**2*mu0*md)*(va0/omega_h)**3*md*va0**2
  write(*,*) 'MEGA energy unit (kJ): ', mega_energy_unit2/1000._p_
  write(*,*) 'MEGA energy unit (kJ): ', mega_energy_unit2/1000._p_
  write(*,*) 'MEGA length unit (meter): ', va0/omega_h
  mega_pressure_unit=b_axis**2/mu0  !not b0**2/(2*mu0)

  stored_eng_desired=stored_eng_desired_in_mega_unit*mega_energy_unit1 !to SI unit J

  dtheta=pi/(m-1)
  do i=1,m !pitch angle array
     theta(i)=0._p_+dtheta*(i-1)
  enddo

  dv=(vmax-vmin)/(n-1)
  do j=1,n !velocity array
     v(j)=vmin+dv*(j-1)
  enddo

  dr=(right-left)/(nx-1)
  do ii=1,nx
     r(ii)=left+dr*(ii-1)
  enddo
  dz=(top-bottom)/(nz-1)
  do jj=1,nz
     z(jj)=bottom+dz*(jj-1)
  enddo

  do ii=1,nx
     do jj=1,nz
        pfn_xz(ii,jj)=(psi_func(r(ii),z(jj))-psi_axis)/(psi_lcfs-psi_axis) !normalized poloidal flux
        bp_xz(ii,jj)=psi_gradient_func(r(ii),z(jj))/r(ii)
     enddo
  enddo

  call spline(pfn,fpsi,nflux,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  do ii=1,nx
     do jj=1,nz
        call splint(pfn,fpsi,tmp_y2,nflux,pfn_xz(ii,jj),fpsi_xz(ii,jj)) 
     enddo
  enddo

  do ii=1,nx !calculate the strength of the magnetic field
     do jj=1,nz
        b(ii,jj)=sqrt(fpsi_xz(ii,jj)**2/r(ii)**2+bp_xz(ii,jj)**2)
     enddo
  enddo

  do ii=1,nx
     !$omp parallel do     
     do jj=1,nz
        call calculate_pressure(pfn_xz(ii,jj),b(ii,jj),ppara(ii,jj),pperp(ii,jj))
     enddo
     !$omp end parallel do
  enddo


  write(*,*) 'anisotropicity of the distribution,max(pperp/ppara),min(pperp/ppara)=', maxval(pperp/ppara),minval(pperp/ppara)

!!$  call calculate_pressure(0._p_,b_axis,ppara0,pperp0) !calculate the parallel and perpendicular pressure at the magnetic axis, ppara0 will be used as a normalzing factor
!!$  write(*,*) 'normalizing factor  ','ppara0=',ppara0, 'pperp0=',pperp0
  !  ppara=ppara/ppara0*beta_a0/two  !in mega pressure is normalized by b0**2/mu0 (not by b0**2/(2*mu0))
  !  pperp=pperp/ppara0*beta_a0/two

  stored_eng=0._p_ !calculate stored energy of fast ions
  do ii=1,nx
     do jj=1,nz
        !stored_eng= stored_eng+0.5_p_*(ppara(ii,jj)+two*pperp(ii,jj))*(b_axis**2/mu0)*dr*dz*twopi*r(ii)
        stored_eng= stored_eng+0.5_p_*(ppara(ii,jj)+two*pperp(ii,jj))*dr*dz*twopi*r(ii)
     enddo
  enddo
  write(*,*) 'before normalizing, the fast ions stored energy (kJ): ', stored_eng/1000._p_

  df_norm_factor=stored_eng_desired/stored_eng

  write(*,*) ' the normalizing factor= ',df_norm_factor
  ppara=ppara*df_norm_factor
  pperp=pperp*df_norm_factor


  ppara_mega=ppara/mega_pressure_unit !normalized to mega unit (b0**2/mu0)
  pperp_mega=pperp/mega_pressure_unit !normalized to mega unit (b0**2/mu0)
  call write_gnuplot_contour_data2(nx,nz,r,z,ppara_mega,'ppara.txt') !to be used by mega code
  call write_gnuplot_contour_data2(nx,nz,r,z,pperp_mega,'pperp.txt') !to be used by mega code

  stored_eng=0._p_ !calculate stored energy of fast ions to check the stored energy
  do ii=1,nx
     do jj=1,nz
        stored_eng= stored_eng+0.5_p_*(ppara(ii,jj)+two*pperp(ii,jj))*dr*dz*twopi*r(ii)
     enddo
  enddo
  write(*,*) 'fast ions stored energy (kJ): ', stored_eng/1000._p_
  write(*,*) 'fast ions stored energy in MEGA unit: ', stored_eng/mega_energy_unit1

  stored_eng=0._p_ !calculate stored energy of fast ions to check the stored energy within boundary flux surface
  do ii=1,nx
     do jj=1,nz
        if(pfn_xz(ii,jj)<(1.0-0.04)) then
           vac=1._p_
        else
           vac=0._p_
        endif
        stored_eng= stored_eng+0.5_p_*(ppara(ii,jj)+two*pperp(ii,jj))*dr*dz*twopi*r(ii)*vac
     enddo
  enddo

  write(*,*) 'fast ions stored energy (kJ) within boundary flux surface: ', stored_eng/1000._p_
  write(*,*) 'fast ions stored energy in MEGA unit: ', stored_eng/mega_energy_unit1

  call calculate_density(0._p_,b_axis,density0)
  write(*,*) 'density of fast ions at magnetic axis (10^19m^-3) is ',density0/(1.d19)

!calculate pressure on magnetic flux coordinate grids
  do jj=1,nflux
     !$omp parallel do     
     do ii=1,mpoloidal
        b_mc(ii,jj)=sqrt(psi_gradient_func(r_mc(ii,jj),z_mc(ii,jj))**2+fpsi(jj)**2)/r_mc(ii,jj)
        call calculate_pressure(pfn(jj),b_mc(ii,jj),ppara_mc(ii,jj),pperp_mc(ii,jj))
     enddo
     !$omp end parallel do
  enddo

  ppara_mc=ppara_mc*df_norm_factor/mega_pressure_unit !to mega unit
  pperp_mc=pperp_mc*df_norm_factor/mega_pressure_unit !to mega unit

  open(123,file='pressure_h.txt')
  do jj=1,nflux
     do ii=1,mpoloidal
        write(123,*) ii,ppara_mc(ii,jj),pperp_mc(ii,jj)
     enddo
        write(123,*)
        write(123,*)
  enddo
  close(123)

call fourier_tansformation(mpoloidal,nflux,pfn,ppara_mc, 'ppara_h00.txt') 
call fourier_tansformation(mpoloidal,nflux,pfn,pperp_mc, 'pperp_h00.txt') 


contains
  function f(psi,v,clambda) result (z) !distribution function of fast ions
    implicit none
    real(p_):: psi,v,clambda,z

    z=exp(-psi/psi_scale) &
         & *1._p_/(v**3 + vcrit**3)*0.50d0*erfc((v-vbirth)/deltav) &
          *exp(-(clambda-clambda0)**2/dclambda**2)
  end function f


  subroutine calculate_pressure(pfn_val,b_val,ppara_val,pperp_val)
    implicit none
    real(p_):: pfn_val,b_val,ppara_val,pperp_val
    real(p_)::clambda,kenergy,mu, dvol
    integer:: i,j

    ppara_val=0._p_ 
    pperp_val=0._p_
    do i=1,m
       do j=1,n
          kenergy=md*v(j)**2/two
          mu=md*(v(j)*sin(theta(i)))**2/(two*b_val)
          clambda=mu/kenergy*b_axis
          dvol=twopi*v(j)**2*sin(theta(i))*dtheta*dv !volume element in spherical coordinates (velocity space)
          ppara_val=ppara_val+md*(v(j)*cos(theta(i)))**2*f(pfn_val,v(j),clambda)*dvol
          pperp_val=pperp_val+md*(v(j)*sin(theta(i)))**2*f(pfn_val,v(j),clambda)*dvol
       enddo
    enddo
    pperp_val=pperp_val/two !dividing by two because the perpendicular pressure is contributed by two directions (for axisymmetrical velcoity distriubtion)

  end subroutine calculate_pressure

  subroutine calculate_density(pfn_val,b_val,density_val)
    implicit none
    real(p_):: pfn_val,b_val,density_val
    real(p_):: kenergy,mu,clambda,dvol
    integer:: i,j

    !    write(*,*) 'df_norm_factor=',df_norm_factor
    !write(*,*) 'pfn_val=',pfn_val ,'b_val=',b_val

    density_val=0._p_ 
    do i=1,m
       do j=1,n
          kenergy=md*v(j)**2/two
          mu=md*(v(j)*sin(theta(i)))**2/(two*b_val)
          clambda=mu/kenergy*b_axis
          dvol=twopi*v(j)**2*sin(theta(i))*dtheta*dv !volume element in velocity space
          density_val=density_val+f(pfn_val,v(j),clambda)*df_norm_factor*dvol
       enddo
    enddo
  end subroutine calculate_density

end subroutine fast_ions_pressure


subroutine fourier_tansformation(mpoloidal,nflux,pfn,kappas_in,filename) 
!modified form filter_theta subroutine
  use precision,only:p_
  use constants,only: zero,one,two,pi,twopi
  implicit none
  integer,intent(in)::mpoloidal,nflux
  real(p_),intent(in)::pfn(nflux), kappas_in(mpoloidal,nflux)
  character(*),intent(in):: filename
!  real(p_),intent(out):: kappas_out(mpoloidal,nflux)

  integer,parameter:: n=5 !number of terms included in the Fourier series
  real(p_):: a(0:n),b(0:n)  !fourier coefficients
  real(p_):: dth,theta(mpoloidal),data(mpoloidal),output(mpoloidal) 
  real(p_):: sum1,sum2
  integer:: i,j,k
integer,parameter:: nrad_mega=101
real(p_):: sqrt_pfn(nflux), pressure_h00(nflux),tmp_y2(nflux)
real(p_)::sqrt_pfn_new(nrad_mega),pressure_h00_new(nrad_mega)

  dth=twopi/(mpoloidal-1)
  do i=1,mpoloidal
     theta(i)=zero+dth*(i-1)
  enddo

open(14,file=filename)

  do j=1,nflux

     do i=1,mpoloidal
        data(i)=kappas_in(i,j)
     enddo

     !calculate the Fourier coefficients
     do k=0,n
        sum1=zero
        sum2=zero
        do i=1,mpoloidal-1 !the integration in the Fourier coefficients
           sum1=sum1+data(i)*cos(k*theta(i))*dth
           sum2=sum2+data(i)*sin(k*theta(i))*dth
        enddo
        a(k)=sum1/pi
        b(k)=sum2/pi
     enddo

!write(14,*) sqrt(pfn(j)), a(0)/2,b(0),(a(k),b(k),k=1,n)

pressure_h00(j)=a(0)/2._p_


!!$     do k=0,6
!!$        write(*,*) k,a(k),b(k)
!!$     enddo

!!$     b(1)=b(1)*(-120.0) !artificially change the coefficient, alter this to get different up-down asymmetricity
!!$
!!$     do i=1,mpoloidal !reconstruct the original function
!!$        output(i)=a(0)/two
!!$        do k=1,n
!!$           output(i)=output(i)+a(k)*cos(k*theta(i))+b(k)*sin(k*theta(i))
!!$        enddo
!!$     enddo
!!$
!!$     do i=1,mpoloidal !store the filtered datas
!!$        kappas_out(i,j)=output(i)
!!$     enddo

  enddo

sqrt_pfn=sqrt(pfn)

do j=1,nrad_mega
sqrt_pfn_new(j)=0.+1._p_/(nrad_mega-1)*(j-1)
enddo

call spline(sqrt_pfn,pressure_h00,nflux,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation



do j=1,nrad_mega
     call splint(sqrt_pfn,pressure_h00,tmp_y2,nflux,sqrt_pfn(j),pressure_h00_new(j)) 
  enddo

do j=1,nrad_mega
write(14,*) sqrt_pfn_new(j), pressure_h00_new(j)
enddo
close(14)

end
