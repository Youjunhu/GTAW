program GTAW
  !-------------------------------------------------------------------------------
  !Code name: GTAW (General Tokamak Alfven Wave)
  !GTAW is a Fortran code calculateing the MHD continuous spectrum and Alfven eigenmodes in general tokamak geometry.
  !GTAW reads a G-EQdsk file (one of the output files of EFIT code) to get the equilibrium information of axisymmetric tokamak plasmas
  !Then GTAW uses the equilibrium information to construct a magnetic surface coordinate system with desired Jacobian
  !In constructing the coordinate system, the code first finds out a series of magnetic surfaces, then construct desired poloidal grids.
  !Using the above equilibrium and flux coordinate system, GTAW calculates the MHD continuous spectrum and find global Alfven gap modes by using the ideal MHD model for the plasmas. 
  !For the details of the model and numerical method used in GTAW, refer to the following documents: 
  !Youjun Hu, et al., Phys. Plasmas 21, 052510 (2014)
  !http://theory.ipp.ac.cn/~yj/research_notes/mhd.pdf

  !GTAW is a free software, which means that you can redistribute it and/or modify it under the terms of the GNU General Public License.
  !The development of GTAW was supported by the National Magnetic Confinement Fusion Science Program of China under Grant No. 2013GB112010.
  !External Library used in GTAW: Lapack (Refer to makefile for compiling and running GTAW).
  !Code author: Youjun Hu (Email: yjhu@ipp.cas.cn), Institute of Plasmas Physics, Chinese Academy of Sciences
  !Input files: a namelist file "gtaw.in", G-file (filename is set in gtaw.in), a file for the electron number density profile (filename is set in gtaw.in)
  !Output files: continua.txt, which contians the Alfven continuous spectrum of the given equilibrium
  !Todo: to replace the spline interpolation with simple linear interpolation to improve the compuational efficiency
  !-----------------------------------------------------------------------------------
  use precision,only:p_
  use constants,only:zero,one,two,three,four,five,twopi,pi,mu0,gamma, slow_sound_approximation !,slow_sound_approximation_for_mode !gamma is polytrope index
  use poloidal_flux_2d,only:xarray,zarray,nx,nz,psi,psi_gradient,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx, &
       & y2a_psi,y2a_gradient,y2a_psi_x,y2a_psi_z,y2a_psi_xx,y2a_psi_zz,y2a_psi_xz,y2a_psi_zx !y2a_psi_* is a tempory array used in 2d cubic spline interpolation of psi_*, psi_x is the partial derivative with respect to x, and similar meaning for psi_x, psi_z etc
  use flux_grids,only: nflux,mpoloidal,r_new,z_new !,theta !, y2a_r_new,y2a_z_new ! r_new and z_new store the magnetic surfaces found in subroutine calculate_contour(), y2a_* for spline interpolation in magnetic surface coordinates
  use radial_module,only: psival_new,ra,dra,dpsidra,pfn,tfn,&
       & rho_normalized,pressure_normalized,pprime_normalized,qpsi,qprime,fpsi,ffprime, & !well known magnetic surface functions (in SI unit)
       & b0_axis,current,psi_axis,psi_lcfs,r_minor,r_axis,z_axis,eps !b0_axis is the strength of magnetic field at the magnetic axis, !(r_axis,z_axis) is the location of the magnetic axis, eps is the inverse aspect ratio, 
  use poloidal_harmonics,only: mhtot !mhtot is the total poloidal Fourier harmonics used to expand perturbations
  use toroidal_harmonics,only: nh !toroidal mode number

  implicit none
  integer,parameter:: nwmax=2001 !assumed maximum length of the arrays, the actual length will be determined by nx, nz, or, np_lcfs
  real(p_):: tmp_psi(nwmax,nwmax)
  real(p_):: tmp_fpsi(nwmax),tmp_qpsi(nwmax),tmp_press(nwmax),tmp_pprime(nwmax),tmp_ffprime(nwmax)
  real(p_):: rleft,zmid,xdim,zdim !specification of the rectangular compuational within which the value of psi is known
  real(p_):: x_lcfs0(nwmax),z_lcfs0(nwmax) !the x and z coordinates of the points on the LCFS, nwmax is the assumed maximum length of the two arrays, the actual length will be determined by np_lcfs
  integer:: np_lcfs ! the actual length of arrays x_lcfs and z_lcfs
  real(p_):: psi_func !the interpolating poloidal flux function of two variables x and z (constructed by 2d cubic spline interpolation)

  real(p_):: pressure(nflux),pprime(nflux),rho(nflux)

  real(p_):: global_shear(nflux), alpha(nflux)  !alpha is the normalized pressure gradient
  real(p_):: circumference(nflux),r_grad_psi_s(nflux) !circumference of the poloidal cross section of every magnetic surface
  real(p_):: rth(mpoloidal,nflux),zth(mpoloidal,nflux),rpsi(mpoloidal,nflux),zpsi(mpoloidal,nflux) !partial derivatives

  real(p_):: av_one_over_rsq(nflux),bp_sq_av(nflux),b_av(nflux),one_over_b_av(nflux),averaged_pressure
  real(p_):: dvdpsi(nflux),dvdpsi2(nflux),safety_factor(nflux)
  real(p_):: psi_gradient_fcg(mpoloidal,nflux),bsq(mpoloidal,nflux),bsq_av(nflux) !bsq is square of the strength of magnetic field
  real(p_):: b0_th(mpoloidal,nflux),b_toroidal(mpoloidal,nflux),bp_boundary_av
  real(p_):: psi_x_fcg(mpoloidal,nflux),psi_z_fcg(mpoloidal,nflux), psi_xx_fcg(mpoloidal,nflux),&  
       &  psi_zz_fcg(mpoloidal,nflux), psi_xz_fcg(mpoloidal,nflux),psi_zx_fcg(mpoloidal,nflux) !various quanties on the flux coordinates grids
  real(p_):: jacobian(mpoloidal,nflux),jacobian2(mpoloidal,nflux),jacobian_av(nflux) !,sign_gradient_psi
  real(p_):: sign_bth,sign_bphi
  real(p_):: kappas(mpoloidal,nflux),kappa_psi(mpoloidal,nflux),kappa_psi2(mpoloidal,nflux),sigma_mu0(mpoloidal,nflux)
  real(p_):: q_local1(mpoloidal,nflux),q_local2(mpoloidal,nflux),q_local3(mpoloidal,nflux)
  real(p_):: local_shear(mpoloidal,nflux),local_shear2(mpoloidal,nflux), local_shear0(mpoloidal,nflux)
  real(p_):: jacobian_th(mpoloidal,nflux),jacobian_th2(mpoloidal,nflux)
  real(p_):: psi_dot_theta(mpoloidal,nflux),psi_dot_theta2(mpoloidal,nflux),psi_dot_theta_psi(mpoloidal,nflux)

  !weight functions appearing in Fourier integral. The weight functions are defined in the document: /home/yj/theory/mhd/mhd.tm.
  real(p_):: w1(mpoloidal,nflux),w2(mpoloidal,nflux),w3(mpoloidal,nflux),w4(mpoloidal,nflux) 
  real(p_):: w5(mpoloidal,nflux),w6(mpoloidal,nflux),w7(mpoloidal,nflux),w8(mpoloidal,nflux)
  real(p_):: w9(mpoloidal,nflux),w10(mpoloidal,nflux),w11(mpoloidal,nflux),w12(mpoloidal,nflux)
  real(p_):: w13(mpoloidal,nflux),w14(mpoloidal,nflux),w15(mpoloidal,nflux),w16(mpoloidal,nflux)
  real(p_):: w17(mpoloidal,nflux),w18(mpoloidal,nflux),w19(mpoloidal,nflux),w20(mpoloidal,nflux)
  real(p_):: w21(mpoloidal,nflux),w22(mpoloidal,nflux),w23(mpoloidal,nflux),w24(mpoloidal,nflux)

  real(p_),dimension(:),allocatable:: psival_nx,pfn_nx,press_nx,pprime_nx,qpsi_nx,fpsi_nx,ffprime_nx,tf_nx,tfn_nx 
  real(p_),dimension(:),allocatable:: x_lcfs,z_lcfs
  real(p_):: theta(mpoloidal)

  real(p_):: psi_gradient_func,psi_x_func,psi_z_func,psi_xx_func,psi_zz_func,psi_xz_func,psi_zx_func !names of interpolating functions
  !real(p_):: rtbis0 !name of a bisection root finder
  real(p_):: wa0,va0  !wa0 is a characteristic Alfven frenquency, which is defined by wa0=VA0/R0, where VA0 is the Alfven speed at the magnetic axis
  real(p_):: total_tf
  real(p_):: dtheta,dpsi_nx,vol,sum,sum_tmp1,sum_tmp2,tmp_x

  logical:: recalculate_flux_coordinates,calculate_continua,calculate_global_mode,calculate_fast_ions_pressure, filter
  integer:: i,j,k
  real(p_):: tarray(2) !store the cpu clock time
  character(100):: poloidal_angle_type
  namelist/poloidal_angle/poloidal_angle_type
  namelist/control_parameters/nh,gamma,&
       & calculate_continua, recalculate_flux_coordinates,&
       & slow_sound_approximation,filter,calculate_global_mode,calculate_fast_ions_pressure !slow_sound_approximation_for_mode
  INTEGER :: count1,count2, count_rate, count_max

  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  call cpu_time(tarray(1))    !cpu_time is a f95 intrinsic subroutine

  open(11,file='gtaw.in')
  read(11,control_parameters)
  close(11)
  open(11,file='gtaw.in')
  read(11,poloidal_angle)
  close (11)

  write(*,control_parameters)
  write(*,poloidal_angle)
  write(*,*) 'nflux= ', nflux, 'mpoloidal= ', mpoloidal

  call read_gfile(nwmax,tmp_psi,nx,nz,rleft,zmid,xdim,zdim,psi_axis,psi_lcfs, &
       & r_axis,z_axis,x_lcfs0,z_lcfs0,np_lcfs,tmp_fpsi,tmp_qpsi,tmp_press,tmp_pprime,tmp_ffprime,current)

  write(*,*) '(R,Z) grids in G-file are ', 'nx=',nx,'nz=',nz
  !now the value of nx and nz is known, we allocate the arrays using the actual lenght of the corresponding array:
  !Note that nx in g-file is also used to define the number of radial grid points.
  allocate(psi(nx,nz))
  allocate(xarray(nx))
  allocate(zarray(nz))
  allocate(fpsi_nx(nx))
  allocate(qpsi_nx(nx))
  allocate(press_nx(nx))
  allocate(psival_nx(nx))
  allocate(pfn_nx(nx))
  allocate(pprime_nx(nx))
  allocate(ffprime_nx(nx))
  allocate(tf_nx(nx))
  allocate(tfn_nx(nx))

  !now we use the actual length of the array, instead of the assumed maximum length
  do i=1,nx
     do j=1,nz
        psi(i,j)=tmp_psi(i,j) !use the actual length of the psi array, instead of the assumed maximum length
     enddo
  enddo

  dpsi_nx=(psi_lcfs-psi_axis)/(nx-1)
  do i=1,nx 
     psival_nx(i)=psi_axis+dpsi_nx*(i-1) !uniform psi array, this is the radial coordinator used in G-file for the following magnetic surface function.
  enddo
  pfn_nx=(psival_nx-psival_nx(1))/(psival_nx(nx)-psival_nx(1))


  do i=1,nx !use the actual length of the radial array, instead of the assumed maximum length
     fpsi_nx(i)=tmp_fpsi(i) 
     qpsi_nx(i)=tmp_qpsi(i) 
     press_nx(i)=tmp_press(i) 
     pprime_nx(i)=tmp_pprime(i) 
     ffprime_nx(i)=tmp_ffprime(i) 
     !Actually, tmp_* can be used without the above transformation since they are one-dimension array.
  enddo

!!$  do j=1,nx !pf is the poloidal flux enclosed by a magnetic surface
!!$     pf_nx(j)=twopi*(psival_nx(j)-psi_axis) !if pf>0 then indicate btheta is along anti-clockwise direction
!!$  enddo

  sign_bphi=1._p_
  if(fpsi_nx(1)<0._p_) sign_bphi=-1.0_p_
  !  write(*,*)  'sign_bphi=', sign_bphi
  sign_bth=1._p_
  if(psi_lcfs<psi_axis) sign_bth=-1._p_ !the positive dirctio of theta is anti-clockwise viewed in the grad_phi direction
  !  write(*,*)  'sign_bth=', sign_bth
  qpsi_nx=abs(qpsi_nx)*sign_bth*sign_bphi !the safety-factor given in G-file is the absoulte value while the safety factor used in this code can be negative, which depends on the choice of the positive direction of poloidal angle and toroidal angle.

  tf_nx(1)=0._p_ !tf_nx is the toroidal magnetic flux
  do j=2,nx 
     tf_nx(j)=tf_nx(j-1)+(qpsi_nx(j)+qpsi_nx(j-1))/two*twopi*dpsi_nx  !using the formula dtf=q*dpf=q*d(pf_gs)*twopi
  enddo
  total_tf=tf_nx(nx)
  tfn_nx= tf_nx/tf_nx(nx) !normalized to unit

  do i=1,nx !construct the X array
     xarray(i)=rleft+xdim/(nx-1)*(i-1)
  enddo

  do j=1,nz !construct the Z array
     zarray(j)=(zmid-zdim/two)+zdim/(nz-1)*(j-1)
  enddo

  call write_gnuplot_contour_data2(nx,nz,xarray,zarray,psi,'psi.txt')
  write(*,*) 'psi_max=',maxval(psi),'psi_min=',minval(psi)

  allocate(y2a_psi(nx,nz)) !this is an intermedial array needed in cubic spline interpolation.
  call splie2(xarray,zarray,psi,nx,nz,y2a_psi) !after this call, the function psi_func (defined later in this file) is ready to be used.


  allocate(psi_x(nx,nz))
  allocate(psi_z(nx,nz))
  allocate(psi_xx(nx,nz))
  allocate(psi_zz(nx,nz))
  allocate(psi_xz(nx,nz))
  allocate(psi_zx(nx,nz))
  allocate(psi_gradient(nx,nz))
  call calculate_poloidal_flux_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
  allocate(y2a_gradient(nx,nz)) ! an array in spline interpolation to store 2nd derivatives
  allocate(y2a_psi_x(nx,nz))
  allocate(y2a_psi_z(nx,nz))
  allocate(y2a_psi_xx(nx,nz))
  allocate(y2a_psi_zz(nx,nz))
  allocate(y2a_psi_xz(nx,nz))
  allocate(y2a_psi_zx(nx,nz))
  call splie2(xarray,zarray,psi_gradient,nx,nz,y2a_gradient) !after this call, the function psi_gradient_func is ready to be used.
  call splie2(xarray,zarray,psi_x,nx,nz,y2a_psi_x) !after this call, the function psi_x_func is ready to be used.
  call splie2(xarray,zarray,psi_z,nx,nz,y2a_psi_z) !after this call, the function psi_z_func is ready to be used.
  call splie2(xarray,zarray,psi_xx,nx,nz,y2a_psi_xx) !after this call, the function psi_xx_func is ready to be used.
  call splie2(xarray,zarray,psi_zz,nx,nz,y2a_psi_zz) !after this call, the function psi_zz_func is ready to be used.
  call splie2(xarray,zarray,psi_xz,nx,nz,y2a_psi_xz) !after this call, the function psi_xz_func is ready to be used.
  call splie2(xarray,zarray,psi_zx,nx,nz,y2a_psi_zx) !after this call, the function psi_zx_func is ready to be used.


  !test whether the interpolating function psi_func is implemented correctly by comparing its results with the original data
!!$  open(111,file='psicontour')
!!$  do i=1,nx
!!$     do j=1,nz
!!$        write(111,*) xarray(i),zarray(j), psi_func(xarray(i),zarray(j))
!!$     enddo
!!$     write(111,*)
!!$  enddo
!!$  close(111)
!!$
!!$ open(111,file='psicontour2')
!!$  do i=1,nx
!!$     do j=1,nz
!!$        write(111,*) xarray(i),zarray(j), psi(i,j)
!!$     enddo
!!$     write(111,*)
!!$  enddo
!!$  close(111)

  allocate(x_lcfs(np_lcfs))
  allocate(z_lcfs(np_lcfs))
  do i=1,np_lcfs
     x_lcfs(i)=x_lcfs0(i)
     z_lcfs(i)=z_lcfs0(i)
  enddo

  !  write(*,*) 'minimum value of psi determined directly from the discrete psi data', minval(psi) !,  'min value location=',minloc(psi)
  write(*,*) 'value of psi at magnetic axis calculated from interpolating function',psi_func(r_axis,z_axis)
  write(*,*) 'value of psi at magnetic axis specified in g-file, psi_axis=',psi_axis

  call choose_radial_coordinate(psi_axis,psi_lcfs,psival_nx,tfn_nx,qpsi_nx,nx,total_tf,nflux,psival_new,ra,dra,dpsidra)

  do j=1,nflux !this is the normalized poloidal flux array corresponding to psival_new array
     pfn(j)=(psival_new(j)-psi_axis)/(psi_lcfs-psi_axis)
  enddo

  call interpolate_flux_function(nx,pfn_nx,psival_nx,fpsi_nx,qpsi_nx,press_nx,pprime_nx,ffprime_nx,tfn_nx, &
       & nflux,pfn,fpsi,qpsi,pressure,pprime,ffprime,qprime,tfn) !interpolate to psival_new grids

  call arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis) !this is to guarantee that (x_lcfs_new(1),z_lcfs_new(1)) is exactly on the low-field-side of the midplane

  call calculate_contours_wrapper(recalculate_flux_coordinates,poloidal_angle_type,mpoloidal,nflux,psival_new,&
       & fpsi,r_axis,z_axis,x_lcfs, z_lcfs,np_lcfs,r_new,z_new,&
       & av_one_over_rsq,dvdpsi2,safety_factor,circumference,r_grad_psi_s,bp_sq_av,b_av,one_over_b_av)
  call midplane() !verify that the values of poloidal angle on the low-field side of the midplane is zero, and calculate the major radius of the low-field-side point on every flux surface

  !call artifial_mode(r_new,z_new,ra,nflux,mpoloidal) !this is to construct a mode using analytical expressions, for test purpose
  !----test begin
  !  do j=2,nflux !use the relationship between dVdPsi and toroidal field function F and safety factor to calculate dVdPsi
  do j=1,nflux !use the relationship of VdPsi with toroidal field function F and safety factor to calculate dVdPsi
     !here V is the volume within a magnetic surface
     dvdpsi(j)= qpsi(j)*twopi**2/(fpsi(j)*av_one_over_rsq(j))
  enddo
  !to check that two methods of calculating dvdpsi give the same result. The agreement between the results indicates the constructed interpolating function for the gradient of Psi is correct.
  open(11,file='dvdpsi.txt')
  do j=1,nflux
     write(11,*) j,dvdpsi(j),dvdpsi2(j)
  enddo
  close(11)
  !  call calculate_fpsi_sunist(nflux,dvdpsi2,av_one_over_rsq,qpsi,fpsi)

  !----test end
  do i=1,mpoloidal !calculate the gradient of psi and its six partial derivatives on the flux coordinate grid points by using the intepolating function
     do j=1,nflux !here fcg stands for flux coordinate grids
        psi_gradient_fcg(i,j)=psi_gradient_func(r_new(i,j),z_new(i,j))
        psi_x_fcg(i,j)=psi_x_func(r_new(i,j),z_new(i,j)) 
        psi_z_fcg(i,j)=psi_z_func(r_new(i,j),z_new(i,j))
        psi_xx_fcg(i,j)=psi_xx_func(r_new(i,j),z_new(i,j))
        psi_zz_fcg(i,j)=psi_zz_func(r_new(i,j),z_new(i,j))
        psi_xz_fcg(i,j)=psi_xz_func(r_new(i,j),z_new(i,j)) 
        psi_zx_fcg(i,j)=psi_zx_func(r_new(i,j),z_new(i,j))
     enddo
  enddo
  call plot_poloidal(mpoloidal,nflux,psi_gradient_fcg,sqrt(psi_x_fcg**2+psi_z_fcg**2),'psi_gradient.txt')
  call plot_poloidal(mpoloidal,nflux,psi_gradient_fcg/r_new,sqrt(psi_x_fcg**2+psi_z_fcg**2),'bp.txt') !poloidal magnetic field
  call plot_poloidal(mpoloidal,nflux,psi_x_fcg,psi_z_fcg,'psi_x_z.txt')
  call plot_poloidal(mpoloidal,nflux,psi_xx_fcg,psi_zz_fcg,'psi_xx_zz.txt')
  call plot_poloidal(mpoloidal,nflux,psi_xz_fcg,psi_zx_fcg,'psi_xz_zx.txt') !the results indicate psi_xz_fcg and psi_zx_fcg are equal to each other

  ! sign_gradient_psi=one
  !  if((psi_lcfs-psi_axis).lt.0.)  sign_gradient_psi=-one !this quantity is used to determine the sign of Jacobian

  call calculate_jacobian(mpoloidal,nflux,r_new,z_new,ra,dpsidra,poloidal_angle_type,&
       & psi_gradient_fcg,circumference,r_grad_psi_s,jacobian,jacobian_th)

  !--test the accuray of jacobian for a analytical coordinate transformation
!!$  do i=1,mpoloidal 
!!$     theta(i)=0.0_p_+twopi/(mpoloidal-1)*(i-1)
!!$  enddo
!!$  do i=1,mpoloidal
!!$     do j=1,nflux
!!$        r_new(i,j)=2._p_+ra(j)*cos(theta(i))
!!$        z_new(i,j)=ra(j)*sin(theta(i))
!!$     enddo
!!$  enddo
!!$
!!$  !--test end
!!$
!!$  do i=1,mpoloidal
!!$     do j=1,nflux
!!$jacobian(i,j)=-r_new(i,j)*ra(j)
!!$enddo
!!$enddo

  !Verify the correctness of the jacobian by directly calculating the numerical jacobian
  !and calculate a new quantity psi_dot_theta, the scalar product between grad_psi and grad_theta)/grad_psi^2
  dtheta=twopi/(mpoloidal-1) !grid interval, uniform theta grid is assumed
  call calculate_jacobian_directly(mpoloidal,nflux,r_new,z_new,dtheta,dra,jacobian2, &
       & psi_dot_theta2,rpsi, rth, zpsi,zth,jacobian_th2) 

  call plot_poloidal(mpoloidal,nflux,jacobian,jacobian2,'jacobian.txt')

  if (sign(one,jacobian(mpoloidal/2,nflux/2)) .ne. sign(one,jacobian2(mpoloidal/2,nflux/2))) then
     stop 'The signs of Jacobian1 and Jacobian2 are different, please check!'
  endif

  do i=1,mpoloidal
     do j=1,nflux
        q_local1(i,j)=-jacobian(i,j)/r_new(i,j)**2*fpsi(j)/dpsidra(j)
        q_local2(i,j)=-jacobian2(i,j)/r_new(i,j)**2*fpsi(j)/dpsidra(j)
        q_local3(i,j)=qpsi(j)
     enddo
  enddo

  call plot_poloidal(mpoloidal,nflux,q_local1,q_local3,'q_local.txt')

  !----------------
  call plot_poloidal(mpoloidal,nflux,rpsi,zpsi,'rpsi_zpsi.txt')
  call plot_radial(mpoloidal,nflux,rpsi,zpsi,'rpsi_zpsi_radial.txt')
  call plot_radial(mpoloidal,nflux,r_new,z_new,'r_z_radial.txt')
  call safety_factor2(mpoloidal,nflux,r_new,z_new,jacobian,fpsi,dpsidra) !test, calculating q using jacobian
  call calculate_psi_dot_theta(mpoloidal,nflux,rpsi,rth,zpsi,zth,jacobian,r_new,psi_dot_theta,psi_dot_theta_psi) !psi_dot_theta is the scalar product between grad_psi and grad_theta
  !call plot_poloidal(mpoloidal,nflux,-(zth*zpsi+rth*rpsi),psi_dot_theta,'psi_dot_theta.txt')
  call plot_poloidal(mpoloidal,nflux,psi_dot_theta,psi_dot_theta2,'psi_dot_theta.txt')

  !safety_factor is calculated numerically in "calculate_contours" subroutine, but a sign factor is dropped, now take back this:
  safety_factor=safety_factor*sign(one,jacobian(mpoloidal/2,nflux/2))*sign(one,dpsidra(nflux/2))

  !compare the safety_factor() with that given in G-file qpsi() (note that the safety-factor given in G-file is the absoulte value while the safety factor calculated in this code can be negative, which depends on the choice of the positive direction of poloidal angle and toroidal angle).

  open(1234,file='profile.txt') 
  do j=1,nflux
     write(1234,*) ra(j),sqrt(pfn(j)),sqrt(tfn(j)),qpsi(j),safety_factor(j),pressure(j),dpsidra(j) &
          &  , fpsi(j),ffprime(j)/fpsi(j),qprime(j),qprime(j)*dpsidra(j),pprime_normalized(j)
     !     write(*,*) j, sqrt(pfn(j)),sqrt(tfn(j)),qpsi(j),qprime(j), qprime(j)*dpsidra(j)
  enddo
  close(1234)

  
  if(safety_factor(nflux/2)>0) then
     !q>0 for left-hand helix
     write(*,*) '----------safety_factor is positive---------------------'
  else
     write(*,*) '----------safety_factor is negative---------------------'
  endif
  if(sign(one,safety_factor(nflux/2))*qpsi(1)<0) stop '****the sign of safety factor needs checking********'
  !  qpsi=abs(qpsi)*sign(one,safety_factor(nflux/2)) !this is to ensure qpsi takes correct sign (the qpsi given in G-file is usually the absolute value)
!!$  do j=13,nflux-2
!!$     qpsi(j)=safety_factor(j)
!!$  enddo
  !  write(*,*) 'q(0) in gfile',qpsi_nx(1),'q(0) from interpolating',safety_factor(1)
  write(*,*) 'q(0) in gfile',qpsi_nx(1),'q(0) from interpolating',&
       safety_factor(1)+(safety_factor(2)-safety_factor(1))/(psival_new(2)-psival_new(1))*(psi_axis-psival_new(1))
  !write(*,*) 'q at LCFS from gfile',qpsi_nx(nx),'q at LCFS from interpolating',safety_factor(nflux)
  write(*,*) 'q at LCFS from gfile',qpsi_nx(nx),'q at LCFS from interpolating',safety_factor(nflux)+&
       & (safety_factor(nflux)-safety_factor(nflux-1))/(psival_new(nflux)-psival_new(nflux-1))*(psi_lcfs-psival_new(nflux))
  if(minval(abs(qpsi_nx))<one)  call find_q_equal_one_surface(tfn_nx,psival_nx,pfn_nx,qpsi_nx,nx,x_lcfs,z_lcfs,&
       & np_lcfs,r_axis,z_axis)
  call calculate_bsq(mpoloidal,nflux,r_new,z_new,psi_gradient_fcg,fpsi,bsq,b0_th,b_toroidal,bp_boundary_av)
  b0_axis=abs(fpsi(1)/r_axis) !the strength of magnetic field at magnetic axis
  pprime_normalized=pprime*two*mu0/b0_axis**2 !the gradient of the thermal pressure normalized to magnetic pressure at magnetic axis
  call plot_poloidal(mpoloidal,nflux,b_toroidal,b_toroidal,'b_toroidal.txt')

  call normalized_parallel_current_density(mpoloidal,nflux,fpsi,ffprime,b0_axis,pprime_normalized,bsq,sigma_mu0) !defined by sigma_mu0=mu0*J_dot_B/B^2
  write(*,*) 'toroidal current given in G-file is (kA)', current/(1d3)
  call toroidal_current(mpoloidal,nflux,dra,dtheta,jacobian,ffprime,pprime,r_new,current)
  call toroidal_current2(mpoloidal,nflux,dra,dtheta,jacobian,psi_x_fcg,psi_xx_fcg,psi_zz_fcg,r_new,current)
  !call force_analysis(mpoloidal,nflux,psi_x_fcg,psi_xx_fcg,psi_zz_fcg,r_new,z_new,ffprime,pprime)
  call check_directions(current,psi_axis,psi_lcfs,fpsi_nx(1),fpsi_nx(nx))

  call plot_poloidal(mpoloidal,nflux,sigma_mu0,sigma_mu0/jacobian,'sigma_mu0.txt')

  open(11,file='av_j_dot_b.txt') 
  do j=2,nflux-1
     bsq_av(j)=bp_sq_av(j)+av_one_over_rsq(j)*fpsi(j)**2
     write(11,*)  pfn(j), mu0*fpsi(j)*pprime(j), ffprime(j)/fpsi(j)*bsq_av(j), &
          &  (mu0*fpsi(j)*pprime(j)+ffprime(j)/fpsi(j)*bsq_av(j)), & !give mu0*<J.B>, which has two terms, one of which is proportional to pprime, one to fprime
          &  (mu0*fpsi(j)*pprime(j)*one_over_b_av(j)+ffprime(j)/fpsi(j)*b_av(j)) !gvie mu0<J.B/B>
  enddo
  close(11)

  call plot_poloidal(mpoloidal,nflux,jacobian_th,jacobian_th2,'jacobian_th.txt')
  call plot_poloidal(mpoloidal,nflux,r_new,z_new,'r_z.txt')
  !call plot_poloidal2(mpoloidal,nflux,rth,zth,rpsi,zpsi,r_new,'r_z_p.txt')
  call plot_poloidal(mpoloidal,nflux,jacobian**2*psi_gradient_fcg**2/r_new**2 &
       &,zth*zpsi+rth*rpsi,'jacobian_zth_rth.txt')

  call curvature(mpoloidal,nflux,r_new,z_new,fpsi,dpsidra,jacobian,jacobian_th,rth,zth,bsq,b0_th,circumference,kappas,kappa_psi)
  !  call plot_poloidal(mpoloidal,nflux,kappa_psi,kappa_psi,'normal_curvature.txt')

  if (filter.eqv. .true.) then
     call filter_theta(mpoloidal,nflux,kappa_psi,kappa_psi2) !filter out some Fourier components to generate up-down asymmetric kappas
     call plot_poloidal(mpoloidal,nflux,kappa_psi,kappa_psi2,'normal_curvature.txt')
     kappa_psi=kappa_psi2 !use the filtered data
  endif


  call local_magnetic_shear(mpoloidal,nflux,r_new,fpsi,ffprime,psi_x_fcg,psi_z_fcg,&
       & psi_xx_fcg,psi_zz_fcg,psi_xz_fcg,local_shear) !calculate the local shear in the cylindrical coordinates
  call local_magnetic_shear2(mpoloidal,nflux,dra,fpsi,jacobian,r_new,psi_dot_theta_psi,dpsidra,local_shear2) !calculate the local shear in the flux coordinates

  call plot_poloidal(mpoloidal,nflux,local_shear,local_shear2,'local_shear.txt')
  call plot_radial  (mpoloidal,nflux,local_shear,local_shear2,'local_shear_radial.txt')
  call write_gnuplot_contour_data(mpoloidal,nflux,r_new,z_new,local_shear,'shear_contour.txt')
  call benchmark_local_gloabal_shear(mpoloidal,nflux,jacobian,qprime,dpsidra,local_shear,'global_shear.txt')
  call benchmark_local_gloabal_shear(mpoloidal,nflux,jacobian,qprime,dpsidra,local_shear2,'global_shear2.txt')

  call local_magnetic_shear0(mpoloidal,nflux,r_new,fpsi,ffprime,psi_x_fcg,psi_z_fcg,&
       & psi_xx_fcg,psi_zz_fcg,psi_xz_fcg,local_shear0) !local_shear0 is local_shear*psi_gradient_sq
  call plot_poloidal(mpoloidal,nflux,local_shear*psi_gradient_fcg**2,local_shear0,'s_grad_psi_sq.txt')

  write(*,*) 'Magnetic field (Telsa) at magentic axis, B_phi (can be negative), is', fpsi(1)/r_axis !,b0_axis
  call vacuum_toroidal_magnetic_field(r_axis)
  write(*,*) 'Pressure at the magnetic axis is (kPa)', pressure(1)/1000._p_
  pressure_normalized=pressure/(b0_axis**2/(two*mu0)) !thermal pressure normalized to the magnetic pressure at magnetic axis, i.e. beta
  write(*,*) 'Beta value at magnetic_axis is ', pressure_normalized(1)
  call calculate_volume_averaged_beta(nflux,mpoloidal,pressure_normalized,jacobian,dra,dtheta,averaged_pressure)
  call calculate_surface_averaged_beta(nflux,mpoloidal,pressure_normalized,jacobian,r_new,dra,dtheta)
  write(*,*) 'Normalized beta (betaN)=',averaged_pressure*10**8*(r_minor*b0_axis)/abs(current)
  write(*,*) 'Poloidal beta= ',averaged_pressure*(b0_axis/bp_boundary_av)**2
  call plasma_volume(nflux,mpoloidal,jacobian,dra,dtheta,dvdpsi,psival_new,vol) !plasma volume within LCFS
  write(*,*) 'Store energy is (kJ)',three/two*averaged_pressure*b0_axis**2/(two*mu0)*vol/1000._p_
  write(*,*) '--------finish setting equilibrium-------'

  !check whether the magnetic shear at the tip of the gap is weak enough for upper TAEs (TAEs with frequency near the upper tip of the gapxsa) to exist:
  !eps=0.45/1.9
  !eps=(maxval(x_lcfs)-minval(x_lcfs))/2._p_/r_axis
  !write(*,*) 'eps=',eps
  open(11,file='s_alpha.txt')
  do j=1,nflux
     global_shear(j)=one/(two*qpsi(j))*qprime(j)*(psi_lcfs-psi_axis)
     alpha(j)=-one/eps*qpsi(j)**2*pprime_normalized(j)/(two*sqrt(pfn(j)))*(psi_lcfs-psi_axis)
     write(11,*) j, sqrt(pfn(j)),global_shear(j),alpha(j),eps*sqrt(pfn(j))
  enddo
  close(11)

  !  j=31 !the grid cooresponding to the radial location of the gap tip
!!$  j=116
!!$  write(*,*) 'shear=', global_shear(j)
!!$  write(*,*) 'alpha=', alpha(j)
!!$  write(*,*) alpha(j), -global_shear(j)**2+eps

  write(*,*) 'The Greenwalt density for this equilibrium is (10^20m^-3):',abs(current)/(pi*r_minor**2)/10**6
  !  call set_density(nflux,psival_new,psi_axis,psi_lcfs,rho) !set radial mass density profile
  call set_density(nflux,pfn,psi_axis,psi_lcfs,rho) !set radial mass density profile
  rho_normalized=rho/rho(1) !mass density normalized by its value at magentic axis
  va0=b0_axis/sqrt(mu0*rho(1)) !Alfven speed at the magnetic axis
  wa0=va0/r_axis !wa0 is the characteristic Alfven frequency, where VA0 is the Alfven speed at the magnetic axis
  write(*,*) 'Alfven velocity at magnetic_axis is (10^6m/s):  ', va0/10**6
  write(*,*) 'Sound speed at the magnetic axis is (10^6m/s): ',  sqrt(gamma*pressure(1)/rho(1))/10**6
  write(*,*) 'Characteristic Alfven frequency, Va0/r_axis, is (KHZ):',wa0/1000._p_
  !  call set_rotation(rotation_file_name,prof_grid_type,nflux,psival_new,psi_axis,psi_lcfs,rotation_new)
  !  call fast_ions_pressure(nflux,mpoloidal,pfn,sqrt(bsq),b0_axis)
  !  call fast_ions_pressure(psival_nx,fpsi_nx,nx,b0_axis,psi_axis,psi_lcfs,va0) !for mega code
  if(calculate_fast_ions_pressure.eqv..true.)  call fast_ions_pressure(pfn,fpsi,nflux,mpoloidal,&
       & r_new,z_new,b0_axis,psi_axis,psi_lcfs,va0) !for mega code

!   call provide_equ_to_GEM_code(nflux,mpoloidal,r_new,z_new,ra,qpsi,fpsi)
  
  if ((calculate_continua .eqv. .false.) .and. (calculate_global_mode .eqv. .false.)) goto 111 !end the program
  call cylindrical_continua(mpoloidal,nflux,pfn,tfn,qpsi,rho,pressure,jacobian,bsq,dpsidra)  !obtain the corresponding cylindrical continuous spectrum, which will be compared with continua in toroidal geometry

  call calculate_weight_functions(mpoloidal,nflux,r_new,circumference,jacobian, &
       & psi_gradient_fcg,psi_x_fcg,psi_z_fcg,psi_xx_fcg,psi_zz_fcg,psi_xz_fcg, &
       & psi_dot_theta, bsq,kappas,kappa_psi,sigma_mu0,local_shear0, &
       & fpsi,ffprime,qpsi,qprime,dra,dpsidra, &
       & w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,w23,w24) !local_shear0 is local_shear*psi_gradient_sq

  call plot_poloidal24(mpoloidal,nflux,w1,w2,w3,w4,w5,w6,w7,w8,w9,&
       & w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,w23,w24,'weight_poloidal.txt') !all these data can be examined in one page by using gnuplot 

  call calculate_kernel(mpoloidal,nflux,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,w23,w24)

  if(calculate_continua .eqv. .true.) call calculate_continous_spectrum(nflux,nh,r_axis,b0_axis,wa0)

  if(calculate_global_mode .eqv. .true.) then
     call calculate_radial_matrix(mpoloidal,nflux,mhtot)
     call find_global_mode(nh,mhtot,nflux,wa0)
  endif

111 call cpu_time(tarray(2))
  write (*,*) 'CPU time used (seconds)', tarray(2)-tarray(1)
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  write(*,*) 'Wall time used (seconds) ', (count2-count1)/count_rate
end program GTAW

!!$function r_func(theta0,psi0)
!!$  use precision,only:p_
!!$  use flux_grids,only: theta,r_new,y2a_r_new,mpoloidal,nflux
!!$  use radial_module,only:pfn
!!$  implicit none
!!$  real(p_):: r_func,theta0,psi0,rval
!!$  call splin2(theta,pfn,r_new,y2a_r_new,mpoloidal,nflux,theta0,psi0,rval)
!!$  r_func=rval
!!$end function 
!!$
!!$function z_func(theta0,psi0)
!!$  use precision,only:p_
!!$  use flux_grids,only: theta,z_new,y2a_z_new,mpoloidal,nflux
!!$  use radial_module,only:pfn
!!$  implicit none
!!$  real(p_):: z_func,theta0,psi0,zval
!!$  call splin2(theta,pfn,z_new,y2a_z_new,mpoloidal,nflux,theta0,psi0,zval)
!!$  z_func=zval
!!$end function 


subroutine calculate_fpsi_sunist(nflux,dvdpsi,av_one_over_rsq,qpsi,fpsi)
  !for the study of AEs on sunist tokamak 2015-11-29, to callculate fpsi consistent with the given q profile
 use precision,only:p_
  use constants,only:twopi
  implicit none
  integer,intent(in):: nflux
  real(p_),intent(in):: dvdpsi(nflux),av_one_over_rsq(nflux),qpsi(nflux),fpsi(nflux)
  real(p_)::fpsi_new(nflux)
  integer:: j
  open(52,file='poloidal_current_function.txt')
  do j=1,nflux
     fpsi_new(j)=qpsi(j)*twopi**2/(dvdpsi(j)*av_one_over_rsq(j))
     write(52,*) j,fpsi_new(j),fpsi(j)
  enddo
  close(52)
!     dvdpsi(j)= qpsi(j)*twopi**2/(fpsi(j)*av_one_over_rsq(j))
end subroutine calculate_fpsi_sunist
