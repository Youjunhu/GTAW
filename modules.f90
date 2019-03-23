module  precision
  implicit none
  !integer,parameter:: p_=kind(1.0)
  integer,parameter:: p_=kind(1.0d0)
end module precision

module  poloidal_flux_2d
!poloidal flux and its partial derivatives on (R,Z) plane, where (R,fai,Z) is the cylindrical coordinate system
  use precision,only:p_
  implicit none
  integer:: nx,nz !nx,nz are respectively the numbers of grids in R and Z directions (specified in G-file)
  real(p_),dimension(:),allocatable ::xarray,zarray ! R and Z array
  real(p_),dimension(:,:),allocatable ::psi,y2a_psi !psi, 2D array on (R,Z) plane, is the poloidal flux function appearing in Grad-Shafranov equation. 
  !psi is related to poloidal magnetic field by Bp=\nabal{psi}\times\nabla{fai},where fai is the usual toroidal angle;  y2a_psi is the 2nd order erivative array used in the 2D spline interpolation of psi
  real(p_),dimension(:,:),allocatable ::psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx !psi_x is the partial derivative with respect to x, similar meaning for others
  real(p_),dimension(:,:),allocatable ::psi_gradient,y2a_gradient !strength of the gradient of the poloidal flux, y2a_gradient is the second derivative array used in 2D spline interpolation
  real(p_),dimension(:,:),allocatable ::y2a_psi_x,y2a_psi_z,y2a_psi_xx,y2a_psi_zz,y2a_psi_xz,y2a_psi_zx ! y2a_* is the second derivative array used in 2D spline interpolation

end module poloidal_flux_2d

module flux_grids
  use precision,only:p_
  implicit none
  !integer,parameter:: nflux=66*2.5+30 !-39 !radial grid number
!   integer,parameter:: nflux=250 
!  integer,parameter:: nflux=180
!  integer,parameter:: nflux=98 !for mega
!  integer,parameter:: nflux=129
  integer,parameter:: nflux=129
  !integer,parameter:: nflux=66*2.5+50
!  integer,parameter:: nflux=256
  integer,parameter:: mpoloidal=129 
!  integer,parameter:: mpoloidal=2**7+1 !usually enough
!  integer,parameter:: mpoloidal=2**8+1 !high resolution, used for DIIID case (high q, q95=10.5)
  !integer,parameter:: mpoloidal=2**7+3 !equilibrium poloidal grid points number
  !integer,parameter:: mpoloidal=2**7-2**4+1 !equilibrium poloidal grid points number
  real(p_):: r_new(mpoloidal,nflux),z_new(mpoloidal,nflux) !store the magnetic surfaces found in subroutine calculate_contour()
  real(p_):: delta_q(mpoloidal,nflux) !delta_q is the toroidal angle shift used in defining the generalized toroidal angle
  real(p_):: midplane_theta(nflux) !the values of poloidal angle on the low-field side of the midplane
  real(p_):: x0(nflux) !values of R of magnetic surfaces on low-field-side of the midplane 
 ! real(p_):: theta(mpoloidal)
!  real(p_):: y2a_r_new(mpoloidal,nflux),y2a_z_new(mpoloidal,nflux) !for spline interpolation in magnetic surface coordinates
end module flux_grids

module radial_module
  use precision,only:p_
  use flux_grids,only: nflux
  implicit none
  real(p_):: psival_new(nflux),pfn(nflux),tfn(nflux) !psival_new is Aphi*R, tfn is the normalized toroidal flux, pfn is the normalized poloidal flux.
  real(p_):: ra(nflux) !ra is the radial coordinate, 
  real(p_):: dpsidra(nflux) !d(Aphi*R)/dra, the derivative of Aphi*R (GS poloidal flux) with respect to the radial coordinate ra
  real(p_):: dra !dra is the radial grid interval
!  real(p_):: sqrt_toroidal_flux2(nflux),poloidal_flux2(nflux) !radial grids corresponding to uniform sqrt(toroidal_flux)
  real(p_):: rho_normalized(nflux) !normalized mass density
  real(p_):: pressure_normalized(nflux)
  real(p_):: pprime_normalized(nflux)
  real(p_):: qpsi(nflux),qprime(nflux),fpsi(nflux),ffprime(nflux) !well known flux functions
  real(p_):: r_axis,z_axis   !location of magnetic axis
  real(p_):: r_major,r_minor,eps !eps is the inverse aspect ratio, alpha is the normalized pressure gradient
  real(p_):: b0_axis !magnetic field strength at magnetic axis
  real(p_):: current !the total toroidal current within the LCFS
  real(p_):: psi_axis,psi_lcfs
  integer:: starting_surface_number,ending_surface_number
end module radial_module

!!$module radial_interpolation
!!$  use precision,only:p_
!!$ implicit none
!!$
!!$end module radial_interpolation

!---------------perturbation relavant variables---------------------------------
module poloidal_harmonics
  implicit none
  !integer,parameter:: mh_upp=40 !mh_upp is the upper mode number
!  integer,parameter:: mh_upp=15 !mh_upp is the upper mode number
  integer,parameter:: mh_upp=30 !mh_upp is the upper mode number
!  integer,parameter:: mh_upp=50 !mh_upp is the upper mode number, 50 is used for DIIID case (high q, q95=10.5)
!  integer,parameter:: mh_upp=27 !mh_upp is the upper mode number
!  integer,parameter:: mh_low=-40 !mh_low is the lower mode number (can be a negative integer)
!  integer,parameter:: mh_low=-8 !mh_low is the lower mode number (can be a negative integer)
!  integer,parameter:: mh_low=-4 !mh_low is the lower mode number (can be a negative integer)
  integer,parameter:: mh_low=-30 !mh_low is the lower mode number (can be a negative integer)
!integer,parameter:: mh_low=-mh_upp !
! integer,parameter:: mh_low=-3 !mh_low is the lower mode number (can be a negative integer) 
  integer,parameter:: mhtot=mh_upp-mh_low+1 !the total number of the poloidal harmonics
end module poloidal_harmonics

module toroidal_harmonics
  implicit none
  integer:: nh
end module toroidal_harmonics


module kernel_module
  use precision,only:p_
  use flux_grids,only:nflux
  use poloidal_harmonics,only: mhtot
  implicit none

  complex(p_):: full_kernel0(nflux,mhtot,mhtot),full_kernel1(nflux,mhtot,mhtot),full_kernel2(nflux,mhtot,mhtot)
  complex(p_):: full_kernel3(nflux,mhtot,mhtot),full_kernel4(nflux,mhtot,mhtot),full_kernel5(nflux,mhtot,mhtot)
  complex(p_):: full_kernel6(nflux,mhtot,mhtot),full_kernel7(nflux,mhtot,mhtot),full_kernel8(nflux,mhtot,mhtot)
  complex(p_):: full_kernel9(nflux,mhtot,mhtot)
  complex(p_):: full_kernel10(nflux,mhtot,mhtot),full_kernel11(nflux,mhtot,mhtot),full_kernel12(nflux,mhtot,mhtot)
  complex(p_):: full_kernel13(nflux,mhtot,mhtot),full_kernel14(nflux,mhtot,mhtot),full_kernel15(nflux,mhtot,mhtot)
  complex(p_):: full_kernel16(nflux,mhtot,mhtot),full_kernel17(nflux,mhtot,mhtot),full_kernel18(nflux,mhtot,mhtot)
  complex(p_):: full_kernel19(nflux,mhtot,mhtot),full_kernel20(nflux,mhtot,mhtot),full_kernel21(nflux,mhtot,mhtot)
  complex(p_):: full_kernel22(nflux,mhtot,mhtot),full_kernel23(nflux,mhtot,mhtot),full_kernel24(nflux,mhtot,mhtot)

  complex(p_):: full_radial_matrix(nflux,2*mhtot,2*mhtot)

end module kernel_module


module constants
  use precision,only:p_
  implicit none
  logical:: slow_sound_approximation
  logical:: slow_sound_approximation_for_mode
  real(p_):: gamma !polytrope index
  real(p_),parameter:: pi=3.1415926_p_
  real(p_),parameter:: twopi=pi*2.0_p_
  real(p_),parameter:: fourpi=pi*4.0_p_
  real(p_),parameter:: mu0=fourpi*1.0d-7 !permeability in SI unit
  real(p_),parameter:: zero=0.0_p_
  real(p_),parameter:: one=1.0_p_
  real(p_),parameter:: two=2.0_p_
  real(p_),parameter:: three=3.0_p_
  real(p_),parameter:: four=4.0_p_
  real(p_),parameter:: five=5.0_p_
  real(p_),parameter:: six=6.0_p_
  real(p_),parameter:: seven=7.0_p_
  real(p_),parameter:: eight=8.0_p_
  real(p_),parameter:: nine=9.0_p_
  real(p_),parameter:: ten=10.0_p_
  real(p_),parameter:: eleven=11.0_p_
  real(p_),parameter:: twelve=12.0_p_
  real(p_),parameter:: thirteen=13.0_p_
  real(p_),parameter:: fourteen=14.0_p_
  real(p_),parameter:: fifteen=15.0_p_

  real(p_),parameter:: one_half=0.5_p_
  real(p_),parameter:: one_third=one/three
  real(p_),parameter:: one_fifth=0.2_p_
  real(p_),parameter:: three_halfs=1.5_p_
  complex(p_),parameter:: ii=(0._p_,1.0_p_)
end module constants

module fast_ions
  use precision,only:p_
  implicit none
!  real(p_),parameter:: md=3.3452d-27 !particle mass in kg (mass=2._p_*1.6726d-27 is for Deuterium, mass=9.1094d-31 is for electron)
  real(p_),parameter:: md=3.343d-27 !particle mass in kg (mass=2._p_*1.6726d-27 is for Deuterium, mass=9.1094d-31 is for electron)
  real(p_),parameter:: me=9.1094d-31 !me is electron mass
  real(p_),parameter:: kev=1.6022d-16   !unit J
  real(p_),parameter:: charge=1.6022d-19 
end module fast_ions
