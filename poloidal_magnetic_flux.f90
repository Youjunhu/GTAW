

function psi_func(xval,zval)
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi,y2a_psi,nx,nz
  implicit none
  real(p_):: psi_func
  real(p_)::xval,zval,psival
  call splin2(xarray,zarray,psi,y2a_psi,nx,nz,xval,zval,psival)

  psi_func=psival
end function psi_func


function psi_gradient_func(xval,zval) result (z)
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi_gradient,y2a_gradient,nx,nz
  implicit none
  real(p_):: z,xval,zval
  call splin2(xarray,zarray,psi_gradient,y2a_gradient,nx,nz,xval,zval,z)
end function psi_gradient_func

function psi_x_func(xval,zval) result (z)
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi_x,y2a_psi_x,nx,nz
  implicit none
  real(p_):: z,xval,zval
  call splin2(xarray,zarray,psi_x,y2a_psi_x,nx,nz,xval,zval,z)
end function 

function psi_z_func(xval,zval) result (z)
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi_z,y2a_psi_z,nx,nz
  implicit none
  real(p_):: z,xval,zval
  call splin2(xarray,zarray,psi_z,y2a_psi_z,nx,nz,xval,zval,z)
end function 

function psi_xx_func(xval,zval) result (z)
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi_xx,y2a_psi_xx,nx,nz
  implicit none
  real(p_):: z,xval,zval
  call splin2(xarray,zarray,psi_xx,y2a_psi_xx,nx,nz,xval,zval,z)
end function 

function psi_zz_func(xval,zval) result (z)
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi_zz,y2a_psi_zz,nx,nz
  implicit none
  real(p_):: z,xval,zval
  call splin2(xarray,zarray,psi_zz,y2a_psi_zz,nx,nz,xval,zval,z)
end function 


function psi_xz_func(xval,zval) result (z)
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi_xz,y2a_psi_xz,nx,nz
  implicit none
  real(p_):: z,xval,zval
  call splin2(xarray,zarray,psi_xz,y2a_psi_xz,nx,nz,xval,zval,z)
end function 

function psi_zx_func(xval,zval) result (z)
  use precision,only:p_
  use poloidal_flux_2d,only: xarray,zarray,psi_zx,y2a_psi_zx,nx,nz
  implicit none
  real(p_):: z,xval,zval
  call splin2(xarray,zarray,psi_zx,y2a_psi_zx,nx,nz,xval,zval,z)
end function 




