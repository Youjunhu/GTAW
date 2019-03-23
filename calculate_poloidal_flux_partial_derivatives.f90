subroutine calculate_poloidal_flux_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z,&
     & psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
  !subroutine calculate_poloidal_flux_gradient(nx,nz,xarray,zarray,psi,psi_gradient)
  use precision,only:p_
  use constants,only:one,two,twopi
  implicit none
  integer,intent(in)::nx,nz
  real(p_),intent(in):: psi(nx,nz)
  real(p_),intent(in):: xarray(nx),zarray(nz)
  real(p_),intent(out):: psi_x(nx,nz),psi_z(nx,nz),psi_xx(nx,nz),psi_zz(nx,nz),psi_xz(nx,nz),psi_zx(nx,nz)
  real(p_),intent(out):: psi_gradient(nx,nz)

  integer:: i,j,i1,i2,j1,j2

  !first-order partial derivatives
  do i=1,nx
     do j=1,nz
        i2=i+1
        i1=i-1
        j2=j+1
        j1=j-1
        if(i.eq.1) i1=i
        if(j.eq.1) j1=j
        if(i.eq.nx) i2=i
        if(j.eq.nz) j2=j
        psi_x(i,j)=(psi(i2,j)-psi(i1,j))/(xarray(i2)-xarray(i1))
        psi_z(i,j)=(psi(i,j2)-psi(i,j1))/(zarray(j2)-zarray(j1))
     enddo
  enddo

  !second-order partial derivatives
  do i=1,nx
     do j=1,nz
        i2=i+1
        i1=i-1
        j2=j+1
        j1=j-1
        if(i.eq.1) i1=i
        if(j.eq.1) j1=j
        if(i.eq.nx) i2=i
        if(j.eq.nz) j2=j
        psi_xx(i,j)=(psi_x(i2,j)-psi_x(i1,j))/(xarray(i2)-xarray(i1))
        psi_zz(i,j)=(psi_z(i,j2)-psi_z(i,j1))/(zarray(j2)-zarray(j1))
        psi_xz(i,j)=(psi_x(i,j2)-psi_x(i,j1))/(zarray(j2)-zarray(j1))
        psi_zx(i,j)=(psi_z(i2,j)-psi_z(i1,j))/(xarray(i2)-xarray(i1))
     enddo
  enddo

!!$  !bottom left corner
!!$  i=1
!!$  j=1
!!$  psi_x(i,j)=(psi(i+1,j)-psi(i,j))/(xarray(i+1)-xarray(i))
!!$  psi_z(i,j)=(psi(i,j+1)-psi(i,j))/(zarray(j+1)-zarray(j))
!!$
!!$  !bottom right corner  
!!$  i=nx
!!$  j=1
!!$  psi_x(i,j)=(psi(i,j)-psi(i-1,j))/(xarray(i)-xarray(i-1))
!!$  psi_z(i,j)=(psi(i,j+1)-psi(i,j))/(zarray(j+1)-zarray(j))
!!$
!!$  !upper left corner
!!$  i=1
!!$  j=nz
!!$  psi_x(i,j)=(psi(i+1,j)-psi(i,j))/(xarray(i+1)-xarray(i))
!!$  psi_z(i,j)=(psi(i,j)-psi(i,j-1))/(zarray(j)-zarray(j-1))
!!$
!!$  !upper right corner
!!$  i=nx
!!$  j=nz
!!$  psi_x(i,j)=(psi(i,j)-psi(i-1,j))/(xarray(i)-xarray(i-1))
!!$  psi_z(i,j)=(psi(i,j)-psi(i,j-1))/(zarray(j)-zarray(j-1))
!!$
!!$  i=1
!!$  do j=2,nz-1
!!$     psi_x(i,j)=(psi(i+1,j)-psi(i,j))/(xarray(i+1)-xarray(i))
!!$     psi_z(i,j)=(psi(i,j+1)-psi(i,j-1))/(zarray(j+1)-zarray(j-1))
!!$  enddo
!!$
!!$  i=nx
!!$  do j=2,nz-1
!!$     psi_x(i,j)=(psi(i,j)-psi(i-1,j))/(xarray(i)-xarray(i-1))
!!$     psi_z(i,j)=(psi(i,j+1)-psi(i,j-1))/(zarray(j+1)-zarray(j-1))
!!$  enddo
!!$
!!$
!!$  j=1
!!$  do i=2,nx-1
!!$     psi_x(i,j)=(psi(i+1,j)-psi(i-1,j))/(xarray(i+1)-xarray(i-1))
!!$     psi_z(i,j)=(psi(i,j+1)-psi(i,j))/(zarray(j+1)-zarray(j))
!!$  enddo
!!$
!!$
!!$  j=nz
!!$  do i=2,nx-1
!!$     psi_x(i,j)=(psi(i+1,j)-psi(i-1,j))/(xarray(i+1)-xarray(i-1))
!!$     psi_z(i,j)=(psi(i,j)-psi(i,j-1))/(zarray(j)-zarray(j-1))
!!$  enddo

  psi_gradient=sqrt(psi_x**2+psi_z**2)


end subroutine calculate_poloidal_flux_partial_derivatives
