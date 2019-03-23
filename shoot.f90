!U    SUBROUTINE shoot(n2,v,f) is named "funcv" for use with "newt"
SUBROUTINE funcv(n2,v,f)
!SUBROUTINE funcv_shoot(n2,v,f) !rename it to funcv if you want to use it in shooting process
  use precision,only: p_
  use radial_module, only: ra,starting_surface_number,ending_surface_number
  use poloidal_harmonics,only: mhtot
  implicit none

  INTEGER,intent(in):: n2
  complex(p_),intent(in):: v(n2)
  complex(p_),intent(out):: f(n2)

  real(p_):: x1,x2
  complex(p_):: y(2*mhtot+1) !2*mhtot+1 is the number of functions, which is also the number of ordinary difference equations

  !x1=psival_new(starting_surface_number)
  !x2=psival_new(ending_surface_number)
  x1=ra(starting_surface_number)
  x2=ra(ending_surface_number)

  call load(x1,v,y) !set value of y at the left boundary point by using v(n2) array
  call odeint_yj(y,2*mhtot+1,x1,x2) !advancing
  call score(x2,y,f)  !calculate the derivation of y from the given boundary condition at the right boundary

END SUBROUTINE


SUBROUTINE load(x1,v,y) 
!This routine sets the values of 2*mhtot+1 (complex) unknown fucntions y at x=x1 through the inital guess (stored in v array of length mhtot) and the known boundary conditions for the mhtot poloidal harmonics of plasma displacement. Therefore, there remains (2*mhtot+1-mhtot-mhtot)=1 function whose values at x1 need to be set. In this routine, the values of the last poloidal harmonic of the perturbed pressure at x1 are set to be a very small number.
  use precision,only:p_
  use poloidal_harmonics,only: mhtot
  implicit none
  !  real(p_),intent(in):: x1, v(2*mhtot) !real number array of 2*mhtot length, which will represent 1*mhtot complex variables
  real(p_),intent(in):: x1
  complex(p_),intent(in)::  v(mhtot) !real number array of 2*mhtot length, which will represent 1*mhtot complex variables
  complex(p_),intent(out):: y(2*mhtot+1) !values of 2*mhtot+1 unknown functions (the last unknown function is wsq) at the left boundary point
  !complex(p_),parameter:: ii=(0.0,1.0_p_)
  integer:: i

  !the first mhtot functions are poloidal harmonics of P1 (P1 is the sum of perturbed thermal pressure and  magnetic pressure), the second mhtot functions are  poloidal harmonics of radial displacement, the last function is wsq.
  !Note that both the P1 and the radial displacement are complex functions, wsq is also asumed to be a complex number although it is actually a real number in ideal MHD case.

!!$  k=1
!!$  do i=1,mhtot-1 !note the store arrangement of v arrray: v(k) and v(k+1) store respectively the real and imaginary part of a complex number
!!$     y(i)=v(k)+ii*v(k+1)
!!$     k=k+2
!!$  enddo
  do i=1,mhtot-1
     y(i)=v(i) !initial guess
  enddo

  !The following is the known boundary conditions
  y(mhtot)=(1.0d-6,1.0d-6) ! the value of the last poloidal harmonics of P1 at innermost magnetic surface is set to be a small number
  !y(mhtot)=v(mhtot-1)*0.8_p_ ! the value of the last poloidal harmonics of P1 at innermost magnetic surface is set to be equal to  the  value of the second last harmonic
  !y(mhtot)=(1.0d-16,1.0d-16) ! the value of the last poloidal harmonics of P1 at innermost magnetic surface is set to be a small number
  !y(mhtot)=(1.0d-8,1.0d-8) ! the value of the last poloidal harmonics of P1 at innermost magnetic surface is set to be a small number
  do i=mhtot+1,2*mhtot !the values of all the poloidal harmonics of radial displacement at innermost magnetic surface are set to be zero
     y(i)=(1.0d-12,1.0d-12)
     ! y(i)=(1.0d-6,1.0d-6)
     !y(i)=(0.0_p_,0.0_p_)
     ! y(i)=(1.0d-7,1.0d-7)
  enddo

  !y(2*mhtot+1)=v(2*mhtot-1)+ii*v(2*mhtot) !the last unknown function is wsq
  y(2*mhtot+1)=v(mhtot) !the last unknown function is wsq
  !write(*,*) (real(y(i)),i=1,2*mhtot,3)
  !write(*,*) 'Y=',Y
end SUBROUTINE load


SUBROUTINE score(x2,y,f)
  use precision,only:p_
  use poloidal_harmonics,only: mhtot
  implicit none

  INTEGER i
  real(p_),intent(in):: x2
  complex(p_),intent(in):: y(2*mhtot+1)
!  real(p_),intent(out):: f(2*mhtot)
  complex(p_),intent(out):: f(mhtot)

!!$  do i=1,mhtot !the goal is to make all the mhtot poloidal harmonics of the plasma displacement be zero.
!!$     f(i)=real(y(i+mhtot))-0.0
!!$     f(i+mhtot)=imag(y(i+mhtot))-0.0
!!$  enddo

  do i=1,mhtot !the goal is to make all the mhtot poloidal harmonics of the plasma displacement be zero.
     f(i)=y(i+mhtot)-0.0
  enddo

!the values of fi are the derivation of values of y at the right boundary point (obtained by Runk-Kutta solver) from the values specified by the known boundary conditions at the right boundary point (The right boundary condition is that some of the values of yi are zero at the point).
end SUBROUTINE score


SUBROUTINE derivs(x,y,dydx)
  use precision,only:p_
  use constants,only: two
  use poloidal_harmonics,only: mhtot
  use toroidal_harmonics,only: nh
  use flux_grids,only: nflux
  use radial_module,only: dra,ra,pressure_normalized,&
       & rho_normalized,pprime_normalized,qpsi,fpsi,ffprime, b0_axis,r_axis,dpsidra
  use kernel_module,only: full_kernel0,full_kernel1,full_kernel2,full_kernel3,full_kernel4,full_kernel5, &
       & full_kernel6,full_kernel7,full_kernel8,full_kernel9,full_kernel10, full_kernel11,full_kernel12, &
       & full_kernel13,full_kernel14,full_kernel15,full_kernel16,full_kernel17,full_kernel18,full_kernel19,&
       & full_kernel20,full_kernel21,full_kernel22,full_kernel23,full_kernel24,full_radial_matrix

  implicit none
  real(p_),intent(in):: x !x is the value of PSi, which determines the radial location
  complex(p_),intent(in):: y(2*mhtot+1) !the unknown functions: poloidal harmonics of P1 and radial displacement, wsq
  complex(p_),intent(out):: dydx(2*mhtot+1) !derivative of y with respect to x
  complex(p_):: wsq
  complex(p_):: kernel0(mhtot,mhtot),kernel1(mhtot,mhtot),kernel2(mhtot,mhtot),kernel3(mhtot,mhtot), &
       & kernel4(mhtot,mhtot),kernel5(mhtot,mhtot),kernel6(mhtot,mhtot),kernel7(mhtot,mhtot), &
       & kernel8(mhtot,mhtot),kernel9(mhtot,mhtot),kernel10(mhtot,mhtot),kernel11(mhtot,mhtot), &
       & kernel12(mhtot,mhtot),kernel13(mhtot,mhtot),kernel14(mhtot,mhtot),kernel15(mhtot,mhtot), &
       & kernel16(mhtot,mhtot),kernel17(mhtot,mhtot),kernel18(mhtot,mhtot),kernel19(mhtot,mhtot), &
       & kernel20(mhtot,mhtot),kernel21(mhtot,mhtot),kernel22(mhtot,mhtot),kernel23(mhtot,mhtot), &
       & kernel24(mhtot,mhtot)
  complex(p_):: radial_matrix(2*mhtot,2*mhtot)

  complex(p_):: matrix_c(2*mhtot,2*mhtot),matrix_d(2*mhtot,2*mhtot),matrix_e(2*mhtot,2*mhtot),matrix_f(2*mhtot,2*mhtot)
  complex(p_):: matrix_e_inverse(2*mhtot,2*mhtot)

  complex(p_):: C_Y(2*mhtot),F_Y(2*mhtot),ER_F_Y(2*mhtot),D_ER_F_Y(2*mhtot)
  complex(p_):: dydx0(2*mhtot),sum
!  real(p_):: interval
  integer:: k,i,j
  !complex(p_):: test_matrix(2*mhtot,2*mhtot)

  !infering from the value of x which magnetic surface we are on.
  !interval=(psival_new(nflux)-psival_new(1))/(nflux-1)
  !k=nint(1+(x-psival_new(1))/interval) !nint() is a fortran intrinsic function, which returns the nearest integer to its argument
   k=nint(1+(x-ra(1))/dra) !nint() is a fortran intrinsic function, which returns the nearest integer to its argument
  !write(*,*) 'magnetic surface No.', k
  !then fetch the kernel matrixes of this magnetic surface
  do i=1,mhtot
     do j=1,mhtot
        kernel0(i,j)=full_kernel0(k,i,j)
        kernel1(i,j)=full_kernel1(k,i,j)
        kernel2(i,j)=full_kernel2(k,i,j)
        kernel3(i,j)=full_kernel3(k,i,j)
        kernel4(i,j)=full_kernel4(k,i,j)
        kernel5(i,j)=full_kernel5(k,i,j)
        kernel6(i,j)=full_kernel6(k,i,j)
        kernel7(i,j)=full_kernel7(k,i,j)
        kernel8(i,j)=full_kernel8(k,i,j)
        kernel9(i,j)=full_kernel9(k,i,j)
        kernel10(i,j)=full_kernel10(k,i,j)
        kernel11(i,j)=full_kernel11(k,i,j)
        kernel12(i,j)=full_kernel12(k,i,j)
        kernel13(i,j)=full_kernel13(k,i,j)
        kernel14(i,j)=full_kernel14(k,i,j)
        kernel15(i,j)=full_kernel15(k,i,j)
        kernel16(i,j)=full_kernel16(k,i,j)
        kernel17(i,j)=full_kernel17(k,i,j)
        kernel18(i,j)=full_kernel18(k,i,j)
        kernel19(i,j)=full_kernel19(k,i,j)
        kernel20(i,j)=full_kernel20(k,i,j)
        kernel21(i,j)=full_kernel21(k,i,j)
        kernel22(i,j)=full_kernel22(k,i,j)
        kernel23(i,j)=full_kernel23(k,i,j)
        kernel24(i,j)=full_kernel24(k,i,j)
     enddo
  enddo

  do i=1,2*mhtot
     do j=1,2*mhtot
        radial_matrix(i,j)=full_radial_matrix(k,i,j)
     enddo
  enddo

  wsq=y(2*mhtot+1) !the No. 2*mhtot+1 unknown function is wsq

  call calculate_matrix_elements(wsq,nh,r_axis,b0_axis, &
       & rho_normalized(k),pressure_normalized(k),pprime_normalized(k),fpsi(k),ffprime(k),qpsi(k),dpsidra(k),&
       & kernel0,kernel1,kernel2,kernel3,kernel4,kernel5,kernel6,kernel7,kernel8,kernel9,kernel10, &
       & kernel11,kernel12,kernel13, kernel14, kernel15, kernel16, kernel17, kernel18, kernel19, &
       & kernel20,kernel21,kernel22,kernel23,kernel24, matrix_c,matrix_d,matrix_e,matrix_f)

  !invert the matrix E, output is matrix_e_inverse
  call invert_matrix(matrix_e,2*mhtot,matrix_e_inverse)

  !check whether matrix_e*matrix_e_inverse equals unit matrix
!!$  do i=1,2*mhtot
!!$     do j=1,2*mhtot
!!$        sum=0.
!!$        do k=1,2*mhtot
!!$           sum=sum+matrix_e(i,k)*matrix_e_inverse(k,j)
!!$        enddo
!!$        test_matrix(i,j)=sum
!!$     enddo
!!$  enddo
!!$  do i=2,2*mhtot
!!$     write(*,*) test_matrix(i,i-1)
!!$  enddo
!!$  sum=0.
!!$  do i=1,2*mhtot
!!$     do j=1,2*mhtot
!!$        sum=sum+abs(test_matrix(i,j))
!!$     enddo
!!$     enddo
!!$  
!!$  write(*,*) "sum of absolute value of matrix_e*matrix_e_inverse", sum/(2*mhtot)
!!$  stop 'test finished'


  !after this, matrix f,er,d, and c are ready to be used.
  !next, calculate the right-hand side of the eigenmode equations
  do i=1,2*mhtot
     sum=0.
     do j=1,2*mhtot
        sum=sum+matrix_f(i,j)*y(j)
     enddo
     F_Y(i)=sum !F_Y is the product of matrix F with collomn vector Y
  enddo

  do i=1,2*mhtot
     sum=0.
     do j=1,2*mhtot
        sum=sum+matrix_e_inverse(i,j)*F_Y(j)
     enddo
     ER_F_Y(i)=sum
  enddo


  do i=1,2*mhtot
     sum=0.
     do j=1,2*mhtot
        sum=sum+matrix_d(i,j)*ER_F_Y(j)
     enddo
     D_ER_F_Y(i)=sum
  enddo

!!$  write(*,*) 'D_ER_F_Y='
!!$  do i=1,2*mhtot
!!$     write(*,*) D_ER_F_Y(i)
!!$  enddo

  ! write(*,*) 'd22=',d22
  do i=1,2*mhtot
     sum=0.
     do j=1,2*mhtot
        sum=sum+matrix_c(i,j)*y(j)
     enddo
     C_Y(i)=sum
  enddo

  do i=1,2*mhtot
     dydx0(i)=C_Y(i)+D_ER_F_Y(i)
  enddo

  !radial matrix
  do i=1,2*mhtot
     sum=0.
     do j=1,2*mhtot
        sum=sum+radial_matrix(i,j)*dydx0(j)
     enddo
     dydx(i)=sum
  enddo

  dydx(2*mhtot+1)=0.0_p_ !the equation for wsq

!!$  do i=mhtot+1,2*mhtot
!!$     write(*,*) C_Y(i),D_ER_F_Y(i)
!!$enddo
!!$write(*,*) 'dydx='
!!$if (k.eq.5) then
!!$   write(*,*) (dydx(i),i=1,2*mhtot)
!!$ endif
end SUBROUTINE derivs

subroutine odeint_yj(y,nvar,x1,x2)
  use precision,only:p_
  use radial_module, only: starting_surface_number,ending_surface_number,dra
  implicit none
  integer,intent(in):: nvar
  complex(p_),intent(inout)::y(nvar) !as input, y provides the values of fuctions at x1, as output, y provides the values at x2
  real(p_),intent(in):: x1,x2
  real(p_):: x,x_old
  complex(p_)::dydx(nvar)
  integer:: j
  complex(p_):: y_old(nvar)

  x=x1
  do j=starting_surface_number,ending_surface_number-1 !advance 
!1st order Euler
     call derivs(x,y,dydx) !prepare derivatives for advance
     y=y+dydx*dra
     x=x+dra

!2nd Runge-Kutta scheme, wrong because the quantities at the middle-point is not calculated
!!$     y_old=y
!!$     x_old=x
!!$     call derivs(x,y,dydx) !prepare derivatives for advance
!!$     y=y+dydx*dra*0.5_p_
!!$     x=x+dra*0.5_p_
!!$
!!$     call derivs(x,y,dydx) !prepare derivatives for advance
!!$     y=y_old+dydx*dra
!!$     x=x_old+dra

  end do
end subroutine odeint_yj


subroutine odeint_yj2(y,nvar,x1,x2,xpath,ypath)
  use precision,only:p_
  use constants,only:one,one_half
  use flux_grids,only: nflux
  use radial_module, only: starting_surface_number,ending_surface_number,dra
  !use path,only: xpath,ypath_real,ypath_imag !record the path (eigenfucntion)
  !use path,only: xpath,ypath !record the path (eigenfucntion)
  implicit none
  integer,intent(in):: nvar
  complex(p_),intent(inout)::y(nvar) !as input, y provides the values of fuctions at x1, as output, y provides the values at x2
  real(p_),intent(in):: x1,x2
  real(p_),intent(out):: xpath(nflux)
  complex(p_),intent(out):: ypath(nvar,nflux)
  real(p_):: x,x_old
  complex(p_)::dydx(nvar)
  complex(p_):: k1(nvar),k2(nvar)
  complex(p_):: y_old(nvar)

  integer:: i,j
  !external:: derivs
  x=x1
  !record the initial location and values of functions at this location
  xpath(starting_surface_number)=x1
  do i=1,nvar
     ypath(i,starting_surface_number)=y(i)
  enddo

  do j=starting_surface_number,ending_surface_number-1 !advance 
!1st order Euler scheme
     call derivs(x,y,dydx) !prepare derivatives for advance
     y=y+dydx*dra
     x=x+dra

!2nd Runge-Kutta scheme,wrong!
!!$     y_old=y
!!$     x_old=x
!!$     call derivs(x,y,dydx) !prepare derivatives for advance
!!$     y=y+dydx*dra*0.5_p_
!!$     x=x+dra*0.5_p_
!!$
!!$     call derivs(x,y,dydx) !prepare derivatives for advance
!!$     y=y_old+dydx*dra
!!$     x=x_old+dra

     !record the radial location and values of functions at this location
     xpath(j+1)=x 
     do i=1,nvar 
        ypath(i,j+1)=y(i)
     enddo

  end do

end subroutine odeint_yj2

