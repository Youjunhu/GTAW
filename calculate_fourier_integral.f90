subroutine calculate_fourier_integral(mpoloidal,ws1,ws2,ws3,ws4,ws5,ws6,ws7, ws8,ws9,&
     & ws10,ws11,ws12,ws13,ws14,ws15,ws16,ws17,ws18,ws19,ws20,ws21,ws22,ws23,ws24,&
     & kernel0,kernel1,kernel2,kernel3,kernel4,kernel5,kernel6,kernel7,&
     & kernel8,kernel9,kernel10,kernel11,kernel12,kernel13,kernel14,&
     & kernel15,kernel16,kernel17,kernel18,kernel19,kernel20,kernel21,kernel22,kernel23,kernel24)
  !this subroutine calculate the fourier integration (in the direction of theta) on a single magnetic surface
  use precision,only:p_
  use constants,only: zero,one,two,twopi
  use poloidal_harmonics,only: mhtot,mh_low,mh_upp
  implicit none
  integer,intent(in):: mpoloidal
  real(p_),intent(in):: ws1(mpoloidal),ws2(mpoloidal),ws3(mpoloidal),ws4(mpoloidal), ws5(mpoloidal)
  real(p_),intent(in):: ws6(mpoloidal),ws7(mpoloidal),ws8(mpoloidal),ws9(mpoloidal), ws10(mpoloidal)
  real(p_),intent(in):: ws11(mpoloidal),ws12(mpoloidal),ws13(mpoloidal),ws14(mpoloidal), ws15(mpoloidal)
  real(p_),intent(in):: ws16(mpoloidal),ws17(mpoloidal),ws18(mpoloidal),ws19(mpoloidal), ws20(mpoloidal)
  real(p_),intent(in):: ws21(mpoloidal),ws22(mpoloidal),ws23(mpoloidal),ws24(mpoloidal)
  complex(p_),intent(out):: kernel0(mhtot,mhtot),kernel1(mhtot,mhtot),kernel2(mhtot,mhtot),kernel3(mhtot,mhtot)
  complex(p_),intent(out):: kernel4(mhtot,mhtot),kernel5(mhtot,mhtot),kernel6(mhtot,mhtot),kernel7(mhtot,mhtot)
  complex(p_),intent(out):: kernel8(mhtot,mhtot),kernel9(mhtot,mhtot),kernel10(mhtot,mhtot),kernel11(mhtot,mhtot)
  complex(p_),intent(out):: kernel12(mhtot,mhtot),kernel13(mhtot,mhtot),kernel14(mhtot,mhtot),kernel15(mhtot,mhtot)
  complex(p_),intent(out):: kernel16(mhtot,mhtot),kernel17(mhtot,mhtot),kernel18(mhtot,mhtot),kernel19(mhtot,mhtot)
  complex(p_),intent(out):: kernel20(mhtot,mhtot),kernel21(mhtot,mhtot),kernel22(mhtot,mhtot)
  complex(p_),intent(out):: kernel23(mhtot,mhtot),kernel24(mhtot,mhtot)
  integer:: mhi,mhj,mh
  integer:: i,j

  !real(p_):: ws0(mpoloidal),kernel00(mhtot,mhtot)

  call fouier_integral(ws1,kernel1) !output kernel1(mhtot,mhtot)
  call fouier_integral(ws2,kernel2) !output kernel2(mhtot,mhtot)
  call fouier_integral(ws3,kernel3) !output kernel3(mhtot,mhtot)
  call fouier_integral(ws4,kernel4) !output kernel4(mhtot,mhtot)
  call fouier_integral(ws5,kernel5) !output kernel5(mhtot,mhtot)
  call fouier_integral(ws6,kernel6) !output kernel6(mhtot,mhtot)
  call fouier_integral(ws7,kernel7) !output kernel7(mhtot,mhtot)
  call fouier_integral(ws8,kernel8) 
  call fouier_integral(ws9,kernel9) 
  call fouier_integral(ws10,kernel10) 
  call fouier_integral(ws11,kernel11) 
  call fouier_integral(ws12,kernel12) 
  call fouier_integral(ws13,kernel13) 
  call fouier_integral(ws14,kernel14) 
  call fouier_integral(ws15,kernel15) 
  call fouier_integral(ws16,kernel16) 
  call fouier_integral(ws17,kernel17) 
  call fouier_integral(ws18,kernel18) 
  call fouier_integral(ws19,kernel19) 
  call fouier_integral(ws20,kernel20) 
  call fouier_integral(ws21,kernel21) 
  call fouier_integral(ws22,kernel22) 
  call fouier_integral(ws23,kernel23) 
  call fouier_integral(ws24,kernel24) 


!!$  do i=1,mpoloidal
!!$     ws0(i)=1.0_p_
!!$  enddo
!!$  call fouier_integral(ws0,kernel0)

  kernel1=kernel1/twopi
  kernel2=kernel2/twopi
  kernel3=kernel3/twopi
  kernel4=kernel4/twopi
  kernel5=kernel5/twopi
  kernel6=kernel6/twopi
  kernel7=kernel7/twopi
  kernel8=kernel8/twopi
  kernel9=kernel9/twopi
  kernel10=kernel10/twopi
  kernel11=kernel11/twopi
  kernel12=kernel12/twopi
  kernel13=kernel13/twopi
  kernel14=kernel14/twopi
  kernel15=kernel15/twopi
  kernel16=kernel16/twopi
  kernel17=kernel17/twopi
  kernel18=kernel18/twopi
  kernel19=kernel19/twopi
  kernel20=kernel20/twopi
  kernel21=kernel21/twopi
  kernel22=kernel22/twopi
  kernel23=kernel23/twopi
  kernel24=kernel24/twopi

!!$  kernel0=kernel0/twopi

  !kernel0 is the inner product with the wight function w=1. This can be done analytically
  do i=1,mhtot
     mhi=mh_low+(i-1)
     do j=1,mhtot
        mhj=mh_low+(j-1)
        !mh=mhi+mhj
        mh=-mhi+mhj
        if (mh .eq.0) then
           kernel0(i,j)=1.0_p_
        else
           kernel0(i,j)=0._p_
        endif
     enddo
  enddo

end subroutine calculate_fourier_integral


subroutine fouier_integral_origin(weight_array,kernel)
!This subroutine makes use of "dfint" subroutine given in Numerical recipe book to calculate the fourier integration.
!when I calculateing the continua for DIIID high q equilibrium, I need to include more poloidal harmonics and thus need to increase the value of mh_upp.
!However, when mh_upp=50, the dftcor subroutine generates the error message: "STOP bad arguments to dftcor".
!To make GTAW able to deal with larger mh_upp, I decide to writing a new subroutine to replace this one.
!The new subroutine (at the end of this file) uses simple numerical integration formula to calculate the Fourier integratiion.
!The results indicate that there is no obvious difference between the results of continua given by "dftint" and the new subroutine.
!This indicates that the very tricky methods used in "dftint" may be not necessary for the present case.
!the cpu time used by the two subroutines are similar
  use precision,only: p_
  use constants,only: zero,one,twopi,ii
  use flux_grids,only: mpoloidal
  use poloidal_harmonics,only:mh_low,mh_upp,mhtot
  implicit none
  real(p_),intent(in):: weight_array(mpoloidal)
  complex(p_),intent(out):: kernel(mhtot,mhtot)
  integer:: iii,m,i,j,mhi,mhj,mh
  real(p_):: mh_real
  real(p_):: sinint,cosint,sign
  real(p_):: sinint_array(0:mh_upp-mh_low),cosint_array(0:mh_upp-mh_low)
  logical:: recalculate

  m=mpoloidal-1
  recalculate=.true.
  do iii=0,mh_upp-mh_low !mh_upp-mh_low is the largest possible difference between two poloidal mode numbers
     mh_real=real(iii)
     call dftint(recalculate,m,weight_array,zero,twopi,mh_real,cosint,sinint)
     !write(*,*) 'cosint= ',cosint,'sinint= ',sinint
     cosint_array(iii)=cosint
     sinint_array(iii)=sinint
  enddo

  do i=1,mhtot
     mhi=mh_low+(i-1)
     do j=1,mhtot
        mhj=mh_low+(j-1)
        !mh=mhi+mhj
        mh=-mhi+mhj
        sign=+one
        if (mh .lt. 0) then
           mh=-mh
           sign=-one
        endif
        kernel(i,j)=cosint_array(mh)+ii*sinint_array(mh)*sign
     enddo
  enddo
end subroutine



subroutine fouier_integral(weight_array,kernel)
!see the comments in subroutine "fouier_integral_origin"
  use precision,only: p_
  use constants,only: zero,one,two,twopi,ii
  use flux_grids,only: mpoloidal
  use poloidal_harmonics,only:mh_low,mh_upp,mhtot
  implicit none
  real(p_),intent(in):: weight_array(mpoloidal)
  complex(p_),intent(out):: kernel(mhtot,mhtot)
  integer:: iii,m,i,j,mhi,mhj,mh
  real(p_):: mh_real
  real(p_):: sinint,cosint,sign,dtheta,theta(mpoloidal)
  real(p_):: sinint_array(0:mh_upp-mh_low),cosint_array(0:mh_upp-mh_low)
!  logical:: recalculate

!  m=mpoloidal-1
  dtheta=twopi/(mpoloidal-1)
  do i=1,mpoloidal
     theta(i)=0._p_+dtheta*(i-1)
  enddo

!  recalculate=.true.
  do iii=0,mh_upp-mh_low !mh_upp-mh_low is the largest possible difference between two poloidal mode numbers
     mh_real=real(iii)

     cosint=0.
     sinint=0.
     do i=1,mpoloidal-1
!        cosint=cosint+weight_array(i)*cos(mh_real*theta(i))*dtheta !rectangular formula
!        sinint=sinint+weight_array(i)*sin(mh_real*theta(i))*dtheta !rectangular formula
        cosint=cosint+(weight_array(i)+weight_array(i+1))/two*cos(mh_real*(theta(i)+theta(i+1))/two)*dtheta !trapezoid formula
        sinint=sinint+(weight_array(i)+weight_array(i+1))/two*sin(mh_real*(theta(i)+theta(i+1))/two)*dtheta !trapezoid formula
     enddo

     !call dftint(recalculate,m,weight_array,zero,twopi,mh_real,cosint,sinint)
     !write(*,*) 'cosint= ',cosint,'sinint= ',sinint
     cosint_array(iii)=cosint
     sinint_array(iii)=sinint
  enddo

  do i=1,mhtot
     mhi=mh_low+(i-1)
     do j=1,mhtot
        mhj=mh_low+(j-1)
        !mh=mhi+mhj
        mh=-mhi+mhj
        sign=+one
        if (mh .lt. 0) then
           mh=-mh
           sign=-one
        endif
        kernel(i,j)=cosint_array(mh)+ii*sinint_array(mh)*sign
     enddo
  enddo
end subroutine
