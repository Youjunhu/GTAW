subroutine calculate_matrix_elements(wsq,nh,r_axis,b0_axis, &
     & rho_normalized,pressure_normalized,pprime_normalized,fpsi_val,ffprime_val,q_val,dpsidra,&
     & kernel0,kernel1,kernel2,kernel3,kernel4,kernel5,kernel6,kernel7,kernel8,kernel9,kernel10, &
     & kernel11,kernel12,kernel13, kernel14, kernel15, kernel16, kernel17, kernel18, kernel19, &
     & kernel20,kernel21,kernel22, kernel23,kernel24,matrix_c,matrix_d,matrix_e,matrix_f)
  !this subroutine calculates the discrete form of the surface operator matrixes (on one magnetic surface)
  use precision,only:p_
  use constants,only: zero,one,two,four,pi,twopi,gamma,ii !,slow_sound_approximation_for_mode
  use poloidal_harmonics,only: mhtot,mh_low !mh_low is the smallest poloidal mode number (can be negative), mhtot is the total poloidal numbers included in the calculation
  implicit none
  complex(p_),intent(in):: wsq ! wsq=(wave frequency)^2/wa0^2, where wa0 is the characteristic Alfven frequency defined by wa=VA/R0, where VA is the Alfven speed at the magnetic axis, 
  integer,intent(in):: nh !toroidal mode number
  real(p_),intent(in):: r_axis,b0_axis
  real(p_),intent(in):: rho_normalized,pressure_normalized,pprime_normalized,fpsi_val,ffprime_val,q_val,dpsidra

  complex(p_),intent(in):: kernel0(mhtot,mhtot),kernel1(mhtot,mhtot),kernel2(mhtot,mhtot),kernel3(mhtot,mhtot), &
       & kernel4(mhtot,mhtot),kernel5(mhtot,mhtot),kernel6(mhtot,mhtot),kernel7(mhtot,mhtot), &
       & kernel8(mhtot,mhtot),kernel9(mhtot,mhtot),kernel10(mhtot,mhtot),kernel11(mhtot,mhtot), &
       & kernel12(mhtot,mhtot),kernel13(mhtot,mhtot),kernel14(mhtot,mhtot),kernel15(mhtot,mhtot), &
       & kernel16(mhtot,mhtot),kernel17(mhtot,mhtot),kernel18(mhtot,mhtot),kernel19(mhtot,mhtot), &
       & kernel20(mhtot,mhtot),kernel21(mhtot,mhtot),kernel22(mhtot,mhtot),kernel23(mhtot,mhtot),&
       & kernel24(mhtot,mhtot)

  complex(p_),intent(out):: matrix_c(2*mhtot,2*mhtot),matrix_d(2*mhtot,2*mhtot),&
       & matrix_e(2*mhtot,2*mhtot),matrix_f(2*mhtot,2*mhtot) !matrix_e_inverse(2*mhtot,2*mhtot)
  complex(p_)::matrix_ea(2*mhtot,2*mhtot),matrix_eb(2*mhtot,2*mhtot)

  complex(p_):: c11(mhtot,mhtot),c12(mhtot,mhtot),c21(mhtot,mhtot),c22(mhtot,mhtot)
  complex(p_):: d11(mhtot,mhtot),d12(mhtot,mhtot),d21(mhtot,mhtot),d22(mhtot,mhtot)
  !complex(p_):: e11(mhtot,mhtot),e12(mhtot,mhtot),e21(mhtot,mhtot),e22(mhtot,mhtot)
  complex(p_):: f11(mhtot,mhtot),f12(mhtot,mhtot),f21(mhtot,mhtot),f22(mhtot,mhtot)
  complex(p_):: h11(mhtot,mhtot),h12(mhtot,mhtot),h21(mhtot,mhtot),h22(mhtot,mhtot)
  integer:: i,j,mhj
  real(p_):: factor

  do i=1,mhtot
     do j=1,mhtot
        mhj=mh_low+(j-1)
        factor=mhj-nh*q_val !resonant term

        c11(i,j)=two*kernel15(i,j)
        c12(i,j)=two/r_axis**2*wsq*rho_normalized*kernel0(i,j) &
             & + two*dpsidra**2/b0_axis**2*(ii*factor*kernel21(i,j)-factor**2*kernel22(i,j)) &
             &  -two/b0_axis**2*kernel17(i,j)+two*pprime_normalized*kernel15(i,j)
        c21(i,j)=zero
        !c22(i,j)=pprime_normalized*b0_axis**2/two*kernel13(i,j)+ffprime_val*kernel14(i,j)-kernel15(i,j)-kernel16(i,j)
        c22(i,j)=-kernel16(i,j)

        d11(i,j)=-ii*factor*two*dpsidra/b0_axis**2*kernel18(i,j)
        d12(i,j)=two*gamma*pressure_normalized*kernel15(i,j)
        d21(i,j)=-ii*factor*dpsidra*fpsi_val*kernel23(i,j)+ii*nh*kernel13(i,j)+two*kernel24(i,j)
        d22(i,j)=kernel13(i,j) +dpsidra**2*gamma*pressure_normalized/(two*wsq*rho_normalized)*r_axis**2 &
             & *(ii*factor*kernel14(i,j)-factor**2*kernel3(i,j))


        h11(i,j)=-ii*mhj*dpsidra*kernel19(i,j)+ii*nh*dpsidra*kernel20(i,j)
        h12(i,j)=zero
        h21(i,j)=zero
        h22(i,j)=h11(i,j)


!!$        e11(i,j)= -two*wsq*rho_normalized/r_axis**2*kernel5(i,j) &
!!$             -two/b0_axis**2*ii*factor*kernel4(i,j)+two/b0_axis**2*factor**2*kernel3(i,j) 
!!$
!!$        e12(i,j)=-two*gamma*pressure_normalized*kernel6(i,j)
!!$
!!$        e21(i,j)=four/b0_axis**2*kernel6(i,j)
!!$
!!$        if(slow_sound_approximation_for_mode.eqv..true.) then
!!$           e22(i,j)=two/b0_axis**2*kernel0(i,j)+ gamma*pressure_normalized*kernel7(i,j)
!!$        else
!!$           e22(i,j)= two/b0_axis**2*kernel0(i,j)+ gamma*pressure_normalized*kernel7(i,j) &
!!$                +r_axis**2/b0_axis**2*gamma*pressure_normalized/(wsq*rho_normalized)* &
!!$                & (ii*factor*kernel2(i,j)-factor**2*kernel1(i,j)) 
!!$        endif
        f11(i,j)=two*kernel6(i,j)-ii*dpsidra*fpsi_val*factor*kernel8(i,j)+ii*nh*kernel0(i,j)
        f12(i,j)=-two*dpsidra/b0_axis**2*ii*factor*kernel9(i,j) &
             & +two*dpsidra/b0_axis**2*(ii*factor*kernel10(i,j)+ kernel11(i,j)) &
             & +two*pprime_normalized*kernel6(i,j)

        f21(i,j)=-kernel7(i,j)
        f22(i,j)=-four/b0_axis**2*kernel12(i,j)

     enddo
  enddo

  !merge matrix h with matrix c and the resulting matrix is still denoted by c:
  c11=c11+h11
  c12=c12+h12
  c21=c21+h21
  c22=c22+h22

  !Next, contructe a (2*mhtot x 2*mhtot) matrix C, which is composed of c11,c12,c21,c22, as follows
  do i=1,mhtot
     do j=1,mhtot
        matrix_c(i,j)=c11(i,j)
        matrix_c(i,j+mhtot)=c12(i,j)
        matrix_c(i+mhtot,j)=c21(i,j)
        matrix_c(i+mhtot,j+mhtot)=c22(i,j)
     enddo
  enddo
  !Next, contructe a (2*mhtot x 2*mhtot) matrix d, which is composed of d11,d12,d21,d22, as follows
  do i=1,mhtot
     do j=1,mhtot
        matrix_d(i,j)=d11(i,j)
        matrix_d(i,j+mhtot)=d12(i,j)
        matrix_d(i+mhtot,j)=d21(i,j)
        matrix_d(i+mhtot,j+mhtot)=d22(i,j)
     enddo
  enddo
  !Next, contructe a (2*mhtot x 2*mhtot) matrix f, which is composed of f11,f12,f21,f22, as follows
  do i=1,mhtot
     do j=1,mhtot
        matrix_f(i,j)=f11(i,j)
        matrix_f(i,j+mhtot)=f12(i,j)
        matrix_f(i+mhtot,j)=wsq*f21(i,j)
        matrix_f(i+mhtot,j+mhtot)=wsq*f22(i,j)
     enddo
  enddo

  !Next, contructe a (2*mhtot x 2*mhtot) matrix E, which is composed of e11,e12,e21,e22, as follows
!!$  do i=1,mhtot
!!$     do j=1,mhtot
!!$        matrix_e(i,j)=e11(i,j)
!!$        matrix_e(i,j+mhtot)=e12(i,j)
!!$        matrix_e(i+mhtot,j)=e21(i,j)
!!$        matrix_e(i+mhtot,j+mhtot)=e22(i,j)
!!$     enddo
!!$  enddo

 !a subroutine for calculating the matrix E is available, thus call it, instead of using the above code:
  call calculate_matrix_e(nh,r_axis,b0_axis,rho_normalized,pressure_normalized,q_val, dpsidra,&
       & kernel0,kernel1,kernel2,kernel3,kernel4,kernel5,kernel6,kernel7,matrix_ea,matrix_eb) !calculate matrix E 
  matrix_e=matrix_ea+wsq*matrix_eb
  !invert the matrix E, output is matrix_e_inverse
  !call invert_matrix(matrix_e,2*mhtot,matrix_e_inverse)

end subroutine calculate_matrix_elements



subroutine calculate_matrix_e(nh,r_axis,b0_axis,rho_normalized,pressure_normalized,q_val, dpsidra,&
     & kernel0,kernel1,kernel2,kernel3,kernel4,kernel5,kernel6,kernel7,matrix_a,matrix_b)
  !finding the continuous spectrum reduces to solving a generalized eigenvale problem, which takes the form matrix_a*x=-wsq**2*matrix_b*x
  !this routine is to construct the matrix_a and matrix_b
  use precision,only:p_
  use constants,only: zero,one,two,four,pi,twopi,gamma,slow_sound_approximation
  use poloidal_harmonics,only: mhtot,mh_low

  implicit none
  integer,intent(in):: nh !toroidal mode number
  real(p_),intent(in):: r_axis,b0_axis !normalizing constants
  real(p_),intent(in):: rho_normalized,pressure_normalized,q_val,dpsidra
  complex(p_),intent(in):: kernel0(mhtot,mhtot),kernel1(mhtot,mhtot),kernel2(mhtot,mhtot),kernel3(mhtot,mhtot)
  complex(p_),intent(in):: kernel4(mhtot,mhtot),kernel5(mhtot,mhtot),kernel6(mhtot,mhtot),kernel7(mhtot,mhtot)
  complex(p_),intent(out):: matrix_a(2*mhtot,2*mhtot),matrix_b(2*mhtot,2*mhtot) !the two matrixes that form the generalized eigenvale equation

  !here two matrixs are returned, the final matrix we need is that matrix=matrix_a+wsq*matrix_b
  !complex(p_),intent(out):: matrix_a(2*mhtot,2*mhtot),matrix_b(2*mhtot,2*mhtot)
  complex(p_),parameter:: ii=(0.0_p_,1.0_p_)
  integer:: i,j,mhj
  real(p_):: factor
  !  character*80:: format_string,format_string2
  complex(p_):: e11_a(mhtot,mhtot),e12_a(mhtot,mhtot),e21_a(mhtot,mhtot),e22_a(mhtot,mhtot)
  complex(p_):: e11_b(mhtot,mhtot),e12_b(mhtot,mhtot),e21_b(mhtot,mhtot),e22_b(mhtot,mhtot)


  do i=1,mhtot
     do j=1,mhtot
        mhj=mh_low+(j-1)
        factor=mhj-nh*q_val
        e11_a(i,j)= -two*dpsidra**2/b0_axis**2*(ii*factor*kernel4(i,j)-factor**2*kernel3(i,j))
        e11_b(i,j)= -two*rho_normalized/r_axis**2*kernel5(i,j)
        e12_a(i,j)=-two*gamma*pressure_normalized*kernel6(i,j)
        e12_b(i,j)=0.
        e21_a(i,j)=0.
        e21_b(i,j)=four/b0_axis**2*kernel6(i,j)
        if(slow_sound_approximation.eqv..true.) then
           e22_a(i,j)=0._p_
        else
           e22_a(i,j)= r_axis**2/b0_axis**2*gamma*pressure_normalized/rho_normalized*dpsidra**2* &
                & (ii*factor*kernel2(i,j)-factor**2*kernel1(i,j))
        endif
        e22_b(i,j)=two/b0_axis**2*kernel0(i,j)+ gamma*pressure_normalized*kernel7(i,j)
     enddo
  enddo

  !store them in a bigger matrix
  do i=1,mhtot
     do j=1,mhtot
        matrix_a(i,j)=e11_a(i,j)
        matrix_a(i,j+mhtot)=e12_a(i,j)
        matrix_a(i+mhtot,j)=e21_a(i,j)
        matrix_a(i+mhtot,j+mhtot)=e22_a(i,j)
     enddo
  enddo

  do i=1,mhtot
     do j=1,mhtot
        matrix_b(i,j)=e11_b(i,j)
        matrix_b(i,j+mhtot)=e12_b(i,j)
        matrix_b(i+mhtot,j)=e21_b(i,j)
        matrix_b(i+mhtot,j+mhtot)=e22_b(i,j)
     enddo
  enddo

end subroutine calculate_matrix_e
