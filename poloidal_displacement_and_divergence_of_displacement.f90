subroutine poloidal_displacement_and_divergence_of_displacement(ypath,ER_F_Y)
  use precision,only:p_
  use constants,only: ii,one
  use poloidal_harmonics,only: mhtot
  use toroidal_harmonics,only: nh
  use flux_grids,only: nflux
  use radial_module,only: psival_new,pressure_normalized,&
       & rho_normalized,pprime_normalized,qpsi,dpsidra,fpsi,ffprime, b0_axis,r_axis,starting_surface_number,ending_surface_number
  use kernel_module,only: full_kernel0,full_kernel1,full_kernel2,full_kernel3,full_kernel4,full_kernel5, &
       & full_kernel6,full_kernel7,full_kernel8,full_kernel9,full_kernel10, full_kernel11,full_kernel12, &
       & full_kernel13,full_kernel14,full_kernel15,full_kernel16,full_kernel17,full_kernel18,full_kernel19,&
       & full_kernel20,full_kernel21,full_kernel22,full_kernel23,full_kernel24,full_radial_matrix

  implicit none
  complex(p_),intent(in):: ypath(2*mhtot+1,nflux) !poloidal harmonics of P1 and radial displacement, wsq
  complex(p_),intent(out):: ER_F_Y(2*mhtot,nflux)

  complex(p_):: kernel0(mhtot,mhtot),kernel1(mhtot,mhtot),kernel2(mhtot,mhtot),kernel3(mhtot,mhtot), &
       & kernel4(mhtot,mhtot),kernel5(mhtot,mhtot),kernel6(mhtot,mhtot),kernel7(mhtot,mhtot), &
       & kernel8(mhtot,mhtot),kernel9(mhtot,mhtot),kernel10(mhtot,mhtot),kernel11(mhtot,mhtot), &
       & kernel12(mhtot,mhtot),kernel13(mhtot,mhtot),kernel14(mhtot,mhtot),kernel15(mhtot,mhtot), &
       & kernel16(mhtot,mhtot),kernel17(mhtot,mhtot),kernel18(mhtot,mhtot),kernel19(mhtot,mhtot), &
       & kernel20(mhtot,mhtot),kernel21(mhtot,mhtot),kernel22(mhtot,mhtot),kernel23(mhtot,mhtot), &
       & kernel24(mhtot,mhtot)
  complex(p_):: matrix_c(2*mhtot,2*mhtot),matrix_d(2*mhtot,2*mhtot),matrix_e(2*mhtot,2*mhtot),matrix_f(2*mhtot,2*mhtot)
  complex(p_):: matrix_e_inverse(2*mhtot,2*mhtot)
  complex(p_):: y(2*mhtot) !the array consisting of harmonics of P1 and Xir
  complex(p_):: F_Y(2*mhtot),wsq,sum
  integer:: k,i,j

  do k=starting_surface_number,ending_surface_number

     do i=1,2*mhtot
        y(i)=ypath(i,k) !fetch the values of eigen functions on this flux surface
     enddo
     wsq=ypath(2*mhtot+1,k) !the No. 2*mhtot+1 unknown function is wsq

     !then fetch the kernel matrixes on this magnetic surface
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



     call calculate_matrix_elements(wsq,nh,r_axis,b0_axis, &
          & rho_normalized(k),pressure_normalized(k),pprime_normalized(k),fpsi(k),ffprime(k),qpsi(k),dpsidra(k),&
          & kernel0,kernel1,kernel2,kernel3,kernel4,kernel5,kernel6,kernel7,kernel8,kernel9,kernel10, &
          & kernel11,kernel12,kernel13, kernel14, kernel15, kernel16, kernel17, kernel18, kernel19, &
          & kernel20,kernel21,kernel22,kernel23,kernel24, matrix_c,matrix_d,matrix_e,matrix_f)

     !invert the matrix E, output is matrix_e_inverse
     call invert_matrix(matrix_e,2*mhtot,matrix_e_inverse)

     !after this, matrix F and ER are known. we can calculate ER*F*Y:
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
        ER_F_Y(i,k)=sum !record the value on every flux surface
     enddo

  end do
end subroutine poloidal_displacement_and_divergence_of_displacement
