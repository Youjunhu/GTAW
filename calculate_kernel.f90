subroutine calculate_kernel(mpoloidal,nflux,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,w23,w24)
  use kernel_module,only: full_kernel0,full_kernel1,full_kernel2,full_kernel3,full_kernel4,full_kernel5 !Fourier integrations on every flux surface
  use kernel_module,only: full_kernel6,full_kernel7,full_kernel8,full_kernel9,full_kernel10,full_kernel11
  use kernel_module,only: full_kernel12,full_kernel13,full_kernel14,full_kernel15, full_kernel16,full_kernel17
  use kernel_module,only: full_kernel18,full_kernel19,full_kernel20,full_kernel21,full_kernel22
  use kernel_module,only: full_kernel23,full_kernel24 !as output
  use poloidal_harmonics,only: mhtot

  use precision,only:p_
  use constants,only:zero,one,two
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: w1(mpoloidal,nflux),w2(mpoloidal,nflux),w3(mpoloidal,nflux),w4(mpoloidal,nflux) 
  real(p_),intent(in):: w5(mpoloidal,nflux),w6(mpoloidal,nflux),w7(mpoloidal,nflux),w8(mpoloidal,nflux)
  real(p_),intent(in):: w9(mpoloidal,nflux),w10(mpoloidal,nflux),w11(mpoloidal,nflux),w12(mpoloidal,nflux)
  real(p_),intent(in):: w13(mpoloidal,nflux),w14(mpoloidal,nflux),w15(mpoloidal,nflux),w16(mpoloidal,nflux)
  real(p_),intent(in):: w17(mpoloidal,nflux),w18(mpoloidal,nflux),w19(mpoloidal,nflux),w20(mpoloidal,nflux)
  real(p_),intent(in):: w21(mpoloidal,nflux),w22(mpoloidal,nflux),w23(mpoloidal,nflux),w24(mpoloidal,nflux)

  !weight functions on a single magnetic surface
  real(p_):: ws1(mpoloidal),ws2(mpoloidal),ws3(mpoloidal),ws4(mpoloidal),ws5(mpoloidal),ws6(mpoloidal),ws7(mpoloidal) 
  real(p_):: ws8(mpoloidal),ws9(mpoloidal),ws10(mpoloidal),ws11(mpoloidal),ws12(mpoloidal),ws13(mpoloidal),ws14(mpoloidal) 
  real(p_):: ws15(mpoloidal),ws16(mpoloidal),ws17(mpoloidal),ws18(mpoloidal),ws19(mpoloidal),ws20(mpoloidal),ws21(mpoloidal) 
  real(p_):: ws22(mpoloidal),ws23(mpoloidal),ws24(mpoloidal)
  !inner products between two fourier harmonics, The different numbers correspond to different weight functions
  complex(p_):: kernel0(mhtot,mhtot),kernel1(mhtot,mhtot),kernel2(mhtot,mhtot),kernel3(mhtot,mhtot)
  complex(p_):: kernel4(mhtot,mhtot),kernel5(mhtot,mhtot),kernel6(mhtot,mhtot),kernel7(mhtot,mhtot)
  complex(p_):: kernel8(mhtot,mhtot),kernel9(mhtot,mhtot),kernel10(mhtot,mhtot),kernel11(mhtot,mhtot)
  complex(p_):: kernel12(mhtot,mhtot),kernel13(mhtot,mhtot),kernel14(mhtot,mhtot),kernel15(mhtot,mhtot)
  complex(p_):: kernel16(mhtot,mhtot),kernel17(mhtot,mhtot),kernel18(mhtot,mhtot),kernel19(mhtot,mhtot)
  complex(p_):: kernel20(mhtot,mhtot),kernel21(mhtot,mhtot),kernel22(mhtot,mhtot),kernel23(mhtot,mhtot)
  complex(p_):: kernel24(mhtot,mhtot)
  integer:: k,i,j

  do k=1,nflux
     do i=1,mpoloidal
        ws1(i)=w1(i,k) !here "s" stands for single surface
        ws2(i)=w2(i,k)
        ws3(i)=w3(i,k)
        ws4(i)=w4(i,k)
        ws5(i)=w5(i,k)
        ws6(i)=w6(i,k)
        ws7(i)=w7(i,k)
        ws8(i)=w8(i,k)
        ws9(i)=w9(i,k)
        ws10(i)=w10(i,k)
        ws11(i)=w11(i,k)
        ws12(i)=w12(i,k)
        ws13(i)=w13(i,k)
        ws14(i)=w14(i,k)
        ws15(i)=w15(i,k)
        ws16(i)=w16(i,k)
        ws17(i)=w17(i,k)
        ws18(i)=w18(i,k)
        ws19(i)=w19(i,k)
        ws20(i)=w20(i,k)
        ws21(i)=w21(i,k)
        ws22(i)=w22(i,k)
        ws23(i)=w23(i,k)
        ws24(i)=w24(i,k)
     enddo

     !on a single flux surface
     call calculate_fourier_integral(mpoloidal,ws1,ws2,ws3,ws4,ws5,ws6,ws7, ws8,ws9,&
          & ws10,ws11,ws12,ws13,ws14,ws15,ws16,ws17,ws18,ws19,ws20,ws21,ws22,ws23,ws24,&
          & kernel0,kernel1,kernel2,kernel3,kernel4,kernel5,kernel6,kernel7,&
          & kernel8,kernel9,kernel10,kernel11,kernel12,kernel13,kernel14,&
          & kernel15,kernel16,kernel17,kernel18,kernel19,kernel20,kernel21,kernel22, &
          & kernel23,kernel24)

     !record matrix on all magnetic surfaces
     do i=1,mhtot
        do j=1,mhtot
           full_kernel0(k,i,j)=kernel0(i,j)
           full_kernel1(k,i,j)=kernel1(i,j)
           full_kernel2(k,i,j)=kernel2(i,j)
           full_kernel3(k,i,j)=kernel3(i,j)
           full_kernel4(k,i,j)=kernel4(i,j)
           full_kernel5(k,i,j)=kernel5(i,j)
           full_kernel6(k,i,j)=kernel6(i,j)
           full_kernel7(k,i,j)=kernel7(i,j)
           full_kernel8(k,i,j)=kernel8(i,j)
           full_kernel9(k,i,j)=kernel9(i,j)
           full_kernel10(k,i,j)=kernel10(i,j)
           full_kernel11(k,i,j)=kernel11(i,j)
           full_kernel12(k,i,j)=kernel12(i,j)
           full_kernel13(k,i,j)=kernel13(i,j)
           full_kernel14(k,i,j)=kernel14(i,j)
           full_kernel15(k,i,j)=kernel15(i,j)
           full_kernel16(k,i,j)=kernel16(i,j)
           full_kernel17(k,i,j)=kernel17(i,j)
           full_kernel18(k,i,j)=kernel18(i,j)
           full_kernel19(k,i,j)=kernel19(i,j)
           full_kernel20(k,i,j)=kernel20(i,j)
           full_kernel21(k,i,j)=kernel21(i,j)
           full_kernel22(k,i,j)=kernel22(i,j)
           full_kernel23(k,i,j)=kernel23(i,j)
           full_kernel24(k,i,j)=kernel24(i,j)
        enddo
     enddo

  enddo
end subroutine calculate_kernel
