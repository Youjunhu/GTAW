subroutine calculate_continous_spectrum(nflux,nh,r_axis,b0_axis,wa0)
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  use poloidal_harmonics,only: mhtot
  use flux_grids,only: mpoloidal
  use radial_module,only:pfn,psival_new,rho_normalized,pressure_normalized,tfn,qpsi,dpsidra
  use kernel_module,only: full_kernel0,full_kernel1,full_kernel2,full_kernel3,full_kernel4,full_kernel5 !Fourier integrations on every flux surface
  use kernel_module,only: full_kernel6,full_kernel7 !as input
  implicit none
  integer,intent(in):: nflux,nh
  real(p_),intent(in):: r_axis,b0_axis,wa0
  
  complex(p_):: kernel0(mhtot,mhtot),kernel1(mhtot,mhtot),kernel2(mhtot,mhtot),kernel3(mhtot,mhtot)
  complex(p_):: kernel4(mhtot,mhtot),kernel5(mhtot,mhtot),kernel6(mhtot,mhtot),kernel7(mhtot,mhtot)
  complex(p_):: matrix_ea(2*mhtot,2*mhtot),matrix_eb(2*mhtot,2*mhtot)
  integer,parameter::  nn=2*mhtot
  complex(p_):: omega_sq(nn) !omega_sq=(w/wa0)**2, where w is the frequency of the continua  
  integer:: k,i,j

  open(15,file='continua.txt') !record the continua
  !open(16,file='continua1.txt') !record the continua in the form that I can plot continua as lines, instead of dots. However this turns out to be not feasible because the sequency of the array real(omega_sq(i=1,nn)) on different magnetic surfaces are different.

  do k=1,nflux
     !fetch the kernels
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
        enddo
     enddo

     call calculate_matrix_e(nh,r_axis,b0_axis,rho_normalized(k),pressure_normalized(k),qpsi(k), dpsidra(k),&
          & kernel0,kernel1,kernel2,kernel3,kernel4,kernel5,kernel6,kernel7,matrix_ea,matrix_eb) !calculate matrix E 

     call solve_generalized_eigen_value_problem(nn,matrix_ea,matrix_eb,omega_sq)

     do i=1,nn !record the continua. Note that there are nn eigenvalues on every flux surface
        write(15,*) sqrt(pfn(k)),psival_new(k),sqrt(tfn(k)), & !three different choices of the radial coordinate, psival_new is the poloidal magnetic flux (Tm^2) given in gfiel
             & real(omega_sq(i)),imag(omega_sq(i)), & !this is (w*R0/VA)**2
             & sqrt(abs(real(omega_sq(i))))*wa0/twopi/1000._p_, & !this is the frequency in kHz
             & sqrt(abs(real(omega_sq(i)))) !this is w*R0/VA
     enddo
     !write(16,*) pfn(k),  (sqrt(abs(real(omega_sq(i))))*wa0/twopi/1000._p_,i=1,nn)  !this is frequency in kHz
  enddo

  close(15)
  write(*,*) '>>>>>Finish calculating MHD continua.'

end subroutine calculate_continous_spectrum



subroutine solve_generalized_eigen_value_problem(nn,matrix_a,matrix_b,omega_sq)
  use constants,only:p_
  implicit none
  integer,intent(in):: nn
  complex(p_),intent(in)::matrix_a(nn,nn),matrix_b(nn,nn)
  complex(p_),intent(out):: omega_sq(nn)
  !--work array for lapack routine--
  integer,parameter:: lwork=10000
  complex(p_):: alfa(nn),beta(nn)
  complex(p_):: vl(nn,nn),vr(nn,nn)
  complex(p_):: work(lwork),rwork(8*nn)
  integer:: info,i

  !--the following is the professional way of calculating the continua, i.e. put the problem as a generalized eigenvalue problem matrix_a*x=-wsq**2*matrix_b*x and solve it by using the robust lapack subroutine:
  if (p_.eq.kind(1.d0)) then
     call zggev( 'n', 'n', nn,matrix_a,nn,matrix_b,nn,alfa,beta,vl,nn,vr,nn,work,lwork,rwork,info)
  else
     call cggev( 'n', 'n', nn,matrix_a,nn,matrix_b,nn,alfa,beta,vl,nn,vr,nn,work,lwork,rwork,info)
  endif

  do i=1,nn !record the eigenvalues, which is (w/(VA0/R0))**2
     omega_sq(i)=-alfa(i)/beta(i) !which is a complex number, however the results indicate that the imaginary part is small, which is consistent with the analytic conclusion that omegasq must be a real number
  enddo
  !-----------lapack routine end------------
end subroutine solve_generalized_eigen_value_problem

