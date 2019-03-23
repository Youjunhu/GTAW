subroutine write_gnuplot_contour_data(m,n,x,y,u,filename)
  !write data in gnuplot's 'grid data format', so that gnuplot can read it to draw contours of the data
  use precision,only:p_
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: x(m,n),y(m,n),u(m,n)
  character(*),intent(in):: filename
  integer:: i,j
  open(40,file=filename)
  do i=1,m
     do j=1,n
        write(40,*) x(i,j),y(i,j),u(i,j)
     enddo
     ! write a blank line to the file, this blank is needed by gnuplot to recongnize the data as grid data
     write(40,*)          
  enddo
  close(40)
end subroutine write_gnuplot_contour_data


subroutine write_gnuplot_contour_data2(m,n,x,y,u,filename)
  !write data in gnuplot's 'grid data format', so that gnuplot can read it to draw contours of the data
  use precision,only:p_
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: x(m),y(n),u(m,n)
  character(*),intent(in):: filename
  integer:: i,j
  open(40,file=filename)
  do i=1,m
     do j=1,n
        write(40,*) x(i),y(j),u(i,j)
     enddo
     ! write a blank line to the file, this blank is needed by gnuplot to recongnize the data as grid data
     write(40,*)          
  enddo
  close(40)
end subroutine write_gnuplot_contour_data2



subroutine interpolate(x, y, n, new_x, new_y, new_n)
  !linear interpolation
  use precision,only:p_
  use constants,only:zero,one
  implicit none
  integer,intent(in):: n, new_n
  real(p_),intent(in):: x(n),y(n),new_x(new_n)
  real(p_),intent(out):: new_y(new_n)

  real(p_):: weight,tmp
  integer:: flag
  integer:: i,j

  do i=1,new_n

     tmp=new_x(i)

     flag=0
     do j=1,n-1
        if((x(j)-tmp)*(x(j+1)-tmp).le.zero) then
           flag=j
           exit
        endif
     enddo
     if(flag.eq.0) stop 'warning***an element in new_x array is not in the range of x array'
     weight=(tmp-x(flag))/(x(flag+1)-x(flag))
     new_y(i)=(one-weight)*y(flag)+weight*y(flag+1)

  enddo


end subroutine interpolate


subroutine calculate_determinant(matrix_original,nd,z)
  use precision,only:p_
  use constants,only: zero,one,two,twopi
  implicit none
  integer,intent(in):: nd
  real(p_),intent(in):: matrix_original(nd,nd)
  real(p_),intent(out):: z
  
  real(p_)::d
  integer::info,j,ipiv(nd)
  real(p_):: matrix(nd,nd)

  matrix=matrix_original
  ! Compute LU Factorization
  if ( p_ == kind(1.0e1) ) then   !sgetrf and dgetrf are programs in Lapack
     call sgetrf(nd,nd,matrix,nd,ipiv,info)
  else
     call dgetrf(nd,nd,matrix,nd,ipiv,info)
  endif

!  compute determinant
    if ( info .ge. 0 ) then
    d = 1.0_p_
    do j=1,nd
       d = d*matrix(j,j)
       if ( ipiv(j) /= j ) d=-d
    enddo
    z=d
    else
       print *," *** Error in computing determinant"
       print *," info =",info
    endif
end subroutine calculate_determinant

subroutine calculate_determinant_complex(a,n,d)
  use precision ,only: p_
  implicit none
  integer:: n
  complex(p_)::a(n,n),d,s  
  integer::info,j,ipiv(n)

  ! Compute LU Factorization
  if ( p_ == kind(1.0e1) ) then   !cgetrf and zgetrf are programs in Lapack
     call cgetrf(n,n,a,n,ipiv,info)
  else
     call zgetrf(n,n,a,n,ipiv,info)
  endif

  !  compute determinant
  if ( info .ge. 0 ) then
     d = (1.0,0.0)
     do j=1,n
        d = d*a(j,j)
        if ( ipiv(j) /= j ) d=-d
     enddo
  else
     print *," *** Error in computing determinant"
     print *," info =",info
  endif
end subroutine calculate_determinant_complex


subroutine plot_poloidal(mpoloidal,nflux,jacobian,jacobian2,filename)
  use precision,only:p_
  implicit none
  integer,intent(in):: mpoloidal,nflux
  character(len=*),intent(in):: filename
  real(p_),intent(in):: jacobian(mpoloidal,nflux),jacobian2(mpoloidal,nflux)
  integer:: i,j

  open(717,file=filename)
  do j=1,nflux
     do i=1,mpoloidal
        write(717,*) i,jacobian(i,j),jacobian2(i,j)
     enddo
     write(717,*)
     write(717,*)
  enddo
  close(717)

end subroutine plot_poloidal


subroutine plot_radial(mpoloidal,nflux,jacobian,jacobian2,filename)
  use precision,only:p_
  implicit none
  integer,intent(in):: mpoloidal,nflux
  character(len=*),intent(in):: filename
  real(p_),intent(in):: jacobian(mpoloidal,nflux),jacobian2(mpoloidal,nflux)
  integer:: i,j

  open(717,file=filename)
  do i=1,mpoloidal
     do j=1,nflux
        write(717,*) j,jacobian(i,j),jacobian2(i,j)
     enddo
     write(717,*)
     write(717,*)
  enddo
  close(717)

end subroutine plot_radial


subroutine plot_poloidal2(mpoloidal,nflux,jacobian,jacobian2,a,b,c,filename)
  use precision,only:p_
  implicit none
  integer,intent(in):: mpoloidal,nflux
  character(len=*),intent(in):: filename
  real(p_),intent(in):: jacobian(mpoloidal,nflux),jacobian2(mpoloidal,nflux)
  real(p_),intent(in):: a(mpoloidal,nflux),b(mpoloidal,nflux),c(mpoloidal,nflux)
  integer:: i,j

  open(717,file=filename)
  do j=1,nflux
     do i=1,mpoloidal
        write(717,*) i,jacobian(i,j),jacobian2(i,j),a(i,j),b(i,j),c(i,j)
     enddo
     write(717,*)
     write(717,*)
  enddo
  close(717)

end subroutine plot_poloidal2



subroutine plot_poloidal24(mpoloidal,nflux,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,&
     & w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,w23,w24,filename)
  use precision,only:p_
  implicit none
  integer,intent(in):: mpoloidal,nflux
  character(len=*),intent(in):: filename
  real(p_),intent(in):: w1(mpoloidal,nflux),w2(mpoloidal,nflux),w3(mpoloidal,nflux),w4(mpoloidal,nflux),&
       & w5(mpoloidal,nflux),w6(mpoloidal,nflux),w7(mpoloidal,nflux),w8(mpoloidal,nflux),&
       & w9(mpoloidal,nflux),w10(mpoloidal,nflux),w11(mpoloidal,nflux),w12(mpoloidal,nflux),&
       & w13(mpoloidal,nflux),w14(mpoloidal,nflux),w15(mpoloidal,nflux),w16(mpoloidal,nflux),&
       & w17(mpoloidal,nflux),w18(mpoloidal,nflux),w19(mpoloidal,nflux),w20(mpoloidal,nflux),&
       & w21(mpoloidal,nflux),w22(mpoloidal,nflux),w23(mpoloidal,nflux),w24(mpoloidal,nflux)
  integer:: i,j

  open(717,file=filename)
  do j=1,nflux
     do i=1,mpoloidal
        write(717,*) i,w1(i,j),w2(i,j),w3(i,j),w4(i,j),w5(i,j),w6(i,j),w7(i,j),w8(i,j),&
             & w9(i,j),w10(i,j),w11(i,j),w12(i,j),w13(i,j),w14(i,j),w15(i,j),w16(i,j),&
             & w17(i,j),w18(i,j),w19(i,j),w20(i,j),w21(i,j),w22(i,j),w23(i,j),w24(i,j)
     enddo
     write(717,*)
     write(717,*)
  enddo
  close(717)

end subroutine plot_poloidal24



subroutine plot_radial22(mpoloidal,nflux,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,&
     & w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,filename)
  use precision,only:p_
  implicit none
  integer,intent(in):: mpoloidal,nflux
  character(len=*),intent(in):: filename
  real(p_),intent(in):: w1(mpoloidal,nflux),w2(mpoloidal,nflux),w3(mpoloidal,nflux),w4(mpoloidal,nflux),&
       & w5(mpoloidal,nflux),w6(mpoloidal,nflux),w7(mpoloidal,nflux),w8(mpoloidal,nflux),&
       & w9(mpoloidal,nflux),w10(mpoloidal,nflux),w11(mpoloidal,nflux),w12(mpoloidal,nflux),&
       & w13(mpoloidal,nflux),w14(mpoloidal,nflux),w15(mpoloidal,nflux),w16(mpoloidal,nflux),&
       & w17(mpoloidal,nflux),w18(mpoloidal,nflux),w19(mpoloidal,nflux),w20(mpoloidal,nflux),&
       & w21(mpoloidal,nflux),w22(mpoloidal,nflux)
  integer:: i,j

  open(717,file=filename)
  do i=1,mpoloidal
     do j=1,nflux
        write(717,*) j,w1(i,j),w2(i,j),w3(i,j),w4(i,j),w5(i,j),w6(i,j),w7(i,j),w8(i,j),&
             & w9(i,j),w10(i,j),w11(i,j),w12(i,j),w13(i,j),w14(i,j),w15(i,j),w16(i,j),&
             & w17(i,j),w18(i,j),w19(i,j),w20(i,j),w21(i,j),w22(i,j)
     enddo
     write(717,*)
     write(717,*)
  enddo
  close(717)

end subroutine plot_radial22



subroutine wrt_r_z(m,n,r,z,filename1,filename2)
  !plot curves correspond to poloidal_flux=constant and theta=constant.
  !The former corresponds to  magnetic surfaces, the latter corresponds to equal-theta lines
  use precision,only:p_
  implicit none
  integer,intent(in)::m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  character(len=*),intent(in):: filename1,filename2
  integer:: i,j
  open(33,file=filename1)
  do j=1,n
     do i=1,m
        write(33,*) r(i,j),z(i,j)
     enddo
     write(33,*)
     write(33,*)
  enddo
  close(33)

  open(33,file=filename2)
  do i=1,m
     do j=1,n
        write(33,*) r(i,j),z(i,j)
     enddo
     write(33,*)
     write(33,*)
  enddo
  close(33)
end subroutine wrt_r_z

