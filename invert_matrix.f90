!!$program main
!! this main program is for test purpose
!!$  use precision,only:p_
!!$implicit none
!!$
!!$real(p_):: a(2,2),b(2,2)
!!$a(1,1)=1
!!$a(1,2)=2
!!$a(2,1)=3
!!$a(2,2)=4
!!$
!!$call invert_matrix(a,2,b)
!!$write (*,*) b
!!$end program main


subroutine invert_matrix(matrix0,nd,matrix_inverse)
!invert complex matrix0
  use precision,only:p_
  implicit none
  integer,intent(in):: nd
  complex(p_),intent(in):: matrix0(nd,nd)
  complex(p_),intent(out):: matrix_inverse(nd,nd)
  integer::info,j,ipiv(nd)
  complex(p_):: work(10*nd)

  matrix_inverse=matrix0
  !Compute LU Factorization
  if ( p_ == kind(1.0e1) ) then   !sgetrf and dgetrf are programs in Lapack
     call cgetrf(nd,nd,matrix_inverse,nd,ipiv,info)
     call cgetri(nd,matrix_inverse,nd,IPIV, WORK,10*nd,INFO ) 		
  else
     call zgetrf(nd,nd,matrix_inverse,nd,ipiv,info)
     call zgetri(nd,matrix_inverse,nd,IPIV, wORK, 10*nd,INFO ) 		
  endif

end subroutine invert_matrix


subroutine invert_hermite_matrix(matrix0,nd,matrix_inverse)
!invert complex matrix0
  use precision,only:p_
  implicit none
  integer,intent(in):: nd
  complex(p_),intent(in):: matrix0(nd,nd)
  complex(p_),intent(out):: matrix_inverse(nd,nd)
  integer::info,j,ipiv(nd)
  complex(p_):: work(10*nd)

  matrix_inverse=matrix0
  !Compute LU Factorization
  if ( p_ == kind(1.0e1) ) then   !sgetrf and dgetrf are programs in Lapack
     call chetrf(nd,nd,matrix_inverse,nd,ipiv,info)
     call chetri(nd,matrix_inverse,nd,IPIV, WORK,10*nd,INFO ) 		
  else
     call zhetrf(nd,nd,matrix_inverse,nd,ipiv,info)
     call zhetri(nd,matrix_inverse,nd,IPIV, wORK, 10*nd,INFO ) 		
  endif

end subroutine invert_hermite_matrix






