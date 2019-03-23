subroutine calculate_radial_matrix(mpoloidal,nflux,mhtot)
  use precision,only:p_
  use radial_module,only: dpsidra
  use kernel_module,only:full_kernel13 !as input
  use kernel_module,only: full_radial_matrix !as outpu
  implicit none
  integer,intent(in):: mpoloidal,nflux,mhtot
  complex(p_):: matrix(mhtot,mhtot),matrix_inverse(mhtot,mhtot)
  !complex(p_):: full_matrix(2*mhtot,2*mhtot),full_matrix_inverse(2*mhtot,2*mhtot)
  integer:: i,j,k

  do k=1,nflux

     do i=1,mhtot
        do j=1,mhtot
           matrix(i,j)=full_kernel13(k,i,j)/dpsidra(k)
        enddo
     enddo

     call invert_matrix(matrix,mhtot,matrix_inverse)

     do i=1,mhtot
        do j=1,mhtot
           full_radial_matrix(k,i,j)=matrix_inverse(i,j)
           full_radial_matrix(k,i,j+mhtot)=0.
           full_radial_matrix(k,i+mhtot,j)=0.
           full_radial_matrix(k,i+mhtot,j+mhtot)=matrix_inverse(i,j)
        enddo
     enddo
!-----another way---

!!$   do i=1,mhtot
!!$        do j=1,mhtot
!!$           full_matrix(i,j)=matrix(i,j)
!!$           full_matrix(i,j+mhtot)=0.
!!$           full_matrix(i+mhtot,j)=0.
!!$           full_matrix(i+mhtot,j+mhtot)=matrix(i,j)
!!$        enddo
!!$     enddo
!!$
!!$   call invert_matrix(full_matrix,2*mhtot,full_matrix_inverse)
!!$
!!$   do i=1,2*mhtot
!!$        do j=1,2*mhtot
!!$           full_radial_matrix(k,i,j)=full_matrix_inverse(i,j)
!!$        enddo
!!$     enddo

  enddo

end subroutine calculate_radial_matrix
