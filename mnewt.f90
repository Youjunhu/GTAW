SUBROUTINE mnewt(ntrial,x,n,tolx,tolf)
!modified from codes in Numerical Recipes book
  use precision,only: p_
  implicit none
  integer,intent(in):: ntrial,n
  complex(p_),intent(inout):: x(n)
  real(p_),intent(in):: tolx,tolf
  !USES lubksb,ludcmp,usrfun
  INTEGER i,k,indx(n),info,ipiv(n)
  complex(p_):: fjac(n,n),fvec(n),p(n)
  real(p_)::errf,errx

    write(*,'(a10,6a16)') 'root=', 'real(x(1))','imag(x(1))', 'real(x(n))','imag(x(n))', 'errf','errx'
  do   k=1,ntrial
     call usrfun(x,n,fvec,fjac) !to calculate values of functions and Jacobian matrix
     write(*,'(a10,4E16.8)',advance='no') 'root=', real(x(1)),imag(x(1)), real(x(n)),imag(x(n))

     errf=0._p_
     do  i=1,n
        errf=errf+abs(fvec(i))
     enddo
     write(*,'(E16.8)',advance='no') errf
     if(errf.le.tolf) return

     do  i=1,n
        p(i)=-fvec(i)
     enddo
     !call ludcmp(fjac,n,n,indx,d)
     !call lubksb(fjac,n,n,indx,p)
     !call dgesv(n,1,fjac,n,ipiv,p,n,info) !lapack routine to slove linear equations system
     call zgesv(n,1,fjac,n,ipiv,p,n,info) !lapack routine to slove linear equations system
     errx=0._p_
     do  i=1,n
        errx=errx+abs(p(i))
        x(i)=x(i)+p(i)
     enddo
     write(*,'(E16.8)')  errx
     if(errx.le.tolx) return
  enddo

  write(*,*) '******maximum iteration number exceeded********'
END SUBROUTINE mnewt

 

subroutine usrfun(x,n,fvec,fjac)
  use precision,only: p_
  implicit none
  integer,intent(in):: n
  complex(p_),intent(in):: x(n)
  complex(p_),intent(out):: fvec(n),fjac(n,n)
  call funcv(n,x,fvec)
  call fdjac(n,x,fvec,fjac) !calculate the Jacobian matrix of func
end subroutine usrfun


SUBROUTINE fdjac(n,x,fvec,fjac) !calculate the Jacobian matrix of func
  !modified from codes in Numerical Recipes book
  use precision,only:p_
  integer,intent(in):: n
  complex(p_),intent(in):: fvec(n)
  complex(p_),intent(out):: fjac(n,n)
  complex(p_):: x(n)
!complex(p_):: EPS=(1.d-5,1.d-5)
  complex(p_):: EPS=(1.d-4,1.d-4)
  !complex(p_):: EPS=(1.d-3,1.d-3)
  !U    USES funcv
  INTEGER i,j
  complex(p_) h,temp,f(n)

  do j=1,n
     temp=x(j)
     h=EPS*abs(temp)
     if(abs(h).eq.0.) h=EPS
     x(j)=temp+h
     h=x(j)-temp !determine the actuall interval
     call funcv(n,x,f)
     x(j)=temp
     do i=1,n
        fjac(i,j)=(f(i)-fvec(i))/h
     enddo
  enddo
END SUBROUTINE fdjac


