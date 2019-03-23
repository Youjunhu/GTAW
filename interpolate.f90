subroutine linear_1d_interpolation(n,x,y,xval,yval)  
  use precision,only:p_
  use constants,only:one
  implicit none
  integer,intent(in):: n
  real(p_),intent(in):: x(n),y(n)
  real(p_),intent(in):: xval
  real(p_),intent(out):: yval

  real(p_):: slope
  integer:: i

  !dx=x(2)-x(1)
  !i=floor(one+(xval-x(1))/dx) !this for uniform x, otherwise we need to call location() subroutine to locate xval
  call location(n,x,xval,i)
  if(i.ge.n) i=n-1

  slope=(y(i+1)-y(i))/(x(i+1)-x(i))
  yval=y(i)+slope*(xval-x(i))

end subroutine linear_1d_interpolation



subroutine location(n,x,xval,k) !use bisection method to locate xval in an array
  use precision,only:p_
  implicit none
  integer,intent(in):: n
  real(p_),intent(in):: x(n),xval
  integer,intent(out)::k
  integer:: kl,ku,km

  kl=1
  ku=n
30 if(ku-kl .gt. 1) then  !use bisection method to search location of theta
     km=(ku+kl)/2
     if((x(n).ge.x(1)).eqv.(xval.ge.x(km))) then
        kl=km
     else
        ku=km
     endif
     goto 30
  endif
  k=kl
end subroutine location
