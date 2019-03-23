subroutine normalized_parallel_current_density(mpoloidal,nflux,fpsi,ffprime,b0_axis,pprime_normalized,bsq,sigma_mu0)
  !calculate sigma_mu0, which is defined by sigma_mu0=mu0*J_dot_B/B^2, where mu0 the permeability in SI unit
  use precision,only:p_
  use constants,only: two !,mu0
  implicit none
  integer,intent(in):: mpoloidal,nflux
  real(p_),intent(in):: fpsi(nflux),ffprime(nflux),pprime_normalized(nflux),bsq(mpoloidal,nflux),b0_axis
  real(p_),intent(out):: sigma_mu0(mpoloidal,nflux)
  integer:: i,j
  real(p_):: fprime(nflux)

  fprime=ffprime/fpsi

  do i=1,mpoloidal
     do j=1,nflux
        sigma_mu0(i,j)=fpsi(j)*pprime_normalized(j)/bsq(i,j)*b0_axis**2/two+fprime(j) !refer to the formula in /home/yj/theory/tokamak_equilibrium.tm
     enddo
  enddo
end subroutine normalized_parallel_current_density

