subroutine set_radial_range(rad_left,rad_right)
  use precision,only:p_
  use constants,only: one,two,three
  use flux_grids,only: nflux
  use radial_module,only: tfn
  use radial_module,only:starting_surface_number,ending_surface_number !radial grid points, as output
  implicit none
  real(p_),intent(in):: rad_left,rad_right
  real(p_)::tf_sq(nflux)
  integer:: k
   
  tf_sq=sqrt(tfn)
  call location(nflux,tf_sq,rad_left,k) !use bisection method to locate xval in an array
  starting_surface_number=k+1
  call location(nflux,tf_sq,rad_right,k) !use bisection method to locate xval in an array
  ending_surface_number=k

  !starting_surface_number=1 !test
  write(*,*) 'radial region used in searching for global modes: ', 'starting_surface_number ',&
       & starting_surface_number,'ending_surface_number',ending_surface_number
end subroutine set_radial_range
