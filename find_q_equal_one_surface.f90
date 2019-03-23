subroutine find_q_equal_one_surface(tfn_nx,psival_nx,pfn_nx,qpsi_nx,nx,x_lcfs,z_lcfs,np_lcfs,x_axis,z_axis)
  use precision,only:p_
  use constants,only:zero,one,two,twopi
  implicit none

  integer,intent(in):: nx,np_lcfs
  real(p_),intent(in):: psival_nx(nx),pfn_nx(nx),qpsi_nx(nx),tfn_nx(nx)
  real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs),x_axis,z_axis
  real(p_)::  tmp_y2(nx)
  real(p_),parameter:: xacc=1.0d-6 !tolerance used in bi-section root-finder
  real(p_):: pfn_val,psival,tfn_val
  real(p_):: rtbis2,q_func !function name
  external:: q_minus_one_func
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  integer:: i

  call spline(pfn_nx,qpsi_nx,nx,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation, after this call, the function q_func is ready to be used.
  pfn_val=rtbis2(q_minus_one_func,pfn_nx(1),pfn_nx(nx),xacc,nx,pfn_nx,qpsi_nx,tmp_y2)

  write(*,*) 'The normalized poloidal magnetic flux at q=1 surface is', pfn_val

  !find the corresponding normalized toroidal flux
  call spline(pfn_nx,tfn_nx,nx,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
  call splint(pfn_nx,tfn_nx,tmp_y2,nx,pfn_val,tfn_val) 
  write(*,*) 'the normalized toroidal flux at q=1 surface is', tfn_val, 'sqrt of normalized toroidal flux is ', sqrt(tfn_val)

  psival=psival_nx(1)+pfn_val*(psival_nx(nx)-psival_nx(1))
  call contour(x_lcfs,z_lcfs,np_lcfs,x_axis,z_axis,psival,x_contour,z_contour)
  open(398,file='q1_surface')
  do i=1,np_lcfs
     write(398,*) x_contour(i),z_contour(i)
  enddo
  close(398)

end subroutine find_q_equal_one_surface



function q_minus_one_func(pfn_val,nx,pfn_nx,qpsi_nx,tmp_y2)
  use precision,only:p_
  implicit none
  real(p_):: q_minus_one_func
  integer,intent(in):: nx
  real(p_),intent(in):: pfn_val,pfn_nx(nx),qpsi_nx(nx), tmp_y2(nx)
  real(p_):: tmp_y

  call splint(pfn_nx,qpsi_nx,tmp_y2,nx,pfn_val,tmp_y) 
  !q_minus_one_func=tmp_y-1.5_p_
  if(qpsi_nx(nx/2)<0._p_) then
     q_minus_one_func=tmp_y-(-1._p_)
  else
     q_minus_one_func=tmp_y-1._p_
  endif

end function q_minus_one_func




