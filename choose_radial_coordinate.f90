subroutine choose_radial_coordinate(psi_axis,psi_lcfs,psival_nx,tfn_nx,qpsi_nx,nx,total_tf,nflux,psival_new,ra,dra,dpsidra)
  use constants,only: one,two,twopi
  use precision,only:p_
  implicit none
  real(p_),intent(in)::psi_axis,psi_lcfs
  integer,intent(in)::nx,nflux
  real(p_),intent(in)::psival_nx(nx),qpsi_nx(nx),tfn_nx(nx),total_tf
  real(p_),intent(out):: psival_new(nflux),dpsidra(nflux),ra(nflux),dra !!dpsidr is the differential of GS psi (Aphi*R) with respct to the radial coordinate
  real(p_):: qpsi(nflux),psival_nx_t(nx),qpsi_nx_t(nx)
  real(p_)::inner_pfn !to avoid singularity at the magnetic axis, the first flux surface is chosen at psi=psi_axis+delta_psi, where delta_psi=(psi_lcfs-psi_axis)*inner_pfn
  real(p_):: bdry_pfn !the radial boundary is choosen at psi_bdry= psi_axis+(psi_lcfs-psi_axis)*bdry_pfn
  real(p_):: sqrt_ntf_in !sqrt_ntf_in is the value of the sqrt_normalized_toroidal_flux of the first radial location
  real(p_):: sqrt_ntf_bdry !sqrt_ntf_bdry is the value of the sqrt_normalized_toroidal_flux of the boundary radial location
  character(100):: radial_coordinate_type
  ! real(p_):: psi_in,psi_bdry !psi_axis and psi_lcfs are respectively the value of psi at magnetic axis and last-colsed-flux-surface (LCFS), !psi_in is the value of psi on the flux surface that is very near the magnetic axis
  real(p_):: sqrt_tfn_nx(nx),tmp_y2(nx),q_y2(nx),y_tmp
  !  real(p_):: uniform_sqrt_tfn(nflux)
  integer:: j

  namelist/radial_coordinate/radial_coordinate_type,inner_pfn,bdry_pfn,sqrt_ntf_in,sqrt_ntf_bdry

  open(111,file='gtaw.in')
  read(111,radial_coordinate)
  close (111)
  write(*,radial_coordinate)

  if(trim(radial_coordinate_type).eq.'gs_pf') then
     do j=1,nflux 
        ra(j)=inner_pfn+(bdry_pfn-inner_pfn)/(nflux-1)*(j-1)
     enddo
     ra=psi_axis+ra*(psi_lcfs-psi_axis) !conerted to GS poloidal flux funtion psi
     dra=ra(2)-ra(1) !interval of psi between two neighbour magnetic surfaces, grid interval, uniform grid is assumed
     psival_new=ra 
     dpsidra=one !dpsidr is the differential of GS psi (Aphi*R) with respct to the radial coordinate

  else if(trim(radial_coordinate_type).eq.'pfn') then
     do j=1,nflux 
        ra(j)=inner_pfn+(bdry_pfn-inner_pfn)/(nflux-1)*(j-1)
     enddo
     dra=ra(2)-ra(1) !interval of psi between two neighbour magnetic surfaces, grid interval, uniform grid is assumed
     psival_new=psi_axis+ra*(psi_lcfs-psi_axis) !set the GS poloidal flux vlaue
     do j=1,nflux
        dpsidra(j)=psi_lcfs-psi_axis !dpsidr is the differential of GS psi (Aphi*R) with respct to the radial coordinate
     enddo

!!$do j=1,nflux
!!$write(*,*) j, dpsidra(j), 'choose'
!!$enddo

  else if  (trim(radial_coordinate_type).eq.'tfn_sqrt') then
     do j=1,nflux !uniform grid in sqrt(toroidal_flux), used in analysing the data of MEGA code
        ra(j)=sqrt_ntf_in+(sqrt_ntf_bdry-sqrt_ntf_in)/(nflux-1)*(j-1)
     enddo
     dra=(ra(nflux)-ra(1))/(nflux-1)
     call spline(sqrt(tfn_nx),psival_nx,nx,2.d30,2.d30,tmp_y2) !prepare the second order derivative needed in the cubic spline interpolation
     do j=1,nflux 
        call splint(sqrt(tfn_nx),psival_nx,tmp_y2,nx,ra(j),psival_new(j)) !to get the value of gs_psi on the new radial grids
     enddo

     !the spline interpolating works only for sequence x1<x2<x3...<xn. The following step is to gurantee that psival_nx_t is in this order.
     if(psival_nx(nx)<psival_nx(1)) then
        do j=1,nx
           psival_nx_t(j)=psival_nx(nx+1-j)
           qpsi_nx_t(j)=qpsi_nx(nx+1-j)
        enddo
     else
        psival_nx_t=psival_nx
        qpsi_nx_t=qpsi_nx
     endif

     call spline(psival_nx_t,qpsi_nx_t,nx,2.d30,2.d30,q_y2) !prepare the second order derivative needed in the cubic spline interpolation
     do j=1,nflux
        call splint(psival_nx_t,qpsi_nx_t,q_y2,nx,psival_new(j),y_tmp) !to get qpsi array corresponding to the new radial grid points
        qpsi(j)=y_tmp 
        !write(*,*) j, qpsi(j), 'yj'
     enddo

     do j=1,nflux
        dpsidra(j)=one/(twopi*qpsi(j))*total_tf*two*ra(j)
     enddo
  else
     stop 'please choose correct radial_coordinate_type: gs_pf, pfn, or tfn_sqrt'
  endif
  ! stop 'ddddddddddd'
end subroutine choose_radial_coordinate
