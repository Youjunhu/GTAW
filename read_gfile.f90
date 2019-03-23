subroutine read_gfile(nwmax,psi,nx,nz,rleft,zmid,xdim,zdim,psi_axis,psi_lcfs,xmaxis,zmaxis, &
     & x_lcfs,z_lcfs,np_lcfs,fpsi,qpsi,presspsi,pprime,ffprime,current)
  use precision,only:p_
  implicit none
  integer,intent(in)::nwmax
  real(p_),intent(out):: psi(nwmax,nwmax)
  integer,intent(out):: nx,nz
  character*8:: ntitle(5),vid
  integer:: neq,ipestg
  real(p_)::xdim, zdim, rmajor_mk, rleft, zmid
  real(p_)::xmaxis, zmaxis, psi_axis, psi_lcfs, btorus
  real(p_):: dumaraya4(4),dumarayb5(5) 
  real(p_),intent(out):: fpsi(nwmax),qpsi(nwmax),presspsi(nwmax)
  real(p_),intent(out):: ffprime(nwmax),pprime(nwmax)
  real(p_),intent(out):: current !current the is the total toroidal current within the LCFS
  integer::np_lcfs, nlim_eqd
  real(p_):: x_lcfs(nwmax),z_lcfs(nwmax),rlim_eqd(nwmax), zlim_eqd(nwmax)
  integer:: i,j
  character(80):: gfile_name,gfile_name_out

  namelist/gfile_namelist/gfile_name
  open(11,file='gtaw.in')
  read(11,gfile_namelist)
  close(11)
  write(*,gfile_namelist)

  neq=123
  !open and read in eqdsk file
  open(neq,file=gfile_name,status='old')

  ipestg = 4
  read (neq, '(6a8, 3i4)') (ntitle(i), i=1,5), vid, ipestg, nx, nz
  read (neq,300) xdim, zdim, rmajor_mk, rleft, zmid
  read (neq,300) xmaxis, zmaxis, psi_axis, psi_lcfs, btorus
  read (neq,300) current,dumaraya4
  read (neq,300) dumarayb5
  read (neq ,300) (fpsi(i), i=1,nx)
  read (neq ,300) (presspsi(i), i=1,nx)
  read (neq ,300) (ffprime(i), i=1,nx)
  read (neq ,300) (pprime(i), i=1,nx)
  read (neq ,300) ((psi(i,j), i=1,nx), j=1,nz)
  read (neq ,300) (qpsi(i), i=1,nx)
  read (neq ,'(2i5)') np_lcfs, nlim_eqd
  read (neq ,300) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
  read (neq ,300) (rlim_eqd(i), zlim_eqd(i), i=1,nlim_eqd)
  close(neq)

  !write(*,*) 'toroidal current given in G-file is (kA)', current/(1d3)
  !to verify that I have read the eqdsk file correctly, I write the data read to a new file.
  !After the program finished, I compare this file with the original file using diff command
  !the output of diff command indicates that the two files are idential, which shows I have read the eqdsk file correctly
  !Somtimes, I alter some quantities (e.g. increase the pressure by a constant),in this case, the out gfile is different from the original one
!!$  do i=1,nx
!!$     presspsi(i)=presspsi(i)+0.5*presspsi(1) !increase the presure
!!$  enddo
 

 !fpsi=-fpsi !revert the toroidal magnetic field


!!$psi=-psi !revert direction of the torodial current
!!$psi_axis=-psi_axis
!!$psi_lcfs=-psi_lcfs
!!$ffprime=-ffprime
!!$pprime=-pprime


! fpsi=sign(1._p_,fpsi(1))*sqrt(fpsi**2-0.5*abs(fpsi(1))) !change the toroidal magnetic field, to change the q profile

  neq=111
  gfile_name_out=trim(gfile_name)//'_GTAW_output'
  open(neq,file=gfile_name_out)
  write (neq, '(6a8, 3i4)') (ntitle(i), i=1,5), vid, ipestg, nx, nz
  write (neq,300) xdim, zdim, rmajor_mk, rleft, zmid
  write (neq,300) xmaxis, zmaxis, psi_axis, psi_lcfs, btorus !note that btorus is the vacuum magnetic field at rmajor_mk, which is not at the magnetic axis and is usually different the actual magnetic field in plasma ! Futher, I find that this value given in gfile of EAST48916@4.5s is incorrect
  write (neq,300) current, dumaraya4
  write (neq,300) dumarayb5
  write (neq ,300) (fpsi(i), i=1,nx)
  write (neq ,300) (presspsi(i), i=1,nx)
  write (neq ,300) (ffprime(i), i=1,nx)
  write (neq ,300) (pprime(i), i=1,nx)
  write (neq ,300) ((psi(i,j), i=1,nx), j=1,nz)
  write (neq ,300) (qpsi(i), i=1,nx)
  write (neq ,'(2i5)') np_lcfs, nlim_eqd
  write (neq ,300) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
  write (neq ,300) (rlim_eqd(i), zlim_eqd(i), i=1,nlim_eqd)
  close(neq)

300 format (5e16.9)

  write(*,*) 'Computational box used in g-file ','Rleft=',rleft,'Rlength=',xdim, 'Zlength=',zdim, 'Zmid=',zmid
  write(*,*) 'dx (m)=',xdim/(nx-1),'dz (m)=',zdim/(nz-1)
  write(*,*) 'Magnetic location ', 'r_axis=',xmaxis, 'z_axis=',zmaxis
  write(*,*) 'psi_axis=',psi_axis,'psi_lcfs=',psi_lcfs
  write(*,*) 'number of points used in specifying the shape of LCFS, np_lcfs=',np_lcfs
  !write(*,*) 'x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)=', x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)


!stop 'stop in read gfile'
end subroutine read_gfile
