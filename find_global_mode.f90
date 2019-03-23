subroutine find_global_mode(nh,mhtot,nflux,wa0)
  !use shooting method to find eigenvalues and the corresponding eigenfunctions
  use precision,only:p_
  use constants,only: one,two,three,twopi,ii !ii=(0._p_,1.0_p_)
  use radial_module,only: starting_surface_number,ending_surface_number,qpsi!radial grid points
  use flux_grids,only: mpoloidal
  implicit none
  integer,intent(in):: nh,mhtot,nflux
  real(p_),intent(in):: wa0 !wa0 is the normalizing frequency

  complex(p_)::ypath(2*mhtot+1,nflux),ypath2(2*mhtot,nflux),xir2d(mpoloidal,nflux)  !eigenfunctions
  real(p_):: phase(mhtot,nflux),phase_total(mpoloidal,nflux),phase_midplane(nflux) !radial phase for poloidal harmonics of the radial displacement
  complex(p_):: xir_midplane(mhtot,nflux)
  real(p_):: xpath(nflux)

  integer:: n2
  integer:: i,j,k
  !real(p_):: v_tmp(2*mhtot),f_tmp(2*mhtot) !independent variable and function used in shooting method
  complex(p_):: v_tmp(mhtot),f_tmp(mhtot) !independent variable and function used in shooting method
  logical:: check
  real(p_):: guess_freq,frequency_khz,growth_khz,guess_freqval
  logical:: restart_from_last
  logical:: cond1,cond2,cond3,cond4
  integer:: interation_n !number of interation used in the newton root finder
  integer:: mh1,mh2
  real(p_):: rad_left,rad_right !the radial region used in searching global modes, in unit of the sqrt(normalized_toroidal_flux)
  namelist/global_mode_parameters/rad_left,rad_right,restart_from_last,guess_freq,interation_n !rad_left is in unit of sqrt_normlaized_toroidal_flux

  open(11,file='gtaw.in')
  read(11,global_mode_parameters)
  close(11)  

  !---for n=1, m=1,2 TAEs 38300------
!!$  starting_surface_number=3
!!$  ending_surface_number=80
  !for g048916.04500
  ! starting_surface_number=3
  !  ending_surface_number=nflux-50
  call set_radial_range(rad_left,rad_right) !determine the sequence number of the innermost and outmost surface
  !use the full radial range, In this case, I use less poloidal harmonics (mh_low=-3 and mh_upp=8), so that the results are not seriously blured
!!$  starting_surface_number=2
!!$  ending_surface_number=nflux
  !----end-----------------
  !---for n=1, m=2,3 TAEs------
!!$  starting_surface_number=61
!!$  ending_surface_number=nflux-25
  !-------end--------------
  !--for n=1,m=3,4 TAEs--
!!$  starting_surface_number=58
!!$  ending_surface_number=nflux
  !----end-------

  !n2=2*mhtot !number of independent variables (and is also the number of nonlinear algebra equations), the 2 factor is due to the variables are complex numbers, which need two real numbers to represent one variable.
  n2=mhtot !using complex numbers

!!$  do i=1,2*(mhtot-1) !this sets the guess values of the first mhtot-1 harmonics of perturbed pressure at the left boundary. Since two real numbers are needed to determine one complex number, the total number of real numbers needed is 2*(mhtot-1)
!!$     v_tmp(i)= 2.0d-3
!!$  enddo
!!$
!!$  v_tmp(2*mhtot-1)=25._p_*1000._p_*twopi/wa0  !guess value of real part of wsq
!!$  v_tmp(2*mhtot)=3.0790701597008d-8    !guess value of imaginary part of wsq
!!$
!!$  do k=1,40
!!$     do i=1,n2
!!$        v_tmp(i)=i*0.5+6.0_p_*k*10
!!$     enddo
!!$     call funcv(n2,v_tmp,f_tmp)
!!$
!!$     write(*,*) "f_tmp=",(f_tmp(j),j=1,5)
!!$  enddo

!!$  call newt(v_tmp,n2,check)
!!$  if(check)then
!!$     write(*,*)'shoot failed; bad initial guess'
!!$  else
!!$     write(*,*) 'shoot succeed,',' eigenvalue real part',v_tmp(2*mhtot-1),'imaginary part',v_tmp(2*mhtot)
!!$     write(*,*)  'Mode frequency (KHz)=', sqrt(v_tmp(2*mhtot-1))*wa0/twopi/1000._p_
!!$  endif
  write(*,*) 'begin to search eigenmodes...'
  open(1238,file='possible_modes.txt') !select modes we are interested in

  do k=1,50
     write(*,*) 'restart_from_last=',restart_from_last
     if(restart_from_last.eqv. .true.) then
        open(111,file='root_guess.txt')
        do i=1,mhtot
           read(111,*) v_tmp(i)
        enddo
        close(111)
     else
        do i=1,mhtot-1 !this sets the guess values of the first mhtot-1 harmonics of perturbed pressure at the left boundary. 
           v_tmp(i)= (1.0d-6,1.0d-6)
        enddo
        !v_tmp(mhtot)=((101.+k*1)*1000._p_*twopi/wa0)**2+ii*(0.5_p_*1000._p_*twopi/wa0)**2 !this is the guess value of w^2, which is a complex number
        !v_tmp(mhtot)=((296+k*1)*1000._p_*twopi/wa0)**2+ii*(0.5_p_*1000._p_*twopi/wa0)**2 !this is the guess value of w^2, which is a complex number
        guess_freqval=guess_freq+(k-1)*0.1_p_ !in kHz
        v_tmp(mhtot)=(guess_freqval*1000._p_*twopi/wa0)**2+ii*(0.1_p_*1000._p_*twopi/wa0)**2 !this is the guess value of w^2/wa0^2, which is a complex number
        write(*,*) 'guess frequency=', guess_freqval
     endif
     call mnewt(interation_n,v_tmp,n2,1.0d-12,1.0d-12)

     open(111,file='root_guess.txt') !store the roots, which will be used as intial guess of roots when refining the roots
     do i=1,mhtot
        write(111,*) v_tmp(i)
     enddo
     close(111)

     if(real(v_tmp(mhtot))<0.) write(*,*) '***warning**** wsq<0'
     frequency_khz= sqrt(real(v_tmp(mhtot)))*wa0/twopi/1000._p_
     growth_khz= sqrt(abs(imag(v_tmp(mhtot))))*wa0/twopi/1000._p_
     write(*,*)'----------------','k=',k,'---------------------------------------'
     write(*,*) 'Mode frequency (KHz)=', frequency_khz, 'Imaginary part of the root (KHz)=', growth_khz

     !     write(*,*) '(Imaginary part of a reasonable root should be much smaller than the real part)'
     write(*,*) '(q(a)*w/wa)**2=',real(v_tmp(mhtot))*qpsi(nflux)**2 !this is the normalized wsq output by NOVA code, for the purpose of comparison

     call calculate_eigenfunction(n2,v_tmp,xpath,ypath,ypath2,xir2d,phase,phase_total,phase_midplane,xir_midplane) !calculate the eigen function and store them in ypath arrays
     !call determine_dominant_modes(mhtot,nflux,xpath,ypath,starting_surface_number,ending_surface_number,nh,wa0) !determine the dominant poloidal harmonics
     call determine_dominant_modes(mhtot,nflux,xpath,ypath,nh,wa0,mh1,mh2) !determine the dominant poloidal harmonics

cond1=((mh1.eq.1) .and. (mh2.eq.2)).or. ((mh1.eq.2) .and. (mh2.eq.1))
cond2=((mh1.eq.2) .and. (mh2.eq.3)).or. ((mh1.eq.3) .and. (mh2.eq.2))
cond3=((mh1.eq.3) .and. (mh2.eq.4)).or. ((mh1.eq.4) .and. (mh2.eq.3))


     if(cond1 .or. cond2 .or. cond3) then
        write(1238,*) 'k=',k, 'mode freq.=', frequency_khz , 'mh1=',mh1,'mh2=',mh2
     endif
     call record_eigenfunctions(k,mhtot,nflux,frequency_khz,growth_khz,xpath,ypath,ypath2,xir2d,phase,&
          & phase_total,phase_midplane,xir_midplane) !write the eigenfunctions in data files

  enddo
  close(1238)
end subroutine find_global_mode


subroutine determine_dominant_modes(mhtot,nflux,xpath,ypath,nh,wa0,mh1,mh2)
  use precision,only:p_
  use constants,only: two,one,three,twopi
  use poloidal_harmonics,only: mh_low
  use radial_module,only:tfn,starting_surface_number,ending_surface_number
  implicit none
  integer,intent(in):: mhtot,nflux,nh
  real(p_),intent(in):: xpath(nflux),wa0
  complex(p_),intent(in):: ypath(2*mhtot+1,nflux)
  integer:: mharmonic
  real(p_)::max_displacement
  integer:: max_displacement_loc(2)
  integer:: i,j
  integer,intent(out):: mh1,mh2
  real(p_):: displacement(mhtot,nflux)

 displacement=0.

  do i=1,mhtot
     do j=starting_surface_number,ending_surface_number
        displacement(i,j)=abs(ypath(i+mhtot,j))
     enddo
  enddo

  max_displacement=maxval(displacement) !maxval is Fortran intrinsic function
  !write(*,*) 'max_displacement= ', max_displacement
  max_displacement_loc=maxloc(displacement) !maxloc is Fortran intrinsic function, which returns the index of the largerst elements of a array
  mh1=mh_low+max_displacement_loc(1)-1
  write(*,*) 'Largest poloidal mode number determined from the amplitude of plasma displacement:'
  write(*,*) 'polodial_mode_number=',mh1
  write(*,*) 'peak_location (sqrt(toroidal_flux_normalized))=',sqrt(tfn(max_displacement_loc(2)))

  do j=1,nflux !zero the largest poloidal harmonics so that then I can select the second largest poloidal harmonics
     displacement(max_displacement_loc(1),j)=0.
  enddo

  max_displacement=maxval(displacement) !maxval is Fortran intrinsic function
  !write(*,*) 'max_displacement= ', max_displacement

  max_displacement_loc=maxloc(displacement) !maxloc is Fortran intrinsic function, which returns the index of the largerst elements of a array
  mh2=mh_low+max_displacement_loc(1)-1
  write(*,*) 'Second largest poloidal mode number determined from the amplitude of plasma displacement:'
  write(*,*) 'poloidal_mode_number=',mh2
  write(*,*) 'peak_location (sqrt(toroidal_flux_normalized))=',sqrt(tfn(max_displacement_loc(2)))

  mharmonic=min(mh1,mh2) !return the smaller one of the two domenant mode numbers

  write(*,*) 'approximate centre frequency (KHZ) of the TAE gap ', wa0*abs(nh/    (two*mharmonic+one))/twopi/1000._p_
  write(*,*) 'approximate centre frequency (KHZ) of the EAE gap ', wa0*two*abs(nh/(two*mharmonic+two))/twopi/1000._p_
  write(*,*) 'approximate centre frequency (KHZ) of the NAE gap ', wa0*three*abs(nh/(two*mharmonic+three))/twopi/1000._p_

end subroutine determine_dominant_modes



subroutine record_eigenfunctions(k,mhtot,nflux,frequency_khz,growth_khz,xpath, &
     & ypath,ypath2,xir2d,phase,phase_total,phase_midplane,xir_midplane)
  use precision,only:p_
  use constants,only: one,two,three,pi,twopi,ii,zero
  !use path,only: xpath,ypath
  use radial_module,only: pfn,tfn,starting_surface_number,ending_surface_number !radial grid points
  use poloidal_harmonics,only: mh_low
  use toroidal_harmonics,only:nh
  use flux_grids,only: mpoloidal,x0
  use flux_grids,only: delta_q,r_new,z_new !r_new and z_new contains the magnetic surfaces found in subroutine calculate_contour()
  implicit none
  integer,intent(in):: k,mhtot,nflux
  complex(p_),intent(in)::ypath(2*mhtot+1,nflux),ypath2(2*mhtot,nflux),xir2d(mpoloidal,nflux)
  real(p_),intent(in)::frequency_khz,growth_khz, phase(mhtot,nflux),phase_total(mpoloidal,nflux),phase_midplane(nflux)
  real(p_):: theta(mpoloidal)
  complex(p_),intent(in):: xir_midplane(mhtot,nflux)
  real(p_),intent(in):: xpath(nflux)
  character(80):: filename1,filename2  !dynamic file name. 
  integer:: i,j


  !--record the eigen-functions in data files with dynamic filenames
  write(filename1,'(a,i3.3,a)') 'eigenpressure',k,'.txt'  !dynamic file name. Note that the editor descriptor iw.m , where w specifies the width, m specifies the minimum number of digits to output.
  write(filename2,'(a,i3.3,a)') 'eigendisplacement',k,'.txt'  !dynamic file name. 

  open(113,file=filename1) !record the eigenfunctions
  open(114,file=filename2) !record the eigenfunctions

  do j=starting_surface_number,ending_surface_number
     write(113,*) sqrt(pfn(j)),sqrt(tfn(j)),&
          & (real(ypath(i,j)),i=1,mhtot),&
          & (imag(ypath(i,j)),i=1,mhtot),&
          & (abs(ypath(i,j)), i=1,mhtot)

     write(114,*) sqrt(pfn(j)),sqrt(tfn(j)), &
          & (real(ypath(i,j)),i=mhtot+1,2*mhtot), &
          & (imag(ypath(i,j)),i=mhtot+1,2*mhtot), &
          & (abs(ypath(i,j)), i=mhtot+1,2*mhtot)
  enddo
  close(113)
  close(114)

  write(filename1,'(a,i3.3,a)') 'xis',k,'.txt' 
  write(filename2,'(a,i3.3,a)') 'divxi',k,'.txt'  

  open(113,file=filename1) !record the eigenfunctions
  open(114,file=filename2) !record the eigenfunctions
  do j=starting_surface_number,ending_surface_number
     write(113,*) sqrt(pfn(j)), sqrt(tfn(j)), &
          & (real(ypath2(i,j)),i=1,mhtot),&
          & (imag(ypath2(i,j)),i=1,mhtot),&
          & (abs(ypath2(i,j)), i=1,mhtot)

     write(114,*) sqrt(pfn(j)), sqrt(tfn(j)), &
          & (real(ypath2(i,j)),i=mhtot+1,2*mhtot), &
          & (imag(ypath2(i,j)),i=mhtot+1,2*mhtot), &
          & (abs(ypath2(i,j)), i=mhtot+1,2*mhtot)
  enddo
  close(113)
  close(114)

  write(filename1,'(a,i3.3,a)') 'frequency',k,'.txt' 
  open(113,file=filename1) !record the eigenfunctions
  do j=starting_surface_number,ending_surface_number
     write(113,*) sqrt(pfn(j)), sqrt(tfn(j)),frequency_khz,growth_khz
  enddo
  close(113)

  open(11,file='toroidal_angle_shift.txt')
  do i=1,mpoloidal
     do j=1,nflux
        write(11,*) r_new(i,j),z_new(i,j),delta_q(i,j)
        write(11,*) ! write a blank line to tell gnuplot that this is a 2d grid file
     enddo
  enddo
  close(11)
  call plot_poloidal(mpoloidal,nflux,delta_q,delta_q,'toroidal_angle_shif2.txt')


  do i=1,mpoloidal
     theta(i)=zero+twopi/(mpoloidal-1)*(i-1)
     if (theta(i).ge.pi) theta(i)=theta(i)-twopi
  enddo


  write(filename1,'(a,i3.3,a)') 'xir2d',k,'.txt'  !dynamic file name. 
  open(113,file=filename1) !record the radial phase
  do j=1,nflux
!  do j=101,101
     do i=1,mpoloidal
        !do j=starting_surface_number,ending_surface_number
        !  write(113,*) r_new(i,j),z_new(i,j),imag(xir2d(i,j))
        write(113,*) r_new(i,j),z_new(i,j),real(xir2d(i,j)),imag(xir2d(i,j)),abs(xir2d(i,j)),theta(i)
     enddo
     write(113,*) ! write a blank line to tell gnuplot that this is a 2d grid file
!     write(113,*) ! write two blank lines to indicate the data block
  enddo
  close(113)

  ! call plot_poloidal(mpoloidal,nflux,imag(xir2d),real(xir2d),'xir2d.txt')

  write(*,*) 'maximum of real_xir=',maxval(real(xir2d)),'minmum of real_xir=',minval(real(xir2d))

  write(filename1,'(a,i3.3,a)') 'phase',k,'.txt'  !dynamic file name. 
  open(113,file=filename1) !record the radial phase
  do j=starting_surface_number,ending_surface_number
     !   write(113,*) sqrt(pfn(j)), (phase(i,j),i=1,mhtot),phase_total(j)
     write(113,*) sqrt(pfn(j)),phase_midplane(j),x0(j) ,(phase_total(i,j),i=1,mpoloidal)
  enddo
  close(113)

  write(filename1,'(a,i3.3,a)') 'xir_midplane',k,'.txt' 
  open(113,file=filename1) !record the eigenfunctions
  do j=starting_surface_number,ending_surface_number
     write(113,*) sqrt(pfn(j)), & 
          & (real(xir_midplane(i,j)),i=1,mhtot), &
          & (imag(xir_midplane(i,j)),i=1,mhtot), &
          & (abs(xir_midplane(i,j)), i=1,mhtot)
  enddo
  close(113)


end subroutine record_eigenfunctions


subroutine calculate_eigenfunction(n2,v,xpath,ypath,ypath2,xir2d,phase,phase_total,phase_midplane,xir_midplane)
  use precision,only: p_
  use constants,only:one,pi,ii,twopi
  use flux_grids,only: nflux,mpoloidal,delta_q
  use radial_module, only: ra,starting_surface_number,ending_surface_number
  use poloidal_harmonics,only: mh_low,mhtot
  use toroidal_harmonics,only: nh
  implicit none

  INTEGER,intent(in):: n2
  complex(p_),intent(in):: v(n2)
  real(p_),intent(out):: xpath(nflux)
  complex(p_),intent(out):: ypath(2*mhtot+1,nflux),ypath2(2*mhtot,nflux),xir2d(mpoloidal,nflux)
  real(p_),intent(out):: phase(mhtot,nflux),phase_total(mpoloidal,nflux),phase_midplane(nflux)
  complex(p_),intent(out):: xir_midplane(mhtot,nflux)
  complex(p_):: sum

  real(p_):: x1,x2
  complex(p_):: y(2*mhtot+1) !2*mhtot+1 is the number of functions, which is also the number of ordinary difference equations
  integer:: i,j,k,mhk
  real(p_):: theta(mpoloidal)
  !x1=psival_new(starting_surface_number)
  !x2=psival_new(ending_surface_number)
x1=ra(starting_surface_number)
x2=ra(ending_surface_number)
  call load(x1,v,y) !set the values of functions at the starting point
  call odeint_yj2(y,2*mhtot+1,x1,x2,xpath,ypath) !advanc from the innermost surface to the outmost surface, and record the values of P1 and Xi on every surface

  do i=1,mpoloidal
     theta(i)=0.+(i-1)*twopi/(mpoloidal-1)
  enddo

  do j=starting_surface_number,ending_surface_number
     do i=1,2*mhtot !exclude the wsq
        if(i.le.mhtot) then !determine the poloidal mode numbers of the eigenfunctions
           mhk=mh_low+(i-1)
        else 
           mhk=mh_low+(i-mhtot-1)
        endif
        !ypath(i,j)=ypath(i,j)*exp(ii*mhk*pi) !measured at theta=pi, !this factor is -1 or 1 depending on whether mhk is old or even
        !We usually evaluate the poloidal harmonics at theta=0. In GTAW, theta=0 is at the high-field side of the midplane, while, in NOVA, the location theta=0 is at the low-field side of the midplane. Thus, to compare GTAW's results with NOVA's, the above factor should be included.
     enddo
  enddo

  call poloidal_displacement_and_divergence_of_displacement(ypath,ypath2) !ypath2 records Xis (the poloidal displacement) and div(Xi) (the divergence of Xi) 



  !calculate the 2d structure of the radial displacement on the fai=0 plane in (R,fai,Z) coordinate system:
  xir2d=0._p_
  do i=1,mpoloidal
     do j=starting_surface_number,ending_surface_number
        sum=0.
        do k=1,mhtot
           mhk=mh_low+k-1
           sum=sum+ypath(k+mhtot,j)*exp(ii*mhk*theta(i))*exp(ii*nh*delta_q(i,j))
        enddo
        xir2d(i,j)=sum
     enddo
  enddo

  call calculate_phase(ypath,phase,phase_total,phase_midplane, xir_midplane)!calculate the phase of the radial displacement across the minor radius
  

end subroutine calculate_eigenfunction
