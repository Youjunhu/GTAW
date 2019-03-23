      FUNCTION rtbis(func,x1,x2,xacc,xmaxis,zmaxis,slope,psival)
      !modified by Youjun Hu
      use precision,only: p_
      implicit none
      INTEGER JMAX
      REAL(p_) rtbis,x1,x2,xacc,func,xmaxis,zmaxis,slope,psival
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      REAL(p_) dx,f,fmid,xmid
      fmid=func(xmaxis,zmaxis,slope,psival,x2)
      f=   func(xmaxis,zmaxis,slope,psival,x1)
!      write(*,*) 'f1=', f, 'f2=',fmid
      if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        fmid=func(xmaxis,zmaxis,slope,psival,xmid)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      stop 'too many bisections in rtbis'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


      FUNCTION rtbis2(func,x1,x2,xacc,nx,psival_nx,qpsi_nx,tmp_y2)
      use precision,only:p_
      implicit none
      INTEGER JMAX
      REAL(p_) rtbis2,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=100)
      INTEGER j
      REAL(p_) dx,f,fmid,xmid
      integer,intent(in):: nx
      real(p_),intent(in):: psival_nx(nx),qpsi_nx(nx), tmp_y2(nx)
      fmid=func(x2,nx,psival_nx,qpsi_nx,tmp_y2)
      f=func(x1,nx,psival_nx,qpsi_nx,tmp_y2)
      if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis2'
      if(f.lt.0.)then
        rtbis2=x1
        dx=x2-x1
      else
        rtbis2=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis2+dx
        fmid=func(xmid,nx,psival_nx,qpsi_nx,tmp_y2)
        if(fmid.le.0.)rtbis2=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      stop 'too many bisections in rtbis2'
      END


      FUNCTION rtbis0(func,x1,x2,xacc)
      use precision,only:p_
      implicit none
      INTEGER JMAX
      REAL(p_) rtbis0,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=100)
      INTEGER j
      REAL(p_) dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis0'
      if(f.lt.0.)then
        rtbis0=x1
        dx=x2-x1
      else
        rtbis0=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis0+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbis0=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      stop 'too many bisections in rtbis0'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
