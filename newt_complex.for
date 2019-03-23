      SUBROUTINE newt(x,n,check)
      use precision,only: p_
      implicit none
      INTEGER n,nn,NP,MAXITS
      LOGICAL check
cyj      REAL(p_) x(n),fvec,
      complex(p_) x(n),fvec
      real(p_)::TOLF,TOLMIN,TOLX,STPMX
cyj      PARAMETER (NP=40,MAXITS=200,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7,
      PARAMETER (NP=80,MAXITS=300,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-6,
     *     STPMX=100.)
      COMMON /newtv/ fvec(NP),nn
      SAVE /newtv/
C     U    USES fdjac,fmin,lnsrch,lubksb,ludcmp
      INTEGER i,its,j,indx(NP)
cyj      REAL(p_) d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),
      complex(p_) sum,fjac(NP,NP),g(NP),p(NP),xold(NP)
      real(p_)::den,f,fold,fmin, stpmax,temp,test
      EXTERNAL fmin
      INTEGER ::info,ipiv(n)
      nn=n
      f=fmin(x)
      test=0.
      do 11 i=1,n
         if(abs(fvec(i)).gt.test)test=abs(fvec(i))
 11   continue
      if(test.lt..01*TOLF)  then
         write(*,*) 'return from location 0'         
         return
      endif
      sum=0.
      do 12 i=1,n
         sum=sum+x(i)**2
 12   continue
      stpmax=STPMX*max(sqrt(abs(sum)),float(n))
      do 21 its=1,MAXITS
cyj         call fdjac(n,x,fvec,NP,fjac)
         call fdjac(n,x,fvec,fjac)
         do 14 i=1,n
            sum=0.
            do 13 j=1,n
               sum=sum+fjac(j,i)*fvec(j)
 13         continue
            g(i)=sum
 14      continue
         do 15 i=1,n
            xold(i)=x(i)
 15      continue
         fold=f
         do 16 i=1,n
            p(i)=-fvec(i)
 16      continue
cyj         call ludcmp(fjac,n,NP,indx,d)
cyj         call lubksb(fjac,n,NP,indx,p)
         call zgesv(n,1,fjac,n,ipiv,p,n,info) !lapack routine to slove linear equations system
         call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
cyj
         write(*,*) 'the new root is ',x(n-1),x(n) !, x(n-3),x(n-2)
cyj         
         test=0.
         do 17 i=1,n
            if(abs(fvec(i)).gt.test)test=abs(fvec(i))
 17      continue
         if(test.lt.TOLF)then
            check=.false.
            write(*,*) 'return from location 1'
            return
         endif
         if(check)then
            test=0.
            den=max(f,.5*n)
            do 18 i=1,n
               temp=abs(g(i))*max(abs(x(i)),1.)/den
               if(temp.gt.test)test=temp
 18         continue
            if(test.lt.TOLMIN)then
               check=.true.
            else
               check=.false.
            endif
            write(*,*) 'return from location 2'
            return
         endif
         test=0.
         do 19 i=1,n
            temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
            if(temp.gt.test)test=temp
 19      continue
         if(test.lt.TOLX) then
            write(*,*) 'return from location 3'
            return
         endif
 21   continue
c     pause 'MAXITS exceeded in newt'
      stop 'MAXITS exceeded in newt'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
