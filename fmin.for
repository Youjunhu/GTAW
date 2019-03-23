      FUNCTION fmin(x)
      use precision,only:p_
      INTEGER n,NP
cyj      REAL(p_) fmin,x(*),fvec
        REAL(p_) fmin
        complex(p_):: x(*),fvec
cyj      PARAMETER (NP=40)
      PARAMETER (NP=80)
      COMMON /newtv/ fvec(NP),n
      SAVE /newtv/
CU    USES funcv
      INTEGER i
      REAL(p_) sum
      call funcv(n,x,fvec)
      sum=0.
      do 11 i=1,n
        sum=sum+abs(fvec(i))**2
11    continue
      fmin=0.5*sum
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
