      SUBROUTINE rsolv(a,n,np,d,b)
         use precision,only:p_
      implicit none
      INTEGER n,np
      REAL(p_) a(np,np),b(n),d(n)
      INTEGER i,j
      REAL(p_) sum
      b(n)=b(n)/d(n)
      do 12 i=n-1,1,-1
        sum=0.
        do 11 j=i+1,n
          sum=sum+a(i,j)*b(j)
11      continue
        b(i)=(b(i)-sum)/d(i)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
