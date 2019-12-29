      SUBROUTINE slvbd
C
C ---------------------------------------------------------------------
C I                             SLVBD                                 I
C I								      I
C I   FELEMA vs 1.0		  DJB		      14 / Mar / 1994 I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      COMMON /av/ a(maxav), jgf, jgsm, np, nbw
      INTEGER i,j,k,l,n,ij,il,mj,nj,np1
      REAL sum
C
C ---------------------------------------------------------------------
C I            DECOMPOSITION OF THE GLOBAL FORCE VECTOR               I
C ---------------------------------------------------------------------
C
      np1=np-1
      DO 791 i=1,np1
	 mj=i+nbw-1
	 IF (mj.GT.np) mj=np
	 nj=i+1
	 l=1
	 DO 790 j=nj,mj
	    l=l+1
	    il=jgsm+(l-1)*np+i-(l-1)*(l-2)/2
	    a(jgf+j)=a(jgf+j)-a(il)*a(jgf+i)/a(jgsm+i)
790	 CONTINUE
791   CONTINUE
C
C ---------------------------------------------------------------------
C I       BACKWARD SUBSTITUTION TO CALCULATE THE NODAL VALUES         I
C ---------------------------------------------------------------------
C
      a(np)=a(jgf+np)/a(jgsm+np)
      DO 793 k=1,np1
	 i=np-k
	 mj=nbw
	 IF ((i+nbw-1).gt.np) mj=np-i+1
	 sum=0.0
	 DO 792 j=2,mj
	    n=i+j-1
	    ij=jgsm+(j-1)*np+i-(j-1)*(j-2)/2
	    sum=sum+a(ij)*a(n)
792	 CONTINUE
	 a(i)=(a(jgf+i)-sum)/a(jgsm+i)
793   CONTINUE
      RETURN
      END
