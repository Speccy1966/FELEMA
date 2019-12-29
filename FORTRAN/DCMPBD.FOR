      SUBROUTINE dcmpbd
C
C ---------------------------------------------------------------------
C I                                                                   I
C I            MATRIX SOLVER FOR SYSTEM OF LINEAR EQUATIONS,          I
C I                       GAUSSIAN ELIMINATION.                       I
C I                                                                   I
C I								      I
C I   FELEMA vs 1.0		  DJB		      14 / Mar / 1994 I
C ---------------------------------------------------------------------
C I                             DCMPBD                                I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      COMMON /av/ a(maxav), jgf, jgsm, np, nbw
      INTEGER i,j,k,ii,jk,mj,mk,nd,nl,nk,inl,ink,np1
C
C ---------------------------------------------------------------------
C I        DECOMPOSITION OF A BANDED MATRIX INTO AN UPPER             I
C I          TRIANGULAR FORM USING GAUSSIAN ELIMINATION               I
C ---------------------------------------------------------------------
C
      np1=np-1
      DO 789 i=1,np1
	 mj=i+nbw-1
	 IF (mj.GT.np) mj=np
	 nj=i+1
	 mk=nbw
	 IF((np-i+1).LT.nbw) mk=np-i+1
	 nd=0
	 DO 788 j=nj,mj
	    mk=mk-1
	    nd=nd+1
	    nl=nd+1
	    DO 787 k=1,mk
	       nk=nd+k
	       jk=jgsm+(k-1)*np+j-(k-1)*(k-2)/2
	       inl=jgsm+(nl-1)*np+i-(nl-1)*(nl-2)/2
	       ink=jgsm+(nk-1)*np+i-(nk-1)*(nk-2)/2
	       ii=jgsm+i
	       a(jk)=a(jk)-a(inl)*a(ink)/a(ii)
787	    CONTINUE
788	 CONTINUE
789   CONTINUE
      RETURN
      END
