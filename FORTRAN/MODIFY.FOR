      SUBROUTINE modify
C
C ---------------------------------------------------------------------
C I                             MODIFY                                I
C I								      I
C I   FELEMA vs 1.0		  DJB		      14 / Mar / 1994 I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      INTEGER j,k,m,ij,kj,kk
C
C ---------------------------------------------------------------------
C I          MODIFICATION OF THE GLOBAL STIFFNESS MATRIX              I
C I          AND THE GLOBAL FORCE VECTOR USING THE METHOD             I
C I               OF DELETION OF ROWS AND COLUMNS                     I
C ---------------------------------------------------------------------
C
      DO 786 kk=1,nbc
	 k=ib(kk)-1
	 DO 785 j=2,nbw
	    m=ib(kk)+j-1
	    IF (m.LE.np) THEN
	       ij=jgsm+(j-1)*np+ib(kk)-(j-1)*(j-2)/2
	       a(jgf+m)=a(jgf+m)-a(ij)*bv(kk)
	       a(ij)=0.0
            END IF
   	    IF (k.GT.0) THEN
	       kj=jgsm+(j-1)*np+k-(j-1)*(j-2)/2
	       a(jgf+k)=a(jgf+k)-a(kj)*bv(kk)
	       a(kj)=0.0
	       k=k-1
            END IF
785	 CONTINUE
	 a(jgf+ib(kk))=a(jgsm+ib(kk))*bv(kk)
786   CONTINUE
      RETURN
      END
