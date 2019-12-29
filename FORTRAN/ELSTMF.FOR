      SUBROUTINE elstmf(kk,iplvl,iprint)
C
C ---------------------------------------------------------------------
C I                              ELSTMF                               I
C I								      I
C I   FELEMA vs 1.0		  DJB		      7 / July / 1994 I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      INTEGER i,j
      REAL ba(3),ca(3),ar2
C
C ---------------------------------------------------------------------
C I                    LINEAR TRIANGULAR ELEMENT                      I
C ---------------------------------------------------------------------
C
      CALL triangle(x,y,ba,ca,ar2)
      IF (ABS(ar2).LT.1.0E-8) CALL warning(iprint,kk,noutp)
      DO 781 i=1,3
	 ef(kk,i)=qe*ar2/6
	 DO 780 j=1,3
	    esmuni(kk,i,j)=(ba(i)*ba(j)+ca(i)*ca(j))/(ar2*2.)
	    esm(i,j)=dxe*esmuni(kk,i,j)
780	 CONTINUE
781   CONTINUE
C
C ---------------------------------------------------------------------
C I          OUTPUT OF THE ELEMENT MATRICES IF IPLVL=1                I
C ---------------------------------------------------------------------
C
      IF (iplvl.EQ.0) RETURN
      WRITE(noutp,490) kk
490   FORMAT (/10X,'Element',I4/10X,'Force Vector',10X,
     + 'Stiffness Matrix')
      DO 784 i=1,3
	 WRITE (noutp,491) ef(kk,i), (esm(i,j), j=1,3)
784   CONTINUE
491   FORMAT (10X,E12.5,10X,4E13.5)
      RETURN
      END



