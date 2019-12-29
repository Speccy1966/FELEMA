      SUBROUTINE alstmf(kk,iplvl,iprint)
C
C ---------------------------------------------------------------------
C I				 ALSTMF				      I
C I								      I
C I   Subroutine added to include an axisymmetric linear triangular   I
C I   element for magnetic and electric analyses.		      I
C I                                                                   I
C I   FELEMA vs 1.0		  DJB		      7 / July / 1994 I
C ---------------------------------------------------------------------
C
      PARAMETER (PI=3.141592654)
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      INTEGER i,j
      REAL ba(3),ca(3),ar2,eg(3,3),rbar,fcons,econs
      DATA eg/2.,1.,1.,1.,2.,1.,1.,1.,2./
C
C ---------------------------------------------------------------------
C I		 AXIYSMMETRIC LINEAR TRIANGULAR ELEMENT		      I
C ---------------------------------------------------------------------
C
      rbar=(x(1)+x(2)+x(3))/3
      CALL triangle(x,y,ba,ca,ar2)
      IF (ABS(ar2).LT.1.0E-8) CALL warning(iprint,kk,noutp)
      fcons=PI*qe*ar2/12
      econs=PI*rbar/ar2
      DO 783 i=1,3
	 ef(kk,i)=0.0
	 DO 782 j=1,3
	    ef(kk,i)=ef(kk,i)+fcons*eg(i,j)*x(i)
	    esmuni(kk,i,j)=(ba(i)*ba(j)+ca(i)*ca(j))*econs
	    esm(i,j)=dxe*esmuni(kk,i,j)
782	 CONTINUE
783   CONTINUE
      RETURN
      END







