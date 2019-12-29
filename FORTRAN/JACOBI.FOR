      SUBROUTINE jacobi(kk,iplvl,matl)
C
C ---------------------------------------------------------------------
C I                              JACOBI                               I
C I								      I
C I   FELEMA vs 1.0		  DJB		      14 / Mar / 1994 I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      REAL aux,ep,relold,drdb2,auxe(3)
      INTEGER i,j,matl
C
C ---------------------------------------------------------------------
C I      CALCULATE ELEMENT JACOBI MATRIX AND RESIDUAL VECTOR          I
C ---------------------------------------------------------------------
C I                    PRELIMINARY CALCULATIONS                       I
C ---------------------------------------------------------------------
C
      relold=erel(kk)
      IF (matl.EQ.2) THEN
	 drdb2=const(1)*const(2)*EXP(const(2)*EBSQ(KK))
      ELSE
	 drdb2=0.0
      END IF
C
C ---------------------------------------------------------------------
C I                    ELEMENT STIFFNESS MATRIX                       I
C ---------------------------------------------------------------------
C
      DO 802 i=1,3
	 DO 800 j=1,3
	    ESM(i,j)=relold*esmuni(kk,i,j)
800	 CONTINUE
802   CONTINUE
C
C ---------------------------------------------------------------------
C I                 RESIDUAL AND AUXILLIARY VECTORS                   I
C ---------------------------------------------------------------------
C
      DO 806 i=1,3
	 ep=0.0
	 auxe(i)=0.0
	 DO 804 j=1,3
	    ep=ep+esm(i,j)*phi(J)
	    auxe(i)=auxe(i)+esmuni(kk,i,j)*phi(j)
804	 CONTINUE
	 resv(i)=ep-ef(kk,i)
806   CONTINUE
C
C ---------------------------------------------------------------------
C I                  THE ELEMENT JACOBIAN MATRIX                      I
C ---------------------------------------------------------------------
C
      DO 810 i=1,3
	 DO 808 j=1,3
	    aux=drdb2*auxe(i)*auxe(j)
	    ajac(i,j)=esm(i,j)+(4/ar2)*aux
808	 CONTINUE
810   CONTINUE
C
C ---------------------------------------------------------------------
C I      OUTPUT ELEMENT JACOBIAN AND RESIDUAL VECTOR IF IPLVL=1       I
C ---------------------------------------------------------------------
C
      IF (iplvl.EQ.0) RETURN
      WRITE (noutp,501) kk
501   FORMAT (/10X,'Element',I4/10X,'Residual Vector',10X,
     + 'Jacobi Matrix')
      DO 812 i=1,3
	 WRITE (noutp,502) resv(i), (ajac(i,j), j=1,3)
812   CONTINUE
502   FORMAT (10X,E12.5,10X,4E13.5)
      RETURN
      END


