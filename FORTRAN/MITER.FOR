      SUBROUTINE miter(ne,iter,niter,iterc,iplvl,rhoc,rhonr,iflag)
C
C ---------------------------------------------------------------------
C I                            MITER                                  I
C I                                                                   I
C I PROCESS SOLUTION FOR THE ELEMENT GRADIENTS (B VALUES), USE THESE  I
C I   TO DETERMINE THE MATERIAL PROPERTIES FOR THE NEXT ITERATION.    I
C I CALCULATE THE ERROR BETWEEN ITERATIONS AND CHECK FOR CONVERGENCE. I
C I                                                                   I
C I   FELEMA vs 1.0		  DJB		      7 / July / 1994 I
C ---------------------------------------------------------------------
C I          CALCULATION OF THE ELEMENT GRADIENTS, B=CURL A           I
C ---------------------------------------------------------------------
C I                     LOOP AROUND KK ELEMENTS                       I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      REAL eb1
      DO 722 kk=1,ne
C
C ---------------------------------------------------------------------
C I                CALL SUBROUTINE FOR B VALUES                       I
C ---------------------------------------------------------------------
C
         CALL calccurl(kk)
C
C ---------------------------------------------------------------------
C I         AFTER PRELIMINARY ITERATIONS USE NEWTON RAPHSON           I
C ---------------------------------------------------------------------
C
     	 IF ((iter.LE.iterc).AND.(nmtl(kk).EQ.2))
     +	 CALL calcnu(kk,rhoc)
	 IF (iter.GT.iterc) THEN
	    CALL jacobi(kk,iplvl,nmtl(kk))
	    CALL dirstif(ns,resv,ajac)
	 END IF
722   CONTINUE
C
C ---------------------------------------------------------------------
C I      MODIFICATION AND SOLUTION OF THE JACOBIAN SYSTEM OF          I
C I      EQUATIONS AND OUTPUT OF THE CALCULATED ERROR VALUES          I
C ---------------------------------------------------------------------
C
      IF (iter.GT.iterc) THEN
	 CALL modify
	 CALL dcmpbd
	 CALL slvbd
C
C ---------------------------------------------------------------------
C I          ADD THE ERROR TERMS TO THE ORIGINAL SOLUTION             I
C ---------------------------------------------------------------------
C
	 DO 842 i=1,np
	    soln(i)=soln(i)-a(i)
842	 CONTINUE
C
C ---------------------------------------------------------------------
C I      CALCULATE NEW ELEMENT GRADIENTS AND RELUCTIVITY VALUES       I
C I                    FOR STEEL ELEMENTS ONLY.                       I
C ---------------------------------------------------------------------
C
	 DO 844 kk=1,ne
	    IF (nmtl(kk).EQ.2) THEN
	       CALL calccurl(kk)
	       CALL calcnu(kk,rhonr)
	    END IF
844	 CONTINUE
      END IF
C
C ---------------------------------------------------------------------
C I      WRITE OUT THE ITERATION NUMBER AND CHECK VALUE SUMWTD        I
C ---------------------------------------------------------------------
C
      eb1=sqrt(ebsq(1))
      WRITE (nscr,420) iter,sumwtd,eb1
      WRITE (noutp,420) iter,sumwtd,eb1
420   FORMAT (10X,I6,E15.5,3X,E12.4)
C
C ---------------------------------------------------------------------
C I                     CHECK FOR CONVERGENCE                         I
C ---------------------------------------------------------------------
C
      IF (iter.EQ.niter) THEN
          WRITE (nscr,324) iter
          WRITE (noutp,324) iter
      END IF    
      IF (sumwtd.LT.1e-1) THEN
         WRITE (nscr,430) iter
         WRITE (noutp,430) iter
         iflag=0
      END IF
430   FORMAT (/10X,'PROCESS HAS CONVERGED AFTER ',I4,' ITERATIONS')
324   FORMAT (/10X,'PROCESS HAS NOT CONVERGED AFTER ',I4,' ITERATIONS'
     + //10X,'Modify input file and run again.')
      RETURN
      END




