      SUBROUTINE warning(iprint,kk,noutp)
C
C ---------------------------------------------------------------------
C I                              WARNING                              I
C I								      I
C I   FELEMA vs 1.0		  DJB		      7 / July / 1994 I
C ---------------------------------------------------------------------
C I                         DIAGNOSTIC OUTPUT                         I
C ---------------------------------------------------------------------
C
      IF (iprint.EQ.1) WRITE (noutp,492) kk
492   FORMAT (//10X,'*** WARNING ***'//
     + 10X,'THE AREA OF ELEMENT ',I4,
     + ' IS LESS THAN 1E-8'/
     + 10X,'OR THE NODE NUMBERS ARE IN THE WRONG ORDER'/
     + 10X,'OR THE NODES FORM A STRAIGHT LINE')
      RETURN
      END
