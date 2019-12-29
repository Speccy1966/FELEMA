      SUBROUTINE magnet(nlchk)
C
C ---------------------------------------------------------------------
C I				 MAGNET				      I
C I								      I
C I   Subroutine to set up the magnetic properties of material 2.     I
C I								      I
C I   FELEMA vs 1.0		  DJB		      17 / Mar / 1994 I
C ---------------------------------------------------------------------
C
C ---------------------------------------------------------------------
C I             OUTPUT OF THE NON LINEAR CHARACTERISTIC               I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      WRITE (noutp,85)
85    FORMAT (//10X,'Non Linear Material Characteristic'
     + //10X,'H (A/m)',7X,'B (Tesla)')
      WRITE (noutp,86) (h(i),b(i),i=1,maxbh)
86    FORMAT (10X,2F10.4)
C
C ---------------------------------------------------------------------
C I               CONVERT B-H DATA TO REL-BSQ DATA                    I
C ---------------------------------------------------------------------
C
      DO 775 i=1,maxbh
	 rel(i)=h(i)/b(i)
	 bsq(i)=b(i)**2
775   CONTINUE
      WRITE (noutp,434)
434   FORMAT (/10X,'BSq(T^2)',6X,'Rel(m/H)')
      WRITE (noutp,435) (bsq(i),rel(i),i=1,maxbh)
435   FORMAT (10X,F8.6,2X,F10.4)
      WRITE (17,205) (bsq(i),rel(i),i=1,maxbh)
205   FORMAT (2E15.5)
      IF (nlchk.NE.1) CALL convbh
      WRITE (noutp,485)
485   FORMAT (/10X,'Constants that give the best fit curve for'
     + //10X, 'Rel = C1 * EXP(C2 * BSq) + C3')
      WRITE (noutp,207) const(1)
      WRITE (noutp,208) const(2)
      WRITE (noutp,209) const(3)
207   FORMAT (/10X,'C1 :',F10.4)
208   FORMAT (10X,'C2 :',F10.4)
209   FORMAT (10X,'C3 :',F10.4)
      WRITE (17,206) (const(i),i=1,3)
206   FORMAT (E15.5)
C
C ---------------------------------------------------------------------
C I       CALCULATE THE MAXIMUM GRADIENT OF THE REL-BSQ CURVE         I
C I                   AT REL(maxbh), BSQ(maxbh)                       I
C ---------------------------------------------------------------------
C
      grdmax=const(1)*const(2)*EXP(const(2)*bsq(maxbh))
      WRITE (noutp,486) grdmax
486   FORMAT (/10X,'Maximum gradient of Rel - BSq curve '
     + /10X,'grdmax = ',F10.4)
      RETURN
      END











