      SUBROUTINE output(title,ne,ityp)
C
C ---------------------------------------------------------------------
C I                              OUTPUT                               I
C I                                                                   I
C I   Subroutine to print out results to output file                  I
C I                                                                   I
C I   FELEMA vs 1.0               DJB                 17 / Mar / 1994 I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      CHARACTER title(20)
      INTEGER iline, kk
      WRITE (noutp,220) title
220   FORMAT (//10X,20A)
C
C ---------------------------------------------------------------------
C I                  OUTPUT OF THE POTENTIAL HEADINGS                 I
C ---------------------------------------------------------------------
C
      WRITE (noutp,230)
      WRITE (noutp,240) (i,soln(i), i=1,np)
230   FORMAT (//10X,'Calculated quantities'/
     +12X,'Nodal values for potential ')
240   FORMAT (12X,I3,E14.5,3X,I3,E14.5,3X,I3,E14.5)
C
C-----Write POTENTIAL files
C
      WRITE (13,320) (i,soln(i),i=1,np)
      WRITE (14,320) (i,cur(i),i=1,np)
320   FORMAT (I5,E15.5)
C
C ---------------------------------------------------------------------
C I                  OUTPUT OF THE GRADIENT HEADINGS                  I
C ---------------------------------------------------------------------
C
      iline=0
      DO 846 kk=1,ne
C
C ---------------------------------------------------------------------
C I       CONVERT GRADIENTS TO A/cm^2 FOR ITYP 2 & 3                  I
C ---------------------------------------------------------------------
C
	 IF ((ityp.EQ.3).OR.(ityp.EQ.4)) THEN         
	    gradx(kk)=gradx(kk)/10000
	    grady(kk)=grady(kk)/10000
	    vsum(kk)=vsum(kk)/10000
	 END IF
	 IF ((iline.EQ.0).AND.(ityp.EQ.1)) THEN 
	    WRITE (noutp,250) title
	    WRITE (noutp,300)
	 END IF
	 IF ((iline.EQ.0).AND.(ityp.GT.1)) THEN 
	    WRITE (noutp,250) title
	    WRITE (noutp,302)
	 END IF
C
C ---------------------------------------------------------------------
C I                      INCREASE THE LINE COUNT                      I
C ---------------------------------------------------------------------
C
	 iline=iline+2
	 IF (iline.GT.50) iline=0
	 WRITE (noutp,311) kk,gradx(kk),grady(kk),vsum(kk)
C
C-----Write GRADIENT files
C
	 WRITE (15,322) kk,gradx(kk),grady(kk),vsum(kk)
	 WRITE (16,322) kk,gradx(kk),grady(kk),vsum(kk)
846   CONTINUE
250   FORMAT (//10X,20A)
300   FORMAT (//10X,'Magnetic Flux Density ''B'' Field.'/
     +10X,'Element',4X,'Location',10X,'B(X)',11X,'B(Y)',
     +10X,'Mod B')
 302  FORMAT (//10X,'Electric Current Density ''J'' Field.'/
     +10X,'Element',4X,'Location',10X,'J(X)',11X,'J(Y)',
     +10X,'Mod J')
311   FORMAT (13X,I3,5X,'CENTRE',2X,3E15.5)
322   FORMAT (I5,3E15.5)
      RETURN
      END










