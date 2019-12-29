      SUBROUTINE convbh
C
C ---------------------------------------------------------------------
C I                                                                   I
C I                            SUBROUTINES                            I
C I                                                                   I
C ---------------------------------------------------------------------
C I				 CONVBH                               I
C I								      I
C I   FELEMA vs 1.0		  DJB		      14 / Mar / 1994 I
C ---------------------------------------------------------------------
C I   CONVERTS DATA FROM A STANDARD B-H CURVE INTO AN EQUATION FOR    I
C I   RELUCTIVITY VS FLUX DENSITY SQUARED. THE FORM OF THE EQUATION   I
C I                                                                   I
C I		    REL = C1 * EXP(C2 * BSQ) + C3                     I
C I                                                                   I
C I     WAS PROPOSED BY BRAUER (1975). THE NEWTON ALGORITHM FOR       I
C I   FINDING THE COEFFICIENTS C1, C2 AND C3 IS DERIVED FROM HOOLE    I
C I                          AND HOOLE (1986).                        I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      DIMENSION relfit(maxbh), c2(2)
      REAL c1,c3,midrel,midbsq,err,toterr,minerr,test,aa,bb,cc,dd,ee,
     + f,dfdc2
      INTEGER i,k,lend
      minerr=10000
      WRITE (noutp,440)
440   FORMAT (/36X,'CONVBH',/35X,'--------'//
     + 10X,'Newton method for calculating Nu-BSq curve constants')
C
C ---------------------------------------------------------------------
C I             THE FIRST AND LAST DATA POINTS ARE FIXED,             I
C I        THE ROUTINE LOOPS ABOUT THE REMAINING DATA POINTS          I
C I        TO FIND THE MID POINT WHICH GIVES THE LEAST ERROR          I
C ---------------------------------------------------------------------
C
      lend=maxbh-1
      DO 778 k=2,lend
	 midrel=rel(k)
	 midbsq=bsq(k)
	 toterr=0.0
	 err=0.0
C
C ---------------------------------------------------------------------
C I         PERFORM 20 NEWTON ITERATIONS TO CONVERGE ON C2            I
C ---------------------------------------------------------------------
C
	 j=1
	 test=1.0
	 c2(1)=1.0
	 WRITE (noutp,450)
450	 FORMAT (/10X,'I',5X,'C2(I)',3X,'F(C2)',6X,'F''C(2)',8X,
     +    'C2(I+1)')
	 DO 776 i=1,20
C
C ---------------------------------------------------------------------
C I                 CALCULATE THE FUNCTION OF C2                      I
C ---------------------------------------------------------------------
C
	    aa=EXP(c2(1)*midbsq)
	    bb=EXP(c2(1)*bsq(maxbh))
	    f=(aa-1)/(bb-1)-(midrel-rel(1))/(rel(maxbh)-rel(1))
C
C ---------------------------------------------------------------------
C I          CALCULATE THE DIFFERENTIAL OF THE FUNCTION               I
C ---------------------------------------------------------------------
C
	    cc=(bb-1)*midbsq*aa
	    dd=(aa-1)*bsq(maxbh)*bb
	    ee=(bb-1)**2
	    dfdc2=(cc-dd)/ee
C
C ---------------------------------------------------------------------
C I                 CALCULATE THE NEW VALUE OF C2                     I
C ---------------------------------------------------------------------
C
	    c2(2)=c2(1)-f/dfdc2
	    WRITE (noutp,460) i,c2(1),f,dfdc2,c2(2)
460	    FORMAT (8X,I3,F10.4,2E12.4,F10.4)
C
C ---------------------------------------------------------------------
C I                     CHECK FOR CONVERGENCE                         I
C ---------------------------------------------------------------------
C
	    test=ABS(c2(2)-c2(1))
	    IF (test.LT.0.0001) GOTO 625
	    c2(1)=c2(2)
776	 CONTINUE
C
C ---------------------------------------------------------------------
C I              IF NOT CONVERGED THEN HALT PROGRAM                   I
C ---------------------------------------------------------------------
C
	 WRITE (noutp,465)
465	 FORMAT (/10X,'Incorrect B-H data'/
     +	  /10X,'** INPUT TERMINATED **')
	 STOP
C
C ---------------------------------------------------------------------
C I                           SUCCESS !                               I
C ---------------------------------------------------------------------
C
625	 WRITE (noutp,470) i
470	 FORMAT (/10X,'Converged after ',I2,' iterations.')
C
C ---------------------------------------------------------------------
C I            CALCULATE REMAINING CONSTANTS C1 AND C3                I
C ---------------------------------------------------------------------
C
	 c1=(midrel-rel(1))/(aa-1)
	 c3=rel(1)-c1
	 WRITE (noutp,475) c1
	 WRITE (noutp,476) c2(1)
	 WRITE (noutp,477) c3
475      FORMAT (/10X,'C1 :',F10.4)
476	 FORMAT (10X,'C2 :',F10.4)
477	 FORMAT (10X,'C3 :',F10.4)
C
C ---------------------------------------------------------------------
C I   CALCULATE THE FITTED CURVE AND CHECK THE ERROR BETWEEN KNOWN    I
C I   AND CALCULATED VALUES. THE MIDPOINT WHICH GIVES THE SMALLEST    I
C I      ERROR FIXES THE THREE EQUATION CONSTANTS C1, C2 AND C3       I
C ---------------------------------------------------------------------
C
	 DO 777 i=1,maxbh
	    relfit(i)=c1*EXP(c2(1)*bsq(i))+c3
	    err=ABS((rel(i)-relfit(i))/rel(i))
	    toterr=toterr+err
777	 CONTINUE
	 IF (toterr.LT.minerr) THEN
	    const(1)=c1
	    const(2)=c2(1)
	    const(3)=c3
	    minerr=toterr
	 END IF
	 WRITE (noutp,480) k,toterr
480      FORMAT (/10X,'Data point : ',I3,'  Error : ',F10.4)
778   CONTINUE
      RETURN
      END



