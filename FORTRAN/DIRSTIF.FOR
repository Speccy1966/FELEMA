      SUBROUTINE dirstif(ns,fv,asm)
C
C ---------------------------------------------------------------------
C I				 DIRSTIF                              I
C I								      I
C I   FELEMA vs 1.0		  DJB		      14 / Mar / 1994 I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      COMMON /av/ a(maxav), jgf, jgsm, np, nbw
      DIMENSION ns(3), fv(3), asm(3,3)
      INTEGER i,j,ii,jj,j1
C
C ---------------------------------------------------------------------
C I                   DIRECT STIFFNESS PROCEDURE                      I
C ---------------------------------------------------------------------
C
      DO 713 i=1,3
	 ii=ns(i)
	 a(jgf+ii)=a(jgf+ii)+fv(i)
	 DO 712 j=1,3
	    jj=ns(j)+1-ii
	    IF (jj.GT.0) THEN
	       j1=jgsm+(jj-1)*np+ii-(jj-1)*(jj-2)/2
	       a(j1)=a(j1)+asm(i,j)
            END IF
712	 CONTINUE
713   CONTINUE
      RETURN
      END
