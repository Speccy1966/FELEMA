      SUBROUTINE calccurl(kk)
C
C ---------------------------------------------------------------------
C I                           CALCCURL                                I
C I								      I
C I   FELEMA vs 1.0		  DJB		      7 / July / 1994 I
C I                                                                   I
C I   This version of CALCCURL for 2D cartesian co-ordinates only     I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      REAL xx(3),yy(3),bb(3),cc(3),vsold,dbwt
C
C ---------------------------------------------------------------------
C I       RETRIEVAL OF THE NODAL COORDINATES, THE NODE NUMBERS        I
C I                   AND THE NODAL VALUES OF A                       I
C ---------------------------------------------------------------------
C
      vsold=sqrt(ebsq(kk))
      DO 720 i=1,3
	 ns(i)=nel(kk,i)
	 j=ns(i)
	 xx(i)=xc(j)
	 yy(i)=yc(j)
	 phi(i)=soln(j)
720   CONTINUE
C
C ---------------------------------------------------------------------
C I             EVALUATION OF THE ELEMENT CURL                        I
C ---------------------------------------------------------------------
C
      CALL triangle(xx,yy,bb,cc,ar2)
C
C-----The i component of Curl is the Y potential gradient
C 
      gradx(kk)=(cc(1)*phi(1)+cc(2)*phi(2)+cc(3)*phi(3))/ar2
C
C-----The j component of Curl is the negative X potential gradient
C
      grady(kk)=-(bb(1)*phi(1)+bb(2)*phi(2)+bb(3)*phi(3))/ar2
C
C-----The magnitude is the vector sum
C
      ebsq(kk)=gradx(kk)**2+grady(kk)**2
      vsum(kk)=sqrt(ebsq(kk))
      IF (vsum(kk).GT.0.0) THEN
	 dbwt=ABS((vsum(kk)-vsold)/vsum(kk))
	 ELSE
	 dbwt=0.0
      END IF
      sumwtd=sumwtd+dbwt
      RETURN
      END




