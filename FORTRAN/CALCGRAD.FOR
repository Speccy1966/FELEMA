      SUBROUTINE calcgrad(kk)
C
C ---------------------------------------------------------------------
C I                           CALCGRAD                                I
C I								      I
C I   FELEMA vs 1.0		  DJB		      7 / July / 1994 I
C I								      I
C I   CALCGRAD work for both 2D cartesian and 2D cylindrical          I
C I   co-ordinate systems			                      I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      REAL xx(3),yy(3),bb(3),cc(3)
C
C ---------------------------------------------------------------------
C I       RETRIEVAL OF THE NODAL COORDINATES, THE NODE NUMBERS        I
C I                   AND THE NODAL VALUES OF A                       I
C ---------------------------------------------------------------------
C
      DO 720 i=1,3
	 ns(i)=nel(kk,i)
	 j=ns(i)
	 xx(i)=xc(j)
	 yy(i)=yc(j)
	 phi(i)=soln(j)
720   CONTINUE
C
C ---------------------------------------------------------------------
C I             EVALUATION OF THE ELEMENT GRADIENTS                   I
C ---------------------------------------------------------------------
C
      CALL triangle(xx,yy,bb,cc,ar2)
C
C-----The i component of Grad is the X potential gradient
C
      gradx(kk)=(bb(1)*phi(1)+bb(2)*phi(2)+bb(3)*phi(3))/ar2
C
C-----The j component of Grad is the Y potential gradient
C
      grady(kk)=(cc(1)*phi(1)+cc(2)*phi(2)+cc(3)*phi(3))/ar2
C
C-----The magnitude is the vector sum
C
      ebsq(kk)=gradx(kk)**2+grady(kk)**2
      vsum(kk)=sqrt(ebsq(kk))
      RETURN
      END



