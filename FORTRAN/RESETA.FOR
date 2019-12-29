      SUBROUTINE reseta(a,jend)
C
C ---------------------------------------------------------------------
C I				 RESETA                               I
C I								      I
C I   FELEMA vs 1.0		  DJB		      14 / Mar / 1994 I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      DIMENSION a(maxav)
      DO 710 i=1,jend
	 a(i)=0.0
710   CONTINUE
      RETURN
      END
