      SUBROUTINE calcnu(kk,rho)
C
C ---------------------------------------------------------------------
C I                            CALCNU                                 I
C I								      I
C I   FELEMA vs 1.0		  DJB		      14 / Mar / 1994 I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      REAL relold,relnew
C
C ---------------------------------------------------------------------
C I       CALCULATE ELEMENT RELUCTIVITY USING FITTED EQUATION         I
C I           OF FORM :- Rel = C1 * EXP(C2 * BSq) + C3                I
C I             FOR ELEMENTS OF MATERIAL 2 (STEEL).                   I
C I           AT HIGH BSQ A LINEAR RELATIONSHIP IS USED               I
C ---------------------------------------------------------------------
C I     A WEIGHTING FACTOR 'RHO' IS USED TO IMPROVE CONVERGENCE       I
C ---------------------------------------------------------------------
C
      relold=erel(kk)
      IF (ebsq(kk).GT.bsq(maxbh)) THEN
	 relnew=rel(maxbh)+(ebsq(kk)-bsq(maxbh))*grdmax
      ELSE
	 relnew=const(1)*EXP(const(2)*ebsq(kk))+const(3)
      END IF
C
C ---------------------------------------------------------------------
C I     THE MODIFIED RELUCTIVITY VALUE IS A WEIGHTED COMBINATION      I
C I                      OF RELOLD AND RELNEW.                        I
C ---------------------------------------------------------------------
C
      erel(kk)=(1-rho)*relold+rho*relnew
      RETURN
      END
