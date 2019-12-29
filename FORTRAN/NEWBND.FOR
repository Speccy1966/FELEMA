      SUBROUTINE newbnd(kk,nbound,cdens,wgt)
C
C ---------------------------------------------------------------------
C I                              NEWBND                               I
C I								      I
C I   FELEMA vs 1.0		  DJB		     11 / July / 1994 I
C ---------------------------------------------------------------------
C I                    RECALCULATE BOUNDARY VOLTAGE                   I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      COMMON /echem/ nsbe, necbe(maxbc), necb1(maxbc), necb2(maxbc),
     + ecba, ecbb, ecbc, ecmc, eccc, volts, nbel(5), monp
      REAL pot
C
C ---------------------------------------------------------------------
C I  For the cathodic surface potential is calculated from Tafels law I
C ---------------------------------------------------------------------
C 
      pot=(LOG(cdens)-eccc)*(1/ecmc)
      IF (volts.LT.0.0) THEN
         pot=volts+pot
      ELSE
         pot=volts-pot
      END IF
      ib(kk+nbound)=necb1(kk)
      ib(kk+nbound+1)=necb2(kk)
      bv(kk+nbound)=(1-wgt)*bv(kk+nbound)+wgt*pot
      bv(kk+nbound+1)=(1-wgt)*bv(kk+nbound+1)+wgt*pot
      nbc=nbc+1
      RETURN
      END
