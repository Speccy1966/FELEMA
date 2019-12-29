      SUBROUTINE eliter(ne,iter,niter,iflag,sum,sumold,wgt)
C
C ---------------------------------------------------------------------
C I                              ELITER                               I
C I                                                                   I
C I PROCESS SOLUTION FOR THE ELEMENT GRADIENTS (J VALUES), USE THESE  I
C I TO DETERMINE THE DIRICHLET BOUNDARIES FOR THE NEXT ITERATION.     I
C I CALCULATE THE ERROR BETWEEN ITERATIONS AND CHECK FOR CONVERGENCE. I
C I                                                                   I
C I                                                                   I
C I   FELEMA vs 1.0               DJB                11 / July / 1994 I
C ---------------------------------------------------------------------
C I              ITERATES FOR ELECTROCHEMICAL SYSTEM                  I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      COMMON /echem/ nsbe, necbe(maxbc), necb1(maxbc), necb2(maxbc),
     + ecba, ecbb, ecbc, ecmc, eccc, volts, nbel(5), monp
      REAL cdens,pot,wgt,emon,diff
      INTEGER nbound, nflag, nlhs, ntop, nrhs, nbot, nanod
      cdens=0.0
      pot=0.0
      nbound=nbc-(nsbe+1)
      diff=0.0
      nflag=1
C
C ---------------------------------------------------------------------
C     The monitor point is the value of potential at node 'monp'
C ---------------------------------------------------------------------
C
      emon=soln(monp)
      DO 200 kk=1,ne
	 CALL calcgrad(kk)
C
C ---------------------------------------------------------------------
C I         Gradients are converted to Current density                I
C ---------------------------------------------------------------------
C
	 gradx(kk)=-erel(kk)*gradx(kk)
	 grady(kk)=-erel(kk)*grady(kk)
	 vsum(kk)=SQRT(gradx(kk)**2+grady(kk)**2)
 200  CONTINUE
C
C ---------------------------------------------------------------------
C    Calculate current and new potential for electrode surface elements
C ---------------------------------------------------------------------
C     The number of dirichlet boundaries is initially set not to 
C     include the reacting surface
C ---------------------------------------------------------------------
C
      nbc=nbound
C
C ---------------------------------------------------------------------
C I         Calculate pointers for electrode surfaces                 I
C ---------------------------------------------------------------------
C
      nlhs=nbel(1)
      ntop=nlhs+nbel(2)
      nrhs=ntop+nbel(3)
      nbot=nrhs+nbel(4)
      nanod=nbot+nbel(5)
C
C----------------------------------------------------------------------
C     Do the cathodic surface first
C----------------------------------------------------------------------
C
C-----LHS of cathode
C
      DO 205 kk=1,nlhs
	 cdens=ABS(gradx(necbe(kk)))
	 CALL newbnd(kk,nbound,cdens,wgt)
 205  CONTINUE
C
C-----TOP of cathode
C
      DO 206 kk=(nlhs+1),ntop
	 cdens=ABS(grady(necbe(kk)))
	 CALL newbnd(kk,nbound,cdens,wgt)
 206  CONTINUE
C
C-----RHS of cathode
C
      DO 207 kk=(ntop+1),nrhs
	 cdens=ABS(gradx(necbe(kk)))
	 CALL newbnd(kk,nbound,cdens,wgt)
 207  CONTINUE
C
C-----BOTTOM of cathode
C
      DO 208 kk=(nrhs+1),(nbot-1)
	 cdens=ABS(grady(necbe(kk)))
	 CALL newbnd(kk,nbound,cdens,wgt)
 208  CONTINUE
C
C----------------------------------------------------------------------
C     Last node of the cathode bottom joins up with first node of LHS
C----------------------------------------------------------------------
C
      kk=nbot
      cdens=ABS(grady(necbe(kk)))
      pot=(LOG(cdens)-eccc)*(1/ecmc)
      IF (volts.LT.0.0) THEN
	 pot=volts+pot
      ELSE
	 pot=volts-pot
      END IF
      ib(kk+nbound)=necb1(kk)
      ib(nbound+1)=necb2(kk)
      bv(kk+nbound)=(1-wgt)*bv(kk+nbound)+wgt*pot
      bv(nbound+1)=(1-wgt)*bv(nbound+1)+wgt*pot
      nbc=nbc+1
      write(noutp,*) nbound, nbc
C
C----------------------------------------------------------------------
C     Do the anodic surface
C----------------------------------------------------------------------
C
C-----Start the count at first anodic boundary element
C
      jj=nbot+1
      DO 210 kk=jj,nanod
C
C--------Calculate the normal component of current density
C
	 cdens=ABS(gradx(necbe(kk)))
C
C--------Use empirical quadratic function to get potential
C
	 pot=(-ecbb+SQRT(ecbb**2-4*ecba*(ecbc-cdens)))/(2*ecba)
	 IF (volts.LT.0.0) pot=-pot
C
C--------Only if the element potential is greater than 0.5
C        do the nodal potentials change, othaerwise they 
C        revert to the Neumann condition
C 
C        IF (ABS(pot).GT.0.5) THEN
	     ib(jj+nbound)=necb1(jj)
	     ib(jj+nbound+1)=necb2(jj)
	     bv(jj+nbound)=(1-wgt)*bv(jj+nbound)+wgt*pot
	     bv(jj+nbound+1)=(1-wgt)*bv(jj+nbound+1)+wgt*pot
C
C--------Incrementing NBC by 1 
C        means only the bottom node is changed
C
	     nbc=nbc+1
	     jj=jj+1
C        END IF
  210 CONTINUE
C
C-----Include the top node only if the potential in the
C     top element is greater than 0.5
C
C     IF (ABS(pot).GT.0.5) 
      nbc=nbc+1
      diff=ABS(sum-sumold)
      DO 220 i=1,nbc
	 write(noutp,*) i, ib(i), bv(i)
 220  CONTINUE
C
C ---------------------------------------------------------------------
C I      WRITE OUT THE ITERATION NUMBER AND CHECK VALUE OF DIFF       I
C ---------------------------------------------------------------------
C
      WRITE (nscr,420) iter,nbc,diff,emon
      WRITE (noutp,420) iter,nbc,diff,emon
420   FORMAT (10X,2I6,E15.5,E12.4)
C
C ---------------------------------------------------------------------
C I                     CHECK FOR CONVERGENCE                         I
C ---------------------------------------------------------------------
C
      IF (iter.EQ.niter) THEN
	  WRITE (nscr,324) iter
	  WRITE (noutp,324) iter
      END IF    
      IF (diff.LT.1.0e-2) THEN
	 WRITE (nscr,430) iter
	 WRITE (noutp,430) iter
	 iflag=0
      END IF
430   FORMAT (/10X,'PROCESS HAS CONVERGED AFTER ',I4,' ITERATIONS')
324   FORMAT (/10X,'PROCESS HAS NOT CONVERGED AFTER ',I4,' ITERATIONS'
     + //10X,'Modify input file and run again.')
      RETURN
      END











