      PROGRAM MAIN
C
C ---------------------------------------------------------------------
C I                                                                   I
C I    FFFFF    EEEEE    L        EEEEE    M   M     AAA              I
C I    F        E        L        E        MM MM    A   A             I
C I    FFFF     EEEE     L        EEEE     M M M    AAAAA             I
C I    F        E        L        E        M   M    A   A             I
C I    F        EEEEE    LLLLL    EEEEE    M   M    A   A             I
C I                                                                   I
C ---------------------------------------------------------------------
C
C =====================================================================
C I                                                                   I
C I                             FELEMA                                I
C I                           ----------                              I
C I                                                                   I
C I       A finite element solver for 2D and 2D axisymmetric          I
C I           non linear electromagnetic field problems               I
C I                                                                   I
C I                     Version 1 : 7 July 1994                       I
C I                                                                   I
C I                 Solver is Gaussian Elimination                    I
C I                                                                   I
C I  Four types of Field analysis are possible                        I
C I                                                                   I
C I  1) Electromagnetic field analysis 1/u . del^2 A + J = 0          I
C I     The non linear relationship between u and A is handled        I
C I     by simple substitution and Newton Raphson                     I
C I                                                                   I
C I  2) Electrostatic analysis           e . del^2 PHI + RHO = 0      I
C I                                                                   I
C I  3) Electric current analysis    SIGMA . del^2 PHI = 0            I
C I                                                                   I
C I  4) Electrochemical systems      SIGMA . del^2 PHI = 0            I
C I     The Tafel equation is used to balance the over potential      I
C I     at the electrode surfaces with the current density.           I
C I                                                                   I
C I                                                                   I
C I     Modified from the TDFIELD & MAG2DNL programs by D.J.Burt      I
C I                                                                   I
C I             Risley, Engineering Technology Services.              I
C I                                                                   I
C I                 Rutherford House R201, Ext 2694                   I
C =====================================================================
C
C     Standard i/o for the SGI is stream 5 for READ and 6 for WRITE   
C
C ---------------------------------------------------------------------
C I                                                                   I
C I        DEFINE VARIABLE TYPES, ARRAYS AND COMMON BLOCKS            I
C I                                                                   I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      COMMON /echem/ nsbe, necbe(maxbc), necb1(maxbc), necb2(maxbc),
     + ecba, ecbb, ecbc, ecmc, eccc, volts, nbel(5), monp
      CHARACTER title(20), datin*20, datout*20
      DIMENSION dx(maxld), q(maxld), ick(maxnp), fv(3)
      INTEGER i,j,k,n,ii,jj,kk,ij,nb,ne,ncoef,iplvl,niter,istop,iterc,
     + nid,inbw,jend,iter,iprint,iaxis,ityp,iflag,nbcold,monp
      REAL rhoc,rhonr,sum,sumold
C
C ---------------------------------------------------------------------
C I                                                                   I
C I               DEFINITION OF THE INPUT PARAMETERS                  I
C I                                                                   I
C ---------------------------------------------------------------------
C
C   Title and Parameters
C
C         title  - A descriptive statement of the
C                  problem being solved (max :- 20 chars)
C
C         np     - Number of equations (also number of nodes)
C         ne     - Number of elements
C         ncoef  - Number of sets of load cases
C                  maximum of five
C         nbc    - Number of known nodal values, used for applying
C                  Dirichlet boundary condition
C         nsbe   - Number of special boundary elements, used for 
C                  electrochemical Dirichlet boundaries
C         iplvl  - Print level
C                  0 - Don't write the element matrices
C                  1 - Write the element matrices
C                  Set to zero if ne > 20
C         iprint - Iteration print level
C                  0 - Don't write solution for this iteration
C                  1 - Write solution for this iteration
C
C         niter  - Maximum number of iterations
C
C         iaxis  - Geometry configuration
C                  0 - Use 2d linear triangular elements
C                  1 - Use 2d axisymmetric triangular elements
C
C         ityp   - Type of field analysis
C                  1 - Magnetostatic
C                  2 - Electrostatic
C                  3 - Electric Current Analysis
C                  4 - Electrochemical system
C
C   Material Set Loading
C
C         The number of sets must equal ncoef
C
C         q(i)       - Material loads
C         For ITYP 1 - Current density in the Z direction A/m^2
C                  2 - Electric charge density C/m^2 / unit depth
C                      (unit depth could be metre or radian)
C                  3 - Should be zero for Laplacian formulation
C                  4 - Should be zero for Laplacian formulation
C         dx(i)      - Initial material properties
C         For ITYP 1 - Reluctivity m/H
C                  2 - Permittivity F/m
C                  3 - Conductivity S/m
C                  4 - Conductivity S/m
C
C   Non Linear Material Data
C
C         So far this has only been used for ITYP 1, 17/1/93
C
C           h(i) - Magnetic field intensity  (A/m)
C           b(i) - Magnetic flux density     (Tesla)
C
C          nlchk - Switch for data
C              0 - Use the B-H data and the CONVBH routine
C              1 - Use the input constants (const(3))
C
C       const(i) - Curve fit constants can be calculated
C                  from experimental data with EASYPLOT
C
C   Iteration Controls
C
C         iterc  - Number of preliminary iterations
C         rhoc   - Weighting used for preliminary iterations
C         rhonr  - Weighting used by Newton Raphson iteration
C         iflag  - Set to 1 at start, set to zero 0 when converged
C
C   Nodal Coordinate Values
C
C         xc(i)  - X Coordinates of the nodes
C         yc(i)  - Y Coordinates of the nodes
C                  The coordinates must be in numerical sequence
C                  relative to the node numbers
C
C   Elemant Data
C
C         n        - Element number
C         nmtl     - Integer specifying the load case set
C         nel(n,1) - Numerical value of node I
C         nel(n,2) - Numerical value of node J
C         nel(n,3) - Numerical value of node K
C                    Element connectivity must form a right 
C                    right handed set to work with PPMAG2D  
C
C   Boundary Condition Data
C
C         ib(i) - Node number
C         bv(i) - Known value of potential
C                 (zero for Dirichlet boundary)
C      necbe(i) - Elements with face on electrochemical
C                 reacting surface
C      necb1(i) - Node 1 of electrochemical element face
C      necb2(i) - Node 2 of electrochemical element face
C
C ---------------------------------------------------------------------
C I                                                                   I
C I                     INITIALISE ARRAYS                             I
C I                                                                   I
C ---------------------------------------------------------------------
C
      DO 751 kk=1,maxne
	 erel(kk)=0.0
	 ebsq(kk)=0.0
	 nmtl(kk)=0
	 DO 750 jj=1,3
	    nel(kk,jj)=0
	    ef(kk,jj)=0.0
	    DO 754 ii=1,3
	       esmuni(kk,jj,ii)=0.0
754         CONTINUE
750      CONTINUE
751   CONTINUE
      DO 752 kk=1,maxnp
	 xc(kk)=0.0
	 yc(kk)=0.0
	 ick(kk)=0
	 soln(kk)=0.0
	 cur(kk)=0.0
752   CONTINUE
      DO 753 kk=1,maxbh
	 b(kk)=0.0
	 h(kk)=0.0
	 rel(kk)=0.0
	 bsq(kk)=0.0
753   CONTINUE
      DO 755 kk=1,maxbc
	 ib(kk)=0
	 bv(kk)=0.0
	 necbe(kk)=0
	 necb1(kk)=0
	 necb2(kk)=0         
755   CONTINUE
      DO 756 kk=1,maxld
	 dx(kk)=0.0
	 q(kk)=0.0
756   CONTINUE
      DO 758 kk=1,3
	 DO 757 jj=1,3
	    esm(kk,jj)=0.0
	    ajac(kk,jj)=0.0
757      CONTINUE
	 x(kk)=0.0
	 y(kk)=0.0
	 const(kk)=0.0
	 ns(kk)=0
	 resv(kk)=0.0
	 fv(kk)=0.0
758   CONTINUE
      CALL reseta(a,maxav)
C
C ---------------------------------------------------------------------
C I                                                                   I
C I                DATA INPUT SECTION OF THE PROGRAM                  I
C I                                                                   I
C ---------------------------------------------------------------------
C I           INPUT OF THE TITLE AND CONTROL PARAMETERS               I
C ---------------------------------------------------------------------
C
      WRITE (nscr,5)
5     FORMAT (//36X,'FELEMA'/35X,'--------'//15X,
     + 'A Finite Element Solver for 2D Non Linear'/15X,
     + 'problems taking form of the Poisson equation.')
C
C ---------------------------------------------------------------------
C I                        GET FILE NAMES                             I
C ---------------------------------------------------------------------
C
      WRITE (nscr,7)
7     FORMAT (/'Input file name : ')
      READ (5,8) datin
8     FORMAT (14A)
      WRITE (nscr,9)
9     FORMAT ('Output file name : ')
      READ (nkbd,8) datout
C
C ---------------------------------------------------------------------
C I                        OPEN DATA FILES                            I
C ---------------------------------------------------------------------
C
      OPEN (UNIT=ninp, FILE=datin)
      OPEN (UNIT=noutp, FILE=datout)
      OPEN (UNIT=13, FILE='POT1.INP')
      OPEN (UNIT=14, FILE='POT2.INP')
      OPEN (UNIT=15, FILE='VEC1.INP')
      OPEN (UNIT=16, FILE='VEC2.INP')
      OPEN (UNIT=17, FILE='NL.INP')
      OPEN (UNIT=18, FILE='avs.inp')
      WRITE (noutp,5)
C
C ---------------------------------------------------------------------
C I                    READ TITLE AND PARAMETERS                      I
C ---------------------------------------------------------------------
C
      READ (ninp,6) title
6     FORMAT (20A)
      READ (ninp,*) np,ne,ncoef,nbc,iplvl,iprint,niter,iaxis,ityp
      IF (ne.GT.20) iplvl=0
      IF (iplvl.GT.1) iplvl=0
      IF (iprint.GT.1) iprint=0
C     iplvl=1
C
C ---------------------------------------------------------------------
C I                  PUT HEADINGS INTO PPMAG2D FILES                  I
C ---------------------------------------------------------------------
C
      IF (ityp.EQ.1) THEN
	       WRITE (13,*) 'Magnetic Vector Potential'
	       WRITE (13,*) '(Wb/m)'
	       WRITE (14,*) 'Current Density'
	       WRITE (14,*) '(A/m^2)'
	       WRITE (15,*) 'Magnetic Flux Density'
	       WRITE (15,*) '(Tesla)'
	       WRITE (16,*) 'Not used for Magnetostatics'
	       WRITE (17,*) 'Non Linear Rel-BSq Data'
      END IF
      IF (ityp.EQ.2) THEN
	       WRITE (13,*) 'Electric Potential'
	       WRITE (13,*) '(Volts)'
	       WRITE (14,*) 'Electric Charge Density'
	       WRITE (14,*) '(C/m^2/m)'
	       WRITE (15,*) 'Electric Flux Density'
	       WRITE (15,*) '(C/m^2)'
	       WRITE (16,*) 'Electric Field Strength'
	       WRITE (16,*) '(V/m)'
	       WRITE (17,*) 'Not used for Electrostatics'
      END IF
      IF ((ityp.EQ.3).OR.(ityp.EQ.4)) THEN
	       WRITE (13,*) 'Electric Potential'
	       WRITE (13,*) '(Volts)'
	       WRITE (14,*) 'Not used for Electric Current Analysis'
	       WRITE (15,*) 'Current Density'
	       WRITE (15,*) '(A/cm^2)'
	       WRITE (16,*) 'Not used for Electric Current Analysis'
	       WRITE (17,*) 'Not used for Electric Current Analysis'
      END IF
C
C ---------------------------------------------------------------------
C I             CHECK OF NP, NE, NCOEF AND NBC AGAINST                I
C I             THE VALUES IN THE PARAMETER STATEMENT                 I
C ---------------------------------------------------------------------
C
      istop=0
      IF (np.GT.maxnp) THEN
	     WRITE (noutp,10) maxnp
	     istop=1
      END IF
10    FORMAT (10X,'Number of nodes exceeds ',I4/
     +/10X,'** INPUT TERMINATED **')
      IF (ne.GT.maxne) THEN
	     WRITE (noutp,20) maxne
	     istop=1
      END IF
20    FORMAT (10X,'Number of elements exceeds ',I4/
     +/10X,'** INPUT TERMINATED **')
      IF (nbc.GT.maxbc) THEN
	     WRITE (noutp,30) maxbc
	     istop=1
      END IF
30    FORMAT (10X,'Boundary condition data'/
     + 'exceeds dimension of ',I4/10X,'** INPUT TERMINATED **')
      IF (ncoef.GT.maxld) THEN
	     WRITE (noutp,50) maxld
	     istop=1
      END IF
50    FORMAT (10X,'ncoef exceeds ',I3/
     +/10X,'** INPUT TERMINATED **')
      IF (istop.EQ.1) STOP
C
C ---------------------------------------------------------------------
C I    INPUT OF EQUATION COEFFICIENTS AND THE NODAL COORDINATES       I
C ---------------------------------------------------------------------
C
      READ (ninp,*) (dx(i),q(i), i=1,ncoef)
      READ (ninp,*) (h(i),b(i), i=1,maxbh)
      READ (ninp,*) nlchk
      READ (ninp,*) (const(i), i=1,3)
      READ (ninp,*) iterc,rhoc,rhonr
      READ (ninp,*) (xc(i), i=1,np)
      READ (ninp,*) (yc(i), i=1,np)
C
C ---------------------------------------------------------------------
C I               OUTPUT OF THE TITLE AND PARAMETERS                  I
C ---------------------------------------------------------------------
C
      WRITE (noutp,60) title,np,ne,iplvl,iprint,niter,iaxis,ityp
60    FORMAT (//10X,20A//10X,'np = ',I5/10X,'ne = ',I5/
     + 10X,'iplvl = ',I2/10X,'iprint = ',I2/10X,'niter = ',I5/
     + 10X,'iaxis = ',I2/10X,'ityp = ',I2)
C
C ---------------------------------------------------------------------
C I               OUTPUT OF THE MATERIAL SET LOADING                  I
C ---------------------------------------------------------------------
C
      WRITE (noutp,70)
70    FORMAT (//10X,'Material Set Loading'
     + //13X,'Set',8X,'dx',9X,'q')
      WRITE (noutp,80) (i,dx(i),q(i),i=1,ncoef)
80    FORMAT (14X,I2,2E15.5)
C
C ---------------------------------------------------------------------
C I           IF MAGNETIC PROBLEM SET UP NON LINEAR DATA              I
C ---------------------------------------------------------------------
C
      IF (ityp.eq.1) CALL magnet(nlchk)
C
C ---------------------------------------------------------------------
C I                  OUTPUT OF ITERATION CONTROLS                     I
C ---------------------------------------------------------------------
C
      WRITE (noutp,87) iterc,rhoc,rhonr
87    FORMAT (//10X,'Iteration Controls'//10X,'iterc = ',I5
     + /10X,'rhoc = ',F5.2/10X,'rhonr= ',F5.2)
C
C ---------------------------------------------------------------------
C I                OUTPUT OF THE NODAL COORDINATES                    I
C ---------------------------------------------------------------------
C
      WRITE (noutp,90)
90    FORMAT (//10X,'Nodal Coordinates'/10X,
     + 'NODE',5X,'x',14X,'y')
      WRITE (noutp,100) (i,xc(i),yc(i),i=1,np)
100   FORMAT (10X,I4,2E15.5)
C
C ---------------------------------------------------------------------
C I        INPUT AND ECHO PRINT OF THE ELEMENT NODAL DATA             I
C I   AND SET UP INITIAL VALUES FOR ELEMENT MATERIAL PROPERTIES       I
C ---------------------------------------------------------------------
C
      WRITE (noutp,110) title
110   FORMAT (//10X,20A//10X,'Element Data'/
     + 15X,'nel',4X,'nmtl',4X,'Node Numbers')
      nid=0
      DO 701 kk=1,ne
	 READ (ninp,*) n,nmtl(kk),(nel(n,i),i=1,3)
	 IF ((n-1).NE.nid) WRITE (noutp,120) n
	 nid=n
	 WRITE (noutp,130) n,nmtl(kk),(nel(n,i),i=1,3)
	 ii=nmtl(kk)
	 erel(kk)=dx(ii)
701   CONTINUE
120   FORMAT (10X,'ELEMENT',I4,' NOT IN SEQUENCE')
130   FORMAT (15X,I3,5X,I3,2X,3I4)
C
C ---------------------------------------------------------------------
C I      INPUT AND ECHO PRINT OF THE BOUNDARY CONDITION DATA          I
C ---------------------------------------------------------------------
C
      READ (ninp,*) (ib(i), i=1,nbc)
      READ (ninp,*) (bv(i), i=1,nbc)
      WRITE (noutp,140)
140   FORMAT (//10X,'Dirichlet Boundary',/10X,'Known values of A')
      WRITE (noutp,150) (ib(i), bv(i), i=1,nbc)
150   FORMAT (10X,I3,E15.5)
C
C ---------------------------------------------------------------------
C I      INPUT ANY ELCTROCHEMICAL BOUNDARY DATA                       I
C ---------------------------------------------------------------------
C
      IF (ityp.EQ.4) THEN
	 READ (ninp,*) nsbe
	 READ (ninp,*) volts,ecmc,eccc
	 READ (ninp,*) ecba,ecbb,ecbc
	 READ (ninp,*) (nbel(i), i=1,5)
	 READ (ninp,*) monp
	 READ (ninp,*) (necbe(i),necb1(i),necb2(i),i=1,nsbe)
      END IF
      CLOSE (UNIT=ninp)
C
C ---------------------------------------------------------------------
C I                                                                   I
C I                  ANALYSIS OF THE NODE NUMBERS                     I
C I                                                                   I
C ---------------------------------------------------------------------
C I          CHECK TO SEE IF ANY NODE NUMBERS EXCEED NP               I
C ---------------------------------------------------------------------
C
      DO 704 i=1,ne
	 DO 703 j=1,3
	    k=nel(i,j)
	    ick(k)=1
	    IF (k.GT.np) WRITE (noutp,160) j,i,np
703      CONTINUE
704   CONTINUE
160   FORMAT (/10X,'NODE',I4,' OF ELEMENT',I4,
     + ' EXCEEDS NP =',I4)
C
C ---------------------------------------------------------------------
C I    CHECK TO SEE IF ALL NODE NUMBERS THROUGH NP ARE INCLUDED       I
C ---------------------------------------------------------------------
C
      DO 705 i=1,np
	 IF(ick(i).EQ.0) WRITE (noutp,170) i
705   CONTINUE
170   FORMAT(/10X,'NODE',I4,' DOES NOT EXIST')
C
C ---------------------------------------------------------------------
C I                                                                   I
C I          CREATION AND INITIALISATION OF THE A VECTOR              I
C I                                                                   I
C ---------------------------------------------------------------------
C I                 CALCULATION OF THE BAND WIDTH                     I
C ---------------------------------------------------------------------
C
      inbw=0
      nbw=0
      DO 709 kk=1,ne
	 DO 706 i=1,3
	    ns(i)=nel(kk,i)
706      CONTINUE
	 DO 708 i=1,2
	    ij=i+1
	    DO 707 j=ij,3
	       nb=iabs(ns(i)-ns(j))
	       IF (nb.EQ.0) WRITE (noutp,180) kk
	       IF (nb.GT.nbw) THEN
		  inbw=kk
		  nbw=nb
	       END IF
707         CONTINUE
708      CONTINUE
709   CONTINUE
180   FORMAT (/10X,'ELEMENT',I3,' HAS TWO NODES WITH'/
     + 10X,'THE SAME NODE NUMBER')
      nbw=nbw+1
      WRITE (noutp,190) nbw,inbw
190   FORMAT (//10X,'Bandwidth is',I4,' in element',I4)
C
C ---------------------------------------------------------------------
C I                   CALCULATION OF POINTERS                         I
C ---------------------------------------------------------------------
C
      jgf=np
      jgsm=jgf+np
      jend=jgsm+np*nbw
      IF (jend.GT.maxav) THEN
	      WRITE (noutp,200)
	      STOP
      END IF
200   FORMAT (10X,'DIMENSION OF A VECTOR EXCEEDED')

C
C ---------------------------------------------------------------------
C I                                                                   I
C I                           MAIN LOOP                               I
C I                                                                   I
C ---------------------------------------------------------------------
C I  PERFORM A MAXIMUM OF 'NITER' ITERATIONS                          I
C ---------------------------------------------------------------------
C 
      nbcold=nbc
      iflag=1
      IF (ityp.EQ.1) THEN
	 WRITE (nscr,430)
	 WRITE (noutp,430)
      END IF
      IF (ityp.EQ.4) THEN
	 WRITE (nscr,440)
	 WRITE (noutp,440)
      END IF
 430  FORMAT (/12X,'ITER',9X,'SUMWTD',2X,'MONITOR POINT')
 440  FORMAT (/12X,'ITER',3X,'NBC',7X,'RESIDUAL',2X,
     + 'MONITOR POINT')
      DO 770 iter=1,niter
       IF (iflag.EQ.1) THEN 
C
C ---------------------------------------------------------------------
C I              INITIALISATION OF THE COLUMN VECTOR A                I
C ---------------------------------------------------------------------
C
	 CALL reseta(a,jend)
	 sumwtd=0.0
	 sum=0.0
C
C ---------------------------------------------------------------------
C I             GENERATION OF THE SYSTEM OF EQUATIONS                 I
C ---------------------------------------------------------------------
C
	 DO 714 kk=1,ne
C
C ---------------------------------------------------------------------
C I        RETRIEVAL OF NODAL COORDINATES AND NODE NUMBERS            I
C ---------------------------------------------------------------------
C
	    IF (iter.EQ.1) THEN
	       DO 711 i=1,3
		  ns(i)=nel(kk,i)
		  j=ns(i)
		  x(i)=xc(j)
		  y(i)=yc(j)
711            CONTINUE
C
C ---------------------------------------------------------------------
C I             SET MATERIAL COEFFICIENT AND LOAD CASE                I
C ---------------------------------------------------------------------
C
	       dxe=erel(kk)
	       ii=nmtl(kk)
	       qe=q(ii)
C
C ---------------------------------------------------------------------
C I          CALCULATION OF THE ELEMENT STIFFNESS MATRIX              I
C I                 AND THE ELEMENT FORCE VECTOR                      I
C ---------------------------------------------------------------------
C

	       IF (iaxis.EQ.0) CALL elstmf(kk,iplvl,iprint)
	       IF (iaxis.EQ.1) CALL alstmf(kk,iplvl,iprint)
	    END IF
C
C ---------------------------------------------------------------------
C I     RECALCULATE ELEMENT STIFFNESS MATRIX WITH NEW RELUCTIVITY     I
C I     ONLY FOR MAGNETIC SYSTEMS ityp=1                              I     
C ---------------------------------------------------------------------
C
	    IF ((iter.GT.1).AND.(ityp.EQ.1)) THEN
	       DO 832 i=1,3
		  DO 830 j=1,3
		     esm(i,j)=erel(kk)*esmuni(kk,i,j)
830               CONTINUE
832            CONTINUE
	    END IF
	    DO 833 i=1,3
	       ns(i)=nel(kk,i)
	       fv(i)=ef(kk,i)
833         CONTINUE
	    CALL dirstif(ns,fv,esm)
714      CONTINUE
C
C ---------------------------------------------------------------------
C I      MODIFICATION AND SOLUTION OF THE SYSTEM OF EQUATIONS         I
C ---------------------------------------------------------------------
C
	 CALL modify
	 CALL dcmpbd
	 CALL slvbd
C
C ---------------------------------------------------------------------
C I      STORE NODAL POTENTIALS AND FORCING FUNCTION VALUES           I
C I      THEN RESET STORAGE VECTOR                                    I
C ---------------------------------------------------------------------
C
	 DO 840 i=1,np
	    soln(i)=a(i)
	    j=i+np
	    cur(i)=a(j)
	    sum=sum+a(i)
	    IF (iter.EQ.1) sumold=sum+1
840      CONTINUE
	 CALL reseta(a,jend)
C
C ---------------------------------------------------------------------
C I      IF MODEL IS MAGNETIC ITERATE TO UPDATE RELUCTIVITIES         I
C ---------------------------------------------------------------------
C
	 IF (ityp.EQ.1) THEN 
	    CALL miter(ne,iter,niter,iterc,iplvl,rhoc,rhonr,iflag)
	 END IF
C
C ---------------------------------------------------------------------
C I      IF MODEL IS ELECTROCHEMICAL ITERATE TO UPDATE DIRICHLET      I
C I      BOUNDARY CONDITIONS                                          I
C ---------------------------------------------------------------------
C
	 IF (ityp.EQ.4) THEN
	    nbc=nbcold
	    CALL eliter(ne,iter,niter,iflag,sum,sumold,rhoc)
	 END IF
C
C ---------------------------------------------------------------------
C I      FOR ELECTROSTATIC AND ELECTRIC CURRENT ANALYSIS ONLY         I
C I      ONE ITERATION IS REQUIRED                                    I
C ---------------------------------------------------------------------
C
	 IF ((ityp.EQ.2).OR.(ityp.EQ.3)) THEN
	    iflag=0
	    DO 775 kk=1,ne
	       CALL calcgrad(kk)
C
C ---------------------------------------------------------------------
C I      Gradients are converted to current density                   I
C ---------------------------------------------------------------------
C
	       gradx(kk)=-erel(kk)*gradx(kk)
	       grady(kk)=-erel(kk)*grady(kk)
 775        CONTINUE
	 END IF
       END IF
       sumold=sum
 770  CONTINUE
C
C ---------------------------------------------------------------------
C I                     SEND RESULTS TO OUTPUT                        I
C ---------------------------------------------------------------------
C
      CALL output(title,ne,ityp)
C
C ---------------------------------------------------------------------
C I                   CALL SUBROUTINE FOR AVS OUTPUT                  I
C ---------------------------------------------------------------------
C
      CALL avsinp(ne,ityp)
C
C ---------------------------------------------------------------------
C I                   CLOSE UP AND FINISH                             I
C ---------------------------------------------------------------------
C
      DO 740 i=12,18
	 CLOSE (UNIT=i)
740   CONTINUE
      WRITE (nscr,390)
390   FORMAT (/10X,'RUN COMPLETE..')
      STOP
      END












