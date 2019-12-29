      SUBROUTINE avsinp(ne,ityp)
C
C ---------------------------------------------------------------------
C I                              AVSINP                               I
C I								      I
C I   FELEMA vs 1.0		  DJB		      14 / Mar / 1994 I
C ---------------------------------------------------------------------
C
      INCLUDE 'MAG2DV5.CFG'
      INCLUDE 'common.f'
      REAL gz,zc
      INTEGER ndata,ncdata,mdata,ncom,nscal,nvec,kk
C
C ---------------------------------------------------------------------
C I      WRITE OUT FIRST LINE; NUMBER OF NODES, ELEMENTS & DATA       I
C I      There are 4 data components, 1 scalar at nodes               I
C I      and 3 vector components at cell centres                      I
C ---------------------------------------------------------------------
C
      ndata=1 
      ncdata=3
      mdata=0
      WRITE(18,180) np,ne,ndata,ncdata,mdata
 180  FORMAT (5I8)
C
C ---------------------------------------------------------------------
C I                   WRITE NODAL COORDINATES                         I
C ---------------------------------------------------------------------
C
      zc=0.0
      WRITE(18,190) (i,xc(i),yc(i),zc,i=1,np) 
 190  FORMAT (I4,3E15.5)
C
C ---------------------------------------------------------------------
C I                   WRITE ELEMENT DATA                              I
C ---------------------------------------------------------------------
C
      DO 100 i=1,ne
         WRITE (18,200) i,nmtl(i),(nel(i,j),j=1,3)
 100  CONTINUE
 200  FORMAT (I4,I3,' tri ',3I4)
C
C ---------------------------------------------------------------------
C I                   WRITE DATA SETS                                 I
C I                                                                   I
C ---------------------------------------------------------------------
C
      ncom=1
      nscal=1
      WRITE(18,205) ncom,nscal
 205  FORMAT (2I3)
      IF (ityp.EQ.1) THEN
         WRITE(18,210)
      END IF
 210  FORMAT ('A_Field, Wb/m')
      IF ((ityp.EQ.2).OR.(ityp.EQ.3).OR.(ityp.EQ.4)) THEN
         WRITE(18,216)
      END IF
 216  FORMAT ('Potential, Volts')
C
C ---------------------------------------------------------------------
C I                   WRITE NODAL FIELD DATA                          I
C ---------------------------------------------------------------------
C
      WRITE(18,220) (i,soln(i),i=1,np)
 220  FORMAT (I4,E15.5)
C
C ---------------------------------------------------------------------
C I                   WRITE CELL CENTRE FIELD DATA                    I
C ---------------------------------------------------------------------
C
      ncom=1
      nvec=3
      WRITE(18,225) ncom,nvec
      gz=0.0
 225  FORMAT (2I3)
      IF (ityp.EQ.1) THEN
         WRITE(18,230)
      END IF
 230  FORMAT ('B_Field, Tesla')
      IF (ityp.EQ.2) THEN
         WRITE(18,233)
      END IF
 233  FORMAT ('Flux_Dens, C/m^2')
      IF ((ityp.EQ.3).OR.(ityp.EQ.4)) THEN
         WRITE(18,235)
      END IF
 235  FORMAT ('Curr_Dens, A/m^2')
      DO 1000 kk=1,ne
         WRITE (18,240) kk,gradx(kk),grady(kk),gz
 1000 CONTINUE
 240  FORMAT (I5,3E15.5)
      RETURN
      END








