PROGRAM FESTIVE (INPUT,OUTPUT);
{ ---------------------------------------------------------------------
  I                                                                   I
  I    FFFFF    EEEEE     SSSS    TTTTT     III    V     V   EEEEE    I
  I    F        E        S          T        I     V     V   E        I
  I    FFFF     EEEE      SSS       T        I      V   V    EEEE     I
  I    F        E            S      T        I       V V     E        I
  I    F        EEEEE    SSSS       T       III       V      EEEEE    I
  I                                                                   I
  ---------------------------------------------------------------------

  =====================================================================

                                 FESTIVE
                                ---------

               A simple pre and post processor facility for

                     the Finite Element code FELEMA


                             14 / July / 1994

                            Written by D.J.Burt

        Salford University, Department of Pure and Applied Physics.


  =====================================================================}

USES PPUTILS,
     PPGRAPH,
     TITLES,
     CRT,
     GRAPH,
     GLOBALS;

PROCEDURE ReadNL (VAR c: name); {Read any non linear data}
          BEGIN
          Assign (TxtFile,c);
          Reset (TxtFile);
          FOR ii:=1 TO maxbh DO
                    ReadLn (TxtFile,BSq[ii],Rel[ii]);
          FOR ii:=1 TO 3 DO
                    ReadLn (TxtFile,con[ii]);
          Close (TxtFile);
          END; {ReadNL}

PROCEDURE ReadData;

          VAR
             Fini       :BOOLEAN;

PROCEDURE ReadINP;
          BEGIN
          DFile (Datin);
          Assign (TxtFile,Datin);
          Reset (TxtFile);
          ReadLn (TxtFile,title);
          WriteLn;
          WriteLn ('     ',title);
          ReadLn (TxtFile,np,ne,ncoef,nbc,iplvl,iprint,niter,iaxis,ityp);
          FOR ii:=1 TO NCoef DO
                    ReadLn (TxtFile,dx[ii],CurrDens[ii]);
          FOR ii:=1 TO maxbh DO
                    ReadLn (TxtFile,H[ii],B[ii]);
          ReadLn (TxtFile,nlchk);
          ReadLn (TxtFile,con[1],con[2],con[3]);
          ReadLn (TxtFile,iterc,rhoc,rhonl);
          FOR ii:=1 TO np DO
                    Read (TxtFile,xc[ii]);
          FOR ii:=1 TO np DO
                    Read (TxtFile,yc[ii]);
          FOR ii:=1 TO NE DO
                    BEGIN
                    Read (TxtFile,n,nmtl[ii]);
                         FOR jj:=1 TO 3 DO
                                   Read (TxtFile,NEL[ii,jj]);
                    END; {DO ii}
          FOR ii:=1 TO NBC DO
                    Read (TxtFile,IB[ii]);
          FOR ii:=1 TO NBC DO
                    Read (TxtFile,BV[ii]);
          WriteLn ('     Input data read successfully.');
          Close (TxtFile);

          END; {ReadINP}

          BEGIN
          Fini:=False;
          ClrScr;
          Banner;
          WriteLn;
          REPEAT
                WriteLn ('     1)  PRE - PROCESS ');
                WriteLn;
                WriteLn ('     2)  POST - PROCESS ');
                GotoXY(1,8);
                Write ('     Choose option (1-2) : ');
                ReadLn (i);
                WriteLn;
                IF (i>0) AND (i<=2) THEN
                     CASE i OF
                         1: BEGIN
                            ReadINP;
                            PREP:=true;
                            Fini:=true;
                            InKey(cr);
                            END;
                         2: BEGIN
                            ReadINP;
                            Fini:=true;
                            InKey(cr);
                            END;
                     END {Case}
                     ELSE
                            BEGIN
                            GotoXY(1,8);
                            WriteLn ('                           ');
                            END; {Else}
          UNTIL Fini
          END; {ReadData}

PROCEDURE ChooseColour (i: integer);
          BEGIN
          CASE i OF
              1: Col:=Blue;
              2: Col:=LightBlue;
              3: Col:=Cyan;
              4: Col:=LightCyan;
              5: Col:=Green;
              6: Col:=LightGreen;
              7: Col:=Yellow;
              8: Col:=Brown;
              9: Col:=Red;
              10: Col:=LightRed;
              END; {Case}
          END; {ChooseColour}

PROCEDURE KeyScale (Max :REAL);

          VAR
          DelX,DelY,X1,Y1,X2,Y2                 :INTEGER;
          MinStr,MidStr,MaxStr                  :STRING[8];

PROCEDURE NumberScale(NameStr :STRING; StringLen :INTEGER);

          BEGIN
          CASE StringLen OF
               1: INSERT('0.000',NameStr,1);
               2: INSERT('0.00',NameStr,1);
               3: INSERT('0.0',NameStr,1);
               4: INSERT('0.',NameStr,1);
               5: INSERT('.',NameStr,2);
          END; {Case}
          OutTextXY(X,Y,NameStr);
          END; {NumberScale}

          BEGIN
          DelX:=ROUND((TitleX-OldOX)/10);
          DelY:=ROUND((DelX*Aspect)/2);
          X1:=OldOX;
          Y1:=ROUND((YMax-0.06*YMax)-(DelY/2));
          Y2:=Y1+DelY;
          Y:=ROUND((Y1+Y2)/2);
          SetColor(DarkGrey);
          FOR i:= 1 TO 10 DO
              BEGIN
              X2:=X1+DelX;
              Rectangle (X1,Y1,X2,Y2);
              X:=ROUND((X1+X2)/2);
              ChooseColour(i);
              SetFillStyle(1,Col);
              FloodFill(X,Y,DarkGrey);
              X1:=X2;
              END; {DO}
          RestCol;
          MinStr:='0.0';
          X:=OldOX;
          Y:=ROUND(Y1+(2*DelY));
          OutTextXY(X,Y,'Units : '+Units);
          Y:=ROUND(Y1-DelY);
          OutTextXY(X,Y,MinStr);
          X:=TitleX;
          i:=ROUND(Max*10000);
          MaxStr:=IntToStr(i);
          j:=Length(MaxStr);
          NumberScale (MaxStr,j);
          X:=ROUND((TitleX+OldOX)/2);
          i:=ROUND(Max*10000/2);
          MidStr:=IntToStr(i);
          j:=Length(MidStr);
          NumberScale (MidStr,j);
          END; {KeyScale}

PROCEDURE Bndplot;
          CONST
               maxnbl   = 50;
          TYPE
               lines=ARRAY [1..maxnbl] OF INTEGER;
          VAR
          nlin,X1,Y1,X2,Y2    :INTEGER;
          nod1,nod2           :lines;
          BEGIN
          FOR ii:=1 TO maxnbl DO
                    BEGIN
                    nod1[ii]:=0;
                    nod2[ii]:=0;
                    END; {DO}
          Assign (TxtFile,'bounds.inp');
          Reset (TxtFile);
          ReadLn (TxtFile,nlin);
          FOR ii:=1 TO nlin DO
                    BEGIN
                    ReadLn (TxtFile,nod1[ii],nod2[ii])
                    END; {DO}
          FOR ii:=1 TO nlin DO
                    BEGIN
                    X1:=Round(OrigX+XC[nod1[ii]]/ScalerX);
                    Y1:=Round(OrigY-YC[nod1[ii]]/ScalerY);
                    X2:=Round(OrigX+XC[nod2[ii]]/ScalerX);
                    Y2:=Round(OrigY-YC[nod2[ii]]/ScalerY);
                    SetColor(LightGrey);
                    LINE (X1,Y1,X2,Y2);
                    RestCol;
                    END; {DO}
          END; {Bndplot}

PROCEDURE Nplot;
          VAR
             NodeNum       :BOOLEAN;
          BEGIN
          NodeNum:=FALSE;
          Write ('     Node numbering on ');
          Check (ans);
          If ch='Y' THEN NodeNum:=TRUE;
          HeadLine:='Nodal Plot.';
          OpenWindow (HeadLine);
          FOR ii:=1 TO NP DO
                    BEGIN
                    X:=Round(OrigX+XC[ii]/ScalerX);
                    Y:=Round(OrigY-YC[ii]/ScalerY);
                    PutPixel(X,Y,LightRed);
                    IF NodeNum THEN
                              BEGIN
                              X:=X+8;
                              OutTextXY(X,Y,IntToStr(ii))
                              END {IF}
                    END; {DO}
          CloseWindow;
          END; {Nplot}

PROCEDURE Eplot;
          VAR
             ElemNum,ElemColour,Shade            :BOOLEAN;
          BEGIN
          FOR ii:= 1 TO 3 DO
              BEGIN
              EX[ii]:=0;
              EY[ii]:=0;
              END; {DO}
          ElemNum:=FALSE;
          ElemColour:=FALSE;
          Shade:=FALSE;
          Write ('     Element numbering on ');
          Check (ans);
          If ch='Y' THEN ElemNum:=TRUE;
          Write ('     Element colour on ');
          Check (ans);
          If ch='Y' THEN ElemColour:=TRUE;
          If (ch='Y') AND GraMon THEN
             BEGIN
             ElemColour:=FALSE;
             Shade:=TRUE;
             END; {IF}
          HeadLine:='Element Plot.';
          OpenWindow (HeadLine);
          FOR ii:=1 TO NE DO
                    BEGIN
                    Triangle;
                    IF ElemColour THEN BEGIN
                       j:=NMTL[ii];
                       IF (j>0) AND (j<=4) THEN
                          CASE j OF
                               1:   SetFillStyle(1,Blue);
                               2:   SetFillStyle(1,DarkGrey);
                               3:   SetFillStyle(9,Red);
                               4:   SetFillStyle(9,LightRed);
                          END; {Case}
                          FloodFill(X,Y,LightGrey);
                       END; {IF}
                    IF Shade THEN BEGIN
                       j:=NMTL[ii];
                       IF (j>0) AND (j<=4) THEN
                          CASE j OF
                               1:   SetFillStyle(0,White);
                               2:   SetFillStyle(10,Blue);
                               3:   SetFillStyle(11,Blue);
                               4:   SetFillStyle(11,Blue);
                          END; {Case}
                          FloodFill(X,Y,LightGrey);
                       END; {IF}
                    RestCol;
                    IF ElemNum THEN
                       OutTextXY(X,Y,IntToStr(ii))
                    END; {DO}
          CloseWindow;
          END; {Eplot}


{ ---------------------------------------------------------------------
  I                                                                   I
  I                        PLOT 'A' FIELD.                            I
  I                                                                   I
  ---------------------------------------------------------------------}


PROCEDURE PlPot;

          TYPE
             VectInt3 = ARRAY [1..3] OF INTEGER;
             Vect3 = ARRAY [1..3] OF REAL;
             Vect2 = ARRAY [1..2] OF REAL;
          VAR
             Phi,SF,TA,TB,TC,Zeta                : Vect3;
             xFE,yFE,xSE,ySE,PhiFE,PhiSE         : Vect2;
             num                                 : VectInt3;
             AMax,step,Pot,temp,sum,Delta        : REAL;
             contour                             : INTEGER;
             Fini                                : BOOLEAN;

PROCEDURE TriShape(i,j,k :INTEGER);

          VAR
             s,t,u                               : REAL;
          BEGIN

          s:=EX[j];
          u:=EY[k];
          s:=s*u;
          t:=EX[k];
          u:=EY[j];
          t:=t*u;
          TA[i]:=s-t;
          TB[i]:=EY[j]-EY[k];
          TC[i]:=EX[k]-EX[j];
          END; {TriShape}

PROCEDURE AScreen;
          BEGIN
          OpenWindow (HeadLine);
          KeyScale (AMax);
          END; {AScreen}

PROCEDURE SetTriA;
          BEGIN
          CoordTri;
          FOR jj:= 1 TO 3 DO
              Phi[jj]:=ABS(PotA[NEL[ii,jj]]);
          END; {SetTriA}

PROCEDURE PlotEdgeToEdge;

          VAR
            xx1,yy1,xx2,yy2           : INTEGER;
            ratio                     : REAL;

          BEGIN
          IF (PhiFE[2] < PhiFE[1]) THEN
             BEGIN
             temp:=PhiFE[1];
             PhiFE[1]:=PhiFE[2];
             PhiFE[2]:=temp;
             temp:=xFE[1];
             xFE[1]:=xFE[2];
             xFE[2]:=temp;
             temp:=yFE[1];
             yFE[1]:=yFE[2];
             yFE[2]:=temp;
          END; {IF}
          IF (PhiSE[2] < PhiSE[1]) THEN
             BEGIN
             temp:=PhiSE[1];
             PhiSE[1]:=PhiSE[2];
             PhiSE[2]:=temp;
             temp:=xSE[1];
             xSE[1]:=xSE[2];
             xSE[2]:=temp;
             temp:=ySE[1];
             ySE[1]:=ySE[2];
             ySE[2]:=temp;
          END; {IF}
          Pot:=0.0;
          WHILE (Pot < PhiFE[1]) DO
             Pot:=Pot+step;
          WHILE (Pot < PhiFE[2]) DO
             BEGIN
             n:=ROUND(10*((Pot/AMax)+0.05));
             ChooseColour(n);
             SetColor(Col);
             IF ((Pot >= PhiSE[1]) AND (Pot <= PhiSE[2]) AND
                (PhiFE[2]<>PhiFE[1]) AND (PhiSE[1]<>PhiSE[2])) THEN
                   BEGIN
                   ratio:=(Pot-PhiFE[1])/(PhiFE[2]-PhiFE[1]);
                   xx1:=ROUND(xFE[1]+(xFE[2]-xFE[1])*ratio);
                   yy1:=ROUND(yFE[1]+(yFE[2]-yFE[1])*ratio);
                   ratio:=(Pot-PhiSE[1])/(PhiSE[2]-PhiSE[1]);
                   xx2:=ROUND(xSE[1]+(xSE[2]-xSE[1])*ratio);
                   yy2:=ROUND(ySE[1]+(ySE[2]-ySE[1])*ratio);
                   LINE (xx1,yy1,xx2,yy2);
                   END; {IF}
             Pot:=Pot+step;
             END; {DO}
          RestCol;
          END; {PlotEdgeToEdge}

PROCEDURE PlotTri;

          VAR
             fe,se              :INTEGER;

          BEGIN
          FOR fe:=1 TO 2 DO
                BEGIN
                xFE[1]:=EX[fe];
                yFE[1]:=EY[fe];
                PhiFE[1]:=Phi[fe];
                xFE[2]:=EX[fe+1];
                yFE[2]:=EY[fe+1];
                PhiFE[2]:=Phi[fe+1];
                FOR se:=fe+1 TO 3 DO
                      BEGIN
                      xSE[1]:=EX[se];
                      ySE[1]:=EY[se];
                      PhiSE[1]:=Phi[se];
                      xSE[2]:=EX[se MOD 3 + 1];
                      ySE[2]:=EY[se MOD 3 + 1];
                      PhiSE[2]:=Phi[se MOD 3 + 1];
                      PlotEdgeToEdge;
                      END; {DO}
                END; {DO}
          END; {PlotTri}

PROCEDURE AContPlot;
          VAR
             Boundary, Elems      :BOOLEAN;
          BEGIN
          Boundary:=FALSE;
          Elems:=FALSE;
	  Write ('     Number of contours [',rings,'] ? : ');
          ReadLn (rings);
          Write ('     Plot element boundaries ');
          Check(ans);
          IF ch='Y' THEN Elems:=TRUE;
{         Write ('     Plot material boundaries ');
          Check (ans);
          IF ch='Y' THEN Boundary:=TRUE;  }
          step:=AMax/rings;
          AScreen;
          FOR ii:= 1 TO NE DO
                BEGIN
                SetTriA;
                IF Elems THEN DrawTri;
                PlotTri;
                END; {DO}
          IF Boundary THEN Bndplot;
          CloseWindow;
          END; {AContPlot}

PROCEDURE GetAFile;

          VAR
              infile       : name;

PROCEDURE ReadPot (VAR a: name);  {Reads in a potential field}
          BEGIN
          Assign (TxtFile,a);
          Reset (TxtFile);
          ReadLn (TxtFile,HeadLine);
          ReadLn (TxtFile,Units);
          FOR ii:=1 TO np DO
                    ReadLn (TxtFile,n,PotA[ii]);
          Close (TxtFile);
          GotoXY(1,15);
          WriteLn ('     File read successfully. ')
          END; {ReadPot}

PROCEDURE Choice;
          BEGIN
          GotoXY(1,17);
          Write ('     Choose Option (1-2) : ');
          ReadLn(i);
          WriteLn;

          IF (i>0) AND (i<=2) THEN

                   CASE i OF

                        1: BEGIN
                           infile:='POT1.INP';
                           ReadPot (infile);
                           END;
                        2: BEGIN
                           infile:='POT2.INP';
                           ReadPot (infile);
                           END;

                   END {Case}
                   ELSE
                           BEGIN
                           ScrubLine;
                           Choice;
                           END; {Else}
          END; {Choice}

          BEGIN
          TopHead;
             CASE ityp OF
             1: BEGIN
                WriteLn ('     1) MAGNETIC VECTOR POTENTIAL. ');
                WriteLn ('     2) CURRENT DENSITY. ');
                Choice;
                END;
             2: BEGIN
                WriteLn ('     1) ELECTRIC POTENTIAL, VOLTAGE. ');
                WriteLn ('     2) ELECTRIC CHARGE DENSITY. ');
                Choice;
                END;
             3: BEGIN
                WriteLn ('     1) ELECTRIC POTENTIAL, VOLTAGE. ');
                InKey (cr);
                infile:='POT1.INP';
                ReadPot(infile);
                END;
             4: BEGIN
                WriteLn ('     1) ELECTRIC POTENTIAL, VOLTAGE. ');
                InKey (cr);
                infile:='POT1.INP';
                ReadPot(infile);
                END;
             END; {case}

          AMax:=0.0;
          FOR ii:=1 TO NP DO
                    BEGIN
                    IF (ABS(PotA[ii])>AMax) THEN
                             AMax:=ABS(PotA[ii]);
                    END; {DO}

          END; {GetAFile}

          BEGIN
          Fini:=False;
          REPEAT
                TopHead;
                WriteLn ('     1) CONTOUR PLOT. ');
                WriteLn ('     2) RETURN TO MAIN MENU.');
                GotoXY(1,17);
                Write ('     Choose Option (1-2) : ');
                ReadLn(i);
                WriteLn;
                IF (i>0) AND (i<=2) THEN

                   CASE i OF

                   1: BEGIN
                      GetAFile;
                      AContPlot;
                      END;

                   2: BEGIN
                      Fini:=True;
                      END;

                END
                ELSE
                      BEGIN
                      ScrubLine;
                      PlPot;
                      END;
          UNTIL Fini
          END; {PlPot}

{ ---------------------------------------------------------------------
  I                                                                   I
  I                       PLOT 'B' FIELD.                             I
  I                                                                   I
  ---------------------------------------------------------------------}

PROCEDURE PlVec;

          CONST
             PI=3.141592654;
          VAR
             BMax,Theta,Beta,A                  :REAL;
             X1,Y1,X2,Y2,ArrLength              :INTEGER;
             Fini                               :BOOLEAN;

PROCEDURE Arrow(Angle :Real);
          BEGIN
          MoveTo(X1,Y1);
          X2:=ROUND(ArrLength*(COS(Angle)));
          Y2:=ROUND(ArrLength*(SIN(Angle)));
                    IF (X1>X) THEN
                              LineRel(-X2,-Y2)
                              ELSE
                              LineRel(X2,Y2);
          END; {Arrow}

PROCEDURE BScreen;
          BEGIN
          OpenWindow (HeadLine);
          KeyScale (BMax);
          END; {BScreen}

PROCEDURE BColour;
          BEGIN
          i:=ROUND(10*((FluxB[ii]/BMax)+0.05));
          ChooseColour(i);
          END; {BColour}

PROCEDURE BVect;
          VAR
             Boundary, Elems      :BOOLEAN;
          BEGIN
          Boundary:=FALSE;
          Elems:=FALSE;
          GotoXY(1,17);
          Write ('     Plot element boundaries ');
          Check(ans);
          IF ch='Y' THEN Elems:=TRUE;
{          Write ('     Plot material boundaries ');
          Check (ans);
          IF ch='Y' THEN Boundary:=TRUE;  }
          ArrLength:=10;
          Theta:=0.0;
          Beta:=20*PI/360;
          Write ('     Arrow scale  [',ascale:3:2,'] ? : ');
          ReadLn (ascale);
          BScreen;
          FOR ii:=1 TO NE DO
              BEGIN
              CoordTri;
              IF Elems THEN DrawTri;
              CenTri;
              Y1:=Y-ROUND(ascale*By[ii]/ScalerY);
              X1:=X+ROUND(ascale*Bx[ii]/ScalerX);
              IF (X1=X) AND (Y1>Y) THEN
              Theta:= 3*PI/2;
              IF (X1=X) AND (Y1<Y) THEN
              Theta:= PI/2;
              IF (X1>X) OR (X1<X) THEN
              Theta:=ArcTan((Y1-Y)/(X1-X));
              BColour;
              SetColor(Col);
              Line(X,Y,X1,Y1);
              A:=Theta+Beta;
              Arrow(A);
              A:=Theta-Beta;
              Arrow(A);
              RestCol;
              END; {Do}
          IF Boundary THEN Bndplot;
          CloseWindow;
          END; {BVect}

PROCEDURE BModCol;
          BEGIN
          BScreen;
          FOR ii:=1 TO NE DO
              BEGIN
              Triangle;
              BColour;
              SetFillStyle(1,Col);
              FloodFill(X,Y,LightGrey);
              END; {DO}
          CloseWindow;
          END; {BModCol}

PROCEDURE GetBFile;

          VAR
               infile        :name;
               
PROCEDURE ReadVec (VAR b: name);  {Reads in a vector field}
          BEGIN
          Assign (TxtFile,b);
          Reset (TxtFile);
          ReadLn (TxtFile,HeadLine);
          ReadLn (TxtFile,Units);
          FOR ii:=1 TO ne DO
                    ReadLn (TxtFile,n,Bx[ii],By[ii],FluxB[ii]);
          Close (TxtFile);
          GotoXY(1,15);
          WriteLn ('     File read successfully. ')
          END; {ReadVec}

PROCEDURE Choice;
          BEGIN
          GotoXY(1,17);
          Write ('     Choose Option (1-2) : ');
          ReadLn(i);
          WriteLn;

          IF (i>0) AND (i<=2) THEN

                   CASE i OF

                        1: BEGIN
                           infile:='VEC1.INP';
                           ReadVec (infile);
                           END;
                        2: BEGIN
                           infile:='VEC2.INP';
                           ReadVec (infile);
                           END;

                   END {Case}
                   ELSE
                           BEGIN
                           ScrubLine;
                           Choice;
                           END; {Else}
          END; {Choice}

          BEGIN
          TopHead;
             CASE ityp OF
             1: BEGIN
                WriteLn ('     1) MAGNETIC FLUX DENSITY. ');
                InKey (cr);
                infile:='VEC1.INP';
                ReadVec(infile);
                END;
             2: BEGIN
                WriteLn ('     1) ELECTRIC FIELD INTENSITY. ');
                WriteLn ('     2) ELECTRIC FLUX DENSITY. ');
                Choice;
                END;
             3: BEGIN
                WriteLn ('     1) ELECTRIC CURRENT. ');
                InKey (cr);
                infile:='VEC1.INP';
                ReadVec(infile);
                END;
             4: BEGIN
                WriteLn ('     1) ELECTRIC CURRENT. ');
                InKey (cr);
                infile:='VEC1.INP';
                ReadVec(infile);
                END;
             END; {case}

          BMax:=0.0;
          FOR ii:=1 TO NE DO
                    BEGIN
                    IF (ABS(FluxB[ii])>BMax) THEN
                             BMax:=ABS(FluxB[ii]);
                    END; {DO}
{          WriteLn (BMAX);}
          InKey (cr);
          END; {GetBFile}

          BEGIN
          Fini:=False;
          REPEAT
                TopHead;
                WriteLn ('     1) VECTOR PLOT. ');
                WriteLn ('     2) COLOUR ELEMENTS WITH VECTOR MAGNITUDE. ');
                WriteLn ('     3) RETURN TO MAIN MENU.');
                GotoXY(1,17);
                Write ('     Choose Option (1-3) : ');
                ReadLn(i);
                WriteLn;
                IF (i>0) AND (i<=3) THEN

                   CASE i OF

                   1: BEGIN
                      GetBFile;
                      BVect;
                      END;

                   2: BEGIN
                      GetBFile;
                      BModCol;
                      END;

                   3: BEGIN
                      Fini:=True;
                      END;

                END {Case}
                ELSE
                      BEGIN
                      ScrubLine;
                      PlVec;
                      END; {Else}
          UNTIL Fini
          END; {PlVec}

PROCEDURE OptList;

          BEGIN
          TopHead;
          WriteLn ('     1)  PLOT NODES.');
          WriteLn ('     2)  PLOT ELEMENTS.');
          WriteLn ('     3)  POTENTIAL FIELD RESULTS.');
          WriteLn ('     4)  VECTOR FIELD RESULTS.');
          WriteLn ('     5)  ADJUST X AND Y SCALING.');
          WriteLn ('     6)  END FESTIVE.');
          END; {OptList}

PROCEDURE ZOOM;
          BEGIN
          ScalerX:=sxorig;
          ScalerY:=syorig;
          GotoXY(1,18);
          Write ('     Scaler X  [',modx:3:2,'] ? : ');
          ReadLn (modx);
          Write ('     Scaler Y  [',mody:3:2,'] ? : ');
          ReadLn (mody);
          ScalerX:=ScalerX/modx;
          ScalerY:=ScalerY/mody;
          Write ('     Move node to centre screen ');
          Check(ans);
          IF ch='Y' THEN BEGIN
            Write ('     Node number for centre screen [',nodo,'] ? : ');
            ReadLn (nodo);
               IF (nodo>=1) THEN BEGIN
                  OrigY:=ROUND((YMax/2)+(yc[nodo]/ScalerY));
                  OrigX:=ROUND((XMax/2)-(xc[nodo]/ScalerX))
               END;
               IF (nodo=0) THEN BEGIN
                  OrigY:=OldOY;
                  OrigX:=OldOX
               END;
            END; {IF ch}
          END; {zoom}

PROCEDURE NoCanDo;
          BEGIN
          GotoXY(1,18);
          WriteLn('     OPTION NOT AVAILABLE IN PRE-PROCESSOR !');
          InKey(cr);
          END; {NoCanDo}

PROCEDURE Option;

          VAR
             Fini       :BOOLEAN;

          BEGIN
          Fini:=False;
          REPEAT
                OptList;
                GotoXY(1,17);
                Write ('     Choose option (1-6) : ');
                ReadLn (i);
                IF (i>0) AND (i<=6) THEN

                     CASE i OF

                         1: BEGIN
                            NPlot;
                            END;

                         2: BEGIN
                            EPlot;
                            END;

                         3: BEGIN
                            IF PREP THEN
                               NoCanDo
                            ELSE
                               PlPot;
                            END;

                         4: BEGIN
                            IF PREP THEN
                               NoCanDo
                            ELSE
                               PlVec;
                            END;

                         5: BEGIN
                            ZOOM;
                            END;

                         6: BEGIN
                            WriteLn ('     END OF PROGRAM PPMAG2D');
                            Fini:=True;
                            END;

                     END {Case}
                     ELSE
                            BEGIN
                            ScrubLine;
                            Option;
                            END; {Else}
          UNTIL Fini
          END; {Option}



{********** START ********** START ********** START ********** START ******}

BEGIN
PREP:=False;
modx:=1;
mody:=1;
ascale:=0.01;
rings:=10;
nodo:=0;
ZeroVar; {in unit GLOBALS}
ClrScr;
Write ('Do you want b/w reversed ');
Check (ans);
If ch='Y' THEN
   BEGIN
   GraMon:=True;
   END; {IF}
TitleScr; {in unit TITLES}
ReadData;  {in unit READ}
SetUpScreen; {in unit PPGRAPH}
Option;
END. {PPMAG2D}

{*********** END ************ END ************ END ************ END *******}
