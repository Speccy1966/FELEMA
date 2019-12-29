UNIT PPGRAPH;

{Graph utilities used by PPMAG2D,
 UNIT compiled by DJB, 16/1/93}

INTERFACE

USES CRT,
     GLOBALS,
     GRAPH;

PROCEDURE SetUpScreen;
PROCEDURE Box;
PROCEDURE RestCol;
PROCEDURE OpenWindow (HL :name);
PROCEDURE CloseWindow;
PROCEDURE DrawTri;
PROCEDURE CenTri;
PROCEDURE CoordTri;
PROCEDURE Triangle;

IMPLEMENTATION

PROCEDURE SetUpScreen;
          VAR
             XYRatio,RealXMax,RealYMax          :REAL;
             X1,Y1                              :INTEGER;
          BEGIN
          RealXMax:=0;
          RealYMax:=0;

          { Auto detect graphics card information. }

          GraphDriver:=Detect;
          InitGraph (GraphDriver,GraphMode,'C:\TP\BGI');
          ErrorCode:=GraphResult;
          IF ErrorCode <> grOk THEN
             BEGIN
             WriteLn ('Graphics Error : ', GraphErrorMsg(ErrorCode));
             WriteLn;
             WriteLn ('PROGRAM ABORTED...');
             Halt(1);
             END; {IF}

          { Pixel aspect ratio and real screen ratio. }

          GetAspectRatio (Xasp,Yasp);
          Aspect:=Xasp/Yasp;
          XMax:=GetMaxX;
          YMax:=GetMaxY;
          XYRatio:=XMax*Aspect/YMax;

          { Area of screen reserved for drawing. }

          X1:=Round(XMax-XMax*0.2);
          Y1:=Round(YMax-YMax*0.2);

          { Co-ordinates for model origin. }

          OrigX:=Round((XMax-X1)/2);
          OrigY:=Round(YMax - 0.12*YMax);
          OldOX:=OrigX;
          OldOY:=OrigY;

          { Determine Maximum values of X and Y for the model. }

          FOR ii:=1 TO np DO
                    BEGIN
                    IF xc[ii]>RealXMax THEN
                         RealXMax:=xc[ii];
                    IF yc[ii]>RealYMax THEN
                         RealYMax:=yc[ii];
                    END; {DO ii}

          { Calculate scaling factors for X and Y dimensions. }

          IF ((RealYMax*XYRatio) > RealXMax) THEN
               BEGIN
               ScalerY:=RealYMax/Y1;
               ScalerX:=RealXMax/((RealXMax/RealYMax)*(Y1/Aspect));
          END {IF}
          ELSE
               BEGIN
               ScalerX:=RealXMax/X1;
               ScalerY:=RealYMax/((RealYMax/RealXMax)*(X1/Aspect));
          END; {ELSE}
          sxorig:=ScalerX;
          syorig:=ScalerY;

          { Calculate coordinate position for screen headings. }

          titleX:=ROUND(XMax/2);
          titleY:=ROUND(0.025*YMax);

          RestoreCrtMode;
          END; {SetUpScreen}

PROCEDURE Box;
          BEGIN
          Rectangle (0,0,XMax,YMax);
          END; {Box}

PROCEDURE RestCol;
          BEGIN
          IF GraMon THEN
                    SetColor(Blue)
                    ELSE
                    SetColor(White);
          END; {RestCol}

PROCEDURE OpenWindow (HL :name);
          BEGIN
          SetGraphMode(GraphMode);
          IF GraMon THEN
             BEGIN
             SetBkColor(White);
             SetColor(Blue);
             END; {IF}
          Box;
          SetTextJustify (CenterText,CenterText);
          SetTextStyle (SmallFont,HorizDir,6);
          X:=TitleX;
          Y:=TitleY;
          OutTextXY (X,Y,HL);
          SetTextStyle (SmallFont,HorizDir,4);
          END; {OpenWindow}

PROCEDURE CloseWindow;
          BEGIN
          SetTextJustify (CenterText,CenterText);
          SetTextStyle (DefaultFont,HorizDir,NormSize);
          X:=XMax-100;
          Y:=YMax-5;
          OutTextXY (X,Y,'Press [CR] to continue..');
          ReadLn;
          RestoreCrtMode;
          END; {CloseWindow}

PROCEDURE DrawTri;
          BEGIN
          SetColor(LightGrey);
          LINE (EX[1],EY[1],EX[2],EY[2]);
          LINE (EX[2],EY[2],EX[3],EY[3]);
          LINE (EX[3],EY[3],EX[1],EY[1]);
          RestCol;
          END; {DrawTri}

PROCEDURE CenTri;
          BEGIN
          X:=ROUND((EX[1]+EX[2]+EX[3])/3);
          Y:=ROUND((EY[1]+EY[2]+EY[3])/3);
          END; {CenTri}

PROCEDURE CoordTri;
          BEGIN
          FOR jj:= 1 TO 3 DO
              BEGIN
              EX[jj]:=ROUND(OrigX+XC[NEL[ii,jj]]/ScalerX);
              EY[jj]:=ROUND(OrigY-YC[NEL[ii,jj]]/ScalerY);
              END; {DO}
          END; {CoordTri}

PROCEDURE Triangle;
          BEGIN
          CoordTri;
          DrawTri;
          CenTri;
          END; {Triangle}

END. {PPGRAPH}