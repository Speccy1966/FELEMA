UNIT TITLES;

{Titles used by PPMAG2D,
 UNIT compiled by DJB, 14/7/94}

INTERFACE

USES CRT,
     GLOBALS,
     PPUTILS;

PROCEDURE Banner;
PROCEDURE TitleScr;
PROCEDURE TopHead;
PROCEDURE ScrubLine;

IMPLEMENTATION

PROCEDURE Banner;
          BEGIN
          TextBackground(White);
          TextColor(Red);
          Write ('    FESTIVE :  A Graphical interface for the');
          WriteLn (' Finite Element code,  FELEMA.      ');
          TextBackground(Blue);
          TextColor(White);
          END; {Banner}

PROCEDURE TitleScr;
          BEGIN
          TextBackground(Blue);
          TextColor(White);
          ClrScr;
          WriteLn;
          WriteLn;
          WriteLn ('                                FESTIVE.');
          WriteLn;
          WriteLn ('    *******************************************************************');
          WriteLn ('    *         Finite Element Store Topology and View Everthing        *');
          WriteLn ('    *                                                                 *');
          WriteLn ('    *              A simple post processor facility for               *');
          WriteLn ('    *                                                                 *');
          WriteLn ('    *                the Finite Element code FELEMA                   *');
          WriteLn ('    *                                                                 *');
          WriteLn ('    *                           JULY 1994                             *');
          WriteLn ('    *                                                                 *');
          WriteLn ('    *                      Written by D.J.Burt                        *');
          WriteLn ('    *                                                                 *');
          WriteLn ('    *       BNFL, Engineering Technology Section, R201, Risley.       *');
          WriteLn ('    *******************************************************************');
          InKey (cr);
          END; {Title}

PROCEDURE TopHead;
          BEGIN
          ClrScr;
          Banner;
          WriteLn;
          WriteLn ('     MODEL     :',title);
          WriteLn ('     NODES     :',np);
          WriteLn ('     ELEMENTS  :',ne);
          WriteLn;
          WriteLn;
          END; {TopHead}

PROCEDURE ScrubLine;
          BEGIN
          GotoXY(1,17);
          WriteLn ('                           ');
          END; {ScrubLine}

END. {TITLES}