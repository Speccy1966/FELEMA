UNIT PPUTILS;

{Utilities used by PPMAG2D,
 UNIT compiled by DJB, 16/1/93}

INTERFACE

USES CRT,GLOBALS;

PROCEDURE InKey (VAR cr :char);
PROCEDURE Check (VAR ans :char);
FUNCTION  IntToStr(i :integer) :STRING;
PROCEDURE Dfile (VAR a :name);

IMPLEMENTATION

PROCEDURE InKey (VAR cr :char);
          BEGIN
          GotoXY (1,25);
          Write ('Press Any Key....');
          REPEAT
             cr:=ReadKey;
          UNTIL cr<>#0;
          END;

PROCEDURE Check (VAR ans :char);
          BEGIN
          REPEAT
             WRITE ('(YES or NO ?) : ');
             ReadLn (ans);
             ch:=UpCase(ans);
          UNTIL (ch='Y') or (ch='N')
          END;

FUNCTION  IntToStr(i :integer) :STRING;
          VAR
            S : STRING[8];
          BEGIN
            Str (i,s);
            IntToStr := s;
          END; {IntToStr}

PROCEDURE Dfile (Var a :name);
          BEGIN
          Write ('     Name of input data file ? : ');
          ReadLn (a);
          END; {DFile}

END. {PPUTILS}
