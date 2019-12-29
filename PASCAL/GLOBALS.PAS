UNIT GLOBALS;

{Global variables used by PPMAG2D,
 UNIT compiled by DJB, 12/1/93}

INTERFACE

CONST  {VGA Colours}

       Black              = 0;
       Blue               = 1;
       Green              = 2;
       Cyan               = 3;
       Red                = 4;
       Magenta            = 5;
       Brown              = 6;
       LightGrey          = 7;
       DarkGrey           = 8;
       LightBlue          = 9;
       LightGreen         = 10;
       LightCyan          = 11;
       LightRed           = 12;
       LightMagenta       = 13;
       Yellow             = 14;
       White              = 15;

       {Set/GetTextStyle constants}

       DefaultFont        = 0;
       TriplexFont        = 1;
       SmallFont          = 2;
       SansSerifFont      = 3;
       GothicFont         = 4;

       HorizDir           = 0;
       VertDir            = 1;
       NormSize           = 1;

       {Model size parameters modified
        to max size for DOS (640k) 1/1/93 DJB}

       maxnp		  = 900;
       maxne		  = 1600;
       maxbc		  = 250;
       maxbh              = 12;
       maxld              = 5;
       maxcolors          = 15;
TYPE
       vectint=ARRAY [1..maxne] OF INTEGER;
       vectorA=ARRAY [1..maxnp] OF REAL;
       vectorB=ARRAY [1..maxne] OF REAL;
       name=STRING [30];
       ematrix=ARRAY [1..maxne,1..3] OF INTEGER;
       ldata=ARRAY [1..maxld] OF REAL;
       bcdata=ARRAY [1..maxbc] OF REAL;
       bhdata=ARRAY [1..maxbh] OF REAL;
       bcnod=ARRAY [1..maxbc] OF INTEGER;
       vect=ARRAY [1..3] OF INTEGER;

VAR
       n,i,j,k,ii,jj,kk,
       GraphDriver,GraphMode,ErrorCode,nlchk,
       np,ne,ncoef,nbc,x,y,Xmax,Ymax,origX,origY,
       iplvl,iprint,niter,titleX,titleY,col,iterc,
       rings,OldOX,OldOY,nodo,iaxis,ityp
                                                     :INTEGER;

       ScalerX,ScalerY,Aspect,rhoc,rhonl,ascale,
       modx,mody,sxorig,syorig                       :REAL;

       CurrDens,dx,con                               :ldata;
       xc,yc,PotA                                    :vectorA;
       Bx,By,FluxB                                   :vectorB;
       B,H,Rel,BSq                                   :bhdata;
       nmtl                                          :vectint;
       nel                                           :ematrix;
       bv                                            :bcdata;
       ib                                            :bcnod;
       TxtFile                                       :TEXT;
       datin,title,HeadLine,units                    :name;
       ans,ch,cr                                     :CHAR;
       Xasp,Yasp                                     :WORD;
       ex,ey                                         :vect;
       PREP,GraMon                                   :BOOLEAN;

PROCEDURE ZeroVar;

IMPLEMENTATION

PROCEDURE ZeroVar;
          BEGIN
          GraMon:=False;
          FOR jj := 1 TO maxne DO
                 BEGIN
                 FluxB[jj]:=0.0;
                 Bx[jj]:=0.0;
                 By[jj]:=0.0;
                 nmtl[jj]:=0;
                 FOR kk := 1 TO 4 DO
                        nel[jj,kk]:=0;
                 END; {jj}
          FOR jj := 1 TO maxbc DO
                 BEGIN
                 ib[jj]:=0;
                 bv[jj]:=0.0;
                 END; {jj}
          FOR jj := 1 TO maxld DO
                 CurrDens[jj]:=0.0;
          FOR jj := 1 TO maxnp DO
                 BEGIN
                 PotA[jj]:=0.0;
                 xc[jj]:=0.0;
                 yc[jj]:=0.0;
                 END; {jj}
          FOR jj := 1 TO maxbh DO
                 BEGIN
                 B[jj]:=0.0;
                 H[jj]:=0.0;
                 Rel[jj]:=0.0;
                 BSq[jj]:=0.0;
                 END; {jj}
          END; {ZeroVar}

END. {of UNIT Globals}