#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.1416
#define maxnp 50                   /* max no of nodes around perim */
#define maxnr 20                   /* max no of nodes along radius */
#define minnp 10                   /* min no of nodes around perim */
#define minnr 4                    /* min no of nodes along radius */
#define definpdeck "FERRET.DCK"    /* default input deck file name */
#define deflogfile "FERRET.LOG"    /* default log file name */
#define defavsfile "FERRET.INP"    /* default AVS file name */
#define defgrdfile "FERRET.GRD"    /* default grid file name */
#define asmth 1.5                  /* used in elliptic smoother */
#define mat 1                      /* material number, used for AVS file */

/* This is a new version of Ferret which calls the elliptic smoother
   four times, once in each quadrant of the grid. This ensures that
   the main grid diagonals stay where they are and are not swirled
   around the electrode. The elliptic smoother has been put into a
   subroutine instead of in the main code */

/* function prototypes */
void esmooth(int istart, int iend, int esits);
float volemax(float x, float y);
void cmderr(char c);
/* global variables */
FILE *input, *log_f, *grd, *avs;
int np, nr;
float x[maxnp][maxnr], y[maxnp][maxnr];

main(int argc, char *argv[])
{
  int i, ii, j, k1, k2, nnp[4], n1, n2, n3, n4, nelem,
  maxit, elem[2*maxnp*(maxnr-1)][3];
  char c, inpdeck[30], logfile[30], line[100], avsfile[30], grdfile[30];
  float cx[4], cy[4], d, w, l, theta, p[4], perim, ex, ey, r, sx, sy, tx,
  ty, alpha, beta, s, sum, rhoc, volt, cond;
  
  /* Variable usage
     inpdeck         name of file containing input deck
     input           handle for input from inpdeck
     logfile         name of file to use for log
     log_f           handle for output to log
     grdfile         name of grid file
     grd             handle for grdfile
     avsfile         name for AVS file
     avs             handle for avsfile
     line            dummy variable to hold input line during checking
     (cx[i],cy[j])   co-ords of four corners of model/electrode
     d               distance from l.h. edge of model to centre of electrode
     w               width of electrode
     l               length of electrode
     r               half-diagonal length of electrode
     sx, sy, tx, ty  components of half-diagonal of electrode
     theta           angle between electrode and vertical
     s[nr]           s values for radius
     np              number of nodes around perimeter
     nr              number of nodes along radius
     p[i]            length of section i of perimeter (0..3)
     perim           total length of perimeter
     nnp[i]          no of nodes on perimeter section i (0..3)
     x[i][j], y[i][j] co-ords of node (i,j)
     elem[i][j]      node number of jth corner of element i
     esits           no of iterations for elliptic smoother
     rhoc            underrelaxation factor
     volt            measured voltage difference
     cond            electrolyte conductivity
     */
  
  strcpy(inpdeck,definpdeck);
  strcpy(logfile,deflogfile);
  strcpy(grdfile,defgrdfile);
  strcpy(avsfile,defavsfile);
  for (i=1; i<argc; ++i)
    {
      if (argv[i][0] == '-')
	{
	  c = argv[i][1];
	  switch (c) {
	  case 'h':
	    printf("Usage : %s [-h] [-i inpfile] [-l logfile] [-g grdfile] [-a AVSfile]\n",argv[0]);
	    printf("Where : -h prints this help message\n");
	    printf("        -i specifies inpfile as the input deck\n");
	    printf("        -l specifies logfile as the output log file\n");
	    printf("        -g specifies grdfile as the grid file name\n");
	    printf("        -a specifies AVSfile as the AVS UCD file\n");
	    printf("and there MUST be a space between the option and the\n");
	    printf("file name (ie -i input.dck NOT -iinput.dck)\n");
	    printf("\nDefault file names are\n");
	    printf("\tInput deck\t%s\n",definpdeck);
	    printf("\tLog file\t%s\n",deflogfile);
	    printf("\tGrid file\t%s\n",defgrdfile);
	    printf("\tAVS file\t%s\n",defavsfile);
	    exit(0);
	    break;
	  case 'i':
	    if ((++i >= argc) || (argv[i][0] == '-'))
	      cmderr(argv[--i][1]);
	    else
	      strcpy(inpdeck,argv[i]);
	    break;
	  case 'l':
	    if ((++i >= argc) || (argv[i][0] == '-'))
	      cmderr(argv[--i][1]);
	    else
	      strcpy(logfile,argv[i]);
	    break;
	  case 'g':
	    if ((++i >= argc) || (argv[i][0] == '-'))
	      cmderr(argv[--i][1]);
	    else
	      strcpy(grdfile,argv[i]);
	    break;
	  case 'a':
	    if ((++i >= argc) || (argv[i][0] == '-'))
	      cmderr(argv[--i][1]);
	    else
	      strcpy(avsfile,argv[i]);
	    break;
	  default:
	    printf("Unknown option -%c (ignored)\n",c);
	    argc = 0;
	  }
	}
    }

  /* open input deck file */
  if (!(input = fopen(inpdeck,"r")))
    {
      printf("WARNING : File %s does not exist\n",inpdeck);
      exit(1);
    }
  
  /* open log file */
  log_f = fopen(logfile,"w");
  /* write headers to log file */
  fprintf(log_f,"##############################################\n");
  fprintf(log_f,"# F      E       R     R        E       T    #\n");
  fprintf(log_f,"# Finite Element Rapid Response Editing Tool #\n");
  fprintf(log_f,"# Grid generator                             #\n");
  fprintf(log_f,"# Version 1.0                                #\n");
  fprintf(log_f,"# Alan Silver, R201                          #\n");
  fprintf(log_f,"##############################################\n");
  fprintf(log_f,"\nUsing input file %s\n\n",inpdeck);
  
  /* open grid file */
  grd = fopen(grdfile,"w");

  /* open avs file */
  avs = fopen(avsfile,"w");
  fprintf(avs,"# AVS file from finite element grid generator FERRET\n");
  
  /* now read parameters in from input deck */
  /* Layout of input deck :-
     # comment line(s) - OPTIONAL, MUST START WITH A #
     x0  y0   co-ords of point 0, lower l.h. corner
     x1  y1   co-ords of point 1, lower r.h. corner
     x2  y2   co-ords of point 2, upper r.h. corner
     x3  y3   co-ords of point 3, upper l.h. corner
     d        distance from middle of l.h. edge to centre of electrode
     w        width of electrode
     l        length of electrode
     theta    angle of inclination to the vertical (degrees)
     volt     measured voltage difference
     cond     electrolyte conductivity
     rhoc     underrelaxation factor
     np       no of points around perimeter
     nr       no of points along radius
     maxit    maximum number of iterations for solver
     */
  /* read first line and check for comment */
  fscanf(input,"%[^\n]\n",line);
  while (line[0] == '#')
    {
      fprintf(log_f,"%s\n",line);
      fprintf(avs,"%s\n",line);
      fscanf(input,"%[^\n]\n",line);
    }
  
  /* parse line for first co-ordinates */
  sscanf(line,"%f %f\n",&cx[0],&cy[0]);
  fprintf(log_f,"\nINPUT PARAMETERS\n");
  fprintf(log_f,"Co-ordinates of model\nNo\tx\ty\n");
  fprintf(log_f,"%d\t%.3f\t%.3f\n",i,cx[0],cy[0]);
  
  /* now read positions of other corner points */
  for (i=1; i<4; ++i)
    {
      fscanf(input,"%f %f\n",&cx[i],&cy[i]);
      fprintf(log_f,"%d\t%.3f\t%.3f\n",i,cx[i],cy[i]);
    }
  
  /* Error checking on the points */
  if ((cx[0] >= cx[1]) || (cx[3] >= cx[1])
      || (cy[0] >= cy[3]) || (cy[1] >= cy[2]))
    {
      printf("WARNING : Error detected in model co-ordinates.\n");
      printf("Please check that all points are distinct and that\n");
      printf("the points are specified in the following order :-\n");
      printf("\tlower left corner\n");
      printf("\tlower right corner\n");
      printf("\tupper right corner\n");
      printf("\tupper left corner\n");
      printf("Please correct the input file and try again\n");
      fprintf(log_f,"WARNING : Error detected in model co-ordinates.\n");
      fprintf(log_f,"Please check that all points are distinct and that\n");
      fprintf(log_f,"the points are specified in the following order :-\n");
      fprintf(log_f,"\tlower left corner\n");
      fprintf(log_f,"\tlower right corner\n");
      fprintf(log_f,"\tupper right corner\n");
      fprintf(log_f,"\tupper left corner\n");
      fprintf(log_f,"Please correct the input file and try again\n");
      exit(1);
    }
  
  /* next read d, w, l, theta and rhoc */
  fscanf(input,"%f\n%f\n%f\n%f\n%f\n%f\n%f\n",&d,&w,&l,&theta,&volt,&cond,&rhoc);
  fprintf(log_f,"\nElectrode parameters\n");
  fprintf(log_f,"d     = %.3f m\nw     = %.3f m\nl     = %.3f m\n",d,w,l);
  fprintf(log_f,"theta = %.3f degrees\nvolt  = %.3f\ncond  = %.3f\nrhoc  = %.3f\n",
	  theta,volt,cond,rhoc);
  fprintf(log_f,"\nGrid parameters\n");
  /* convert theta to radians */
  theta = theta*pi/180;
  /* finally, read np, nr and maxit */
  fscanf(input,"%d\n%d\n%d\n",&np,&nr,&maxit);
  fprintf(log_f,"np    = %d\nnr    = %d\n",np,nr);
  if (np < minnp)
    {
      fprintf(log_f,"WARNING : number of perimeter nodes too small\n");
      fprintf(log_f,"Minimum allowed = %d, number set = %d\n",minnp,np);
      fprintf(log_f,"Please correct and rerun FERRET\n");
      printf("WARNING : number of perimeter nodes too small\n");
      printf("Minimum allowed = %d, number set = %d\n",minnp,np);
      printf("Please correct and rerun FERRET\n");
      exit(1);
    }
  if (np>maxnp)
    {
      fprintf(log_f,"WARNING : number of perimeter nodes too large\n");
      fprintf(log_f,"Maximum allowed = %d, number set = %d\n",maxnp,np);
      fprintf(log_f,"Please correct and rerun FERRET\n");
      printf("WARNING : number of perimeter nodes too large\n");
      printf("Maximum allowed = %d, number set = %d\n",maxnp,np);
      printf("Please correct and rerun FERRET\n");
      exit(1);
    }
  if (nr < minnr)
    {
      fprintf(log_f,"WARNING : number of radius nodes too small\n");
      fprintf(log_f,"Minimum allowed = %d, number set = %d\n",minnr,nr);
      fprintf(log_f,"Please correct and rerun FERRET\n");
      printf("WARNING : number of radius nodes too small\n");
      printf("Minimum allowed = %d, number set = %d\n",minnr,nr);
      printf("Please correct and rerun FERRET\n");
      exit(1);
    }
  if (nr>maxnr)
    {
      fprintf(log_f,"WARNING : number of radius nodes too large\n");
      fprintf(log_f,"Maximum allowed = %d, number set = %d\n",maxnr,nr);
      fprintf(log_f,"Please correct and rerun FERRET\n");
      printf("WARNING : number of radius nodes too large\n");
      printf("Maximum allowed = %d, number set = %d\n",maxnr,nr);
      printf("Please correct and rerun FERRET\n");
      exit(1);
    }
  fprintf(log_f,"\nMaximum number of iterations for solver = %d\n",maxit);
  
  /* calculate side lengths */
  p[0] = sqrt((cx[0]-cx[1])*(cx[0]-cx[1])+(cy[0]-cy[1])*(cy[0]-cy[1]));
  p[1] = sqrt((cx[1]-cx[2])*(cx[1]-cx[2])+(cy[1]-cy[2])*(cy[1]-cy[2]));
  p[2] = sqrt((cx[2]-cx[3])*(cx[2]-cx[3])+(cy[2]-cy[3])*(cy[2]-cy[3]));
  p[3] = sqrt((cx[3]-cx[0])*(cx[3]-cx[0])+(cy[3]-cy[0])*(cy[3]-cy[0]));
  perim = p[0] + p[1] + p[2] + p[3];
  
  /* now calculate number of nodes on each section of perimeter */
  /* a neater grid will be formed if perim is an integer and
     divides exactly into np */
  for (i=0; i<3; ++i)
    nnp[i] = p[i]*(float)np/perim;
  /* calculate nnp[3] from others to avoid rounding errors */
  nnp[3] = np - (nnp[0] + nnp[1] + nnp[2]);
  fprintf(log_f,"\nCALCULATED QUANTITIES\n");
  fprintf(log_f,"Number of nodes on outside edges\n");
  fprintf(log_f,"edge\tnnodes\n");
  for (i=0; i<4; ++i)
    fprintf(log_f,"%d\t%d\n",i,nnp[i]);
  
  /******************************************************/
  /* now calculate grid nodes around the perimeter      */
  fprintf(log_f,"\nCalculating node positions around the perimeter ...\n");
  j = 0;
  for (k1=0; k1<4; ++k1)
    {
      fprintf(log_f,"\tPerimeter edge %d\n",k1);
      for (i=0; i<nnp[k1]; ++i)
	{
	  ii = i;
	  for (k2=0; k2<k1; ++k2)
	    ii += nnp[k2];
	  k2 = k1 + 1;
	  if (k2 == 4)
	    k2 = 0;
	  x[ii][j] = cx[k1] + (i*(cx[k2]-cx[k1])/nnp[k1]);
	  y[ii][j] = cy[k1] + (i*(cy[k2]-cy[k1])/nnp[k1]);
	}
    }
  
  /******************************************************/
  /* so now the outside perimeter nodes are set, do the */
  /* inside perimeter                                   */
  fprintf(log_f,"\nCalculating co-ordinates of the electrode ...\n");
  /* first calculate ex and ey (centre of the electrode)*/
  ex = cx[3] - cx[0] + d;
  ey = (cy[3] - cy[0])/2;
  fprintf(log_f,"\tcentre of electrode = (%.3f,%.3f)\n",ex,ey);
  
  /* now calc various lengths and angles needed */
  r = sqrt(w*w + l*l)/2;
  alpha = atan(w/l);
  beta = pi/2 - theta;
  
  sx = r*cos(alpha + beta);
  sy = r*sin(alpha + beta);
  tx = r*cos(beta - alpha);
  ty = r*sin(beta - alpha);
  
  /* now recalculate cx and cy - now the corners of electrode */
  cx[0] = ex - tx;
  cy[0] = ey - ty;
  cx[1] = ex - sx;
  cy[1] = ey - sy;
  cx[2] = ex + tx;
  cy[2] = ey + ty;
  cx[3] = ex + sx;
  cy[3] = ey + sy;
  fprintf(log_f,"\tCo-ordinates of corners of electrode\n");
  fprintf(log_f,"\tNo\tx\ty\n");
  for (i=0; i<4; ++i)
    fprintf(log_f,"\t%d\t%.3f\t%.3f\n",i,cx[i],cy[i]);
  
  /* now calculate grid nodes around the electrode */
  fprintf(log_f,"\nCalculating node positions around the electrode ...\n");
  j = nr-1;
  for (k1=0; k1<4; ++k1)
    {
      fprintf(log_f,"\tElectrode edge %d\n",k1);
      for (i=0; i<nnp[k1]; ++i)
	{
	  ii = i;
	  for (k2=0; k2<k1; ++k2)
	    ii += nnp[k2];
	  k2 = k1 + 1;
	  if (k2 == 4)
	    k2 = 0;
	  x[ii][j] = cx[k1] + (i*(cx[k2]-cx[k1])/nnp[k1]);
	  y[ii][j] = cy[k1] + (i*(cy[k2]-cy[k1])/nnp[k1]);
	}
    }
  
  /******************************************************/
  /* now do radius. Note that the nodes on the edge j=0
     are identical with those on the edge j=nr-1 and that we have
     already calculated the co-ords of points (0,0), (0,nr-1),
     (np-1,0) and (np-1,nr-1) */
  fprintf(log_f,"\nCalculating node positions along the radius\n");
  /* now set up node co-ords along the radius */
  /* don't bother resetting end points (see comment above) */
  s = 1/((float)nr-1);
  sum = s;
  fprintf(log_f,"\tsetting node co-ordinates along the radius\n");
  for (j=1; j<nr-1; ++j)
    {
      x[0][j] = x[0][0] + sum*(cx[0]-x[0][0]);
      y[0][j] = y[0][0] + sum*(cy[0]-y[0][0]);
      x[nr-1][j] = x[0][j];
      y[nr-1][j] = y[0][j];
      sum += s;
    }
  
  /*******************************************************/
  /* Now use eliptic smoother to generate interior nodes */
  /* first set all interior node co-ords to simple value */
  for (i=1; i<np; ++i)
    {
      sum = s;
      for (j=1; j<nr-1; ++j)
	{
	  x[i][j] = x[i][0] + sum*(x[i][nr-1]-x[i][0]);
	  y[i][j] = y[i][0] + sum*(y[i][nr-1]-y[i][0]);
	  sum += s;
	}
    }
  
  fprintf(log_f,"\nSmoothing interior node co-ordinates\n");
  /* call elliptic smoother */
  /* first smooth each quadrant individualy */
  esmooth(1,                     nnp[0],20);
  esmooth(nnp[0]+1,              nnp[0]+nnp[1],20);
  esmooth(nnp[0]+nnp[1]+1,       nnp[0]+nnp[1]+nnp[2],20);
  esmooth(nnp[0]+nnp[1]+nnp[2]+1,nnp[0]+nnp[1]+nnp[2]+nnp[3],20);
  /* now smooth whole grid a little bit */
  esmooth(0,np,5);
  /****************************************************************/
  
  /* so now we have all of the points. Next we create the fe grid */
  /* write the header to the AVS file */
  fprintf(avs,"%d\t%d\t1\t0\t0\n",np*nr,2*np*(nr-1));
  

  /* write the nodes out to the AVS file */
  for (j=0; j<nr; ++j)
    for (i=0; i<np; ++i)
      fprintf(avs,"%d\t%.8e\t%.8e\t0\n",(j*np)+i+1,x[i][j],y[i][j]);
  
  /* no of elements will be twice number of square cells */
  /* nelem = 2*(np-1)*(nr-1) */
  for (j=0; j<nr-1; ++j)
    for (i=0; i<np; ++i)
      {
	/* compute number of each node in this square */
	/* (i,j) is n1, then go anti-clockwise */
	if (i == np-1)
	  {
	    n1 = j*np + i;
	    n2 = j*np;
	    n3 = (j+1)*np;
	    n4 = (j+1)*np + i;
	  }
	else
	  {
	    n1 = j*np + i;
	    n2 = j*np + i+1;
	    n3 = (j+1)*np + i+1;
	    n4 = (j+1)*np + i;
	  }
	/* First tri in square - compute element number */
	nelem = 2*(i+j*np);
	elem[nelem][1] = n1;
	elem[nelem][2] = n3;
	elem[nelem][3] = n4;
	fprintf(avs,"%d %d tri %d %d %d\n",nelem+1,mat,n1+1,n3+1,n4+1);
	/* Second tri in square - compute element number */
	nelem = 2*(i+j*np)+1;
	elem[nelem][1] = n1;
	elem[nelem][2] = n2;
	elem[nelem][3] = n3;
	fprintf(avs,"%d %d tri %d %d %d\n",nelem+1,mat,n1+1,n2+1,n3+1);
      }

  /* DEBUGGER - print data values (node number) at nodes */
  /* IMPORTANT - REMEMBER TO CHANGE HEADER LINE TO NO DATA */
  fprintf(avs,"1 1\n");
  fprintf(avs,"ferret, ins\n");
  for (j=0; j<nr; ++j)
    for (i=0; i<np; ++i)
      fprintf(avs,"%d %d\n",(j*np)+i,(j*np)+i);
  
  /****************************************************************/
  /* write out grid file */
  fprintf(log_f,"\nWriting grid file ...\n");
  fprintf(log_f,"\twriting headers\n");
  fprintf(grd,"Gibbon\n");
  fprintf(grd,"%d %d 1 %d 0 0 %d 1 4\n",np*nr,2*np*(nr-1),np+nnp[3]+1,maxit);
  fprintf(grd,"%.3f 0.0\n",cond);
  for (i=0; i<12; ++i)
    fprintf(grd,"0.0 0.0\n");
  fprintf(grd,"1\n1.0 1.0 1.0\n");
  fprintf(grd,"1 %.3f 1.0\n",rhoc);
  /* write out x-co-ords of nodes - go in node number order */
  fprintf(log_f,"\twriting x co-ordinates\n");
  k1 = -1;
  for (j=0; j<nr; ++j)
    for (i=0; i<np; ++i)
      {
	fprintf(grd,"%.3e\t",x[i][j]);
	if (++k1%8 == 7)
	  fprintf(grd,"\n");
      }
  /* write out y-co-ords of nodes - go in node number order */
  fprintf(log_f,"\twriting y co-ordinates\n");
  for (j=0; j<nr; ++j)
    for (i=0; i<np; ++i)
      {
	fprintf(grd,"%.3e\t",y[i][j]);
	if (++k1%8 == 7)
	  fprintf(grd,"\n");
      }
  if (k1%8 != 7)
    fprintf(grd,"\n");

  /* next write out elements */
  fprintf(log_f,"\twriting elements\n");
  nelem = 2*np*(nr-1);
  for (i=0; i<nelem; ++i)
    fprintf(grd,"%3d %d %3d %3d %3d\n",i+1,mat,elem[i][1]+1,elem[i][2]+1,elem[i][3]+1);

  /* write out boundary elements */
  fprintf(log_f,"\twriting boundary nodes\n");
  /* write electrode first */
  k1 = -1;
  for (i=0; i<np; ++i)
    {
      fprintf(grd,"%3d ",(nr-1)*np+i+1);
      if (++k1%10 == 9)
	fprintf(grd,"\n");
    }
  /* now anodic nodes */
  ++k1;
  fprintf(grd,"1   ");
  for (i=nnp[0]+nnp[1]+nnp[2]; i<np; ++i)
    {
      fprintf(grd,"%3d ",i+1);
      if (++k1%10 == 9)
	fprintf(grd,"\n");
    }
  if (k1%10 != 9)
    fprintf(grd,"\n");

  /* now write out potential on electrode */
  fprintf(log_f,"\twriting potential at boundary nodes\n");
  k1 = -1;
  for (i=0; i<np; ++i)
    {
      fprintf(grd,"%6.3f ",volt);
      if (++k1%6 == 5)
	fprintf(grd,"\n");
    }
  /* and potential on the anode */
  for (i=0; i<nnp[3]+1; ++i)
    {
      fprintf(grd,"0.0 ");
      if (++k1%6 == 5)
	fprintf(grd,"\n");
    }
  if (k1%6 != 5)
    fprintf(grd,"\n");

  /* now write out the electrochemical data */
  fprintf(log_f,"\twriting electrochemical data\n");
  fprintf(grd,"%d\n",np+nnp[3]);
  fprintf(grd,"%.3f  5.3   -3.12\n3712.0 -255.0 -467.0\n",volt);
  fprintf(grd,"%3d %3d %3d %3d %3d\n",nnp[3],nnp[2],nnp[1],nnp[0],nnp[3]);
  fprintf(grd,"%d\n",nr*np/2 + nnp[0]+nnp[1]+nnp[2]+nnp[3]/2);

  /* now write out element numbers and nodes on boundaries */
  fprintf(log_f,"\twriting special boundary information\n");
  /* first do the electrode */
  /* do case i=np-1 first, then loop from np-2 down to i=0 */
  fprintf(grd,"%3d %3d %3d\n",2*(np-1+(nr-2)*np)+1,(nr-1)*np+1,(nr-1)*np+np-1+1);
  for (i=np-2; i>=0; --i)
    fprintf(grd,"%3d %3d %3d\n",2*(i+(nr-2)*np)+1,(nr-1)*np+i+2,(nr-1)*np+i+1);
  /* now do the anode */
  fprintf(grd,"%3d %3d %3d\n",2*(np-1)+2,1,np);
  for (i=np-2; i>=np-nnp[3]; --i)
    fprintf(grd,"%3d %3d %3d\n",2*i+2,i+2,i+1);

  /****************************************************************/
  fprintf(log_f,"\nFerret execution ended normally\n");
  exit(0);
}

/* maximum function */
float volemax(float x, float y)
{
  if (x>y)
    return x;
  else
    return y;
}

/* error in command line arguments function */
void cmderr(char c)
{
  printf("Error in command line arguments\n");
  printf("\teither : file name missing after -%c\n",c);
  printf("\t    or : no space between -%c and file name\n",c);
  exit(1);
}

/* elliptic smoother. Always sweeps across all interior nodes in
   the j-direction. The nodes to be swept in the i-direction are
   passed as parameters. */
void esmooth(int istart, int iend, int esits)
{
  int i, j, its, k1, k2;
  float xxi, yxi, xxixi, yxixi, xeta, yeta, xxieta, yxieta, xeta2,
  yeta2, g11, g22, g12, g11p22, xtmp, ytmp, dx, dy;

  fprintf(log_f,"\tsmoothing interior nodes %d to %d\n",istart,iend);
  fprintf(log_f,"\titer\tdx\t\tdy\n");
  for (its=0; its<esits; ++its)
    {
      dx = 0;
      dy = 0;
      for (j=1; j<nr-1; ++j)
	for (i=istart; i<iend; ++i)
	  {
	    /* since we are smoothing nodes on the i=0 and i=np-1
	       boundaries, we need to check for these cases as i-1
	       or i+1 would take us out of the array */
	    if (i-1 < 0)
	      k1 = np-1;
	    else
	      k1 = i-1;
	    if (i+1 == np)
	      k2 = 0;
	    else
	      k2 = i+1;
	    /* so now we use k1 for i-1 and k2 for i+1 */
	    /* the following bit of code is taken from GRIDGEN (c) CFDS, */
	    /* so the weird variable names are not my fault !  */
	    xxi   = 0.5*(x[k1][j]  - x[k2][j]);
	    yxi   = 0.5*(y[k1][j]  - y[k2][j]);
	    
	    xxixi =      x[k2][j]  + x[k1][j];
	    yxixi =      y[k2][j]  + y[k1][j];
	    
	    xeta  = 0.5*(x[i][j+1] - x[i][j-1]);
	    yeta  = 0.5*(y[i][j+1] - y[i][j-1]);
	    
	    xxieta = 0.25*(x[k1][j+1] - x[k1][j-1]
			   - x[k2][j+1] + x[k2][j-1]);
	    yxieta = 0.25*(y[k1][j+1] - y[k1][j-1]
			   - y[k2][j+1] + y[k2][j-1]);
	    
	    xeta2 = x[i][j+1] + x[i][j-1];
	    yeta2 = y[i][j+1] + y[i][j-1];
	    
	    g11 = xxi*xxi + yxi*yxi;
	    g22 = xeta*xeta + yeta*yeta;
	    g12 = xxi*xeta + yxi*yeta;
	    g11p22 = g11 + g22;
	    if (g11p22 == 0)
	      g11p22 = 1e-7;
	    
	    xtmp = 0.5*(g22*xxixi + g11*xeta2 - 2*g12*xxieta)/g11p22;
	    ytmp = 0.5*(g22*yxixi + g11*yeta2 - 2*g12*yxieta)/g11p22;
	    xtmp = asmth*xtmp + (1-asmth)*x[i][j];
	    ytmp = asmth*ytmp + (1-asmth)*y[i][j];
	    
	    dx = volemax(dx,fabs(xtmp-x[i][j]));
	    dy = volemax(dy,fabs(ytmp-y[i][j]));
	    
	    x[i][j] = xtmp;
	    y[i][j] = ytmp;
	  }
      fprintf(log_f,"\t%d\t%.6f\t%.6f\n",its,dx,dy);
    }
}
