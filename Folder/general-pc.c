/* The original version of this code was written by
   Prof. Istvan Szalai (Eötvös University, Budapest, Hungary). */

#include <stdio.h>
#include <math.h>

/* Header files with a description of contents used */

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_band.h>
#include <sundials/sundials_band.h>
#include <cvode/cvode_spgmr.h>
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <cvode/cvode_bandpre.h>

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NC   4                /* number of components */
#define NX   800              /* number of space points */
#define NEQ  3204             /* number of equations (NC*(NX+1) */
#define BAND 4                /* banded Jacobian */


/* Functions Called by the Solver */

static int f_osfr(realtype t, N_Vector y_osfr, N_Vector y_osfrdot, void *user_data);

/* Private functions to output results */

static void PrintOutput(realtype t, N_Vector y_osfr, int scan_par);
static void PrintOutput2(realtype t, N_Vector y_osfr, int scan_par);

static void PrintFinal(N_Vector y_osfr);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);

int NX1=NX+1;
int NEQ1=NEQ+1;
realtype par[50];  
realtype num_par[10]; 

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  realtype reltol, abstol, t, tout, delay;
  N_Vector y_osfr;
  void *cvode_mem;
  int flag, flagr, iout;
  realtype t0,delta,tend;
  int i,j,N;

  realtype cp0,dcp;
  int sweep, scan_par;

  y_osfr = NULL;
  abstol = reltol = 0;
  cvode_mem = NULL;

  /* Create serial vector of length NEQ for I.C. and abstol */
  y_osfr = N_VNew_Serial(NEQ);
  if (check_flag((void *)y_osfr, "N_VNew_Serial", 0)) return(1);

  /* Initialize y */
  input_osfr (y_osfr);  
  input_par (par);
  input_numpar (num_par);
  t0=num_par[0];
  tend=num_par[1];
  delta=num_par[2];
  delay=num_par[3];
  reltol=num_par[4];
  abstol=num_par[5];
  N= (int) num_par[6];


  t=t0;
  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  
  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f_osfr, t0, y_osfr);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and absolute tolerances */
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

  /* Call CVBand to specify the CVBAND band linear solver
   * flag = CVBand(cvode_mem, NEQ, BAND, BAND);
   * if(check_flag(&flag, "CVBand", 1)) return(1); */

  flag = CVSpgmr(cvode_mem, PREC_LEFT, 0);
  if(check_flag(&flag, "CVSpgmr", 1)) return(1);

  flag = CVBandPrecInit(cvode_mem, NEQ, BAND, BAND);
  if(check_flag(&flag, "CVBandPrecInit", 0)) return(1);


  /* In loop, call CVode, print results, and test for error. */

  tout = t0+delta;
  if (t>=delay) {
  PrintOutput(t, y_osfr,scan_par);
  PrintOutput2(t, y_osfr,scan_par);
  }
  j=1;
  while(1) {
    flag = CVode(cvode_mem, tout, y_osfr, &t, CV_NORMAL);
    if (j % N == 0 && t>=delay) {
    PrintOutput(t, y_osfr,scan_par);
    PrintOutput2(t, y_osfr,scan_par);
    }
    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS) {
      tout += delta;
    }
    if (tout > tend) break;
  ++j;
  }


  PrintFinal(y_osfr);
  /* Print some final statistics */
  /* PrintFinalStats(cvode_mem); */

  /* Free y and abstol vectors */
  N_VDestroy_Serial(y_osfr);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y).
 */

static int f_osfr(realtype t, N_Vector y_osfr, N_Vector y_osfrdot, void *user_data)
{

/*
A + B --> C
   4C --> D

A: freely diffusing strand that causes plasticization
B: low mobility strand in the PC
C: low mobility complex of AB in the glassy state
D: low mobility complex of AB in the plasticized state */


/* Chemical components */

  realtype A[NX1],B[NX1],C[NX1],D[NX1];
  realtype dA[NX1],dB[NX1],dC[NX1],dD[NX1];
  realtype v1[NX1],v2[NX1]; 
  int i;


/* Chemical parameters	D0, D1: for component A in the glassy and plasticized regimes
			DS: for the non-diffusing components
			Q: threshold for glassy/plasticized transition */

  realtype k1=par[1];
  realtype k2=par[2];
  realtype Ai=par[3];
  realtype Bi=par[4];
  realtype D0=par[5];
  realtype D1=par[6];
  realtype DS=par[7];
  realtype Q=par[8];
  realtype x=num_par[7];
  realtype d_x=1/(x*x);


/* CSTR: X[0], 0D RD-domain X[1] */

  for (i=0;i<NX1;i++) {
  A[i] = Ith(y_osfr,NC*i+1);
  B[i] = Ith(y_osfr,NC*i+2);
  C[i] = Ith(y_osfr,NC*i+3);
  D[i] = Ith(y_osfr,NC*i+4);
  v1[i]=k1*A[i]*B[i];
  v2[i]=k2*C[i];
  }


/* 0D RD-domain ODEs with Dirichlet and Neumann BC */

/* Source */
  dA[0] = Ith(y_osfrdot,1) = Ai-A[0];
  dB[0] = Ith(y_osfrdot,2) = 0-B[0];
  dC[0] = Ith(y_osfrdot,3) = 0-C[0];
  dD[0] = Ith(y_osfrdot,4) = 0-D[0];

  for (i=1;i<NX;i++) {

/* Glassy */
  if (D[i]<Q) {
  if (D[i+1]<Q) {
  dA[i] = Ith(y_osfrdot,NC*i+1) = -v1[i]	+D0*(A[i+1]+A[i-1]-2*A[i])*d_x;
  dB[i] = Ith(y_osfrdot,NC*i+2) = -v1[i]	+DS*(B[i+1]+B[i-1]-2*B[i])*d_x;
  dC[i] = Ith(y_osfrdot,NC*i+3) = v1[i]-4*v2[i]	+DS*(C[i+1]+C[i-1]-2*C[i])*d_x;
  dD[i] = Ith(y_osfrdot,NC*i+4) = v2[i]		+DS*(D[i+1]+D[i-1]-2*D[i])*d_x;
  }
  }

/* Plasticized */
  if (D[i]>=Q) {
  if (D[i+1]>=Q) {
  dA[i] = Ith(y_osfrdot,NC*i+1) = -v1[i]	+D1*(A[i+1]+A[i-1]-2*A[i])*d_x;
  dB[i] = Ith(y_osfrdot,NC*i+2) = -v1[i]	+DS*(B[i+1]+B[i-1]-2*B[i])*d_x;
  dC[i] = Ith(y_osfrdot,NC*i+3) = v1[i]-4*v2[i]	+DS*(C[i+1]+C[i-1]-2*C[i])*d_x;
  dD[i] = Ith(y_osfrdot,NC*i+4) = v2[i]		+DS*(D[i+1]+D[i-1]-2*D[i])*d_x;
  }
  }

/* Phase-boundary */
  if (D[i]>=Q) {
  if (D[i+1]<Q) {
  dA[i] = Ith(y_osfrdot,NC*i+1) = -v1[i]	+D0*(A[i+1]-A[i])*d_x	-D1*(A[i]-A[i-1])*d_x;
  dB[i] = Ith(y_osfrdot,NC*i+2) = -v1[i]	+DS*(B[i+1]-B[i])*d_x	-DS*(B[i]-B[i-1])*d_x;
  dC[i] = Ith(y_osfrdot,NC*i+3) = v1[i]-4*v2[i]	+DS*(C[i+1]-C[i])*d_x	-DS*(C[i]-C[i-1])*d_x;
  dD[i] = Ith(y_osfrdot,NC*i+4) = v2[i]		+DS*(D[i+1]-D[i])*d_x	-DS*(D[i]-D[i-1])*d_x;
  }
  }

  }

/* Wall */
  dA[NX] = Ith(y_osfrdot,NC*NX+1) = -v1[NX]		-D1*(A[NX]-A[NX-1])*d_x;
  dB[NX] = Ith(y_osfrdot,NC*NX+2) = -v1[NX]		-DS*(B[NX]-B[NX-1])*d_x;
  dC[NX] = Ith(y_osfrdot,NC*NX+3) = v1[NX]-4*v2[NX]	-DS*(C[NX]-C[NX-1])*d_x;
  dD[NX] = Ith(y_osfrdot,NC*NX+4) = v2[NX]		-DS*(D[NX]-D[NX-1])*d_x;


  return(0);
}


/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(realtype t, N_Vector y_osfr, int scan_par)
{
  int i,j;
  realtype x=num_par[7];

/*
    printf("%0.4le", t);
    for (i=1;i<limit_neq;i++) {
      printf("%14.4le", Ith(y_osfr,i));
      }
    printf("\n");
*/

  for (i=0;i<NX1;i++) {
    printf("%14.4le", t);
    printf("%14.3le", i*x);
    for (j=(1+i*NC);j<=(NC*(1+i));j++) {
      printf("%14.4le", Ith(y_osfr,j));
      }
    printf("\n");
    }
  printf("\n");

  return;
}

static void PrintOutput2(realtype t, N_Vector y_osfr, int scan_par)
{
  FILE *out;
  int i,j;
  realtype x=num_par[7];
  out = fopen("out1", "a");

/*  printf("%0.4le", t);
    for (i=1;i<limit_neq;i++) {
      printf("%14.4le", Ith(y_osfr,i));
      }
    printf("\n");

*/


  i=133;
  fprintf(out,"%14.4le", t);
  fprintf(out,"%14.3le", i*x);
  for (j=(1+i*NC);j<=(NC*(1+i));j++) {
     fprintf(out,"%14.4le", Ith(y_osfr,j));
     }
  fprintf(out,"\n");
  fclose(out);


  return;
}


static void PrintFinal(N_Vector y_osfr)
{
  FILE *out;
  int i,j;
  realtype y=num_par[7];
  out = fopen("final_osfr", "w");

  for (i=1;i<=(NX+1);i++) {
  for (j=1;j<=NC;j++){
    fprintf(out,"%14.4le", Ith(y_osfr,(i-1)*NC+j));
    }
  fprintf(out,"\n");
  }

  fclose(out);
  return;
}


/* Input */
int input_osfr (N_Vector y_osfr)
{
  FILE *inf;
  realtype q;
  int i;
  
  if ((inf = fopen("input_osfr", "r")) == NULL)
	{printf("The input_osfr file cannot be opened\n");
	return 1;}
  else {for (i=1;i<NEQ1;i++) {
	fscanf(inf, "%le", &q);
	Ith(y_osfr,i)=q;
        }}
	fclose(inf);
  return 0;
}


int input_numpar (double num_par[10])
{
  FILE *inf;
  realtype q;
  int i, npar;
  char par_name[10];
  
  if ((inf = fopen("input_numpar", "r")) == NULL)
	{printf("The input_numpar file cannot be opened\n");
	return 1;}
  else {fscanf(inf, "%d\n", &npar);
	for (i=0;i<npar;i++) {
	fscanf(inf, "%s%le\n", &par_name, &q);
	num_par[i]=q;
        printf("%s = %0.4le\n", par_name, num_par[i]);
        }}
	fclose(inf);
  return 0;
}


int input_par (double par[50])
{
  FILE *inf;
  realtype q;
  int i, npar;
  char par_name[10];
  
  if ((inf = fopen("input_par", "r")) == NULL)
	{printf("The input_par file cannot be opened\n");
	return 1;}
  else {fscanf(inf, "%d\n", &npar);
	for (i=1;i<npar+1;i++) {
	fscanf(inf, "%s%le\n", &par_name, &q);
	par[i]=q;
        printf("%s = %0.4le\n", par_name, par[i]);
        }}
	fclose(inf);
  return 0;
}


/*
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);

}
