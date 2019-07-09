
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "hessianApprox.hpp"
#include "braid.h"

#define M_PI 3.14159265358979323846

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    sstart;
   double    sstop;
   int       ntime;     /* Number of uniform time points */
   int       ndisc;     /* Number of switching events */
   int       rank;
   // Optim
   double*   design;    /* Stores the switching times */
   double*   gradient;  /* dJ/dDesign */
   double    path_penalty_param;   /* Param for state  penalty */
   // Model
   double    exciter_param;  /* G0 parameter in exciter model */
   double    Vmax;           /* Max limit of exciter */
   double    Vmin;           /* Min limit of exciter */
   double    V0;             /* Initial condition */
   double    tau;            /* Time constant */

} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double volt;      // Exciter voltage
} my_Vector;


/* analytic solution for inputVoltage 2*sin x + 2 */
double exactSolution(braid_App app,
                     double t)
{

   double ratio = app->exciter_param / app->tau;

   double exact = exp(-ratio*t) * app->V0;
   // InputVoltage = 2*sin(t-pi/2) + 2
   // exact += 2. * ratio * app->exciter_param / (pow(ratio,3) + ratio) * ( - pow(ratio,2) * cos(t) + pow(ratio,2) - ratio * sin(t) + 1. - exp(-ratio*t));
   // InputVoltage = sin(t) 
   exact += ratio * app->exciter_param / (pow(ratio,2) + 1) * ( ratio * sin(t) - cos(t) + exp(-ratio*t));


   return exact;
}

/* Dynamic input voltage */
double inputVoltage(double t)
{
   // Put some oscillation here. TODO: ASK PHILIPP!
   // double volt = 2.*sin(t) + 2.;
   // double volt = 2.*sin(t - M_PI/2.0) + 2.;
   double volt = sin(t);

   return volt;
}

/* a(t) = sum_k (-1)^k*design(k)*indicatorfunction_[k,k+1)(t) */
double getA(double* design,
            int     ndesign,
            double  t)
{
   double a = .0;

   /* Find interval */
   for (int i=0; i < ndesign; i++)
   {
      if (i <= t &&  t < i+1)
      {
         a += pow(-1.0,(double) i) * design[i];
      }
   }

   return a;
}

/* Derivative of getA times a_bar */
void getA_diff(double* gradient,
               int     ndesign,
               double  t,
               double  a_bar)
{
   for (int i=0; i < ndesign; i++)
   {
      if (i <= t &&  t < i+1)
      {
         double tmp = pow(-1.0, (double) i) * a_bar;
         gradient[i] += tmp;
      }
   }
}

/* Get original time t(s) */
double getOriginalTime(double *design, 
                    int     ndesign, 
                    double  s)
{
   double ts = 0.0;

   int i=0;
   for (i = 0; i < ndesign-1; i++)
   {
      /* Find interval such that s \in [k,k+1]*/
      if (i <= s && s < i+1)
      {
         break;
      }
      else
      {
         ts += design[i];
      }
   }

   /* Add remaining time (k_0,s)*/
   ts += (s - i) * design[i];

   return ts;
}

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double sstart;             /* current time */
   double sstop;              /* evolve to this time*/
   int    istart;             /* time point index value corresponding to tstop on (global) fine grid */
   int    level, iter;
   double ds;
   double control;

   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetIter(status, &iter);
   braid_StepStatusGetTstartTstop(status, &sstart, &sstop);
   braid_StepStatusGetTIndex(status, &istart);

   /* Find the number of switches in between sstop and sstart */
   int nswitches = (int) ceil(sstop) - (int) floor(sstart) - 1;


   double u_curr = u->volt;
   // if (nswitches > 0)
   // {
   //    // if (nswitches > 1) printf("%f -> %f, %d switches: \n", sstart, sstop, nswitches);
   //    // printf("%f -> %f, %d switches: \n", sstart, sstop, nswitches);

   //    /* Step from sstart to first switch */
   //    ds = ceil(sstart) - sstart;
   //    control = getA(app->design, app->ndisc+1, sstart);
   //    u_curr = 1./(1. - control * ds) * u_curr;
   //    // printf("first to %f: %f %f, ", ceil(sstart), ds, control);

   //    /* Step through all intermediate switches */
   //    for (int iswitch = 0; iswitch < nswitches-1; iswitch++)
   //    {
   //       ds = 1.0; // since switches happen at integer times 
   //       control = getA(app->design, app->ndisc+1, ceil(sstart) + (double) iswitch);
   //       u_curr = 1./(1. - control * ds) * u_curr;
   //       // printf("%d %f %f, ", iswitch, ds, control);
   //    }

   //    /* Step from last switch to sstop */
   //    ds = sstop - floor(sstop);
   //    control = getA(app->design, app->ndisc+1, floor(sstop));
   //    u_curr = 1./(1. - control * ds) * u_curr;
   //    // printf("last to %f: %f %f\n, ", sstop, ds, control);
   // }
   // else
   // {
   //    /* Step from sstart to sstop */
   //    control = getA(app->design, app->ndisc+1, sstart);
   //    ds = sstop - sstart;
   //    u_curr = 1./(1. - control * ds) * u_curr;
   // }

   // /* Set u */
   // u->volt = u_curr;


   /* Step from sstart to sstop */
   /* NO LIMITS */
   ds = sstop - sstart;
   double ratio  = ds / app->tau;
   double unext = 1./(1. + ratio) * (u_curr + ratio * app->exciter_param * inputVoltage(sstop, app->sstop) );
   // printf("step %f->%f, ds=%f, u: %f -> %f\n", sstart, sstop, ds, u_curr, u->volt);

   /* Enforce the limits: If limits are exceeded, set dV/dt = 0 */
   if (unext > app->Vmax)
   {
      // unext = u_curr;
      unext = app->Vmax;
   } 
   else if (unext < app->Vmin)
   {
      // unext = u_curr;
      unext = app->Vmin;
   }

   /* Set the output voltage */
   u->volt = unext;   

   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}


int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   
   if (t == 0.0) /* Initial condition */
   {
      (u->volt) = app->V0;
   }
   else /* All other time points set to arbitrary value */
   {
      (u->volt) = app->V0;
   }
   *u_ptr = u;

   return 0;
}

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;

   v = (my_Vector *) malloc(sizeof(my_Vector));
   (v->volt) = (u->volt);
   *v_ptr = v;

   return 0;
}

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u);

   return 0;
}

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{
   (y->volt) = alpha*(x->volt) + beta*(y->volt);
   
   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot = (u->volt)*(u->volt);
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        index;
   int        iter;
   int        level;
   double     s, ts;
   char       filename[255];
   FILE      *file;
   
   braid_AccessStatusGetLevel(astatus, &level);
   braid_AccessStatusGetIter(astatus, &iter);
   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetT(astatus, &s);
   // sprintf(filename, "%s.%04d.%03d", "ex-01-expanded.out", index, app->rank);
   // file = fopen(filename, "w");
   // fprintf(file, "%1.4f  %.14e\n", t, (u->volt));
   // fflush(file);
   // fclose(file);

   /* Get original time */
   ts = getOriginalTime(app->design, app->ndisc+1, s);

   /* Append all into one file. */
   if (level == 0)
   {
      // sprintf(filename, "%s.iter%03d.%03d", "exciter-model-optim.out", iter, app->rank);
      sprintf(filename, "%s.%03d", "exciter-model-optim.out", app->rank);
      file = fopen(filename, "a");
      fprintf(file, "%04d %1.4f %1.4f %.14e %1.14e\n", index, s, ts, (u->volt), exactSolution(app, s));
      fflush(file);
      fclose(file);

   }

   // if (index == app->ntime) printf("Last: %d  %4f  %1.14e\n", index, t, u->volt);

   return 0;
}


int
my_BufSize(braid_App          app,
           int                *size_ptr,
           braid_BufferStatus bstatus)
{
   *size_ptr = sizeof(double);
   return 0;
}

int
my_BufPack(braid_App          app,
           braid_Vector       u,
           void               *buffer,
           braid_BufferStatus bstatus)
{
   double *dbuffer = (double*)buffer;

   dbuffer[0] = (u->volt);
   braid_BufferStatusSetSize( bstatus, sizeof(double) );

   return 0;
}

int
my_BufUnpack(braid_App          app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus bstatus)
{
   double    *dbuffer = (double*)buffer;
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->volt) = dbuffer[0];
   *u_ptr = u;

   return 0;
}

int
my_ObjectiveT(braid_App app,
              braid_Vector u,
              braid_ObjectiveStatus ostatus,
              double *objectiveT_ptr)
{
   double pathconstraint = 0.0;
   double objT = 0.0;

   /* Get current time */
   double t;
   braid_ObjectiveStatusGetT(ostatus, &t);

   /* --- y(2k)=1 & y(2k+1) = 2 -- */

   // /* if t is an integer and not first or last time step */
   // if ( (t == floor(t)) && (t > 0) && (t < app->ndisc+1)) 
   // {
   //    if ( (int) t % 2 == 0 )
   //    {
   //       // if t is even: y = 1 !
   //       pathconstraint = 0.5 * pow(u->volt-1,2);
   //    }
   //    else
   //    {
   //       // if t is odd :  y = 2 !
   //       pathconstraint = 0.5 * pow(u->volt-2,2);
   //    }
   // }
   // objT += app->path_penalty_param * pathconstraint;


   /* set return value */
   *objectiveT_ptr =  objT;

   return 0;
}

/* Transposed partial derivatives of objectiveT times F_bar */
int
my_ObjectiveT_diff(braid_App            app,
                  braid_Vector          u,
                  braid_Vector          u_bar,
                  braid_Real            F_bar,
                  braid_ObjectiveStatus ostatus)
{
   double fd_diff             = 0.0;
   double oneoverdt           = 0.0;
   double ddpathconstraint    = 0.0;
   double diff, fd;

   /* Get current time */
   double t;
   braid_ObjectiveStatusGetT(ostatus, &t);

   /* --- y(2k)=1 & y(2k+1) = 2 -- */

   // /* if t is an integer and not first or last time step */
   // if ( (t == floor(t)) && (t > 0) && (t < app->ndisc+1)) 
   // {
   //    if ( (int) t % 2 == 0 )
   //    {
   //       // if t is even: y = 1 !
   //       ddpathconstraint += u->volt - 1.0;
   //    }
   //    else
   //    {
   //       // if t is odd :  y = 2 !
   //       ddpathconstraint += u->volt - 2.0;
   //    }
   // }
   // u_bar->volt = app->path_penalty_param * ddpathconstraint * F_bar;



   return 0;
}

int
my_Step_diff(braid_App           app,
             braid_Vector        ustop,
             braid_Vector        u,
             braid_Vector        ustop_bar,
             braid_Vector        u_bar,
             braid_StepStatus    status)
{
   double ddu;      /* Derivative wrt u */
   double ddc;      /* Derivative wrt control */
   double ddesign;  /* Derivative wrt design */
   double sstart;             /* current time */
   double sstop;              /* evolve to this time*/
   double dsigm;
   int istart;
   double ds;
   double control;
   double u_curr;

   braid_StepStatusGetTstartTstop(status, &sstart, &sstop);
   braid_StepStatusGetTIndex(status, &istart);
   double deltat = sstop - sstart;

   /* Find the number of switches in between sstop and sstart */
   int nswitches = (int) ceil(sstop) - (int) floor(sstart) - 1;

   u_curr = u->volt;
   ddu    = u_bar->volt;

   // if (nswitches > 0)
   // {

   //    /* Move forward, storing u */

   //    double* utmp = new double[nswitches];
   //    /* First step */
   //    ds = ceil(sstart) - sstart;
   //    control = getA(app->design, app->ndisc+1, sstart);
   //    u_curr = 1./(1. - control * ds) * u_curr;
   //    utmp[0] = u_curr;
   //    /* Intermediate steps */
   //    for (int iswitch = 0; iswitch < nswitches-1; iswitch++)
   //    {
   //       ds = 1.0; 
   //       control = getA(app->design, app->ndisc+1, ceil(sstart) + (double) iswitch);
   //       u_curr = 1./(1. - control * ds) * u_curr;
   //       utmp[iswitch+1] = u_curr;
   //    }
   //    /* last step */
   //    ds = sstop - floor(sstop);
   //    control = getA(app->design, app->ndisc+1, ceil(sstart) + nswitches-1);
   //    u_curr = 1./(1. - control * ds) * u_curr;


   //    /* Step backwards */

   //    /* Last step */
   //    ds = sstop - floor(sstop);
   //    control = getA(app->design, app->ndisc+1, floor(sstop));
   //    ddu = 1./(1. - control * ds) * u_bar->volt;
   //    ddc = (ds* (u_curr)) / pow(1. - ds*control,2) * (u_bar->volt);
   //    getA_diff(app->gradient, app->ndisc+1, floor(sstop), ddc); 
   //    /* Intermediate switches */
   //    for (int iswitch = nswitches - 2; iswitch > -1; iswitch--)
   //    {
   //       ds = 1.0; 
   //       control = getA(app->design, app->ndisc+1, ceil(sstart) + (double) iswitch);
   //       ddc = (ds* (utmp[iswitch])) / pow(1. - ds*control,2) * (ddu);
   //       ddu  = 1./(1. - control * ds) * ddu;
   //       getA_diff(app->gradient, app->ndisc+1, ceil(sstart) + (double) iswitch, ddc); 
   //    }
   //    /* first step */
   //    ds = ceil(sstart) - sstart;
   //    control = getA(app->design, app->ndisc+1, sstart);
   //    ddc = (ds* (utmp[0])) / pow(1. - ds*control,2) * (ddu);
   //    ddu = 1./(1. - control * ds) * ddu;
   //    getA_diff(app->gradient, app->ndisc+1, sstart, ddc); 

      // delete [] utmp;

   // }
   // else
   // {
   //    /* Step from sstart to sstop */
   //    control = getA(app->design, app->ndisc+1, sstart);
   //    ds = sstop - sstart;
   //    ddc = (ds* (u->volt)) / pow(1. - ds*control,2) * (ddu);
   //    ddu = 1./(1. - control * ds) * ddu;
   //    getA_diff(app->gradient, app->ndisc+1, sstart, ddc); 
   // }


   /* Update u_bar */
   u_bar->volt = ddu;              

   return 0;
}


/* Set the gradient to zero */
int 
my_ResetGradient(braid_App app)
{
   for (int i = 0; i < app->ndisc+1; i++)
   {
      app->gradient[i] = 0.0;
   }

   return 0;
}

/* Function to allow for the computation of the gradient */
int
gradient_allreduce(braid_App app)
{
   /* Copy gradient */
   double *mygradient = (double*) malloc((app->ndisc+1)*sizeof(double)); 
   for (int i = 0; i<app->ndisc+1; i++)
   {
      mygradient[i] = app->gradient[i];
   }

   /* Collect sensitivities from all processors and broadcast it */
   MPI_Allreduce(mygradient, app->gradient, app->ndisc+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   /* Clean up */
   free(mygradient);

   return 0;
}


void write_vector(char   *filename,
                  double * var, 
                  int      dimN)
{
   FILE *file;
   int i;

   /* open file */
   file = fopen(filename, "w");
   if (file == NULL)
   {
      printf("Can't open %s \n", filename);
      exit(1);
   }
   /* Write data */
   printf("Writing file %s\n", filename);
   for ( i = 0; i < dimN; i++)
   {
      fprintf(file, "%04d  %1.14e\n", i, var[i]);
   }

   /* close file */
   fclose(file);
}            


/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   braid_Core    core;
   MPI_Comm      comm;
   my_App       *app;
   double        objective;

   double        gnorm;
   int           optimiter;
   int           arg_index;
   int           rank, size;
   double        mygradnorm;
   FILE         *optimfile;

   double ls_objective;
   double ls_param    = 1e-4;
   double ls_maxiter  = 40;

   /* Setting defaults */
   // Braid
   int         max_levels = 5;
   int         nrelax     = 1;
   int         nrelax0    = -1;
   double      tol        = 1.0e-06;
   int         cfactor    = 2;
   int         print_level= 1;
   int         max_iter   = 100;
   int         fmg        = 0;
   int         storage    = -1;
   // Optim
   double      path_penalty_param   = 1;    /* Param for path constraint */
   int         maxoptimiter         = 100;  /* Maximum optimization iterations */
   double      stepsize_init        = 1.0;  /* Set initial Step size for design updates */
   double      gtol                 = 1e-6; /* Stopping criterion on the gradient norm */
   // Exciter model
   double      exciter_param   = 2.0;   /* G0 parameter in exciter model */
   double      Vmax            = 1.0;   /* Max limit of exciter */
   double      Vmin            = -1.0;   /* Min limit of exciter */
   double      V0              = 0.0;   /* Initial condition */
   double      tau             = .1;    /* Time constant */

   /* Default time domain */
   int ntime  = 400;          /* Number of time steps */
   int ndisc  = 3;            /* Number of discontinuities / switches */

   /* Initialize MPI */
   comm   = MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &size);

   /* Hack: Open and close output file. */
   char       filename[255];
   FILE      *file;
   sprintf(filename, "%s.%03d", "exciter-model-optim.out", rank);
   file = fopen(filename, "w");
   fclose(file);

   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( rank == 0 )
         {
            printf("\nExample 1: Solve a scalar ODE \n\n");
            printf("  -ntime <ntime>         : set num time points\n");
            printf("  -ndisc <ndisc>         : set num of switching points\n");
            printf("  -ml  <max_levels>      : set max levels\n");
            printf("  -nu  <nrelax>          : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>          : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>             : set stopping tolerance\n");
            printf("  -cf  <cfactor>         : set coarsening factor\n");
            printf("  -pl  <printlevel>      : set print level\n");
            printf("  -fmg                   : use FMG cycling\n");
            printf("  -storage <level>       : full storage on levels >= level\n");
            printf("  -mi  <max_iter>        : max braid iterations\n");
            printf("  -moi <int>             : max optimization iter\n");
            printf("  -dstep <double>        : step size for design updates\n");
            printf("  -gtol <double>         : optimization stopping tolerance\n");
            printf("  -ppp <double>          : path constraint penalty parameter \n");
         }
         exit(1);
      } 
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ndisc") == 0 )
      {
         arg_index++;
         ndisc = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu0") == 0 )
      {
         arg_index++;
         nrelax0 = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 ) 
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-gtol") == 0 ) 
      {
         arg_index++;
         gtol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-dstep") == 0 ) 
      {
         arg_index++;
         stepsize_init = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-pl") == 0 )
      {
         arg_index++;
         print_level = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         max_iter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-moi") == 0 )
      {
         arg_index++;
         maxoptimiter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
      }
      else if ( strcmp(argv[arg_index], "-storage") == 0 ){
         arg_index++;
         storage = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ppp") == 0 ) 
      {
         arg_index++;
         path_penalty_param = atof(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   /* Transformed time domain */
   double sstart = 0.0;
   // double sstop  = (double) (ndisc + 1);  
   double sstop  = 10;  
   sstop = 320;
   ntime = 6400;
   /* Sanity check: ntime / (ndisc + 1) must be integer for objective function evaluation */
   if (ntime % (ndisc + 1) != 0)
   {
      printf("\nError: Choose ntime, ndisc such that ntime / (ndisc+1) is integer.\n");
      printf("This is required for evaluating the objective function.\n\n");
      exit(1);
   }


   /* initialize optimization */
   double* design;
   double* gradient;
   double* design0;
   double* gradient0;
   double* ascentdir;
   design    = (double*) malloc((ndisc+1) * sizeof(double));
   gradient  = (double*) malloc((ndisc+1) * sizeof(double));
   design0   = (double*) malloc((ndisc+1) * sizeof(double));
   gradient0 = (double*) malloc((ndisc+1) * sizeof(double));
   ascentdir = (double*) malloc((ndisc+1) * sizeof(double));
   for (int i = 0; i < ndisc+1; i++)
   {
      design[i]    =  .5; // Initial segment length in original time t.
      gradient[i]  = 0.0;
      ascentdir[i] = 0.0;
   }
   HessianApprox* hessian = new L_BFGS(MPI_COMM_WORLD, ndisc+1, 50);

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm)     = comm;
   (app->sstart)   = sstart;
   (app->sstop)    = sstop;
   (app->ntime)    = ntime;
   (app->ndisc)    = ndisc;
   (app->rank)     = rank;
   // Optim
   (app->design)   = design;
   (app->gradient) = gradient;
   (app->path_penalty_param)   = path_penalty_param;
   // Model
   app->exciter_param = exciter_param;
   app->Vmax          = Vmax;
   app->Vmin          = Vmin;
   app->V0            = V0;  
   app->tau           = tau;

   /* initialize XBraid */
   braid_Init(comm, comm, sstart, sstop, ntime, app,
            my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
            my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

  /* Initialize adjoint-based gradient computation */
   braid_InitAdjoint( my_ObjectiveT, my_ObjectiveT_diff, my_Step_diff, my_ResetGradient, &core);


   /* Set options */
   braid_SetPrintLevel( core, print_level);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   if (fmg)
   {
      braid_SetFMG(core);
   }
   if (storage >= -2)
   {
      braid_SetStorage(core, storage);
   }

   /* Prepare output */
   sprintf(filename, "%s", "optim.out");
   optimfile = fopen(filename, "w");
   if (rank == 0) 
   {
      printf("#Iter Objective            rel.obj     ||grad||             rel.||grad||    stepsize\n"); 
      fprintf(optimfile, "#Iter Objective            rel.obj      ||grad||           rel.||grad||    stepsize\n"); 
   }


   /* Optimization iteration */
   double obj_init    = 1.0;
   double gnorm_init  = 1.0;
   double ls_stepsize = 0.002;  // HACK!
   double StartTime = MPI_Wtime();
   for (optimiter = 0; optimiter < maxoptimiter; optimiter++)

   {

      /* Run adjoint XBraid to compute objective function and gradient */
      braid_SetAccessLevel(core, 0);
      braid_SetObjectiveOnly(core, 0);
      braid_Drive(core);

      /* Get the objective function value */
      braid_GetObjective(core, &objective);
      if (optimiter == 0) obj_init = objective;

      /* Allreduce gradient and compute its norm */
      gradient_allreduce(app);
      gnorm = 0.0;
      for (int id = 0; id < ndisc+1; id++)
      {
         gnorm += pow(app->gradient[id], 2);
      }
      gnorm = sqrt(gnorm);
      if (optimiter == 0) gnorm_init = gnorm;

      /* Output */
      if (rank == 0) 
      {
         printf("%3d: %1.14e  %1.4e  %1.14e  %1.4e  %.8f\n", optimiter, objective, objective/obj_init, gnorm, gnorm/gnorm_init, ls_stepsize);
         fprintf(optimfile, "%3d: %1.14e  %1.4e  %1.14e  %1.4e  %.8f\n", optimiter, objective, objective/obj_init, gnorm, gnorm/gnorm_init, ls_stepsize);
      }

      /* Check optimization convergence */
      if (gnorm/gnorm_init < gtol)
      {
         break;
      }

      /* Store design and gradient */
      for (int idx = 0; idx < ndisc; idx++)
      {
         gradient0[idx] = app->gradient[idx];
         design0[idx]   = app->design[idx];
      }

      /* Update hessian approximation */
      hessian->updateMemory(optimiter, app->design, app->gradient);
      hessian->computeAscentDir(optimiter, app->gradient, ascentdir);

      /* Design update using simple steepest descent method */
      // printf("Init design update stepsize %f\n", ls_stepsize);
      for (int idx = 0; idx < ndisc; idx++)
      {
         app->design[idx] = design0[idx] - ls_stepsize * ascentdir[idx];
      }

      /* --- Backtracking linesearch -- */

      /* Compute wolfe condition */
      double wolfe = 0.0;
      for (int idx = 0; idx < ndisc; idx++)
      {
         wolfe += pow(app->gradient[idx], 2);
      } 
      
      /* Perform linesearch updates */
      for (int ls_iter = 0; ls_iter < ls_maxiter; ls_iter++)
      {
         /* Get objective */
         braid_SetObjectiveOnly(core, 1);
         braid_Drive(core);
         braid_GetObjective(core, &ls_objective);

         double test_obj = objective - ls_param * ls_stepsize * wolfe;
         // printf("ls %d: new step %1.14e < %1.14e ?\n", ls_iter, ls_objective, test_obj);
         if (ls_objective <= test_obj)
         {
            break; // line search success, use this design 
         }
         else
         {
            /* Test for line-search failure */
            if (ls_iter == ls_maxiter - 1)
            {
               if (rank == 0) printf("\n WARNING: LINESEARCH FAILED! \n");
               if (rank == 0) fprintf(optimfile,"# WARNING: LINESEARCH FAILED! \n");
            }

            /* Go back half of the step */
            ls_stepsize = 0.5 * ls_stepsize;
            // printf("LS design update stepsize %f\n", ls_stepsize);
            for (int idx = 0; idx < ndisc; idx++)
            {
               app->design[idx] = design0[idx] - ls_stepsize * ascentdir[idx];
            }
         }
         
      }

      /* Reset stepsize */
      ls_stepsize = stepsize_init;

   }

   /* Timing */
   double myUsedTime = MPI_Wtime() - StartTime;
   double UsedTime;
   MPI_Allreduce(&myUsedTime, &UsedTime, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   UsedTime = UsedTime / ((double) size);
   if (rank == 0) printf("\n Used time: %f sec\n", UsedTime);

   /* Print some statistics about the optimization run */
   if (rank == 0)
   {
      if (optimiter == maxoptimiter)
      {
         printf("\n Max. number of iterations reached! \n\n"); 
         fprintf(optimfile, "#Max. number of iterations reached! \n\n"); 
      }
      else
      {
         printf("\n");
         printf("  Optimization has converged.\n");
         printf("  Be happy and go home!      \n");
         printf("\n"); 
         printf("  Objective function value = %1.8e\n", objective);
         printf("  Gradient norm            = %1.8e\n", gnorm);
         printf("\n");
         printf("  optimization iterations  = %d\n", optimiter);
         printf("  max optim iterations     = %d\n", maxoptimiter);
         printf("  rel. gradient norm tolerance = %1.1e", gtol);
         printf("\n");
      }
   }

   /* Get final access */
   braid_SetAccessLevel(core, 1);
   braid_SetObjectiveOnly(core, 1);
   braid_Drive(core);

   /* Print XBraid statistics */
   braid_PrintStats(core);


   /* Finish braid */
   braid_Destroy(core);

   /* print */
   if (rank == 0)
   {
      sprintf(filename, "%s.%03d", "design.out", rank);
      write_vector(filename, app->design, ndisc+1);
      sprintf(filename, "%s.%03d", "gradient.out", rank);
      write_vector(filename, app->gradient, ndisc+1);
   }

   /* Close optimization output file */
   if (rank == 0) fclose(optimfile);

#if 0
   /* --- Finite differences test --- */
   printf("\n\n --- FINITE DIFFERENCE TESTING ---\n\n");

   double objective_orig, objective_perturb;
   double findiff, relerror;
   double errornorm = 0.0;


   /* Compute original objective and gradient */ 
   braid_Init(comm, comm, sstart, sstop, ntime, app,
            my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
            my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
   braid_InitAdjoint( my_ObjectiveT, my_ObjectiveT_diff, my_Step_diff, my_ResetGradient, &core);
   braid_SetPrintLevel( core, 1);
   braid_SetMaxLevels(core, max_levels);
   braid_SetNRelax(core, -1, nrelax);
   if (nrelax0 > -1)
   {
      braid_SetNRelax(core,  0, nrelax0);
   }
   braid_SetAbsTol(core, tol);
   braid_SetCFactor(core, -1, cfactor);
   braid_SetMaxIter(core, max_iter);
   if (fmg)
   {
      braid_SetFMG(core);
   }
   if (storage >= -2)
   {
      braid_SetStorage(core, storage);
   }
      
   /* Run simulation */
   my_ResetGradient(app);
   braid_Drive(core);

   /* Get the perturbed objective function value */
   braid_GetObjective(core, &objective_orig);
   printf("original Objective: %1.14e\n", objective_orig);

   /* Destroy new braid instance */
   braid_Destroy(core);


   /* Store original design and gradient */
   for (int idx = 0; idx < ndisc+1; idx++)
   {
      design0[idx]   = app->design[idx];  
      gradient0[idx] = app->gradient[idx];
   }

   /* FD step size */
   double EPS = 1e-6;
   /* Iterate over all design elements */
   // int idx = 0;  
   for (int idx = 0; idx < ndisc+1; idx ++)
   {
      
      /* perturb design */
      app->design[idx] += EPS;

      /* Create and run another braid instance */ 
      braid_Init(comm, comm, sstart, sstop, ntime, app,
               my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
               my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
      braid_InitAdjoint( my_ObjectiveT, my_ObjectiveT_diff, my_Step_diff, my_ResetGradient, &core);
      braid_SetPrintLevel( core, 1);
      braid_SetMaxLevels(core, max_levels);
      braid_SetNRelax(core, -1, nrelax);
      if (nrelax0 > -1)
      {
         braid_SetNRelax(core,  0, nrelax0);
      }
      braid_SetAbsTol(core, tol);
      braid_SetCFactor(core, -1, cfactor);
      braid_SetMaxIter(core, max_iter);
      if (fmg)
      {
         braid_SetFMG(core);
      }
      if (storage >= -2)
      {
         braid_SetStorage(core, storage);
      }
      /* Run simulation */
      my_ResetGradient(app);
      braid_Drive(core);
      /* Get the perturbed objective function value */
      braid_GetObjective(core, &objective_perturb);
      printf("perturbed Objective: %1.14e\n", objective_perturb);

      /* Destroy new braid instance */
      braid_Destroy(core);

      /* FD */
      findiff = (objective_perturb - objective_orig) / EPS;
      relerror = (gradient0[idx] - findiff) / findiff;
      printf("Finite Difference Test: \n");
      printf("%03d: gradient %1.14e findiff %1.14e \n", idx, gradient0[idx], findiff);
      printf("%03d: rel. error %1.14e\n", idx, relerror);
      errornorm += pow(relerror, 2);

      /* Reset the design */
      app->design[idx] = design0[idx];
   }

   errornorm = sqrt(errornorm);
   printf("\n Global errornorm %1.14e\n", errornorm);

#endif

   /* Clean up */
   delete hessian;
   free(app);
   free(design);
   free(gradient);
   free(design0);
   free(gradient0);
   free(ascentdir);
   MPI_Finalize();

   return (0);
}
