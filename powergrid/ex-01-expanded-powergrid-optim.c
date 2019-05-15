
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

/* On/Off switch for finite difference testing */
#define FINDEF 0

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
   MPI_Comm  comm;
   double    tstart;
   double    tstop;
   int       ntime;     /* Number of uniform time points */
   int       rank;
   double*   design;    /* Stores a(t) forall t */
   double*   gradient;  /* dJ/dDesign */

   double    gamma;     /* Regularization parameter */
   double    penalty_p; /* parameter of the penalty curve */ 

} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double value;  // PDE solution value
} my_Vector;


/* This drives a(t) towards -1 or 1, param defines steepness */
double penalty(int param, double a)
{
   /* Identity */
   // return a;

   /* a\in[0,1], exponential, then linear transform */
   // return 2 * pow(a, param) - 1;

   /* Sigmoid */
   return tanh(param * a);
}

double penalty_diff(int param, double a)
{
   /* Identity */
   // return 1.0;

   /* a\in[0,1], exponential, then linear transform */
   // return 2 * param * pow(a, param-1);

   /* Sigmoid */
   return param / pow(cosh(param * a), 2);
}

int
my_Step(braid_App        app,
        braid_Vector     ustop,
        braid_Vector     fstop,
        braid_Vector     u,
        braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   int    istart;             /* time point index value corresponding to tstop on (global) fine grid */
   int    level, iter;

   braid_StepStatusGetLevel(status, &level);
   braid_StepStatusGetIter(status, &iter);
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &istart);

   /* Get the design variable from the app */
   double design = app->design[istart];
   
   /* Push towards -1 or 1 */
   design = penalty(app->penalty_p, design); 
   
   /* Use backward Euler and current design to propagate solution forward */
   (u->value) = 1./(1. + (-design)*(tstop-tstart))*(u->value);

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
      (u->value) = 1.5;
   }
   else /* All other time points set to arbitrary value */
   {
      (u->value) = 1.5;
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
   (v->value) = (u->value);
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
   (y->value) = alpha*(x->value) + beta*(y->value);
   
   return 0;
}

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot;

   dot = (u->value)*(u->value);
   *norm_ptr = sqrt(dot);

   return 0;
}

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int        index;
   double     t;
   char       filename[255];
   FILE      *file;
   
   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetT(astatus, &t);
   // sprintf(filename, "%s.%04d.%03d", "ex-01-expanded.out", index, app->rank);
   // file = fopen(filename, "w");
   // fprintf(file, "%1.4f  %.14e\n", t, (u->value));
   // fflush(file);
   // fclose(file);

   /* Append all into one file. */
   sprintf(filename, "%s.%03d", "ex-01-expanded.out", app->rank);
   file = fopen(filename, "a");
   fprintf(file, "%04d %1.4f  %.14e\n", index, t, (u->value));
   fflush(file);
   fclose(file);


   // if (index == app->ntime) printf("Last: %d  %4f  %1.14e\n", index, t, u->value);

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
   double *dbuffer = buffer;

   dbuffer[0] = (u->value);
   braid_BufferStatusSetSize( bstatus, sizeof(double) );

   return 0;
}

int
my_BufUnpack(braid_App          app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus bstatus)
{
   double    *dbuffer = buffer;
   my_Vector *u;

   u = (my_Vector *) malloc(sizeof(my_Vector));
   (u->value) = dbuffer[0];
   *u_ptr = u;

   return 0;
}

int
my_ObjectiveT(braid_App app,
              braid_Vector u,
              braid_ObjectiveStatus ostatus,
              double *objectiveT_ptr)
{
   double objT;

   /* one norm */
   objT = fabs(u->value - 2) + fabs(u->value - 1) - 1.0;

   /* Regularization */
   int    idx;
   braid_ObjectiveStatusGetTIndex(ostatus, &idx);
   objT += app->gamma * fabs(app->design[idx] - 1.0);

   /* two norm */
   // objT  = ( u->value - 2.0 ) * ( u->value - 2.0 );
   // objT += ( 1.0 - u->value ) * ( 1.0 - u->value );
   // objT = 0.5 * objT;

   *objectiveT_ptr = objT;

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

   /* Partial wrt u times F_bar */
   /* one norm */
   if      (u->value < 1) u_bar->value = -1.0 * F_bar;
   else if (u->value > 2) u_bar->value =  1.0 * F_bar;
   else                   u_bar->value =  0.0;


   /* two norm */
   // u_bar->value = ( 2.0 * u->value - 3.0 ) * F_bar;

   /* Partial wrt design times F_bar -> Regularization . */
   int    idx;
   braid_ObjectiveStatusGetTIndex(ostatus, &idx);
   if (app->design[idx] < 1) app->gradient[idx] += -1.0 * app->gamma * F_bar;
   if (app->design[idx] > 1) app->gradient[idx] +=  1.0 * app->gamma * F_bar;

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
   double ddesign;  /* Derivative wrt design */
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   double dsigm;
   int istart;

   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
   braid_StepStatusGetTIndex(status, &istart);
   double deltat = tstop - tstart;

   /* Get the design from the app */
   double design = app->design[istart];

   /* Push towards -1 or 1 */
   design = penalty(app->penalty_p, design); 

   /* Get gradient of penalty */
   dsigm  = penalty_diff(app->penalty_p, app->design[istart]);

   /* Transposed derivative of step wrt u times u_bar */
   ddu = 1./(1. + (-design)*deltat)*(u_bar->value);

   /* Transposed derivative of step wrt design times u_bar */
   ddesign = dsigm * (deltat * (u->value)) / pow(1. - deltat*design,2) * (u_bar->value);

   /* Update u_bar and add to gradient */
   u_bar->value        = ddu;              
   app->gradient[istart] += ddesign;


   return 0;
}


/* Set the gradient to zero */
int 
my_ResetGradient(braid_App app)
{
   for (int i = 0; i < app->ntime; i++)
   {
      app->gradient[i] = 0.0;
   }

   return 0;
}

// * Function to allow for the computation of the gradient */
// int
// gradient_allreduce(braid_App app)
// {
//    double mygradient = app->gradient;
//    double gradient; 

//    /* Collect sensitivities from all processors and broadcast it */
//    MPI_Allreduce(&mygradient, &gradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    app->gradient = gradient;

//    return 0;
// }


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
   double        tstart, tstop;
   int           ntime;
   double        objective;

   double        gnorm;
   int           optimiter;
   int           arg_index;
   int           rank;
   double        mygradnorm;

   /* Setting defaults */
   int         max_levels = 2;
   int         nrelax     = 1;
   int         nrelax0    = -1;
   double      tol        = 1.0e-06;
   int         cfactor    = 2;
   int         max_iter   = 100;
   int         fmg        = 0;
   int         storage    = -1;

   double      penalty_p    = 10;     /* Parameter for penalty function */
   int         maxoptimiter = 100;    /* Maximum number of optimization iterations */
   double      stepsize     = 1.0;    /* Step size for design updates */
   double      gtol         = 1e-6;   /* Stopping criterion on the gradient norm */
   double      gamma        = 1e-4;   /* Regularization parameter */

   /* Define time domain: ntime intervals */
   ntime  = 100;
   tstart = 0.0;
   tstop  = 2.0; 

   /* Initialize MPI */
   comm   = MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(comm, &rank);

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
            printf("  -ml  <max_levels>      : set max levels\n");
            printf("  -nu  <nrelax>          : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>          : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>             : set stopping tolerance\n");
            printf("  -cf  <cfactor>         : set coarsening factor\n");
            printf("  -fmg                   : use FMG cycling\n");
            printf("  -storage <level>       : full storage on levels >= level\n");
            printf("  -mi  <max_iter>        : set max braid iterations\n");
            printf("  -moi <max_optim_iter>  : set max optimization iter\n");
            printf("  -dstep <stepsize>      : set step size for design updates\n");
            printf("  -gtol <gtol>           : set optimization stopping tolerance\n");
            printf("  -pp <penalty_param>    : set parameter for penalty function \n");
            printf("  -regul <regularization param> : set the regularization parameter \n");
         }
         exit(1);
      } 
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
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
         stepsize = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
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
      else if ( strcmp(argv[arg_index], "-pp") == 0 ) 
      {
         arg_index++;
         penalty_p = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-regul") == 0 ) 
      {
         arg_index++;
         gamma = atof(argv[arg_index++]);
      }
      else
      {
         arg_index++;
      }
   }

   /* Hack: Open and close output file. */
   char       filename[255];
   FILE      *file;
   sprintf(filename, "%s.%03d", "ex-01-expanded.out", rank);
   file = fopen(filename, "w");
   fclose(file);

   /* initialize design and gradient */
   double* design;
   double* gradient;
   design   = (double*) malloc(ntime * sizeof(double));
   gradient = (double*) malloc(ntime * sizeof(double));
   for (int i = 0; i < ntime; i++)
   {
      design[i]   = 1.0;  // Initial guess
      gradient[i] = 0.0;
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm)   = comm;
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;
   (app->rank)   = rank;
   (app->design) = design;
   (app->gradient) = gradient;
   (app->penalty_p) = penalty_p;
   (app->gamma)     = gamma;


   /* initialize XBraid */
   
   braid_Init(comm, comm, tstart, tstop, ntime, app,
            my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
            my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

  /* Initialize adjoint-based gradient computation */
   braid_InitAdjoint( my_ObjectiveT, my_ObjectiveT_diff, my_Step_diff, my_ResetGradient, &core);


   /* Set options */
   braid_SetPrintLevel( core, 0);
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

   /* Optimization iteration */
   for (optimiter = 0; optimiter < maxoptimiter; optimiter++)
   {
      /* Run adjoint XBraid to compute objective function and gradient */
      braid_SetAccessLevel(core, 0);
      braid_Drive(core);

      /* Get the objective function value */
      braid_GetObjective(core, &objective);

      /* Compute gradient norm */
      mygradnorm = 0.0;
      for (int id = 0; id < app->ntime; id++)
      {
         mygradnorm += pow(app->gradient[id], 2);
      }
      MPI_Allreduce(&mygradnorm, &gnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      gnorm = sqrt(gnorm);

      /* Output */
      if (rank == 0) 
      {
         printf("\n %3d: Objective = %1.8e,  || Gradient || = %1.8e\n", optimiter, objective, gnorm);
      }

      /* Check optimization convergence */
      if (gnorm < gtol)
      {
         break;
      }

      /* Design update using simple steepest descent method with fixed stepsize */
      for (int idx = 0; idx < app->ntime; idx++)
      {
         app->design[idx] -= stepsize * app->gradient[idx];
      }
   }


   /* Print some statistics about the optimization run */
   if (rank == 0)
   {
      if (optimiter == maxoptimiter)
      {
         printf("\n Max. number of iterations reached! \n\n"); 
      }
      else
      {
         printf("\n");
         printf("  Optimization has converged.\n");
         printf("\n"); 
         printf("  Objective function value = %1.8e\n", objective);
         printf("  Gradient norm            = %1.8e\n", gnorm);
         printf("\n");
         printf("  optimization iterations  = %d\n", optimiter);
         printf("  max optim iterations     = %d\n", maxoptimiter);
         printf("  gradient norm tolerance  = %1.1e", gtol);
         printf("\n");
      }
   }

   /* Print XBraid statistics */
   braid_PrintStats(core);

   /* Get final access */
   braid_SetAccessLevel(core, 1);
   braid_Drive(core);

   /* Finish braid */
   braid_Destroy(core);

   /* print */
   sprintf(filename, "%s.%03d", "design.out", rank);
   write_vector(filename, app->design, app->ntime);
   sprintf(filename, "%s.%03d", "gradient.out", rank);
   write_vector(filename, app->gradient, app->ntime);

   /* Print Penalty(design) */
   double* Sa = (double*) malloc(app->ntime * sizeof(double));
   for (int i=0; i<app->ntime; i++)
   {
      Sa[i] = penalty(app->penalty_p, app->design[i]);
   }
   sprintf(filename, "%s.%03d", "Sa.out", rank);
   write_vector(filename, Sa, app->ntime);



#if FINDEF
   /* --- Finite differences test --- */
   double objective_orig, objective_perturb;
   double findiff, relerror;
   double errornorm = 0.0;

   double EPS = 1e-6;   // FD step size

   /* Store original design and gradient */
   double* design0   = (double*) malloc(app->ntime * sizeof(double));
   double* gradient0 = (double*) malloc(app->ntime * sizeof(double));
   for (int idx = 0; idx < app->ntime; idx++)
   {
      design0[idx]   = app->design[idx];  
      gradient0[idx] = app->gradient[idx];
   }

   /* Iterate over all design elements */
   // int idx = 96;         // design element id
   for (int idx = 0; idx < app->ntime; idx ++)
   {
      /* Reset the design */
      app->design[idx] = design0[idx];

      /* Run first braid instance */ 
      braid_Init(comm, comm, tstart, tstop, ntime, app,
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
      braid_Drive(core);

      /* Get the perturbed objective function value */
      braid_GetObjective(core, &objective_orig);

      /* Destroy new braid instance */
      braid_Destroy(core);

      
      /* perturb design */
      app->design[idx] += EPS;

      /* Create and run another braid instance */ 
      braid_Init(comm, comm, tstart, tstop, ntime, app,
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
   }

   errornorm = sqrt(errornorm);
   printf("\n Global errornorm %1.14e\n", errornorm);

   /* Clean up Finite Differences */
   free(design0);
   free(gradient0);

#endif

   /* Clean up */
   free(app);
   free(design);
   free(gradient);
   MPI_Finalize();

   return (0);
}
