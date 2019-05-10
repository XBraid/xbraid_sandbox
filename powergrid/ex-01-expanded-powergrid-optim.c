
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"

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

} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double value;  // PDE solution value
   double coeff;  // PDE coefficient
} my_Vector;



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
   
   /* Use backward Euler and current coeff to propagate solution forward */
   (u->value) = 1./(1. + (-u->coeff)*(tstop-tstart))*(u->value);

   /* Compute new nonlinear PDE coefficient (based on ustop, or u?) */
   if (u->value >= 2.0) 
   {
      u->coeff = -1.0;
   }
   else if (u->value <= 1.0)
   {
      u->coeff = 1.0;
   }
   /* else do nothing, keep existing value of u->coeff */


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
   (u->coeff) = 1.0;
   
   if (t == 0.0) /* Initial condition */
   {
      (u->value) = 1.5;
   }
   else /* All other time points set to arbitrary value */
   {
      (u->value) = 0.0;
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
   (v->coeff) = (u->coeff);
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
   
   /* if(x->level < y->level)
   {
      y->coeff = x->coeff;
   } */
  
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
   *size_ptr = 2*sizeof(double);
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
   dbuffer[1] = (u->coeff);
   braid_BufferStatusSetSize( bstatus, 2*sizeof(double) );

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
   (u->coeff) = dbuffer[1];
   *u_ptr = u;

   return 0;
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

   int           max_levels = 2;
   int           nrelax     = 1;
   int           nrelax0    = -1;
   double        tol        = 1.0e-06;
   int           cfactor    = 2;
   int           max_iter   = 100;
   int           fmg        = 0;
   int           storage    = -1;

   int           arg_index;
   int           rank;

   /* Initialize MPI */
   comm   = MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(comm, &rank);

   /* Define time domain: ntime intervals */
   ntime  = 100;
   tstart = 0.0;
   tstop  = 2.0; 
   
   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         if ( rank == 0 )
         {
            printf("\nExample 1: Solve a scalar ODE \n\n");
            printf("  -ntime <ntime>    : set num time points\n");
            printf("  -ml  <max_levels> : set max levels\n");
            printf("  -nu  <nrelax>     : set num F-C relaxations\n");
            printf("  -nu0 <nrelax>     : set num F-C relaxations on level 0\n");
            printf("  -tol <tol>        : set stopping tolerance\n");
            printf("  -cf  <cfactor>    : set coarsening factor\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -fmg              : use FMG cycling\n");
            printf("  -storage <level>  : full storage on levels >= level\n");
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
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
      }
      else if ( strcmp(argv[arg_index], "-storage") == 0 ){
         arg_index++;
         storage = atoi(argv[arg_index++]);
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

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->comm)   = comm;
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;
   (app->rank)   = rank;


   /* initialize XBraid and set options */
   braid_Init(comm, comm, tstart, tstop, ntime, app,
            my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
            my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);
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
   braid_SetMaxLevels(core, 1); //serial timestepping only. 
   braid_Drive(core);


   /* Clean up */
   braid_Destroy(core);
   free(app);
   MPI_Finalize();

   return (0);
}
