/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. Email xbraid-support@llnl.gov for support.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER*/

/**
 * Example:       ex-01-expanded.c
 *
 * Interface:     C
 * 
 * Requires:      only C-language support     
 *
 * Compile with:  make ex-01-expanded
 *
 * Help with:     ex-01-expanded -help
 *
 * Sample run:    mpirun -np 2 ex-01-expanded
 *
 * Description:   solve the scalar ODE 
 *                   u' = lambda u, 
 *                   y(0) = 1.0
 *                   lambda(0) = 1.0
 *
 *                   Nonlinear lambda to implement the thermostat problem
 *                   lambda(t) <--  -1   if  u(t)  >= 2
 *                   lambda(t) <--   1   if  u(t)  <= 1
 *                   lambda(t) unchanged, otherwise (equal to limit value of lambda from left)
 *
 *                Same as ex-01, only implements more advanced XBraid features.
 *                
 *                When run with the default 10 time steps, the solution is:
 *                $ ./ex-01-expanded
 *                $ cat ex-01-expanded.out.00*
 *                  1.00000000000000e+00
 *                  6.66666666666667e-01
 *                  4.44444444444444e-01
 *                  2.96296296296296e-01
 *                  1.97530864197531e-01
 *                  1.31687242798354e-01
 *                  8.77914951989026e-02
 *                  5.85276634659351e-02
 *                  3.90184423106234e-02
 *                  2.60122948737489e-02
 *                  1.73415299158326e-02
 **/


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
   double   *t;
   int       mydt;
   int       ntime;     /* Number of uniform time points */
   int       ndisc;     /* Number of discontinuities */
   int       rank;

} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
   double value;  // PDE solution value
   double coeff;  // PDE coefficient
   int    iter, level;
} my_Vector;

/* Generate time grid, take union of uniform time grid with discontinuity locations */ 
int
init_TimeSteps(braid_App  app)
{
   
   int ntime = app->ntime;
   double dt = (app->tstop - app->tstart) / ntime;
   double logTwo = log(2.0); 
   int i;
   int counter = 1;
   (app->t) = (double*) malloc( (ntime + app->ndisc + 1)*sizeof(double));
   
   app->t[0] = 0.0;
   for (i = 1; i <= ntime + app->ndisc; i++)
   {
      if( (app->t[i-1] + dt) > counter*logTwo)
      {
         app->t[i] = counter*logTwo;
         counter += 1;
         i += 1;
         app->t[i] = app->t[i-2] + dt;
      }
      else
      {
         app->t[i] = app->t[i-1] + dt;
      }
   }
   
   return 0;
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
   
   /* Store the level and iter for this Step, so that the correct precedence
    * rules for inheriting u->coeff in sum can be used */
   u->level = level;
   u->iter = iter;

   /* Account for XBraid right-hand-side */
   if (fstop != NULL)
   {
      (u->value) += (fstop->value);
   }

   
   /* Compute nonlinear PDE coefficient (based on ustop, or u?) */
   if (u->value >= 2.0) 
   {
      u->coeff = -1.0;
   }
   else if (u->value <= 1.0)
   {
      u->coeff = 1.0;
   }
   /* else do nothing, keep existing value of u->coeff */

   /* Use backward Euler to propagate solution */
   (u->value) = 1./(1. + (-u->coeff)*(tstop-tstart))*(u->value);


   /* no refinement */
   braid_StepStatusSetRFactor(status, 1);

   return 0;
}

int
my_Residual(braid_App        app,
            braid_Vector     ustop,
            braid_Vector     r,
            braid_StepStatus status)
{
   double tstart;             /* current time */
   double tstop;              /* evolve to this time*/
   braid_StepStatusGetTstartTstop(status, &tstart, &tstop);

   /* compute A(u_i, u_{i-1}) */
   (r->value) = (1. + (-ustop->coeff)*(tstop-tstart))*(ustop->value) - (r->value);

   return 0;
}

int
print_my_timegrid(braid_App        app,
                  braid_Real      *ta,
                  braid_Int       *ilower,
                  braid_Int       *iupper)
{
   int   i, lower, upper;
   char  filename[255];
   FILE *file;
   
   lower = *ilower;
   upper = *iupper;

   /* filename could be anything that helps you track the current time grid */
   sprintf(filename, "timegrid.%d.info", app->rank);
   file = fopen(filename, "w");
   if (file != NULL) {
      for (i = lower; i <= upper; i++)
      {
         fprintf(file, "%d %f\n", i, ta[i-lower]);
      }
   }
   else
   {
      return 1;
   }
   fclose(file);

   return 0;
}

int
my_timegrid(braid_App        app,
            braid_Real      *ta,
            braid_Int       *ilower,
            braid_Int       *iupper)
{
   int i, lower, upper;
   
   lower = *ilower;
   upper = *iupper;
   
   /* Assign time point values for local time point index values lower:upper */
   for (i = lower; i <= upper; i++)
   {
      ta[i-lower]  = app->t[i-lower];
   }
   print_my_timegrid(app, ta, &lower, &upper);

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
   (u->level) = 999.0;
   (u->iter)  = 999.0;
   
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
   int        index, iter, level;
   double     t;
   char       filename[255];
   FILE      *file;
   
   braid_AccessStatusGetTIndex(astatus, &index);
   braid_AccessStatusGetT(astatus, &t);
   braid_AccessStatusGetIter(astatus, &iter);
   braid_AccessStatusGetLevel(astatus, &level);
   // sprintf(filename, "%s.%04d.%03d", "ex-01-expanded.out", index, app->rank);
   // file = fopen(filename, "w");
   // fprintf(file, "%.14e\n", (u->value));
   // fflush(file);
   // fclose(file);

   /* Append all into one file. */
   if (level == 0)
   {
      sprintf(filename, "%s.%03d.%03d", "ex-01-expanded.out", iter, app->rank);
      file = fopen(filename, "a");
      fprintf(file, "%04d %1.4f  %.14e\n", index, t, (u->value));
      fflush(file);
      fclose(file);
   }
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
   double *dbuffer = (double*) buffer;

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
   double    *dbuffer = (double*) buffer;
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
   int           res        = 0;
   int           mydt       = 0;
   int           storage    = -1;
   int           print_level= 1;

   int           arg_index;
   int           rank;

   /* Initialize MPI */
   comm   = MPI_COMM_WORLD;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(comm, &rank);

   /* Define time domain: ntime intervals */
   ntime  = 400;
   tstart = 0.0;
   tstop  = 5.0; 
   
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
            printf("  -pl  <printlevel> : set print level\n");
            printf("  -mi  <max_iter>   : set max iterations\n");
            printf("  -fmg              : use FMG cycling\n");
            printf("  -res              : use my residual\n");
            printf("  -tg               : use a uniform time grid unioned with the discontinuity locations\n");
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
      else if ( strcmp(argv[arg_index], "-fmg") == 0 )
      {
         arg_index++;
         fmg = 1;
      }
      else if ( strcmp(argv[arg_index], "-res") == 0 )
      {
         arg_index++;
         res = 1;
      }
      else if ( strcmp(argv[arg_index], "-storage") == 0 ){
         arg_index++;
         storage = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tg") == 0 )
      {
         arg_index++;
         mydt = 1;
      }
      else
      {
         arg_index++;
      }
   }

   /* set up app structure */
   app = (my_App *) malloc(sizeof(my_App));
   (app->t)      = NULL;
   (app->mydt)   = mydt;
   (app->comm)   = comm;
   (app->tstart) = tstart;
   (app->tstop)  = tstop;
   (app->ntime)  = ntime;
   (app->rank)   = rank;
   if (app->mydt > 0)
   {
      (app->ndisc)  = floor(tstop / log(2.0));
   }
   else
   {
      (app->ndisc)  = 0; 
   }
   /* initialize XBraid and set options */
   braid_Init(comm, comm, tstart, tstop, ntime + app->ndisc, app,
             my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
             my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

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
   if (res)
   {
      braid_SetResidual(core, my_Residual);
   }
   if (app->mydt > 0)
   {
      init_TimeSteps(app);
      braid_SetTimeGrid(core, my_timegrid);
   }

   /* Run simulation, and then clean up */
   braid_SetAccessLevel(core, 2);
   braid_Drive(core);

   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}
