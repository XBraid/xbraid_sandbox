#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.h"
#include "braid_test.h"
#include "xbraid-SWE-lib.c"

/*--------------------------------------------------------------------------
 * User-defined routines and structures
 *--------------------------------------------------------------------------*/

/* App structure can contain anything, and be named anything as well */
typedef struct _braid_App_struct
{
	MPI_Comm  comm;
	double    tstart;       /* Define the temporal domain */
	double    tstop;
	int       ntime;
	double    xstart;       /* Define the spatial domain */
	double    xstop;
	int       nspace;
	double*   sc_info;
} my_App;

/* Vector structure can contain anything, and be name anything as well */
typedef struct _braid_Vector_struct
{
	int size;
	double* h;
	double* uh;
} my_Vector;

	void
create_vector(my_Vector **u, 
		int size)
{
	(*u) = (my_Vector *) malloc(sizeof(my_Vector));
	((*u)->size) = size;
	((*u)->h) = (double *) malloc(size*sizeof(double));
	((*u)->uh) = (double *) malloc(size*sizeof(double));
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
	my_Vector* prev;
	double g = 9.8;
	int level;
	braid_StepStatusGetLevel(status,&level);
	braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
	int size = (u->size);
	int i;
	double deltaX,deltaT;
	deltaT = tstop - tstart;
	deltaX = (app->xstop - app->xstart) / (ustop->size - 1.0);

	create_vector(&prev,size);

	double dx = (app->xstop - app->xstart)/(ustop->size-1.0);
	double dt = (tstop - tstart);
	my_Clone(app,u,&prev);

	for (i = 1; i < size-1; i++)
	{
		(u->h)[i] = (prev->h)[i] - dt*((prev->uh)[i+1] - (prev->uh)[i-1])/(2.0 * dx);
		(u->uh)[i] = (prev->uh)[i] - dt*((((prev->uh)[i+1]*(prev->uh)[i+1])/(prev->h)[i+1] +
					g*(prev->h)[i+1]*(prev->h)[i+1]*0.5) - 
				(((prev->uh)[i-1]*(prev->uh)[i-1])/(prev->h)[i-1] +
				 g*(prev->h)[i-1]*(prev->h)[i-1]* 0.5))/(2.0*dx);
	//	if ((abs(u->uh[i]/u->h[i]) + sqrt(g*u->h[i]))*dt/dx >= 1.) printf("Courant number large %d : %f\n",i,(abs(u->uh[i]/u->h[i]) + sqrt(g*u->h[i]))*dt/dx);
	//	else if ((abs(u->uh[i]/u->h[i]) + sqrt(g*u->h[i]))*dt/dx <= 0.1) printf("Courant number small %d : %f\n",i,(abs(u->uh[i]/u->h[i]) + sqrt(g*u->h[i]))*dt/dx);
	//	else printf("OK\n");
	}

	/*
	   Update the boundary conditions.
	   */
	// Periodic BC:
	//(u->h)[0] = (u->h)[size-2];
	//(u->h)[size-1] = (u->h)[1];
	//(u->uh)[0] = (u->uh)[size-2];
	//(u->uh)[size-1] = (u->uh)[1];

	/* Store info on space-time grids visited during the simulation */
	(app->sc_info)[ (2*level) ] = deltaX;
	(app->sc_info)[ (2*level) + 1] = deltaT;

	braid_StepStatusSetRFactor(status,1);


	// Free BC:
	(u->h)[0] = (prev->h)[1];
	(u->h)[size-1] = (prev->h)[size-2];
	(u->uh)[0] = (prev->uh)[1];
	(u->uh)[size-1] = (prev->uh)[size-2];
	return 0;
}

	int
my_Init(braid_App     app,
		double        t,
		braid_Vector *u_ptr)
{
	my_Vector *u;
	int nspace = (app->nspace); 
	double pi = 3.141592653589793;
	double* x = (double*) malloc(nspace*sizeof(double));
	// Make a linspace vector x
	if (nspace==1) {
		x[0] = ((app->xstop) + (app->xstart)) / 2.0;
	}
	else {
		for (int i = 0; i < nspace; i++) {
			x[i] = ((double)(nspace-1-i)*(app->xstart)
					+ (double)(i)*(app->xstop))
				/ (double)(nspace-1);
		}
	}
	create_vector(&u,nspace);
	if (t == 0.0) {
		for (int i = 0; i < nspace; i++) {
			(u->h)[i] = 2.0 + sin(.05*pi*x[i]); 
			(u->uh)[i] = u->h[i]*2.0;
		}
		// Periodic BC:
		//(u->h)[0] = (u->h)[nspace-2];
		//(u->h)[nspace-1] = (u->h)[1];
		//(u->uh)[0] = (u->uh)[nspace-2];
		//(u->uh)[nspace-1] = (u->uh)[1];
		// Free BC:
		(u->h)[0] = (u->h)[1];
		(u->h)[nspace-1] = (u->h)[nspace-2];
		(u->uh)[0] = (u->uh)[1];
		(u->uh)[nspace-1] = (u->uh)[nspace-2];
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
	int size = (u->size);
	int i;

	create_vector(&v, size);
	for (i = 0; i < size; i++)
	{
		(v->h)[i] = (u->h)[i];
		(v->uh)[i] = (u->uh)[i];
	}
	*v_ptr = v;

	return 0;
}

	int
my_Free(braid_App    app,
		braid_Vector u)
{
	free(u->h);
	free(u->uh);
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
	int size = (y->size);
	for (int i = 0; i < size; i++) {
		(y->h)[i] = alpha*(x->h)[i] + beta*(y->h)[i];
		(y->uh)[i] = alpha*(x->uh)[i] + beta*(y->uh)[i];
	}
	return 0;
}

	int
my_SpatialNorm(braid_App     app,
		braid_Vector  u,
		double       *norm_ptr)
{
	int size = (u->size);
	double dot = 0.;
	for (int i = 0; i < size; i++) {
		dot = dot + (u->h)[i]*(u->h)[i]; //+ (u->uh)[i]*(u->uh)[i]; 
	}
	*norm_ptr = sqrt(dot);
	return 0;
}

	int
my_Access(braid_App          app,
		braid_Vector       u,
		braid_AccessStatus astatus)
{
	int        index;
	char       filename[255];
	FILE      *file;

	braid_AccessStatusGetTIndex(astatus, &index);
	sprintf(filename, "000%d", index);
	file = fopen(filename, "w");

	int size = (u->size);
	for (int i = 0; i < size; i++) {
	    fprintf(file, "%.16e \n", u->h[i]);
	}
	fflush(file);
	fclose(file);
	return 0;
}

	int
my_BufSize(braid_App          app,
		int                *size_ptr,
		braid_BufferStatus bstatus)
{
	int size = (app->nspace);
	*size_ptr = (2*size+1)*sizeof(double);
	return 0;
}

	int
my_BufPack(braid_App          app,
		braid_Vector       u,
		void               *buffer,
		braid_BufferStatus bstatus)
{
	double *dbuffer = buffer;
	int size = (u->size);
	dbuffer[0] = size;
	for (int i = 0; i < size; i++) {
		dbuffer[i+1] = (u->h)[i];
		dbuffer[size+i+1] = (u->uh)[i];
	}
	braid_BufferStatusSetSize( bstatus, (2*size+1)*sizeof(double) );

	return 0;
}

	int
my_BufUnpack(braid_App          app,
		void               *buffer,
		braid_Vector       *u_ptr,
		braid_BufferStatus bstatus)
{
	double    *dbuffer = buffer;
	my_Vector *u = NULL;
	int size;

	size = dbuffer[0];
	create_vector(&u, size);

	for (int i = 0; i < size; i++)
	{
		(u->h)[i] = dbuffer[i+1];
		(u->uh)[i] = dbuffer[size+i+1];
	}
	*u_ptr = u;

	return 0;
}

/* Bilinear Coarsening */
	int
my_Coarsen(braid_App              app,           
		braid_Vector           fu,
		braid_Vector          *cu_ptr,
		braid_CoarsenRefStatus status)
{

	int i,fidx,csize,level;
	double* fh = fu->h;
	double* fuh = fu->uh;
	my_Vector *v;
	braid_CoarsenRefStatusGetLevel(status, &level);
	if( level < floor(log2(app->nspace)) - 1 )
	{
		csize = (fu->size - 1)/2 + 1;
		v = (my_Vector *) malloc(sizeof(my_Vector));
		(v->size)   = csize;
		(v->h) = (double *) malloc(csize*sizeof(double));
		(v->uh) = (double *) malloc(csize*sizeof(double));
		for (i = 1; i < csize-1; i++)
		{
			fidx = 2*i;
			(v->h)[i] = 0.5*fh[fidx] + 0.25*fh[fidx+1] + 0.25*fh[fidx-1];
			(v->uh)[i] = 0.5*fuh[fidx] + 0.25*fuh[fidx+1] + 0.25*fuh[fidx-1];
		}
		// Boundary conditions
		// Periodic BC
		//(v->h)[0]=.5*fh[0]+.25*fh[1]+.25*fh[(fu->size - 1)];
		//(v->uh)[0]=.5*fuh[0]+.25*fuh[1]+.25*fuh[(fu->size - 1)];
		//(v->h)[csize-1]=.5*fh[fu->size-1]+.25*fh[0]+.25*fh[(fu->size - 2)];
		//(v->uh)[csize-1]=.5*fuh[fu->size-1]+.25*fuh[0]+.25*fuh[(fu->size - 2)];
		//Free BC
		(v->h)[0]=fh[0];
		(v->uh)[0]=fuh[0];
		(v->h)[csize-1]=fh[fu->size-1];
		(v->uh)[csize-1]=fuh[fu->size-1];
	}
	else
	{
		/* No coarsening, clone the vector */
		my_Clone(app, fu, &v);
	}
	*cu_ptr = v;
	return 0;
}

/* Bilinear interpolation */
	int
my_Interp(braid_App              app,           
		braid_Vector           cu,
		braid_Vector          *fu_ptr,
		braid_CoarsenRefStatus status)
{

	int i, fsize, level;
	my_Vector *v;
	double* ch = cu->h;
	double* cuh = cu->uh;
	braid_CoarsenRefStatusGetLevel(status, &level);
	if( level < floor(log2(app->nspace)) - 1 )
	{
		fsize = (cu->size - 1)*2 + 1;
		v = (my_Vector *) malloc(sizeof(my_Vector));
		(v->size)   = fsize;
		(v->h) = (double *) malloc(fsize*sizeof(double));
		(v->uh) = (double *) malloc(fsize*sizeof(double));
		for (i = 1; i < fsize-1; i++)
		{
			if(i%2 == 1)
			{
				(v->h)[i] = 0.5*ch[i/2] + 0.5*ch[(i+1)/2];
				(v->uh)[i] = 0.5*cuh[i/2] + 0.5*cuh[(i+1)/2];
			}
			else
			{
				(v->h)[i] = ch[i/2];
				(v->uh)[i] = cuh[i/2];
			}

		}

		/* Boundary Conditions */
		(v->uh)[0] = cuh[0];
		(v->uh)[fsize-1] = cuh[cu->size-1];
		(v->h)[0] = ch[0];
		(v->h)[fsize-1] = ch[cu->size-1];
	}
	else
	{
		/* No refinement, clone the vector */
		my_Clone(app, cu, &v);
	}
	*fu_ptr = v;
	return 0;
}

	int
my_Residual(braid_App        app,
		braid_Vector     ustop,
		braid_Vector     r,
		braid_StepStatus status)
{
	return 0; }


	/*--------------------------------------------------------------------------
	 * Main driver
	 *--------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
	braid_Core    core;
	my_App       *app;
	double        tstart, tstop, loglevels;
	MPI_Comm      comm, comm_x, comm_t;
	int           ntime, rank;
	int           scoarsen = 0;

	/* Define time domain: ntime intervals */
	ntime  = 1000;
	tstart = 0.0;
	tstop  = 1.0;
	double    xstart        =  0.0;
	double    xstop         =  75.0;
	int      nspace        = 257;
	int print_level = 2;
	int max_levels = 2;


	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);

	/* set up app structure */
	app = (my_App *) malloc(sizeof(my_App));
	(app->comm)          = comm;
	(app->tstart)        = tstart;
	(app->tstop)         = tstop;
	(app->ntime)         = ntime;
	(app->xstart)        = xstart;
	(app->xstop)         = xstop;
	(app->nspace)        = nspace;


	/* Initialize storage for sc_info, for tracking space-time grids visited during the simulation */
	app->sc_info = (double*) malloc( 2*max_levels*sizeof(double) );
	for(int i = 0; i < 2*max_levels; i++) {
		app->sc_info[i] = -1.0;
	}


	/* initialize XBraid and set options */
	braid_Init(MPI_COMM_WORLD, comm, tstart, tstop, ntime, app,
			my_Step, my_Init, my_Clone, my_Free, my_Sum, my_SpatialNorm, 
			my_Access, my_BufSize, my_BufPack, my_BufUnpack, &core);

	/* Create spatial communicator for wrapper-tests */
	braid_SplitCommworld(&comm, 1, &comm_x, &comm_t);

	braid_TestAll(app, comm_x, stdout, 0.0, (tstop-tstart)/ntime,
			2*(tstop-tstart)/ntime, my_Init, my_Free, my_Clone,
			my_Sum, my_SpatialNorm, my_BufSize, my_BufPack,
			my_BufUnpack, my_Coarsen, my_Interp, my_Residual, my_Step);



	/* Scale tol by domain */
	double      tol = 1.0e-5;
	tol = tol/( sqrt((tstop - tstart)/(ntime-1))*sqrt((xstop - xstart)/(nspace-1)) );

	/* Set some typical Braid parameters */
	braid_SetPrintLevel(core, print_level);
	braid_SetMaxLevels(core, max_levels);
	braid_SetAbsTol(core, tol);
	//   braid_SetSeqSoln(core, 1);
	braid_SetCFactor(core, -1, 2);
	braid_SetNRelax(core, -1, 1);
	braid_SetMaxIter(core, 300);
	braid_SetMinCoarse(core,3);

	/* Spatial Coarsening */
	loglevels = log2(nspace - 1.0);
	if ( scoarsen && ( fabs(loglevels - round(loglevels)) > 1e-10 ))
	{
		if(rank == 0)
		{
			fprintf(stderr, "\nWarning!\nFor spatial coarsening, spatial grids must be a "
					"power of 2 + 1, \ni.e., nspace = 2^k + 1.  Your spatial grid is of size"
					" %d.  Spatial \ncoarsening is therefore being ignored.\n\n", nspace);
		}
	}
	else if (scoarsen)
	{
		braid_SetSpatialCoarsen(core, my_Coarsen);
		braid_SetSpatialRefine(core,  my_Interp);
	}

	/* Print accumulated info on space-time grids visited during the simulation */

	/* Run simulation, and then clean up */
	braid_Drive(core);

	print_sc_info(app->sc_info, max_levels);

	braid_Destroy(core);
	free(app->sc_info);
	free(app);
	MPI_Finalize();

	return (0);
}

