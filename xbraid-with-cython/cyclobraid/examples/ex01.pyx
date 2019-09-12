from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from math import sqrt
from mpi4py import MPI
cimport mpi4py.MPI as MPI
cimport mpi4py.libmpi as libmpi
import numpy as np 
#cimport numpy as np


##
# Define Braid Vector
cdef struct _braid_Vector_struct:
    double* value
##
ctypedef _braid_Vector_struct* braid_Vector
ctypedef _braid_Vector_struct my_Vector

##
# Define Braid App
cdef struct _braid_App_struct:
    int rank
ctypedef _braid_App_struct my_App
ctypedef _braid_App_struct *braid_App

##
# Import Cython Braid Wrappers 
# (only do after declaration of Braid Vector and App)
include "braid.pyx"

##
# Begin Defining Step, Init, ...
##

##
# Define "Global" Python Objects that would normally go in App here.
##
# O = some object ... 
# Then, you can use O below in my_step

cdef int my_step(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status):
    cdef double tstart
    cdef double tstop
    tstart = 0.0
    tstop = 5.0
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop)
    
    # Convert u.value into numpy array using memory views
    # - Note, you have to cast u.value as a length-1 double array
    # - Note, this is ugly synatax, but val_arr and u.value share the same memory
    cdef double [:] val_view = <double[:1]> u.value
    val_arr = np.asarray(val_view)
    
    # Do time-stepping using numpy array
    # - Note how val_arr is written "in-place" with the colon, very important
    val_arr[:] = 1./(1. + tstop-tstart)*val_arr

    return 0

cdef int my_init(braid_App app, double t, braid_Vector *u_ptr):
    cdef my_Vector* u
    u = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))
    u.value = <double*>PyMem_Malloc(1*sizeof(double))
    if (t == 0.0):
        u.value[0] = 1.0
    else:
        u.value[0] = 0.456
    u_ptr[0] = u
    return 0

cdef int my_clone(braid_App app, braid_Vector u, braid_Vector *v_ptr):
    cdef my_Vector* v
    v = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))
    v.value = <double*>PyMem_Malloc(1*sizeof(double))
    v.value[0] = u.value[0]
    v_ptr[0] = v
    return 0

cdef int my_free(braid_App app, braid_Vector u):
    PyMem_Free(u.value)
    PyMem_Free(u)
    return 0

cdef int my_sum(braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y):
    y.value[0] = alpha*x.value[0] + beta*y.value[0]
    return 0

cdef int my_norm(braid_App app, braid_Vector u, double *norm_ptr):
    cdef double dot
    dot = u.value[0]*u.value[0]
    norm_ptr[0] = sqrt(dot)
    return 0

cdef int my_access(braid_App app, braid_Vector u, braid_AccessStatus status):
    cdef double t
    braid_AccessStatusGetT(status, &t)
    print(u.value[0], t)
    return 0

cdef int my_bufsize(braid_App app, int *size_ptr, braid_BufferStatus status):
    size_ptr[0] = sizeof(double)
    return 0

cdef int my_bufpack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status):
    cdef double *dbuffer = <double*> buffer
    dbuffer[0] = u.value[0]
    braid_BufferStatusSetSize(status, sizeof(double))
    return 0

cdef int my_bufunpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status):
    cdef double *dbuffer = <double*> buffer
    cdef my_Vector *u
    u = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))
    u.value = <double*>PyMem_Malloc(1*sizeof(double))
    u.value[0] = dbuffer[0]
    u_ptr[0] = u
    return 0

def braid_drive_py():
    cdef braid_Core core
    cdef my_App *app
    cdef double tstart
    cdef double tstop
    cdef MPI.Comm comm = MPI.COMM_WORLD
    cdef int ntime
    cdef int rank
 
    ntime = 10
    tstart = 0.0
    tstop = tstart + ntime/2.
    rank = 0
    app = <my_App*>PyMem_Malloc(sizeof(my_App))
    app.rank = rank

    braid_Init(comm.ob_mpi, comm.ob_mpi, tstart, tstop, ntime, app, my_step, my_init, my_clone, my_free, my_sum, my_norm, my_access, my_bufsize, my_bufpack, my_bufunpack, &core)
    
    braid_Drive(core)
    


