from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from math import sqrt
from mpi4py import MPI
cimport mpi4py.MPI as MPI
cimport mpi4py.libmpi as libmpi

import time
import sys
import os
import pickle

from spectral_toolbox_1D import SpectralToolbox
from RSWE_exponential_integrator import ExponentialIntegrator
from scipy.signal import gaussian
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import RSWE_direct
import cyclops_base
import numpy as np
import cyclops_control

#setup github to push this & 'sequentialized' cyclops lite to xbraid sandbox

cdef struct _braid_Vector_struct:
    double* v1
    double* v2
    double* h
ctypedef _braid_Vector_struct my_Vector
ctypedef _braid_Vector_struct *braid_Vector

cdef struct _braid_App_struct:
    int rank
ctypedef _braid_App_struct my_App
ctypedef _braid_App_struct *braid_App

include "braid.pyx"

control = cyclops_control.setup_control(None)
expInt = ExponentialIntegrator(control)
st = SpectralToolbox(control['Nx'], control['Lx'])
ICs = cyclops_base.h_init(control)

cdef int my_step(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status):
    cdef double tstart
    cdef double tstop
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop)
    cdef double[:] v1_view = <double[:control['Nx']]> u.v1
    cdef double[:] v2_view = <double[:control['Nx']]> u.v2
    cdef double[:] h_view = <double[:control['Nx']]> u.h

    v1_arr = np.asarray(v1_view)
    v2_arr = np.asarray(v2_view)
    h_arr = np.asarray(h_view)

    output = RSWE_direct.solve(control, expInt, st, np.array([v1_arr,v2_arr,h_arr]), solver = 'fine_propagator', realspace = False)
    
    v1_arr[:] = output[0,:]
    v2_arr[:] = output[1,:]
    h_arr[:] = output[2,:]

    return 0

cdef int my_init(braid_App app, double t, braid_Vector *u_ptr):
    cdef my_Vector* u
    u = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))

    cdef double[:] v1_view = <double[:control['Nx']]> u.v1
    cdef double[:] v2_view = <double[:control['Nx']]> u.v2
    cdef double[:] h_view = <double[:control['Nx']]> u.h

    v1_arr = np.asarray(v1_view)
    v2_arr = np.asarray(v2_view)
    h_arr = np.asarray(h_view)

    if (t == 0.0):
        v1_arr[:] = ICs[0,:]
        v2_arr[:] = ICs[1,:]
        h_arr[:] = ICs[2,:]
    else:
        v1_arr[:] = np.random.rand(control['Nx'])
        v2_arr[:] = np.random.rand(control['Nx'])
        h_arr[:] = np.random.rand(control['Nx'])

    #Copy ICs to v1[:] v2[:] h[:] for time 0, else skip option/randomness

    return 0

cdef int my_clone(braid_App app, braid_Vector u, braid_Vector *v_ptr):
    cdef my_Vector* v
    v = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))

    cdef double[:] u_v1_view = <double[:control['Nx']]> u.v1
    cdef double[:] u_v2_view = <double[:control['Nx']]> u.v2
    cdef double[:] u_h_view = <double[:control['Nx']]> u.h

    u_v1_arr = np.asarray(u_v1_view)
    u_v2_arr = np.asarray(u_v2_view)
    u_h_arr = np.asarray(u_h_view)
        
    cdef double[:] v_v1_view = <double[:control['Nx']]> v.v1
    cdef double[:] v_v2_view = <double[:control['Nx']]> v.v2
    cdef double[:] v_h_view = <double[:control['Nx']]> v.h

    v_v1_arr = np.asarray(v_v1_view)
    v_v2_arr = np.asarray(v_v2_view)
    v_h_arr = np.asarray(v_h_view)
   
    v_v1_arr[:] = u_v1_arr[:]
    v_v2_arr[:] = u_v2_arr[:]
    v_h_arr[:] = u_h_arr[:]

    # copy over 3 vectors using np.asarray numpy-style assignments
    v_ptr[0] = v
    return 0

cdef int my_free(braid_App app, braid_Vector u):
    PyMem_Free(u.v1)
    PyMem_Free(u.v2)
    PyMem_Free(u.h)
    PyMem_Free(u)
    return 0

cdef int my_sum(braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y):
    #y.value = alpha*x.value + beta*y.value
    cdef double[:] x_v1_view = <double[:control['Nx']]> x.v1
    cdef double[:] x_v2_view = <double[:control['Nx']]> x.v2
    cdef double[:] x_h_view = <double[:control['Nx']]> x.h

    x_v1_arr = np.asarray(x_v1_view)
    x_v2_arr = np.asarray(x_v2_view)
    x_h_arr = np.asarray(x_h_view)
        
    cdef double[:] y_v1_view = <double[:control['Nx']]> y.v1
    cdef double[:] y_v2_view = <double[:control['Nx']]> y.v2
    cdef double[:] y_h_view = <double[:control['Nx']]> y.h

    y_v1_arr = np.asarray(y_v1_view)
    y_v2_arr = np.asarray(y_v2_view)
    y_h_arr = np.asarray(y_h_view)

    y_v1_arr[:] = alpha*x_v1_arr[:] + beta*y_v1_arr[:]
    y_v2_arr[:] = alpha*x_v2_arr[:] + beta*y_v2_arr[:]
    y_h_arr[:] = alpha*x_h_arr[:] + beta*y_h_arr[:]

    #use np.asarray to do np-style sum
    return 0

cdef int my_norm(braid_App app, braid_Vector u, double *norm_ptr):
    cdef double dot

    cdef double[:] v1_view = <double[:control['Nx']]> u.v1
    cdef double[:] v2_view = <double[:control['Nx']]> u.v2
    cdef double[:] h_view = <double[:control['Nx']]> u.h

    v1_arr = np.asarray(v1_view)
    v2_arr = np.asarray(v2_view)
    h_arr = np.asarray(h_view)
    
    dot = np.linalg.norm(np.hstack([v1_arr,v2_arr,h_arr]))

    norm_ptr[0] = dot
    return 0

cdef int my_access(braid_App app, braid_Vector u, braid_AccessStatus status):
    print('my_access placeholder')
    #np.savetext to dump vectors
    return 0

cdef int my_bufsize(braid_App app, int *size_ptr, braid_BufferStatus status):
    size_ptr[0] = sizeof(double)*3*control['Nx'] # Need Nx (Control)
    return 0

cdef int my_bufpack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status):
    cdef double *dbuffer = <double*> buffer

    #cdef double[:] v1_view = <double[:control['Nx']]> u.v1
    #cdef double[:] v2_view = <double[:control['Nx']]> u.v2
    #cdef double[:] h_view = <double[:control['Nx']]> u.h

    #v1_arr = np.asarray(v1_view)
    #v2_arr = np.asarray(v2_view)
    #h_arr = np.asarray(h_view)

    #dbuffer[0:control['Nx']] = v1_arr[:]
    #dbuffer[control['Nx']:2*control['Nx']] = v2_arr[:]
    #dbuffer[2*control['Nx']:3*control['Nx']] = h_arr[:]

    for i in range(control['Nx']):
        dbuffer[i] = u.v1[i]
        dbuffer[control['Nx']+i] = u.v2[i]
        dbuffer[2*control['Nx']+i] = u.h[i]


    braid_BufferStatusSetSize(status, sizeof(double))
    return 0

cdef int my_bufunpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status):
    cdef double *dbuffer = <double*> buffer
    cdef my_Vector *u
    u = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))

    for i in range(control['Nx']):
        u.v1[i] = dbuffer[i]
        u.v2[i] = dbuffer[control['Nx']+i]
        u.h[i] = dbuffer[2*control['Nx']+i]

    #cdef double[:] v1_view = <double[:control['Nx']]> u.v1
    #cdef double[:] v2_view = <double[:control['Nx']]> u.v2
    #cdef double[:] h_view = <double[:control['Nx']]> u.h

    #v1_arr = np.asarray(v1_view)
    #v2_arr = np.asarray(v2_view)
    #h_arr = np.asarray(h_view)

    #v1_arr[:] = dbuffer[0:control['Nx']] 
    #v2_arr[:] = dbuffer[control['Nx']:2*control['Nx']] 
    #h_arr[:] = dbuffer[2*control['Nx']:3*control['Nx']]

    u_ptr[0] = u
    return 0

def braid_init_py():
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
    
    braid_SetMaxLevels(core, 1)

    braid_Drive(core)
    
