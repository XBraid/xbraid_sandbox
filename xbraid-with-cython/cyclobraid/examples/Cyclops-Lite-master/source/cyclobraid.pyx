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

    double* v1_imag
    double* v2_imag
    double* h_imag

ctypedef _braid_Vector_struct *braid_Vector
ctypedef _braid_Vector_struct my_Vector

cdef struct _braid_App_struct:
    int rank

ctypedef _braid_App_struct my_App
ctypedef _braid_App_struct *braid_App

include "braid.pyx"

# Set up global objects that traditionally go into "app", but in Python can be
# easily declared here.
control = cyclops_control.setup_control(None)
expInt = ExponentialIntegrator(control)
st = SpectralToolbox(control['Nx'], control['Lx'])
ICs = cyclops_base.h_init(control)

if control['HMM_T0'] is None:
    control['HMM_T0'] = control['coarse_timestep']/(control['epsilon']**0.2)
control['HMM_M_bar'] = max(25, int(80*control['HMM_T0']))
 
# Switch to turn on/off the time-averaged coarse-grid correction due to Terry and Beth
HMM_METHOD = 1 #On
#HMM_METHOD = 0 #Off

cdef int my_step(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status):
    cdef double tstart
    cdef double tstop
    cdef int level
    
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop)
    braid_StepStatusGetLevel(status, &level)

    #print("tstart = "+str(tstart))
    #print("app.rank = "+str(app.rank))
    #print('\n')
    
    ##
    # Create Python data type "memory views" of the double * arrays in braid_Vector
    cdef double[:] v1_view = <double[:control['Nx']]> u.v1
    cdef double[:] v2_view = <double[:control['Nx']]> u.v2
    cdef double[:] h_view = <double[:control['Nx']]> u.h

    v1_arr = np.asarray(v1_view)
    v2_arr = np.asarray(v2_view)
    h_arr = np.asarray(h_view)

    cdef double[:] v1_i_view = <double[:control['Nx']]> u.v1_imag
    cdef double[:] v2_i_view = <double[:control['Nx']]> u.v2_imag
    cdef double[:] h_i_view = <double[:control['Nx']]> u.h_imag

    v1_i_arr = np.asarray(v1_i_view)
    v2_i_arr = np.asarray(v2_i_view)
    h_i_arr = np.asarray(h_i_view)

    state = np.zeros((3,control['Nx']), dtype='complex')
    
    ##
    # Copy braid_Vector into state
    state[0,:].real = v1_arr[:] 
    state[1,:].real = v2_arr[:] 
    state[2,:].real = h_arr[:]
    state[0,:].imag = v1_i_arr[:] 
    state[1,:].imag = v2_i_arr[:] 
    state[2,:].imag = h_i_arr[:] 

    control['final_time'] = tstop
    control['start_time'] = tstart

    if (HMM_METHOD):
        if level == 0:
            control['fine_timestep'] = tstop - tstart
            output = RSWE_direct.solve(control, expInt, st, state, solver = 'fine_propagator', realspace = False)
        else:
            control['coarse_timestep'] = tstop - tstart
            output = RSWE_direct.solve(control, expInt, st, state, solver = 'coarse_propagator', realspace = False)
    else:
        control['fine_timestep'] = tstop - tstart
        output = RSWE_direct.solve(control, expInt, st, state, solver = 'fine_propagator', realspace = False)


    ##
    # Copy output into braid_Vector
    v1_arr[:] = output[0,:].real
    v2_arr[:] = output[1,:].real
    h_arr[:] = output[2,:].real
    v1_i_arr[:] = output[0,:].imag
    v2_i_arr[:] = output[1,:].imag
    h_i_arr[:] = output[2,:].imag

    return 0

cdef int my_init(braid_App app, double t, braid_Vector *u_ptr):
    cdef my_Vector* u
    u = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))

    u.v1 = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.v2 = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.h = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.v1_imag = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.v2_imag = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.h_imag = <double*>PyMem_Malloc(control['Nx']*sizeof(double))

    #for k in range(control['Nx']):
    #    print(k)
    #    u.v1[k] = 0. 
    #    u.v2[k] = 0. 
    #    u.h[k] = 0. 
    #    u.v1_imag[k] = 0.
    #    u.v2_imag[k] = 0.
    #    u.h_imag[k] = 0. 

    cdef double[:] v1_view = <double[:control['Nx']]> u.v1
    cdef double[:] v2_view = <double[:control['Nx']]> u.v2
    cdef double[:] h_view = <double[:control['Nx']]> u.h

    v1_arr = np.asarray(v1_view)
    v2_arr = np.asarray(v2_view)
    h_arr = np.asarray(h_view)

    cdef double[:] v1_i_view = <double[:control['Nx']]> u.v1_imag
    cdef double[:] v2_i_view = <double[:control['Nx']]> u.v2_imag
    cdef double[:] h_i_view = <double[:control['Nx']]> u.h_imag

    v1_i_arr = np.asarray(v1_i_view)
    v2_i_arr = np.asarray(v2_i_view)
    h_i_arr = np.asarray(h_i_view)

    state = np.zeros((3, control['Nx']), dtype='complex')
    for k in range(3):
        state[k,:] = st.forward_fft(ICs[k,:])

    if (t == 0.0):
        v1_arr[:] = state[0,:].real
        v2_arr[:] = state[1,:].real
        h_arr[:] = state[2,:].real
    else:
        v1_arr[:] = np.random.rand(control['Nx'])
        v2_arr[:] = np.random.rand(control['Nx'])
        h_arr[:] = np.random.rand(control['Nx'])

    if (t == 0.0):
        v1_i_arr[:] = state[0,:].imag
        v2_i_arr[:] = state[1,:].imag
        h_i_arr[:] = state[2,:].imag
    else:
        v1_i_arr[:] = np.random.rand(control['Nx'])
        v2_i_arr[:] = np.random.rand(control['Nx'])
        h_i_arr[:] = np.random.rand(control['Nx'])

    #Copy ICs to v1[:] v2[:] h[:] for time 0, else skip option/randomness

    u_ptr[0] = u
    return 0

cdef int my_clone(braid_App app, braid_Vector u, braid_Vector *v_ptr):
    cdef my_Vector* v
    v = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))

    v.v1 = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    v.v2 = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    v.h = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    v.v1_imag = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    v.v2_imag = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    v.h_imag = <double*>PyMem_Malloc(control['Nx']*sizeof(double))


    cdef double[:] u_v1_view = <double[:control['Nx']]> u.v1
    cdef double[:] u_v2_view = <double[:control['Nx']]> u.v2
    cdef double[:] u_h_view = <double[:control['Nx']]> u.h

    u_v1_arr = np.asarray(u_v1_view)
    u_v2_arr = np.asarray(u_v2_view)
    u_h_arr = np.asarray(u_h_view)
 
    cdef double[:] u_v1_i_view = <double[:control['Nx']]> u.v1_imag
    cdef double[:] u_v2_i_view = <double[:control['Nx']]> u.v2_imag
    cdef double[:] u_h_i_view = <double[:control['Nx']]> u.h_imag

    u_v1_i_arr = np.asarray(u_v1_i_view)
    u_v2_i_arr = np.asarray(u_v2_i_view)
    u_h_i_arr = np.asarray(u_h_i_view)
         
    cdef double[:] v_v1_view = <double[:control['Nx']]> v.v1
    cdef double[:] v_v2_view = <double[:control['Nx']]> v.v2
    cdef double[:] v_h_view = <double[:control['Nx']]> v.h

    v_v1_arr = np.asarray(v_v1_view)
    v_v2_arr = np.asarray(v_v2_view)
    v_h_arr = np.asarray(v_h_view)
        
    cdef double[:] v_v1_i_view = <double[:control['Nx']]> v.v1_imag
    cdef double[:] v_v2_i_view = <double[:control['Nx']]> v.v2_imag
    cdef double[:] v_h_i_view = <double[:control['Nx']]> v.h_imag

    v_v1_i_arr = np.asarray(v_v1_i_view)
    v_v2_i_arr = np.asarray(v_v2_i_view)
    v_h_i_arr = np.asarray(v_h_i_view)
   
    v_v1_arr[:] = u_v1_arr[:]
    v_v2_arr[:] = u_v2_arr[:]
    v_h_arr[:] = u_h_arr[:]
   
    v_v1_i_arr[:] = u_v1_i_arr[:]
    v_v2_i_arr[:] = u_v2_i_arr[:]
    v_h_i_arr[:] = u_h_i_arr[:]

    # copy over 3 vectors using np.asarray numpy-style assignments
    v_ptr[0] = v
    return 0

cdef int my_free(braid_App app, braid_Vector u):
    PyMem_Free(u.v1)
    PyMem_Free(u.v2)
    PyMem_Free(u.h)
    PyMem_Free(u.v1_imag)
    PyMem_Free(u.v2_imag)
    PyMem_Free(u.h_imag)
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

    cdef double[:] x_v1_i_view = <double[:control['Nx']]> x.v1_imag
    cdef double[:] x_v2_i_view = <double[:control['Nx']]> x.v2_imag
    cdef double[:] x_h_i_view = <double[:control['Nx']]> x.h_imag

    x_v1_i_arr = np.asarray(x_v1_i_view)
    x_v2_i_arr = np.asarray(x_v2_i_view)
    x_h_i_arr = np.asarray(x_h_i_view)
        
    cdef double[:] y_v1_i_view = <double[:control['Nx']]> y.v1_imag
    cdef double[:] y_v2_i_view = <double[:control['Nx']]> y.v2_imag
    cdef double[:] y_h_i_view = <double[:control['Nx']]> y.h_imag

    y_v1_i_arr = np.asarray(y_v1_i_view)
    y_v2_i_arr = np.asarray(y_v2_i_view)
    y_h_i_arr = np.asarray(y_h_i_view)

    y_v1_i_arr[:] = alpha*x_v1_i_arr[:] + beta*y_v1_i_arr[:]
    y_v2_i_arr[:] = alpha*x_v2_i_arr[:] + beta*y_v2_i_arr[:]
    y_h_i_arr[:] = alpha*x_h_i_arr[:] + beta*y_h_i_arr[:]

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
    cdef int index
    braid_AccessStatusGetTIndex(status,&index);
    cdef double[:] v1_view = <double[:control['Nx']]> u.v1
    cdef double[:] v2_view = <double[:control['Nx']]> u.v2
    cdef double[:] h_view = <double[:control['Nx']]> u.h

    v1_arr = np.asarray(v1_view)
    v2_arr = np.asarray(v2_view)
    h_arr = np.asarray(h_view)

    np.savetxt('cyclobraid.out-v1-'+str(index),st.inverse_fft(v1_arr[:]).real)
    np.savetxt('cyclobraid.out-v2-'+str(index),st.inverse_fft(v2_arr[:]).real) 
    np.savetxt('cyclobraid.out-h-'+str(index), st.inverse_fft(h_arr[:]).real)
    return 0

cdef int my_bufsize(braid_App app, int *size_ptr, braid_BufferStatus status):
    size_ptr[0] = sizeof(double)*6*control['Nx'] 
    return 0

cdef int my_bufpack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status):
    cdef double *dbuffer = <double*> buffer

    for i in range(control['Nx']):
        dbuffer[i] = u.v1[i]
        dbuffer[control['Nx']+i] = u.v2[i]
        dbuffer[2*control['Nx']+i] = u.h[i]
        dbuffer[3*control['Nx']+i] = u.v1_imag[i]
        dbuffer[4*control['Nx']+i] = u.v2_imag[i]
        dbuffer[5*control['Nx']+i] = u.h_imag[i]

    braid_BufferStatusSetSize(status, sizeof(double)*6*control['Nx'])
    return 0

cdef int my_bufunpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status):
    cdef double *dbuffer = <double*> buffer
    cdef my_Vector *u
    u = <my_Vector*>PyMem_Malloc(sizeof(my_Vector))
    u.v1 = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.v2 = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.h = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.v1_imag = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.v2_imag = <double*>PyMem_Malloc(control['Nx']*sizeof(double))
    u.h_imag = <double*>PyMem_Malloc(control['Nx']*sizeof(double))

    for i in range(control['Nx']):
        u.v1[i] = dbuffer[i]
        u.v2[i] = dbuffer[control['Nx']+i]
        u.h[i] = dbuffer[2*control['Nx']+i]
        u.v1_imag[i] = dbuffer[3*control['Nx']+i]
        u.v2_imag[i] = dbuffer[4*control['Nx']+i]
        u.h_imag[i] = dbuffer[5*control['Nx']+i]

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
 
    ntime = 50
    tstart = 0.0
    tstop = 5.0
    MPI_Comm_rank(comm.ob_mpi, &rank)

    app = <my_App*>PyMem_Malloc(sizeof(my_App))
    app.rank = rank
    
    # Update final_time value so Cyclops and Braid agree
    # Not sure that it's safe to overwrite other control values here 
    control['final_time'] = tstop
   
    braid_Init(comm.ob_mpi, comm.ob_mpi, tstart, tstop, ntime, app, my_step, my_init, my_clone, my_free, my_sum, my_norm, my_access, my_bufsize, my_bufpack, my_bufunpack, &core)
    
    braid_SetMaxLevels(core,2)
    braid_SetCFactor(core,-1,2)
    #braid_SetAbsTol(core,0.0)
    #braid_SetSeqSoln(core,1)
    #braid_SetMaxIter(core,3)

    ## Below options apply to multi-level runs
    braid_SetNRelax(core,-1,0)
    braid_SetSkip(core,1)

    braid_Drive(core)
    
