from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython.ref cimport PyObject, Py_XINCREF, Py_INCREF, Py_DECREF
from mpi4py import MPI
cimport mpi4py.MPI as MPI
cimport mpi4py.libmpi as libmpi
import numpy as np 

##
# See ex_01-setup.py for notes on installing and running
##

##
# TODO for this example
#
# - can anything go into library routines, so that we can avoid non-Python language here, like
#    (1) dereferencing a pointer uptr[0] with a macro? 
#    (2) any casting between C and python types?
# - remove the braid app (?)
#

cdef class py_braid_Vector:
    '''
    Use a Python class to define your braid_Vector
    '''
    cdef object value

    def __cinit__(self): 
        self.value = np.array((11.1,), dtype='float') 

##
# Make braid_Vector just a PyObject* pointer 
# https://cython.readthedocs.io/en/latest/src/userguide/language_basics.html
ctypedef PyObject* braid_Vector


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
# Define user App 
#
# The Python example does not use app, like the C and Fortran codes.  Instead,
# define "Global" Python Objects that would normally go in App here, at the
# start of your program.  Then, you can use Variable below in my_step
##

##
# Define user functions like Step, Init, ...
##

cdef int my_step(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status):
    tstart = 0.0 # declare tstart and tstop as doubles
    tstop = 0.0
    braid_StepStatusGetTstartTstop(status, &tstart, &tstop)
    
    # Cast u as a py_braid_Vector, and do time-stepping.  Note
    # how the value array is written in-place
    pyU = <py_braid_Vector> u
    #print("Step " + str(pyU.value))
    pyU.value[0] = 1./(1. + tstop-tstart)*pyU.value[0]
    #print("   Step Done" + str(pyU.value))

    return 0


cdef int my_init(braid_App app, double t, braid_Vector *u_ptr):
    
    # Allocate new vector
    pyU = py_braid_Vector() 
    
    if (t == 0.0):
        pyU.value[0] = 1.0
    else:
        pyU.value[0] = 0.456
    Py_INCREF(pyU)
    
    #print("Init " + str(pyU.value))
    u_ptr[0] = <braid_Vector> pyU
    #Py_INCREF(<object>u_ptr[0])
    
    return 0

cdef int my_clone(braid_App app, braid_Vector u, braid_Vector *v_ptr):
    
    # Allocate new vector
    pyV = py_braid_Vector()
    
    # Cast u as a py_braid_Vector
    pyU = <py_braid_Vector> u

    pyV.value[0] = pyU.value[0]
    
    #print("Clone  " + str(pyU.value) + "  " + str(pyV.value))
    
    Py_INCREF(pyV)
    v_ptr[0] = <braid_Vector> pyV
    #Py_INCREF(<object>v_ptr[0])
    
    return 0
 
cdef int my_free(braid_App app, braid_Vector u):
    # Cast u as a py_braid_Vector
    #pyU = <py_braid_Vector> u
    #Py_DECREF(pyU)

    #del pyU 
    
    return 0

cdef int my_sum(braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y):
    # Cast x and y as a py_braid_Vector
    pyX = <py_braid_Vector> x
    pyY = <py_braid_Vector> y

    pyY.value[0] = alpha*pyX.value[0] + beta*pyY.value[0]
    #print("Sum " + str(pyY.value[0]) + ",  " + str(pyX.value[0]))
    return 0

cdef int my_norm(braid_App app, braid_Vector u, double *norm_ptr):
    # Cast u as a py_braid_Vector
    pyU = <py_braid_Vector> u

    norm_ptr[0] = np.sqrt(np.dot(pyU.value, pyU.value))
    #print("Norm " + str(norm_ptr[0]) + ",  " + str(pyU.value[0]))
    return 0

cdef int my_access(braid_App app, braid_Vector u, braid_AccessStatus status):
    # Cast u as a py_braid_Vector
    pyU = <py_braid_Vector> u

    # Declare t as a double, and the fill it in with a time
    t = 1.0
    braid_AccessStatusGetT(status, &t)

    print("Access " + str(pyU.value[0]) + ",  " + str(t))
    return 0

cdef int my_bufsize(braid_App app, int *size_ptr, braid_BufferStatus status):
    size_ptr[0] = sizeof(double)
    return 0

cdef int my_bufpack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status):
    # Cast u as a py_braid_Vector
    pyU = <py_braid_Vector> u

    
    # use the convert to array magic routine here...
    cdef double *dbuffer = <double*> buffer
    dbuffer[0] = pyU.value[0]

    braid_BufferStatusSetSize(status, sizeof(double))
    return 0

cdef int my_bufunpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status):
    
   # use the convert to array magic routine here...
    cdef double *dbuffer = <double*> buffer
    
    # Allocate new vector
    pyU = py_braid_Vector()
    Py_INCREF(pyU)
    #print("BufUnpack  " + str(pyU.value)) 

    pyU.value[0] = dbuffer[0]
    u_ptr[0] = <braid_Vector> pyU
    #Py_INCREF(<object>u_ptr[0])
    
    return 0

def InitCore():
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

    # Destroy Braid
    braid_Init(comm.ob_mpi, comm.ob_mpi, tstart, tstop, ntime, 
               app, my_step, my_init, my_clone, my_free, my_sum, 
               my_norm, my_access, my_bufsize, my_bufpack, 
               my_bufunpack, &core)

    # Store the Braid core inside of a Python-compatible object for return 
    pyCore = PyBraidCore()
    pyCore.setCore(core)

    return pyCore

def run_Braid(PyBraidCore pyCore):
    
    # Set Braid options
    braid_SetMaxLevels(pyCore.getCore(), 2)
    braid_SetMaxIter(pyCore.getCore(), 10)

    # Run Braid
    braid_Drive(pyCore.getCore())
    
    # Destroy Braid
    braid_Destroy(pyCore.getCore())


