cdef extern from "braid_status.h":

    cdef struct _braid_Status_struct:
        pass
    ctypedef _braid_Status_struct *braid_Status

    cdef struct _braid_AccessStatus_struct:
        pass
    ctypedef _braid_AccessStatus_struct *braid_AccessStatus 

    cdef struct _braid_StepStatus_struct:
        pass
    ctypedef _braid_StepStatus_struct *braid_StepStatus

    cdef struct _braid_BufferStatus_struct:
        pass
    ctypedef _braid_BufferStatus_struct *braid_BufferStatus

    int braid_StepStatusGetTstartTstop(braid_StepStatus status, double *tstart_ptr, double* tstop_ptr)
    int braid_AccessStatusGetT(braid_AccessStatus status, double *t_ptr) 
    int braid_BufferStatusSetSize( braid_BufferStatus status, int size);

cdef extern from "braid.h":
    
    cdef struct _braid_Core_struct:
        pass
    ctypedef _braid_Core_struct *braid_Core

    ctypedef int (*braid_PtFcnStep)(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status)

    ctypedef int (*braid_PtFcnInit)(braid_App app, double t, braid_Vector *u_ptr)
    
    ctypedef int (*braid_PtFcnClone)(braid_App app, braid_Vector u, braid_Vector *v_ptr)

    ctypedef int (*braid_PtFcnFree)(braid_App app, braid_Vector u)

    ctypedef int (*braid_PtFcnSum)(braid_App app, double alpha, braid_Vector x, double beta, braid_Vector y)

    ctypedef int (*braid_PtFcnSpatialNorm)(braid_App app, braid_Vector u, double *norm_ptr)

    ctypedef int (*braid_PtFcnAccess)(braid_App app, braid_Vector u, braid_AccessStatus status)

    ctypedef int (*braid_PtFcnBufSize)(braid_App app, int *size_ptr, braid_BufferStatus status)      

    ctypedef int (*braid_PtFcnBufPack)(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus status)

    ctypedef int (*braid_PtFcnBufUnpack)(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus status)

    int braid_Init(libmpi.MPI_Comm comm_world, libmpi.MPI_Comm comm, double tstart, double tstop, int ntime, braid_App app, braid_PtFcnStep step, braid_PtFcnInit init, braid_PtFcnClone clone, braid_PtFcnFree free, braid_PtFcnSum sum, braid_PtFcnSpatialNorm spatialnorm, braid_PtFcnAccess access,  braid_PtFcnBufSize bufsize, braid_PtFcnBufPack bufpack, braid_PtFcnBufUnpack bufunpack, braid_Core *core_ptr)
    int braid_SetMaxLevels(braid_Core core, int max_levels)

    int braid_Drive(braid_Core core)
