cdef extern from "braid/braid.h":
    ctypedef struct _braid_App_struct
    ctypedef struct* braid_App:
        pass

    struct _braid_Vector_struct
    ctypedef struct* braid_Vector:
        pass

    ctypedef braid_Int (*braid_PtFcnStep)(braid_App app, braid_Vector ustop, braid_Vector     fstop, braid_Vector u, braid_StepStatus status)

    braid_Int
braid_Init(MPI_Comm comm_world, MPI_Comm               comm, braid_Real             tstart,braid_Real tstop, braid_Int              ntime,braid_App              app,braid_PtFcnStep step, braid_PtFcnInit        init,braid_PtFcnClone       clone,    braid_PtFcnFree        free,   braid_PtFcnSum         sum,      braid_PtFcnSpatialNorm spatialnorm,braid_PtFcnAccess      access, braid_PtFcnBufSize     bufsize, braid_PtFcnBufPack     bufpack, braid_PtFcnBufUnpack   bufunpack, braid_Core            *core_ptr)
