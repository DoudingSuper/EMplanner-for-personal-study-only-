/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Errhandler_get */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Errhandler_get = PMPI_Errhandler_get
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Errhandler_get  MPI_Errhandler_get
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Errhandler_get as PMPI_Errhandler_get
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler)
     __attribute__ ((weak, alias("PMPI_Errhandler_get")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Errhandler_get
#define MPI_Errhandler_get PMPI_Errhandler_get
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Comm_get_errhandler(comm, errhandler);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Errhandler_get(QMPI_Context context, int tool_id, MPI_Comm comm,
                        MPI_Errhandler *errhandler)
{
    return internal_Errhandler_get(comm, errhandler);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler)
{
    QMPI_Context context;
    QMPI_Errhandler_get_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Errhandler_get(context, 0, comm, errhandler);

    fn_ptr = (QMPI_Errhandler_get_t *) MPIR_QMPI_first_fn_ptrs[MPI_ERRHANDLER_GET_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_ERRHANDLER_GET_T], comm, errhandler);
}
#else /* ENABLE_QMPI */

int MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler)
{
    return internal_Errhandler_get(comm, errhandler);
}
#endif /* ENABLE_QMPI */
