/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Errhandler_create */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Errhandler_create = PMPI_Errhandler_create
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Errhandler_create  MPI_Errhandler_create
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Errhandler_create as PMPI_Errhandler_create
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Errhandler_create(MPI_Comm_errhandler_function *comm_errhandler_fn,
                          MPI_Errhandler *errhandler)
                           __attribute__ ((weak, alias("PMPI_Errhandler_create")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Errhandler_create
#define MPI_Errhandler_create PMPI_Errhandler_create
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Errhandler_create(MPI_Comm_errhandler_function *comm_errhandler_fn,
                                      MPI_Errhandler *errhandler)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Comm_create_errhandler(comm_errhandler_fn, errhandler);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Errhandler_create(QMPI_Context context, int tool_id,
                           MPI_Comm_errhandler_function *comm_errhandler_fn,
                           MPI_Errhandler *errhandler)
{
    return internal_Errhandler_create(comm_errhandler_fn, errhandler);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Errhandler_create(MPI_Comm_errhandler_function *comm_errhandler_fn,
                          MPI_Errhandler *errhandler)
{
    QMPI_Context context;
    QMPI_Errhandler_create_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Errhandler_create(context, 0, comm_errhandler_fn, errhandler);

    fn_ptr = (QMPI_Errhandler_create_t *) MPIR_QMPI_first_fn_ptrs[MPI_ERRHANDLER_CREATE_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_ERRHANDLER_CREATE_T],
            comm_errhandler_fn, errhandler);
}
#else /* ENABLE_QMPI */

int MPI_Errhandler_create(MPI_Comm_errhandler_function *comm_errhandler_fn,
                          MPI_Errhandler *errhandler)
{
    return internal_Errhandler_create(comm_errhandler_fn, errhandler);
}
#endif /* ENABLE_QMPI */
