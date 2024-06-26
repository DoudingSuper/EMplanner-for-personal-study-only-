/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Errhandler_set */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Errhandler_set = PMPI_Errhandler_set
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Errhandler_set  MPI_Errhandler_set
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Errhandler_set as PMPI_Errhandler_set
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler)
     __attribute__ ((weak, alias("PMPI_Errhandler_set")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Errhandler_set
#define MPI_Errhandler_set PMPI_Errhandler_set
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Comm_set_errhandler(comm, errhandler);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Errhandler_set(QMPI_Context context, int tool_id, MPI_Comm comm,
                        MPI_Errhandler errhandler)
{
    return internal_Errhandler_set(comm, errhandler);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler)
{
    QMPI_Context context;
    QMPI_Errhandler_set_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Errhandler_set(context, 0, comm, errhandler);

    fn_ptr = (QMPI_Errhandler_set_t *) MPIR_QMPI_first_fn_ptrs[MPI_ERRHANDLER_SET_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_ERRHANDLER_SET_T], comm, errhandler);
}
#else /* ENABLE_QMPI */

int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler)
{
    return internal_Errhandler_set(comm, errhandler);
}
#endif /* ENABLE_QMPI */
