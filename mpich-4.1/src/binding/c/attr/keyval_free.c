/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Keyval_free */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Keyval_free = PMPI_Keyval_free
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Keyval_free  MPI_Keyval_free
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Keyval_free as PMPI_Keyval_free
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Keyval_free(int *keyval)  __attribute__ ((weak, alias("PMPI_Keyval_free")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Keyval_free
#define MPI_Keyval_free PMPI_Keyval_free
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Keyval_free(int *keyval)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Comm_free_keyval(keyval);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Keyval_free(QMPI_Context context, int tool_id, int *keyval)
{
    return internal_Keyval_free(keyval);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Keyval_free(int *keyval)
{
    QMPI_Context context;
    QMPI_Keyval_free_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Keyval_free(context, 0, keyval);

    fn_ptr = (QMPI_Keyval_free_t *) MPIR_QMPI_first_fn_ptrs[MPI_KEYVAL_FREE_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_KEYVAL_FREE_T], keyval);
}
#else /* ENABLE_QMPI */

int MPI_Keyval_free(int *keyval)
{
    return internal_Keyval_free(keyval);
}
#endif /* ENABLE_QMPI */
