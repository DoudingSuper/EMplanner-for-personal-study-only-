/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Keyval_create */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Keyval_create = PMPI_Keyval_create
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Keyval_create  MPI_Keyval_create
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Keyval_create as PMPI_Keyval_create
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn, int *keyval,
                      void *extra_state)  __attribute__ ((weak, alias("PMPI_Keyval_create")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Keyval_create
#define MPI_Keyval_create PMPI_Keyval_create
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn,
                                  int *keyval, void *extra_state)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Comm_create_keyval(copy_fn, delete_fn, keyval, extra_state);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Keyval_create(QMPI_Context context, int tool_id, MPI_Copy_function *copy_fn,
                       MPI_Delete_function *delete_fn, int *keyval, void *extra_state)
{
    return internal_Keyval_create(copy_fn, delete_fn, keyval, extra_state);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn, int *keyval,
                      void *extra_state)
{
    QMPI_Context context;
    QMPI_Keyval_create_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Keyval_create(context, 0, copy_fn, delete_fn, keyval, extra_state);

    fn_ptr = (QMPI_Keyval_create_t *) MPIR_QMPI_first_fn_ptrs[MPI_KEYVAL_CREATE_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_KEYVAL_CREATE_T], copy_fn, delete_fn,
            keyval, extra_state);
}
#else /* ENABLE_QMPI */

int MPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn, int *keyval,
                      void *extra_state)
{
    return internal_Keyval_create(copy_fn, delete_fn, keyval, extra_state);
}
#endif /* ENABLE_QMPI */
