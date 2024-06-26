/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Attr_delete */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Attr_delete = PMPI_Attr_delete
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Attr_delete  MPI_Attr_delete
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Attr_delete as PMPI_Attr_delete
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Attr_delete(MPI_Comm comm, int keyval)  __attribute__ ((weak, alias("PMPI_Attr_delete")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Attr_delete
#define MPI_Attr_delete PMPI_Attr_delete
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Attr_delete(MPI_Comm comm, int keyval)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Comm_delete_attr(comm, keyval);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Attr_delete(QMPI_Context context, int tool_id, MPI_Comm comm, int keyval)
{
    return internal_Attr_delete(comm, keyval);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Attr_delete(MPI_Comm comm, int keyval)
{
    QMPI_Context context;
    QMPI_Attr_delete_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Attr_delete(context, 0, comm, keyval);

    fn_ptr = (QMPI_Attr_delete_t *) MPIR_QMPI_first_fn_ptrs[MPI_ATTR_DELETE_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_ATTR_DELETE_T], comm, keyval);
}
#else /* ENABLE_QMPI */

int MPI_Attr_delete(MPI_Comm comm, int keyval)
{
    return internal_Attr_delete(comm, keyval);
}
#endif /* ENABLE_QMPI */
