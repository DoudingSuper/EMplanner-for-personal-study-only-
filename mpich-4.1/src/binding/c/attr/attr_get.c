/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Attr_get */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Attr_get = PMPI_Attr_get
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Attr_get  MPI_Attr_get
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Attr_get as PMPI_Attr_get
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag)
     __attribute__ ((weak, alias("PMPI_Attr_get")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Attr_get
#define MPI_Attr_get PMPI_Attr_get
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Comm_get_attr(comm, keyval, attribute_val, flag);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Attr_get(QMPI_Context context, int tool_id, MPI_Comm comm, int keyval, void *attribute_val,
                  int *flag)
{
    return internal_Attr_get(comm, keyval, attribute_val, flag);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag)
{
    QMPI_Context context;
    QMPI_Attr_get_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Attr_get(context, 0, comm, keyval, attribute_val, flag);

    fn_ptr = (QMPI_Attr_get_t *) MPIR_QMPI_first_fn_ptrs[MPI_ATTR_GET_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_ATTR_GET_T], comm, keyval,
            attribute_val, flag);
}
#else /* ENABLE_QMPI */

int MPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag)
{
    return internal_Attr_get(comm, keyval, attribute_val, flag);
}
#endif /* ENABLE_QMPI */
