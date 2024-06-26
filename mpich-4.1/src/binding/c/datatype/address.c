/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Address */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Address = PMPI_Address
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Address  MPI_Address
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Address as PMPI_Address
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Address(void *location, MPI_Aint *address)  __attribute__ ((weak, alias("PMPI_Address")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Address
#define MPI_Address PMPI_Address
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Address(void *location, MPI_Aint *address)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Get_address(location, address);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Address(QMPI_Context context, int tool_id, void *location, MPI_Aint *address)
{
    return internal_Address(location, address);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Address(void *location, MPI_Aint *address)
{
    QMPI_Context context;
    QMPI_Address_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Address(context, 0, location, address);

    fn_ptr = (QMPI_Address_t *) MPIR_QMPI_first_fn_ptrs[MPI_ADDRESS_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_ADDRESS_T], location, address);
}
#else /* ENABLE_QMPI */

int MPI_Address(void *location, MPI_Aint *address)
{
    return internal_Address(location, address);
}
#endif /* ENABLE_QMPI */
