/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Type_hvector */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Type_hvector = PMPI_Type_hvector
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Type_hvector  MPI_Type_hvector
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Type_hvector as PMPI_Type_hvector
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype)  __attribute__ ((weak, alias("PMPI_Type_hvector")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Type_hvector
#define MPI_Type_hvector PMPI_Type_hvector
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                                 MPI_Datatype *newtype)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Type_hvector(QMPI_Context context, int tool_id, int count, int blocklength,
                      MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
{
    return internal_Type_hvector(count, blocklength, stride, oldtype, newtype);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype)
{
    QMPI_Context context;
    QMPI_Type_hvector_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Type_hvector(context, 0, count, blocklength, stride, oldtype, newtype);

    fn_ptr = (QMPI_Type_hvector_t *) MPIR_QMPI_first_fn_ptrs[MPI_TYPE_HVECTOR_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_TYPE_HVECTOR_T], count, blocklength,
            stride, oldtype, newtype);
}
#else /* ENABLE_QMPI */

int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype)
{
    return internal_Type_hvector(count, blocklength, stride, oldtype, newtype);
}
#endif /* ENABLE_QMPI */
