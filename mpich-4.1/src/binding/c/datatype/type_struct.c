/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Type_struct */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Type_struct = PMPI_Type_struct
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Type_struct  MPI_Type_struct
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Type_struct as PMPI_Type_struct
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Type_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                    MPI_Datatype array_of_types[], MPI_Datatype *newtype)
                     __attribute__ ((weak, alias("PMPI_Type_struct")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Type_struct
#define MPI_Type_struct PMPI_Type_struct
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Type_struct(int count, int array_of_blocklengths[],
                                MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[],
                                MPI_Datatype *newtype)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_TERSE_ENTER;

    mpi_errno = PMPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);

    MPIR_FUNC_TERSE_EXIT;
    return mpi_errno;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Type_struct(QMPI_Context context, int tool_id, int count, int array_of_blocklengths[],
                     MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[],
                     MPI_Datatype *newtype)
{
    return internal_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Type_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                    MPI_Datatype array_of_types[], MPI_Datatype *newtype)
{
    QMPI_Context context;
    QMPI_Type_struct_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Type_struct(context, 0, count, array_of_blocklengths, array_of_displacements,
                                array_of_types, newtype);

    fn_ptr = (QMPI_Type_struct_t *) MPIR_QMPI_first_fn_ptrs[MPI_TYPE_STRUCT_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_TYPE_STRUCT_T], count,
            array_of_blocklengths, array_of_displacements, array_of_types, newtype);
}
#else /* ENABLE_QMPI */

int MPI_Type_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[],
                    MPI_Datatype array_of_types[], MPI_Datatype *newtype)
{
    return internal_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
}
#endif /* ENABLE_QMPI */
