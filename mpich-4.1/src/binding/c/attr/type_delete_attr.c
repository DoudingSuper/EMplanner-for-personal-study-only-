/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Type_delete_attr */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Type_delete_attr = PMPI_Type_delete_attr
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Type_delete_attr  MPI_Type_delete_attr
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Type_delete_attr as PMPI_Type_delete_attr
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval)
     __attribute__ ((weak, alias("PMPI_Type_delete_attr")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Type_delete_attr
#define MPI_Type_delete_attr PMPI_Type_delete_attr
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Type_delete_attr(MPI_Datatype datatype, int type_keyval)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Datatype *datatype_ptr ATTRIBUTE((unused)) = NULL;
    MPII_Keyval *type_keyval_ptr ATTRIBUTE((unused)) = NULL;

    MPIR_ERRTEST_INITIALIZED_ORDIE();

    MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    MPIR_FUNC_TERSE_ENTER;

#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_DATATYPE(datatype, "datatype", mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    MPIR_Datatype_get_ptr(datatype, datatype_ptr);
    MPII_Keyval_get_ptr(type_keyval, type_keyval_ptr);

#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_Datatype_valid_ptr(datatype_ptr, mpi_errno);
            if (mpi_errno) {
                goto fn_fail;
            }
            MPII_Keyval_valid_ptr(type_keyval_ptr, mpi_errno);
            if (mpi_errno) {
                goto fn_fail;
            }
            MPIR_ERRTEST_KEYVAL(type_keyval, MPIR_DATATYPE, "type_keyval", mpi_errno);
            MPIR_ERRTEST_KEYVAL_PERM(type_keyval, mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ... */
    mpi_errno = MPIR_Type_delete_attr_impl(datatype_ptr, type_keyval_ptr);
    if (mpi_errno) {
        goto fn_fail;
    }
    /* ... end of body of routine ... */

  fn_exit:
    MPIR_FUNC_TERSE_EXIT;
    MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    return mpi_errno;

  fn_fail:
    /* --BEGIN ERROR HANDLINE-- */
#ifdef HAVE_ERROR_CHECKING
    mpi_errno = MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE, __func__, __LINE__, MPI_ERR_OTHER,
                                     "**mpi_type_delete_attr", "**mpi_type_delete_attr %D %K", datatype,
                                     type_keyval);
#endif
    mpi_errno = MPIR_Err_return_comm(0, __func__, mpi_errno);
    /* --END ERROR HANDLING-- */
    goto fn_exit;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Type_delete_attr(QMPI_Context context, int tool_id, MPI_Datatype datatype,
                          int type_keyval)
{
    return internal_Type_delete_attr(datatype, type_keyval);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval)
{
    QMPI_Context context;
    QMPI_Type_delete_attr_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Type_delete_attr(context, 0, datatype, type_keyval);

    fn_ptr = (QMPI_Type_delete_attr_t *) MPIR_QMPI_first_fn_ptrs[MPI_TYPE_DELETE_ATTR_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_TYPE_DELETE_ATTR_T], datatype,
            type_keyval);
}
#else /* ENABLE_QMPI */

int MPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval)
{
    return internal_Type_delete_attr(datatype, type_keyval);
}
#endif /* ENABLE_QMPI */
