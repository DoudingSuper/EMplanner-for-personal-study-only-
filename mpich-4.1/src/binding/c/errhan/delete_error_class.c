/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPIX_Delete_error_class */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPIX_Delete_error_class = PMPIX_Delete_error_class
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPIX_Delete_error_class  MPIX_Delete_error_class
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPIX_Delete_error_class as PMPIX_Delete_error_class
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPIX_Delete_error_class(int errorclass)
     __attribute__ ((weak, alias("PMPIX_Delete_error_class")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPIX_Delete_error_class
#define MPIX_Delete_error_class PMPIX_Delete_error_class
#endif /* MPICH_MPI_FROM_PMPI */

static int internalX_Delete_error_class(int errorclass)
{
    int mpi_errno = MPI_SUCCESS;

    MPIR_ERRTEST_INITIALIZED_ORDIE();

    MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    MPIR_FUNC_TERSE_ENTER;

    /* ... body of routine ... */
    mpi_errno = MPIR_Delete_error_class_impl(errorclass);
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
                                     "**mpix_delete_error_class", "**mpix_delete_error_class %d",
                                     errorclass);
#endif
    mpi_errno = MPIR_Err_return_comm(0, __func__, mpi_errno);
    /* --END ERROR HANDLING-- */
    goto fn_exit;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPIX_Delete_error_class(QMPI_Context context, int tool_id, int errorclass)
{
    return internalX_Delete_error_class(errorclass);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPIX_Delete_error_class(int errorclass)
{
    QMPI_Context context;
    QMPIX_Delete_error_class_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPIX_Delete_error_class(context, 0, errorclass);

    fn_ptr = (QMPIX_Delete_error_class_t *) MPIR_QMPI_first_fn_ptrs[MPIX_DELETE_ERROR_CLASS_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPIX_DELETE_ERROR_CLASS_T], errorclass);
}
#else /* ENABLE_QMPI */

int MPIX_Delete_error_class(int errorclass)
{
    return internalX_Delete_error_class(errorclass);
}
#endif /* ENABLE_QMPI */
