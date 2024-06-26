/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Testall */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Testall = PMPI_Testall
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Testall  MPI_Testall
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Testall as PMPI_Testall
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                MPI_Status *array_of_statuses)  __attribute__ ((weak, alias("PMPI_Testall")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Testall
#define MPI_Testall PMPI_Testall
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Testall(int count, MPI_Request array_of_requests[], int *flag,
                            MPI_Status *array_of_statuses)
{
    int mpi_errno = MPI_SUCCESS;

    MPIR_ERRTEST_INITIALIZED_ORDIE();

    MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    MPIR_FUNC_TERSE_ENTER;

#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            if (count > 0) {
                MPIR_ERRTEST_ARGNULL(array_of_requests, "array_of_requests", mpi_errno);
                for (int i = 0; i < count; i++) {
                    MPIR_ERRTEST_ARRAYREQUEST_OR_NULL(array_of_requests[i], i, mpi_errno);
                }
            }
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_COUNT(count, mpi_errno);
            MPIR_ERRTEST_ARGNULL(flag, "flag", mpi_errno);
            if (count > 0) {
                MPIR_ERRTEST_ARGNULL(array_of_statuses, "array_of_statuses", mpi_errno);
            }
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ... */
    mpi_errno = MPIR_Testall(count, array_of_requests, flag, array_of_statuses);
    if (mpi_errno != MPI_SUCCESS)
        goto fn_fail;
    /* ... end of body of routine ... */

  fn_exit:
    MPIR_FUNC_TERSE_EXIT;
    MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    return mpi_errno;

  fn_fail:
    /* --BEGIN ERROR HANDLINE-- */
#ifdef HAVE_ERROR_CHECKING
    mpi_errno = MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE, __func__, __LINE__, MPI_ERR_OTHER,
                                     "**mpi_testall", "**mpi_testall %d %p %p %p", count,
                                     array_of_requests, flag, array_of_statuses);
#endif
    mpi_errno = MPIR_Err_return_comm(0, __func__, mpi_errno);
    /* --END ERROR HANDLING-- */
    goto fn_exit;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Testall(QMPI_Context context, int tool_id, int count, MPI_Request array_of_requests[],
                 int *flag, MPI_Status *array_of_statuses)
{
    return internal_Testall(count, array_of_requests, flag, array_of_statuses);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                MPI_Status *array_of_statuses)
{
    QMPI_Context context;
    QMPI_Testall_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Testall(context, 0, count, array_of_requests, flag, array_of_statuses);

    fn_ptr = (QMPI_Testall_t *) MPIR_QMPI_first_fn_ptrs[MPI_TESTALL_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_TESTALL_T], count, array_of_requests,
            flag, array_of_statuses);
}
#else /* ENABLE_QMPI */

int MPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                MPI_Status *array_of_statuses)
{
    return internal_Testall(count, array_of_requests, flag, array_of_statuses);
}
#endif /* ENABLE_QMPI */
