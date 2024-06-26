/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_T_event_get_num */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_T_event_get_num = PMPI_T_event_get_num
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_T_event_get_num  MPI_T_event_get_num
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_T_event_get_num as PMPI_T_event_get_num
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_T_event_get_num(int *num_events)  __attribute__ ((weak, alias("PMPI_T_event_get_num")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_T_event_get_num
#define MPI_T_event_get_num PMPI_T_event_get_num
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_T_event_get_num(int *num_events)
{
    int mpi_errno = MPI_SUCCESS;

    MPIT_ERRTEST_MPIT_INITIALIZED();

    MPIR_T_THREAD_CS_ENTER();
    MPIR_FUNC_TERSE_ENTER;

#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIT_ERRTEST_ARGNULL(num_events);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ... */
    mpi_errno = MPIR_T_event_get_num_impl(num_events);
    if (mpi_errno) {
        goto fn_fail;
    }
    /* ... end of body of routine ... */

  fn_exit:
    MPIR_FUNC_TERSE_EXIT;
    MPIR_T_THREAD_CS_EXIT();
    return mpi_errno;

  fn_fail:
    goto fn_exit;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_T_event_get_num(QMPI_Context context, int tool_id, int *num_events)
{
    return internal_T_event_get_num(num_events);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_T_event_get_num(int *num_events)
{
    QMPI_Context context;
    QMPI_T_event_get_num_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_T_event_get_num(context, 0, num_events);

    fn_ptr = (QMPI_T_event_get_num_t *) MPIR_QMPI_first_fn_ptrs[MPI_T_EVENT_GET_NUM_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_T_EVENT_GET_NUM_T], num_events);
}
#else /* ENABLE_QMPI */

int MPI_T_event_get_num(int *num_events)
{
    return internal_T_event_get_num(num_events);
}
#endif /* ENABLE_QMPI */
