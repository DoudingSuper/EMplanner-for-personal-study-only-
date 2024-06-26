/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* -- THIS FILE IS AUTO-GENERATED -- */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Neighbor_alltoall_init */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Neighbor_alltoall_init = PMPI_Neighbor_alltoall_init
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Neighbor_alltoall_init  MPI_Neighbor_alltoall_init
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Neighbor_alltoall_init as PMPI_Neighbor_alltoall_init
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Neighbor_alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                               void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Info info, MPI_Request *request)
                                __attribute__ ((weak, alias("PMPI_Neighbor_alltoall_init")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Neighbor_alltoall_init
#define MPI_Neighbor_alltoall_init PMPI_Neighbor_alltoall_init
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Neighbor_alltoall_init(const void *sendbuf, int sendcount,
                                           MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                           MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                           MPI_Request *request)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Comm *comm_ptr ATTRIBUTE((unused)) = NULL;
    MPIR_Info *info_ptr ATTRIBUTE((unused)) = NULL;

    MPIR_ERRTEST_INITIALIZED_ORDIE();

    MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    MPIR_FUNC_TERSE_ENTER;

#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_COMM(comm, mpi_errno);
            MPIR_ERRTEST_INFO_OR_NULL(info, mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    MPIR_Comm_get_ptr(comm, comm_ptr);
    if (info != MPI_INFO_NULL) {
        MPIR_Info_get_ptr(info, info_ptr);
    }

#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_Comm_valid_ptr(comm_ptr, mpi_errno, FALSE);
            if (mpi_errno) {
                goto fn_fail;
            }
            if (info != MPI_INFO_NULL) {
                MPIR_Info_valid_ptr(info_ptr, mpi_errno);
                if (mpi_errno) {
                    goto fn_fail;
                }
            }
            MPIR_ERRTEST_COUNT(sendcount, mpi_errno);
            if (sendcount > 0) {
                MPIR_ERRTEST_DATATYPE(sendtype, "datatype", mpi_errno);
                if (!HANDLE_IS_BUILTIN(sendtype)) {
                    MPIR_Datatype *datatype_ptr ATTRIBUTE((unused)) = NULL;
                    MPIR_Datatype_get_ptr(sendtype, datatype_ptr);
                    MPIR_Datatype_valid_ptr(datatype_ptr, mpi_errno);
                    if (mpi_errno) {
                        goto fn_fail;
                    }
                    MPIR_Datatype_committed_ptr(datatype_ptr, mpi_errno);
                    if (mpi_errno) {
                        goto fn_fail;
                    }
                }
                MPIR_ERRTEST_USERBUFFER(sendbuf, sendcount, sendtype, mpi_errno);
            }
            MPIR_ERRTEST_COUNT(recvcount, mpi_errno);
            if (recvcount > 0) {
                MPIR_ERRTEST_DATATYPE(recvtype, "datatype", mpi_errno);
                if (!HANDLE_IS_BUILTIN(recvtype)) {
                    MPIR_Datatype *datatype_ptr ATTRIBUTE((unused)) = NULL;
                    MPIR_Datatype_get_ptr(recvtype, datatype_ptr);
                    MPIR_Datatype_valid_ptr(datatype_ptr, mpi_errno);
                    if (mpi_errno) {
                        goto fn_fail;
                    }
                    MPIR_Datatype_committed_ptr(datatype_ptr, mpi_errno);
                    if (mpi_errno) {
                        goto fn_fail;
                    }
                }
                MPIR_ERRTEST_USERBUFFER(recvbuf, recvcount, recvtype, mpi_errno);
            }
            MPIR_ERRTEST_ARGNULL(request, "request", mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ... */
    MPIR_Request *request_ptr = NULL;
    mpi_errno = MPIR_Neighbor_alltoall_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype,
                                            comm_ptr, info_ptr, &request_ptr);
    if (mpi_errno) {
        goto fn_fail;
    }
    if (!request_ptr) {
        request_ptr = MPIR_Request_create_complete(MPIR_REQUEST_KIND__COLL);
    }
    *request = request_ptr->handle;
    /* ... end of body of routine ... */

  fn_exit:
    MPIR_FUNC_TERSE_EXIT;
    MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    return mpi_errno;

  fn_fail:
    /* --BEGIN ERROR HANDLINE-- */
#ifdef HAVE_ERROR_CHECKING
    mpi_errno = MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE, __func__, __LINE__, MPI_ERR_OTHER,
                                     "**mpi_neighbor_alltoall_init",
                                     "**mpi_neighbor_alltoall_init %p %d %D %p %d %D %C %I %p", sendbuf,
                                     sendcount, sendtype, recvbuf, recvcount, recvtype, comm, info,
                                     request);
#endif
    mpi_errno = MPIR_Err_return_comm(comm_ptr, __func__, mpi_errno);
    /* --END ERROR HANDLING-- */
    goto fn_exit;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Neighbor_alltoall_init(QMPI_Context context, int tool_id, const void *sendbuf,
                                int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                MPI_Datatype recvtype, MPI_Comm comm, MPI_Info info,
                                MPI_Request *request)
{
    return internal_Neighbor_alltoall_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, info, request);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Neighbor_alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                               void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Info info, MPI_Request *request)
{
    QMPI_Context context;
    QMPI_Neighbor_alltoall_init_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Neighbor_alltoall_init(context, 0, sendbuf, sendcount, sendtype, recvbuf,
                                           recvcount, recvtype, comm, info, request);

    fn_ptr = (QMPI_Neighbor_alltoall_init_t *) MPIR_QMPI_first_fn_ptrs[MPI_NEIGHBOR_ALLTOALL_INIT_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_NEIGHBOR_ALLTOALL_INIT_T], sendbuf,
            sendcount, sendtype, recvbuf, recvcount, recvtype, comm, info, request);
}
#else /* ENABLE_QMPI */

int MPI_Neighbor_alltoall_init(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                               void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                               MPI_Info info, MPI_Request *request)
{
    return internal_Neighbor_alltoall_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, info, request);
}
#endif /* ENABLE_QMPI */

/* -- Begin Profiling Symbol Block for routine MPI_Neighbor_alltoall_init_c */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Neighbor_alltoall_init_c = PMPI_Neighbor_alltoall_init_c
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Neighbor_alltoall_init_c  MPI_Neighbor_alltoall_init_c
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Neighbor_alltoall_init_c as PMPI_Neighbor_alltoall_init_c
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Neighbor_alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                 MPI_Comm comm, MPI_Info info, MPI_Request *request)
                                  __attribute__ ((weak, alias("PMPI_Neighbor_alltoall_init_c")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Neighbor_alltoall_init_c
#define MPI_Neighbor_alltoall_init_c PMPI_Neighbor_alltoall_init_c
#endif /* MPICH_MPI_FROM_PMPI */

static int internal_Neighbor_alltoall_init_c(const void *sendbuf, MPI_Count sendcount,
                                             MPI_Datatype sendtype, void *recvbuf,
                                             MPI_Count recvcount, MPI_Datatype recvtype,
                                             MPI_Comm comm, MPI_Info info, MPI_Request *request)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Comm *comm_ptr ATTRIBUTE((unused)) = NULL;
    MPIR_Info *info_ptr ATTRIBUTE((unused)) = NULL;

    MPIR_ERRTEST_INITIALIZED_ORDIE();

    MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    MPIR_FUNC_TERSE_ENTER;

#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_COMM(comm, mpi_errno);
            MPIR_ERRTEST_INFO_OR_NULL(info, mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    MPIR_Comm_get_ptr(comm, comm_ptr);
    if (info != MPI_INFO_NULL) {
        MPIR_Info_get_ptr(info, info_ptr);
    }

#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_Comm_valid_ptr(comm_ptr, mpi_errno, FALSE);
            if (mpi_errno) {
                goto fn_fail;
            }
            if (info != MPI_INFO_NULL) {
                MPIR_Info_valid_ptr(info_ptr, mpi_errno);
                if (mpi_errno) {
                    goto fn_fail;
                }
            }
            MPIR_ERRTEST_COUNT(sendcount, mpi_errno);
            if (sendcount > 0) {
                MPIR_ERRTEST_DATATYPE(sendtype, "datatype", mpi_errno);
                if (!HANDLE_IS_BUILTIN(sendtype)) {
                    MPIR_Datatype *datatype_ptr ATTRIBUTE((unused)) = NULL;
                    MPIR_Datatype_get_ptr(sendtype, datatype_ptr);
                    MPIR_Datatype_valid_ptr(datatype_ptr, mpi_errno);
                    if (mpi_errno) {
                        goto fn_fail;
                    }
                    MPIR_Datatype_committed_ptr(datatype_ptr, mpi_errno);
                    if (mpi_errno) {
                        goto fn_fail;
                    }
                }
                MPIR_ERRTEST_USERBUFFER(sendbuf, sendcount, sendtype, mpi_errno);
            }
            MPIR_ERRTEST_COUNT(recvcount, mpi_errno);
            if (recvcount > 0) {
                MPIR_ERRTEST_DATATYPE(recvtype, "datatype", mpi_errno);
                if (!HANDLE_IS_BUILTIN(recvtype)) {
                    MPIR_Datatype *datatype_ptr ATTRIBUTE((unused)) = NULL;
                    MPIR_Datatype_get_ptr(recvtype, datatype_ptr);
                    MPIR_Datatype_valid_ptr(datatype_ptr, mpi_errno);
                    if (mpi_errno) {
                        goto fn_fail;
                    }
                    MPIR_Datatype_committed_ptr(datatype_ptr, mpi_errno);
                    if (mpi_errno) {
                        goto fn_fail;
                    }
                }
                MPIR_ERRTEST_USERBUFFER(recvbuf, recvcount, recvtype, mpi_errno);
            }
            MPIR_ERRTEST_ARGNULL(request, "request", mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ... */
    if (sizeof(MPI_Count) == sizeof(MPI_Aint)) {
        MPIR_Request *request_ptr = NULL;
        mpi_errno = MPIR_Neighbor_alltoall_init(sendbuf, (MPI_Aint) sendcount, sendtype, recvbuf,
                                                (MPI_Aint) recvcount, recvtype, comm_ptr, info_ptr,
                                                &request_ptr);
        if (mpi_errno) {
            goto fn_fail;
        }
        if (!request_ptr) {
            request_ptr = MPIR_Request_create_complete(MPIR_REQUEST_KIND__COLL);
        }
        *request = request_ptr->handle;
    } else {
        /* MPI_Count is bigger than MPI_Aint */
        if (sendcount > MPIR_AINT_MAX) {
            mpi_errno = MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE,
                                             __func__, __LINE__, MPI_ERR_OTHER,
                                             "**too_big_for_input",
                                             "**too_big_for_input %s", "sendcount");
            goto fn_fail;
        }
        if (recvcount > MPIR_AINT_MAX) {
            mpi_errno = MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE,
                                             __func__, __LINE__, MPI_ERR_OTHER,
                                             "**too_big_for_input",
                                             "**too_big_for_input %s", "recvcount");
            goto fn_fail;
        }
        MPIR_Request *request_ptr = NULL;
        mpi_errno = MPIR_Neighbor_alltoall_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype,
                                                comm_ptr, info_ptr, &request_ptr);
        if (mpi_errno) {
            goto fn_fail;
        }
        if (!request_ptr) {
            request_ptr = MPIR_Request_create_complete(MPIR_REQUEST_KIND__COLL);
        }
        *request = request_ptr->handle;
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
                                     "**mpi_neighbor_alltoall_init_c",
                                     "**mpi_neighbor_alltoall_init_c %p %c %D %p %c %D %C %I %p",
                                     sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm,
                                     info, request);
#endif
    mpi_errno = MPIR_Err_return_comm(comm_ptr, __func__, mpi_errno);
    /* --END ERROR HANDLING-- */
    goto fn_exit;
}

#ifdef ENABLE_QMPI
#ifndef MPICH_MPI_FROM_PMPI
int QMPI_Neighbor_alltoall_init_c(QMPI_Context context, int tool_id, const void *sendbuf,
                                  MPI_Count sendcount, MPI_Datatype sendtype, void *recvbuf,
                                  MPI_Count recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                                  MPI_Info info, MPI_Request *request)
{
    return internal_Neighbor_alltoall_init_c(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, info, request);
}
#endif /* MPICH_MPI_FROM_PMPI */
int MPI_Neighbor_alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                 MPI_Comm comm, MPI_Info info, MPI_Request *request)
{
    QMPI_Context context;
    QMPI_Neighbor_alltoall_init_c_t *fn_ptr;

    context.storage_stack = NULL;

    if (MPIR_QMPI_num_tools == 0)
        return QMPI_Neighbor_alltoall_init_c(context, 0, sendbuf, sendcount, sendtype, recvbuf,
                                             recvcount, recvtype, comm, info, request);

    fn_ptr = (QMPI_Neighbor_alltoall_init_c_t *) MPIR_QMPI_first_fn_ptrs[MPI_NEIGHBOR_ALLTOALL_INIT_C_T];

    return (*fn_ptr) (context, MPIR_QMPI_first_tool_ids[MPI_NEIGHBOR_ALLTOALL_INIT_C_T], sendbuf,
            sendcount, sendtype, recvbuf, recvcount, recvtype, comm, info, request);
}
#else /* ENABLE_QMPI */

int MPI_Neighbor_alltoall_init_c(const void *sendbuf, MPI_Count sendcount, MPI_Datatype sendtype,
                                 void *recvbuf, MPI_Count recvcount, MPI_Datatype recvtype,
                                 MPI_Comm comm, MPI_Info info, MPI_Request *request)
{
    return internal_Neighbor_alltoall_init_c(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, info, request);
}
#endif /* ENABLE_QMPI */
