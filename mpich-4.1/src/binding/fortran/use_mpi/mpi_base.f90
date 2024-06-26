!
! Copyright (C) by Argonne National Laboratory
!     See COPYRIGHT in top-level directory
!

! -- THIS FILE IS AUTO-GENERATED -- 

MODULE mpi_base
    IMPLICIT NONE
    INTERFACE
    
    SUBROUTINE MPI_Comm_create_keyval(comm_copy_attr_fn, comm_delete_attr_fn, comm_keyval, extra_state, &
                                      ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        EXTERNAL :: comm_copy_attr_fn
        EXTERNAL :: comm_delete_attr_fn
        INTEGER :: comm_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_create_keyval
    
    SUBROUTINE MPI_Comm_delete_attr(comm, comm_keyval, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: comm_keyval
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_delete_attr
    
    SUBROUTINE MPI_Comm_free_keyval(comm_keyval, ierr)
        IMPLICIT NONE
        INTEGER :: comm_keyval
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_free_keyval
    
    SUBROUTINE MPI_Comm_get_attr(comm, comm_keyval, attribute_val, flag, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: comm_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_get_attr
    
    SUBROUTINE MPI_Comm_set_attr(comm, comm_keyval, attribute_val, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: comm_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_set_attr
    
    SUBROUTINE MPI_Keyval_create(copy_fn, delete_fn, keyval, extra_state, ierr)
        IMPLICIT NONE
        EXTERNAL :: copy_fn
        EXTERNAL :: delete_fn
        INTEGER :: keyval
        INTEGER :: extra_state
        INTEGER :: ierr
    END SUBROUTINE MPI_Keyval_create
    
    SUBROUTINE MPI_Keyval_free(keyval, ierr)
        IMPLICIT NONE
        INTEGER :: keyval
        INTEGER :: ierr
    END SUBROUTINE MPI_Keyval_free
    
    SUBROUTINE MPI_Type_create_keyval(type_copy_attr_fn, type_delete_attr_fn, type_keyval, extra_state, &
                                      ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        EXTERNAL :: type_copy_attr_fn
        EXTERNAL :: type_delete_attr_fn
        INTEGER :: type_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_keyval
    
    SUBROUTINE MPI_Type_delete_attr(datatype, type_keyval, ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: type_keyval
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_delete_attr
    
    SUBROUTINE MPI_Type_free_keyval(type_keyval, ierr)
        IMPLICIT NONE
        INTEGER :: type_keyval
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_free_keyval
    
    SUBROUTINE MPI_Type_get_attr(datatype, type_keyval, attribute_val, flag, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: type_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_get_attr
    
    SUBROUTINE MPI_Type_set_attr(datatype, type_keyval, attribute_val, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: type_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_set_attr
    
    SUBROUTINE MPI_Win_create_keyval(win_copy_attr_fn, win_delete_attr_fn, win_keyval, extra_state, &
                                     ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        EXTERNAL :: win_copy_attr_fn
        EXTERNAL :: win_delete_attr_fn
        INTEGER :: win_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_create_keyval
    
    SUBROUTINE MPI_Win_delete_attr(win, win_keyval, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: win_keyval
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_delete_attr
    
    SUBROUTINE MPI_Win_free_keyval(win_keyval, ierr)
        IMPLICIT NONE
        INTEGER :: win_keyval
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_free_keyval
    
    SUBROUTINE MPI_Win_get_attr(win, win_keyval, attribute_val, flag, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: win_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_get_attr
    
    SUBROUTINE MPI_Win_set_attr(win, win_keyval, attribute_val, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: win_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_set_attr
    
    SUBROUTINE MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Allgather
    
    SUBROUTINE MPI_Allgather_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, &
                                  info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Allgather_init
    
    SUBROUTINE MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm, &
                              ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: displs(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Allgatherv
    
    SUBROUTINE MPI_Allgatherv_init(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, &
                                   comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: displs(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Allgatherv_init
    
    SUBROUTINE MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Allreduce
    
    SUBROUTINE MPI_Allreduce_init(sendbuf, recvbuf, count, datatype, op, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Allreduce_init
    
    SUBROUTINE MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Alltoall
    
    SUBROUTINE MPI_Alltoall_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, info, &
                                 request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Alltoall_init
    
    SUBROUTINE MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, &
                             recvtype, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: sdispls(*)
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: rdispls(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Alltoallv
    
    SUBROUTINE MPI_Alltoallv_init(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, &
                                  recvtype, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: sdispls(*)
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: rdispls(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Alltoallv_init
    
    SUBROUTINE MPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, &
                             recvtypes, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: sdispls(*)
        INTEGER :: sendtypes(*)
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: rdispls(*)
        INTEGER :: recvtypes(*)
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Alltoallw
    
    SUBROUTINE MPI_Alltoallw_init(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, &
                                  recvtypes, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: sdispls(*)
        INTEGER :: sendtypes(*)
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: rdispls(*)
        INTEGER :: recvtypes(*)
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Alltoallw_init
    
    SUBROUTINE MPI_Barrier(comm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Barrier
    
    SUBROUTINE MPI_Barrier_init(comm, info, request, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Barrier_init
    
    SUBROUTINE MPI_Bcast(buffer, count, datatype, root, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buffer
        REAL :: buffer
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Bcast
    
    SUBROUTINE MPI_Bcast_init(buffer, count, datatype, root, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buffer
        REAL :: buffer
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Bcast_init
    
    SUBROUTINE MPI_Exscan(sendbuf, recvbuf, count, datatype, op, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Exscan
    
    SUBROUTINE MPI_Exscan_init(sendbuf, recvbuf, count, datatype, op, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Exscan_init
    
    SUBROUTINE MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Gather
    
    SUBROUTINE MPI_Gather_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, &
                               info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Gather_init
    
    SUBROUTINE MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, &
                           comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: displs(*)
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Gatherv
    
    SUBROUTINE MPI_Gatherv_init(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, &
                                root, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: displs(*)
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Gatherv_init
    
    SUBROUTINE MPI_Iallgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request, &
                              ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Iallgather
    
    SUBROUTINE MPI_Iallgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, &
                               comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: displs(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Iallgatherv
    
    SUBROUTINE MPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Iallreduce
    
    SUBROUTINE MPI_Ialltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, request, &
                             ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ialltoall
    
    SUBROUTINE MPI_Ialltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, &
                              recvtype, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: sdispls(*)
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: rdispls(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ialltoallv
    
    SUBROUTINE MPI_Ialltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, &
                              recvtypes, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: sdispls(*)
        INTEGER :: sendtypes(*)
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: rdispls(*)
        INTEGER :: recvtypes(*)
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ialltoallw
    
    SUBROUTINE MPI_Ibarrier(comm, request, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ibarrier
    
    SUBROUTINE MPI_Ibcast(buffer, count, datatype, root, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buffer
        REAL :: buffer
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ibcast
    
    SUBROUTINE MPI_Iexscan(sendbuf, recvbuf, count, datatype, op, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Iexscan
    
    SUBROUTINE MPI_Igather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, &
                           request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Igather
    
    SUBROUTINE MPI_Igatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, &
                            comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: displs(*)
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Igatherv
    
    SUBROUTINE MPI_Ineighbor_allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, &
                                       request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ineighbor_allgather
    
    SUBROUTINE MPI_Ineighbor_allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, &
                                        recvtype, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: displs(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ineighbor_allgatherv
    
    SUBROUTINE MPI_Ineighbor_alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, &
                                      request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ineighbor_alltoall
    
    SUBROUTINE MPI_Ineighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, &
                                       rdispls, recvtype, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: sdispls(*)
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: rdispls(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ineighbor_alltoallv
    
    SUBROUTINE MPI_Ineighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, &
                                       rdispls, recvtypes, comm, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: sdispls(*)
        INTEGER :: sendtypes(*)
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: rdispls(*)
        INTEGER :: recvtypes(*)
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ineighbor_alltoallw
    
    SUBROUTINE MPI_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ireduce
    
    SUBROUTINE MPI_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ireduce_scatter
    
    SUBROUTINE MPI_Ireduce_scatter_block(sendbuf, recvbuf, recvcount, datatype, op, comm, request, &
                                         ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ireduce_scatter_block
    
    SUBROUTINE MPI_Iscan(sendbuf, recvbuf, count, datatype, op, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Iscan
    
    SUBROUTINE MPI_Iscatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, &
                            request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Iscatter
    
    SUBROUTINE MPI_Iscatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, &
                             comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: displs(*)
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Iscatterv
    
    SUBROUTINE MPI_Neighbor_allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, &
                                      ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_allgather
    
    SUBROUTINE MPI_Neighbor_allgather_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, &
                                           comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_allgather_init
    
    SUBROUTINE MPI_Neighbor_allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, &
                                       recvtype, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: displs(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_allgatherv
    
    SUBROUTINE MPI_Neighbor_allgatherv_init(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, &
                                            recvtype, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: displs(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_allgatherv_init
    
    SUBROUTINE MPI_Neighbor_alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, &
                                     ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_alltoall
    
    SUBROUTINE MPI_Neighbor_alltoall_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, &
                                          comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_alltoall_init
    
    SUBROUTINE MPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, &
                                      rdispls, recvtype, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: sdispls(*)
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: rdispls(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_alltoallv
    
    SUBROUTINE MPI_Neighbor_alltoallv_init(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, &
                                           rdispls, recvtype, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: sdispls(*)
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: rdispls(*)
        INTEGER :: recvtype
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_alltoallv_init
    
    SUBROUTINE MPI_Neighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, &
                                      rdispls, recvtypes, comm, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: sdispls(*)
        INTEGER :: sendtypes(*)
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: rdispls(*)
        INTEGER :: recvtypes(*)
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_alltoallw
    
    SUBROUTINE MPI_Neighbor_alltoallw_init(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, &
                                           rdispls, recvtypes, comm, info, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: sdispls(*)
        INTEGER :: sendtypes(*)
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: rdispls(*)
        INTEGER :: recvtypes(*)
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Neighbor_alltoallw_init
    
    SUBROUTINE MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Reduce
    
    SUBROUTINE MPI_Reduce_init(sendbuf, recvbuf, count, datatype, op, root, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Reduce_init
    
    SUBROUTINE MPI_Reduce_local(inbuf, inoutbuf, count, datatype, op, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: inbuf, inoutbuf
        REAL :: inbuf
        REAL :: inoutbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: ierr
    END SUBROUTINE MPI_Reduce_local
    
    SUBROUTINE MPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Reduce_scatter
    
    SUBROUTINE MPI_Reduce_scatter_block(sendbuf, recvbuf, recvcount, datatype, op, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Reduce_scatter_block
    
    SUBROUTINE MPI_Reduce_scatter_block_init(sendbuf, recvbuf, recvcount, datatype, op, comm, info, &
                                             request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Reduce_scatter_block_init
    
    SUBROUTINE MPI_Reduce_scatter_init(sendbuf, recvbuf, recvcounts, datatype, op, comm, info, request, &
                                       ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: recvcounts(*)
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Reduce_scatter_init
    
    SUBROUTINE MPI_Scan(sendbuf, recvbuf, count, datatype, op, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Scan
    
    SUBROUTINE MPI_Scan_init(sendbuf, recvbuf, count, datatype, op, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Scan_init
    
    SUBROUTINE MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, &
                           ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Scatter
    
    SUBROUTINE MPI_Scatter_init(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, &
                                info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Scatter_init
    
    SUBROUTINE MPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, &
                            comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: displs(*)
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Scatterv
    
    SUBROUTINE MPI_Scatterv_init(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, &
                                 root, comm, info, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcounts(*)
        INTEGER :: displs(*)
        INTEGER :: sendtype
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Scatterv_init
    
    SUBROUTINE MPI_Comm_compare(comm1, comm2, result, ierr)
        IMPLICIT NONE
        INTEGER :: comm1
        INTEGER :: comm2
        INTEGER :: result
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_compare
    
    SUBROUTINE MPI_Comm_create(comm, group, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: group
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_create
    
    SUBROUTINE MPI_Comm_create_group(comm, group, tag, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: group
        INTEGER :: tag
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_create_group
    
    SUBROUTINE MPI_Comm_dup(comm, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_dup
    
    SUBROUTINE MPI_Comm_dup_with_info(comm, info, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_dup_with_info
    
    SUBROUTINE MPI_Comm_free(comm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_free
    
    SUBROUTINE MPI_Comm_get_info(comm, info_used, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: info_used
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_get_info
    
    SUBROUTINE MPI_Comm_get_name(comm, comm_name, resultlen, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        CHARACTER*(*) :: comm_name
        INTEGER :: resultlen
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_get_name
    
    SUBROUTINE MPI_Comm_group(comm, group, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: group
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_group
    
    SUBROUTINE MPI_Comm_idup(comm, newcomm, request, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: newcomm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_idup
    
    SUBROUTINE MPI_Comm_idup_with_info(comm, info, newcomm, request, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: newcomm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_idup_with_info
    
    SUBROUTINE MPI_Comm_rank(comm, rank, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: rank
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_rank
    
    SUBROUTINE MPI_Comm_remote_group(comm, group, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: group
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_remote_group
    
    SUBROUTINE MPI_Comm_remote_size(comm, size, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_remote_size
    
    SUBROUTINE MPI_Comm_set_info(comm, info, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_set_info
    
    SUBROUTINE MPI_Comm_set_name(comm, comm_name, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        CHARACTER*(*) :: comm_name
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_set_name
    
    SUBROUTINE MPI_Comm_size(comm, size, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_size
    
    SUBROUTINE MPI_Comm_split(comm, color, key, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: color
        INTEGER :: key
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_split
    
    SUBROUTINE MPI_Comm_split_type(comm, split_type, key, info, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: split_type
        INTEGER :: key
        INTEGER :: info
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_split_type
    
    SUBROUTINE MPI_Comm_test_inter(comm, flag, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_test_inter
    
    SUBROUTINE MPI_Intercomm_create(local_comm, local_leader, peer_comm, remote_leader, tag, &
                                    newintercomm, ierr)
        IMPLICIT NONE
        INTEGER :: local_comm
        INTEGER :: local_leader
        INTEGER :: peer_comm
        INTEGER :: remote_leader
        INTEGER :: tag
        INTEGER :: newintercomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Intercomm_create
    
    SUBROUTINE MPI_Intercomm_create_from_groups(local_group, local_leader, remote_group, remote_leader, &
                                                stringtag, info, errhandler, newintercomm, ierr)
        IMPLICIT NONE
        INTEGER :: local_group
        INTEGER :: local_leader
        INTEGER :: remote_group
        INTEGER :: remote_leader
        CHARACTER*(*) :: stringtag
        INTEGER :: info
        INTEGER :: errhandler
        INTEGER :: newintercomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Intercomm_create_from_groups
    
    SUBROUTINE MPI_Intercomm_merge(intercomm, high, newintracomm, ierr)
        IMPLICIT NONE
        INTEGER :: intercomm
        LOGICAL :: high
        INTEGER :: newintracomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Intercomm_merge
    
    SUBROUTINE MPIX_Comm_revoke(comm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPIX_Comm_revoke
    
    SUBROUTINE MPIX_Comm_shrink(comm, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPIX_Comm_shrink
    
    SUBROUTINE MPIX_Comm_failure_ack(comm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPIX_Comm_failure_ack
    
    SUBROUTINE MPIX_Comm_failure_get_acked(comm, failedgrp, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: failedgrp
        INTEGER :: ierr
    END SUBROUTINE MPIX_Comm_failure_get_acked
    
    SUBROUTINE MPIX_Comm_agree(comm, flag, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPIX_Comm_agree
    
    SUBROUTINE MPIX_Comm_get_failed(comm, failedgrp, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: failedgrp
        INTEGER :: ierr
    END SUBROUTINE MPIX_Comm_get_failed
    
    SUBROUTINE MPI_Get_address(location, address, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: location
        REAL :: location
        INTEGER(KIND=MPI_ADDRESS_KIND) :: address
        INTEGER :: ierr
    END SUBROUTINE MPI_Get_address
    
    SUBROUTINE MPI_Get_count(status, datatype, count, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: datatype
        INTEGER :: count
        INTEGER :: ierr
    END SUBROUTINE MPI_Get_count
    
    SUBROUTINE MPI_Get_elements(status, datatype, count, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: datatype
        INTEGER :: count
        INTEGER :: ierr
    END SUBROUTINE MPI_Get_elements
    
    SUBROUTINE MPI_Get_elements_x(status, datatype, count, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE, MPI_COUNT_KIND
        IMPLICIT NONE
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: datatype
        INTEGER(KIND=MPI_COUNT_KIND) :: count
        INTEGER :: ierr
    END SUBROUTINE MPI_Get_elements_x
    
    SUBROUTINE MPI_Pack(inbuf, incount, datatype, outbuf, outsize, position, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: inbuf, outbuf
        REAL :: inbuf
        INTEGER :: incount
        INTEGER :: datatype
        REAL :: outbuf
        INTEGER :: outsize
        INTEGER :: position
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Pack
    
    SUBROUTINE MPI_Pack_external(datarep, inbuf, incount, datatype, outbuf, outsize, position, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: inbuf, outbuf
        CHARACTER*(*) :: datarep
        REAL :: inbuf
        INTEGER :: incount
        INTEGER :: datatype
        REAL :: outbuf
        INTEGER(KIND=MPI_ADDRESS_KIND) :: outsize
        INTEGER(KIND=MPI_ADDRESS_KIND) :: position
        INTEGER :: ierr
    END SUBROUTINE MPI_Pack_external
    
    SUBROUTINE MPI_Pack_external_size(datarep, incount, datatype, size, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        CHARACTER*(*) :: datarep
        INTEGER :: incount
        INTEGER :: datatype
        INTEGER(KIND=MPI_ADDRESS_KIND) :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Pack_external_size
    
    SUBROUTINE MPI_Pack_size(incount, datatype, comm, size, ierr)
        IMPLICIT NONE
        INTEGER :: incount
        INTEGER :: datatype
        INTEGER :: comm
        INTEGER :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Pack_size
    
    SUBROUTINE MPI_Status_set_elements(status, datatype, count, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: datatype
        INTEGER :: count
        INTEGER :: ierr
    END SUBROUTINE MPI_Status_set_elements
    
    SUBROUTINE MPI_Status_set_elements_x(status, datatype, count, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE, MPI_COUNT_KIND
        IMPLICIT NONE
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: datatype
        INTEGER(KIND=MPI_COUNT_KIND) :: count
        INTEGER :: ierr
    END SUBROUTINE MPI_Status_set_elements_x
    
    SUBROUTINE MPI_Type_commit(datatype, ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_commit
    
    SUBROUTINE MPI_Type_contiguous(count, oldtype, newtype, ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_contiguous
    
    SUBROUTINE MPI_Type_create_darray(size, rank, ndims, array_of_gsizes, array_of_distribs, &
                                      array_of_dargs, array_of_psizes, order, oldtype, newtype, ierr)
        IMPLICIT NONE
        INTEGER :: size
        INTEGER :: rank
        INTEGER :: ndims
        INTEGER :: array_of_gsizes(ndims)
        INTEGER :: array_of_distribs(ndims)
        INTEGER :: array_of_dargs(ndims)
        INTEGER :: array_of_psizes(ndims)
        INTEGER :: order
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_darray
    
    SUBROUTINE MPI_Type_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, &
                                        newtype, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_blocklengths(count)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: array_of_displacements(count)
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_hindexed
    
    SUBROUTINE MPI_Type_create_hindexed_block(count, blocklength, array_of_displacements, oldtype, &
                                              newtype, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: blocklength
        INTEGER(KIND=MPI_ADDRESS_KIND) :: array_of_displacements(count)
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_hindexed_block
    
    SUBROUTINE MPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: blocklength
        INTEGER(KIND=MPI_ADDRESS_KIND) :: stride
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_hvector
    
    SUBROUTINE MPI_Type_create_indexed_block(count, blocklength, array_of_displacements, oldtype, &
                                             newtype, ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: blocklength
        INTEGER :: array_of_displacements(count)
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_indexed_block
    
    SUBROUTINE MPI_Type_create_resized(oldtype, lb, extent, newtype, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: oldtype
        INTEGER(KIND=MPI_ADDRESS_KIND) :: lb
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extent
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_resized
    
    SUBROUTINE MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, &
                                      array_of_types, newtype, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_blocklengths(count)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: array_of_displacements(count)
        INTEGER :: array_of_types(count)
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_struct
    
    SUBROUTINE MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, &
                                        order, oldtype, newtype, ierr)
        IMPLICIT NONE
        INTEGER :: ndims
        INTEGER :: array_of_sizes(ndims)
        INTEGER :: array_of_subsizes(ndims)
        INTEGER :: array_of_starts(ndims)
        INTEGER :: order
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_subarray
    
    SUBROUTINE MPI_Type_dup(oldtype, newtype, ierr)
        IMPLICIT NONE
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_dup
    
    SUBROUTINE MPI_Type_free(datatype, ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_free
    
    SUBROUTINE MPI_Type_get_contents(datatype, max_integers, max_addresses, max_datatypes, &
                                     array_of_integers, array_of_addresses, array_of_datatypes, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: max_integers
        INTEGER :: max_addresses
        INTEGER :: max_datatypes
        INTEGER :: array_of_integers(max_integers)
        INTEGER(KIND=MPI_ADDRESS_KIND) :: array_of_addresses(max_addresses)
        INTEGER :: array_of_datatypes(max_datatypes)
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_get_contents
    
    SUBROUTINE MPI_Type_get_envelope(datatype, num_integers, num_addresses, num_datatypes, combiner, &
                                     ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: num_integers
        INTEGER :: num_addresses
        INTEGER :: num_datatypes
        INTEGER :: combiner
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_get_envelope
    
    SUBROUTINE MPI_Type_get_extent(datatype, lb, extent, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER(KIND=MPI_ADDRESS_KIND) :: lb
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extent
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_get_extent
    
    SUBROUTINE MPI_Type_get_extent_x(datatype, lb, extent, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_COUNT_KIND
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER(KIND=MPI_COUNT_KIND) :: lb
        INTEGER(KIND=MPI_COUNT_KIND) :: extent
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_get_extent_x
    
    SUBROUTINE MPI_Type_get_name(datatype, type_name, resultlen, ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        CHARACTER*(*) :: type_name
        INTEGER :: resultlen
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_get_name
    
    SUBROUTINE MPI_Type_get_true_extent(datatype, true_lb, true_extent, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER(KIND=MPI_ADDRESS_KIND) :: true_lb
        INTEGER(KIND=MPI_ADDRESS_KIND) :: true_extent
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_get_true_extent
    
    SUBROUTINE MPI_Type_get_true_extent_x(datatype, true_lb, true_extent, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_COUNT_KIND
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER(KIND=MPI_COUNT_KIND) :: true_lb
        INTEGER(KIND=MPI_COUNT_KIND) :: true_extent
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_get_true_extent_x
    
    SUBROUTINE MPI_Type_indexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype, &
                                ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_blocklengths(count)
        INTEGER :: array_of_displacements(count)
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_indexed
    
    SUBROUTINE MPI_Type_match_size(typeclass, size, datatype, ierr)
        IMPLICIT NONE
        INTEGER :: typeclass
        INTEGER :: size
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_match_size
    
    SUBROUTINE MPI_Type_set_name(datatype, type_name, ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        CHARACTER*(*) :: type_name
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_set_name
    
    SUBROUTINE MPI_Type_size(datatype, size, ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_size
    
    SUBROUTINE MPI_Type_size_x(datatype, size, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_COUNT_KIND
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER(KIND=MPI_COUNT_KIND) :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_size_x
    
    SUBROUTINE MPI_Type_vector(count, blocklength, stride, oldtype, newtype, ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: blocklength
        INTEGER :: stride
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_vector
    
    SUBROUTINE MPI_Unpack(inbuf, insize, position, outbuf, outcount, datatype, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: inbuf, outbuf
        REAL :: inbuf
        INTEGER :: insize
        INTEGER :: position
        REAL :: outbuf
        INTEGER :: outcount
        INTEGER :: datatype
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Unpack
    
    SUBROUTINE MPI_Unpack_external(datarep, inbuf, insize, position, outbuf, outcount, datatype, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: inbuf, outbuf
        CHARACTER*(*) :: datarep
        REAL :: inbuf
        INTEGER(KIND=MPI_ADDRESS_KIND) :: insize
        INTEGER(KIND=MPI_ADDRESS_KIND) :: position
        REAL :: outbuf
        INTEGER :: outcount
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_Unpack_external
    
    SUBROUTINE MPI_Address(location, address, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: location
        REAL :: location
        INTEGER(KIND=MPI_ADDRESS_KIND) :: address
        INTEGER :: ierr
    END SUBROUTINE MPI_Address
    
    SUBROUTINE MPI_Type_extent(datatype, extent, ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: extent
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_extent
    
    SUBROUTINE MPI_Type_lb(datatype, displacement, ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: displacement
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_lb
    
    SUBROUTINE MPI_Type_ub(datatype, displacement, ierr)
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: displacement
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_ub
    
    SUBROUTINE MPI_Type_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype, &
                                 ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_blocklengths(count)
        INTEGER :: array_of_displacements(count)
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_hindexed
    
    SUBROUTINE MPI_Type_hvector(count, blocklength, stride, oldtype, newtype, ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: blocklength
        INTEGER :: stride
        INTEGER :: oldtype
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_hvector
    
    SUBROUTINE MPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, &
                               newtype, ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_blocklengths(count)
        INTEGER :: array_of_displacements(count)
        INTEGER :: array_of_types(count)
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_struct
    
    SUBROUTINE MPI_Add_error_class(errorclass, ierr)
        IMPLICIT NONE
        INTEGER :: errorclass
        INTEGER :: ierr
    END SUBROUTINE MPI_Add_error_class
    
    SUBROUTINE MPI_Add_error_code(errorclass, errorcode, ierr)
        IMPLICIT NONE
        INTEGER :: errorclass
        INTEGER :: errorcode
        INTEGER :: ierr
    END SUBROUTINE MPI_Add_error_code
    
    SUBROUTINE MPI_Add_error_string(errorcode, string, ierr)
        IMPLICIT NONE
        INTEGER :: errorcode
        CHARACTER*(*) :: string
        INTEGER :: ierr
    END SUBROUTINE MPI_Add_error_string
    
    SUBROUTINE MPI_Comm_call_errhandler(comm, errorcode, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: errorcode
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_call_errhandler
    
    SUBROUTINE MPI_Comm_create_errhandler(comm_errhandler_fn, errhandler, ierr)
        IMPLICIT NONE
        EXTERNAL :: comm_errhandler_fn
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_create_errhandler
    
    SUBROUTINE MPI_Comm_get_errhandler(comm, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_get_errhandler
    
    SUBROUTINE MPI_Comm_set_errhandler(comm, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_set_errhandler
    
    SUBROUTINE MPI_Errhandler_free(errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Errhandler_free
    
    SUBROUTINE MPI_Error_class(errorcode, errorclass, ierr)
        IMPLICIT NONE
        INTEGER :: errorcode
        INTEGER :: errorclass
        INTEGER :: ierr
    END SUBROUTINE MPI_Error_class
    
    SUBROUTINE MPI_Error_string(errorcode, string, resultlen, ierr)
        IMPLICIT NONE
        INTEGER :: errorcode
        CHARACTER*(*) :: string
        INTEGER :: resultlen
        INTEGER :: ierr
    END SUBROUTINE MPI_Error_string
    
    SUBROUTINE MPI_File_call_errhandler(fh, errorcode, ierr)
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER :: errorcode
        INTEGER :: ierr
    END SUBROUTINE MPI_File_call_errhandler
    
    SUBROUTINE MPI_File_create_errhandler(file_errhandler_fn, errhandler, ierr)
        IMPLICIT NONE
        EXTERNAL :: file_errhandler_fn
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_File_create_errhandler
    
    SUBROUTINE MPI_File_get_errhandler(file, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: file
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_errhandler
    
    SUBROUTINE MPI_File_set_errhandler(file, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: file
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_File_set_errhandler
    
    SUBROUTINE MPI_Session_call_errhandler(session, errorcode, ierr)
        IMPLICIT NONE
        INTEGER :: session
        INTEGER :: errorcode
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_call_errhandler
    
    SUBROUTINE MPI_Session_create_errhandler(session_errhandler_fn, errhandler, ierr)
        IMPLICIT NONE
        EXTERNAL :: session_errhandler_fn
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_create_errhandler
    
    SUBROUTINE MPI_Session_get_errhandler(session, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: session
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_get_errhandler
    
    SUBROUTINE MPI_Session_set_errhandler(session, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: session
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_set_errhandler
    
    SUBROUTINE MPI_Win_call_errhandler(win, errorcode, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: errorcode
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_call_errhandler
    
    SUBROUTINE MPI_Win_create_errhandler(win_errhandler_fn, errhandler, ierr)
        IMPLICIT NONE
        EXTERNAL :: win_errhandler_fn
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_create_errhandler
    
    SUBROUTINE MPI_Win_get_errhandler(win, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_get_errhandler
    
    SUBROUTINE MPI_Win_set_errhandler(win, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_set_errhandler
    
    SUBROUTINE MPIX_Delete_error_class(errorclass, ierr)
        IMPLICIT NONE
        INTEGER :: errorclass
        INTEGER :: ierr
    END SUBROUTINE MPIX_Delete_error_class
    
    SUBROUTINE MPIX_Delete_error_code(errorcode, ierr)
        IMPLICIT NONE
        INTEGER :: errorcode
        INTEGER :: ierr
    END SUBROUTINE MPIX_Delete_error_code
    
    SUBROUTINE MPIX_Delete_error_string(errorcode, ierr)
        IMPLICIT NONE
        INTEGER :: errorcode
        INTEGER :: ierr
    END SUBROUTINE MPIX_Delete_error_string
    
    SUBROUTINE MPI_Errhandler_create(comm_errhandler_fn, errhandler, ierr)
        IMPLICIT NONE
        EXTERNAL :: comm_errhandler_fn
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Errhandler_create
    
    SUBROUTINE MPI_Errhandler_get(comm, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Errhandler_get
    
    SUBROUTINE MPI_Errhandler_set(comm, errhandler, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: errhandler
        INTEGER :: ierr
    END SUBROUTINE MPI_Errhandler_set
    
    SUBROUTINE MPI_Group_compare(group1, group2, result, ierr)
        IMPLICIT NONE
        INTEGER :: group1
        INTEGER :: group2
        INTEGER :: result
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_compare
    
    SUBROUTINE MPI_Group_difference(group1, group2, newgroup, ierr)
        IMPLICIT NONE
        INTEGER :: group1
        INTEGER :: group2
        INTEGER :: newgroup
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_difference
    
    SUBROUTINE MPI_Group_excl(group, n, ranks, newgroup, ierr)
        IMPLICIT NONE
        INTEGER :: group
        INTEGER :: n
        INTEGER :: ranks(n)
        INTEGER :: newgroup
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_excl
    
    SUBROUTINE MPI_Group_free(group, ierr)
        IMPLICIT NONE
        INTEGER :: group
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_free
    
    SUBROUTINE MPI_Group_incl(group, n, ranks, newgroup, ierr)
        IMPLICIT NONE
        INTEGER :: group
        INTEGER :: n
        INTEGER :: ranks(n)
        INTEGER :: newgroup
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_incl
    
    SUBROUTINE MPI_Group_intersection(group1, group2, newgroup, ierr)
        IMPLICIT NONE
        INTEGER :: group1
        INTEGER :: group2
        INTEGER :: newgroup
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_intersection
    
    SUBROUTINE MPI_Group_range_excl(group, n, ranges, newgroup, ierr)
        IMPLICIT NONE
        INTEGER :: group
        INTEGER :: n
        INTEGER :: ranges(3, *)
        INTEGER :: newgroup
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_range_excl
    
    SUBROUTINE MPI_Group_range_incl(group, n, ranges, newgroup, ierr)
        IMPLICIT NONE
        INTEGER :: group
        INTEGER :: n
        INTEGER :: ranges(3, *)
        INTEGER :: newgroup
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_range_incl
    
    SUBROUTINE MPI_Group_rank(group, rank, ierr)
        IMPLICIT NONE
        INTEGER :: group
        INTEGER :: rank
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_rank
    
    SUBROUTINE MPI_Group_size(group, size, ierr)
        IMPLICIT NONE
        INTEGER :: group
        INTEGER :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_size
    
    SUBROUTINE MPI_Group_translate_ranks(group1, n, ranks1, group2, ranks2, ierr)
        IMPLICIT NONE
        INTEGER :: group1
        INTEGER :: n
        INTEGER :: ranks1(n)
        INTEGER :: group2
        INTEGER :: ranks2(n)
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_translate_ranks
    
    SUBROUTINE MPI_Group_union(group1, group2, newgroup, ierr)
        IMPLICIT NONE
        INTEGER :: group1
        INTEGER :: group2
        INTEGER :: newgroup
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_union
    
    SUBROUTINE MPI_Info_create(info, ierr)
        IMPLICIT NONE
        INTEGER :: info
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_create
    
    SUBROUTINE MPI_Info_create_env(info, ierr)
        IMPLICIT NONE
        INTEGER :: info
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_create_env
    
    SUBROUTINE MPI_Info_delete(info, key, ierr)
        IMPLICIT NONE
        INTEGER :: info
        CHARACTER*(*) :: key
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_delete
    
    SUBROUTINE MPI_Info_dup(info, newinfo, ierr)
        IMPLICIT NONE
        INTEGER :: info
        INTEGER :: newinfo
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_dup
    
    SUBROUTINE MPI_Info_free(info, ierr)
        IMPLICIT NONE
        INTEGER :: info
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_free
    
    SUBROUTINE MPI_Info_get(info, key, valuelen, value, flag, ierr)
        IMPLICIT NONE
        INTEGER :: info
        CHARACTER*(*) :: key
        INTEGER :: valuelen
        CHARACTER*(*) :: value
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_get
    
    SUBROUTINE MPI_Info_get_nkeys(info, nkeys, ierr)
        IMPLICIT NONE
        INTEGER :: info
        INTEGER :: nkeys
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_get_nkeys
    
    SUBROUTINE MPI_Info_get_nthkey(info, n, key, ierr)
        IMPLICIT NONE
        INTEGER :: info
        INTEGER :: n
        CHARACTER*(*) :: key
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_get_nthkey
    
    SUBROUTINE MPI_Info_get_string(info, key, buflen, value, flag, ierr)
        IMPLICIT NONE
        INTEGER :: info
        CHARACTER*(*) :: key
        INTEGER :: buflen
        CHARACTER*(*) :: value
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_get_string
    
    SUBROUTINE MPI_Info_get_valuelen(info, key, valuelen, flag, ierr)
        IMPLICIT NONE
        INTEGER :: info
        CHARACTER*(*) :: key
        INTEGER :: valuelen
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_get_valuelen
    
    SUBROUTINE MPI_Info_set(info, key, value, ierr)
        IMPLICIT NONE
        INTEGER :: info
        CHARACTER*(*) :: key
        CHARACTER*(*) :: value
        INTEGER :: ierr
    END SUBROUTINE MPI_Info_set
    
    SUBROUTINE MPIX_Info_set_hex(info, key, value, value_size, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: value
        INTEGER :: info
        CHARACTER*(*) :: key
        REAL :: value
        INTEGER :: value_size
        INTEGER :: ierr
    END SUBROUTINE MPIX_Info_set_hex
    
    SUBROUTINE MPI_Abort(comm, errorcode, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: errorcode
        INTEGER :: ierr
    END SUBROUTINE MPI_Abort
    
    SUBROUTINE MPI_Comm_create_from_group(group, stringtag, info, errhandler, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: group
        CHARACTER*(*) :: stringtag
        INTEGER :: info
        INTEGER :: errhandler
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_create_from_group
    
    SUBROUTINE MPI_Finalize(ierr)
        IMPLICIT NONE
        INTEGER :: ierr
    END SUBROUTINE MPI_Finalize
    
    SUBROUTINE MPI_Finalized(flag, ierr)
        IMPLICIT NONE
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Finalized
    
    SUBROUTINE MPI_Group_from_session_pset(session, pset_name, newgroup, ierr)
        IMPLICIT NONE
        INTEGER :: session
        CHARACTER*(*) :: pset_name
        INTEGER :: newgroup
        INTEGER :: ierr
    END SUBROUTINE MPI_Group_from_session_pset
    
    SUBROUTINE MPI_Init(ierr)
        IMPLICIT NONE
        INTEGER :: ierr
    END SUBROUTINE MPI_Init
    
    SUBROUTINE MPI_Init_thread(required, provided, ierr)
        IMPLICIT NONE
        INTEGER :: required
        INTEGER :: provided
        INTEGER :: ierr
    END SUBROUTINE MPI_Init_thread
    
    SUBROUTINE MPI_Initialized(flag, ierr)
        IMPLICIT NONE
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Initialized
    
    SUBROUTINE MPI_Is_thread_main(flag, ierr)
        IMPLICIT NONE
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Is_thread_main
    
    SUBROUTINE MPI_Query_thread(provided, ierr)
        IMPLICIT NONE
        INTEGER :: provided
        INTEGER :: ierr
    END SUBROUTINE MPI_Query_thread
    
    SUBROUTINE MPI_Session_finalize(session, ierr)
        IMPLICIT NONE
        INTEGER :: session
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_finalize
    
    SUBROUTINE MPI_Session_get_info(session, info_used, ierr)
        IMPLICIT NONE
        INTEGER :: session
        INTEGER :: info_used
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_get_info
    
    SUBROUTINE MPI_Session_get_nth_pset(session, info, n, pset_len, pset_name, ierr)
        IMPLICIT NONE
        INTEGER :: session
        INTEGER :: info
        INTEGER :: n
        INTEGER :: pset_len
        CHARACTER*(*) :: pset_name
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_get_nth_pset
    
    SUBROUTINE MPI_Session_get_num_psets(session, info, npset_names, ierr)
        IMPLICIT NONE
        INTEGER :: session
        INTEGER :: info
        INTEGER :: npset_names
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_get_num_psets
    
    SUBROUTINE MPI_Session_get_pset_info(session, pset_name, info, ierr)
        IMPLICIT NONE
        INTEGER :: session
        CHARACTER*(*) :: pset_name
        INTEGER :: info
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_get_pset_info
    
    SUBROUTINE MPI_Session_init(info, errhandler, session, ierr)
        IMPLICIT NONE
        INTEGER :: info
        INTEGER :: errhandler
        INTEGER :: session
        INTEGER :: ierr
    END SUBROUTINE MPI_Session_init
    
    FUNCTION MPI_Aint_add(base, disp) result(res)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER(KIND=MPI_ADDRESS_KIND) :: base
        INTEGER(KIND=MPI_ADDRESS_KIND) :: disp
        INTEGER(KIND=MPI_ADDRESS_KIND) :: res
    END FUNCTION MPI_Aint_add
    
    FUNCTION MPI_Aint_diff(addr1, addr2) result(res)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER(KIND=MPI_ADDRESS_KIND) :: addr1
        INTEGER(KIND=MPI_ADDRESS_KIND) :: addr2
        INTEGER(KIND=MPI_ADDRESS_KIND) :: res
    END FUNCTION MPI_Aint_diff
    
    SUBROUTINE MPI_Get_library_version(version, resultlen, ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: version
        INTEGER :: resultlen
        INTEGER :: ierr
    END SUBROUTINE MPI_Get_library_version
    
    SUBROUTINE MPI_Get_processor_name(name, resultlen, ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: name
        INTEGER :: resultlen
        INTEGER :: ierr
    END SUBROUTINE MPI_Get_processor_name
    
    SUBROUTINE MPI_Get_version(version, subversion, ierr)
        IMPLICIT NONE
        INTEGER :: version
        INTEGER :: subversion
        INTEGER :: ierr
    END SUBROUTINE MPI_Get_version
    
    SUBROUTINE MPIX_GPU_query_support(gpu_type, is_supported, ierr)
        IMPLICIT NONE
        INTEGER :: gpu_type
        LOGICAL :: is_supported
        INTEGER :: ierr
    END SUBROUTINE MPIX_GPU_query_support
    
    SUBROUTINE MPIX_Query_cuda_support(ierr)
        IMPLICIT NONE
        INTEGER :: ierr
    END SUBROUTINE MPIX_Query_cuda_support
    
    SUBROUTINE MPIX_Query_ze_support(ierr)
        IMPLICIT NONE
        INTEGER :: ierr
    END SUBROUTINE MPIX_Query_ze_support
    
    SUBROUTINE MPIX_Query_hip_support(ierr)
        IMPLICIT NONE
        INTEGER :: ierr
    END SUBROUTINE MPIX_Query_hip_support
    
    SUBROUTINE MPI_Op_commutative(op, commute, ierr)
        IMPLICIT NONE
        INTEGER :: op
        LOGICAL :: commute
        INTEGER :: ierr
    END SUBROUTINE MPI_Op_commutative
    
    SUBROUTINE MPI_Op_create(user_fn, commute, op, ierr)
        IMPLICIT NONE
        EXTERNAL :: user_fn
        LOGICAL :: commute
        INTEGER :: op
        INTEGER :: ierr
    END SUBROUTINE MPI_Op_create
    
    SUBROUTINE MPI_Op_free(op, ierr)
        IMPLICIT NONE
        INTEGER :: op
        INTEGER :: ierr
    END SUBROUTINE MPI_Op_free
    
    SUBROUTINE MPI_Parrived(request, partition, flag, ierr)
        IMPLICIT NONE
        INTEGER :: request
        INTEGER :: partition
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Parrived
    
    SUBROUTINE MPI_Pready(partition, request, ierr)
        IMPLICIT NONE
        INTEGER :: partition
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Pready
    
    SUBROUTINE MPI_Pready_list(length, array_of_partitions, request, ierr)
        IMPLICIT NONE
        INTEGER :: length
        INTEGER :: array_of_partitions(length)
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Pready_list
    
    SUBROUTINE MPI_Pready_range(partition_low, partition_high, request, ierr)
        IMPLICIT NONE
        INTEGER :: partition_low
        INTEGER :: partition_high
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Pready_range
    
    SUBROUTINE MPI_Precv_init(buf, partitions, count, datatype, dest, tag, comm, info, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_COUNT_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: partitions
        INTEGER(KIND=MPI_COUNT_KIND) :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Precv_init
    
    SUBROUTINE MPI_Psend_init(buf, partitions, count, datatype, dest, tag, comm, info, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_COUNT_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: partitions
        INTEGER(KIND=MPI_COUNT_KIND) :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: info
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Psend_init
    
    SUBROUTINE MPI_Bsend(buf, count, datatype, dest, tag, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Bsend
    
    SUBROUTINE MPI_Bsend_init(buf, count, datatype, dest, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Bsend_init
    
    SUBROUTINE MPI_Buffer_attach(buffer, size, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buffer
        REAL :: buffer
        INTEGER :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Buffer_attach
    
    SUBROUTINE MPI_Buffer_detach(buffer_addr, size, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buffer_addr
        REAL :: buffer_addr
        INTEGER :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Buffer_detach
    
    SUBROUTINE MPI_Ibsend(buf, count, datatype, dest, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ibsend
    
    SUBROUTINE MPI_Improbe(source, tag, comm, flag, message, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        LOGICAL :: flag
        INTEGER :: message
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Improbe
    
    SUBROUTINE MPI_Imrecv(buf, count, datatype, message, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: message
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Imrecv
    
    SUBROUTINE MPI_Iprobe(source, tag, comm, flag, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        LOGICAL :: flag
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Iprobe
    
    SUBROUTINE MPI_Irecv(buf, count, datatype, source, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Irecv
    
    SUBROUTINE MPI_Irsend(buf, count, datatype, dest, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Irsend
    
    SUBROUTINE MPI_Isend(buf, count, datatype, dest, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Isend
    
    SUBROUTINE MPI_Isendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, &
                             source, recvtag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        INTEGER :: dest
        INTEGER :: sendtag
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: source
        INTEGER :: recvtag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Isendrecv
    
    SUBROUTINE MPI_Isendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, &
                                     request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: sendtag
        INTEGER :: source
        INTEGER :: recvtag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Isendrecv_replace
    
    SUBROUTINE MPI_Issend(buf, count, datatype, dest, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Issend
    
    SUBROUTINE MPI_Mprobe(source, tag, comm, message, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: message
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Mprobe
    
    SUBROUTINE MPI_Mrecv(buf, count, datatype, message, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: message
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Mrecv
    
    SUBROUTINE MPI_Probe(source, tag, comm, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Probe
    
    SUBROUTINE MPI_Recv(buf, count, datatype, source, tag, comm, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Recv
    
    SUBROUTINE MPI_Recv_init(buf, count, datatype, source, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Recv_init
    
    SUBROUTINE MPI_Rsend(buf, count, datatype, dest, tag, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Rsend
    
    SUBROUTINE MPI_Rsend_init(buf, count, datatype, dest, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Rsend_init
    
    SUBROUTINE MPI_Send(buf, count, datatype, dest, tag, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Send
    
    SUBROUTINE MPI_Send_init(buf, count, datatype, dest, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Send_init
    
    SUBROUTINE MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, &
                            source, recvtag, comm, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        INTEGER :: sendcount
        INTEGER :: sendtype
        INTEGER :: dest
        INTEGER :: sendtag
        REAL :: recvbuf
        INTEGER :: recvcount
        INTEGER :: recvtype
        INTEGER :: source
        INTEGER :: recvtag
        INTEGER :: comm
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Sendrecv
    
    SUBROUTINE MPI_Sendrecv_replace(buf, count, datatype, dest, sendtag, source, recvtag, comm, status, &
                                    ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: sendtag
        INTEGER :: source
        INTEGER :: recvtag
        INTEGER :: comm
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Sendrecv_replace
    
    SUBROUTINE MPI_Ssend(buf, count, datatype, dest, tag, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Ssend
    
    SUBROUTINE MPI_Ssend_init(buf, count, datatype, dest, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Ssend_init
    
    SUBROUTINE MPI_Cancel(request, ierr)
        IMPLICIT NONE
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Cancel
    
    SUBROUTINE MPI_Grequest_complete(request, ierr)
        IMPLICIT NONE
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Grequest_complete
    
    SUBROUTINE MPI_Grequest_start(query_fn, free_fn, cancel_fn, extra_state, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        EXTERNAL :: query_fn
        EXTERNAL :: free_fn
        EXTERNAL :: cancel_fn
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Grequest_start
    
    SUBROUTINE MPI_Request_free(request, ierr)
        IMPLICIT NONE
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Request_free
    
    SUBROUTINE MPI_Request_get_status(request, flag, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: request
        LOGICAL :: flag
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Request_get_status
    
    SUBROUTINE MPI_Start(request, ierr)
        IMPLICIT NONE
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Start
    
    SUBROUTINE MPI_Startall(count, array_of_requests, ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_requests(count)
        INTEGER :: ierr
    END SUBROUTINE MPI_Startall
    
    SUBROUTINE MPI_Status_set_cancelled(status, flag, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: status(MPI_STATUS_SIZE)
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Status_set_cancelled
    
    SUBROUTINE MPI_Test(request, flag, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: request
        LOGICAL :: flag
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Test
    
    SUBROUTINE MPI_Test_cancelled(status, flag, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: status(MPI_STATUS_SIZE)
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Test_cancelled
    
    SUBROUTINE MPI_Testall(count, array_of_requests, flag, array_of_statuses, ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_requests(count)
        LOGICAL :: flag
        INTEGER :: array_of_statuses(*)
        INTEGER :: ierr
    END SUBROUTINE MPI_Testall
    
    SUBROUTINE MPI_Testany(count, array_of_requests, indx, flag, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_requests(count)
        INTEGER :: indx
        LOGICAL :: flag
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Testany
    
    SUBROUTINE MPI_Testsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses, &
                            ierr)
        IMPLICIT NONE
        INTEGER :: incount
        INTEGER :: array_of_requests(incount)
        INTEGER :: outcount
        INTEGER :: array_of_indices(*)
        INTEGER :: array_of_statuses(*)
        INTEGER :: ierr
    END SUBROUTINE MPI_Testsome
    
    SUBROUTINE MPI_Wait(request, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: request
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Wait
    
    SUBROUTINE MPI_Waitall(count, array_of_requests, array_of_statuses, ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_requests(count)
        INTEGER :: array_of_statuses(*)
        INTEGER :: ierr
    END SUBROUTINE MPI_Waitall
    
    SUBROUTINE MPI_Waitany(count, array_of_requests, indx, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_requests(count)
        INTEGER :: indx
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_Waitany
    
    SUBROUTINE MPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses, &
                            ierr)
        IMPLICIT NONE
        INTEGER :: incount
        INTEGER :: array_of_requests(incount)
        INTEGER :: outcount
        INTEGER :: array_of_indices(*)
        INTEGER :: array_of_statuses(*)
        INTEGER :: ierr
    END SUBROUTINE MPI_Waitsome
    
    SUBROUTINE MPI_Accumulate(origin_addr, origin_count, origin_datatype, target_rank, target_disp, &
                              target_count, target_datatype, op, win, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr
        REAL :: origin_addr
        INTEGER :: origin_count
        INTEGER :: origin_datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: target_count
        INTEGER :: target_datatype
        INTEGER :: op
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Accumulate
    
    SUBROUTINE MPI_Alloc_mem(size, info, baseptr, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER(KIND=MPI_ADDRESS_KIND) :: size
        INTEGER :: info
        INTEGER(KIND=MPI_ADDRESS_KIND) :: baseptr
        INTEGER :: ierr
    END SUBROUTINE MPI_Alloc_mem
    
    SUBROUTINE MPI_Compare_and_swap(origin_addr, compare_addr, result_addr, datatype, target_rank, &
                                    target_disp, win, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr, compare_addr, result_addr
        REAL :: origin_addr
        REAL :: compare_addr
        REAL :: result_addr
        INTEGER :: datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Compare_and_swap
    
    SUBROUTINE MPI_Fetch_and_op(origin_addr, result_addr, datatype, target_rank, target_disp, op, win, &
                                ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr, result_addr
        REAL :: origin_addr
        REAL :: result_addr
        INTEGER :: datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: op
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Fetch_and_op
    
    SUBROUTINE MPI_Free_mem(base, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: base
        REAL :: base
        INTEGER :: ierr
    END SUBROUTINE MPI_Free_mem
    
    SUBROUTINE MPI_Get(origin_addr, origin_count, origin_datatype, target_rank, target_disp, &
                       target_count, target_datatype, win, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr
        REAL :: origin_addr
        INTEGER :: origin_count
        INTEGER :: origin_datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: target_count
        INTEGER :: target_datatype
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Get
    
    SUBROUTINE MPI_Get_accumulate(origin_addr, origin_count, origin_datatype, result_addr, result_count, &
                                  result_datatype, target_rank, target_disp, target_count, &
                                  target_datatype, op, win, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr, result_addr
        REAL :: origin_addr
        INTEGER :: origin_count
        INTEGER :: origin_datatype
        REAL :: result_addr
        INTEGER :: result_count
        INTEGER :: result_datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: target_count
        INTEGER :: target_datatype
        INTEGER :: op
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Get_accumulate
    
    SUBROUTINE MPI_Put(origin_addr, origin_count, origin_datatype, target_rank, target_disp, &
                       target_count, target_datatype, win, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr
        REAL :: origin_addr
        INTEGER :: origin_count
        INTEGER :: origin_datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: target_count
        INTEGER :: target_datatype
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Put
    
    SUBROUTINE MPI_Raccumulate(origin_addr, origin_count, origin_datatype, target_rank, target_disp, &
                               target_count, target_datatype, op, win, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr
        REAL :: origin_addr
        INTEGER :: origin_count
        INTEGER :: origin_datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: target_count
        INTEGER :: target_datatype
        INTEGER :: op
        INTEGER :: win
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Raccumulate
    
    SUBROUTINE MPI_Rget(origin_addr, origin_count, origin_datatype, target_rank, target_disp, &
                        target_count, target_datatype, win, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr
        REAL :: origin_addr
        INTEGER :: origin_count
        INTEGER :: origin_datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: target_count
        INTEGER :: target_datatype
        INTEGER :: win
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Rget
    
    SUBROUTINE MPI_Rget_accumulate(origin_addr, origin_count, origin_datatype, result_addr, &
                                   result_count, result_datatype, target_rank, target_disp, &
                                   target_count, target_datatype, op, win, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr, result_addr
        REAL :: origin_addr
        INTEGER :: origin_count
        INTEGER :: origin_datatype
        REAL :: result_addr
        INTEGER :: result_count
        INTEGER :: result_datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: target_count
        INTEGER :: target_datatype
        INTEGER :: op
        INTEGER :: win
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Rget_accumulate
    
    SUBROUTINE MPI_Rput(origin_addr, origin_count, origin_datatype, target_rank, target_disp, &
                        target_count, target_datatype, win, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: origin_addr
        REAL :: origin_addr
        INTEGER :: origin_count
        INTEGER :: origin_datatype
        INTEGER :: target_rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: target_disp
        INTEGER :: target_count
        INTEGER :: target_datatype
        INTEGER :: win
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_Rput
    
    SUBROUTINE MPI_Win_allocate(size, disp_unit, info, comm, baseptr, win, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER(KIND=MPI_ADDRESS_KIND) :: size
        INTEGER :: disp_unit
        INTEGER :: info
        INTEGER :: comm
        INTEGER(KIND=MPI_ADDRESS_KIND) :: baseptr
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_allocate
    
    SUBROUTINE MPI_Win_allocate_shared(size, disp_unit, info, comm, baseptr, win, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER(KIND=MPI_ADDRESS_KIND) :: size
        INTEGER :: disp_unit
        INTEGER :: info
        INTEGER :: comm
        INTEGER(KIND=MPI_ADDRESS_KIND) :: baseptr
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_allocate_shared
    
    SUBROUTINE MPI_Win_attach(win, base, size, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: base
        INTEGER :: win
        REAL :: base
        INTEGER(KIND=MPI_ADDRESS_KIND) :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_attach
    
    SUBROUTINE MPI_Win_complete(win, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_complete
    
    SUBROUTINE MPI_Win_create(base, size, disp_unit, info, comm, win, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: base
        REAL :: base
        INTEGER(KIND=MPI_ADDRESS_KIND) :: size
        INTEGER :: disp_unit
        INTEGER :: info
        INTEGER :: comm
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_create
    
    SUBROUTINE MPI_Win_create_dynamic(info, comm, win, ierr)
        IMPLICIT NONE
        INTEGER :: info
        INTEGER :: comm
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_create_dynamic
    
    SUBROUTINE MPI_Win_detach(win, base, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: base
        INTEGER :: win
        REAL :: base
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_detach
    
    SUBROUTINE MPI_Win_fence(assert, win, ierr)
        IMPLICIT NONE
        INTEGER :: assert
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_fence
    
    SUBROUTINE MPI_Win_flush(rank, win, ierr)
        IMPLICIT NONE
        INTEGER :: rank
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_flush
    
    SUBROUTINE MPI_Win_flush_all(win, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_flush_all
    
    SUBROUTINE MPI_Win_flush_local(rank, win, ierr)
        IMPLICIT NONE
        INTEGER :: rank
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_flush_local
    
    SUBROUTINE MPI_Win_flush_local_all(win, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_flush_local_all
    
    SUBROUTINE MPI_Win_free(win, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_free
    
    SUBROUTINE MPI_Win_get_group(win, group, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: group
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_get_group
    
    SUBROUTINE MPI_Win_get_info(win, info_used, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: info_used
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_get_info
    
    SUBROUTINE MPI_Win_get_name(win, win_name, resultlen, ierr)
        IMPLICIT NONE
        INTEGER :: win
        CHARACTER*(*) :: win_name
        INTEGER :: resultlen
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_get_name
    
    SUBROUTINE MPI_Win_lock(lock_type, rank, assert, win, ierr)
        IMPLICIT NONE
        INTEGER :: lock_type
        INTEGER :: rank
        INTEGER :: assert
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_lock
    
    SUBROUTINE MPI_Win_lock_all(assert, win, ierr)
        IMPLICIT NONE
        INTEGER :: assert
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_lock_all
    
    SUBROUTINE MPI_Win_post(group, assert, win, ierr)
        IMPLICIT NONE
        INTEGER :: group
        INTEGER :: assert
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_post
    
    SUBROUTINE MPI_Win_set_info(win, info, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: info
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_set_info
    
    SUBROUTINE MPI_Win_set_name(win, win_name, ierr)
        IMPLICIT NONE
        INTEGER :: win
        CHARACTER*(*) :: win_name
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_set_name
    
    SUBROUTINE MPI_Win_shared_query(win, rank, size, disp_unit, baseptr, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: rank
        INTEGER(KIND=MPI_ADDRESS_KIND) :: size
        INTEGER :: disp_unit
        INTEGER(KIND=MPI_ADDRESS_KIND) :: baseptr
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_shared_query
    
    SUBROUTINE MPI_Win_start(group, assert, win, ierr)
        IMPLICIT NONE
        INTEGER :: group
        INTEGER :: assert
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_start
    
    SUBROUTINE MPI_Win_sync(win, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_sync
    
    SUBROUTINE MPI_Win_test(win, flag, ierr)
        IMPLICIT NONE
        INTEGER :: win
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_test
    
    SUBROUTINE MPI_Win_unlock(rank, win, ierr)
        IMPLICIT NONE
        INTEGER :: rank
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_unlock
    
    SUBROUTINE MPI_Win_unlock_all(win, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_unlock_all
    
    SUBROUTINE MPI_Win_wait(win, ierr)
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: ierr
    END SUBROUTINE MPI_Win_wait
    
    SUBROUTINE MPI_Close_port(port_name, ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: port_name
        INTEGER :: ierr
    END SUBROUTINE MPI_Close_port
    
    SUBROUTINE MPI_Comm_accept(port_name, info, root, comm, newcomm, ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: port_name
        INTEGER :: info
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_accept
    
    SUBROUTINE MPI_Comm_connect(port_name, info, root, comm, newcomm, ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: port_name
        INTEGER :: info
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_connect
    
    SUBROUTINE MPI_Comm_disconnect(comm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_disconnect
    
    SUBROUTINE MPI_Comm_get_parent(parent, ierr)
        IMPLICIT NONE
        INTEGER :: parent
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_get_parent
    
    SUBROUTINE MPI_Comm_join(fd, intercomm, ierr)
        IMPLICIT NONE
        INTEGER :: fd
        INTEGER :: intercomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_join
    
    SUBROUTINE MPI_Comm_spawn(command, argv, maxprocs, info, root, comm, intercomm, array_of_errcodes, &
                              ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: command
        CHARACTER*(*) :: argv(*)
        INTEGER :: maxprocs
        INTEGER :: info
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: intercomm
        INTEGER :: array_of_errcodes(*)
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_spawn
    
    SUBROUTINE MPI_Comm_spawn_multiple(count, array_of_commands, array_of_argv, array_of_maxprocs, &
                                       array_of_info, root, comm, intercomm, array_of_errcodes, ierr)
        IMPLICIT NONE
        INTEGER :: count
        CHARACTER*(*) :: array_of_commands(*)
        CHARACTER*(*) :: array_of_argv(count, *)
        INTEGER :: array_of_maxprocs(count)
        INTEGER :: array_of_info(count)
        INTEGER :: root
        INTEGER :: comm
        INTEGER :: intercomm
        INTEGER :: array_of_errcodes(*)
        INTEGER :: ierr
    END SUBROUTINE MPI_Comm_spawn_multiple
    
    SUBROUTINE MPI_Lookup_name(service_name, info, port_name, ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: service_name
        INTEGER :: info
        CHARACTER*(*) :: port_name
        INTEGER :: ierr
    END SUBROUTINE MPI_Lookup_name
    
    SUBROUTINE MPI_Open_port(info, port_name, ierr)
        IMPLICIT NONE
        INTEGER :: info
        CHARACTER*(*) :: port_name
        INTEGER :: ierr
    END SUBROUTINE MPI_Open_port
    
    SUBROUTINE MPI_Publish_name(service_name, info, port_name, ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: service_name
        INTEGER :: info
        CHARACTER*(*) :: port_name
        INTEGER :: ierr
    END SUBROUTINE MPI_Publish_name
    
    SUBROUTINE MPI_Unpublish_name(service_name, info, port_name, ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: service_name
        INTEGER :: info
        CHARACTER*(*) :: port_name
        INTEGER :: ierr
    END SUBROUTINE MPI_Unpublish_name
    
    SUBROUTINE MPIX_Stream_create(info, stream, ierr)
        IMPLICIT NONE
        INTEGER :: info
        INTEGER :: stream
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stream_create
    
    SUBROUTINE MPIX_Stream_free(stream, ierr)
        IMPLICIT NONE
        INTEGER :: stream
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stream_free
    
    SUBROUTINE MPIX_Stream_comm_create(comm, stream, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: stream
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stream_comm_create
    
    SUBROUTINE MPIX_Stream_comm_create_multiplex(comm, count, array_of_streams, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: count
        INTEGER :: array_of_streams(count)
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stream_comm_create_multiplex
    
    SUBROUTINE MPIX_Comm_get_stream(comm, idx, stream, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: idx
        INTEGER :: stream
        INTEGER :: ierr
    END SUBROUTINE MPIX_Comm_get_stream
    
    SUBROUTINE MPIX_Stream_progress(stream, ierr)
        IMPLICIT NONE
        INTEGER :: stream
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stream_progress
    
    SUBROUTINE MPIX_Start_progress_thread(stream, ierr)
        IMPLICIT NONE
        INTEGER :: stream
        INTEGER :: ierr
    END SUBROUTINE MPIX_Start_progress_thread
    
    SUBROUTINE MPIX_Stop_progress_thread(stream, ierr)
        IMPLICIT NONE
        INTEGER :: stream
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stop_progress_thread
    
    SUBROUTINE MPIX_Stream_send(buf, count, datatype, dest, tag, comm, source_stream_index, &
                                dest_stream_index, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: source_stream_index
        INTEGER :: dest_stream_index
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stream_send
    
    SUBROUTINE MPIX_Stream_isend(buf, count, datatype, dest, tag, comm, source_stream_index, &
                                 dest_stream_index, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: source_stream_index
        INTEGER :: dest_stream_index
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stream_isend
    
    SUBROUTINE MPIX_Stream_recv(buf, count, datatype, source, tag, comm, source_stream_index, &
                                dest_stream_index, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: source_stream_index
        INTEGER :: dest_stream_index
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stream_recv
    
    SUBROUTINE MPIX_Stream_irecv(buf, count, datatype, source, tag, comm, source_stream_index, &
                                 dest_stream_index, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: source_stream_index
        INTEGER :: dest_stream_index
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPIX_Stream_irecv
    
    SUBROUTINE MPIX_Send_enqueue(buf, count, datatype, dest, tag, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPIX_Send_enqueue
    
    SUBROUTINE MPIX_Recv_enqueue(buf, count, datatype, source, tag, comm, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPIX_Recv_enqueue
    
    SUBROUTINE MPIX_Isend_enqueue(buf, count, datatype, dest, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: dest
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPIX_Isend_enqueue
    
    SUBROUTINE MPIX_Irecv_enqueue(buf, count, datatype, source, tag, comm, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: source
        INTEGER :: tag
        INTEGER :: comm
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPIX_Irecv_enqueue
    
    SUBROUTINE MPIX_Wait_enqueue(request, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        INTEGER :: request
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPIX_Wait_enqueue
    
    SUBROUTINE MPIX_Waitall_enqueue(count, array_of_requests, array_of_statuses, ierr)
        IMPLICIT NONE
        INTEGER :: count
        INTEGER :: array_of_requests(count)
        INTEGER :: array_of_statuses(*)
        INTEGER :: ierr
    END SUBROUTINE MPIX_Waitall_enqueue
    
    SUBROUTINE MPIX_Allreduce_enqueue(sendbuf, recvbuf, count, datatype, op, comm, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
        REAL :: sendbuf
        REAL :: recvbuf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: op
        INTEGER :: comm
        INTEGER :: ierr
    END SUBROUTINE MPIX_Allreduce_enqueue
    
    FUNCTION MPI_Wtick() result(res)
        IMPLICIT NONE
        DOUBLE PRECISION :: res
    END FUNCTION MPI_Wtick
    
    FUNCTION PMPI_Wtick() result(res)
        IMPLICIT NONE
        DOUBLE PRECISION :: res
    END FUNCTION PMPI_Wtick
    
    FUNCTION MPI_Wtime() result(res)
        IMPLICIT NONE
        DOUBLE PRECISION :: res
    END FUNCTION MPI_Wtime
    
    FUNCTION PMPI_Wtime() result(res)
        IMPLICIT NONE
        DOUBLE PRECISION :: res
    END FUNCTION PMPI_Wtime
    
    SUBROUTINE MPI_Cart_coords(comm, rank, maxdims, coords, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: rank
        INTEGER :: maxdims
        INTEGER :: coords(maxdims)
        INTEGER :: ierr
    END SUBROUTINE MPI_Cart_coords
    
    SUBROUTINE MPI_Cart_create(comm_old, ndims, dims, periods, reorder, comm_cart, ierr)
        IMPLICIT NONE
        INTEGER :: comm_old
        INTEGER :: ndims
        INTEGER :: dims(ndims)
        LOGICAL :: periods(ndims)
        LOGICAL :: reorder
        INTEGER :: comm_cart
        INTEGER :: ierr
    END SUBROUTINE MPI_Cart_create
    
    SUBROUTINE MPI_Cart_get(comm, maxdims, dims, periods, coords, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: maxdims
        INTEGER :: dims(maxdims)
        LOGICAL :: periods(maxdims)
        INTEGER :: coords(maxdims)
        INTEGER :: ierr
    END SUBROUTINE MPI_Cart_get
    
    SUBROUTINE MPI_Cart_map(comm, ndims, dims, periods, newrank, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: ndims
        INTEGER :: dims(ndims)
        LOGICAL :: periods(ndims)
        INTEGER :: newrank
        INTEGER :: ierr
    END SUBROUTINE MPI_Cart_map
    
    SUBROUTINE MPI_Cart_rank(comm, coords, rank, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: coords(*)
        INTEGER :: rank
        INTEGER :: ierr
    END SUBROUTINE MPI_Cart_rank
    
    SUBROUTINE MPI_Cart_shift(comm, direction, disp, rank_source, rank_dest, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: direction
        INTEGER :: disp
        INTEGER :: rank_source
        INTEGER :: rank_dest
        INTEGER :: ierr
    END SUBROUTINE MPI_Cart_shift
    
    SUBROUTINE MPI_Cart_sub(comm, remain_dims, newcomm, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        LOGICAL :: remain_dims(*)
        INTEGER :: newcomm
        INTEGER :: ierr
    END SUBROUTINE MPI_Cart_sub
    
    SUBROUTINE MPI_Cartdim_get(comm, ndims, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: ndims
        INTEGER :: ierr
    END SUBROUTINE MPI_Cartdim_get
    
    SUBROUTINE MPI_Dims_create(nnodes, ndims, dims, ierr)
        IMPLICIT NONE
        INTEGER :: nnodes
        INTEGER :: ndims
        INTEGER :: dims(ndims)
        INTEGER :: ierr
    END SUBROUTINE MPI_Dims_create
    
    SUBROUTINE MPI_Dist_graph_create(comm_old, n, sources, degrees, destinations, weights, info, &
                                     reorder, comm_dist_graph, ierr)
        IMPLICIT NONE
        INTEGER :: comm_old
        INTEGER :: n
        INTEGER :: sources(n)
        INTEGER :: degrees(n)
        INTEGER :: destinations(*)
        INTEGER :: weights(*)
        INTEGER :: info
        LOGICAL :: reorder
        INTEGER :: comm_dist_graph
        INTEGER :: ierr
    END SUBROUTINE MPI_Dist_graph_create
    
    SUBROUTINE MPI_Dist_graph_create_adjacent(comm_old, indegree, sources, sourceweights, outdegree, &
                                              destinations, destweights, info, reorder, comm_dist_graph, &
                                              ierr)
        IMPLICIT NONE
        INTEGER :: comm_old
        INTEGER :: indegree
        INTEGER :: sources(indegree)
        INTEGER :: sourceweights(*)
        INTEGER :: outdegree
        INTEGER :: destinations(outdegree)
        INTEGER :: destweights(*)
        INTEGER :: info
        LOGICAL :: reorder
        INTEGER :: comm_dist_graph
        INTEGER :: ierr
    END SUBROUTINE MPI_Dist_graph_create_adjacent
    
    SUBROUTINE MPI_Dist_graph_neighbors(comm, maxindegree, sources, sourceweights, maxoutdegree, &
                                        destinations, destweights, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: maxindegree
        INTEGER :: sources(maxindegree)
        INTEGER :: sourceweights(*)
        INTEGER :: maxoutdegree
        INTEGER :: destinations(maxoutdegree)
        INTEGER :: destweights(*)
        INTEGER :: ierr
    END SUBROUTINE MPI_Dist_graph_neighbors
    
    SUBROUTINE MPI_Dist_graph_neighbors_count(comm, indegree, outdegree, weighted, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: indegree
        INTEGER :: outdegree
        LOGICAL :: weighted
        INTEGER :: ierr
    END SUBROUTINE MPI_Dist_graph_neighbors_count
    
    SUBROUTINE MPI_Graph_create(comm_old, nnodes, indx, edges, reorder, comm_graph, ierr)
        IMPLICIT NONE
        INTEGER :: comm_old
        INTEGER :: nnodes
        INTEGER :: indx(nnodes)
        INTEGER :: edges(*)
        LOGICAL :: reorder
        INTEGER :: comm_graph
        INTEGER :: ierr
    END SUBROUTINE MPI_Graph_create
    
    SUBROUTINE MPI_Graph_get(comm, maxindex, maxedges, indx, edges, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: maxindex
        INTEGER :: maxedges
        INTEGER :: indx(maxindex)
        INTEGER :: edges(maxedges)
        INTEGER :: ierr
    END SUBROUTINE MPI_Graph_get
    
    SUBROUTINE MPI_Graph_map(comm, nnodes, indx, edges, newrank, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: nnodes
        INTEGER :: indx(nnodes)
        INTEGER :: edges(*)
        INTEGER :: newrank
        INTEGER :: ierr
    END SUBROUTINE MPI_Graph_map
    
    SUBROUTINE MPI_Graph_neighbors(comm, rank, maxneighbors, neighbors, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: rank
        INTEGER :: maxneighbors
        INTEGER :: neighbors(maxneighbors)
        INTEGER :: ierr
    END SUBROUTINE MPI_Graph_neighbors
    
    SUBROUTINE MPI_Graph_neighbors_count(comm, rank, nneighbors, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: rank
        INTEGER :: nneighbors
        INTEGER :: ierr
    END SUBROUTINE MPI_Graph_neighbors_count
    
    SUBROUTINE MPI_Graphdims_get(comm, nnodes, nedges, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: nnodes
        INTEGER :: nedges
        INTEGER :: ierr
    END SUBROUTINE MPI_Graphdims_get
    
    SUBROUTINE MPI_Topo_test(comm, status, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: status
        INTEGER :: ierr
    END SUBROUTINE MPI_Topo_test
    
    SUBROUTINE MPI_File_close(fh, ierr)
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER :: ierr
    END SUBROUTINE MPI_File_close
    
    SUBROUTINE MPI_File_delete(filename, info, ierr)
        IMPLICIT NONE
        CHARACTER*(*) :: filename
        INTEGER :: info
        INTEGER :: ierr
    END SUBROUTINE MPI_File_delete
    
    SUBROUTINE MPI_File_get_amode(fh, amode, ierr)
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER :: amode
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_amode
    
    SUBROUTINE MPI_File_get_atomicity(fh, flag, ierr)
        IMPLICIT NONE
        INTEGER :: fh
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_atomicity
    
    SUBROUTINE MPI_File_get_byte_offset(fh, offset, disp, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        INTEGER(KIND=MPI_OFFSET_KIND) :: disp
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_byte_offset
    
    SUBROUTINE MPI_File_get_group(fh, group, ierr)
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER :: group
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_group
    
    SUBROUTINE MPI_File_get_info(fh, info_used, ierr)
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER :: info_used
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_info
    
    SUBROUTINE MPI_File_get_position(fh, offset, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_position
    
    SUBROUTINE MPI_File_get_position_shared(fh, offset, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_position_shared
    
    SUBROUTINE MPI_File_get_size(fh, size, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_size
    
    SUBROUTINE MPI_File_get_type_extent(fh, datatype, extent, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER :: datatype
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extent
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_type_extent
    
    SUBROUTINE MPI_File_get_view(fh, disp, etype, filetype, datarep, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: disp
        INTEGER :: etype
        INTEGER :: filetype
        CHARACTER*(*) :: datarep
        INTEGER :: ierr
    END SUBROUTINE MPI_File_get_view
    
    SUBROUTINE MPI_File_iread(fh, buf, count, datatype, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iread
    
    SUBROUTINE MPI_File_iread_all(fh, buf, count, datatype, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iread_all
    
    SUBROUTINE MPI_File_iread_at(fh, offset, buf, count, datatype, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iread_at
    
    SUBROUTINE MPI_File_iread_at_all(fh, offset, buf, count, datatype, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iread_at_all
    
    SUBROUTINE MPI_File_iread_shared(fh, buf, count, datatype, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iread_shared
    
    SUBROUTINE MPI_File_iwrite(fh, buf, count, datatype, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iwrite
    
    SUBROUTINE MPI_File_iwrite_all(fh, buf, count, datatype, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iwrite_all
    
    SUBROUTINE MPI_File_iwrite_at(fh, offset, buf, count, datatype, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iwrite_at
    
    SUBROUTINE MPI_File_iwrite_at_all(fh, offset, buf, count, datatype, request, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iwrite_at_all
    
    SUBROUTINE MPI_File_iwrite_shared(fh, buf, count, datatype, request, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: request
        INTEGER :: ierr
    END SUBROUTINE MPI_File_iwrite_shared
    
    SUBROUTINE MPI_File_open(comm, filename, amode, info, fh, ierr)
        IMPLICIT NONE
        INTEGER :: comm
        CHARACTER*(*) :: filename
        INTEGER :: amode
        INTEGER :: info
        INTEGER :: fh
        INTEGER :: ierr
    END SUBROUTINE MPI_File_open
    
    SUBROUTINE MPI_File_preallocate(fh, size, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_File_preallocate
    
    SUBROUTINE MPI_File_read(fh, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read
    
    SUBROUTINE MPI_File_read_all(fh, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_all
    
    SUBROUTINE MPI_File_read_all_begin(fh, buf, count, datatype, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_all_begin
    
    SUBROUTINE MPI_File_read_all_end(fh, buf, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_all_end
    
    SUBROUTINE MPI_File_read_at(fh, offset, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND, MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_at
    
    SUBROUTINE MPI_File_read_at_all(fh, offset, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND, MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_at_all
    
    SUBROUTINE MPI_File_read_at_all_begin(fh, offset, buf, count, datatype, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_at_all_begin
    
    SUBROUTINE MPI_File_read_at_all_end(fh, buf, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_at_all_end
    
    SUBROUTINE MPI_File_read_ordered(fh, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_ordered
    
    SUBROUTINE MPI_File_read_ordered_begin(fh, buf, count, datatype, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_ordered_begin
    
    SUBROUTINE MPI_File_read_ordered_end(fh, buf, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_ordered_end
    
    SUBROUTINE MPI_File_read_shared(fh, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_read_shared
    
    SUBROUTINE MPI_File_seek(fh, offset, whence, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        INTEGER :: whence
        INTEGER :: ierr
    END SUBROUTINE MPI_File_seek
    
    SUBROUTINE MPI_File_seek_shared(fh, offset, whence, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        INTEGER :: whence
        INTEGER :: ierr
    END SUBROUTINE MPI_File_seek_shared
    
    SUBROUTINE MPI_File_set_atomicity(fh, flag, ierr)
        IMPLICIT NONE
        INTEGER :: fh
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_File_set_atomicity
    
    SUBROUTINE MPI_File_set_info(fh, info, ierr)
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER :: info
        INTEGER :: ierr
    END SUBROUTINE MPI_File_set_info
    
    SUBROUTINE MPI_File_set_size(fh, size, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: size
        INTEGER :: ierr
    END SUBROUTINE MPI_File_set_size
    
    SUBROUTINE MPI_File_set_view(fh, disp, etype, filetype, datarep, info, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: disp
        INTEGER :: etype
        INTEGER :: filetype
        CHARACTER*(*) :: datarep
        INTEGER :: info
        INTEGER :: ierr
    END SUBROUTINE MPI_File_set_view
    
    SUBROUTINE MPI_File_sync(fh, ierr)
        IMPLICIT NONE
        INTEGER :: fh
        INTEGER :: ierr
    END SUBROUTINE MPI_File_sync
    
    SUBROUTINE MPI_File_write(fh, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write
    
    SUBROUTINE MPI_File_write_all(fh, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_all
    
    SUBROUTINE MPI_File_write_all_begin(fh, buf, count, datatype, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_all_begin
    
    SUBROUTINE MPI_File_write_all_end(fh, buf, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_all_end
    
    SUBROUTINE MPI_File_write_at(fh, offset, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND, MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_at
    
    SUBROUTINE MPI_File_write_at_all(fh, offset, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND, MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_at_all
    
    SUBROUTINE MPI_File_write_at_all_begin(fh, offset, buf, count, datatype, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        INTEGER(KIND=MPI_OFFSET_KIND) :: offset
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_at_all_begin
    
    SUBROUTINE MPI_File_write_at_all_end(fh, buf, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_at_all_end
    
    SUBROUTINE MPI_File_write_ordered(fh, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_ordered
    
    SUBROUTINE MPI_File_write_ordered_begin(fh, buf, count, datatype, ierr)
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_ordered_begin
    
    SUBROUTINE MPI_File_write_ordered_end(fh, buf, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_ordered_end
    
    SUBROUTINE MPI_File_write_shared(fh, buf, count, datatype, status, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_STATUS_SIZE
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: buf
        INTEGER :: fh
        REAL :: buf
        INTEGER :: count
        INTEGER :: datatype
        INTEGER :: status(MPI_STATUS_SIZE)
        INTEGER :: ierr
    END SUBROUTINE MPI_File_write_shared
    
    SUBROUTINE MPI_Register_datarep(datarep, read_conversion_fn, write_conversion_fn, &
                                    dtype_file_extent_fn, extra_state, ierr)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        CHARACTER*(*) :: datarep
        EXTERNAL :: read_conversion_fn
        EXTERNAL :: write_conversion_fn
        EXTERNAL :: dtype_file_extent_fn
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER :: ierr
    END SUBROUTINE MPI_Register_datarep
    
    SUBROUTINE MPI_Type_create_f90_integer(r, newtype, ierr)
        IMPLICIT NONE
        INTEGER :: r
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_f90_integer
    
    SUBROUTINE MPI_Type_create_f90_real(p, r, newtype, ierr)
        IMPLICIT NONE
        INTEGER :: p
        INTEGER :: r
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_f90_real
    
    SUBROUTINE MPI_Type_create_f90_complex(p, r, newtype, ierr)
        IMPLICIT NONE
        INTEGER :: p
        INTEGER :: r
        INTEGER :: newtype
        INTEGER :: ierr
    END SUBROUTINE MPI_Type_create_f90_complex
    
    SUBROUTINE MPI_DUP_FN(oldcomm, keyval, extra_state, attribute_val_in, attribute_val_out, flag, &
                          ierr)
        IMPLICIT NONE
        INTEGER :: oldcomm
        INTEGER :: keyval
        INTEGER :: extra_state
        INTEGER :: attribute_val_in
        INTEGER :: attribute_val_out
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_DUP_FN
    
    SUBROUTINE MPI_NULL_COPY_FN(oldcomm, keyval, extra_state, attribute_val_in, attribute_val_out, flag, &
                                ierr)
        IMPLICIT NONE
        INTEGER :: oldcomm
        INTEGER :: keyval
        INTEGER :: extra_state
        INTEGER :: attribute_val_in
        INTEGER :: attribute_val_out
        LOGICAL :: flag
        INTEGER :: ierr
    END SUBROUTINE MPI_NULL_COPY_FN
    
    SUBROUTINE MPI_NULL_DELETE_FN(comm, keyval, attribute_val, extra_state, ierror)
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: keyval
        INTEGER :: attribute_val
        INTEGER :: extra_state
        INTEGER :: ierror
    END SUBROUTINE MPI_NULL_DELETE_FN
    
    SUBROUTINE MPI_COMM_DUP_FN(oldcomm, comm_keyval, extra_state, attribute_val_in, attribute_val_out, &
                               flag, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: oldcomm
        INTEGER :: comm_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_in
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_out
        LOGICAL :: flag
        INTEGER :: ierror
    END SUBROUTINE MPI_COMM_DUP_FN
    
    SUBROUTINE MPI_COMM_NULL_COPY_FN(oldcomm, comm_keyval, extra_state, attribute_val_in, &
                                     attribute_val_out, flag, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: oldcomm
        INTEGER :: comm_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_in
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_out
        LOGICAL :: flag
        INTEGER :: ierror
    END SUBROUTINE MPI_COMM_NULL_COPY_FN
    
    SUBROUTINE MPI_COMM_NULL_DELETE_FN(comm, comm_keyval, attribute_val, extra_state, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: comm
        INTEGER :: comm_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER :: ierror
    END SUBROUTINE MPI_COMM_NULL_DELETE_FN
    
    SUBROUTINE MPI_TYPE_DUP_FN(oldtype, type_keyval, extra_state, attribute_val_in, attribute_val_out, &
                               flag, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: oldtype
        INTEGER :: type_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_in
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_out
        LOGICAL :: flag
        INTEGER :: ierror
    END SUBROUTINE MPI_TYPE_DUP_FN
    
    SUBROUTINE MPI_TYPE_NULL_COPY_FN(oldtype, type_keyval, extra_state, attribute_val_in, &
                                     attribute_val_out, flag, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: oldtype
        INTEGER :: type_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_in
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_out
        LOGICAL :: flag
        INTEGER :: ierror
    END SUBROUTINE MPI_TYPE_NULL_COPY_FN
    
    SUBROUTINE MPI_TYPE_NULL_DELETE_FN(datatype, type_keyval, attribute_val, extra_state, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: datatype
        INTEGER :: type_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER :: ierror
    END SUBROUTINE MPI_TYPE_NULL_DELETE_FN
    
    SUBROUTINE MPI_WIN_DUP_FN(oldwin, win_keyval, extra_state, attribute_val_in, attribute_val_out, &
                              flag, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: oldwin
        INTEGER :: win_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_in
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_out
        LOGICAL :: flag
        INTEGER :: ierror
    END SUBROUTINE MPI_WIN_DUP_FN
    
    SUBROUTINE MPI_WIN_NULL_COPY_FN(oldwin, win_keyval, extra_state, attribute_val_in, &
                                    attribute_val_out, flag, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: oldwin
        INTEGER :: win_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_in
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val_out
        LOGICAL :: flag
        INTEGER :: ierror
    END SUBROUTINE MPI_WIN_NULL_COPY_FN
    
    SUBROUTINE MPI_WIN_NULL_DELETE_FN(win, win_keyval, attribute_val, extra_state, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_ADDRESS_KIND
        IMPLICIT NONE
        INTEGER :: win
        INTEGER :: win_keyval
        INTEGER(KIND=MPI_ADDRESS_KIND) :: attribute_val
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER :: ierror
    END SUBROUTINE MPI_WIN_NULL_DELETE_FN
    
    SUBROUTINE MPI_CONVERSION_FN_NULL(userbuf, datatype, count, filebuf, position, extra_state, ierror)
        USE MPI_CONSTANTS, ONLY: MPI_OFFSET_KIND, MPI_ADDRESS_KIND
        IMPLICIT NONE
        !GCC$ ATTRIBUTES NO_ARG_CHECK :: userbuf, filebuf
        REAL :: userbuf
        INTEGER :: datatype
        INTEGER :: count
        REAL :: filebuf
        INTEGER(KIND=MPI_OFFSET_KIND) :: position
        INTEGER(KIND=MPI_ADDRESS_KIND) :: extra_state
        INTEGER :: ierror
    END SUBROUTINE MPI_CONVERSION_FN_NULL
    
    END INTERFACE
END MODULE mpi_base
