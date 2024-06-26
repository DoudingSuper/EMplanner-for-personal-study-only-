/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

/* This file contains the prototypes for routines that are used with
   "external" modules such as ROMIO.  These allow the different packages to
   hide their internal datatypes from one another */

#ifndef MPIEXT_H_INCLUDED
#define MPIEXT_H_INCLUDED

#include <stdarg.h>

/* This routine, given an MPI_Errhandler (from a file), returns
   a pointer to the user-supplied error function.  The last argument
   is set to an integer indicating that the function is MPI_ERRORS_RETURN
   (value == 1), MPI_ERRORS_ARE_FATAL (value == 0), a valid user-function
   (value == 2), or a valid user-function that is a C++ routine (value == 3)

   This routine is implemented in mpich/src/mpi/errhan/file_set_errhandler.c
*/
void MPIR_Get_file_error_routine(MPI_Errhandler, void (**)(MPI_File *, int *, ...), int *);

/* Invoke the C++ error handler (this invokes a special C++ routine that
 in turn calls the provided function.  That special routine also
 resets the errorcode to MPI_SUCCESS to prevent the MPICH C++ error handling
 code from throwing an exception when the user routine returns.
*/
int MPIR_File_call_cxx_errhandler(MPI_File *, int *, void (*)(MPI_File *, int *, ...));
/*
   These routines provide access to the MPI_Errhandler field within the
   ROMIO MPI_File structure
 */
int MPIR_ROMIO_Get_file_errhand(MPI_File, MPI_Errhandler *);
int MPIR_ROMIO_Set_file_errhand(MPI_File, MPI_Errhandler);

/* FIXME: This routine is also defined in adio.h */
int MPIO_Err_return_file(MPI_File, int);

int MPIR_Err_create_code_valist(int, int, const char[], int, int,
                                const char[], const char[], va_list);
int MPIR_Err_is_fatal(int);

void MPIR_Get_file_error_routine(MPI_Errhandler, void (**)(MPI_File *, int *, ...), int *);
int MPIR_File_call_cxx_errhandler(MPI_File *, int *, void (*)(MPI_File *, int *, ...));

typedef int (*MPIR_Err_get_class_string_func_t) (int error, char *str, int length);
void MPIR_Err_get_string(int, char *, int, MPIR_Err_get_class_string_func_t);

struct MPIR_Comm;
int MPIR_Abort(MPI_Comm comm, int mpi_errno, int exit_code, const char *error_msg);

int MPIR_Ext_assert_fail(const char *cond, const char *file_name, int line_num);

#if (!defined(NDEBUG) && (1))
#define MPIR_Ext_assert(a_)                                \
    do {                                                   \
        if (!(a_)) {                                       \
            MPIR_Ext_assert_fail(#a_, __FILE__, __LINE__); \
        }                                                  \
    } while (0)
#else
#define MPIR_Ext_assert(a_) do {} while (0)
#endif

extern int MPIR_Ext_dbg_romio_terse_enabled;
extern int MPIR_Ext_dbg_romio_typical_enabled;
extern int MPIR_Ext_dbg_romio_verbose_enabled;

/* to be called early by ROMIO's initialization process in order to setup init-time
 * glue code that cannot be initialized statically */
int MPIR_Ext_init(void);

void MPIR_Ext_mutex_init(void);
void MPIR_Ext_mutex_finalize(void);

void MPIR_Ext_cs_enter(void);
void MPIR_Ext_cs_exit(void);
void MPIR_Ext_cs_yield(void);

/* to facilitate error checking */
int MPIR_Ext_datatype_iscommitted(MPI_Datatype datatype);

/* make comm split based on access to a common file system easier */
int MPIR_Get_node_id(MPI_Comm comm, int rank, int *id);

#ifdef ENABLE_QMPI
extern void (**MPIR_QMPI_pointers) (void);
extern void **MPIR_QMPI_storage;
extern void (**MPIR_QMPI_tool_init_callbacks) (int);
extern int MPIR_QMPI_num_tools;
extern char **MPIR_QMPI_tool_names;
extern void (**MPIR_QMPI_first_fn_ptrs) (void);
extern int *MPIR_QMPI_first_tool_ids;
int MPII_qmpi_stash_first_tools(void);
int MPII_qmpi_pre_init(void);
int MPII_qmpi_init(void);
void MPII_qmpi_teardown(void);
#endif /* ENABLE_QMPI */

#endif
