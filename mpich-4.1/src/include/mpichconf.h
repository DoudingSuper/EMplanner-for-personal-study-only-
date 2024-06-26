/* src/include/mpichconf.h.  Generated from mpichconf.h.in by configure.  */
/* src/include/mpichconf.h.in.  Generated from configure.ac by autoheader.  */

/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */
#ifndef MPICHCONF_H_INCLUDED
#define MPICHCONF_H_INCLUDED


/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* The normal alignment of `bool', in bytes. */
#define ALIGNOF_BOOL 0

/* The normal alignment of `char', in bytes. */
#define ALIGNOF_CHAR 1

/* The normal alignment of `double', in bytes. */
#define ALIGNOF_DOUBLE 8

/* The normal alignment of `float', in bytes. */
#define ALIGNOF_FLOAT 4

/* The normal alignment of `int', in bytes. */
#define ALIGNOF_INT 4

/* The normal alignment of `int16_t', in bytes. */
#define ALIGNOF_INT16_T 2

/* The normal alignment of `int32_t', in bytes. */
#define ALIGNOF_INT32_T 4

/* The normal alignment of `int64_t', in bytes. */
#define ALIGNOF_INT64_T 8

/* The normal alignment of `int8_t', in bytes. */
#define ALIGNOF_INT8_T 1

/* The normal alignment of `long', in bytes. */
#define ALIGNOF_LONG 8

/* The normal alignment of `long double', in bytes. */
#define ALIGNOF_LONG_DOUBLE 16

/* The normal alignment of `long long', in bytes. */
#define ALIGNOF_LONG_LONG 8

/* The normal alignment of `max_align_t', in bytes. */
#define ALIGNOF_MAX_ALIGN_T 0

/* The normal alignment of `short', in bytes. */
#define ALIGNOF_SHORT 2

/* The normal alignment of `wchar_t', in bytes. */
#define ALIGNOF_WCHAR_T 4

/* Define the number of CH3_RANK_BITS */
/* #undef CH3_RANK_BITS */

/* Define the number of rank bits used in UCX */
/* #undef CH4_UCX_RANKBITS */

/* Define the search path for machines files */
/* #undef DEFAULT_MACHINES_PATH */

/* Define the default remote shell program to use */
/* #undef DEFAULT_REMOTE_SHELL */

/* Define to workaround interprocess mutex issue on FreeBSD */
/* #undef DELAY_SHM_MUTEX_DESTROY */

/* Define to enable shared-memory collectives */
/* #undef ENABLED_SHM_COLLECTIVES */

/* Application checkpointing enabled */
/* #undef ENABLE_CHECKPOINTING */

/* define to add per-vc function pointers to override send and recv functions
   */
/* #undef ENABLE_COMM_OVERRIDES */

/* Define to skip initializing builtin world comm during MPI_Session_init */
#define ENABLE_LOCAL_SESSION_INIT 1

/* Define to disable shared-memory communication */
/* #undef ENABLE_NO_LOCAL */

/* Define to 1 to enable getdims-related MPI_T performance variables */
#define ENABLE_PVAR_DIMS 0

/* Define to 1 to enable message count transmitted through multiple NICs MPI_T
   performance variables */
#define ENABLE_PVAR_MULTINIC 0

/* Define to 1 to enable nemesis-related MPI_T performance variables */
#define ENABLE_PVAR_NEM 0

/* Define to 1 to enable message receive queue-related MPI_T performance
   variables */
#define ENABLE_PVAR_RECVQ 0

/* Define to 1 to enable rma-related MPI_T performance variables */
#define ENABLE_PVAR_RMA 0

/* Define if QMPI enabled */
/* #undef ENABLE_QMPI */

/* The value of false in Fortran */
#define F77_FALSE_VALUE 0

/* Fortran names are lowercase with no trailing underscore */
/* #undef F77_NAME_LOWER */

/* Fortran names are lowercase with two trailing underscores */
/* #undef F77_NAME_LOWER_2USCORE */

/* Fortran names are lowercase with two trailing underscores in stdcall */
/* #undef F77_NAME_LOWER_2USCORE_STDCALL */

/* Fortran names are lowercase with no trailing underscore in stdcall */
/* #undef F77_NAME_LOWER_STDCALL */

/* Fortran names are lowercase with one trailing underscore */
#define F77_NAME_LOWER_USCORE 1

/* Fortran names are lowercase with one trailing underscore in stdcall */
/* #undef F77_NAME_LOWER_USCORE_STDCALL */

/* Fortran names preserve the original case */
/* #undef F77_NAME_MIXED */

/* Fortran names preserve the original case in stdcall */
/* #undef F77_NAME_MIXED_STDCALL */

/* Fortran names preserve the original case with one trailing underscore */
/* #undef F77_NAME_MIXED_USCORE */

/* Fortran names preserve the original case with one trailing underscore in
   stdcall */
/* #undef F77_NAME_MIXED_USCORE_STDCALL */

/* Fortran names are uppercase */
/* #undef F77_NAME_UPPER */

/* Fortran names are uppercase in stdcall */
/* #undef F77_NAME_UPPER_STDCALL */

/* The value of true in Fortran */
#define F77_TRUE_VALUE 1

/* Define if we know the value of Fortran true and false */
#define F77_TRUE_VALUE_SET 1

/* Define FALSE */
#define FALSE 0

/* Directory to use in namepub */
/* #undef FILE_NAMEPUB_BASEDIR */

/* Define if addresses are a different size than Fortran integers */
#define HAVE_AINT_DIFFERENT_THAN_FINT 1

/* Define if addresses are larger than Fortran integers */
#define HAVE_AINT_LARGER_THAN_FINT 1

/* Define to 1 if you have the `alarm' function. */
#define HAVE_ALARM 1

/* Define if int32_t works with any alignment */
/* #undef HAVE_ANY_INT32_T_ALIGNMENT */

/* Define if int64_t works with any alignment */
/* #undef HAVE_ANY_INT64_T_ALIGNMENT */

/* Define to 1 if you have the <arpa/inet.h> header file. */
/* #undef HAVE_ARPA_INET_H */

/* Define to 1 if you have the <assert.h> header file. */
#define HAVE_ASSERT_H 1

/* Define to 1 if you have the `bindprocessor' function. */
/* #undef HAVE_BINDPROCESSOR */

/* Define to 1 if the compiler supports __builtin_expect. */
#define HAVE_BUILTIN_EXPECT 1

/* Define if C11 _Static_assert is supported. */
#define HAVE_C11__STATIC_ASSERT 1

/* Define to 1 if you have the `CFUUIDCreate' function. */
/* #undef HAVE_CFUUIDCREATE */

/* Define if debugger support is included for CH4 */
/* #undef HAVE_CH4_DEBUGGER_SUPPORT */

/* OFI netmod is built */
#define HAVE_CH4_NETMOD_OFI 1

/* UCX netmod is built */
/* #undef HAVE_CH4_NETMOD_UCX */

/* IQUEUE submodule is built */
#define HAVE_CH4_SHM_EAGER_IQUEUE 1

/* STUB submodule is built */
/* #undef HAVE_CH4_SHM_EAGER_STUB */

/* Define to 1 if you have the <complex.h> header file. */
#define HAVE_COMPLEX_H 1

/* Define if CPU_SET and CPU_ZERO defined */
#define HAVE_CPU_SET_MACROS 1

/* Define if cpu_set_t is defined in sched.h */
#define HAVE_CPU_SET_T 1

/* Define to 1 if you have the <ctype.h> header file. */
/* #undef HAVE_CTYPE_H */

/* Define if C++ is supported */
#define HAVE_CXX_BINDING 1

/* Define if C++ supports bool types */
#define HAVE_CXX_BOOL 1

/* Define if C++ supports complex types */
#define HAVE_CXX_COMPLEX 1

/* define if the compiler supports exceptions */
#define HAVE_CXX_EXCEPTIONS /**/

/* Define if C++ supports long double complex */
#define HAVE_CXX_LONG_DOUBLE_COMPLEX 1

/* Define if multiple __attribute__((alias)) are supported */
/* #undef HAVE_C_MULTI_ATTR_ALIAS */

/* Define if debugger support is included */
/* #undef HAVE_DEBUGGER_SUPPORT */

/* Define to 1 if you have the declaration of `strerror_r', and to 0 if you
   don't. */
#define HAVE_DECL_STRERROR_R 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if the system has the type `double _Complex'. */
#define HAVE_DOUBLE__COMPLEX 1

/* Define to 1 if you have the <endian.h> header file. */
#define HAVE_ENDIAN_H 1

/* Define to 1 if you have the <errno.h> header file. */
#define HAVE_ERRNO_H 1

/* Define to enable error checking */
#define HAVE_ERROR_CHECKING MPID_ERROR_LEVEL_ALL

/* Define to enable extended context id bit space */
/* #undef HAVE_EXTENDED_CONTEXT_BITS */

/* Define if environ extern is available */
/* #undef HAVE_EXTERN_ENVIRON */

/* Define to 1 if we have Fortran 2008 binding */
#define HAVE_F08_BINDING 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define if Fortran 90 type routines available */
#define HAVE_FC_TYPE_ROUTINES 1

/* Define if Fortran integer are the same size as C ints */
#define HAVE_FINT_IS_INT 1

/* Define if __float128 is supported */
#define HAVE_FLOAT128 1

/* Define if _Float16 is supported */
/* #undef HAVE_FLOAT16 */

/* Define to 1 if the system has the type `float _Complex'. */
#define HAVE_FLOAT__COMPLEX 1

/* Define if Fortran is supported */
#define HAVE_FORTRAN_BINDING 1

/* Define if GNU __attribute__ is supported */
#define HAVE_GCC_ATTRIBUTE 1

/* Define to 1 if you have the `gethostname' function. */
#define HAVE_GETHOSTNAME 1

/* Define to 1 if you have the `getsid' function. */
/* #undef HAVE_GETSID */

/* Define if building hcoll */
/* #undef HAVE_HCOLL */

/* Define if hwloc is available */
#define HAVE_HWLOC 1

/* Define to 1 if you have the `inet_pton' function. */
/* #undef HAVE_INET_PTON */

/* Define if int16_t is supported by the C compiler */
#define HAVE_INT16_T 1

/* Define if int32_t is supported by the C compiler */
#define HAVE_INT32_T 1

/* Define if int64_t is supported by the C compiler */
#define HAVE_INT64_T 1

/* Define if int8_t is supported by the C compiler */
#define HAVE_INT8_T 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if struct iovec defined in sys/uio.h */
/* #undef HAVE_IOVEC_DEFINITION */

/* Define to 1 if you have the `isatty' function. */
/* #undef HAVE_ISATTY */

/* Define to 1 if you have the `cr' library (-lcr). */
/* #undef HAVE_LIBCR */

/* Define to 1 if you have the `fabric' library (-lfabric). */
/* #undef HAVE_LIBFABRIC */

/* Define if libfabric library has nic field in fi_info struct */
#define HAVE_LIBFABRIC_NIC 1

/* Define to 1 if you have the `pmi' library (-lpmi). */
/* #undef HAVE_LIBPMI */

/* Define to 1 if you have the `ucp' library (-lucp). */
/* #undef HAVE_LIBUCP */

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define if long double is supported */
#define HAVE_LONG_DOUBLE 1

/* Define to 1 if the system has the type `long double _Complex'. */
#define HAVE_LONG_DOUBLE__COMPLEX 1

/* Define if long long allowed */
#define HAVE_LONG_LONG 1

/* Define if long long is supported */
#define HAVE_LONG_LONG_INT 1

/* Define if C99-style variable argument list macro functionality */
#define HAVE_MACRO_VA_ARGS 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define if 256 bit streaming memcpy is available */
/* #undef HAVE_MM256_STREAM_SI256 */

/* Define so that we can test whether the mpichconf.h file has been included
   */
#define HAVE_MPICHCONF 1

/* Define if MPI_T Events are enabled */
/* #undef HAVE_MPIT_EVENTS */

/* Define if the Fortran init code for MPI works from C programs without
   special libraries */
#define HAVE_MPI_F_INIT_WORKS_WITH_C 1

/* Define if multiple weak symbols may be defined */
#define HAVE_MULTIPLE_PRAGMA_WEAK 1

/* Define if a name publishing service is available */
#define HAVE_NAMEPUB_SERVICE 1

/* define if the compiler implements namespaces */
#define HAVE_NAMESPACES /**/

/* define if the compiler implements namespace std */
#define HAVE_NAMESPACE_STD /**/

/* Define to 1 if you have the <netdb.h> header file. */
/* #undef HAVE_NETDB_H */

/* Define if netinet/in.h exists */
/* #undef HAVE_NETINET_IN_H */

/* Define to 1 if you have the <netinet/tcp.h> header file. */
/* #undef HAVE_NETINET_TCP_H */

/* Define if netloc is available */
/* #undef HAVE_NETLOC */

/* Define to 1 if you have the <net/if.h> header file. */
/* #undef HAVE_NET_IF_H */

/* Define if the Fortran types are not available in C */
/* #undef HAVE_NO_FORTRAN_MPI_TYPES_IN_C */

/* Define if the OSX thread affinity policy macros defined */
/* #undef HAVE_OSX_THREAD_AFFINITY */

/* Define to 1 if you have the <poll.h> header file. */
/* #undef HAVE_POLL_H */

/* Cray style weak pragma */
/* #undef HAVE_PRAGMA_CRI_DUP */

/* HP style weak pragma */
/* #undef HAVE_PRAGMA_HP_SEC_DEF */

/* Supports weak pragma */
#define HAVE_PRAGMA_WEAK 1

/* Define to 1 if you have the `ptrace' function. */
/* #undef HAVE_PTRACE */

/* Define if ptrace parameters available */
/* #undef HAVE_PTRACE_CONT */

/* Define to 1 if you have the `putenv' function. */
#define HAVE_PUTENV 1

/* Define to 1 if you have the `qsort' function. */
#define HAVE_QSORT 1

/* Define to 1 if you have the <random.h> header file. */
/* #undef HAVE_RANDOM_H */

/* Define to 1 if you have the `random_r' function. */
#define HAVE_RANDOM_R 1

/* Define if ROMIO is enabled */
#define HAVE_ROMIO 1

/* Define to 1 if you have the `sched_getaffinity' function. */
#define HAVE_SCHED_GETAFFINITY 1

/* Define to 1 if you have the <sched.h> header file. */
/* #undef HAVE_SCHED_H */

/* Define to 1 if you have the `sched_setaffinity' function. */
#define HAVE_SCHED_SETAFFINITY 1

/* Define to 1 if you have the `select' function. */
/* #undef HAVE_SELECT */

/* Define to 1 if you have the `setitimer' function. */
#define HAVE_SETITIMER 1

/* Define to 1 if you have the `setsid' function. */
/* #undef HAVE_SETSID */

/* Define to 1 if you have the `sigaction' function. */
/* #undef HAVE_SIGACTION */

/* Define to 1 if you have the `signal' function. */
#define HAVE_SIGNAL 1

/* Define to 1 if you have the <signal.h> header file. */
#define HAVE_SIGNAL_H 1

/* Define to 1 if you have the `sigset' function. */
/* #undef HAVE_SIGSET */

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define if socklen_t is available */
/* #undef HAVE_SOCKLEN_T */

/* Define to 1 if you have the <stdarg.h> header file. */
#define HAVE_STDARG_H 1

/* Define to 1 if you have the <stdbool.h> header file. */
#define HAVE_STDBOOL_H 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strdup' function. */
#define HAVE_STRDUP 1

/* Define to 1 if you have the `strerror' function. */
#define HAVE_STRERROR 1

/* Define to 1 if you have the `strerror_r' function. */
#define HAVE_STRERROR_R 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strncasecmp' function. */
#define HAVE_STRNCASECMP 1

/* Define to 1 if you have the `strsignal' function. */
/* #undef HAVE_STRSIGNAL */

/* Define to 1 if the system has the type `struct random_data'. */
#define HAVE_STRUCT_RANDOM_DATA 1

/* Define if sys/bitypes.h exists */
#define HAVE_SYS_BITYPES_H 1

/* Define to 1 if you have the <sys/ioctl.h> header file. */
/* #undef HAVE_SYS_IOCTL_H */

/* Define to 1 if you have the <sys/ipc.h> header file. */
/* #undef HAVE_SYS_IPC_H */

/* Define to 1 if you have the <sys/mman.h> header file. */
#define HAVE_SYS_MMAN_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/poll.h> header file. */
/* #undef HAVE_SYS_POLL_H */

/* Define to 1 if you have the <sys/ptrace.h> header file. */
/* #undef HAVE_SYS_PTRACE_H */

/* Define to 1 if you have the <sys/select.h> header file. */
/* #undef HAVE_SYS_SELECT_H */

/* Define to 1 if you have the <sys/shm.h> header file. */
/* #undef HAVE_SYS_SHM_H */

/* Define to 1 if you have the <sys/socket.h> header file. */
#define HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/sockio.h> header file. */
/* #undef HAVE_SYS_SOCKIO_H */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <sys/uio.h> header file. */
#define HAVE_SYS_UIO_H 1

/* Define to enable tag error bits */
#define HAVE_TAG_ERROR_BITS 1

/* Define to 1 if you have the `thread_policy_set' function. */
/* #undef HAVE_THREAD_POLICY_SET */

/* Define to 1 if you have the `time' function. */
/* #undef HAVE_TIME */

/* Define to 1 if you have the <time.h> header file. */
/* #undef HAVE_TIME_H */

/* Define if uint16_t is supported by the C compiler */
#define HAVE_UINT16_T 1

/* Define if uint32_t is supported by the C compiler */
#define HAVE_UINT32_T 1

/* Define if uint64_t is supported by the C compiler */
#define HAVE_UINT64_T 1

/* Define if uint8_t is supported by the C compiler */
#define HAVE_UINT8_T 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `unsetenv' function. */
/* #undef HAVE_UNSETENV */

/* Define to 1 if you have the `usleep' function. */
/* #undef HAVE_USLEEP */

/* Define to 1 if you have the `uuid_generate' function. */
/* #undef HAVE_UUID_GENERATE */

/* Define to 1 if you have the <uuid/uuid.h> header file. */
/* #undef HAVE_UUID_UUID_H */

/* Whether C compiler supports symbol visibility or not */
#define HAVE_VISIBILITY 1

/* Define to 1 if you have the `vsnprintf' function. */
#define HAVE_VSNPRINTF 1

/* Define to 1 if you have the `vsprintf' function. */
#define HAVE_VSPRINTF 1

/* Define to 1 if you have the <wait.h> header file. */
/* #undef HAVE_WAIT_H */

/* Attribute style weak pragma */
#define HAVE_WEAK_ATTRIBUTE 1

/* Define to 1 if the system has the type `_Bool'. */
#define HAVE__BOOL 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Controls byte alignment of structures (for aligning allocated structures)
   */
#define MAX_ALIGNMENT 16

/* Datatype engine */
#define MPICH_DATATYPE_ENGINE MPICH_DATATYPE_ENGINE_YAKSA

/* Define to enable checking of handles still allocated at MPI_Finalize */
/* #undef MPICH_DEBUG_HANDLEALLOC */

/* Define to enable handle checking */
/* #undef MPICH_DEBUG_HANDLES */

/* Define if each function exit should confirm memory arena correctness */
/* #undef MPICH_DEBUG_MEMARENA */

/* Define to enable preinitialization of memory used by structures and unions
   */
/* #undef MPICH_DEBUG_MEMINIT */

/* Define to enable mutex debugging */
/* #undef MPICH_DEBUG_MUTEX */

/* define to enable error messages */
#define MPICH_ERROR_MSG_LEVEL MPICH_ERROR_MSG__ALL

/* Define as the name of the debugger support library */
/* #undef MPICH_INFODLL_LOC */

/* MPICH is configured to require thread safety */
#define MPICH_IS_THREADED 1

/* Method used to implement atomic updates and access */
#define MPICH_THREAD_GRANULARITY MPICH_THREAD_GRANULARITY__VCI

/* Level of thread support selected at compile time */
#define MPICH_THREAD_LEVEL MPI_THREAD_MULTIPLE

/* Method used to implement refcount updates */
#define MPICH_THREAD_REFCOUNT MPICH_REFCOUNT__LOCKFREE

/* define to disable reference counting predefined objects like MPI_COMM_WORLD
   */
/* #undef MPICH_THREAD_SUPPRESS_PREDEFINED_REFCOUNTS */

/* CH4 should build locality info */
#define MPIDI_BUILD_CH4_LOCALITY_INFO 1

/* CH4 Directly transfers data through the chosen netmode */
/* #undef MPIDI_CH4_DIRECT_NETMOD */

/* Number of VCIs configured in CH4 */
#define MPIDI_CH4_MAX_VCIS 64

/* Define to use bgq capability set */
/* #undef MPIDI_CH4_OFI_USE_SET_BGQ */

/* Define to use cxi capability set */
/* #undef MPIDI_CH4_OFI_USE_SET_CXI */

/* Define to use gni capability set */
/* #undef MPIDI_CH4_OFI_USE_SET_GNI */

/* Define to use PSM2 capability set */
/* #undef MPIDI_CH4_OFI_USE_SET_PSM2 */

/* Define to use PSM3 capability set */
/* #undef MPIDI_CH4_OFI_USE_SET_PSM3 */

/* Define to use runtime capability set */
#define MPIDI_CH4_OFI_USE_SET_RUNTIME 1

/* Define to use sockets capability set */
/* #undef MPIDI_CH4_OFI_USE_SET_SOCKETS */

/* Define to use verbs;ofi_rxm capability set */
/* #undef MPIDI_CH4_OFI_USE_SET_VERBS_RXM */

/* Define if GPU IPC submodule is enabled */
#define MPIDI_CH4_SHM_ENABLE_GPU 1

/* Enable XPMEM shared memory submodule in CH4 */
/* #undef MPIDI_CH4_SHM_ENABLE_XPMEM */

/* Silently disable XPMEM, if it fails at runtime */
/* #undef MPIDI_CH4_SHM_XPMEM_ALLOW_SILENT_FALLBACK */

/* Define to enable direct multi-threading model */
#define MPIDI_CH4_USE_MT_DIRECT 1

/* Define to enable lockless multi-threading model */
/* #undef MPIDI_CH4_USE_MT_LOCKLESS */

/* Define to enable runtime multi-threading model */
/* #undef MPIDI_CH4_USE_MT_RUNTIME */

/* Method used to select vci */
#define MPIDI_CH4_VCI_METHOD MPICH_VCI__COMM

/* Enables AM-only communication */
/* #undef MPIDI_ENABLE_AM_ONLY */

/* CH4/OFI should use domain for vni contexts */
#define MPIDI_OFI_VNI_USE_DOMAIN 1

/* Define to turn on the inlining optimizations in Nemesis code */
/* #undef MPID_NEM_INLINE */

/* Method for local large message transfers. */
/* #undef MPID_NEM_LOCAL_LMT_IMPL */

/* Define if a port may be used to communicate with the processes */
/* #undef MPIEXEC_ALLOW_PORT */

/* limits.h _MAX constant for MPI_Aint */
#define MPIR_AINT_MAX LONG_MAX

/* limits.h _MIN constant for MPI_Aint */
#define MPIR_AINT_MIN LONG_MIN

/* limits.h _MAX constant for MPI_Count */
#define MPIR_COUNT_MAX LLONG_MAX

/* a C type used to compute C++ bool reductions */
#define MPIR_CXX_BOOL_CTYPE _Bool

/* The C type for FORTRAN DOUBLE PRECISION */
#define MPIR_FC_DOUBLE_CTYPE double

/* The C type for FORTRAN REAL */
#define MPIR_FC_REAL_CTYPE float

/* C type to use for MPI_INTEGER16 */
/* #undef MPIR_INTEGER16_CTYPE */

/* C type to use for MPI_INTEGER1 */
#define MPIR_INTEGER1_CTYPE char

/* C type to use for MPI_INTEGER2 */
#define MPIR_INTEGER2_CTYPE short

/* C type to use for MPI_INTEGER4 */
#define MPIR_INTEGER4_CTYPE int

/* C type to use for MPI_INTEGER8 */
#define MPIR_INTEGER8_CTYPE long

/* limits.h _MAX constant for MPI_Offset */
#define MPIR_OFFSET_MAX LLONG_MAX

/* C type to use for MPI_REAL16 */
#define MPIR_REAL16_CTYPE __float128

/* C type to use for MPI_REAL4 */
#define MPIR_REAL4_CTYPE float

/* C type to use for MPI_REAL8 */
#define MPIR_REAL8_CTYPE double

/* MPIR_Ucount is an unsigned MPI_Count-sized integer */
#define MPIR_Ucount unsigned long long

/* Define to enable timing mutexes */
/* #undef MPIU_MUTEX_WAIT_TIME */

/* Define if /bin must be in path */
/* #undef NEEDS_BIN_IN_PATH */

/* Define if environ decl needed */
/* #undef NEEDS_ENVIRON_DECL */

/* Define if gethostname needs a declaration */
/* #undef NEEDS_GETHOSTNAME_DECL */

/* Define if getsid needs a declaration */
/* #undef NEEDS_GETSID_DECL */

/* Define if _POSIX_SOURCE needed to get sigaction */
/* #undef NEEDS_POSIX_FOR_SIGACTION */

/* Define if putenv needs a declaration */
/* #undef NEEDS_PUTENV_DECL */

/* Define if snprintf needs a declaration */
/* #undef NEEDS_SNPRINTF_DECL */

/* Define if strdup needs a declaration */
/* #undef NEEDS_STRDUP_DECL */

/* Define if strerror_r needs a declaration */
/* #undef NEEDS_STRERROR_R_DECL */

/* Define if strict alignment memory access is required */
#define NEEDS_STRICT_ALIGNMENT 1

/* Define if strsignal needs a declaration */
/* #undef NEEDS_STRSIGNAL_DECL */

/* Define if vsnprintf needs a declaration */
/* #undef NEEDS_VSNPRINTF_DECL */

/* The PMI library does not have PMI_Spawn_multiple. */
/* #undef NO_PMI_SPAWN_MULTIPLE */

/* Name of package */
#define PACKAGE "mpich"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "discuss@mpich.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "MPICH"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "MPICH 4.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "mpich"

/* Define to the home page for this package. */
#define PACKAGE_URL "http://www.mpich.org/"

/* Define to the version of this package. */
#define PACKAGE_VERSION "4.1"

/* Define to turn on the prefetching optimization in Nemesis code */
/* #undef PREFETCH_CELL */

/* The size of `bool', as computed by sizeof. */
#define SIZEOF_BOOL 1

/* The size of `char', as computed by sizeof. */
#define SIZEOF_CHAR 1

/* The size of `Complex', as computed by sizeof. */
#define SIZEOF_COMPLEX 8

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE 8

/* The size of `DoubleComplex', as computed by sizeof. */
#define SIZEOF_DOUBLECOMPLEX 16

/* The size of `double_int', as computed by sizeof. */
#define SIZEOF_DOUBLE_INT 16

/* The size of `double _Complex', as computed by sizeof. */
#define SIZEOF_DOUBLE__COMPLEX 16

/* Define size of PAC_TYPE_NAME */
#define SIZEOF_F77_DOUBLE_PRECISION 8

/* Define size of PAC_TYPE_NAME */
#define SIZEOF_F77_INTEGER 4

/* Define size of PAC_TYPE_NAME */
#define SIZEOF_F77_LOGICAL 4

/* Define size of PAC_TYPE_NAME */
#define SIZEOF_F77_REAL 4

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT 4

/* The size of `float_int', as computed by sizeof. */
#define SIZEOF_FLOAT_INT 8

/* The size of `float _Complex', as computed by sizeof. */
#define SIZEOF_FLOAT__COMPLEX 8

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `int16_t', as computed by sizeof. */
#define SIZEOF_INT16_T 2

/* The size of `int32_t', as computed by sizeof. */
#define SIZEOF_INT32_T 4

/* The size of `int64_t', as computed by sizeof. */
#define SIZEOF_INT64_T 8

/* The size of `int8_t', as computed by sizeof. */
#define SIZEOF_INT8_T 1

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `LongDoubleComplex', as computed by sizeof. */
#define SIZEOF_LONGDOUBLECOMPLEX 32

/* The size of `long double', as computed by sizeof. */
#define SIZEOF_LONG_DOUBLE 16

/* The size of `long_double_int', as computed by sizeof. */
#define SIZEOF_LONG_DOUBLE_INT 32

/* The size of `long double _Complex', as computed by sizeof. */
#define SIZEOF_LONG_DOUBLE__COMPLEX 32

/* The size of `long_int', as computed by sizeof. */
#define SIZEOF_LONG_INT 16

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `MPII_Bsend_data_t', as computed by sizeof. */
#define SIZEOF_MPII_BSEND_DATA_T 96

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `short_int', as computed by sizeof. */
#define SIZEOF_SHORT_INT 8

/* The size of `two_int', as computed by sizeof. */
#define SIZEOF_TWO_INT 8

/* The size of `uint16_t', as computed by sizeof. */
#define SIZEOF_UINT16_T 2

/* The size of `uint32_t', as computed by sizeof. */
#define SIZEOF_UINT32_T 4

/* The size of `uint64_t', as computed by sizeof. */
#define SIZEOF_UINT64_T 8

/* The size of `uint8_t', as computed by sizeof. */
#define SIZEOF_UINT8_T 1

/* The size of `unsigned char', as computed by sizeof. */
#define SIZEOF_UNSIGNED_CHAR 1

/* The size of `unsigned int', as computed by sizeof. */
#define SIZEOF_UNSIGNED_INT 4

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 8

/* The size of `unsigned long long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG_LONG 8

/* The size of `unsigned short', as computed by sizeof. */
#define SIZEOF_UNSIGNED_SHORT 2

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P 8

/* The size of `wchar_t', as computed by sizeof. */
#define SIZEOF_WCHAR_T 4

/* The size of `_Bool', as computed by sizeof. */
#define SIZEOF__BOOL 1

/* The size of `_Float16', as computed by sizeof. */
#define SIZEOF__FLOAT16 0

/* The size of `__float128', as computed by sizeof. */
#define SIZEOF___FLOAT128 16

/* Define calling convention */
#define STDCALL 

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if strerror_r returns char *. */
#define STRERROR_R_CHAR_P 1

/* Define TRUE */
#define TRUE 1

/* Define if MPI_Aint should be used instead of void * for storing attribute
   values */
/* #undef USE_AINT_FOR_ATTRVAL */

/* Define if performing coverage tests */
/* #undef USE_COVERAGE */

/* Define to use the fastboxes in Nemesis code */
/* #undef USE_FASTBOX */

/* Define if file should be used for name publisher */
/* #undef USE_FILE_FOR_NAMEPUB */

/* Define if the length of a CHARACTER*(*) string in Fortran should be passed
   as size_t instead of int */
/* #undef USE_FORT_STR_LEN_SIZET */

/* Define to enable memory tracing */
/* #undef USE_MEMORY_TRACING */

/* Define if mpiexec should create a new process group session */
/* #undef USE_NEW_SESSION */

/* Define if using Slurm PMI 1 */
/* #undef USE_PMI1_SLURM */

/* Define if PMI2 API must be used */
/* #undef USE_PMI2_API */

/* Define if using CRAY PMI 2 */
/* #undef USE_PMI2_CRAY */

/* Define if using Slurm PMI 2 */
/* #undef USE_PMI2_SLURM */

/* Define if PMIx API must be used */
/* #undef USE_PMIX_API */

/* Define if sigaction should be used to set signals */
/* #undef USE_SIGACTION */

/* Define if signal should be used to set signals */
/* #undef USE_SIGNAL */

/* Define it the socket verify macros should be enabled */
/* #undef USE_SOCK_VERIFY */

/* Define if we can use a symmetric heap */
#define USE_SYM_HEAP 1

/* Enable extensions on AIX 3, Interix.  */
#ifndef _ALL_SOURCE
# define _ALL_SOURCE 1
#endif
/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif
/* Enable threading extensions on Solaris.  */
#ifndef _POSIX_PTHREAD_SEMANTICS
# define _POSIX_PTHREAD_SEMANTICS 1
#endif
/* Enable extensions on HP NonStop.  */
#ifndef _TANDEM_SOURCE
# define _TANDEM_SOURCE 1
#endif
/* Enable general extensions on Solaris.  */
#ifndef __EXTENSIONS__
# define __EXTENSIONS__ 1
#endif


/* Define if weak symbols should be used */
#define USE_WEAK_SYMBOLS 1

/* Version number of package */
#define VERSION "4.1"

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define if words are little endian */
#define WORDS_LITTLEENDIAN 1

/* Define if configure will not tell us, for universal binaries */
/* #undef WORDS_UNIVERSAL_ENDIAN */

/* Define to 1 if on MINIX. */
/* #undef _MINIX */

/* Define to 2 if the system does not provide POSIX.1 features except with
   this defined. */
/* #undef _POSIX_1_SOURCE */

/* Define to 1 if you need to in order for `stat' and other things to work. */
/* #undef _POSIX_SOURCE */

/* define if bool is a built-in type */
/* #undef bool */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define restrict __restrict
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define if socklen_t is not defined */
/* #undef socklen_t */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */


/* Include nopackage.h to undef autoconf-defined macros that cause conflicts in
 * subpackages.  This should not be necessary, but some packages are too
 * tightly intertwined right now (such as ROMIO and the MPICH core) */
#include "nopackage.h"

#endif /* !defined(MPICHCONF_H_INCLUDED) */

