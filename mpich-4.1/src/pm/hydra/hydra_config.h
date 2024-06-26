/* hydra_config.h.  Generated from hydra_config.h.in by configure.  */
/* hydra_config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if profiling enabled */
#define ENABLE_PROFILING 1

/* Define if warnings are enabled */
/* #undef ENABLE_WARNINGS */

/* Define to 1 if you have the `alarm' function. */
#define HAVE_ALARM 1

/* Define to 1 if you have the `alloca' function. */
/* #undef HAVE_ALLOCA */

/* Define to 1 if you have the <arpa/inet.h> header file. */
#define HAVE_ARPA_INET_H 1

/* Define if external bss is enabled */
#define HAVE_BSS_EXTERNAL 1

/* Define if persist bss is enabled */
#define HAVE_BSS_PERSIST 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define if error checking is enabled */
#define HAVE_ERROR_CHECKING 1

/* Define if environ extern is available */
#define HAVE_EXTERN_ENVIRON 1

/* Define to 1 if you have the `fcntl' function. */
#define HAVE_FCNTL 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define if GNU __attribute__ is supported */
#define HAVE_GCC_ATTRIBUTE 1

/* Define to 1 if you have the `gethostname' function. */
#define HAVE_GETHOSTNAME 1

/* Define to '1' if getifaddrs is present and usable */
#define HAVE_GETIFADDRS 1

/* Define to 1 if you have the `getpgid' function. */
#define HAVE_GETPGID 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the `hstrerror' function. */
#define HAVE_HSTRERROR 1

/* Define if hwloc is available */
#define HAVE_HWLOC 1

/* Define to 1 if you have the <ifaddrs.h> header file. */
#define HAVE_IFADDRS_H 1

/* Define to 1 if you have the `inet_ntop' function. */
#define HAVE_INET_NTOP 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `isatty' function. */
#define HAVE_ISATTY 1

/* Define to 1 if you have the `killpg' function. */
#define HAVE_KILLPG 1

/* Define if C99-style variable argument list macro functionality */
#define HAVE_MACRO_VA_ARGS 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <netdb.h> header file. */
#define HAVE_NETDB_H 1

/* Define to 1 if you have the <netinet/in.h> header file. */
#define HAVE_NETINET_IN_H 1

/* Define to 1 if you have the <netinet/tcp.h> header file. */
#define HAVE_NETINET_TCP_H 1

/* Define to 1 if you have the `poll' function. */
#define HAVE_POLL 1

/* Define to 1 if you have the <poll.h> header file. */
#define HAVE_POLL_H 1

/* Define this if POSIX compliant regcomp()/regexec() */
#define HAVE_POSIX_REGCOMP 1

/* Define to 1 if you have the <sched.h> header file. */
#define HAVE_SCHED_H 1

/* Define to 1 if you have the `select' function. */
#define HAVE_SELECT 1

/* Define to 1 if you have the `setsid' function. */
#define HAVE_SETSID 1

/* Define to 1 if you have the `sigaction' function. */
#define HAVE_SIGACTION 1

/* Define to 1 if you have the `signal' function. */
#define HAVE_SIGNAL 1

/* Define to 1 if you have the `sigset' function. */
#define HAVE_SIGSET 1

/* Define if slurm is available */
/* #undef HAVE_SLURM */

/* Define to 1 if you have the `stat' function. */
#define HAVE_STAT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strsignal' function. */
#define HAVE_STRSIGNAL 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/socket.h> header file. */
#define HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <sys/un.h> header file. */
#define HAVE_SYS_UN_H 1

/* Define to 1 if you have the `time' function. */
#define HAVE_TIME 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define if tm.h and library are usable. */
/* #undef HAVE_TM_H */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `unsetenv' function. */
#define HAVE_UNSETENV 1

/* Define to 1 if you have the `usleep' function. */
#define HAVE_USLEEP 1

/* Define to 1 if you have the <windows.h> header file. */
/* #undef HAVE_WINDOWS_H */

/* Define to 1 if you have the <winsock2.h> header file. */
/* #undef HAVE_WINSOCK2_H */

/* Definition of enabled demux engines */
#define HYDRA_AVAILABLE_DEMUXES "poll select"

/* Definition of enabled launchers */
#define HYDRA_AVAILABLE_LAUNCHERS "ssh rsh fork slurm ll lsf sge manual persist"

/* Definition of enabled RMKS */
#define HYDRA_AVAILABLE_RMKS "user slurm ll lsf sge pbs cobalt"

/* Definition of enabled processor topology libraries */
#define HYDRA_AVAILABLE_TOPOLIBS "hwloc"

/* C compiler */
#define HYDRA_CC "gcc      "

/* Configure arguments */
#define HYDRA_CONFIGURE_ARGS_CLEAN "'--disable-option-checking' '--prefix=/usr/local/mpich-4.1' '--cache-file=/dev/null' '--srcdir=.' 'CC=gcc' 'CFLAGS= -O2' 'LDFLAGS=' 'LIBS=' 'CPPFLAGS= -DNETMOD_INLINE=__netmod_inline_ofi__ -I/home/superdouding/apollo/my_em_planner/mpich-4.1/src/mpl/include -I/home/superdouding/apollo/my_em_planner/mpich-4.1/modules/json-c -D_REENTRANT -I/home/superdouding/apollo/my_em_planner/mpich-4.1/src/mpi/romio/include -I/home/superdouding/apollo/my_em_planner/mpich-4.1/src/pmi/include -I/home/superdouding/apollo/my_em_planner/mpich-4.1/modules/yaksa/src/frontend/include -I/home/superdouding/apollo/my_em_planner/mpich-4.1/modules/libfabric/include'"

/* Default demux engine */
#define HYDRA_DEFAULT_DEMUX "poll"

/* Default launcher */
#define HYDRA_DEFAULT_LAUNCHER "ssh"

/* Default RMK */
#define HYDRA_DEFAULT_RMK "user"

/* Default processor topology library */
#define HYDRA_DEFAULT_TOPOLIB "hwloc"

/* Hydra PMI proxy executable */
#define HYDRA_PMI_PROXY "hydra_pmi_proxy"

/* Hydra release date information */
#define HYDRA_RELEASE_DATE "Fri Jan 27 13:54:44 CST 2023"

/* Hydra version information */
#define HYDRA_VERSION "4.1"

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Define if environ decl needed */
/* #undef MANUAL_EXTERN_ENVIRON */

/* Define if gethostname needs a declaration */
/* #undef NEEDS_GETHOSTNAME_DECL */

/* Define if getpgid needs a declaration */
/* #undef NEEDS_GETPGID_DECL */

/* Define if gettimeofday needs a declaration */
/* #undef NEEDS_GETTIMEOFDAY_DECL */

/* Define if hstrerror needs a declaration */
/* #undef NEEDS_HSTRERROR_DECL */

/* Define if killpg needs a declaration */
/* #undef NEEDS_KILLPG_DECL */

/* Define if _POSIX_SOURCE needed to get sigaction */
/* #undef NEEDS_POSIX_FOR_SIGACTION */

/* Define if strsignal needs a declaration */
/* #undef NEEDS_STRSIGNAL_DECL */

/* Name of package */
#define PACKAGE "hydra"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME "Hydra"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Hydra 4.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "hydra"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "4.1"

/* Define if we should check for PMI key collisions */
/* #undef PMI_KEY_CHECK */

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 8

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define if performing coverage tests */
/* #undef USE_COVERAGE */

/* Define if memory tracing is enabled */
/* #undef USE_MEMORY_TRACING */

/* Define if sigaction should be used to set signals */
#define USE_SIGACTION /**/

/* Define if signal should be used to set signals */
/* #undef USE_SIGNAL */

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


/* Version number of package */
#define VERSION "4.1"

/* Define to 1 if on MINIX. */
/* #undef _MINIX */

/* Define to 2 if the system does not provide POSIX.1 features except with
   this defined. */
/* #undef _POSIX_1_SOURCE */

/* Define to 1 if you need to in order for `stat' and other things to work. */
/* #undef _POSIX_SOURCE */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

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

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */
