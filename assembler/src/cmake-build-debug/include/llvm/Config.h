#ifndef __LLVM_CONFIG_H__
#define __LLVM_CONFIG_H__

#define LLVM_ON_UNIX 1
#define LLVM_ENABLE_CRASH_DUMPS 0
#define LLVM_ENABLE_THREADS 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <execinfo.h> header file. */
#define HAVE_EXECINFO_H 1

/* Define to 1 if you have the <link.h> header file. */
#define HAVE_LINK_H 1

/* Define to 1 if you have the <mach/mach.h> header file. */
/* #undef HAVE_MACH_MACH_H */

/* Define to 1 if you have the <signal.h> header file. */
#define HAVE_SIGNAL_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/resouce.h> header file. */
#define HAVE_SYS_RESOURCE_H 1

/* Define to 1 if you have the `backtrace' function. */
#define HAVE_BACKTRACE TRUE

#define BACKTRACE_HEADER <execinfo.h>

/* Define to 1 if you have the <cxxabi.h> header file. */
/* #undef HAVE_CXXABI_H */

/* Define as the return type of signal handlers (`int' or `void'). */
#define RETSIGTYPE void

/* Define to 1 if you have the `futimens' function. */
#define HAVE_FUTIMENS 1

/* Define to 1 if you have the `futimes' function. */
#define HAVE_FUTIMES 1

/* Define to 1 if you have the <sys/mman.h> header file. */
#define HAVE_SYS_MMAN_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the `getrusage' function. */
#define HAVE_GETRUSAGE 1

/* Define to 1 if you have the `sysconf' function. */
#define HAVE_SYSCONF 1

/* Define to 1 if you have the `mallctl' function. */
/* #undef HAVE_MALLCTL */

/* Define to 1 if you have the `mallinfo' function. */
#define HAVE_MALLINFO 1

/* Define to 1 if you have the <malloc/malloc.h> header file. */
/* #undef HAVE_MALLOC_MALLOC_H */

/* Define to 1 if you have the `malloc_zone_statistics' function. */
/* #undef HAVE_MALLOC_ZONE_STATISTICS */

/* Define to 1 if you have the `pthread' library (-lpthread). */
#define HAVE_LIBPTHREAD 1

/* Define to 1 if you have the `pthread_getname_np' function. */
/* #undef HAVE_PTHREAD_GETNAME_NP */

/* Define to 1 if you have the `pthread_setname_np' function. */
/* #undef HAVE_PTHREAD_SETNAME_NP */

/* Have pthread_getspecific */
#define HAVE_PTHREAD_GETSPECIFIC 1

/* Define to 1 if you have the <pthread.h> header file. */
#define HAVE_PTHREAD_H 1

/* Have pthread_mutex_lock */
#define HAVE_PTHREAD_MUTEX_LOCK 1

/* Have pthread_rwlock_init */
#define HAVE_PTHREAD_RWLOCK_INIT 1

#endif
