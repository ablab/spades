project(llvm-support)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include)

include(CheckIncludeFile)
include(CheckIncludeFiles)
include(CheckIncludeFileCXX)
include(CheckSymbolExists)

check_include_file(execinfo.h  HAVE_EXECINFO_H)
check_include_file(signal.h    HAVE_SIGNAL_H)
check_include_file(sys/stat.h  HAVE_SYS_STAT_H)
check_include_file(cxxabi.h    HAVE_CXXABI_H)
check_include_file(dlfcn.h     HAVE_DLFCN_H)
check_include_file(mach/mach.h HAVE_MACH_MACH_H)
check_include_file(link.h      HAVE_LINK_H)
check_include_file(unistd.h    HAVE_UNISTD_H)
check_include_file(sys/mman.h  HAVE_SYS_MMAN_H)
check_include_file(fcntl.h     HAVE_FCNTL_H)
check_include_file(malloc/malloc.h HAVE_MALLOC_MALLOC_H)
check_include_file(pthread.h   HAVE_PTHREAD_H)
check_include_file(sys/time.h  HAVE_SYS_TIME_H)
check_include_file(sys/resource.h HAVE_SYS_RESOURCE_H)

find_package(Backtrace)
set(HAVE_BACKTRACE ${Backtrace_FOUND})
set(BACKTRACE_HEADER ${Backtrace_HEADER})

check_symbol_exists(futimens sys/stat.h HAVE_FUTIMENS)
check_symbol_exists(futimes sys/time.h HAVE_FUTIMES)
check_symbol_exists(getrusage sys/resource.h HAVE_GETRUSAGE)
check_symbol_exists(getpagesize unistd.h HAVE_GETPAGESIZE)
check_symbol_exists(sysconf unistd.h HAVE_SYSCONF)
check_symbol_exists(mallinfo malloc.h HAVE_MALLINFO)
check_symbol_exists(mallctl malloc_np.h HAVE_MALLCTL)
check_symbol_exists(malloc_zone_statistics malloc/malloc.h
                    HAVE_MALLOC_ZONE_STATISTICS)

find_library(DL_LIB NAMES "dl")

check_library_exists(pthread pthread_create "" HAVE_LIBPTHREAD)
if (HAVE_LIBPTHREAD)
  check_library_exists(pthread pthread_getspecific "" HAVE_PTHREAD_GETSPECIFIC)
  check_library_exists(pthread pthread_rwlock_init "" HAVE_PTHREAD_RWLOCK_INIT)
  check_library_exists(pthread pthread_mutex_lock "" HAVE_PTHREAD_MUTEX_LOCK)
else()
  check_library_exists(c pthread_create "" PTHREAD_IN_LIBC)
  if (PTHREAD_IN_LIBC)
    check_library_exists(c pthread_getspecific "" HAVE_PTHREAD_GETSPECIFIC)
    check_library_exists(c pthread_rwlock_init "" HAVE_PTHREAD_RWLOCK_INIT)
    check_library_exists(c pthread_mutex_lock "" HAVE_PTHREAD_MUTEX_LOCK)
  endif()
endif()
check_library_exists(dl dlopen "" HAVE_LIBDL)
check_library_exists(rt clock_gettime "" HAVE_LIBRT)

if(HAVE_LIBPTHREAD)
  # We want to find pthreads library and at the moment we do want to
  # have it reported as '-l<lib>' instead of '-pthread'.
  # TODO: switch to -pthread once the rest of the build system can deal with it.
  set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
  set(THREADS_HAVE_PTHREAD_ARG Off)
  find_package(Threads REQUIRED)
  set(LLVM_PTHREAD_LIB ${CMAKE_THREAD_LIBS_INIT})
endif()

add_definitions(-D__STDC_CONSTANT_MACROS)
add_definitions(-D__STDC_LIMIT_MACROS)

# FIXME: Signal handler return type, currently hardcoded to 'void'
set(RETSIGTYPE void)

configure_file(${EXT_DIR}/include/llvm/Config.h.in
               ${SPADES_BUILT_INCLUDE_DIR}/llvm/Config.h)

file(GLOB sources "[a-zA-Z]*.cpp" "[a-zA-Z]*.c")

add_library(llvm-support STATIC
  ${sources})
if (DL_LIB)
  target_link_libraries(llvm-support ${DL_LIB})
endif()
