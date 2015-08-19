# -*- cmake -*-

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(jemalloc C CXX)

include(CheckIncludeFiles)
include(CheckSymbolExists)
include(CheckFunctionExists)
include(CheckCSourceCompiles)
include(CheckTypeSize)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include)
include_directories(${SPADES_BUILT_INCLUDE_DIR})

check_include_files(features.h HAVE_FEATURES_H)
check_include_files(sys/cdefs.h HAVE_SYS_CDEFS_H)

set(JEMALLOC_VERSION "3.2.0-0-g87499f6748ebe4817571e817e9f680ccb5bf54a9")
set(JEMALLOC_VERSION_MAJOR "3")
set(JEMALLOC_VERSION_MINOR "2")
set(JEMALLOC_VERSION_BUGFIX "0")
set(JEMALLOC_VERSION_NREV "0")
set(JEMALLOC_VERSION_GID "87499f6748ebe4817571e817e9f680ccb5bf54a9")

set(TARGET_ARCHITECTURE "generic")
set(_vendor_id)
set(_cpu_family)
set(_cpu_model)
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  file(READ "/proc/cpuinfo" _cpuinfo)
  string(REGEX REPLACE ".*vendor_id[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _vendor_id "${_cpuinfo}")
  string(REGEX REPLACE ".*cpu family[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_family "${_cpuinfo}")
  string(REGEX REPLACE ".*model[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_model "${_cpuinfo}")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  exec_program("/usr/sbin/sysctl -n machdep.cpu.vendor" OUTPUT_VARIABLE _vendor_id)
  exec_program("/usr/sbin/sysctl -n machdep.cpu.model"  OUTPUT_VARIABLE _cpu_model)
  exec_program("/usr/sbin/sysctl -n machdep.cpu.family" OUTPUT_VARIABLE _cpu_family)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  get_filename_component(_vendor_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;VendorIdentifier]" NAME CACHE)
  get_filename_component(_cpu_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;Identifier]" NAME CACHE)
  mark_as_advanced(_vendor_id _cpu_id)
  string(REGEX REPLACE ".* Family ([0-9]+) .*" "\\1" _cpu_family "${_cpu_id}")
  string(REGEX REPLACE ".* Model ([0-9]+) .*" "\\1" _cpu_model "${_cpu_id}")
endif()

if(_vendor_id STREQUAL "GenuineIntel" OR
   _vendor_id STREQUAL "AuthenticAMD")
 set(JEMALLOC_CPU_SPINWAIT "__asm__ volatile(\"pause\")")
endif()

check_symbol_exists(atomic_fetchadd_long "sys/types.h;machine/atomic.h;inttypes.h" JEMALLOC_ATOMIC9)
check_symbol_exists(OSAtomicAdd64 "libkern/OSAtomic.h;inttypes.h" JEMALLOC_OSATOMIC)
check_symbol_exists(OSSpinLockLock "libkern/OSAtomic.h;inttypes.h" JEMALLOC_OSSPIN)
check_function_exists(_malloc_thread_cleanup JEMALLOC_MALLOC_THREAD_CLEANUP)

check_function_exists(_pthread_mutex_init_calloc_cb, JEMALLOC_MUTEX_INIT_CB)
check_function_exists(sbrk JEMALLOC_HAVE_SBRK)

set(CMAKE_REQUIRED_FLAGS "-Werror")
check_c_source_compiles("static __thread int
__attribute__((tls_model(\"initial-exec\"))) foo = 0; int main(void) { return 0; }" JEMALLOC_TLS_MODEL_SUPPORTED)
set(CMAKE_REQUIRED_FLAGS)

if (JEMALLOC_TLS_MODEL_SUPPORTED)
  set(JEMALLOC_TLS_MODEL "__attribute__((tls_model(\"initial-exec\")))")
endif()

check_type_size(void* JEMALLOC_SIZEOF_PTR)
if (JEMALLOC_SIZEOF_PTR EQUAL 8)
  set(JEMALLOC_LG_SIZEOF_PTR 3)
else()
  message(FATAL_ERROR "Unsupported pointer size")
endif()

set(JEMALLOC_USABLE_SIZE_CONST "const")
set(JEMALLOC_TCACHE ON)
set(JEMALLOC_MUNMAP ON)
set(JEMALLOC_IVSALLOC OFF)
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(JEMALLOC_PURGE_MADVISE_DONTNEED ON)
  set(JEMALLOC_THREADED_INIT ON)
  set(JEMALLOC_MUNMAP OFF)
  set(JEMALLOC_USABLE_SIZE_CONST "")
  add_definitions(-D_GNU_SOURCE)
elseif(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
  set(JEMALLOC_LAZY_LOCK ON)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(JEMALLOC_PURGE_MADVISE_FREE ON)
  set(JEMALLOC_ZONE ON)

  set(CMAKE_EXTRA_INCLUDE_FILES "malloc/malloc.h")
  check_type_size(malloc_zone_t JEMALLOC_MZONE_SIZE)
  check_type_size(malloc_introspection_t JEMALLOC_MINT_SIZE)
  set(CMAKE_EXTRA_INCLUDE_FILES)

  if (JEMALLOC_MZONE_SIZE EQUAL 112)
    set(JEMALLOC_ZONE_VERSION "3")
  elseif(JEMALLOC_MZONE_SIZE EQUAL 120)
    set(JEMALLOC_ZONE_VERSION "5")
  elseif(JEMALLOC_MZONE_SIZE EQUAL 128)
    if (JEMALLOC_MINT_SIZE EQUAL 72)
      set(JEMALLOC_ZONE_VERSION "6")
    elseif (JEMALLOC_MINT_SIZE EQUAL 112)
      set(JEMALLOC_ZONE_VERSION "7")
    else()
      message(FATAL_ERROR "Unsupported malloc zone version")
    endif()
  elseif(JEMALLOC_MZONE_SIZE EQUAL 136)
    set(JEMALLOC_ZONE_VERSION "8")
  elseif(JEMALLOC_MZONE_SIZE GREATER 136)
    set(JEMALLOC_ZONE_VERSION "8")
  else()
    message(FATAL_ERROR "Unsupported malloc zone version")
  endif()
  set(JEMALLOC_IVSALLOC ON)
endif()

add_definitions(-D_REENTRANT)
add_compile_options(-O3 -funroll-loops)

# TLS on Darwin is "fake". We cannot use it with jemalloc, because
# emulated tls wants to use calloc() and friends
if (APPLE)
else()
  set(CMAKE_REQUIRED_FLAGS "-Werror")
  check_c_source_compiles("__thread int x;
  int main(void) { x = 42; return 0; }" JEMALLOC_TLS)
  set(CMAKE_REQUIRED_FLAGS)
endif()
check_function_exists(memalign JEMALLOC_OVERRIDE_MEMALIGN)
check_function_exists(valloc JEMALLOC_OVERRIDE_VALLOC)

check_type_size(int JEMALLOC_SIZEOF_INT)
if (JEMALLOC_SIZEOF_INT EQUAL 8)
  set(JEMALLOC_LG_SIZEOF_INT 3)
elseif (JEMALLOC_SIZEOF_INT EQUAL 4)
  set(JEMALLOC_LG_SIZEOF_INT 2)
else()
  message(FATAL_ERROR "Unsupported integer size")
endif()

check_type_size(long JEMALLOC_SIZEOF_LONG)
if (JEMALLOC_SIZEOF_LONG EQUAL 8)
  set(JEMALLOC_LG_SIZEOF_LONG 3)
elseif (JEMALLOC_SIZEOF_LONG EQUAL 4)
  set(JEMALLOC_LG_SIZEOF_LONG 2)
else()
  message(FATAL_ERROR "Unsupported long size")
endif()

check_type_size(intmax_t JEMALLOC_SIZEOF_INTMAX_T)
if (JEMALLOC_SIZEOF_INTMAX_T EQUAL 16)
  set(JEMALLOC_LG_SIZEOF_INTMAX_T 4)
elseif (JEMALLOC_SIZEOF_INTMAX_T EQUAL 8)
  set(JEMALLOC_LG_SIZEOF_INTMAX_T 3)
elseif (JEMALLOC_SIZEOF_INTMAX_T EQUAL 4)
  set(JEMALLOC_LG_SIZEOF_INTMAX_T 2)
else()
  message(FATAL_ERROR "Unsupported intmax_t size")
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CheckPageSize.c.in"
               "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/CheckPageSize.c"
                IMMEDIATE @ONLY)
try_run(JEMALLOC_CHECK_PAGE_SIZE JEMALLOC_HAVE_PAGE_SIZE
        ${CMAKE_BINARY_DIR}
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/CheckPageSize.c"
        RUN_OUTPUT_VARIABLE JEMALLOC_STATIC_PAGE_SHIFT)

configure_file(${EXT_DIR}/include/jemalloc/jemalloc.h.in
               ${SPADES_BUILT_INCLUDE_DIR}/jemalloc/jemalloc.h)

configure_file(${EXT_DIR}/include/jemalloc/jemalloc_defs.h.in
               ${SPADES_BUILT_INCLUDE_DIR}/jemalloc/jemalloc_defs.h)


set(JEMALLOC_SRC arena.c
                 base.c
                 chunk.c
                 chunk_mmap.c
                 ctl.c
                 hash.c
                 jemalloc.c
                 mutex.c
                 quarantine.c
                 stats.c
                 tsd.c
                 atomic.c
                 bitmap.c
                 chunk_dss.c
                 ckh.c
                 extent.c
                 huge.c
                 mb.c
                 prof.c
                 rtree.c
                 tcache.c
                 util.c)


if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(JEMALLOC_SRC ${JEMALLOC_SRC} zone.c)
endif()

add_library(jemalloc STATIC
            ${JEMALLOC_SRC})
