############################################################################
# Copyright (c) 2016 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(bwa C)

find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)
  add_definitions(-DHAVE_PTHREAD)
endif()

include(CheckFunctionExists)
include(CheckLibraryExists)
check_function_exists(shm_open HAVE_SHM_OPEN)
if (NOT HAVE_SHM_OPEN)
  check_library_exists(rt shm_open "" HAVE_SHM_OPEN_IN_RT)
endif()
if (NOT HAVE_SHM_OPEN AND NOT HAVE_SHM_OPEN_IN_RT)
  error("Cannot find shm_open()")
endif()

add_library(bwa STATIC
  utils.c kthread.c ksw.c bwt.c bntseq.c bwa.c bwamem.c bwamem_pair.c
  bwamem_extra.c malloc_wrap.c is.c bwtindex.c rope.c rle.c kstring.c)

target_link_libraries(bwa z m)
if (HAVE_SHM_OPEN_IN_RT)
  target_link_libraries(bwa rt)
endif()
