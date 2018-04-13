############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(bwa-spades C)

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

add_executable(spades-bwa
               bamlite.c
               bwa.c
               bwamem.c
               bwase.c
               bwashm.c
               bwtaln.c
               bwtgap.c
               bwtindex.c
               bwt_lite.c
               bwtsw2_aux.c
               bwtsw2_core.c
               bwtsw2_pair.c
               fastmap.c
               kstring.c
               main.c
               utils.c
               bntseq.c
               bwape.c
               bwaseqio.c
               bwt.c
               bwt_gen.c
               bwtsw2_chain.c
               bwtsw2_main.c
               is.c
               ksw.c
               pemerge.c
               QSufSort.c
               kopen.c
               bwamem_extra.c
               bwamem_pair.c
               kthread.c
               kstring.c
               ksw.c
               malloc_wrap.c
               )

target_link_libraries(spades-bwa z m)
if (HAVE_SHM_OPEN_IN_RT)
  target_link_libraries(spades-bwa rt)
endif()

if (SPADES_STATIC_BUILD)
  set_target_properties(spades-bwa PROPERTIES LINK_SEARCH_START_STATIC 1)
endif()

install(TARGETS spades-bwa
        RUNTIME DESTINATION bin)
