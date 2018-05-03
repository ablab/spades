project(easel C)

set(EASEL_MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(EASEL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(EASEL_TOOLS_BINARY_DIR ${EASEL_BINARY_DIR}/bin)
set(EASEL_BUILT_INCLUDE_DIR ${EASEL_BINARY_DIR}/include)

set(EASEL_DATE "July 2016")
set(EASEL_COPYRIGHT "Copyright (C) 2016 Howard Hughes Medical Institute")
set(EASEL_LICENSE "Freely distributed under the BSD open source license.")
set(EASEL_VERSION 0.43)
set(EASEL_URL "http://bioeasel.org/")

include(CheckIncludeFiles)

check_include_files(endian.h HAVE_FEATURES_H)
check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(sys/param.h HAVE_SYS_PARAM_H)
check_include_files(sys/sysctl.h HAVE_SYS_SYSCTL_H)

check_include_files(emmintrin.h HAVE_SSE2)
check_include_files(pmmintrin.h HAVE_SSE3)
check_include_files(xmmintrin.h HAVE_SSE)
if (HAVE_SSE2)
  set(HAVE_SSE2_CAST 1)
endif()

include(CheckFunctionExists)
check_function_exists(chmod HAVE_CHMOD)
check_function_exists(fseeko HAVE_FSEEKO)
check_function_exists(fstat HAVE_FSTAT)
check_function_exists(getcwd HAVE_GETCWD)
check_function_exists(getpid HAVE_GETPID)
check_function_exists(mkstemp HAVE_MKSTEMP)
check_function_exists(popen HAVE_POPEN)
check_function_exists(putenv HAVE_PUTENV)
check_function_exists(stat HAVE_STAT)
check_function_exists(strcasecmp HAVE_STRCASECMP)
check_function_exists(sysconf HAVE_SYSCONF)
check_function_exists(sysctl HAVE_SYSCTL)
check_function_exists(times HAVE_TIMES)
check_function_exists(erfc HAVE_ERFC)

find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)
  set(HAVE_PTHREAD 1)
endif()

set(eslLIBRARY 1)
set(HAVE_FUNC_ATTRIBUTE_NORETURN 1)

configure_file("${EASEL_MAIN_SRC_DIR}/esl_config.h.cmake"
               "${EASEL_BUILT_INCLUDE_DIR}/esl_config.h")

set(easel_src
  easel.c
  esl_alphabet.c esl_avx.c esl_buffer.c esl_cluster.c esl_composition.c esl_dirichlet.c
  esl_distance.c esl_dmatrix.c esl_dsqdata.c esl_exponential.c esl_fileparser.c esl_gamma.c
  esl_gencode.c esl_getopts.c esl_gev.c esl_gumbel.c esl_heap.c esl_histogram.c
  esl_hmm.c esl_hyperexp.c esl_keyhash.c esl_mem.c esl_minimizer.c esl_mixgev.c
  esl_mpi.c esl_msa.c esl_msacluster.c esl_msafile.c esl_msafile2.c esl_msafile_a2m.c
  esl_msafile_afa.c esl_msafile_clustal.c esl_msafile_phylip.c esl_msafile_psiblast.c
  esl_msafile_selex.c esl_msafile_stockholm.c esl_msashuffle.c esl_msaweight.c
  esl_normal.c esl_paml.c esl_random.c esl_randomseq.c
  esl_ratematrix.c esl_recorder.c esl_regexp.c esl_rootfinder.c
  esl_scorematrix.c esl_sq.c esl_sqio.c esl_sqio_ascii.c
  esl_sqio_ncbi.c esl_sse.c esl_ssi.c esl_stack.c esl_stats.c esl_stopwatch.c
  esl_stretchexp.c esl_threads.c esl_tree.c esl_vectorops.c esl_vmx.c esl_weibull.c esl_workqueue.c
  esl_wuss.c)

add_library(easel STATIC ${easel_src})
target_link_libraries(easel m)
target_include_directories(easel SYSTEM PUBLIC ${EASEL_BUILT_INCLUDE_DIR} ${EASEL_MAIN_SRC_DIR})
