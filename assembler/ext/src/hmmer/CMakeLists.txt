project(hmmer C)

set(HMMER_DATE "July 2016")
set(HMMER_COPYRIGHT "Copyright (C) 2016 Howard Hughes Medical Institute.")
set(HMMER_LICENSE "Freely distributed under the BSD open source license.")
set(HMMER_VERSION "3.1b3")
set(HMMER_URL "http://hmmer.org/")
set(HMMER_SADIR "libdivsufsort")

set(HMMER_MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(HMMER_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(HMMER_TOOLS_BINARY_DIR ${HMMER_BINARY_DIR}/bin)
set(HMMER_BUILT_INCLUDE_DIR ${HMMER_BINARY_DIR}/include)

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
  set(p7_IMPL_SSE 1)
endif()

find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)
  set(HAVE_PTHREAD 1)
  set(HMMER_THREADS 1)
endif()

set(eslLIBRARY 1)

include(CheckFunctionExists)
check_function_exists(fstat HAVE_FSTAT)
check_function_exists(getcwd HAVE_GETCWD)
check_function_exists(getpid HAVE_GETPID)
check_function_exists(mkstemp HAVE_MKSTEMP)
check_function_exists(popen HAVE_POPEN)
check_function_exists(stat HAVE_STAT)
check_function_exists(strcasecmp HAVE_STRCASECMP)
check_function_exists(sysconf HAVE_SYSCONF)
check_function_exists(sysctl HAVE_SYSCTL)
check_function_exists(times HAVE_TIMES)

configure_file("${HMMER_MAIN_SRC_DIR}/p7_config.h.cmake"
               "${HMMER_BUILT_INCLUDE_DIR}/p7_config.h")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(hmmer_src
  src/build.c src/cachedb.c src/emit.c src/errors.c src/evalues.c src/eweight.c
  src/generic_decoding.c src/generic_fwdback.c src/generic_fwdback_chk.c
  src/generic_fwdback_banded.c src/generic_null2.c src/generic_msv.c
  src/generic_optacc.c src/generic_stotrace.c src/generic_viterbi.c
  src/generic_vtrace.c src/h2_io.c src/heatmap.c src/hmmlogo.c src/hmmdmstr.c
  src/hmmdwrkr.c src/hmmdutils.c src/hmmer.c src/logsum.c src/modelconfig.c
  src/modelstats.c src/mpisupport.c src/seqmodel.c src/tracealign.c
  src/p7_alidisplay.c src/p7_bg.c src/p7_builder.c src/p7_domaindef.c
  src/p7_gbands.c src/p7_gmx.c src/p7_gmxb.c src/p7_gmxchk.c
  src/p7_hmm.c src/p7_hmmcache.c src/p7_hmmfile.c src/p7_hmmwindow.c
  src/p7_null3.c src/p7_pipeline.c src/p7_prior.c src/p7_profile.c
  src/p7_spensemble.c src/p7_tophits.c src/p7_trace.c src/p7_scoredata.c
  src/hmmpgmd2msa.c src/fm_alphabet.c src/fm_general.c src/fm_sse.c src/fm_ssv.c)

if (HAVE_SSE2)
  set(hmmer_sse_src
    src/impl_sse/decoding.c src/impl_sse/fwdback.c src/impl_sse/io.c src/impl_sse/ssvfilter.c
    src/impl_sse/msvfilter.c src/impl_sse/null2.c src/impl_sse/optacc.c src/impl_sse/stotrace.c
    src/impl_sse/vitfilter.c src/impl_sse/p7_omx.c src/impl_sse/p7_oprofile.c src/impl_sse/mpi.c)
  set(hmmer_src ${hmmer_src} ${hmmer_sse_src})
endif()

add_library(hmmer STATIC ${hmmer_src})
target_include_directories(hmmer SYSTEM PUBLIC ${HMMER_BUILT_INCLUDE_DIR} ${HMMER_MAIN_SRC_DIR})
target_link_libraries(hmmer easel m)

message("${easel_includes}")
