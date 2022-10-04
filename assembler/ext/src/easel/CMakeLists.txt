project(easel C)

set(EASEL_MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(EASEL_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(EASEL_TOOLS_BINARY_DIR ${EASEL_BINARY_DIR}/bin)
set(EASEL_BUILT_INCLUDE_DIR ${EASEL_BINARY_DIR}/include)

set(EASEL_DATE "Nov 2020")
set(EASEL_COPYRIGHT "Copyright (C) 2020 Howard Hughes Medical Institute")
set(EASEL_LICENSE "Freely distributed under the BSD open source license.")
set(EASEL_VERSION 0.48)
set(EASEL_URL "http://bioeasel.org/")

include(CheckIncludeFiles)
check_include_files(endian.h HAVE_ENDIAN_H)
check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(sys/param.h HAVE_SYS_PARAM_H)
check_include_files(sys/sysctl.h HAVE_SYS_SYSCTL_H)
check_include_files(netinet/in.h HAVE_NETINET_IN_H)

include(CheckCSourceCompiles)
check_c_source_compiles("
#include <x86intrin.h>
#include <stdint.h>
int stub_sse(void) {
__m128i v1 = _mm_set1_epi32(42);
__m128i v2 = _mm_set1_epi32(86);
union { __m128i v; int32_t x[4]; } v3;
v3.v = _mm_add_epi32(v1, v2);
return (int) v3.x[0];
}" eslENABLE_SSE)

check_c_source_compiles("
#include <x86intrin.h>
#include <stdint.h>
int stub_sse4(void) {
__m128i v1 = _mm_set1_epi8(-42);
__m128i v2 = _mm_set1_epi8(-86);
union { __m128i v; int8_t x[16]; } v3;
v3.v = _mm_adds_epi8(v1, v2);
v2   = _mm_max_epi8(v1, v1);
return (int) -v3.x[0];
}
int main(void) { if (stub_sse4() != 128) return 1; else return 0; }" eslENABLE_SSE4)

check_c_source_compiles("
#include <x86intrin.h>
#include <stdio.h>
#include <stdint.h>
int main(void) {
  __m256i v1 = _mm256_set1_epi32(42);
  __m256i v2 = _mm256_set1_epi32(214);
  union { __m256i v; int32_t x[8]; } v3;
  v3.v = _mm256_add_epi32(v1, v2);
  printf(\"%d\n\", v3.x[2]);}" eslENABLE_AVX)

check_c_source_compiles("
#include <x86intrin.h>
#include <stdint.h>
int main(void) {
  __m512i v1 = _mm512_set1_epi8(21);
  union { __m256i v; int8_t x[32]; } v2;
  v2.v = _mm512_extracti32x8_epi32(_mm512_adds_epi8(v1, v1), 0x1);
  return (v2.x[0] == 42 ? 0 : 1);}" eslENABLE_AVX512)

check_c_source_compiles("
#include <stdint.h>
#include <arm_neon.h>
int stub_neon(void) {
  int32x4_t v1 = vdupq_n_s32(42);
  int32x4_t v2 = vdupq_n_s32(86);
  union { int32x4_t v; int32_t x[4]; } v3;
  v3.v = vaddq_s32(v1, v2);
  return (int) v3.x[0];
}
int main(void) { if (stub_neon() != 128) return 1; else return 0; }" eslENABLE_NEON)

check_c_source_compiles("
#include <stdint.h>
#include <arm_neon.h>
int main(void) {
  int16x8_t   a1 = { 1, 2, 8, 3, 4, 5, 6, 7 };
  int16_t     r1 = vmaxvq_s16(a1); // = 8 : horizontal max
  uint8x16_t  a2 = { 1, 2, 16, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
  uint8_t     r2 = vmaxvq_u8(a2);  // = 16 : horizontal max
  float32x4_t a3 = { 1.0, 2.0, 3.0, 4.0 };
  float       r3 = vaddvq_f32(a3); // = 10.0 : horizontal sum
  return 0;
}" eslHAVE_NEON_AARCH64)

check_c_source_compiles("
#include <pmmintrin.h>
int main(void) {
  _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON);
  return 0;
}" HAVE_FLUSH_ZERO_MODE)
  
check_c_source_compiles("
#include <xmmintrin.h>
int main(void) {
  _MM_SET_DENORMALS_ZERO_MODE (_MM_DENORMALS_ZERO_ON);
  return 0;
}" HAVE_DENORMALS_ZERO_MODE)

include(CheckFunctionExists)
check_function_exists(aligned_alloc HAVE_ALIGNED_ALLOC)
check_function_exists(erfc HAVE_ERFC)
check_function_exists(getpid HAVE_GETPID)
check_function_exists(_mm_malloc HAVE__MM_MALLOC)
check_function_exists(popen HAVE_POPEN)
check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
check_function_exists(strcasecmp HAVE_STRCASECMP)
check_function_exists(strsep HAVE_STRSEP)
check_function_exists(sysconf HAVE_SYSCONF)
check_function_exists(sysctl HAVE_SYSCTL)
check_function_exists(times HAVE_TIMES)
check_function_exists(fseeko HAVE_FSEEKO)

find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)
  set(HAVE_PTHREAD 1)
endif()

set(eslLIBRARY 1)
set(HAVE_FUNC_ATTRIBUTE_NORETURN 1)
set(HAVE_FUNC_ATTRIBUTE_FORMAT 1)

configure_file("${EASEL_MAIN_SRC_DIR}/esl_config.h.cmake"
               "${EASEL_BUILT_INCLUDE_DIR}/esl_config.h")

set(easel_src
  easel.c
  esl_alloc.c esl_alphabet.c esl_arr2.c esl_arr3.c esl_bitfield.c
  esl_buffer.c esl_cluster.c esl_composition.c esl_cpu.c esl_dirichlet.c
  esl_distance.c esl_dmatrix.c esl_dsqdata.c esl_exponential.c esl_fileparser.c
  esl_gamma.c esl_gencode.c esl_getopts.c esl_gev.c esl_graph.c esl_gumbel.c
  esl_heap.c esl_histogram.c esl_hmm.c esl_huffman.c esl_hyperexp.c
  esl_iset.c esl_json.c esl_keyhash.c esl_lognormal.c esl_matrixops.c
  esl_mem.c esl_minimizer.c esl_mixdchlet.c esl_mixgev.c esl_mpi.c
  esl_msa.c esl_msacluster.c esl_msafile.c esl_msafile2.c esl_msafile_a2m.c
  esl_msafile_afa.c esl_msafile_clustal.c esl_msafile_phylip.c esl_msafile_psiblast.c
  esl_msafile_selex.c esl_msafile_stockholm.c esl_msashuffle.c esl_msaweight.c
  esl_normal.c esl_paml.c esl_quicksort.c esl_random.c esl_rand64.c esl_randomseq.c
  esl_ratematrix.c esl_recorder.c esl_red_black.c esl_regexp.c esl_rootfinder.c
  esl_scorematrix.c esl_sq.c esl_sqio.c esl_sqio_ascii.c esl_sqio_ncbi.c
  esl_ssi.c esl_stack.c esl_stats.c esl_stopwatch.c esl_stretchexp.c
  esl_subcmd.c esl_threads.c esl_tree.c esl_varint.c esl_vectorops.c esl_weibull.c
  esl_workqueue.c esl_wuss.c
  esl_sse.c esl_avx.c esl_avx512.c esl_neon.c
)
    
add_library(easel STATIC ${easel_src})
target_link_libraries(easel m)
target_include_directories(easel SYSTEM PUBLIC ${EASEL_BUILT_INCLUDE_DIR} ${EASEL_MAIN_SRC_DIR})
