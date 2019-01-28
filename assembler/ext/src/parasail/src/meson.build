################
# core sources #
################

parasail_c_core_sources = files([
  'cigar.c',
  'cpuid.c',
  'function_lookup.c',
  'io.c',
  'isastubs.c',
  'matrix_lookup.c',
  'memory.c',
  'parser.c',
  'pssw.c',
  'time.c',
  'nw_dispatch.c',
  'sg_dispatch.c',
  'sw_dispatch.c',
  'dispatch_profile.c',
  'satcheck.c',
  'striped_unwind.c',
  'traceback.c'])

#################
# vector memory #
#################

parasail_c_mem_sse_sources = files(['memory_sse.c'])
parasail_c_mem_avx2_sources = files(['memory_avx2.c'])

##################
# serial methods #
##################

parasail_c_novec_sources = files([
  'nw.c',
  'sg.c',
  'sw.c',

  'nw_scan.c',
  'sg_scan.c',
  'sw_scan.c',

  'nw_banded.c',

  'nw_stats.c',
  'sg_stats.c',
  'sw_stats.c',

  'nw_stats_scan.c',
  'sg_stats_scan.c',
  'sw_stats_scan.c'])

parasail_c_trace_novec_sources = files([
  'nw_trace.c',
  'sg_trace.c',
  'sw_trace.c',

  'nw_trace_scan.c',
  'sg_trace_scan.c',
  'sw_trace_scan.c'])

####################
# parallel methods #
####################

parasail_c_sse2_sources = []
parasail_c_sse41_sources = []
parasail_c_avx2_sources = []

# parallel scan methods
parasail_c_sse2_sources += files([
  'nw_scan_sse2_128_64.c',
  'sg_scan_sse2_128_64.c',
  'sw_scan_sse2_128_64.c',

  'nw_scan_sse2_128_32.c',
  'sg_scan_sse2_128_32.c',
  'sw_scan_sse2_128_32.c',

  'nw_scan_sse2_128_16.c',
  'sg_scan_sse2_128_16.c',
  'sw_scan_sse2_128_16.c',

  'nw_scan_sse2_128_8.c',
  'sg_scan_sse2_128_8.c',
  'sw_scan_sse2_128_8.c'])

parasail_c_sse41_sources += files([
  'nw_scan_sse41_128_64.c',
  'sg_scan_sse41_128_64.c',
  'sw_scan_sse41_128_64.c',

  'nw_scan_sse41_128_32.c',
  'sg_scan_sse41_128_32.c',
  'sw_scan_sse41_128_32.c',

  'nw_scan_sse41_128_16.c',
  'sg_scan_sse41_128_16.c',
  'sw_scan_sse41_128_16.c',

  'nw_scan_sse41_128_8.c',
  'sg_scan_sse41_128_8.c',
  'sw_scan_sse41_128_8.c'])

parasail_c_avx2_sources += files([
  'nw_scan_avx2_256_64.c',
  'sg_scan_avx2_256_64.c',
  'sw_scan_avx2_256_64.c',

  'nw_scan_avx2_256_32.c',
  'sg_scan_avx2_256_32.c',
  'sw_scan_avx2_256_32.c',

  'nw_scan_avx2_256_16.c',
  'sg_scan_avx2_256_16.c',
  'sw_scan_avx2_256_16.c',

  'nw_scan_avx2_256_8.c',
  'sg_scan_avx2_256_8.c',
  'sw_scan_avx2_256_8.c'])

# parallel diag methods
parasail_c_sse2_sources += files([
  'nw_diag_sse2_128_64.c',
  'sg_diag_sse2_128_64.c',
  'sw_diag_sse2_128_64.c',

  'nw_diag_sse2_128_32.c',
  'sg_diag_sse2_128_32.c',
  'sw_diag_sse2_128_32.c',

  'nw_diag_sse2_128_16.c',
  'sg_diag_sse2_128_16.c',
  'sw_diag_sse2_128_16.c',

  'nw_diag_sse2_128_8.c',
  'sg_diag_sse2_128_8.c',
  'sw_diag_sse2_128_8.c'])

parasail_c_sse41_sources += files([
  'nw_diag_sse41_128_64.c',
  'sg_diag_sse41_128_64.c',
  'sw_diag_sse41_128_64.c',

  'nw_diag_sse41_128_32.c',
  'sg_diag_sse41_128_32.c',
  'sw_diag_sse41_128_32.c',

  'nw_diag_sse41_128_16.c',
  'sg_diag_sse41_128_16.c',
  'sw_diag_sse41_128_16.c',

  'nw_diag_sse41_128_8.c',
  'sg_diag_sse41_128_8.c',
  'sw_diag_sse41_128_8.c'])

parasail_c_avx2_sources += files([
  'nw_diag_avx2_256_64.c',
  'sg_diag_avx2_256_64.c',
  'sw_diag_avx2_256_64.c',

  'nw_diag_avx2_256_32.c',
  'sg_diag_avx2_256_32.c',
  'sw_diag_avx2_256_32.c',

  'nw_diag_avx2_256_16.c',
  'sg_diag_avx2_256_16.c',
  'sw_diag_avx2_256_16.c',

  'nw_diag_avx2_256_8.c',
  'sg_diag_avx2_256_8.c',
  'sw_diag_avx2_256_8.c'])

# parallel striped methods
parasail_c_sse2_sources += files([
  'nw_striped_sse2_128_64.c',
  'sg_striped_sse2_128_64.c',
  'sw_striped_sse2_128_64.c',

  'nw_striped_sse2_128_32.c',
  'sg_striped_sse2_128_32.c',
  'sw_striped_sse2_128_32.c',

  'nw_striped_sse2_128_16.c',
  'sg_striped_sse2_128_16.c',
  'sw_striped_sse2_128_16.c',

  'nw_striped_sse2_128_8.c',
  'sg_striped_sse2_128_8.c',
  'sw_striped_sse2_128_8.c'])

parasail_c_sse41_sources += files([
  'nw_striped_sse41_128_64.c',
  'sg_striped_sse41_128_64.c',
  'sw_striped_sse41_128_64.c',

  'nw_striped_sse41_128_32.c',
  'sg_striped_sse41_128_32.c',
  'sw_striped_sse41_128_32.c',

  'nw_striped_sse41_128_16.c',
  'sg_striped_sse41_128_16.c',
  'sw_striped_sse41_128_16.c',

  'nw_striped_sse41_128_8.c',
  'sg_striped_sse41_128_8.c',
  'sw_striped_sse41_128_8.c'])

parasail_c_avx2_sources += files([
  'nw_striped_avx2_256_64.c',
  'sg_striped_avx2_256_64.c',
  'sw_striped_avx2_256_64.c',

  'nw_striped_avx2_256_32.c',
  'sg_striped_avx2_256_32.c',
  'sw_striped_avx2_256_32.c',

  'nw_striped_avx2_256_16.c',
  'sg_striped_avx2_256_16.c',
  'sw_striped_avx2_256_16.c',

  'nw_striped_avx2_256_8.c',
  'sg_striped_avx2_256_8.c',
  'sw_striped_avx2_256_8.c'])

# parallel blocked methods
parasail_c_sse41_sources += files([
  'sw_blocked_sse41_128_32.c',
  'sw_blocked_sse41_128_16.c'])

##########################
# parallel stats methods #
##########################

# parallel scan methods
parasail_c_sse2_sources += files([
  'nw_stats_scan_sse2_128_64.c',
  'sg_stats_scan_sse2_128_64.c',
  'sw_stats_scan_sse2_128_64.c',

  'nw_stats_scan_sse2_128_32.c',
  'sg_stats_scan_sse2_128_32.c',
  'sw_stats_scan_sse2_128_32.c',

  'nw_stats_scan_sse2_128_16.c',
  'sg_stats_scan_sse2_128_16.c',
  'sw_stats_scan_sse2_128_16.c',

  'nw_stats_scan_sse2_128_8.c',
  'sg_stats_scan_sse2_128_8.c',
  'sw_stats_scan_sse2_128_8.c'])

parasail_c_sse41_sources += files([
  'nw_stats_scan_sse41_128_64.c',
  'sg_stats_scan_sse41_128_64.c',
  'sw_stats_scan_sse41_128_64.c',

  'nw_stats_scan_sse41_128_32.c',
  'sg_stats_scan_sse41_128_32.c',
  'sw_stats_scan_sse41_128_32.c',

  'nw_stats_scan_sse41_128_16.c',
  'sg_stats_scan_sse41_128_16.c',
  'sw_stats_scan_sse41_128_16.c',

  'nw_stats_scan_sse41_128_8.c',
  'sg_stats_scan_sse41_128_8.c',
  'sw_stats_scan_sse41_128_8.c'])

parasail_c_avx2_sources += files([
  'nw_stats_scan_avx2_256_64.c',
  'sg_stats_scan_avx2_256_64.c',
  'sw_stats_scan_avx2_256_64.c',

  'nw_stats_scan_avx2_256_32.c',
  'sg_stats_scan_avx2_256_32.c',
  'sw_stats_scan_avx2_256_32.c',

  'nw_stats_scan_avx2_256_16.c',
  'sg_stats_scan_avx2_256_16.c',
  'sw_stats_scan_avx2_256_16.c',

  'nw_stats_scan_avx2_256_8.c',
  'sg_stats_scan_avx2_256_8.c',
  'sw_stats_scan_avx2_256_8.c'])

# parallel diag methods
parasail_c_sse2_sources += files([
  'nw_stats_diag_sse2_128_64.c',
  'sg_stats_diag_sse2_128_64.c',
  'sw_stats_diag_sse2_128_64.c',

  'nw_stats_diag_sse2_128_32.c',
  'sg_stats_diag_sse2_128_32.c',
  'sw_stats_diag_sse2_128_32.c',

  'nw_stats_diag_sse2_128_16.c',
  'sg_stats_diag_sse2_128_16.c',
  'sw_stats_diag_sse2_128_16.c',

  'nw_stats_diag_sse2_128_8.c',
  'sg_stats_diag_sse2_128_8.c',
  'sw_stats_diag_sse2_128_8.c'])

parasail_c_sse41_sources += files([
  'nw_stats_diag_sse41_128_64.c',
  'sg_stats_diag_sse41_128_64.c',
  'sw_stats_diag_sse41_128_64.c',

  'nw_stats_diag_sse41_128_32.c',
  'sg_stats_diag_sse41_128_32.c',
  'sw_stats_diag_sse41_128_32.c',

  'nw_stats_diag_sse41_128_16.c',
  'sg_stats_diag_sse41_128_16.c',
  'sw_stats_diag_sse41_128_16.c',

  'nw_stats_diag_sse41_128_8.c',
  'sg_stats_diag_sse41_128_8.c',
  'sw_stats_diag_sse41_128_8.c'])

parasail_c_avx2_sources += files([
  'nw_stats_diag_avx2_256_64.c',
  'sg_stats_diag_avx2_256_64.c',
  'sw_stats_diag_avx2_256_64.c',

  'nw_stats_diag_avx2_256_32.c',
  'sg_stats_diag_avx2_256_32.c',
  'sw_stats_diag_avx2_256_32.c',

  'nw_stats_diag_avx2_256_16.c',
  'sg_stats_diag_avx2_256_16.c',
  'sw_stats_diag_avx2_256_16.c',

  'nw_stats_diag_avx2_256_8.c',
  'sg_stats_diag_avx2_256_8.c',
  'sw_stats_diag_avx2_256_8.c'])

# parallel striped methods
parasail_c_sse2_sources += files([
  'nw_stats_striped_sse2_128_64.c',
  'sg_stats_striped_sse2_128_64.c',
  'sw_stats_striped_sse2_128_64.c',

  'nw_stats_striped_sse2_128_32.c',
  'sg_stats_striped_sse2_128_32.c',
  'sw_stats_striped_sse2_128_32.c',

  'nw_stats_striped_sse2_128_16.c',
  'sg_stats_striped_sse2_128_16.c',
  'sw_stats_striped_sse2_128_16.c',

  'nw_stats_striped_sse2_128_8.c',
  'sg_stats_striped_sse2_128_8.c',
  'sw_stats_striped_sse2_128_8.c'])

parasail_c_sse41_sources += files([
  'nw_stats_striped_sse41_128_64.c',
  'sg_stats_striped_sse41_128_64.c',
  'sw_stats_striped_sse41_128_64.c',

  'nw_stats_striped_sse41_128_32.c',
  'sg_stats_striped_sse41_128_32.c',
  'sw_stats_striped_sse41_128_32.c',

  'nw_stats_striped_sse41_128_16.c',
  'sg_stats_striped_sse41_128_16.c',
  'sw_stats_striped_sse41_128_16.c',

  'nw_stats_striped_sse41_128_8.c',
  'sg_stats_striped_sse41_128_8.c',
  'sw_stats_striped_sse41_128_8.c'])

parasail_c_avx2_sources += files([
  'nw_stats_striped_avx2_256_64.c',
  'sg_stats_striped_avx2_256_64.c',
  'sw_stats_striped_avx2_256_64.c',

  'nw_stats_striped_avx2_256_32.c',
  'sg_stats_striped_avx2_256_32.c',
  'sw_stats_striped_avx2_256_32.c',

  'nw_stats_striped_avx2_256_16.c',
  'sg_stats_striped_avx2_256_16.c',
  'sw_stats_striped_avx2_256_16.c',

  'nw_stats_striped_avx2_256_8.c',
  'sg_stats_striped_avx2_256_8.c',
  'sw_stats_striped_avx2_256_8.c'])

##########################
# parallel trace methods #
##########################

parasail_c_trace_sse2_sources = []
parasail_c_trace_sse41_sources = []
parasail_c_trace_avx2_sources = []

# parallel scan methods
parasail_c_trace_sse2_sources += files([
  'nw_trace_scan_sse2_128_64.c',
  'sg_trace_scan_sse2_128_64.c',
  'sw_trace_scan_sse2_128_64.c',

  'nw_trace_scan_sse2_128_32.c',
  'sg_trace_scan_sse2_128_32.c',
  'sw_trace_scan_sse2_128_32.c',

  'nw_trace_scan_sse2_128_16.c',
  'sg_trace_scan_sse2_128_16.c',
  'sw_trace_scan_sse2_128_16.c',

  'nw_trace_scan_sse2_128_8.c',
  'sg_trace_scan_sse2_128_8.c',
  'sw_trace_scan_sse2_128_8.c'])

parasail_c_trace_sse41_sources += files([
  'nw_trace_scan_sse41_128_64.c',
  'sg_trace_scan_sse41_128_64.c',
  'sw_trace_scan_sse41_128_64.c',

  'nw_trace_scan_sse41_128_32.c',
  'sg_trace_scan_sse41_128_32.c',
  'sw_trace_scan_sse41_128_32.c',

  'nw_trace_scan_sse41_128_16.c',
  'sg_trace_scan_sse41_128_16.c',
  'sw_trace_scan_sse41_128_16.c',

  'nw_trace_scan_sse41_128_8.c',
  'sg_trace_scan_sse41_128_8.c',
  'sw_trace_scan_sse41_128_8.c'])

parasail_c_trace_avx2_sources += files([
  'nw_trace_scan_avx2_256_64.c',
  'sg_trace_scan_avx2_256_64.c',
  'sw_trace_scan_avx2_256_64.c',

  'nw_trace_scan_avx2_256_32.c',
  'sg_trace_scan_avx2_256_32.c',
  'sw_trace_scan_avx2_256_32.c',

  'nw_trace_scan_avx2_256_16.c',
  'sg_trace_scan_avx2_256_16.c',
  'sw_trace_scan_avx2_256_16.c',

  'nw_trace_scan_avx2_256_8.c',
  'sg_trace_scan_avx2_256_8.c',
  'sw_trace_scan_avx2_256_8.c'])

# parallel diag methods
parasail_c_trace_sse2_sources += files([
  'nw_trace_diag_sse2_128_64.c',
  'sg_trace_diag_sse2_128_64.c',
  'sw_trace_diag_sse2_128_64.c',

  'nw_trace_diag_sse2_128_32.c',
  'sg_trace_diag_sse2_128_32.c',
  'sw_trace_diag_sse2_128_32.c',

  'nw_trace_diag_sse2_128_16.c',
  'sg_trace_diag_sse2_128_16.c',
  'sw_trace_diag_sse2_128_16.c',

  'nw_trace_diag_sse2_128_8.c',
  'sg_trace_diag_sse2_128_8.c',
  'sw_trace_diag_sse2_128_8.c'])

parasail_c_trace_sse41_sources += files([
  'nw_trace_diag_sse41_128_64.c',
  'sg_trace_diag_sse41_128_64.c',
  'sw_trace_diag_sse41_128_64.c',

  'nw_trace_diag_sse41_128_32.c',
  'sg_trace_diag_sse41_128_32.c',
  'sw_trace_diag_sse41_128_32.c',

  'nw_trace_diag_sse41_128_16.c',
  'sg_trace_diag_sse41_128_16.c',
  'sw_trace_diag_sse41_128_16.c',

  'nw_trace_diag_sse41_128_8.c',
  'sg_trace_diag_sse41_128_8.c',
  'sw_trace_diag_sse41_128_8.c'])

parasail_c_trace_avx2_sources += files([
  'nw_trace_diag_avx2_256_64.c',
  'sg_trace_diag_avx2_256_64.c',
  'sw_trace_diag_avx2_256_64.c',

  'nw_trace_diag_avx2_256_32.c',
  'sg_trace_diag_avx2_256_32.c',
  'sw_trace_diag_avx2_256_32.c',

  'nw_trace_diag_avx2_256_16.c',
  'sg_trace_diag_avx2_256_16.c',
  'sw_trace_diag_avx2_256_16.c',

  'nw_trace_diag_avx2_256_8.c',
  'sg_trace_diag_avx2_256_8.c',
  'sw_trace_diag_avx2_256_8.c'])

# parallel striped methods
parasail_c_trace_sse2_sources += files([
  'nw_trace_striped_sse2_128_64.c',
  'sg_trace_striped_sse2_128_64.c',
  'sw_trace_striped_sse2_128_64.c',

  'nw_trace_striped_sse2_128_32.c',
  'sg_trace_striped_sse2_128_32.c',
  'sw_trace_striped_sse2_128_32.c',

  'nw_trace_striped_sse2_128_16.c',
  'sg_trace_striped_sse2_128_16.c',
  'sw_trace_striped_sse2_128_16.c',

  'nw_trace_striped_sse2_128_8.c',
  'sg_trace_striped_sse2_128_8.c',
  'sw_trace_striped_sse2_128_8.c'])

parasail_c_trace_sse41_sources += files([
  'nw_trace_striped_sse41_128_64.c',
  'sg_trace_striped_sse41_128_64.c',
  'sw_trace_striped_sse41_128_64.c',

  'nw_trace_striped_sse41_128_32.c',
  'sg_trace_striped_sse41_128_32.c',
  'sw_trace_striped_sse41_128_32.c',

  'nw_trace_striped_sse41_128_16.c',
  'sg_trace_striped_sse41_128_16.c',
  'sw_trace_striped_sse41_128_16.c',

  'nw_trace_striped_sse41_128_8.c',
  'sg_trace_striped_sse41_128_8.c',
  'sw_trace_striped_sse41_128_8.c'])

parasail_c_trace_avx2_sources += files([
  'nw_trace_striped_avx2_256_64.c',
  'sg_trace_striped_avx2_256_64.c',
  'sw_trace_striped_avx2_256_64.c',

  'nw_trace_striped_avx2_256_32.c',
  'sg_trace_striped_avx2_256_32.c',
  'sw_trace_striped_avx2_256_32.c',

  'nw_trace_striped_avx2_256_16.c',
  'sg_trace_striped_avx2_256_16.c',
  'sw_trace_striped_avx2_256_16.c',

  'nw_trace_striped_avx2_256_8.c',
  'sg_trace_striped_avx2_256_8.c',
  'sw_trace_striped_avx2_256_8.c'])

#####################
# memory allocation #
#####################

parasail_config = configuration_data()
parasail_config.set('SIZEOF_INT', parasail_cc.sizeof('int'))

foreach funcname : [
  '_aligned_malloc',
  'aligned_alloc',
  'memalign',
  'memset',
  'posix_memalign']
    if parasail_cc.has_function(funcname)
      result = 1
    else
      result = false
    endif
    parasail_config.set('HAVE_' + funcname.to_upper(), result)
endforeach

#########################
# convenience libraries #
#########################

parasail_simd = import('unstable_simd')

parasail_simd_variants = [
  'avx2',
  'sse41',
  'sse2',
  'novec']

parasail_simd_sources = [
  parasail_c_avx2_sources,
  parasail_c_sse41_sources,
  parasail_c_sse2_sources,
  parasail_c_novec_sources]

parasail_simd_macros = [
  [
    'HAVE_AVX2',
    'HAVE_AVX2_MM256_EXTRACT_EPI16',
    'HAVE_AVX2_MM256_EXTRACT_EPI32',
    'HAVE_AVX2_MM256_EXTRACT_EPI64',
    'HAVE_AVX2_MM256_EXTRACT_EPI8',
    'HAVE_AVX2_MM256_INSERT_EPI16',
    'HAVE_AVX2_MM256_INSERT_EPI32',
    'HAVE_AVX2_MM256_INSERT_EPI64',
    'HAVE_AVX2_MM256_INSERT_EPI8',
    'HAVE_AVX2_MM256_SET1_EPI64X',
    'HAVE_AVX2_MM256_SET_EPI64X'],
  [
    'HAVE_SSE41',
    'HAVE_SSE41_MM_EXTRACT_EPI64',
    'HAVE_SSE41_MM_INSERT_EPI64'],
  [
    'HAVE_SSE2',
    'HAVE_SSE2_MM_SET1_EPI64X',
    'HAVE_SSE2_MM_SET_EPI64X'],
  []]

parasail_build_variants = [
  ['plain', [], [
    parasail_c_trace_avx2_sources,
    parasail_c_trace_sse41_sources,
    parasail_c_trace_sse2_sources,
    parasail_c_trace_novec_sources]],
  ['table', ['-DPARASAIL_TABLE=1'], []],
  ['rowcol', ['-DPARASAIL_ROWCOL=1'], []]]

parasail_individual_libs = []

foreach j : parasail_build_variants
  simd_source_files = []

  foreach i : [0, 1, 2, 3]
    simd_name = parasail_simd_variants[i]
    extra_sources = []

    if simd_name == 'novec'
      # Non-Vectorized
      if j[0] == 'plain'
        extra_sources = j[2][i]
      endif

      parasail_individual_libs += static_library(
        simd_name + '_' + j[0],
        parasail_simd_sources[i] + extra_sources,
        include_directories : parasail_include_directories,
        c_args : j[1])
    else
      # SIMD
      if get_option(simd_name)
        if j[0] == 'plain'
          extra_sources = j[2][i]

          if simd_name.startswith('sse')
            extra_sources += parasail_c_mem_sse_sources
            parasail_c_mem_sse_sources = []
          endif
          if simd_name.startswith('avx')
            extra_sources += parasail_c_mem_avx2_sources
            parasail_c_mem_avx2_sources = []
          endif
        endif

        simd_source_files += [[parasail_simd_sources[i], extra_sources]]
        result = 1
      else
        # currently Meson has a bug where an empty
        # file list argument to simd.check() generates
        # empty libraries, which in turn causes the
        # macOS linker to trip over.
        # By including a dummy file with a static
        # function, we circumvent this issue until
        # it is fixed in Meson upstream.
        simd_source_files += [files('dummy.c')]
        result = false
      endif
    endif

    foreach k : parasail_simd_macros[i]
      parasail_config.set(k, result)
    endforeach
  endforeach

  message('Configuring ' + j[0] + ' variant')

  parasail_individual_libs += parasail_simd.check(
    simd_name + '_' + j[0],
    avx2 : simd_source_files[0],
    sse41 : simd_source_files[1],
    sse2 : simd_source_files[2],
    include_directories : parasail_include_directories,
    c_args : j[1],
    compiler : parasail_cc)[0]
endforeach

# try to find restrict support in C++
# precedence order taken from AC_C_RESTRICT
parasail_attempt_restrict_variants = [
  '__restrict',
  '__restrict__',
  '_Restrict',
  'restrict']

parasail_working_restrict_variants = []
foreach i : parasail_attempt_restrict_variants
  if parasail_cxx.compiles('int f(int* ' + i + ' p) { return p[0]; }', name : 'C++ compiler supports ' + i + ' keyword')
    parasail_working_restrict_variants+= [i]
  endif
endforeach

if parasail_working_restrict_variants.length() > 0
  parasail_config.set('restrict', parasail_working_restrict_variants[0])
else
  parasail_config.set('restrict', '')
endif

parasail_config_h = configure_file(
  output : 'config.h',
  configuration : parasail_config)

################
# main library #
################

# install library if
# - either running as a proper project
# - or using shared libraries
parasail_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

parasail_lib = library(
  'parasail', [
    parasail_config_h,
    parasail_c_core_sources],
  soversion : 3,
  version : meson.project_version(),
  install : parasail_lib_install,
  dependencies : [parasail_zlib_dep, parasail_libm_dep],
  link_whole : parasail_individual_libs,
  include_directories : parasail_include_directories)
