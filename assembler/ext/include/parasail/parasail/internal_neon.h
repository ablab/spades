/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyb (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_INTERNAL_NEON_H_
#define _PARASAIL_INTERNAL_NEON_H_

#include <stdint.h>

/* we use simde neon equivalents up to sse4.1 */
#include "x86/sse4.1.h"

#ifdef __cplusplus
extern "C" {
#endif

extern simde__m128i * parasail_memalign_simde__m128i(size_t alignment, size_t size);

extern void parasail_memset_simde__m128i(simde__m128i *b, simde__m128i c, size_t len);

extern void parasail_free_simde__m128i(void *ptr);

/* this function is part of sse4.2 */
SIMDE__FUNCTION_ATTRIBUTES
simde__m128i
simde_mm_cmpgt_epi64 (simde__m128i a, simde__m128i b) {
  simde__m128i r;
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.i64) / sizeof(r.i64[0])) ; i++) {
    r.i64[i] = (a.i64[i] > b.i64[i]) ? 0xff : 0x00;
  }
  return r;
}

/* this function does not exist, but is easy to make */
SIMDE__FUNCTION_ATTRIBUTES
simde__m128i
simde_mm_cmplt_epi64 (simde__m128i a, simde__m128i b) {
  simde__m128i r;
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.i64) / sizeof(r.i64[0])) ; i++) {
    r.i64[i] = (a.i64[i] < b.i64[i]) ? 0xff : 0x00;
  }
  return r;
}

/* this function is part of AVX512VL + AVX512F */
SIMDE__FUNCTION_ATTRIBUTES
simde__m128i
simde_mm_max_epi64 (simde__m128i a, simde__m128i b) {
  simde__m128i r;
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.i64) / sizeof(r.i64[0])) ; i++) {
    r.i64[i] = (a.i64[i] > b.i64[i]) ? a.i64[i] : b.i64[i];
  }
  return r;
}

/* this function is part of AVX512VL + AVX512F */
SIMDE__FUNCTION_ATTRIBUTES
simde__m128i
simde_mm_min_epi64 (simde__m128i a, simde__m128i b) {
  simde__m128i r;
  SIMDE__VECTORIZE
  for (size_t i = 0 ; i < (sizeof(r.i64) / sizeof(r.i64[0])) ; i++) {
    r.i64[i] = (a.i64[i] < b.i64[i]) ? a.i64[i] : b.i64[i];
  }
  return r;
}

/* "hmax" is made up, but used */

static inline int8_t simde_mm_hmax_epi8(simde__m128i v)
{
    v = simde_mm_max_epi8(v, simde_mm_srli_si128(v, 8));
    v = simde_mm_max_epi8(v, simde_mm_srli_si128(v, 4));
    v = simde_mm_max_epi8(v, simde_mm_srli_si128(v, 2));
    v = simde_mm_max_epi8(v, simde_mm_srli_si128(v, 1));
    return simde_mm_extract_epi8(v, 0);
}

static inline int16_t simde_mm_hmax_epi16(simde__m128i v)
{
    v = simde_mm_max_epi16(v, simde_mm_srli_si128(v, 8));
    v = simde_mm_max_epi16(v, simde_mm_srli_si128(v, 4));
    v = simde_mm_max_epi16(v, simde_mm_srli_si128(v, 2));
    return simde_mm_extract_epi16(v, 0);
}

static inline int32_t simde_mm_hmax_epi32(simde__m128i v)
{
    v = simde_mm_max_epi32(v, simde_mm_srli_si128(v, 8));
    v = simde_mm_max_epi32(v, simde_mm_srli_si128(v, 4));
    return simde_mm_extract_epi32(v, 0);
}

static inline int64_t simde_mm_hmax_epi64(simde__m128i v)
{
    v = simde_mm_max_epi64(v, simde_mm_srli_si128(v, 8));
    return simde_mm_extract_epi64(v, 0);
}

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_INTERNAL_NEON_H_ */
