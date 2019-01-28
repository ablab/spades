/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_INTERNAL_SSE_H_
#define _PARASAIL_INTERNAL_SSE_H_

#include <stdint.h>

#include <emmintrin.h>

#include "parasail/parasail.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef union __m128i_8 {
    __m128i m;
    int8_t v[16];
} __m128i_8_t;

typedef union __m128i_16 {
    __m128i m;
    int16_t v[8];
} __m128i_16_t;

typedef union __m128i_32 {
    __m128i m;
    int32_t v[4];
} __m128i_32_t;

typedef union __m128i_64 {
    __m128i m;
    int64_t v[2];
} __m128i_64_t;

extern
__m128i * parasail_memalign___m128i(size_t alignment, size_t size);

extern
void parasail_memset___m128i(__m128i *b, __m128i c, size_t len);

extern
void parasail_free___m128i(void *ptr);

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_INTERNAL_SSE_H_ */
