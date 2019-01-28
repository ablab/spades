/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_neon.h"

simde__m128i * parasail_memalign_simde__m128i(size_t alignment, size_t size)
{
    return (simde__m128i *) parasail_memalign(alignment, size*sizeof(simde__m128i));
}

void parasail_memset_simde__m128i(simde__m128i *b, simde__m128i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        simde_mm_store_si128(&b[i], c);
    }
}

parasail_profile_t * parasail_profile_create_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.i8[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], t);
            ++index;
        }
    }

    profile->profile8.score = vProfile;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 8; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.i16[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], t);
            ++index;
        }
    }

    profile->profile16.score = vProfile;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 4; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.i32[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], t);
            ++index;
        }
    }

    profile->profile32.score = vProfile;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 2; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.i64[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], t);
            ++index;
        }
    }

    profile->profile64.score = vProfile;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t* parasail_profile_create_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile8 = parasail_profile_create_neon_128_8(s1, s1Len, matrix);
    parasail_profile_t *profile16 = parasail_profile_create_neon_128_16(s1, s1Len, matrix);
    parasail_profile_t *profile32 = parasail_profile_create_neon_128_32(s1, s1Len, matrix);
    profile8->profile16 = profile16->profile16;
    profile8->profile32 = profile32->profile32;
    free(profile16);
    free(profile32);

    return profile8;
}

parasail_profile_t * parasail_profile_create_stats_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileM = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileS = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i p;
            simde__m128i m;
            simde__m128i s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p.i8[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m.i8[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s.i8[segNum] = p.i8[segNum] > 0;
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], p);
            simde_mm_store_si128(&vProfileM[index], m);
            simde_mm_store_si128(&vProfileS[index], s);
            ++index;
        }
    }

    profile->profile8.score = vProfile;
    profile->profile8.matches = vProfileM;
    profile->profile8.similar = vProfileS;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 8; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileM = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileS = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i p;
            simde__m128i m;
            simde__m128i s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p.i16[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m.i16[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s.i16[segNum] = p.i16[segNum] > 0;
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], p);
            simde_mm_store_si128(&vProfileM[index], m);
            simde_mm_store_si128(&vProfileS[index], s);
            ++index;
        }
    }

    profile->profile16.score = vProfile;
    profile->profile16.matches = vProfileM;
    profile->profile16.similar = vProfileS;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 4; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileM = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileS = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i p;
            simde__m128i m;
            simde__m128i s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p.i32[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m.i32[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s.i32[segNum] = p.i32[segNum] > 0;
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], p);
            simde_mm_store_si128(&vProfileM[index], m);
            simde_mm_store_si128(&vProfileS[index], s);
            ++index;
        }
    }

    profile->profile32.score = vProfile;
    profile->profile32.matches = vProfileM;
    profile->profile32.similar = vProfileS;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 2; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileM = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileS = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i p;
            simde__m128i m;
            simde__m128i s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p.i64[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m.i64[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s.i64[segNum] = p.i64[segNum] > 0;
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], p);
            simde_mm_store_si128(&vProfileM[index], m);
            simde_mm_store_si128(&vProfileS[index], s);
            ++index;
        }
    }

    profile->profile64.score = vProfile;
    profile->profile64.matches = vProfileM;
    profile->profile64.similar = vProfileS;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t* parasail_profile_create_stats_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile8 = parasail_profile_create_stats_neon_128_8(s1, s1Len, matrix);
    parasail_profile_t *profile16 = parasail_profile_create_stats_neon_128_16(s1, s1Len, matrix);
    parasail_profile_t *profile32 = parasail_profile_create_stats_neon_128_32(s1, s1Len, matrix);
    profile8->profile16 = profile16->profile16;
    profile8->profile32 = profile32->profile32;
    free(profile16);
    free(profile32);

    return profile8;
}

void parasail_free_simde__m128i(void *ptr)
{
    parasail_free((simde__m128i*)ptr);
}

