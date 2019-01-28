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
#include "parasail/parasail/internal_altivec.h"

vec128i * parasail_memalign_vec128i(size_t alignment, size_t size)
{
    return (vec128i *) parasail_memalign(alignment, size*sizeof(vec128i));
}

void parasail_memset_vec128i(vec128i *b, vec128i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm_store_si128(&b[i], c);
    }
}

parasail_profile_t * parasail_profile_create_altivec_128_8(
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
    vec128i* const restrict vProfile = parasail_memalign_vec128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            vec128i_8_t t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            _mm_store_si128(&vProfile[index], t.m);
            ++index;
        }
    }

    profile->profile8.score = vProfile;
    profile->free = &parasail_free_vec128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_altivec_128_16(
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
    vec128i* const restrict vProfile = parasail_memalign_vec128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            vec128i_16_t t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            _mm_store_si128(&vProfile[index], t.m);
            ++index;
        }
    }

    profile->profile16.score = vProfile;
    profile->free = &parasail_free_vec128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_altivec_128_32(
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
    vec128i* const restrict vProfile = parasail_memalign_vec128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            vec128i_32_t t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            _mm_store_si128(&vProfile[index], t.m);
            ++index;
        }
    }

    profile->profile32.score = vProfile;
    profile->free = &parasail_free_vec128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_altivec_128_64(
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
    vec128i* const restrict vProfile = parasail_memalign_vec128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            vec128i_64_t t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            _mm_store_si128(&vProfile[index], t.m);
            ++index;
        }
    }

    profile->profile64.score = vProfile;
    profile->free = &parasail_free_vec128i;
    return profile;
}

parasail_profile_t* parasail_profile_create_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile8 = parasail_profile_create_altivec_128_8(s1, s1Len, matrix);
    parasail_profile_t *profile16 = parasail_profile_create_altivec_128_16(s1, s1Len, matrix);
    parasail_profile_t *profile32 = parasail_profile_create_altivec_128_32(s1, s1Len, matrix);
    profile8->profile16 = profile16->profile16;
    profile8->profile32 = profile32->profile32;
    free(profile16);
    free(profile32);

    return profile8;
}

parasail_profile_t * parasail_profile_create_stats_altivec_128_8(
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
    vec128i* const restrict vProfile = parasail_memalign_vec128i(16, n * segLen);
    vec128i* const restrict vProfileM = parasail_memalign_vec128i(16, n * segLen);
    vec128i* const restrict vProfileS = parasail_memalign_vec128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            vec128i_8_t p;
            vec128i_8_t m;
            vec128i_8_t s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m.v[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s.v[segNum] = p.v[segNum] > 0;
                j += segLen;
            }
            _mm_store_si128(&vProfile[index], p.m);
            _mm_store_si128(&vProfileM[index], m.m);
            _mm_store_si128(&vProfileS[index], s.m);
            ++index;
        }
    }

    profile->profile8.score = vProfile;
    profile->profile8.matches = vProfileM;
    profile->profile8.similar = vProfileS;
    profile->free = &parasail_free_vec128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_altivec_128_16(
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
    vec128i* const restrict vProfile = parasail_memalign_vec128i(16, n * segLen);
    vec128i* const restrict vProfileM = parasail_memalign_vec128i(16, n * segLen);
    vec128i* const restrict vProfileS = parasail_memalign_vec128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            vec128i_16_t p;
            vec128i_16_t m;
            vec128i_16_t s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m.v[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s.v[segNum] = p.v[segNum] > 0;
                j += segLen;
            }
            _mm_store_si128(&vProfile[index], p.m);
            _mm_store_si128(&vProfileM[index], m.m);
            _mm_store_si128(&vProfileS[index], s.m);
            ++index;
        }
    }

    profile->profile16.score = vProfile;
    profile->profile16.matches = vProfileM;
    profile->profile16.similar = vProfileS;
    profile->free = &parasail_free_vec128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_altivec_128_32(
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
    vec128i* const restrict vProfile = parasail_memalign_vec128i(16, n * segLen);
    vec128i* const restrict vProfileM = parasail_memalign_vec128i(16, n * segLen);
    vec128i* const restrict vProfileS = parasail_memalign_vec128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            vec128i_32_t p;
            vec128i_32_t m;
            vec128i_32_t s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m.v[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s.v[segNum] = p.v[segNum] > 0;
                j += segLen;
            }
            _mm_store_si128(&vProfile[index], p.m);
            _mm_store_si128(&vProfileM[index], m.m);
            _mm_store_si128(&vProfileS[index], s.m);
            ++index;
        }
    }

    profile->profile32.score = vProfile;
    profile->profile32.matches = vProfileM;
    profile->profile32.similar = vProfileS;
    profile->free = &parasail_free_vec128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_altivec_128_64(
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
    vec128i* const restrict vProfile = parasail_memalign_vec128i(16, n * segLen);
    vec128i* const restrict vProfileM = parasail_memalign_vec128i(16, n * segLen);
    vec128i* const restrict vProfileS = parasail_memalign_vec128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            vec128i_64_t p;
            vec128i_64_t m;
            vec128i_64_t s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m.v[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s.v[segNum] = p.v[segNum] > 0;
                j += segLen;
            }
            _mm_store_si128(&vProfile[index], p.m);
            _mm_store_si128(&vProfileM[index], m.m);
            _mm_store_si128(&vProfileS[index], s.m);
            ++index;
        }
    }

    profile->profile64.score = vProfile;
    profile->profile64.matches = vProfileM;
    profile->profile64.similar = vProfileS;
    profile->free = &parasail_free_vec128i;
    return profile;
}

parasail_profile_t* parasail_profile_create_stats_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile8 = parasail_profile_create_stats_altivec_128_8(s1, s1Len, matrix);
    parasail_profile_t *profile16 = parasail_profile_create_stats_altivec_128_16(s1, s1Len, matrix);
    parasail_profile_t *profile32 = parasail_profile_create_stats_altivec_128_32(s1, s1Len, matrix);
    profile8->profile16 = profile16->profile16;
    profile8->profile32 = profile32->profile32;
    free(profile16);
    free(profile32);

    return profile8;
}

void parasail_free_vec128i(void *ptr)
{
    parasail_free((vec128i*)ptr);
}

