/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

#include <immintrin.h>

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_avx.h"

#define NEG_INF INT8_MIN

#if HAVE_AVX2_MM256_INSERT_EPI8
#define _mm256_insert_epi8_rpl _mm256_insert_epi8
#else
static inline __m256i _mm256_insert_epi8_rpl(__m256i a, int8_t i, int imm) {
    __m256i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI8
#define _mm256_extract_epi8_rpl _mm256_extract_epi8
#else
static inline int8_t _mm256_extract_epi8_rpl(__m256i a, int imm) {
    __m256i_8_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_srli_si256_rpl(a,imm) _mm256_or_si256(_mm256_slli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)), 16-imm), _mm256_srli_si256(a, imm))

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 31);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 30);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 29);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 28);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[1LL*(i+4)*s2Len + (j-4)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 27);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[1LL*(i+5)*s2Len + (j-5)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 26);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[1LL*(i+6)*s2Len + (j-6)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 25);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[1LL*(i+7)*s2Len + (j-7)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 24);
    }
    if (0 <= i+8 && i+8 < s1Len && 0 <= j-8 && j-8 < s2Len) {
        array[1LL*(i+8)*s2Len + (j-8)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 23);
    }
    if (0 <= i+9 && i+9 < s1Len && 0 <= j-9 && j-9 < s2Len) {
        array[1LL*(i+9)*s2Len + (j-9)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 22);
    }
    if (0 <= i+10 && i+10 < s1Len && 0 <= j-10 && j-10 < s2Len) {
        array[1LL*(i+10)*s2Len + (j-10)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 21);
    }
    if (0 <= i+11 && i+11 < s1Len && 0 <= j-11 && j-11 < s2Len) {
        array[1LL*(i+11)*s2Len + (j-11)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 20);
    }
    if (0 <= i+12 && i+12 < s1Len && 0 <= j-12 && j-12 < s2Len) {
        array[1LL*(i+12)*s2Len + (j-12)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 19);
    }
    if (0 <= i+13 && i+13 < s1Len && 0 <= j-13 && j-13 < s2Len) {
        array[1LL*(i+13)*s2Len + (j-13)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 18);
    }
    if (0 <= i+14 && i+14 < s1Len && 0 <= j-14 && j-14 < s2Len) {
        array[1LL*(i+14)*s2Len + (j-14)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 17);
    }
    if (0 <= i+15 && i+15 < s1Len && 0 <= j-15 && j-15 < s2Len) {
        array[1LL*(i+15)*s2Len + (j-15)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 16);
    }
    if (0 <= i+16 && i+16 < s1Len && 0 <= j-16 && j-16 < s2Len) {
        array[1LL*(i+16)*s2Len + (j-16)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 15);
    }
    if (0 <= i+17 && i+17 < s1Len && 0 <= j-17 && j-17 < s2Len) {
        array[1LL*(i+17)*s2Len + (j-17)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 14);
    }
    if (0 <= i+18 && i+18 < s1Len && 0 <= j-18 && j-18 < s2Len) {
        array[1LL*(i+18)*s2Len + (j-18)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 13);
    }
    if (0 <= i+19 && i+19 < s1Len && 0 <= j-19 && j-19 < s2Len) {
        array[1LL*(i+19)*s2Len + (j-19)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 12);
    }
    if (0 <= i+20 && i+20 < s1Len && 0 <= j-20 && j-20 < s2Len) {
        array[1LL*(i+20)*s2Len + (j-20)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 11);
    }
    if (0 <= i+21 && i+21 < s1Len && 0 <= j-21 && j-21 < s2Len) {
        array[1LL*(i+21)*s2Len + (j-21)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 10);
    }
    if (0 <= i+22 && i+22 < s1Len && 0 <= j-22 && j-22 < s2Len) {
        array[1LL*(i+22)*s2Len + (j-22)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 9);
    }
    if (0 <= i+23 && i+23 < s1Len && 0 <= j-23 && j-23 < s2Len) {
        array[1LL*(i+23)*s2Len + (j-23)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 8);
    }
    if (0 <= i+24 && i+24 < s1Len && 0 <= j-24 && j-24 < s2Len) {
        array[1LL*(i+24)*s2Len + (j-24)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 7);
    }
    if (0 <= i+25 && i+25 < s1Len && 0 <= j-25 && j-25 < s2Len) {
        array[1LL*(i+25)*s2Len + (j-25)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 6);
    }
    if (0 <= i+26 && i+26 < s1Len && 0 <= j-26 && j-26 < s2Len) {
        array[1LL*(i+26)*s2Len + (j-26)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 5);
    }
    if (0 <= i+27 && i+27 < s1Len && 0 <= j-27 && j-27 < s2Len) {
        array[1LL*(i+27)*s2Len + (j-27)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 4);
    }
    if (0 <= i+28 && i+28 < s1Len && 0 <= j-28 && j-28 < s2Len) {
        array[1LL*(i+28)*s2Len + (j-28)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 3);
    }
    if (0 <= i+29 && i+29 < s1Len && 0 <= j-29 && j-29 < s2Len) {
        array[1LL*(i+29)*s2Len + (j-29)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 2);
    }
    if (0 <= i+30 && i+30 < s1Len && 0 <= j-30 && j-30 < s2Len) {
        array[1LL*(i+30)*s2Len + (j-30)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 1);
    }
    if (0 <= i+31 && i+31 < s1Len && 0 <= j-31 && j-31 < s2Len) {
        array[1LL*(i+31)*s2Len + (j-31)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        __m256i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int8_t)_mm256_extract_epi8_rpl(vWH, 31);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 31);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int8_t)_mm256_extract_epi8_rpl(vWH, 30);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 30);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int8_t)_mm256_extract_epi8_rpl(vWH, 29);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 29);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int8_t)_mm256_extract_epi8_rpl(vWH, 28);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 28);
    }
    if (i+4 == s1Len-1 && 0 <= j-4 && j-4 < s2Len) {
        row[j-4] = (int8_t)_mm256_extract_epi8_rpl(vWH, 27);
    }
    if (j-4 == s2Len-1 && 0 <= i+4 && i+4 < s1Len) {
        col[(i+4)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 27);
    }
    if (i+5 == s1Len-1 && 0 <= j-5 && j-5 < s2Len) {
        row[j-5] = (int8_t)_mm256_extract_epi8_rpl(vWH, 26);
    }
    if (j-5 == s2Len-1 && 0 <= i+5 && i+5 < s1Len) {
        col[(i+5)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 26);
    }
    if (i+6 == s1Len-1 && 0 <= j-6 && j-6 < s2Len) {
        row[j-6] = (int8_t)_mm256_extract_epi8_rpl(vWH, 25);
    }
    if (j-6 == s2Len-1 && 0 <= i+6 && i+6 < s1Len) {
        col[(i+6)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 25);
    }
    if (i+7 == s1Len-1 && 0 <= j-7 && j-7 < s2Len) {
        row[j-7] = (int8_t)_mm256_extract_epi8_rpl(vWH, 24);
    }
    if (j-7 == s2Len-1 && 0 <= i+7 && i+7 < s1Len) {
        col[(i+7)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 24);
    }
    if (i+8 == s1Len-1 && 0 <= j-8 && j-8 < s2Len) {
        row[j-8] = (int8_t)_mm256_extract_epi8_rpl(vWH, 23);
    }
    if (j-8 == s2Len-1 && 0 <= i+8 && i+8 < s1Len) {
        col[(i+8)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 23);
    }
    if (i+9 == s1Len-1 && 0 <= j-9 && j-9 < s2Len) {
        row[j-9] = (int8_t)_mm256_extract_epi8_rpl(vWH, 22);
    }
    if (j-9 == s2Len-1 && 0 <= i+9 && i+9 < s1Len) {
        col[(i+9)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 22);
    }
    if (i+10 == s1Len-1 && 0 <= j-10 && j-10 < s2Len) {
        row[j-10] = (int8_t)_mm256_extract_epi8_rpl(vWH, 21);
    }
    if (j-10 == s2Len-1 && 0 <= i+10 && i+10 < s1Len) {
        col[(i+10)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 21);
    }
    if (i+11 == s1Len-1 && 0 <= j-11 && j-11 < s2Len) {
        row[j-11] = (int8_t)_mm256_extract_epi8_rpl(vWH, 20);
    }
    if (j-11 == s2Len-1 && 0 <= i+11 && i+11 < s1Len) {
        col[(i+11)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 20);
    }
    if (i+12 == s1Len-1 && 0 <= j-12 && j-12 < s2Len) {
        row[j-12] = (int8_t)_mm256_extract_epi8_rpl(vWH, 19);
    }
    if (j-12 == s2Len-1 && 0 <= i+12 && i+12 < s1Len) {
        col[(i+12)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 19);
    }
    if (i+13 == s1Len-1 && 0 <= j-13 && j-13 < s2Len) {
        row[j-13] = (int8_t)_mm256_extract_epi8_rpl(vWH, 18);
    }
    if (j-13 == s2Len-1 && 0 <= i+13 && i+13 < s1Len) {
        col[(i+13)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 18);
    }
    if (i+14 == s1Len-1 && 0 <= j-14 && j-14 < s2Len) {
        row[j-14] = (int8_t)_mm256_extract_epi8_rpl(vWH, 17);
    }
    if (j-14 == s2Len-1 && 0 <= i+14 && i+14 < s1Len) {
        col[(i+14)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 17);
    }
    if (i+15 == s1Len-1 && 0 <= j-15 && j-15 < s2Len) {
        row[j-15] = (int8_t)_mm256_extract_epi8_rpl(vWH, 16);
    }
    if (j-15 == s2Len-1 && 0 <= i+15 && i+15 < s1Len) {
        col[(i+15)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 16);
    }
    if (i+16 == s1Len-1 && 0 <= j-16 && j-16 < s2Len) {
        row[j-16] = (int8_t)_mm256_extract_epi8_rpl(vWH, 15);
    }
    if (j-16 == s2Len-1 && 0 <= i+16 && i+16 < s1Len) {
        col[(i+16)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 15);
    }
    if (i+17 == s1Len-1 && 0 <= j-17 && j-17 < s2Len) {
        row[j-17] = (int8_t)_mm256_extract_epi8_rpl(vWH, 14);
    }
    if (j-17 == s2Len-1 && 0 <= i+17 && i+17 < s1Len) {
        col[(i+17)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 14);
    }
    if (i+18 == s1Len-1 && 0 <= j-18 && j-18 < s2Len) {
        row[j-18] = (int8_t)_mm256_extract_epi8_rpl(vWH, 13);
    }
    if (j-18 == s2Len-1 && 0 <= i+18 && i+18 < s1Len) {
        col[(i+18)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 13);
    }
    if (i+19 == s1Len-1 && 0 <= j-19 && j-19 < s2Len) {
        row[j-19] = (int8_t)_mm256_extract_epi8_rpl(vWH, 12);
    }
    if (j-19 == s2Len-1 && 0 <= i+19 && i+19 < s1Len) {
        col[(i+19)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 12);
    }
    if (i+20 == s1Len-1 && 0 <= j-20 && j-20 < s2Len) {
        row[j-20] = (int8_t)_mm256_extract_epi8_rpl(vWH, 11);
    }
    if (j-20 == s2Len-1 && 0 <= i+20 && i+20 < s1Len) {
        col[(i+20)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 11);
    }
    if (i+21 == s1Len-1 && 0 <= j-21 && j-21 < s2Len) {
        row[j-21] = (int8_t)_mm256_extract_epi8_rpl(vWH, 10);
    }
    if (j-21 == s2Len-1 && 0 <= i+21 && i+21 < s1Len) {
        col[(i+21)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 10);
    }
    if (i+22 == s1Len-1 && 0 <= j-22 && j-22 < s2Len) {
        row[j-22] = (int8_t)_mm256_extract_epi8_rpl(vWH, 9);
    }
    if (j-22 == s2Len-1 && 0 <= i+22 && i+22 < s1Len) {
        col[(i+22)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 9);
    }
    if (i+23 == s1Len-1 && 0 <= j-23 && j-23 < s2Len) {
        row[j-23] = (int8_t)_mm256_extract_epi8_rpl(vWH, 8);
    }
    if (j-23 == s2Len-1 && 0 <= i+23 && i+23 < s1Len) {
        col[(i+23)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 8);
    }
    if (i+24 == s1Len-1 && 0 <= j-24 && j-24 < s2Len) {
        row[j-24] = (int8_t)_mm256_extract_epi8_rpl(vWH, 7);
    }
    if (j-24 == s2Len-1 && 0 <= i+24 && i+24 < s1Len) {
        col[(i+24)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 7);
    }
    if (i+25 == s1Len-1 && 0 <= j-25 && j-25 < s2Len) {
        row[j-25] = (int8_t)_mm256_extract_epi8_rpl(vWH, 6);
    }
    if (j-25 == s2Len-1 && 0 <= i+25 && i+25 < s1Len) {
        col[(i+25)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 6);
    }
    if (i+26 == s1Len-1 && 0 <= j-26 && j-26 < s2Len) {
        row[j-26] = (int8_t)_mm256_extract_epi8_rpl(vWH, 5);
    }
    if (j-26 == s2Len-1 && 0 <= i+26 && i+26 < s1Len) {
        col[(i+26)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 5);
    }
    if (i+27 == s1Len-1 && 0 <= j-27 && j-27 < s2Len) {
        row[j-27] = (int8_t)_mm256_extract_epi8_rpl(vWH, 4);
    }
    if (j-27 == s2Len-1 && 0 <= i+27 && i+27 < s1Len) {
        col[(i+27)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 4);
    }
    if (i+28 == s1Len-1 && 0 <= j-28 && j-28 < s2Len) {
        row[j-28] = (int8_t)_mm256_extract_epi8_rpl(vWH, 3);
    }
    if (j-28 == s2Len-1 && 0 <= i+28 && i+28 < s1Len) {
        col[(i+28)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 3);
    }
    if (i+29 == s1Len-1 && 0 <= j-29 && j-29 < s2Len) {
        row[j-29] = (int8_t)_mm256_extract_epi8_rpl(vWH, 2);
    }
    if (j-29 == s2Len-1 && 0 <= i+29 && i+29 < s1Len) {
        col[(i+29)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 2);
    }
    if (i+30 == s1Len-1 && 0 <= j-30 && j-30 < s2Len) {
        row[j-30] = (int8_t)_mm256_extract_epi8_rpl(vWH, 1);
    }
    if (j-30 == s2Len-1 && 0 <= i+30 && i+30 < s1Len) {
        col[(i+30)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 1);
    }
    if (i+31 == s1Len-1 && 0 <= j-31 && j-31 < s2Len) {
        row[j-31] = (int8_t)_mm256_extract_epi8_rpl(vWH, 0);
    }
    if (j-31 == s2Len-1 && 0 <= i+31 && i+31 < s1Len) {
        col[(i+31)] = (int8_t)_mm256_extract_epi8_rpl(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_table_diag_avx2_256_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_rowcol_diag_avx2_256_8
#else
#define FNAME parasail_nw_diag_avx2_256_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 32; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int8_t * const restrict s1 = parasail_memalign_int8_t(32, s1Len+PAD);
    int8_t * const restrict s2B= parasail_memalign_int8_t(32, s2Len+PAD2);
    int8_t * const restrict _H_pr = parasail_memalign_int8_t(32, s2Len+PAD2);
    int8_t * const restrict _F_pr = parasail_memalign_int8_t(32, s2Len+PAD2);
    int8_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int8_t * const restrict H_pr = _H_pr+PAD;
    int8_t * const restrict F_pr = _F_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = s1Len-1;
    int32_t end_ref = s2Len-1;
    int8_t score = NEG_INF;
    __m256i vNegInf = _mm256_set1_epi8(NEG_INF);
    __m256i vOpen = _mm256_set1_epi8(open);
    __m256i vGap  = _mm256_set1_epi8(gap);
    __m256i vOne = _mm256_set1_epi8(1);
    __m256i vN = _mm256_set1_epi8(N);
    __m256i vGapN = _mm256_set1_epi8(gap*N);
    __m256i vNegOne = _mm256_set1_epi8(-1);
    __m256i vI = _mm256_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31);
    __m256i vJreset = _mm256_set_epi8(0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,-27,-28,-29,-30,-31);
    __m256i vMax = vNegInf;
    __m256i vILimit = _mm256_set1_epi8(s1Len);
    __m256i vILimit1 = _mm256_subs_epi8(vILimit, vOne);
    __m256i vJLimit = _mm256_set1_epi8(s2Len);
    __m256i vJLimit1 = _mm256_subs_epi8(vJLimit, vOne);
    __m256i vIBoundary = _mm256_set_epi8(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap,
            -open-4*gap,
            -open-5*gap,
            -open-6*gap,
            -open-7*gap,
            -open-8*gap,
            -open-9*gap,
            -open-10*gap,
            -open-11*gap,
            -open-12*gap,
            -open-13*gap,
            -open-14*gap,
            -open-15*gap,
            -open-16*gap,
            -open-17*gap,
            -open-18*gap,
            -open-19*gap,
            -open-20*gap,
            -open-21*gap,
            -open-22*gap,
            -open-23*gap,
            -open-24*gap,
            -open-25*gap,
            -open-26*gap,
            -open-27*gap,
            -open-28*gap,
            -open-29*gap,
            -open-30*gap,
            -open-31*gap
            );
    __m256i vNegLimit = _mm256_set1_epi8(INT8_MIN);
    __m256i vPosLimit = _mm256_set1_epi8(INT8_MAX);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;

    /* convert _s1 from char to int in range 0-23 */
    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    /* pad back of s1 with dummy values */
    for (i=s1Len; i<s1Len_PAD; ++i) {
        s1[i] = 0; /* point to first matrix row because we don't care */
    }

    /* convert _s2 from char to int in range 0-23 */
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }
    /* pad front of s2 with dummy values */
    for (j=-PAD; j<0; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }
    /* pad back of s2 with dummy values */
    for (j=s2Len; j<s2Len_PAD; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }

    /* set initial values for stored row */
    for (j=0; j<s2Len; ++j) {
        H_pr[j] = -open - j*gap;
        F_pr[j] = NEG_INF;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        H_pr[j] = NEG_INF;
        F_pr[j] = NEG_INF;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = NEG_INF;
        F_pr[j] = NEG_INF;
    }
    H_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m256i vNH = vNegInf;
        __m256i vWH = vNegInf;
        __m256i vE = vNegInf;
        __m256i vF = vNegInf;
        __m256i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        const int * const restrict matrow4 = &matrix->matrix[matrix->size*s1[i+4]];
        const int * const restrict matrow5 = &matrix->matrix[matrix->size*s1[i+5]];
        const int * const restrict matrow6 = &matrix->matrix[matrix->size*s1[i+6]];
        const int * const restrict matrow7 = &matrix->matrix[matrix->size*s1[i+7]];
        const int * const restrict matrow8 = &matrix->matrix[matrix->size*s1[i+8]];
        const int * const restrict matrow9 = &matrix->matrix[matrix->size*s1[i+9]];
        const int * const restrict matrow10 = &matrix->matrix[matrix->size*s1[i+10]];
        const int * const restrict matrow11 = &matrix->matrix[matrix->size*s1[i+11]];
        const int * const restrict matrow12 = &matrix->matrix[matrix->size*s1[i+12]];
        const int * const restrict matrow13 = &matrix->matrix[matrix->size*s1[i+13]];
        const int * const restrict matrow14 = &matrix->matrix[matrix->size*s1[i+14]];
        const int * const restrict matrow15 = &matrix->matrix[matrix->size*s1[i+15]];
        const int * const restrict matrow16 = &matrix->matrix[matrix->size*s1[i+16]];
        const int * const restrict matrow17 = &matrix->matrix[matrix->size*s1[i+17]];
        const int * const restrict matrow18 = &matrix->matrix[matrix->size*s1[i+18]];
        const int * const restrict matrow19 = &matrix->matrix[matrix->size*s1[i+19]];
        const int * const restrict matrow20 = &matrix->matrix[matrix->size*s1[i+20]];
        const int * const restrict matrow21 = &matrix->matrix[matrix->size*s1[i+21]];
        const int * const restrict matrow22 = &matrix->matrix[matrix->size*s1[i+22]];
        const int * const restrict matrow23 = &matrix->matrix[matrix->size*s1[i+23]];
        const int * const restrict matrow24 = &matrix->matrix[matrix->size*s1[i+24]];
        const int * const restrict matrow25 = &matrix->matrix[matrix->size*s1[i+25]];
        const int * const restrict matrow26 = &matrix->matrix[matrix->size*s1[i+26]];
        const int * const restrict matrow27 = &matrix->matrix[matrix->size*s1[i+27]];
        const int * const restrict matrow28 = &matrix->matrix[matrix->size*s1[i+28]];
        const int * const restrict matrow29 = &matrix->matrix[matrix->size*s1[i+29]];
        const int * const restrict matrow30 = &matrix->matrix[matrix->size*s1[i+30]];
        const int * const restrict matrow31 = &matrix->matrix[matrix->size*s1[i+31]];
        vNH = _mm256_srli_si256_rpl(vNH, 1);
        vNH = _mm256_insert_epi8_rpl(vNH, H_pr[-1], 31);
        vWH = _mm256_srli_si256_rpl(vWH, 1);
        vWH = _mm256_insert_epi8_rpl(vWH, -open - i*gap, 31);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m256i vMat;
            __m256i vNWH = vNH;
            vNH = _mm256_srli_si256_rpl(vWH, 1);
            vNH = _mm256_insert_epi8_rpl(vNH, H_pr[j], 31);
            vF = _mm256_srli_si256_rpl(vF, 1);
            vF = _mm256_insert_epi8_rpl(vF, F_pr[j], 31);
            vF = _mm256_max_epi8(
                    _mm256_subs_epi8(vNH, vOpen),
                    _mm256_subs_epi8(vF, vGap));
            vE = _mm256_max_epi8(
                    _mm256_subs_epi8(vWH, vOpen),
                    _mm256_subs_epi8(vE, vGap));
            vMat = _mm256_set_epi8(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]],
                    matrow8[s2[j-8]],
                    matrow9[s2[j-9]],
                    matrow10[s2[j-10]],
                    matrow11[s2[j-11]],
                    matrow12[s2[j-12]],
                    matrow13[s2[j-13]],
                    matrow14[s2[j-14]],
                    matrow15[s2[j-15]],
                    matrow16[s2[j-16]],
                    matrow17[s2[j-17]],
                    matrow18[s2[j-18]],
                    matrow19[s2[j-19]],
                    matrow20[s2[j-20]],
                    matrow21[s2[j-21]],
                    matrow22[s2[j-22]],
                    matrow23[s2[j-23]],
                    matrow24[s2[j-24]],
                    matrow25[s2[j-25]],
                    matrow26[s2[j-26]],
                    matrow27[s2[j-27]],
                    matrow28[s2[j-28]],
                    matrow29[s2[j-29]],
                    matrow30[s2[j-30]],
                    matrow31[s2[j-31]]
                    );
            vNWH = _mm256_adds_epi8(vNWH, vMat);
            vWH = _mm256_max_epi8(vNWH, vE);
            vWH = _mm256_max_epi8(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m256i cond = _mm256_cmpeq_epi8(vJ,vNegOne);
                vWH = _mm256_blendv_epi8(vWH, vIBoundary, cond);
                vF = _mm256_blendv_epi8(vF, vNegInf, cond);
                vE = _mm256_blendv_epi8(vE, vNegInf, cond);
            }
            /* check for saturation */
            {
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vWH);
                vSaturationCheckMin = _mm256_min_epi8(vSaturationCheckMin, vWH);
            }
#ifdef PARASAIL_TABLE
            arr_store_si256(result->tables->score_table, vWH, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->rowcols->score_row, result->rowcols->score_col, vWH, i, s1Len, j, s2Len);
#endif
            H_pr[j-31] = (int8_t)_mm256_extract_epi8_rpl(vWH,0);
            F_pr[j-31] = (int8_t)_mm256_extract_epi8_rpl(vF,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                __m256i cond_valid_I = _mm256_cmpeq_epi8(vI, vILimit1);
                __m256i cond_valid_J = _mm256_cmpeq_epi8(vJ, vJLimit1);
                __m256i cond_all = _mm256_and_si256(cond_valid_I, cond_valid_J);
                vMax = _mm256_blendv_epi8(vMax, vWH, cond_all);
            }
            vJ = _mm256_adds_epi8(vJ, vOne);
        }
        vI = _mm256_adds_epi8(vI, vN);
        vIBoundary = _mm256_subs_epi8(vIBoundary, vGapN);
    }

    /* max in vMax */
    for (i=0; i<N; ++i) {
        int8_t value;
        value = (int8_t) _mm256_extract_epi8_rpl(vMax, 31);
        if (value > score) {
            score = value;
        }
        vMax = _mm256_slli_si256_rpl(vMax, 1);
    }

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            _mm256_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_32;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}


