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

#include <immintrin.h>

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_avx.h"

#define NEG_INF INT8_MIN

#if HAVE_AVX2_MM256_EXTRACT_EPI8
#define _mm256_extract_epi8_rpl _mm256_extract_epi8
#else
static inline int8_t _mm256_extract_epi8_rpl(__m256i a, int imm) {
    __m256i_8_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

static inline int8_t _mm256_hmax_epi8_rpl(__m256i a) {
    a = _mm256_max_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 4));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 2));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 1));
    return _mm256_extract_epi8_rpl(a, 31);
}


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen,
        int32_t bias)
{
    array[1LL*( 0*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  0) - bias;
    array[1LL*( 1*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  1) - bias;
    array[1LL*( 2*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  2) - bias;
    array[1LL*( 3*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  3) - bias;
    array[1LL*( 4*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  4) - bias;
    array[1LL*( 5*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  5) - bias;
    array[1LL*( 6*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  6) - bias;
    array[1LL*( 7*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  7) - bias;
    array[1LL*( 8*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  8) - bias;
    array[1LL*( 9*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  9) - bias;
    array[1LL*(10*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 10) - bias;
    array[1LL*(11*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 11) - bias;
    array[1LL*(12*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 12) - bias;
    array[1LL*(13*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 13) - bias;
    array[1LL*(14*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 14) - bias;
    array[1LL*(15*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 15) - bias;
    array[1LL*(16*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 16) - bias;
    array[1LL*(17*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 17) - bias;
    array[1LL*(18*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 18) - bias;
    array[1LL*(19*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 19) - bias;
    array[1LL*(20*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 20) - bias;
    array[1LL*(21*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 21) - bias;
    array[1LL*(22*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 22) - bias;
    array[1LL*(23*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 23) - bias;
    array[1LL*(24*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 24) - bias;
    array[1LL*(25*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 25) - bias;
    array[1LL*(26*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 26) - bias;
    array[1LL*(27*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 27) - bias;
    array[1LL*(28*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 28) - bias;
    array[1LL*(29*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 29) - bias;
    array[1LL*(30*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 30) - bias;
    array[1LL*(31*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 31) - bias;
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t bias)
{
    col[ 0*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  0) - bias;
    col[ 1*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  1) - bias;
    col[ 2*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  2) - bias;
    col[ 3*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  3) - bias;
    col[ 4*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  4) - bias;
    col[ 5*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  5) - bias;
    col[ 6*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  6) - bias;
    col[ 7*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  7) - bias;
    col[ 8*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  8) - bias;
    col[ 9*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  9) - bias;
    col[10*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 10) - bias;
    col[11*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 11) - bias;
    col[12*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 12) - bias;
    col[13*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 13) - bias;
    col[14*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 14) - bias;
    col[15*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 15) - bias;
    col[16*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 16) - bias;
    col[17*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 17) - bias;
    col[18*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 18) - bias;
    col[19*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 19) - bias;
    col[20*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 20) - bias;
    col[21*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 21) - bias;
    col[22*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 22) - bias;
    col[23*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 23) - bias;
    col[24*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 24) - bias;
    col[25*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 25) - bias;
    col[26*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 26) - bias;
    col[27*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 27) - bias;
    col[28*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 28) - bias;
    col[29*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 29) - bias;
    col[30*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 30) - bias;
    col[31*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 31) - bias;
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_table_striped_avx2_256_8
#define PNAME parasail_sw_table_striped_profile_avx2_256_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_rowcol_striped_avx2_256_8
#define PNAME parasail_sw_rowcol_striped_profile_avx2_256_8
#else
#define FNAME parasail_sw_striped_avx2_256_8
#define PNAME parasail_sw_striped_profile_avx2_256_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_avx_256_8(s1, s1Len, matrix);
    parasail_result_t *result = PNAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 32; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m256i* const restrict vProfile = (__m256i*)profile->profile8.score;
    __m256i* restrict pvHStore = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLoad = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMax = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvE = parasail_memalign___m256i(32, segLen);
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    __m256i vZero = _mm256_set1_epi8(0);
    int8_t bias = INT8_MIN;
    int8_t score = bias;
    __m256i vBias = _mm256_set1_epi8(bias);
    __m256i vMaxH = vBias;
    __m256i vMaxHUnit = vBias;
    int8_t maxp = INT8_MAX - (int8_t)(matrix->max+1);
    __m256i insert_mask = _mm256_cmpgt_epi8(
            _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
            vZero);
    /*int8_t stop = profile->stop == INT32_MAX ?  INT8_MAX : (int8_t)profile->stop-bias;*/
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    /* initialize H and E */
    parasail_memset___m256i(pvHStore, vBias, segLen);
    parasail_memset___m256i(pvE, vBias, segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        __m256i vF;
        __m256i vH;
        const __m256i* vP = NULL;
        __m256i* pv = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vF = vBias;

        /* load final segment of pvHStore and shift left by 1 bytes */
        vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 1);
        vH = _mm256_blendv_epi8(vH, vBias, insert_mask);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        if (end_ref == j-2) {
            /* Swap in the max buffer. */
            pv = pvHMax;
            pvHMax = pvHLoad;
            pvHLoad = pvHStore;
            pvHStore = pv;
        }
        else {
            /* Swap the 2 H buffers. */
            pv = pvHLoad;
            pvHLoad = pvHStore;
            pvHStore = pv;
        }

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = _mm256_adds_epi8(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi8(vH, vE);
            vH = _mm256_max_epi8(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
            vMaxH = _mm256_max_epi8(vH, vMaxH);

            /* Update vE value. */
            vH = _mm256_subs_epi8(vH, vGapO);
            vE = _mm256_subs_epi8(vE, vGapE);
            vE = _mm256_max_epi8(vE, vH);
            _mm256_store_si256(pvE + i, vE);

            /* Update vF value. */
            vF = _mm256_subs_epi8(vF, vGapE);
            vF = _mm256_max_epi8(vF, vH);

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            vF = _mm256_slli_si256_rpl(vF, 1);
            vF = _mm256_blendv_epi8(vF, vBias, insert_mask);
            for (i=0; i<segLen; ++i) {
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi8(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
                vMaxH = _mm256_max_epi8(vH, vMaxH);
                vH = _mm256_subs_epi8(vH, vGapO);
                vF = _mm256_subs_epi8(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi8(vF, vH))) goto end;
                /*vF = _mm256_max_epi8(vF, vH);*/
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm256_load_si256(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 1);
            }
            result->rowcols->score_row[j] = (int8_t) _mm256_extract_epi8_rpl (vH, 31) - bias;
        }
#endif

        {
            __m256i vCompare = _mm256_cmpgt_epi8(vMaxH, vMaxHUnit);
            if (_mm256_movemask_epi8(vCompare)) {
                score = _mm256_hmax_epi8_rpl(vMaxH);
                /* if score has potential to overflow, abort early */
                if (score > maxp) {
                    result->flag |= PARASAIL_FLAG_SATURATED;
                    break;
                }
                vMaxHUnit = _mm256_set1_epi8(score);
                end_ref = j;
            }
        }

        /*if (score == stop) break;*/
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m256i vH = _mm256_load_si256(pvHStore+i);
        arr_store_col(result->rowcols->score_col, vH, i, segLen, bias);
    }
#endif

    if (score == INT8_MAX) {
        result->flag |= PARASAIL_FLAG_SATURATED;
    }

    if (parasail_result_is_saturated(result)) {
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }
    else {
        if (end_ref == j-1) {
            /* end_ref was the last store column */
            __m256i *pv = pvHMax;
            pvHMax = pvHStore;
            pvHStore = pv;
        }
        else if (end_ref == j-2) {
            /* end_ref was the last load column */
            __m256i *pv = pvHMax;
            pvHMax = pvHLoad;
            pvHLoad = pv;
        }
        /* Trace the alignment ending position on read. */
        {
            int8_t *t = (int8_t*)pvHMax;
            int32_t column_len = segLen * segWidth;
            end_query = s1Len - 1;
            for (i = 0; i<column_len; ++i, ++t) {
                if (*t == score) {
                    int32_t temp = i / segWidth + i % segWidth * segLen;
                    if (temp < end_query) {
                        end_query = temp;
                    }
                }
            }
        }
    }

    result->score = score - bias;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_32;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(pvE);
    parasail_free(pvHMax);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

