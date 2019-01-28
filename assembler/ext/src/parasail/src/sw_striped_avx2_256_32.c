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

#define NEG_INF (INT32_MIN/(int32_t)(2))

#if HAVE_AVX2_MM256_EXTRACT_EPI32
#define _mm256_extract_epi32_rpl _mm256_extract_epi32
#else
static inline int32_t _mm256_extract_epi32_rpl(__m256i a, int imm) {
    __m256i_32_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

static inline int32_t _mm256_hmax_epi32_rpl(__m256i a) {
    a = _mm256_max_epi32(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi32(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi32(a, _mm256_slli_si256(a, 4));
    return _mm256_extract_epi32_rpl(a, 7);
}


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*(0*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32_rpl(vH, 0);
    array[1LL*(1*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32_rpl(vH, 1);
    array[1LL*(2*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32_rpl(vH, 2);
    array[1LL*(3*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32_rpl(vH, 3);
    array[1LL*(4*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32_rpl(vH, 4);
    array[1LL*(5*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32_rpl(vH, 5);
    array[1LL*(6*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32_rpl(vH, 6);
    array[1LL*(7*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32_rpl(vH, 7);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int32_t)_mm256_extract_epi32_rpl(vH, 0);
    col[1*seglen+t] = (int32_t)_mm256_extract_epi32_rpl(vH, 1);
    col[2*seglen+t] = (int32_t)_mm256_extract_epi32_rpl(vH, 2);
    col[3*seglen+t] = (int32_t)_mm256_extract_epi32_rpl(vH, 3);
    col[4*seglen+t] = (int32_t)_mm256_extract_epi32_rpl(vH, 4);
    col[5*seglen+t] = (int32_t)_mm256_extract_epi32_rpl(vH, 5);
    col[6*seglen+t] = (int32_t)_mm256_extract_epi32_rpl(vH, 6);
    col[7*seglen+t] = (int32_t)_mm256_extract_epi32_rpl(vH, 7);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_table_striped_avx2_256_32
#define PNAME parasail_sw_table_striped_profile_avx2_256_32
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_rowcol_striped_avx2_256_32
#define PNAME parasail_sw_rowcol_striped_profile_avx2_256_32
#else
#define FNAME parasail_sw_striped_avx2_256_32
#define PNAME parasail_sw_striped_profile_avx2_256_32
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_avx_256_32(s1, s1Len, matrix);
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
    const int32_t segWidth = 8; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m256i* const restrict vProfile = (__m256i*)profile->profile32.score;
    __m256i* restrict pvHStore = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLoad =  parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMax = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvE = parasail_memalign___m256i(32, segLen);
    __m256i vGapO = _mm256_set1_epi32(open);
    __m256i vGapE = _mm256_set1_epi32(gap);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vNegInf = _mm256_set1_epi32(NEG_INF);
    int32_t score = NEG_INF;
    __m256i vMaxH = vNegInf;
    __m256i vMaxHUnit = vNegInf;
    int32_t maxp = INT32_MAX - (int32_t)(matrix->max+1);
    /*int32_t stop = profile->stop == INT32_MAX ?  INT32_MAX : (int32_t)profile->stop;*/
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
    parasail_memset___m256i(pvHStore, vZero, segLen);
    parasail_memset___m256i(pvE, _mm256_set1_epi32(-open), segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        __m256i vF;
        __m256i vH;
        const __m256i* vP = NULL;
        __m256i* pv = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = vZero;

        /* load final segment of pvHStore and shift left by 4 bytes */
        vH = _mm256_load_si256(&pvHStore[segLen - 1]);
        vH = _mm256_slli_si256_rpl(vH, 4);

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
            vH = _mm256_add_epi32(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi32(vH, vE);
            vH = _mm256_max_epi32(vH, vF);
            vH = _mm256_max_epi32(vH, vZero);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif
            vMaxH = _mm256_max_epi32(vH, vMaxH);

            /* Update vE value. */
            vH = _mm256_sub_epi32(vH, vGapO);
            vE = _mm256_sub_epi32(vE, vGapE);
            vE = _mm256_max_epi32(vE, vH);
            _mm256_store_si256(pvE + i, vE);

            /* Update vF value. */
            vF = _mm256_sub_epi32(vF, vGapE);
            vF = _mm256_max_epi32(vF, vH);

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            vF = _mm256_slli_si256_rpl(vF, 4);
            for (i=0; i<segLen; ++i) {
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi32(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                vMaxH = _mm256_max_epi32(vH, vMaxH);
                vH = _mm256_sub_epi32(vH, vGapO);
                vF = _mm256_sub_epi32(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi32(vF, vH))) goto end;
                /*vF = _mm256_max_epi32(vF, vH);*/
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
                vH = _mm256_slli_si256_rpl(vH, 4);
            }
            result->rowcols->score_row[j] = (int32_t) _mm256_extract_epi32_rpl (vH, 7);
        }
#endif

        {
            __m256i vCompare = _mm256_cmpgt_epi32(vMaxH, vMaxHUnit);
            if (_mm256_movemask_epi8(vCompare)) {
                score = _mm256_hmax_epi32_rpl(vMaxH);
                /* if score has potential to overflow, abort early */
                if (score > maxp) {
                    result->flag |= PARASAIL_FLAG_SATURATED;
                    break;
                }
                vMaxHUnit = _mm256_set1_epi32(score);
                end_ref = j;
            }
        }

        /*if (score == stop) break;*/
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m256i vH = _mm256_load_si256(pvHStore+i);
        arr_store_col(result->rowcols->score_col, vH, i, segLen);
    }
#endif

    if (score == INT32_MAX) {
        result->flag |= PARASAIL_FLAG_SATURATED;
    }

    if (parasail_result_is_saturated(result)) {
        score = INT32_MAX;
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
            int32_t *t = (int32_t*)pvHMax;
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

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_BITS_32 | PARASAIL_FLAG_LANES_8;
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


