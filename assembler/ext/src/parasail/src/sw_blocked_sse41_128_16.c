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

#include <emmintrin.h>
#include <smmintrin.h>

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_sse.h"

#define NEG_INF_16 (INT16_MIN/(int16_t)(2))

#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t i,
        int32_t segWidth,
        int32_t j,
        int32_t s2Len)
{
    array[(i*segWidth+0)*s2Len + j] = (int16_t)_mm_extract_epi16(vH, 0);
    array[(i*segWidth+1)*s2Len + j] = (int16_t)_mm_extract_epi16(vH, 1);
    array[(i*segWidth+2)*s2Len + j] = (int16_t)_mm_extract_epi16(vH, 2);
    array[(i*segWidth+3)*s2Len + j] = (int16_t)_mm_extract_epi16(vH, 3);
    array[(i*segWidth+4)*s2Len + j] = (int16_t)_mm_extract_epi16(vH, 4);
    array[(i*segWidth+5)*s2Len + j] = (int16_t)_mm_extract_epi16(vH, 5);
    array[(i*segWidth+6)*s2Len + j] = (int16_t)_mm_extract_epi16(vH, 6);
    array[(i*segWidth+7)*s2Len + j] = (int16_t)_mm_extract_epi16(vH, 7);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m128i vH,
        int32_t i,
        int32_t segWidth)
{
    col[i*segWidth+0] = (int16_t)_mm_extract_epi16(vH, 0);
    col[i*segWidth+1] = (int16_t)_mm_extract_epi16(vH, 1);
    col[i*segWidth+2] = (int16_t)_mm_extract_epi16(vH, 2);
    col[i*segWidth+3] = (int16_t)_mm_extract_epi16(vH, 3);
    col[i*segWidth+4] = (int16_t)_mm_extract_epi16(vH, 4);
    col[i*segWidth+5] = (int16_t)_mm_extract_epi16(vH, 5);
    col[i*segWidth+6] = (int16_t)_mm_extract_epi16(vH, 6);
    col[i*segWidth+7] = (int16_t)_mm_extract_epi16(vH, 7);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_table_blocked_sse41_128_16
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_rowcol_blocked_sse41_128_16
#else
#define FNAME parasail_sw_blocked_sse41_128_16
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 8; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict vProfile = parasail_memalign___m128i(16, n * segLen);
    __m128i* restrict pvH = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvE = parasail_memalign___m128i(16, segLen);
    int score = NEG_INF_16;
    __m128i vGapO = _mm_set1_epi16(open);
    __m128i vGapE = _mm_set1_epi16(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vNegInf = _mm_set1_epi16(NEG_INF_16);
    __m128i vMaxH = vZero;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
    const int32_t offset = segLen - 1;
    const int32_t position = s1Len % segWidth;
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m128i_16_t t;
                j = i*segWidth;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix->matrix[matrix->size*matrix->mapper[(unsigned char)s1[j]]];
                    j += 1;
                }
                _mm_store_si128(&vProfile[index], t.m);
                ++index;
            }
        }
    }

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m128i_16_t h;
            __m128i_16_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = -open;
            }
            _mm_store_si128(&pvH[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vX = vZero;
        __m128i vF = vNegInf;
        const __m128i* pvP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;
        for (i=0; i<segLen; ++i) {
            __m128i vP;
            __m128i vH;
            __m128i vE;
            __m128i vT1;

            vH = _mm_load_si128(pvH + i);
            vE = _mm_load_si128(pvE + i);

            vT1 = _mm_srli_si128(vH, 14); /* rshift 3 */
            vH = _mm_slli_si128(vH, 2); /* lshift 1 */
            vH = _mm_or_si128(vH, vX);
            vX = vT1;

            vP = _mm_load_si128(pvP + i);
            vH = _mm_add_epi16(vH, vP);
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vZero);

            vF = _mm_srli_si128(vF, 14);
            vF = _mm_or_si128(vF, _mm_slli_si128(vH, 2));
            vF = _mm_sub_epi16(vF, vGapO);
            if (_mm_movemask_epi8(_mm_cmpgt_epi16(vF, vZero))) {
                __m128i vT2 = vF;
                while (_mm_movemask_epi8(_mm_cmpgt_epi16(vT2, vZero))) {
                    vT2 = _mm_slli_si128(vT2, 2);
                    vT2 = _mm_sub_epi16(vT2, vGapE);
                    vF = _mm_max_epi16(vF, vT2);
                }
                vH = _mm_max_epi16(vH, vF);
                vF = _mm_add_epi16(vF, vGapO);
                vF = _mm_sub_epi16(vF, vGapE);
                vF = _mm_max_epi16(vH, vF);
            }
            else {
                vF = vH;
            }

            _mm_store_si128(pvH + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->tables->score_table, vH, i, segWidth, j, s2Len);
#endif
            vMaxH = _mm_max_epi16(vMaxH, vH);

            vH = _mm_sub_epi16(vH, vGapO);
            vE = _mm_sub_epi16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            _mm_store_si128(pvE + i, vE);

        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            __m128i vH = _mm_load_si128(pvH + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 2);
            }
            result->rowcols->score_row[j] = (int16_t) _mm_extract_epi16 (vH, 7);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m128i vH = _mm_load_si128(pvH + i);
        arr_store_col(result->rowcols->score_col, vH, i, segWidth);
    }
#endif

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int16_t value = (int16_t) _mm_extract_epi16(vMaxH, 7);
        if (value > score) {
            score = value;
        }
        vMaxH = _mm_slli_si128(vMaxH, 2);
    }

    result->score = score;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_BLOCKED
        | PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_8;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(pvE);
    parasail_free(pvH);
    parasail_free(vProfile);

    return result;
}

