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


#if HAVE_AVX2_MM256_INSERT_EPI64
#define _mm256_insert_epi64_rpl _mm256_insert_epi64
#else
static inline __m256i _mm256_insert_epi64_rpl(__m256i a, int64_t i, int imm) {
    __m256i_64_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_SET1_EPI64X
#define _mm256_set1_epi64x_rpl _mm256_set1_epi64x
#else
static inline __m256i _mm256_set1_epi64x_rpl(int64_t i) {
    __m256i_64_t A;
    A.v[0] = i;
    A.v[1] = i;
    A.v[2] = i;
    A.v[3] = i;
    return A.m;
}
#endif

static inline __m256i _mm256_max_epi64_rpl(__m256i a, __m256i b) {
    __m256i_64_t A;
    __m256i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    A.v[2] = (A.v[2]>B.v[2]) ? A.v[2] : B.v[2];
    A.v[3] = (A.v[3]>B.v[3]) ? A.v[3] : B.v[3];
    return A.m;
}

#if HAVE_AVX2_MM256_SET_EPI64X
#define _mm256_set_epi64x_rpl _mm256_set_epi64x
#else
static inline __m256i _mm256_set_epi64x_rpl(int64_t e3, int64_t e2, int64_t e1, int64_t e0) {
    __m256i_64_t A;
    A.v[0] = e0;
    A.v[1] = e1;
    A.v[2] = e2;
    A.v[3] = e3;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI64
#define _mm256_extract_epi64_rpl _mm256_extract_epi64
#else
static inline int64_t _mm256_extract_epi64_rpl(__m256i a, int imm) {
    __m256i_64_t A;
    A.m = a;
    return A.v[imm];
}
#endif

static inline __m256i _mm256_min_epi64_rpl(__m256i a, __m256i b) {
    __m256i_64_t A;
    __m256i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]<B.v[1]) ? A.v[1] : B.v[1];
    A.v[2] = (A.v[2]<B.v[2]) ? A.v[2] : B.v[2];
    A.v[3] = (A.v[3]<B.v[3]) ? A.v[3] : B.v[3];
    return A.m;
}

#define _mm256_cmplt_epi64_rpl(a,b) _mm256_cmpgt_epi64(b,a)

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

static inline int64_t _mm256_hmax_epi64_rpl(__m256i a) {
    a = _mm256_max_epi64_rpl(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi64_rpl(a, _mm256_slli_si256(a, 8));
    return _mm256_extract_epi64_rpl(a, 3);
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
    array[1LL*(0*seglen+t)*dlen + d] = (int64_t)_mm256_extract_epi64_rpl(vH, 0);
    array[1LL*(1*seglen+t)*dlen + d] = (int64_t)_mm256_extract_epi64_rpl(vH, 1);
    array[1LL*(2*seglen+t)*dlen + d] = (int64_t)_mm256_extract_epi64_rpl(vH, 2);
    array[1LL*(3*seglen+t)*dlen + d] = (int64_t)_mm256_extract_epi64_rpl(vH, 3);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int64_t)_mm256_extract_epi64_rpl(vH, 0);
    col[1*seglen+t] = (int64_t)_mm256_extract_epi64_rpl(vH, 1);
    col[2*seglen+t] = (int64_t)_mm256_extract_epi64_rpl(vH, 2);
    col[3*seglen+t] = (int64_t)_mm256_extract_epi64_rpl(vH, 3);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_table_striped_avx2_256_64
#define PNAME parasail_sg_table_striped_profile_avx2_256_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_rowcol_striped_avx2_256_64
#define PNAME parasail_sg_rowcol_striped_profile_avx2_256_64
#else
#define FNAME parasail_sg_striped_avx2_256_64
#define PNAME parasail_sg_striped_profile_avx2_256_64
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_avx_256_64(s1, s1Len, matrix);
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
    const int32_t segWidth = 4; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict vProfile = (__m256i*)profile->profile64.score;
    __m256i* restrict pvHStore = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLoad =  parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvE = parasail_memalign___m256i(32, segLen);
    const __m256i vGapO = _mm256_set1_epi64x_rpl(open);
    const __m256i vGapE = _mm256_set1_epi64x_rpl(gap);
    const int64_t NEG_LIMIT = (-open < matrix->min ?
        INT64_MIN + open : INT64_MIN - matrix->min) + 1;
    const int64_t POS_LIMIT = INT64_MAX - matrix->max - 1;
    int64_t score = NEG_LIMIT;
    __m256i vNegLimit = _mm256_set1_epi64x_rpl(NEG_LIMIT);
    __m256i vPosLimit = _mm256_set1_epi64x_rpl(POS_LIMIT);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
    __m256i vMaxH = vNegLimit;
    __m256i vPosMask = _mm256_cmpeq_epi64(_mm256_set1_epi64x_rpl(position),
            _mm256_set_epi64x_rpl(0,1,2,3));
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    /* initialize H and E */
    parasail_memset___m256i(pvHStore, _mm256_set1_epi64x_rpl(0), segLen);
    parasail_memset___m256i(pvE, _mm256_set1_epi64x_rpl(-open), segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        __m256i vF = vNegLimit;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m256i vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 8);

        /* Correct part of the vProfile */
        const __m256i* vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m256i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = _mm256_add_epi64(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi64_rpl(vH, vE);
            vH = _mm256_max_epi64_rpl(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
            vSaturationCheckMin = _mm256_min_epi64_rpl(vSaturationCheckMin, vH);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vH);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm256_sub_epi64(vH, vGapO);
            vE = _mm256_sub_epi64(vE, vGapE);
            vE = _mm256_max_epi64_rpl(vE, vH);
            _mm256_store_si256(pvE + i, vE);

            /* Update vF value. */
            vF = _mm256_sub_epi64(vF, vGapE);
            vF = _mm256_max_epi64_rpl(vF, vH);

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            vF = _mm256_slli_si256_rpl(vF, 8);
            vF = _mm256_insert_epi64_rpl(vF, -open, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi64_rpl(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                vSaturationCheckMin = _mm256_min_epi64_rpl(vSaturationCheckMin, vH);
                vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vH);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm256_sub_epi64(vH, vGapO);
                vF = _mm256_sub_epi64(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi64(vF, vH))) goto end;
                /*vF = _mm256_max_epi64_rpl(vF, vH);*/
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            __m256i vCompare;
            vH = _mm256_load_si256(pvHStore + offset);
            vCompare = _mm256_and_si256(vPosMask, _mm256_cmpgt_epi64(vH, vMaxH));
            vMaxH = _mm256_max_epi64_rpl(vH, vMaxH);
            if (_mm256_movemask_epi8(vCompare)) {
                end_ref = j;
                end_query = s1Len - 1;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 8);
            }
            result->rowcols->score_row[j] = (int64_t) _mm256_extract_epi64_rpl (vH, 3);
#endif
        }
    }

    /* max last value from all columns */
    {
        for (k=0; k<position; ++k) {
            vMaxH = _mm256_slli_si256_rpl(vMaxH, 8);
        }
        score = (int64_t) _mm256_extract_epi64_rpl(vMaxH, 3);
    }

    /* max of last column */
    {
        int64_t score_last;
        vMaxH = vNegLimit;

        for (i=0; i<segLen; ++i) {
            __m256i vH = _mm256_load_si256(pvHStore + i);
            vMaxH = _mm256_max_epi64_rpl(vH, vMaxH);
#ifdef PARASAIL_ROWCOL
            arr_store_col(result->rowcols->score_col, vH, i, segLen);
#endif
        }

        /* max in vec */
        score_last = _mm256_hmax_epi64_rpl(vMaxH);
        if (score_last > score || (score_last == score && end_ref == s2Len - 1)) {
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                int64_t *t = (int64_t*)pvHStore;
                int32_t column_len = segLen * segWidth;
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
    }

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmplt_epi64_rpl(vSaturationCheckMin, vNegLimit),
            _mm256_cmpgt_epi64(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_4;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

