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
#include <string.h>

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <emmintrin.h>
#include <smmintrin.h>
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_sse.h"


static inline int8_t _mm_hmax_epi8_rpl(__m128i a) {
    a = _mm_max_epi8(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi8(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi8(a, _mm_srli_si128(a, 2));
    a = _mm_max_epi8(a, _mm_srli_si128(a, 1));
    return _mm_extract_epi8(a, 0);
}


static inline void arr_store(
        __m128i *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    _mm_store_si128(array + (1LL*d*seglen+t), vH);
}

static inline __m128i arr_load(
        __m128i *array,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    return _mm_load_si128(array + (1LL*d*seglen+t));
}

#define FNAME parasail_sw_trace_scan_sse41_128_8
#define PNAME parasail_sw_trace_scan_profile_sse41_128_8

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_sse_128_8(s1, s1Len, matrix);
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
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict pvP  = (__m128i*)profile->profile8.score;
    __m128i* const restrict pvE  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHt = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvH  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHMax  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvGapper = parasail_memalign___m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    const int8_t NEG_LIMIT = (-open < matrix->min ?
        INT8_MIN + open : INT8_MIN - matrix->min) + 1;
    const int8_t POS_LIMIT = INT8_MAX - matrix->max - 1;
    __m128i vZero = _mm_setzero_si128();
    int8_t score = NEG_LIMIT;
    __m128i vNegLimit = _mm_set1_epi8(NEG_LIMIT);
    __m128i vPosLimit = _mm_set1_epi8(POS_LIMIT);
    __m128i vSaturationCheckMin = vPosLimit;
    __m128i vSaturationCheckMax = vNegLimit;
    __m128i vMaxH = vNegLimit;
    __m128i vMaxHUnit = vNegLimit;
    __m128i vNegInfFront = vZero;
    __m128i vSegLenXgap;
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(__m128i));
    __m128i vTZero = _mm_set1_epi8(PARASAIL_ZERO);
    __m128i vTIns  = _mm_set1_epi8(PARASAIL_INS);
    __m128i vTDel  = _mm_set1_epi8(PARASAIL_DEL);
    __m128i vTDiag = _mm_set1_epi8(PARASAIL_DIAG);
    __m128i vTDiagE = _mm_set1_epi8(PARASAIL_DIAG_E);
    __m128i vTInsE = _mm_set1_epi8(PARASAIL_INS_E);
    __m128i vTDiagF = _mm_set1_epi8(PARASAIL_DIAG_F);
    __m128i vTDelF = _mm_set1_epi8(PARASAIL_DEL_F);

    vNegInfFront = _mm_insert_epi8(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = _mm_adds_epi8(vNegInfFront,
            _mm_slli_si128(_mm_set1_epi8(-segLen*gap), 1));

    parasail_memset___m128i(pvH, vZero, segLen);
    parasail_memset___m128i(pvE, vNegLimit, segLen);
    {
        __m128i vGapper = _mm_subs_epi8(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            _mm_store_si128(pvGapper+i, vGapper);
            vGapper = _mm_subs_epi8(vGapper, vGapE);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vE_ext;
        __m128i vE_opn;
        __m128i vHt;
        __m128i vF;
        __m128i vF_ext;
        __m128i vF_opn;
        __m128i vH;
        __m128i vHp;
        __m128i *pvW;
        __m128i vW;
        __m128i case1;
        __m128i case2;
        __m128i case0;
        __m128i vGapper;
        __m128i vT;
        __m128i vET;
        __m128i vFT;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = _mm_load_si128(pvH+(segLen-1));
        vHp = _mm_slli_si128(vHp, 1);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vZero;
        vF = vNegLimit;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vGapper = _mm_load_si128(pvGapper+i);
            vE_opn = _mm_subs_epi8(vH, vGapO);
            vE_ext = _mm_subs_epi8(vE, vGapE);
            case1 = _mm_cmpgt_epi8(vE_opn, vE_ext);
            vET = _mm_blendv_epi8(vTInsE, vTDiagE, case1);
            arr_store(result->trace->trace_table, vET, i, segLen, j);
            vE = _mm_max_epi8(vE_opn, vE_ext);
            vGapper = _mm_adds_epi8(vHt, vGapper);
            vF = _mm_max_epi8(vF, vGapper);
            vHp = _mm_adds_epi8(vHp, vW);
            vHt = _mm_max_epi8(vE, vHp);
            _mm_store_si128(pvE+i, vE);
            _mm_store_si128(pvHt+i, vHt);
            _mm_store_si128(pvH+i, vHp);
            vHp = vH;
        }

        /* pseudo prefix scan on F and H */
        vHt = _mm_slli_si128(vHt, 1);
        vGapper = _mm_load_si128(pvGapper);
        vGapper = _mm_adds_epi8(vHt, vGapper);
        vF = _mm_max_epi8(vF, vGapper);
        for (i=0; i<segWidth-2; ++i) {
            __m128i vFt = _mm_slli_si128(vF, 1);
            vFt = _mm_adds_epi8(vFt, vSegLenXgap);
            vF = _mm_max_epi8(vF, vFt);
        }

        /* calculate final H */
        vF = _mm_slli_si128(vF, 1);
        vF = _mm_adds_epi8(vF, vNegInfFront);
        vH = _mm_max_epi8(vF, vHt);
        vH = _mm_max_epi8(vH, vZero);
        for (i=0; i<segLen; ++i) {
            vET = arr_load(result->trace->trace_table, i, segLen, j);
            vHp = _mm_load_si128(pvH+i);
            vHt = _mm_load_si128(pvHt+i);
            vF_opn = _mm_subs_epi8(vH, vGapO);
            vF_ext = _mm_subs_epi8(vF, vGapE);
            vF = _mm_max_epi8(vF_opn, vF_ext);
            case1 = _mm_cmpgt_epi8(vF_opn, vF_ext);
            vFT = _mm_blendv_epi8(vTDelF, vTDiagF, case1);
            vH = _mm_max_epi8(vHt, vF);
            vH = _mm_max_epi8(vH, vZero);
            case0 = _mm_cmpeq_epi8(vH, vZero);
            case1 = _mm_cmpeq_epi8(vH, vHp);
            case2 = _mm_cmpeq_epi8(vH, vF);
            vT = _mm_blendv_epi8(
                    _mm_blendv_epi8(vTIns, vTDel, case2),
                    vTDiag, case1);
            vT = _mm_blendv_epi8(vT, vTZero, case0);
            vT = _mm_or_si128(vT, vET);
            vT = _mm_or_si128(vT, vFT);
            arr_store(result->trace->trace_table, vT, i, segLen, j);
            _mm_store_si128(pvH+i, vH);
            vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vH);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vH);
            {
                vMaxH = _mm_max_epi8(vH, vMaxH);
            }
        } 

        {
            __m128i vCompare = _mm_cmpgt_epi8(vMaxH, vMaxHUnit);
            if (_mm_movemask_epi8(vCompare)) {
                score = _mm_hmax_epi8_rpl(vMaxH);
                vMaxHUnit = _mm_set1_epi8(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvH, sizeof(__m128i)*segLen);
            }
        }
    }

    /* Trace the alignment ending position on read. */
    {
        int8_t *t = (int8_t*)pvHMax;
        int32_t column_len = segLen * segWidth;
        end_query = s1Len;
        for (i = 0; i<column_len; ++i, ++t) {
            if (*t == score) {
                int32_t temp = i / segWidth + i % segWidth * segLen;
                if (temp < end_query) {
                    end_query = temp;
                }
            }
        }
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmplt_epi8(vSaturationCheckMin, vNegLimit),
            _mm_cmpgt_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;

    parasail_free(pvGapper);
    parasail_free(pvHMax);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}


