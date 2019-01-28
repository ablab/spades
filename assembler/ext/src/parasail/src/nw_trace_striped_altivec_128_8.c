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



#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_altivec.h"

#define SWAP(A,B) { vec128i* tmp = A; A = B; B = tmp; }

#define NEG_INF INT8_MIN


static inline void arr_store(
        vec128i *array,
        vec128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    _mm_store_si128(array + (1LL*d*seglen+t), vH);
}

static inline vec128i arr_load(
        vec128i *array,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    return _mm_load_si128(array + (1LL*d*seglen+t));
}

#define FNAME parasail_nw_trace_striped_altivec_128_8
#define PNAME parasail_nw_trace_striped_profile_altivec_128_8

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_altivec_128_8(s1, s1Len, matrix);
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
    const int s1Len = profile->s1Len;
    int32_t end_query = s1Len-1;
    int32_t end_ref = s2Len-1;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    vec128i* const restrict vProfile = (vec128i*)profile->profile8.score;
    vec128i* restrict pvHStore = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvHLoad =  parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvE = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvEaStore = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvEaLoad = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvHT = parasail_memalign_vec128i(16, segLen);
    int8_t* const restrict boundary = parasail_memalign_int8_t(16, s2Len+1);
    vec128i vGapO = _mm_set1_epi8(open);
    vec128i vGapE = _mm_set1_epi8(gap);
    vec128i vNegInf = _mm_set1_epi8(NEG_INF);
    int8_t score = NEG_INF;
    vec128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    vec128i vPosLimit = _mm_set1_epi8(INT8_MAX);
    vec128i vSaturationCheckMin = vPosLimit;
    vec128i vSaturationCheckMax = vNegLimit;
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(vec128i));
    vec128i vTIns  = _mm_set1_epi8(PARASAIL_INS);
    vec128i vTDel  = _mm_set1_epi8(PARASAIL_DEL);
    vec128i vTDiag = _mm_set1_epi8(PARASAIL_DIAG);
    vec128i vTDiagE = _mm_set1_epi8(PARASAIL_DIAG_E);
    vec128i vTInsE = _mm_set1_epi8(PARASAIL_INS_E);
    vec128i vTDiagF = _mm_set1_epi8(PARASAIL_DIAG_F);
    vec128i vTDelF = _mm_set1_epi8(PARASAIL_DEL_F);
    vec128i vTMask = _mm_set1_epi8(PARASAIL_ZERO_MASK);
    vec128i vFTMask = _mm_set1_epi8(PARASAIL_F_MASK);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            vec128i_8_t h;
            vec128i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
            }
            _mm_store_si128(&pvHStore[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            _mm_store_si128(&pvEaStore[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT8_MIN ? INT8_MIN : tmp;
        }
    }

    for (i=0; i<segLen; ++i) {
        arr_store(result->trace->trace_table, vTDiagE, i, segLen, 0);
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        vec128i vEF_opn;
        vec128i vE;
        vec128i vE_ext;
        vec128i vF;
        vec128i vF_ext;
        vec128i vFa;
        vec128i vFa_ext;
        vec128i vH;
        vec128i vH_dag;
        const vec128i* vP = NULL;

        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vF = vNegInf;

        /* load final segment of pvHStore and shift left by 1 bytes */
        vH = _mm_load_si128(&pvHStore[segLen - 1]);
        vH = _mm_slli_si128(vH, 1);

        /* insert upper boundary condition */
        vH = _mm_insert_epi8(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad, pvHStore)
        SWAP(pvEaLoad, pvEaStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = _mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH_dag = _mm_adds_epi8(vH, _mm_load_si128(vP + i));
            vH = _mm_max_epi8(vH_dag, vE);
            vH = _mm_max_epi8(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
            /* check for saturation */
            {
                vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vH);
                vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vE);
                vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vF);
            }

            {
                vec128i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                vec128i case1 = _mm_cmpeq_epi8(vH, vH_dag);
                vec128i case2 = _mm_cmpeq_epi8(vH, vF);
                vec128i vT = _mm_blendv_epi8(
                        _mm_blendv_epi8(vTIns, vTDel, case2),
                        vTDiag, case1);
                _mm_store_si128(pvHT + i, vT);
                vT = _mm_or_si128(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i, segLen, j);
            }

            vEF_opn = _mm_subs_epi8(vH, vGapO);

            /* Update vE value. */
            vE_ext = _mm_subs_epi8(vE, vGapE);
            vE = _mm_max_epi8(vEF_opn, vE_ext);
            _mm_store_si128(pvE + i, vE);
            {
                vec128i vEa = _mm_load_si128(pvEaLoad + i);
                vec128i vEa_ext = _mm_subs_epi8(vEa, vGapE);
                vEa = _mm_max_epi8(vEF_opn, vEa_ext);
                _mm_store_si128(pvEaStore + i, vEa);
                if (j+1<s2Len) {
                    vec128i cond = _mm_cmpgt_epi8(vEF_opn, vEa_ext);
                    vec128i vT = _mm_blendv_epi8(vTInsE, vTDiagE, cond);
                    arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                }
            }

            /* Update vF value. */
            vF_ext = _mm_subs_epi8(vF, vGapE);
            vF = _mm_max_epi8(vEF_opn, vF_ext);
            if (i+1<segLen) {
                vec128i vTAll = arr_load(result->trace->trace_table, i+1, segLen, j);
                vec128i cond = _mm_cmpgt_epi8(vEF_opn, vF_ext);
                vec128i vT = _mm_blendv_epi8(vTDelF, vTDiagF, cond);
                vT = _mm_or_si128(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i+1, segLen, j);
            }

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }


        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            int64_t tmp = boundary[j+1]-open;
            int8_t tmp2 = tmp < INT8_MIN ? INT8_MIN : tmp;
            vec128i vHp = _mm_load_si128(&pvHLoad[segLen - 1]);
            vHp = _mm_slli_si128(vHp, 1);
            vHp = _mm_insert_epi8(vHp, boundary[j], 0);
            vEF_opn = _mm_slli_si128(vEF_opn, 1);
            vEF_opn = _mm_insert_epi8(vEF_opn, tmp2, 0);
            vF_ext = _mm_slli_si128(vF_ext, 1);
            vF_ext = _mm_insert_epi8(vF_ext, NEG_INF, 0);
            vF = _mm_slli_si128(vF, 1);
            vF = _mm_insert_epi8(vF, tmp2, 0);
            vFa_ext = _mm_slli_si128(vFa_ext, 1);
            vFa_ext = _mm_insert_epi8(vFa_ext, NEG_INF, 0);
            vFa = _mm_slli_si128(vFa, 1);
            vFa = _mm_insert_epi8(vFa, tmp2, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi8(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
                /* check for saturation */
            {
                vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vH);
                vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vE);
                vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vF);
            }
                {
                    vec128i vTAll;
                    vec128i vT;
                    vec128i case1;
                    vec128i case2;
                    vec128i cond;
                    vHp = _mm_adds_epi8(vHp, _mm_load_si128(vP + i));
                    case1 = _mm_cmpeq_epi8(vH, vHp);
                    case2 = _mm_cmpeq_epi8(vH, vF);
                    cond = _mm_andnot_si128(case1,case2);
                    vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    vT = _mm_load_si128(pvHT + i);
                    vT = _mm_blendv_epi8(vT, vTDel, cond);
                    _mm_store_si128(pvHT + i, vT);
                    vTAll = _mm_and_si128(vTAll, vTMask);
                    vTAll = _mm_or_si128(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                /* Update vF value. */
                {
                    vec128i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    vec128i cond = _mm_cmpgt_epi8(vEF_opn, vFa_ext);
                    vec128i vT = _mm_blendv_epi8(vTDelF, vTDiagF, cond);
                    vTAll = _mm_and_si128(vTAll, vFTMask);
                    vTAll = _mm_or_si128(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                vEF_opn = _mm_subs_epi8(vH, vGapO);
                vF_ext = _mm_subs_epi8(vF, vGapE);
                {
                    vec128i vEa = _mm_load_si128(pvEaLoad + i);
                    vec128i vEa_ext = _mm_subs_epi8(vEa, vGapE);
                    vEa = _mm_max_epi8(vEF_opn, vEa_ext);
                    _mm_store_si128(pvEaStore + i, vEa);
                    if (j+1<s2Len) {
                        vec128i cond = _mm_cmpgt_epi8(vEF_opn, vEa_ext);
                        vec128i vT = _mm_blendv_epi8(vTInsE, vTDiagE, cond);
                        arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                    }
                }
                if (! _mm_movemask_epi8(
                            _mm_or_si128(
                                _mm_cmpgt_epi8(vF_ext, vEF_opn),
                                _mm_cmpeq_epi8(vF_ext, vEF_opn))))
                    goto end;
                /*vF = _mm_max_epi8(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = _mm_subs_epi8(vFa, vGapE);
                vFa = _mm_max_epi8(vEF_opn, vFa_ext);
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
        }
    }

    /* extract last value from the last column */
    {
        vec128i vH = _mm_load_si128(pvHStore + offset);
        for (k=0; k<position; ++k) {
            vH = _mm_slli_si128 (vH, 1);
        }
        score = (int8_t) _mm_extract_epi8 (vH, 15);
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            _mm_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;

    parasail_free(boundary);
    parasail_free(pvHT);
    parasail_free(pvEaLoad);
    parasail_free(pvEaStore);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}


