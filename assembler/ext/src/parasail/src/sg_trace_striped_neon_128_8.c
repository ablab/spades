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
#include "parasail/parasail/internal_neon.h"

#define SWAP(A,B) { simde__m128i* tmp = A; A = B; B = tmp; }

#define NEG_INF INT8_MIN


static inline void arr_store(
        simde__m128i *array,
        simde__m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    simde_mm_store_si128(array + (1LL*d*seglen+t), vH);
}

static inline simde__m128i arr_load(
        simde__m128i *array,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    return simde_mm_load_si128(array + (1LL*d*seglen+t));
}

#define FNAME parasail_sg_trace_striped_neon_128_8
#define PNAME parasail_sg_trace_striped_profile_neon_128_8

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_neon_128_8(s1, s1Len, matrix);
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
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    simde__m128i* const restrict vProfile = (simde__m128i*)profile->profile8.score;
    simde__m128i* restrict pvHStore = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHLoad = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvE = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvEaStore = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvEaLoad = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHT = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i vGapO = simde_mm_set1_epi8(open);
    simde__m128i vGapE = simde_mm_set1_epi8(gap);
    simde__m128i vNegInf = simde_mm_set1_epi8(NEG_INF);
    int8_t score = NEG_INF;
    simde__m128i vMaxH = vNegInf;
    simde__m128i vPosMask = simde_mm_cmpeq_epi8(simde_mm_set1_epi8(position),
            simde_mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15));
    simde__m128i vNegLimit = simde_mm_set1_epi8(INT8_MIN);
    simde__m128i vPosLimit = simde_mm_set1_epi8(INT8_MAX);
    simde__m128i vSaturationCheckMin = vPosLimit;
    simde__m128i vSaturationCheckMax = vNegLimit;
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(simde__m128i));
    simde__m128i vTIns  = simde_mm_set1_epi8(PARASAIL_INS);
    simde__m128i vTDel  = simde_mm_set1_epi8(PARASAIL_DEL);
    simde__m128i vTDiag = simde_mm_set1_epi8(PARASAIL_DIAG);
    simde__m128i vTDiagE = simde_mm_set1_epi8(PARASAIL_DIAG_E);
    simde__m128i vTInsE = simde_mm_set1_epi8(PARASAIL_INS_E);
    simde__m128i vTDiagF = simde_mm_set1_epi8(PARASAIL_DIAG_F);
    simde__m128i vTDelF = simde_mm_set1_epi8(PARASAIL_DEL_F);
    simde__m128i vTMask = simde_mm_set1_epi8(PARASAIL_ZERO_MASK);
    simde__m128i vFTMask = simde_mm_set1_epi8(PARASAIL_F_MASK);

    /* initialize H and E */
    parasail_memset_simde__m128i(pvHStore, simde_mm_set1_epi8(0), segLen);
    parasail_memset_simde__m128i(pvE, simde_mm_set1_epi8(-open), segLen);
    parasail_memset_simde__m128i(pvEaStore, simde_mm_set1_epi8(-open), segLen);

    for (i=0; i<segLen; ++i) {
        arr_store(result->trace->trace_table, vTDiagE, i, segLen, 0);
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vEF_opn;
        simde__m128i vE;
        simde__m128i vE_ext;
        simde__m128i vF;
        simde__m128i vF_ext;
        simde__m128i vFa;
        simde__m128i vFa_ext;
        simde__m128i vH;
        simde__m128i vH_dag;
        const simde__m128i* vP = NULL;

        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = vNegInf;

        /* load final segment of pvHStore and shift left by 1 bytes */
        vH = simde_mm_load_si128(&pvHStore[segLen - 1]);
        vH = simde_mm_slli_si128(vH, 1);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad, pvHStore)
        SWAP(pvEaLoad, pvEaStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = simde_mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH_dag = simde_mm_adds_epi8(vH, simde_mm_load_si128(vP + i));
            vH = simde_mm_max_epi8(vH_dag, vE);
            vH = simde_mm_max_epi8(vH, vF);
            /* Save vH values. */
            simde_mm_store_si128(pvHStore + i, vH);
            /* check for saturation */
            {
                vSaturationCheckMax = simde_mm_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vH);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vE);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vF);
            }

            {
                simde__m128i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                simde__m128i case1 = simde_mm_cmpeq_epi8(vH, vH_dag);
                simde__m128i case2 = simde_mm_cmpeq_epi8(vH, vF);
                simde__m128i vT = simde_mm_blendv_epi8(
                        simde_mm_blendv_epi8(vTIns, vTDel, case2),
                        vTDiag, case1);
                simde_mm_store_si128(pvHT + i, vT);
                vT = simde_mm_or_si128(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i, segLen, j);
            }

            vEF_opn = simde_mm_subs_epi8(vH, vGapO);

            /* Update vE value. */
            vE_ext = simde_mm_subs_epi8(vE, vGapE);
            vE = simde_mm_max_epi8(vEF_opn, vE_ext);
            simde_mm_store_si128(pvE + i, vE);
            {
                simde__m128i vEa = simde_mm_load_si128(pvEaLoad + i);
                simde__m128i vEa_ext = simde_mm_subs_epi8(vEa, vGapE);
                vEa = simde_mm_max_epi8(vEF_opn, vEa_ext);
                simde_mm_store_si128(pvEaStore + i, vEa);
                if (j+1<s2Len) {
                    simde__m128i cond = simde_mm_cmpgt_epi8(vEF_opn, vEa_ext);
                    simde__m128i vT = simde_mm_blendv_epi8(vTInsE, vTDiagE, cond);
                    arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                }
            }

            /* Update vF value. */
            vF_ext = simde_mm_subs_epi8(vF, vGapE);
            vF = simde_mm_max_epi8(vEF_opn, vF_ext);
            if (i+1<segLen) {
                simde__m128i vTAll = arr_load(result->trace->trace_table, i+1, segLen, j);
                simde__m128i cond = simde_mm_cmpgt_epi8(vEF_opn, vF_ext);
                simde__m128i vT = simde_mm_blendv_epi8(vTDelF, vTDiagF, cond);
                vT = simde_mm_or_si128(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i+1, segLen, j);
            }

            /* Load the next vH. */
            vH = simde_mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            simde__m128i vHp = simde_mm_load_si128(&pvHLoad[segLen - 1]);
            vHp = simde_mm_slli_si128(vHp, 1);
            vEF_opn = simde_mm_slli_si128(vEF_opn, 1);
            vEF_opn = simde_mm_insert_epi8(vEF_opn, -open, 0);
            vF_ext = simde_mm_slli_si128(vF_ext, 1);
            vF_ext = simde_mm_insert_epi8(vF_ext, NEG_INF, 0);
            vF = simde_mm_slli_si128(vF, 1);
            vF = simde_mm_insert_epi8(vF, -open, 0);
            vFa_ext = simde_mm_slli_si128(vFa_ext, 1);
            vFa_ext = simde_mm_insert_epi8(vFa_ext, NEG_INF, 0);
            vFa = simde_mm_slli_si128(vFa, 1);
            vFa = simde_mm_insert_epi8(vFa, -open, 0);
            for (i=0; i<segLen; ++i) {
                vH = simde_mm_load_si128(pvHStore + i);
                vH = simde_mm_max_epi8(vH,vF);
                simde_mm_store_si128(pvHStore + i, vH);
                /* check for saturation */
            {
                vSaturationCheckMax = simde_mm_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vH);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vE);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vF);
            }
                {
                    simde__m128i vTAll;
                    simde__m128i vT;
                    simde__m128i case1;
                    simde__m128i case2;
                    simde__m128i cond;
                    vHp = simde_mm_adds_epi8(vHp, simde_mm_load_si128(vP + i));
                    case1 = simde_mm_cmpeq_epi8(vH, vHp);
                    case2 = simde_mm_cmpeq_epi8(vH, vF);
                    cond = simde_mm_andnot_si128(case1,case2);
                    vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    vT = simde_mm_load_si128(pvHT + i);
                    vT = simde_mm_blendv_epi8(vT, vTDel, cond);
                    simde_mm_store_si128(pvHT + i, vT);
                    vTAll = simde_mm_and_si128(vTAll, vTMask);
                    vTAll = simde_mm_or_si128(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                /* Update vF value. */
                {
                    simde__m128i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    simde__m128i cond = simde_mm_cmpgt_epi8(vEF_opn, vFa_ext);
                    simde__m128i vT = simde_mm_blendv_epi8(vTDelF, vTDiagF, cond);
                    vTAll = simde_mm_and_si128(vTAll, vFTMask);
                    vTAll = simde_mm_or_si128(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                vEF_opn = simde_mm_subs_epi8(vH, vGapO);
                vF_ext = simde_mm_subs_epi8(vF, vGapE);
                {
                    simde__m128i vEa = simde_mm_load_si128(pvEaLoad + i);
                    simde__m128i vEa_ext = simde_mm_subs_epi8(vEa, vGapE);
                    vEa = simde_mm_max_epi8(vEF_opn, vEa_ext);
                    simde_mm_store_si128(pvEaStore + i, vEa);
                    if (j+1<s2Len) {
                        simde__m128i cond = simde_mm_cmpgt_epi8(vEF_opn, vEa_ext);
                        simde__m128i vT = simde_mm_blendv_epi8(vTInsE, vTDiagE, cond);
                        arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                    }
                }
                if (! simde_mm_movemask_epi8(
                            simde_mm_or_si128(
                                simde_mm_cmpgt_epi8(vF_ext, vEF_opn),
                                simde_mm_cmpeq_epi8(vF_ext, vEF_opn))))
                    goto end;
                /*vF = simde_mm_max_epi8(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = simde_mm_subs_epi8(vFa, vGapE);
                vFa = simde_mm_max_epi8(vEF_opn, vFa_ext);
                vHp = simde_mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            simde__m128i vCompare;
            vH = simde_mm_load_si128(pvHStore + offset);
            vCompare = simde_mm_and_si128(vPosMask, simde_mm_cmpgt_epi8(vH, vMaxH));
            vMaxH = simde_mm_max_epi8(vH, vMaxH);
            if (simde_mm_movemask_epi8(vCompare)) {
                end_ref = j;
                end_query = s1Len - 1;
            }
        }
    }

    /* max last value from all columns */
    {
        for (k=0; k<position; ++k) {
            vMaxH = simde_mm_slli_si128(vMaxH, 1);
        }
        score = (int8_t) simde_mm_extract_epi8(vMaxH, 15);
    }

    /* max of last column */
    {
        int8_t score_last;
        vMaxH = vNegInf;

        for (i=0; i<segLen; ++i) {
            simde__m128i vH = simde_mm_load_si128(pvHStore + i);
            vMaxH = simde_mm_max_epi8(vH, vMaxH);
        }

        /* max in vec */
        score_last = simde_mm_hmax_epi8(vMaxH);
        if (score_last > score || (score_last == score && end_ref == s2Len - 1)) {
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                int8_t *t = (int8_t*)pvHStore;
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

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;

    parasail_free(pvHT);
    parasail_free(pvEaLoad);
    parasail_free(pvEaStore);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

