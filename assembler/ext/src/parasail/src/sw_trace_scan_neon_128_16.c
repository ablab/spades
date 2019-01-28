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



#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_neon.h"



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

#define FNAME parasail_sw_trace_scan_neon_128_16
#define PNAME parasail_sw_trace_scan_profile_neon_128_16

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_neon_128_16(s1, s1Len, matrix);
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
    const int32_t segWidth = 8; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict pvP  = (simde__m128i*)profile->profile16.score;
    simde__m128i* const restrict pvE  = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHt = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvH  = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHMax  = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvGapper = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i vGapO = simde_mm_set1_epi16(open);
    simde__m128i vGapE = simde_mm_set1_epi16(gap);
    const int16_t NEG_LIMIT = (-open < matrix->min ?
        INT16_MIN + open : INT16_MIN - matrix->min) + 1;
    const int16_t POS_LIMIT = INT16_MAX - matrix->max - 1;
    simde__m128i vZero = simde_mm_setzero_si128();
    int16_t score = NEG_LIMIT;
    simde__m128i vNegLimit = simde_mm_set1_epi16(NEG_LIMIT);
    simde__m128i vPosLimit = simde_mm_set1_epi16(POS_LIMIT);
    simde__m128i vSaturationCheckMin = vPosLimit;
    simde__m128i vSaturationCheckMax = vNegLimit;
    simde__m128i vMaxH = vNegLimit;
    simde__m128i vMaxHUnit = vNegLimit;
    simde__m128i vNegInfFront = vZero;
    simde__m128i vSegLenXgap;
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(simde__m128i));
    simde__m128i vTZero = simde_mm_set1_epi16(PARASAIL_ZERO);
    simde__m128i vTIns  = simde_mm_set1_epi16(PARASAIL_INS);
    simde__m128i vTDel  = simde_mm_set1_epi16(PARASAIL_DEL);
    simde__m128i vTDiag = simde_mm_set1_epi16(PARASAIL_DIAG);
    simde__m128i vTDiagE = simde_mm_set1_epi16(PARASAIL_DIAG_E);
    simde__m128i vTInsE = simde_mm_set1_epi16(PARASAIL_INS_E);
    simde__m128i vTDiagF = simde_mm_set1_epi16(PARASAIL_DIAG_F);
    simde__m128i vTDelF = simde_mm_set1_epi16(PARASAIL_DEL_F);

    vNegInfFront = simde_mm_insert_epi16(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = simde_mm_add_epi16(vNegInfFront,
            simde_mm_slli_si128(simde_mm_set1_epi16(-segLen*gap), 2));

    parasail_memset_simde__m128i(pvH, vZero, segLen);
    parasail_memset_simde__m128i(pvE, vNegLimit, segLen);
    {
        simde__m128i vGapper = simde_mm_sub_epi16(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            simde_mm_store_si128(pvGapper+i, vGapper);
            vGapper = simde_mm_sub_epi16(vGapper, vGapE);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vE;
        simde__m128i vE_ext;
        simde__m128i vE_opn;
        simde__m128i vHt;
        simde__m128i vF;
        simde__m128i vF_ext;
        simde__m128i vF_opn;
        simde__m128i vH;
        simde__m128i vHp;
        simde__m128i *pvW;
        simde__m128i vW;
        simde__m128i case1;
        simde__m128i case2;
        simde__m128i case0;
        simde__m128i vGapper;
        simde__m128i vT;
        simde__m128i vET;
        simde__m128i vFT;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = simde_mm_load_si128(pvH+(segLen-1));
        vHp = simde_mm_slli_si128(vHp, 2);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vZero;
        vF = vNegLimit;
        for (i=0; i<segLen; ++i) {
            vH = simde_mm_load_si128(pvH+i);
            vE = simde_mm_load_si128(pvE+i);
            vW = simde_mm_load_si128(pvW+i);
            vGapper = simde_mm_load_si128(pvGapper+i);
            vE_opn = simde_mm_sub_epi16(vH, vGapO);
            vE_ext = simde_mm_sub_epi16(vE, vGapE);
            case1 = simde_mm_cmpgt_epi16(vE_opn, vE_ext);
            vET = simde_mm_blendv_epi8(vTInsE, vTDiagE, case1);
            arr_store(result->trace->trace_table, vET, i, segLen, j);
            vE = simde_mm_max_epi16(vE_opn, vE_ext);
            vGapper = simde_mm_add_epi16(vHt, vGapper);
            vF = simde_mm_max_epi16(vF, vGapper);
            vHp = simde_mm_add_epi16(vHp, vW);
            vHt = simde_mm_max_epi16(vE, vHp);
            simde_mm_store_si128(pvE+i, vE);
            simde_mm_store_si128(pvHt+i, vHt);
            simde_mm_store_si128(pvH+i, vHp);
            vHp = vH;
        }

        /* pseudo prefix scan on F and H */
        vHt = simde_mm_slli_si128(vHt, 2);
        vGapper = simde_mm_load_si128(pvGapper);
        vGapper = simde_mm_add_epi16(vHt, vGapper);
        vF = simde_mm_max_epi16(vF, vGapper);
        for (i=0; i<segWidth-2; ++i) {
            simde__m128i vFt = simde_mm_slli_si128(vF, 2);
            vFt = simde_mm_add_epi16(vFt, vSegLenXgap);
            vF = simde_mm_max_epi16(vF, vFt);
        }

        /* calculate final H */
        vF = simde_mm_slli_si128(vF, 2);
        vF = simde_mm_add_epi16(vF, vNegInfFront);
        vH = simde_mm_max_epi16(vF, vHt);
        vH = simde_mm_max_epi16(vH, vZero);
        for (i=0; i<segLen; ++i) {
            vET = arr_load(result->trace->trace_table, i, segLen, j);
            vHp = simde_mm_load_si128(pvH+i);
            vHt = simde_mm_load_si128(pvHt+i);
            vF_opn = simde_mm_sub_epi16(vH, vGapO);
            vF_ext = simde_mm_sub_epi16(vF, vGapE);
            vF = simde_mm_max_epi16(vF_opn, vF_ext);
            case1 = simde_mm_cmpgt_epi16(vF_opn, vF_ext);
            vFT = simde_mm_blendv_epi8(vTDelF, vTDiagF, case1);
            vH = simde_mm_max_epi16(vHt, vF);
            vH = simde_mm_max_epi16(vH, vZero);
            case0 = simde_mm_cmpeq_epi16(vH, vZero);
            case1 = simde_mm_cmpeq_epi16(vH, vHp);
            case2 = simde_mm_cmpeq_epi16(vH, vF);
            vT = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vTIns, vTDel, case2),
                    vTDiag, case1);
            vT = simde_mm_blendv_epi8(vT, vTZero, case0);
            vT = simde_mm_or_si128(vT, vET);
            vT = simde_mm_or_si128(vT, vFT);
            arr_store(result->trace->trace_table, vT, i, segLen, j);
            simde_mm_store_si128(pvH+i, vH);
            vSaturationCheckMin = simde_mm_min_epi16(vSaturationCheckMin, vH);
            vSaturationCheckMax = simde_mm_max_epi16(vSaturationCheckMax, vH);
            {
                vMaxH = simde_mm_max_epi16(vH, vMaxH);
            }
        } 

        {
            simde__m128i vCompare = simde_mm_cmpgt_epi16(vMaxH, vMaxHUnit);
            if (simde_mm_movemask_epi8(vCompare)) {
                score = simde_mm_hmax_epi16(vMaxH);
                vMaxHUnit = simde_mm_set1_epi16(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvH, sizeof(simde__m128i)*segLen);
            }
        }
    }

    /* Trace the alignment ending position on read. */
    {
        int16_t *t = (int16_t*)pvHMax;
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

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmplt_epi16(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpgt_epi16(vSaturationCheckMax, vPosLimit)))) {
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
        | PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_8;

    parasail_free(pvGapper);
    parasail_free(pvHMax);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}


