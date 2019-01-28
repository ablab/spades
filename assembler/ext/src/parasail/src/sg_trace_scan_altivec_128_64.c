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

#define FNAME parasail_sg_trace_scan_altivec_128_64
#define PNAME parasail_sg_trace_scan_profile_altivec_128_64

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_altivec_128_64(s1, s1Len, matrix);
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
    const int32_t segWidth = 2; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    vec128i* const restrict pvP  = (vec128i*)profile->profile64.score;
    vec128i* const restrict pvE  = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvHt = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvH  = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvGapper = parasail_memalign_vec128i(16, segLen);
    vec128i vGapO = _mm_set1_epi64(open);
    vec128i vGapE = _mm_set1_epi64(gap);
    const int64_t NEG_LIMIT = (-open < matrix->min ?
        INT64_MIN + open : INT64_MIN - matrix->min) + 1;
    const int64_t POS_LIMIT = INT64_MAX - matrix->max - 1;
    vec128i vZero = _mm_setzero_si128();
    int64_t score = NEG_LIMIT;
    vec128i vNegLimit = _mm_set1_epi64(NEG_LIMIT);
    vec128i vPosLimit = _mm_set1_epi64(POS_LIMIT);
    vec128i vSaturationCheckMin = vPosLimit;
    vec128i vSaturationCheckMax = vNegLimit;
    vec128i vMaxH = vNegLimit;
    vec128i vPosMask = _mm_cmpeq_epi64(_mm_set1_epi64(position),
            _mm_set_epi64(0,1));
    vec128i vNegInfFront = vZero;
    vec128i vSegLenXgap;
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(vec128i));
    vec128i vTIns  = _mm_set1_epi64(PARASAIL_INS);
    vec128i vTDel  = _mm_set1_epi64(PARASAIL_DEL);
    vec128i vTDiag = _mm_set1_epi64(PARASAIL_DIAG);
    vec128i vTDiagE = _mm_set1_epi64(PARASAIL_DIAG_E);
    vec128i vTInsE = _mm_set1_epi64(PARASAIL_INS_E);
    vec128i vTDiagF = _mm_set1_epi64(PARASAIL_DIAG_F);
    vec128i vTDelF = _mm_set1_epi64(PARASAIL_DEL_F);

    vNegInfFront = _mm_insert_epi64(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = _mm_add_epi64(vNegInfFront,
            _mm_slli_si128(_mm_set1_epi64(-segLen*gap), 8));

    /* initialize H and E */
    parasail_memset_vec128i(pvH, vZero, segLen);
    parasail_memset_vec128i(pvE, vNegLimit, segLen);
    {
        vec128i vGapper = _mm_sub_epi64(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            _mm_store_si128(pvGapper+i, vGapper);
            vGapper = _mm_sub_epi64(vGapper, vGapE);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        vec128i vE;
        vec128i vE_ext;
        vec128i vE_opn;
        vec128i vHt;
        vec128i vF;
        vec128i vF_ext;
        vec128i vF_opn;
        vec128i vH;
        vec128i vHp;
        vec128i *pvW;
        vec128i vW;
        vec128i case1;
        vec128i case2;
        vec128i vGapper;
        vec128i vT;
        vec128i vET;
        vec128i vFT;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = _mm_load_si128(pvH+(segLen-1));
        vHp = _mm_slli_si128(vHp, 8);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = _mm_sub_epi64(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vGapper = _mm_load_si128(pvGapper+i);
            vE_opn = _mm_sub_epi64(vH, vGapO);
            vE_ext = _mm_sub_epi64(vE, vGapE);
            case1 = _mm_cmpgt_epi64(vE_opn, vE_ext);
            vET = _mm_blendv_epi8(vTInsE, vTDiagE, case1);
            arr_store(result->trace->trace_table, vET, i, segLen, j);
            vE = _mm_max_epi64(vE_opn, vE_ext);
            vSaturationCheckMin = _mm_min_epi64(vSaturationCheckMin, vE);
            vGapper = _mm_add_epi64(vHt, vGapper);
            vF = _mm_max_epi64(vF, vGapper);
            vHp = _mm_add_epi64(vHp, vW);
            vHt = _mm_max_epi64(vE, vHp);
            _mm_store_si128(pvE+i, vE);
            _mm_store_si128(pvHt+i, vHt);
            _mm_store_si128(pvH+i, vHp);
            vHp = vH;
        }

        /* pseudo prefix scan on F and H */
        vHt = _mm_slli_si128(vHt, 8);
        vGapper = _mm_load_si128(pvGapper);
        vGapper = _mm_add_epi64(vHt, vGapper);
        vF = _mm_max_epi64(vF, vGapper);
        for (i=0; i<segWidth-2; ++i) {
            vec128i vFt = _mm_slli_si128(vF, 8);
            vFt = _mm_add_epi64(vFt, vSegLenXgap);
            vF = _mm_max_epi64(vF, vFt);
        }

        /* calculate final H */
        vF = _mm_slli_si128(vF, 8);
        vF = _mm_add_epi64(vF, vNegInfFront);
        vH = _mm_max_epi64(vF, vHt);
        for (i=0; i<segLen; ++i) {
            vET = arr_load(result->trace->trace_table, i, segLen, j);
            vHp = _mm_load_si128(pvH+i);
            vHt = _mm_load_si128(pvHt+i);
            vF_opn = _mm_sub_epi64(vH, vGapO);
            vF_ext = _mm_sub_epi64(vF, vGapE);
            vF = _mm_max_epi64(vF_opn, vF_ext);
            case1 = _mm_cmpgt_epi64(vF_opn, vF_ext);
            vFT = _mm_blendv_epi8(vTDelF, vTDiagF, case1);
            vH = _mm_max_epi64(vHt, vF);
            case1 = _mm_cmpeq_epi64(vH, vHp);
            case2 = _mm_cmpeq_epi64(vH, vF);
            vT = _mm_blendv_epi8(
                    _mm_blendv_epi8(vTIns, vTDel, case2),
                    vTDiag, case1);
            vT = _mm_or_si128(vT, vET);
            vT = _mm_or_si128(vT, vFT);
            arr_store(result->trace->trace_table, vT, i, segLen, j);
            _mm_store_si128(pvH+i, vH);
            vSaturationCheckMin = _mm_min_epi64(vSaturationCheckMin, vH);
            vSaturationCheckMin = _mm_min_epi64(vSaturationCheckMin, vF);
            vSaturationCheckMax = _mm_max_epi64(vSaturationCheckMax, vH);
        }

        /* extract vector containing last value from column */
        {
            vec128i vCompare;
            vH = _mm_load_si128(pvH + offset);
            vCompare = _mm_and_si128(vPosMask, _mm_cmpgt_epi64(vH, vMaxH));
            vMaxH = _mm_max_epi64(vH, vMaxH);
            if (_mm_movemask_epi8(vCompare)) {
                end_ref = j;
                end_query = s1Len - 1;
            }
        }
    } 

    /* max last value from all columns */
    {
        int64_t value;
        for (k=0; k<position; ++k) {
            vMaxH = _mm_slli_si128(vMaxH, 8);
        }
        value = (int64_t) _mm_extract_epi64(vMaxH, 1);
        if (value > score) {
            score = value;
        }
    }

    /* max of last column */
    {
        int64_t score_last;
        vMaxH = vNegLimit;

        for (i=0; i<segLen; ++i) {
            vec128i vH = _mm_load_si128(pvH + i);
            vMaxH = _mm_max_epi64(vH, vMaxH);
        }

        /* max in vec */
        score_last = _mm_hmax_epi64(vMaxH);
        if (score_last > score || (score_last == score && end_ref == s2Len - 1)) {
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                int64_t *t = (int64_t*)pvH;
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

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmplt_epi64(vSaturationCheckMin, vNegLimit),
            _mm_cmpgt_epi64(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;

    parasail_free(pvGapper);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}


