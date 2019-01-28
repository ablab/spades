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

#define FNAME parasail_nw_trace_scan_altivec_128_8
#define PNAME parasail_nw_trace_scan_profile_altivec_128_8

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
    vec128i* const restrict pvP  = (vec128i*)profile->profile8.score;
    vec128i* const restrict pvE  = parasail_memalign_vec128i(16, segLen);
    int8_t* const restrict boundary = parasail_memalign_int8_t(16, s2Len+1);
    vec128i* const restrict pvHt = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvH  = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvGapper = parasail_memalign_vec128i(16, segLen);
    vec128i vGapO = _mm_set1_epi8(open);
    vec128i vGapE = _mm_set1_epi8(gap);
    const int8_t NEG_LIMIT = (-open < matrix->min ?
        INT8_MIN + open : INT8_MIN - matrix->min) + 1;
    const int8_t POS_LIMIT = INT8_MAX - matrix->max - 1;
    vec128i vZero = _mm_setzero_si128();
    int8_t score = NEG_LIMIT;
    vec128i vNegLimit = _mm_set1_epi8(NEG_LIMIT);
    vec128i vPosLimit = _mm_set1_epi8(POS_LIMIT);
    vec128i vSaturationCheckMin = vPosLimit;
    vec128i vSaturationCheckMax = vNegLimit;
    vec128i vNegInfFront = vZero;
    vec128i vSegLenXgap;
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(vec128i));
    vec128i vTIns  = _mm_set1_epi8(PARASAIL_INS);
    vec128i vTDel  = _mm_set1_epi8(PARASAIL_DEL);
    vec128i vTDiag = _mm_set1_epi8(PARASAIL_DIAG);
    vec128i vTDiagE = _mm_set1_epi8(PARASAIL_DIAG_E);
    vec128i vTInsE = _mm_set1_epi8(PARASAIL_INS_E);
    vec128i vTDiagF = _mm_set1_epi8(PARASAIL_DIAG_F);
    vec128i vTDelF = _mm_set1_epi8(PARASAIL_DEL_F);

    vNegInfFront = _mm_insert_epi8(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = _mm_adds_epi8(vNegInfFront,
            _mm_slli_si128(_mm_set1_epi8(-segLen*gap), 1));

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
            _mm_store_si128(&pvH[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
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

    {
        vec128i vGapper = _mm_subs_epi8(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            _mm_store_si128(pvGapper+i, vGapper);
            vGapper = _mm_subs_epi8(vGapper, vGapE);
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
        vHp = _mm_slli_si128(vHp, 1);
        vHp = _mm_insert_epi8(vHp, boundary[j], 0);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = _mm_subs_epi8(vNegLimit, pvGapper[0]);
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
            vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vE);
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
        vHt = _mm_insert_epi8(vHt, boundary[j+1], 0);
        vGapper = _mm_load_si128(pvGapper);
        vGapper = _mm_adds_epi8(vHt, vGapper);
        vF = _mm_max_epi8(vF, vGapper);
        for (i=0; i<segWidth-2; ++i) {
            vec128i vFt = _mm_slli_si128(vF, 1);
            vFt = _mm_adds_epi8(vFt, vSegLenXgap);
            vF = _mm_max_epi8(vF, vFt);
        }

        /* calculate final H */
        vF = _mm_slli_si128(vF, 1);
        vF = _mm_adds_epi8(vF, vNegInfFront);
        vH = _mm_max_epi8(vF, vHt);
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
            case1 = _mm_cmpeq_epi8(vH, vHp);
            case2 = _mm_cmpeq_epi8(vH, vF);
            vT = _mm_blendv_epi8(
                    _mm_blendv_epi8(vTIns, vTDel, case2),
                    vTDiag, case1);
            vT = _mm_or_si128(vT, vET);
            vT = _mm_or_si128(vT, vFT);
            arr_store(result->trace->trace_table, vT, i, segLen, j);
            _mm_store_si128(pvH+i, vH);
            vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vH);
            vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vF);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vH);
        }
    }

    /* extract last value from the last column */
    {
        vec128i vH = _mm_load_si128(pvH + offset);
        for (k=0; k<position; ++k) {
            vH = _mm_slli_si128(vH, 1);
        }
        score = (int8_t) _mm_extract_epi8 (vH, 15);
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
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;

    parasail_free(pvGapper);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(boundary);
    parasail_free(pvE);

    return result;
}


