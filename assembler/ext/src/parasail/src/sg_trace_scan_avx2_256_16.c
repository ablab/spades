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


#define _mm256_cmplt_epi16_rpl(a,b) _mm256_cmpgt_epi16(b,a)

#if HAVE_AVX2_MM256_INSERT_EPI16
#define _mm256_insert_epi16_rpl _mm256_insert_epi16
#else
static inline __m256i _mm256_insert_epi16_rpl(__m256i a, int16_t i, int imm) {
    __m256i_16_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI16
#define _mm256_extract_epi16_rpl _mm256_extract_epi16
#else
static inline int16_t _mm256_extract_epi16_rpl(__m256i a, int imm) {
    __m256i_16_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

static inline int16_t _mm256_hmax_epi16_rpl(__m256i a) {
    a = _mm256_max_epi16(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi16(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi16(a, _mm256_slli_si256(a, 4));
    a = _mm256_max_epi16(a, _mm256_slli_si256(a, 2));
    return _mm256_extract_epi16_rpl(a, 15);
}


static inline void arr_store(
        __m256i *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    _mm256_store_si256(array + (1LL*d*seglen+t), vH);
}

static inline __m256i arr_load(
        __m256i *array,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    return _mm256_load_si256(array + (1LL*d*seglen+t));
}

#define FNAME parasail_sg_trace_scan_avx2_256_16
#define PNAME parasail_sg_trace_scan_profile_avx2_256_16

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_avx_256_16(s1, s1Len, matrix);
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
    __m256i* const restrict pvP  = (__m256i*)profile->profile16.score;
    __m256i* const restrict pvE  = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHt = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvH  = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvGapper = parasail_memalign___m256i(32, segLen);
    __m256i vGapO = _mm256_set1_epi16(open);
    __m256i vGapE = _mm256_set1_epi16(gap);
    const int16_t NEG_LIMIT = (-open < matrix->min ?
        INT16_MIN + open : INT16_MIN - matrix->min) + 1;
    const int16_t POS_LIMIT = INT16_MAX - matrix->max - 1;
    __m256i vZero = _mm256_setzero_si256();
    int16_t score = NEG_LIMIT;
    __m256i vNegLimit = _mm256_set1_epi16(NEG_LIMIT);
    __m256i vPosLimit = _mm256_set1_epi16(POS_LIMIT);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
    __m256i vMaxH = vNegLimit;
    __m256i vPosMask = _mm256_cmpeq_epi16(_mm256_set1_epi16(position),
            _mm256_set_epi16(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15));
    __m256i vNegInfFront = vZero;
    __m256i vSegLenXgap;
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, 32, sizeof(__m256i));
    __m256i vTIns  = _mm256_set1_epi16(PARASAIL_INS);
    __m256i vTDel  = _mm256_set1_epi16(PARASAIL_DEL);
    __m256i vTDiag = _mm256_set1_epi16(PARASAIL_DIAG);
    __m256i vTDiagE = _mm256_set1_epi16(PARASAIL_DIAG_E);
    __m256i vTInsE = _mm256_set1_epi16(PARASAIL_INS_E);
    __m256i vTDiagF = _mm256_set1_epi16(PARASAIL_DIAG_F);
    __m256i vTDelF = _mm256_set1_epi16(PARASAIL_DEL_F);

    vNegInfFront = _mm256_insert_epi16_rpl(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = _mm256_add_epi16(vNegInfFront,
            _mm256_slli_si256_rpl(_mm256_set1_epi16(-segLen*gap), 2));

    /* initialize H and E */
    parasail_memset___m256i(pvH, vZero, segLen);
    parasail_memset___m256i(pvE, vNegLimit, segLen);
    {
        __m256i vGapper = _mm256_sub_epi16(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            _mm256_store_si256(pvGapper+i, vGapper);
            vGapper = _mm256_sub_epi16(vGapper, vGapE);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        __m256i vE_ext;
        __m256i vE_opn;
        __m256i vHt;
        __m256i vF;
        __m256i vF_ext;
        __m256i vF_opn;
        __m256i vH;
        __m256i vHp;
        __m256i *pvW;
        __m256i vW;
        __m256i case1;
        __m256i case2;
        __m256i vGapper;
        __m256i vT;
        __m256i vET;
        __m256i vFT;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = _mm256_load_si256(pvH+(segLen-1));
        vHp = _mm256_slli_si256_rpl(vHp, 2);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = _mm256_sub_epi16(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        for (i=0; i<segLen; ++i) {
            vH = _mm256_load_si256(pvH+i);
            vE = _mm256_load_si256(pvE+i);
            vW = _mm256_load_si256(pvW+i);
            vGapper = _mm256_load_si256(pvGapper+i);
            vE_opn = _mm256_sub_epi16(vH, vGapO);
            vE_ext = _mm256_sub_epi16(vE, vGapE);
            case1 = _mm256_cmpgt_epi16(vE_opn, vE_ext);
            vET = _mm256_blendv_epi8(vTInsE, vTDiagE, case1);
            arr_store(result->trace->trace_table, vET, i, segLen, j);
            vE = _mm256_max_epi16(vE_opn, vE_ext);
            vSaturationCheckMin = _mm256_min_epi16(vSaturationCheckMin, vE);
            vGapper = _mm256_add_epi16(vHt, vGapper);
            vF = _mm256_max_epi16(vF, vGapper);
            vHp = _mm256_add_epi16(vHp, vW);
            vHt = _mm256_max_epi16(vE, vHp);
            _mm256_store_si256(pvE+i, vE);
            _mm256_store_si256(pvHt+i, vHt);
            _mm256_store_si256(pvH+i, vHp);
            vHp = vH;
        }

        /* pseudo prefix scan on F and H */
        vHt = _mm256_slli_si256_rpl(vHt, 2);
        vGapper = _mm256_load_si256(pvGapper);
        vGapper = _mm256_add_epi16(vHt, vGapper);
        vF = _mm256_max_epi16(vF, vGapper);
        for (i=0; i<segWidth-2; ++i) {
            __m256i vFt = _mm256_slli_si256_rpl(vF, 2);
            vFt = _mm256_add_epi16(vFt, vSegLenXgap);
            vF = _mm256_max_epi16(vF, vFt);
        }

        /* calculate final H */
        vF = _mm256_slli_si256_rpl(vF, 2);
        vF = _mm256_add_epi16(vF, vNegInfFront);
        vH = _mm256_max_epi16(vF, vHt);
        for (i=0; i<segLen; ++i) {
            vET = arr_load(result->trace->trace_table, i, segLen, j);
            vHp = _mm256_load_si256(pvH+i);
            vHt = _mm256_load_si256(pvHt+i);
            vF_opn = _mm256_sub_epi16(vH, vGapO);
            vF_ext = _mm256_sub_epi16(vF, vGapE);
            vF = _mm256_max_epi16(vF_opn, vF_ext);
            case1 = _mm256_cmpgt_epi16(vF_opn, vF_ext);
            vFT = _mm256_blendv_epi8(vTDelF, vTDiagF, case1);
            vH = _mm256_max_epi16(vHt, vF);
            case1 = _mm256_cmpeq_epi16(vH, vHp);
            case2 = _mm256_cmpeq_epi16(vH, vF);
            vT = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vTIns, vTDel, case2),
                    vTDiag, case1);
            vT = _mm256_or_si256(vT, vET);
            vT = _mm256_or_si256(vT, vFT);
            arr_store(result->trace->trace_table, vT, i, segLen, j);
            _mm256_store_si256(pvH+i, vH);
            vSaturationCheckMin = _mm256_min_epi16(vSaturationCheckMin, vH);
            vSaturationCheckMin = _mm256_min_epi16(vSaturationCheckMin, vF);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vH);
        }

        /* extract vector containing last value from column */
        {
            __m256i vCompare;
            vH = _mm256_load_si256(pvH + offset);
            vCompare = _mm256_and_si256(vPosMask, _mm256_cmpgt_epi16(vH, vMaxH));
            vMaxH = _mm256_max_epi16(vH, vMaxH);
            if (_mm256_movemask_epi8(vCompare)) {
                end_ref = j;
                end_query = s1Len - 1;
            }
        }
    } 

    /* max last value from all columns */
    {
        int16_t value;
        for (k=0; k<position; ++k) {
            vMaxH = _mm256_slli_si256_rpl(vMaxH, 2);
        }
        value = (int16_t) _mm256_extract_epi16_rpl(vMaxH, 15);
        if (value > score) {
            score = value;
        }
    }

    /* max of last column */
    {
        int16_t score_last;
        vMaxH = vNegLimit;

        for (i=0; i<segLen; ++i) {
            __m256i vH = _mm256_load_si256(pvH + i);
            vMaxH = _mm256_max_epi16(vH, vMaxH);
        }

        /* max in vec */
        score_last = _mm256_hmax_epi16_rpl(vMaxH);
        if (score_last > score || (score_last == score && end_ref == s2Len - 1)) {
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                int16_t *t = (int16_t*)pvH;
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
            _mm256_cmplt_epi16_rpl(vSaturationCheckMin, vNegLimit),
            _mm256_cmpgt_epi16(vSaturationCheckMax, vPosLimit)))) {
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
        | PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_16;

    parasail_free(pvGapper);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}


