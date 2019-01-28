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

#define SWAP(A,B) { __m256i* tmp = A; A = B; B = tmp; }

#define NEG_INF (INT64_MIN/(int64_t)(2))

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

#if HAVE_AVX2_MM256_EXTRACT_EPI64
#define _mm256_extract_epi64_rpl _mm256_extract_epi64
#else
static inline int64_t _mm256_extract_epi64_rpl(__m256i a, int imm) {
    __m256i_64_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)


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

#define FNAME parasail_nw_trace_striped_avx2_256_64
#define PNAME parasail_nw_trace_striped_profile_avx2_256_64

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
    const int s1Len = profile->s1Len;
    int32_t end_query = s1Len-1;
    int32_t end_ref = s2Len-1;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 4; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict vProfile = (__m256i*)profile->profile64.score;
    __m256i* restrict pvHStore = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLoad =  parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvE = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvEaStore = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvEaLoad = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHT = parasail_memalign___m256i(32, segLen);
    int64_t* const restrict boundary = parasail_memalign_int64_t(32, s2Len+1);
    __m256i vGapO = _mm256_set1_epi64x_rpl(open);
    __m256i vGapE = _mm256_set1_epi64x_rpl(gap);
    __m256i vNegInf = _mm256_set1_epi64x_rpl(NEG_INF);
    int64_t score = NEG_INF;
    
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, 32, sizeof(__m256i));
    __m256i vTIns  = _mm256_set1_epi64x_rpl(PARASAIL_INS);
    __m256i vTDel  = _mm256_set1_epi64x_rpl(PARASAIL_DEL);
    __m256i vTDiag = _mm256_set1_epi64x_rpl(PARASAIL_DIAG);
    __m256i vTDiagE = _mm256_set1_epi64x_rpl(PARASAIL_DIAG_E);
    __m256i vTInsE = _mm256_set1_epi64x_rpl(PARASAIL_INS_E);
    __m256i vTDiagF = _mm256_set1_epi64x_rpl(PARASAIL_DIAG_F);
    __m256i vTDelF = _mm256_set1_epi64x_rpl(PARASAIL_DEL_F);
    __m256i vTMask = _mm256_set1_epi64x_rpl(PARASAIL_ZERO_MASK);
    __m256i vFTMask = _mm256_set1_epi64x_rpl(PARASAIL_F_MASK);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            __m256i_64_t h;
            __m256i_64_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT64_MIN ? INT64_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT64_MIN ? INT64_MIN : tmp;
            }
            _mm256_store_si256(&pvHStore[index], h.m);
            _mm256_store_si256(&pvE[index], e.m);
            _mm256_store_si256(&pvEaStore[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT64_MIN ? INT64_MIN : tmp;
        }
    }

    for (i=0; i<segLen; ++i) {
        arr_store(result->trace->trace_table, vTDiagE, i, segLen, 0);
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vEF_opn;
        __m256i vE;
        __m256i vE_ext;
        __m256i vF;
        __m256i vF_ext;
        __m256i vFa;
        __m256i vFa_ext;
        __m256i vH;
        __m256i vH_dag;
        const __m256i* vP = NULL;

        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vF = vNegInf;

        /* load final segment of pvHStore and shift left by 8 bytes */
        vH = _mm256_load_si256(&pvHStore[segLen - 1]);
        vH = _mm256_slli_si256_rpl(vH, 8);

        /* insert upper boundary condition */
        vH = _mm256_insert_epi64_rpl(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad, pvHStore)
        SWAP(pvEaLoad, pvEaStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = _mm256_load_si256(pvE + i);

            /* Get max from vH, vE and vF. */
            vH_dag = _mm256_add_epi64(vH, _mm256_load_si256(vP + i));
            vH = _mm256_max_epi64_rpl(vH_dag, vE);
            vH = _mm256_max_epi64_rpl(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
            

            {
                __m256i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                __m256i case1 = _mm256_cmpeq_epi64(vH, vH_dag);
                __m256i case2 = _mm256_cmpeq_epi64(vH, vF);
                __m256i vT = _mm256_blendv_epi8(
                        _mm256_blendv_epi8(vTIns, vTDel, case2),
                        vTDiag, case1);
                _mm256_store_si256(pvHT + i, vT);
                vT = _mm256_or_si256(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i, segLen, j);
            }

            vEF_opn = _mm256_sub_epi64(vH, vGapO);

            /* Update vE value. */
            vE_ext = _mm256_sub_epi64(vE, vGapE);
            vE = _mm256_max_epi64_rpl(vEF_opn, vE_ext);
            _mm256_store_si256(pvE + i, vE);
            {
                __m256i vEa = _mm256_load_si256(pvEaLoad + i);
                __m256i vEa_ext = _mm256_sub_epi64(vEa, vGapE);
                vEa = _mm256_max_epi64_rpl(vEF_opn, vEa_ext);
                _mm256_store_si256(pvEaStore + i, vEa);
                if (j+1<s2Len) {
                    __m256i cond = _mm256_cmpgt_epi64(vEF_opn, vEa_ext);
                    __m256i vT = _mm256_blendv_epi8(vTInsE, vTDiagE, cond);
                    arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                }
            }

            /* Update vF value. */
            vF_ext = _mm256_sub_epi64(vF, vGapE);
            vF = _mm256_max_epi64_rpl(vEF_opn, vF_ext);
            if (i+1<segLen) {
                __m256i vTAll = arr_load(result->trace->trace_table, i+1, segLen, j);
                __m256i cond = _mm256_cmpgt_epi64(vEF_opn, vF_ext);
                __m256i vT = _mm256_blendv_epi8(vTDelF, vTDiagF, cond);
                vT = _mm256_or_si256(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i+1, segLen, j);
            }

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
        }


        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            int64_t tmp = boundary[j+1]-open;
            int64_t tmp2 = tmp < INT64_MIN ? INT64_MIN : tmp;
            __m256i vHp = _mm256_load_si256(&pvHLoad[segLen - 1]);
            vHp = _mm256_slli_si256_rpl(vHp, 8);
            vHp = _mm256_insert_epi64_rpl(vHp, boundary[j], 0);
            vEF_opn = _mm256_slli_si256_rpl(vEF_opn, 8);
            vEF_opn = _mm256_insert_epi64_rpl(vEF_opn, tmp2, 0);
            vF_ext = _mm256_slli_si256_rpl(vF_ext, 8);
            vF_ext = _mm256_insert_epi64_rpl(vF_ext, NEG_INF, 0);
            vF = _mm256_slli_si256_rpl(vF, 8);
            vF = _mm256_insert_epi64_rpl(vF, tmp2, 0);
            vFa_ext = _mm256_slli_si256_rpl(vFa_ext, 8);
            vFa_ext = _mm256_insert_epi64_rpl(vFa_ext, NEG_INF, 0);
            vFa = _mm256_slli_si256_rpl(vFa, 8);
            vFa = _mm256_insert_epi64_rpl(vFa, tmp2, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi64_rpl(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                
                {
                    __m256i vTAll;
                    __m256i vT;
                    __m256i case1;
                    __m256i case2;
                    __m256i cond;
                    vHp = _mm256_add_epi64(vHp, _mm256_load_si256(vP + i));
                    case1 = _mm256_cmpeq_epi64(vH, vHp);
                    case2 = _mm256_cmpeq_epi64(vH, vF);
                    cond = _mm256_andnot_si256(case1,case2);
                    vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    vT = _mm256_load_si256(pvHT + i);
                    vT = _mm256_blendv_epi8(vT, vTDel, cond);
                    _mm256_store_si256(pvHT + i, vT);
                    vTAll = _mm256_and_si256(vTAll, vTMask);
                    vTAll = _mm256_or_si256(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                /* Update vF value. */
                {
                    __m256i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    __m256i cond = _mm256_cmpgt_epi64(vEF_opn, vFa_ext);
                    __m256i vT = _mm256_blendv_epi8(vTDelF, vTDiagF, cond);
                    vTAll = _mm256_and_si256(vTAll, vFTMask);
                    vTAll = _mm256_or_si256(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                vEF_opn = _mm256_sub_epi64(vH, vGapO);
                vF_ext = _mm256_sub_epi64(vF, vGapE);
                {
                    __m256i vEa = _mm256_load_si256(pvEaLoad + i);
                    __m256i vEa_ext = _mm256_sub_epi64(vEa, vGapE);
                    vEa = _mm256_max_epi64_rpl(vEF_opn, vEa_ext);
                    _mm256_store_si256(pvEaStore + i, vEa);
                    if (j+1<s2Len) {
                        __m256i cond = _mm256_cmpgt_epi64(vEF_opn, vEa_ext);
                        __m256i vT = _mm256_blendv_epi8(vTInsE, vTDiagE, cond);
                        arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                    }
                }
                if (! _mm256_movemask_epi8(
                            _mm256_or_si256(
                                _mm256_cmpgt_epi64(vF_ext, vEF_opn),
                                _mm256_cmpeq_epi64(vF_ext, vEF_opn))))
                    goto end;
                /*vF = _mm256_max_epi64_rpl(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = _mm256_sub_epi64(vFa, vGapE);
                vFa = _mm256_max_epi64_rpl(vEF_opn, vFa_ext);
                vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
end:
        {
        }
    }

    /* extract last value from the last column */
    {
        __m256i vH = _mm256_load_si256(pvHStore + offset);
        for (k=0; k<position; ++k) {
            vH = _mm256_slli_si256_rpl (vH, 8);
        }
        score = (int64_t) _mm256_extract_epi64_rpl (vH, 3);
    }

    

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_4;

    parasail_free(boundary);
    parasail_free(pvHT);
    parasail_free(pvEaLoad);
    parasail_free(pvEaStore);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}


