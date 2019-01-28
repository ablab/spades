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

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <emmintrin.h>
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_sse.h"

#define FASTSTATS

#define SWAP(A,B) { __m128i* tmp = A; A = B; B = tmp; }
#define SWAP3(A,B,C) { __m128i* tmp = A; A = B; B = C; C = tmp; }

#define NEG_INF INT8_MIN

static inline __m128i _mm_blendv_epi8_rpl(__m128i a, __m128i b, __m128i mask) {
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
}

static inline __m128i _mm_max_epi8_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi8(a, b);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
}

static inline int8_t _mm_extract_epi8_rpl(__m128i a, const int imm) {
    __m128i_8_t A;
    A.m = a;
    return A.v[imm];
}

static inline int8_t _mm_hmax_epi8_rpl(__m128i a) {
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 2));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 1));
    return _mm_extract_epi8_rpl(a, 0);
}


#ifdef PARASAIL_TABLE
static inline void arr_store(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen,
        int32_t bias)
{
    array[1LL*( 0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  0) - bias;
    array[1LL*( 1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  1) - bias;
    array[1LL*( 2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  2) - bias;
    array[1LL*( 3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  3) - bias;
    array[1LL*( 4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  4) - bias;
    array[1LL*( 5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  5) - bias;
    array[1LL*( 6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  6) - bias;
    array[1LL*( 7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  7) - bias;
    array[1LL*( 8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  8) - bias;
    array[1LL*( 9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  9) - bias;
    array[1LL*(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 10) - bias;
    array[1LL*(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 11) - bias;
    array[1LL*(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 12) - bias;
    array[1LL*(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 13) - bias;
    array[1LL*(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 14) - bias;
    array[1LL*(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 15) - bias;
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t bias)
{
    col[ 0*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  0) - bias;
    col[ 1*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  1) - bias;
    col[ 2*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  2) - bias;
    col[ 3*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  3) - bias;
    col[ 4*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  4) - bias;
    col[ 5*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  5) - bias;
    col[ 6*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  6) - bias;
    col[ 7*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  7) - bias;
    col[ 8*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  8) - bias;
    col[ 9*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  9) - bias;
    col[10*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 10) - bias;
    col[11*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 11) - bias;
    col[12*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 12) - bias;
    col[13*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 13) - bias;
    col[14*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 14) - bias;
    col[15*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 15) - bias;
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_striped_sse2_128_8
#define PNAME parasail_sw_stats_table_striped_profile_sse2_128_8
#define INAME PNAME
#define STATIC
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_striped_sse2_128_8
#define PNAME parasail_sw_stats_rowcol_striped_profile_sse2_128_8
#define INAME PNAME
#define STATIC
#else
#define FNAME parasail_sw_stats_striped_sse2_128_8
#ifdef FASTSTATS
#define PNAME parasail_sw_stats_striped_profile_sse2_128_8_internal
#define INAME parasail_sw_stats_striped_profile_sse2_128_8
#define STATIC static
#else
#define PNAME parasail_sw_stats_striped_profile_sse2_128_8
#define INAME PNAME
#define STATIC
#endif
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_sse_128_8(s1, s1Len, matrix);
    parasail_result_t *result = INAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

STATIC parasail_result_t* PNAME(
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
    __m128i* const restrict vProfile  = (__m128i*)profile->profile8.score;
    __m128i* const restrict vProfileM = (__m128i*)profile->profile8.matches;
    __m128i* const restrict vProfileS = (__m128i*)profile->profile8.similar;
    __m128i* restrict pvHStore        = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHLoad         = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHMStore       = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHMLoad        = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHSStore       = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHSLoad        = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHLStore       = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHLLoad        = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvE       = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvEM      = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvES      = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvEL      = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHMax          = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHMMax         = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHSMax         = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHLMax         = parasail_memalign___m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi8(1);
    int8_t bias = INT8_MIN;
    int8_t score = bias;
    int8_t matches = bias;
    int8_t similar = bias;
    int8_t length = bias;
    __m128i vBias = _mm_set1_epi8(bias);
    __m128i vBias1 = _mm_adds_epi8(vBias,vOne);
    __m128i vMaxH = vBias;
    __m128i vMaxHUnit = vBias;
    __m128i insert_mask = _mm_cmpgt_epi8(
            _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
            vZero);
    __m128i vSaturationCheckMax = vBias;
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
    int8_t maxp = INT8_MAX - (int8_t)(matrix->max+1);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    parasail_memset___m128i(pvHStore, vBias, segLen);
    parasail_memset___m128i(pvHMStore, vBias, segLen);
    parasail_memset___m128i(pvHSStore, vBias, segLen);
    parasail_memset___m128i(pvHLStore, vBias, segLen);
    parasail_memset___m128i(pvE, vBias, segLen);
    parasail_memset___m128i(pvEM, vBias, segLen);
    parasail_memset___m128i(pvES, vBias, segLen);
    parasail_memset___m128i(pvEL, vBias1, segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vEF_opn;
        __m128i vE;
        __m128i vE_ext;
        __m128i vEM;
        __m128i vES;
        __m128i vEL;
        __m128i vF;
        __m128i vF_ext;
        __m128i vFM;
        __m128i vFS;
        __m128i vFL;
        __m128i vH;
        __m128i vH_dag;
        __m128i vHM;
        __m128i vHS;
        __m128i vHL;
        const __m128i* vP = NULL;
        const __m128i* vPM = NULL;
        const __m128i* vPS = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = vBias;
        vFM = vBias;
        vFS = vBias;
        vFL = vBias1;

        /* load final segment of pvHStore and shift left by 1 bytes */
        vH = _mm_load_si128(&pvHStore[segLen - 1]);
        vHM = _mm_load_si128(&pvHMStore[segLen - 1]);
        vHS = _mm_load_si128(&pvHSStore[segLen - 1]);
        vHL = _mm_load_si128(&pvHLStore[segLen - 1]);
        vH = _mm_slli_si128(vH, 1);
        vHM = _mm_slli_si128(vHM, 1);
        vHS = _mm_slli_si128(vHS, 1);
        vHL = _mm_slli_si128(vHL, 1);
        vH = _mm_blendv_epi8_rpl(vH, vBias, insert_mask);
        vHM = _mm_blendv_epi8_rpl(vHM, vBias, insert_mask);
        vHS = _mm_blendv_epi8_rpl(vHS, vBias, insert_mask);
        vHL = _mm_blendv_epi8_rpl(vHL, vBias, insert_mask);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPM = vProfileM + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + matrix->mapper[(unsigned char)s2[j]] * segLen;

        if (end_ref == j-2) {
            /* Swap in the max buffer. */
            SWAP3(pvHMax,  pvHLoad,  pvHStore)
            SWAP3(pvHMMax, pvHMLoad, pvHMStore)
            SWAP3(pvHSMax, pvHSLoad, pvHSStore)
            SWAP3(pvHLMax, pvHLLoad, pvHLStore)
        }
        else {
            /* Swap the 2 H buffers. */
            SWAP(pvHLoad,  pvHStore)
            SWAP(pvHMLoad, pvHMStore)
            SWAP(pvHSLoad, pvHSStore)
            SWAP(pvHLLoad, pvHLStore)
        }

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            __m128i cond_zero;
            __m128i case1;
            __m128i case2;

            vE = _mm_load_si128(pvE+ i);
            vEM = _mm_load_si128(pvEM+ i);
            vES = _mm_load_si128(pvES+ i);
            vEL = _mm_load_si128(pvEL+ i);

            /* Get max from vH, vE and vF. */
            vH_dag = _mm_adds_epi8(vH, _mm_load_si128(vP + i));
            /*vH = _mm_max_epi8_rpl(vH_dag, vBias);*/
            vH = _mm_max_epi8_rpl(vH_dag, vE);
            vH = _mm_max_epi8_rpl(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
            cond_zero = _mm_cmpeq_epi8(vH, vBias);

            case1 = _mm_cmpeq_epi8(vH, vH_dag);
            case2 = _mm_cmpeq_epi8(vH, vF);

            /* calculate vM */
            vHM = _mm_blendv_epi8_rpl(
                    _mm_blendv_epi8_rpl(vEM, vFM, case2),
                    _mm_adds_epi8(vHM, _mm_load_si128(vPM + i)), case1);
            vHM = _mm_blendv_epi8_rpl(vHM, vBias, cond_zero);
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vS */
            vHS = _mm_blendv_epi8_rpl(
                    _mm_blendv_epi8_rpl(vES, vFS, case2),
                    _mm_adds_epi8(vHS, _mm_load_si128(vPS + i)), case1);
            vHS = _mm_blendv_epi8_rpl(vHS, vBias, cond_zero);
            _mm_store_si128(pvHSStore + i, vHS);

            /* calculate vL */
            vHL = _mm_blendv_epi8_rpl(
                    _mm_blendv_epi8_rpl(vEL, vFL, case2),
                    _mm_adds_epi8(vHL, vOne), case1);
            vHL = _mm_blendv_epi8_rpl(vHL, vBias, cond_zero);
            _mm_store_si128(pvHLStore + i, vHL);

            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len, bias);
            arr_store(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len, bias);
            arr_store(result->stats->tables->length_table, vHL, i, segLen, j, s2Len, bias);
            arr_store(result->stats->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
            vMaxH = _mm_max_epi8_rpl(vH, vMaxH);
            vEF_opn = _mm_subs_epi8(vH, vGapO);

            /* Update vE value. */
            vE_ext = _mm_subs_epi8(vE, vGapE);
            vE = _mm_max_epi8_rpl(vEF_opn, vE_ext);
            case1 = _mm_cmpgt_epi8(vEF_opn, vE_ext);
            vEM = _mm_blendv_epi8_rpl(vEM, vHM, case1);
            vES = _mm_blendv_epi8_rpl(vES, vHS, case1);
            vEL = _mm_blendv_epi8_rpl(
                    _mm_adds_epi8(vEL, vOne),
                    _mm_adds_epi8(vHL, vOne),
                    case1);
            _mm_store_si128(pvE + i, vE);
            _mm_store_si128(pvEM + i, vEM);
            _mm_store_si128(pvES + i, vES);
            _mm_store_si128(pvEL + i, vEL);

            /* Update vF value. */
            vF_ext = _mm_subs_epi8(vF, vGapE);
            vF = _mm_max_epi8_rpl(vEF_opn, vF_ext);
            case1 = _mm_cmpgt_epi8(vEF_opn, vF_ext);
            vFM = _mm_blendv_epi8_rpl(vFM, vHM, case1);
            vFS = _mm_blendv_epi8_rpl(vFS, vHS, case1);
            vFL = _mm_blendv_epi8_rpl(
                    _mm_adds_epi8(vFL, vOne),
                    _mm_adds_epi8(vHL, vOne),
                    case1);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
            vHM = _mm_load_si128(pvHMLoad + i);
            vHS = _mm_load_si128(pvHSLoad + i);
            vHL = _mm_load_si128(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            __m128i vHp = _mm_load_si128(&pvHLoad[segLen - 1]);
            vHp = _mm_slli_si128(vHp, 1);
            vF = _mm_slli_si128(vF, 1);
            vFM = _mm_slli_si128(vFM, 1);
            vFS = _mm_slli_si128(vFS, 1);
            vFL = _mm_slli_si128(vFL, 1);
            vHp = _mm_blendv_epi8_rpl(vHp, vBias, insert_mask);
            vF = _mm_blendv_epi8_rpl(vF, vBias, insert_mask);
            vFM = _mm_blendv_epi8_rpl(vFM, vBias, insert_mask);
            vFS = _mm_blendv_epi8_rpl(vFS, vBias, insert_mask);
            vFL = _mm_blendv_epi8_rpl(vFL, vBias1, insert_mask);
            for (i=0; i<segLen; ++i) {
                __m128i case1;
                __m128i case2;
                __m128i cond;

                vHp = _mm_adds_epi8(vHp, _mm_load_si128(vP + i));
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi8_rpl(vH, vF);
                _mm_store_si128(pvHStore + i, vH);
                case1 = _mm_cmpeq_epi8(vH, vHp);
                case2 = _mm_cmpeq_epi8(vH, vF);
                cond = _mm_andnot_si128(case1, case2);

                /* calculate vM */
                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_blendv_epi8_rpl(vHM, vFM, cond);
                _mm_store_si128(pvHMStore + i, vHM);

                /* calculate vS */
                vHS = _mm_load_si128(pvHSStore + i);
                vHS = _mm_blendv_epi8_rpl(vHS, vFS, cond);
                _mm_store_si128(pvHSStore + i, vHS);

                /* calculate vL */
                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_blendv_epi8_rpl(vHL, vFL, cond);
                _mm_store_si128(pvHLStore + i, vHL);

                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHM);
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHS);
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len, bias);
                arr_store(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len, bias);
                arr_store(result->stats->tables->length_table, vHL, i, segLen, j, s2Len, bias);
                arr_store(result->stats->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
                vMaxH = _mm_max_epi8_rpl(vH, vMaxH);
                /* Update vF value. */
                vEF_opn = _mm_subs_epi8(vH, vGapO);
                vF_ext = _mm_subs_epi8(vF, vGapE);
                if (! _mm_movemask_epi8(
                            _mm_or_si128(
                                _mm_cmpgt_epi8(vF_ext, vEF_opn),
                                _mm_and_si128(
                                    _mm_cmpeq_epi8(vF_ext, vEF_opn),
                                    _mm_cmpgt_epi8(vF_ext, vBias)))))
                    goto end;
                /*vF = _mm_max_epi8_rpl(vEF_opn, vF_ext);*/
                vF = vF_ext;
                cond = _mm_cmpgt_epi8(vEF_opn, vF_ext);
                vFM = _mm_blendv_epi8_rpl(vFM, vHM, cond);
                vFS = _mm_blendv_epi8_rpl(vFS, vHS, cond);
                vFL = _mm_blendv_epi8_rpl(
                        _mm_adds_epi8(vFL, vOne),
                        _mm_adds_epi8(vHL, vOne),
                        cond);
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm_load_si128(pvHStore + offset);
            vHM = _mm_load_si128(pvHMStore + offset);
            vHS = _mm_load_si128(pvHSStore + offset);
            vHL = _mm_load_si128(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 1);
                vHM = _mm_slli_si128(vHM, 1);
                vHS = _mm_slli_si128(vHS, 1);
                vHL = _mm_slli_si128(vHL, 1);
            }
            result->stats->rowcols->score_row[j] = (int8_t) _mm_extract_epi8_rpl (vH, 15) - bias;
            result->stats->rowcols->matches_row[j] = (int8_t) _mm_extract_epi8_rpl (vHM, 15) - bias;
            result->stats->rowcols->similar_row[j] = (int8_t) _mm_extract_epi8_rpl (vHS, 15) - bias;
            result->stats->rowcols->length_row[j] = (int8_t) _mm_extract_epi8_rpl (vHL, 15) - bias;
        }
#endif

        {
            __m128i vCompare = _mm_cmpgt_epi8(vMaxH, vMaxHUnit);
            if (_mm_movemask_epi8(vCompare)) {
                score = _mm_hmax_epi8_rpl(vMaxH);
                /* if score has potential to overflow, abort early */
                if (score > maxp) {
                    result->flag |= PARASAIL_FLAG_SATURATED;
                    break;
                }
                vMaxHUnit = _mm_set1_epi8(score);
                end_ref = j;
            }
        }
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m128i vH = _mm_load_si128(pvHStore+i);
        __m128i vHM = _mm_load_si128(pvHMStore+i);
        __m128i vHS = _mm_load_si128(pvHSStore+i);
        __m128i vHL = _mm_load_si128(pvHLStore+i);
        arr_store_col(result->stats->rowcols->score_col, vH, i, segLen, bias);
        arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen, bias);
        arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen, bias);
        arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen, bias);
    }
#endif

    if (score == INT8_MAX
            || _mm_movemask_epi8(_mm_cmpeq_epi8(vSaturationCheckMax,vPosLimit))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
    }

    if (parasail_result_is_saturated(result)) {
        score = 0;
        end_query = 0;
        end_ref = 0;
        matches = 0;
        similar = 0;
        length = 0;
    }
    else {
        if (end_ref == j-1) {
            /* end_ref was the last store column */
            SWAP(pvHMax,  pvHStore)
            SWAP(pvHMMax, pvHMStore)
            SWAP(pvHSMax, pvHSStore)
            SWAP(pvHLMax, pvHLStore)
        }
        else if (end_ref == j-2) {
            /* end_ref was the last load column */
            SWAP(pvHMax,  pvHLoad)
            SWAP(pvHMMax, pvHMLoad)
            SWAP(pvHSMax, pvHSLoad)
            SWAP(pvHLMax, pvHLLoad)
        }
        /* Trace the alignment ending position on read. */
        {
            int8_t *t = (int8_t*)pvHMax;
            int8_t *m = (int8_t*)pvHMMax;
            int8_t *s = (int8_t*)pvHSMax;
            int8_t *l = (int8_t*)pvHLMax;
            int32_t column_len = segLen * segWidth;
            end_query = s1Len;
            for (i = 0; i<column_len; ++i, ++t, ++m, ++s, ++l) {
                if (*t == score) {
                    int32_t temp = i / segWidth + i % segWidth * segLen;
                    if (temp < end_query) {
                        end_query = temp;
                        matches = *m;
                        similar = *s;
                        length = *l;
                    }
                }
            }
        }
    }

    result->score = score - bias;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->stats->matches = matches - bias;
    result->stats->similar = similar - bias;
    result->stats->length = length - bias;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(pvHLMax);
    parasail_free(pvHSMax);
    parasail_free(pvHMMax);
    parasail_free(pvHMax);
    parasail_free(pvEL);
    parasail_free(pvES);
    parasail_free(pvEM);
    parasail_free(pvE);
    parasail_free(pvHLLoad);
    parasail_free(pvHLStore);
    parasail_free(pvHSLoad);
    parasail_free(pvHSStore);
    parasail_free(pvHMLoad);
    parasail_free(pvHMStore);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

#ifdef FASTSTATS
#ifdef PARASAIL_TABLE
#else
#ifdef PARASAIL_ROWCOL
#else
#include <assert.h>
parasail_result_t* INAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    const char *s1 = profile->s1;
    const parasail_matrix_t *matrix = profile->matrix;

    /* find the end loc first with the faster implementation */
    parasail_result_t *result = parasail_sw_striped_profile_sse2_128_8(profile, s2, s2Len, open, gap);
    if (!parasail_result_is_saturated(result)) {
#if 0
        int s1Len_new = 0;
        int s2Len_new = 0;
        char *s1_new = NULL;
        char *s2_new = NULL;
        parasail_profile_t *profile_new = NULL;
        parasail_result_t *result_new = NULL;
        int s1_begin = 0;
        int s2_begin = 0;
        int s1Len_final = 0;
        int s2Len_final = 0;
        parasail_profile_t *profile_final = NULL;
        parasail_result_t *result_final = NULL;

        /* using the end loc and the non-stats version of the function,
         * reverse the inputs and find the beg loc */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;
        s1_new = parasail_reverse(s1, s1Len_new);
        s2_new = parasail_reverse(s2, s2Len_new);
        profile_new = parasail_profile_create_sse_128_8(
                s1_new, s1Len_new, matrix);
        profile_new->stop = result->score;
        result_new = parasail_sw_striped_profile_sse2_128_8(
                profile_new, s2_new, s2Len_new, open, gap);

        /* using both the beg and end loc, call the original stats func */
        s1_begin = s1Len_new - result_new->end_query - 1;
        s2_begin = s2Len_new - result_new->end_ref - 1;
        s1Len_final = s1Len_new - s1_begin;
        s2Len_final = s2Len_new - s2_begin;
        assert(s1_begin >= 0);
        assert(s2_begin >= 0);
        assert(s1Len_new > s1_begin);
        assert(s2Len_new > s2_begin);
        profile_final = parasail_profile_create_stats_sse_128_8(
                &s1[s1_begin], s1Len_final, matrix);
        result_final = PNAME(
                profile_final, &s2[s2_begin], s2Len_final, open, gap);

        /* clean up all the temporary profiles, sequences, and results */
        free(s1_new);
        free(s2_new);
        parasail_profile_free(profile_new);
        parasail_profile_free(profile_final);
        parasail_result_free(result);
        parasail_result_free(result_new);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
#else
        int s1Len_new = 0;
        int s2Len_new = 0;
        parasail_profile_t *profile_final = NULL;
        parasail_result_t *result_final = NULL;

        /* using the end loc, call the original stats function */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;
        profile_final = parasail_profile_create_stats_sse_128_8(
                s1, s1Len_new, matrix);
        result_final = PNAME(
                profile_final, s2, s2Len_new, open, gap);

        /* clean up all the temporary profiles, sequences, and results */
        parasail_profile_free(profile_final);
        parasail_result_free(result);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
#endif
    }
    else {
        return result;
    }
}
#endif
#endif
#endif


