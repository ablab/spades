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


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*( 0*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  0);
    array[1LL*( 1*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  1);
    array[1LL*( 2*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  2);
    array[1LL*( 3*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  3);
    array[1LL*( 4*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  4);
    array[1LL*( 5*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  5);
    array[1LL*( 6*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  6);
    array[1LL*( 7*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  7);
    array[1LL*( 8*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  8);
    array[1LL*( 9*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  9);
    array[1LL*(10*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 10);
    array[1LL*(11*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 11);
    array[1LL*(12*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 12);
    array[1LL*(13*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 13);
    array[1LL*(14*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 14);
    array[1LL*(15*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 15);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  0);
    col[ 1*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  1);
    col[ 2*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  2);
    col[ 3*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  3);
    col[ 4*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  4);
    col[ 5*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  5);
    col[ 6*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  6);
    col[ 7*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  7);
    col[ 8*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  8);
    col[ 9*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  9);
    col[10*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 10);
    col[11*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 11);
    col[12*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 12);
    col[13*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 13);
    col[14*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 14);
    col[15*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_stats_table_striped_avx2_256_16
#define PNAME parasail_nw_stats_table_striped_profile_avx2_256_16
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_stats_rowcol_striped_avx2_256_16
#define PNAME parasail_nw_stats_rowcol_striped_profile_avx2_256_16
#else
#define FNAME parasail_nw_stats_striped_avx2_256_16
#define PNAME parasail_nw_stats_striped_profile_avx2_256_16
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_avx_256_16(s1, s1Len, matrix);
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
    __m256i* const restrict vProfile  = (__m256i*)profile->profile16.score;
    __m256i* const restrict vProfileM = (__m256i*)profile->profile16.matches;
    __m256i* const restrict vProfileS = (__m256i*)profile->profile16.similar;
    __m256i* restrict pvHStore        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLoad         = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMStore       = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMLoad        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHSStore       = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHSLoad        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLStore       = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLLoad        = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvE             = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEM      = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvES      = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEL      = parasail_memalign___m256i(32, segLen);
    int16_t* const restrict boundary  = parasail_memalign_int16_t(32, s2Len+1);
    const __m256i vGapO = _mm256_set1_epi16(open);
    const __m256i vGapE = _mm256_set1_epi16(gap);
    const int16_t NEG_LIMIT = (-open < matrix->min ?
        INT16_MIN + open : INT16_MIN - matrix->min) + 1;
    const int16_t POS_LIMIT = INT16_MAX - matrix->max - 1;
    const __m256i vZero = _mm256_setzero_si256();
    const __m256i vOne = _mm256_set1_epi16(1);
    int16_t score = NEG_LIMIT;
    int16_t matches = 0;
    int16_t similar = 0;
    int16_t length = 0;
    __m256i vNegLimit = _mm256_set1_epi16(NEG_LIMIT);
    __m256i vPosLimit = _mm256_set1_epi16(POS_LIMIT);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    parasail_memset___m256i(pvHMStore, vZero, segLen);
    parasail_memset___m256i(pvHSStore, vZero, segLen);
    parasail_memset___m256i(pvHLStore, vZero, segLen);
    parasail_memset___m256i(pvEM, vZero, segLen);
    parasail_memset___m256i(pvES, vZero, segLen);
    parasail_memset___m256i(pvEL, vOne, segLen);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            __m256i_16_t h;
            __m256i_16_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT16_MIN ? INT16_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT16_MIN ? INT16_MIN : tmp;
            }
            _mm256_store_si256(&pvHStore[index], h.m);
            _mm256_store_si256(&pvE[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT16_MIN ? INT16_MIN : tmp;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vEF_opn;
        __m256i vE;
        __m256i vE_ext;
        __m256i vEM;
        __m256i vES;
        __m256i vEL;
        __m256i vF;
        __m256i vF_ext;
        __m256i vFM;
        __m256i vFS;
        __m256i vFL;
        __m256i vH;
        __m256i vH_dag;
        __m256i vHM;
        __m256i vHS;
        __m256i vHL;
        const __m256i* vP = NULL;
        const __m256i* vPM = NULL;
        const __m256i* vPS = NULL;

        /* Initialize F value to neg inf.  Any errors to vH values will
         * be corrected in the Lazy_F loop.  */
        vF = vNegLimit;
        vFM = vZero;
        vFS = vZero;
        vFL = vOne;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm256_load_si256(&pvHStore[segLen - 1]);
        vHM = _mm256_load_si256(&pvHMStore[segLen - 1]);
        vHS = _mm256_load_si256(&pvHSStore[segLen - 1]);
        vHL = _mm256_load_si256(&pvHLStore[segLen - 1]);
        vH = _mm256_slli_si256_rpl(vH, 2);
        vHM = _mm256_slli_si256_rpl(vHM, 2);
        vHS = _mm256_slli_si256_rpl(vHS, 2);
        vHL = _mm256_slli_si256_rpl(vHL, 2);

        /* insert upper boundary condition */
        vH = _mm256_insert_epi16_rpl(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPM = vProfileM + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad,  pvHStore)
        SWAP(pvHMLoad, pvHMStore)
        SWAP(pvHSLoad, pvHSStore)
        SWAP(pvHLLoad, pvHLStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            __m256i case1;
            __m256i case2;

            vE = _mm256_load_si256(pvE+ i);
            vEM = _mm256_load_si256(pvEM+ i);
            vES = _mm256_load_si256(pvES+ i);
            vEL = _mm256_load_si256(pvEL+ i);

            /* Get max from vH, vE and vF. */
            vH_dag = _mm256_add_epi16(vH, _mm256_load_si256(vP + i));
            vH = _mm256_max_epi16(vH_dag, vE);
            vH = _mm256_max_epi16(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);

            case1 = _mm256_cmpeq_epi16(vH, vH_dag);
            case2 = _mm256_cmpeq_epi16(vH, vF);

            /* calculate vM */
            vHM = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEM, vFM, case2),
                    _mm256_add_epi16(vHM, _mm256_load_si256(vPM + i)),
                    case1);
            _mm256_store_si256(pvHMStore + i, vHM);

            /* calculate vS */
            vHS = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vES, vFS, case2),
                    _mm256_add_epi16(vHS, _mm256_load_si256(vPS + i)),
                    case1);
            _mm256_store_si256(pvHSStore + i, vHS);

            /* calculate vL */
            vHL = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEL, vFL, case2),
                    _mm256_add_epi16(vHL, vOne),
                    case1);
            _mm256_store_si256(pvHLStore + i, vHL);

            vSaturationCheckMin = _mm256_min_epi16(vSaturationCheckMin, vH);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vH);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si256(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si256(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si256(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
            vEF_opn = _mm256_sub_epi16(vH, vGapO);

            /* Update vE value. */
            vE_ext = _mm256_sub_epi16(vE, vGapE);
            vE = _mm256_max_epi16(vEF_opn, vE_ext);
            case1 = _mm256_cmpgt_epi16(vEF_opn, vE_ext);
            vEM = _mm256_blendv_epi8(vEM, vHM, case1);
            vES = _mm256_blendv_epi8(vES, vHS, case1);
            vEL = _mm256_blendv_epi8(
                    _mm256_add_epi16(vEL, vOne),
                    _mm256_add_epi16(vHL, vOne),
                    case1);
            _mm256_store_si256(pvE + i, vE);
            _mm256_store_si256(pvEM + i, vEM);
            _mm256_store_si256(pvES + i, vES);
            _mm256_store_si256(pvEL + i, vEL);

            /* Update vF value. */
            vF_ext = _mm256_sub_epi16(vF, vGapE);
            vF = _mm256_max_epi16(vEF_opn, vF_ext);
            case1 = _mm256_cmpgt_epi16(vEF_opn, vF_ext);
            vFM = _mm256_blendv_epi8(vFM, vHM, case1);
            vFS = _mm256_blendv_epi8(vFS, vHS, case1);
            vFL = _mm256_blendv_epi8(
                    _mm256_add_epi16(vFL, vOne),
                    _mm256_add_epi16(vHL, vOne),
                    case1);

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
            vHM = _mm256_load_si256(pvHMLoad + i);
            vHS = _mm256_load_si256(pvHSLoad + i);
            vHL = _mm256_load_si256(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            __m256i vHp = _mm256_load_si256(&pvHLoad[segLen - 1]);
            int64_t tmp = boundary[j+1]-open;
            int16_t tmp2 = tmp < INT16_MIN ? INT16_MIN : tmp;
            vHp = _mm256_slli_si256_rpl(vHp, 2);
            vF = _mm256_slli_si256_rpl(vF, 2);
            vFM = _mm256_slli_si256_rpl(vFM, 2);
            vFS = _mm256_slli_si256_rpl(vFS, 2);
            vFL = _mm256_slli_si256_rpl(vFL, 2);
            vHp = _mm256_insert_epi16_rpl(vHp, boundary[j], 0);
            vF = _mm256_insert_epi16_rpl(vF, tmp2, 0);
            vFL = _mm256_insert_epi16_rpl(vFL, 1, 0);
            for (i=0; i<segLen; ++i) {
                __m256i case1;
                __m256i case2;
                __m256i cond;

                vHp = _mm256_add_epi16(vHp, _mm256_load_si256(vP + i));
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi16(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                case1 = _mm256_cmpeq_epi16(vH, vHp);
                case2 = _mm256_cmpeq_epi16(vH, vF);
                cond = _mm256_andnot_si256(case1, case2);

                /* calculate vM */
                vHM = _mm256_load_si256(pvHMStore + i);
                vHM = _mm256_blendv_epi8(vHM, vFM, cond);
                _mm256_store_si256(pvHMStore + i, vHM);

                /* calculate vS */
                vHS = _mm256_load_si256(pvHSStore + i);
                vHS = _mm256_blendv_epi8(vHS, vFS, cond);
                _mm256_store_si256(pvHSStore + i, vHS);

                /* calculate vL */
                vHL = _mm256_load_si256(pvHLStore + i);
                vHL = _mm256_blendv_epi8(vHL, vFL, cond);
                _mm256_store_si256(pvHLStore + i, vHL);

                vSaturationCheckMin = _mm256_min_epi16(vSaturationCheckMin, vH);
                vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vH);
                vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHM);
                vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHS);
                vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si256(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si256(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si256(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                /* Update vF value. */
                vEF_opn = _mm256_sub_epi16(vH, vGapO);
                vF_ext = _mm256_sub_epi16(vF, vGapE);
                if (! _mm256_movemask_epi8(
                            _mm256_or_si256(
                                _mm256_cmpgt_epi16(vF_ext, vEF_opn),
                                _mm256_cmpeq_epi16(vF_ext, vEF_opn))))
                    goto end;
                /*vF = _mm256_max_epi16(vEF_opn, vF_ext);*/
                vF = vF_ext;
                cond = _mm256_cmpgt_epi16(vEF_opn, vF_ext);
                vFM = _mm256_blendv_epi8(vFM, vHM, cond);
                vFS = _mm256_blendv_epi8(vFS, vHS, cond);
                vFL = _mm256_blendv_epi8(
                        _mm256_add_epi16(vFL, vOne),
                        _mm256_add_epi16(vHL, vOne),
                        cond);
                vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm256_load_si256(pvHStore + offset);
            vHM = _mm256_load_si256(pvHMStore + offset);
            vHS = _mm256_load_si256(pvHSStore + offset);
            vHL = _mm256_load_si256(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl (vH, 2);
                vHM = _mm256_slli_si256_rpl (vHM, 2);
                vHS = _mm256_slli_si256_rpl (vHS, 2);
                vHL = _mm256_slli_si256_rpl (vHL, 2);
            }
            result->stats->rowcols->score_row[j] = (int16_t) _mm256_extract_epi16_rpl (vH, 15);
            result->stats->rowcols->matches_row[j] = (int16_t) _mm256_extract_epi16_rpl (vHM, 15);
            result->stats->rowcols->similar_row[j] = (int16_t) _mm256_extract_epi16_rpl (vHS, 15);
            result->stats->rowcols->length_row[j] = (int16_t) _mm256_extract_epi16_rpl (vHL, 15);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m256i vH = _mm256_load_si256(pvHStore+i);
        __m256i vHM = _mm256_load_si256(pvHMStore+i);
        __m256i vHS = _mm256_load_si256(pvHSStore+i);
        __m256i vHL = _mm256_load_si256(pvHLStore+i);
        arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
        arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
        arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
        arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        __m256i vH = _mm256_load_si256(pvHStore + offset);
        __m256i vHM = _mm256_load_si256(pvHMStore + offset);
        __m256i vHS = _mm256_load_si256(pvHSStore + offset);
        __m256i vHL = _mm256_load_si256(pvHLStore + offset);
        for (k=0; k<position; ++k) {
            vH = _mm256_slli_si256_rpl (vH, 2);
            vHM = _mm256_slli_si256_rpl (vHM, 2);
            vHS = _mm256_slli_si256_rpl (vHS, 2);
            vHL = _mm256_slli_si256_rpl (vHL, 2);
        }
        score = (int16_t) _mm256_extract_epi16_rpl (vH, 15);
        matches = (int16_t) _mm256_extract_epi16_rpl (vHM, 15);
        similar = (int16_t) _mm256_extract_epi16_rpl (vHS, 15);
        length = (int16_t) _mm256_extract_epi16_rpl (vHL, 15);
    }

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmplt_epi16_rpl(vSaturationCheckMin, vNegLimit),
            _mm256_cmpgt_epi16(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        matches = 0;
        similar = 0;
        length = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->stats->matches = matches;
    result->stats->similar = similar;
    result->stats->length = length;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_16;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(boundary);
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


