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



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        vec128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*(0*seglen+t)*dlen + d] = (int32_t)_mm_extract_epi32(vH, 0);
    array[1LL*(1*seglen+t)*dlen + d] = (int32_t)_mm_extract_epi32(vH, 1);
    array[1LL*(2*seglen+t)*dlen + d] = (int32_t)_mm_extract_epi32(vH, 2);
    array[1LL*(3*seglen+t)*dlen + d] = (int32_t)_mm_extract_epi32(vH, 3);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        vec128i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int32_t)_mm_extract_epi32(vH, 0);
    col[1*seglen+t] = (int32_t)_mm_extract_epi32(vH, 1);
    col[2*seglen+t] = (int32_t)_mm_extract_epi32(vH, 2);
    col[3*seglen+t] = (int32_t)_mm_extract_epi32(vH, 3);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_stats_table_striped_altivec_128_32
#define PNAME parasail_nw_stats_table_striped_profile_altivec_128_32
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_stats_rowcol_striped_altivec_128_32
#define PNAME parasail_nw_stats_rowcol_striped_profile_altivec_128_32
#else
#define FNAME parasail_nw_stats_striped_altivec_128_32
#define PNAME parasail_nw_stats_striped_profile_altivec_128_32
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_altivec_128_32(s1, s1Len, matrix);
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
    vec128i* const restrict vProfile  = (vec128i*)profile->profile32.score;
    vec128i* const restrict vProfileM = (vec128i*)profile->profile32.matches;
    vec128i* const restrict vProfileS = (vec128i*)profile->profile32.similar;
    vec128i* restrict pvHStore        = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvHLoad         = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvHMStore       = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvHMLoad        = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvHSStore       = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvHSLoad        = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvHLStore       = parasail_memalign_vec128i(16, segLen);
    vec128i* restrict pvHLLoad        = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvE             = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvEM      = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvES      = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvEL      = parasail_memalign_vec128i(16, segLen);
    int32_t* const restrict boundary  = parasail_memalign_int32_t(16, s2Len+1);
    const vec128i vGapO = _mm_set1_epi32(open);
    const vec128i vGapE = _mm_set1_epi32(gap);
    const int32_t NEG_LIMIT = (-open < matrix->min ?
        INT32_MIN + open : INT32_MIN - matrix->min) + 1;
    const int32_t POS_LIMIT = INT32_MAX - matrix->max - 1;
    const vec128i vZero = _mm_setzero_si128();
    const vec128i vOne = _mm_set1_epi32(1);
    int32_t score = NEG_LIMIT;
    int32_t matches = 0;
    int32_t similar = 0;
    int32_t length = 0;
    vec128i vNegLimit = _mm_set1_epi32(NEG_LIMIT);
    vec128i vPosLimit = _mm_set1_epi32(POS_LIMIT);
    vec128i vSaturationCheckMin = vPosLimit;
    vec128i vSaturationCheckMax = vNegLimit;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    parasail_memset_vec128i(pvHMStore, vZero, segLen);
    parasail_memset_vec128i(pvHSStore, vZero, segLen);
    parasail_memset_vec128i(pvHLStore, vZero, segLen);
    parasail_memset_vec128i(pvEM, vZero, segLen);
    parasail_memset_vec128i(pvES, vZero, segLen);
    parasail_memset_vec128i(pvEL, vOne, segLen);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            vec128i_32_t h;
            vec128i_32_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT32_MIN ? INT32_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT32_MIN ? INT32_MIN : tmp;
            }
            _mm_store_si128(&pvHStore[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT32_MIN ? INT32_MIN : tmp;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        vec128i vEF_opn;
        vec128i vE;
        vec128i vE_ext;
        vec128i vEM;
        vec128i vES;
        vec128i vEL;
        vec128i vF;
        vec128i vF_ext;
        vec128i vFM;
        vec128i vFS;
        vec128i vFL;
        vec128i vH;
        vec128i vH_dag;
        vec128i vHM;
        vec128i vHS;
        vec128i vHL;
        const vec128i* vP = NULL;
        const vec128i* vPM = NULL;
        const vec128i* vPS = NULL;

        /* Initialize F value to neg inf.  Any errors to vH values will
         * be corrected in the Lazy_F loop.  */
        vF = vNegLimit;
        vFM = vZero;
        vFS = vZero;
        vFL = vOne;

        /* load final segment of pvHStore and shift left by 4 bytes */
        vH = _mm_load_si128(&pvHStore[segLen - 1]);
        vHM = _mm_load_si128(&pvHMStore[segLen - 1]);
        vHS = _mm_load_si128(&pvHSStore[segLen - 1]);
        vHL = _mm_load_si128(&pvHLStore[segLen - 1]);
        vH = _mm_slli_si128(vH, 4);
        vHM = _mm_slli_si128(vHM, 4);
        vHS = _mm_slli_si128(vHS, 4);
        vHL = _mm_slli_si128(vHL, 4);

        /* insert upper boundary condition */
        vH = _mm_insert_epi32(vH, boundary[j], 0);

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
            vec128i case1;
            vec128i case2;

            vE = _mm_load_si128(pvE+ i);
            vEM = _mm_load_si128(pvEM+ i);
            vES = _mm_load_si128(pvES+ i);
            vEL = _mm_load_si128(pvEL+ i);

            /* Get max from vH, vE and vF. */
            vH_dag = _mm_add_epi32(vH, _mm_load_si128(vP + i));
            vH = _mm_max_epi32(vH_dag, vE);
            vH = _mm_max_epi32(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            case1 = _mm_cmpeq_epi32(vH, vH_dag);
            case2 = _mm_cmpeq_epi32(vH, vF);

            /* calculate vM */
            vHM = _mm_blendv_epi8(
                    _mm_blendv_epi8(vEM, vFM, case2),
                    _mm_add_epi32(vHM, _mm_load_si128(vPM + i)),
                    case1);
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vS */
            vHS = _mm_blendv_epi8(
                    _mm_blendv_epi8(vES, vFS, case2),
                    _mm_add_epi32(vHS, _mm_load_si128(vPS + i)),
                    case1);
            _mm_store_si128(pvHSStore + i, vHS);

            /* calculate vL */
            vHL = _mm_blendv_epi8(
                    _mm_blendv_epi8(vEL, vFL, case2),
                    _mm_add_epi32(vHL, vOne),
                    case1);
            _mm_store_si128(pvHLStore + i, vHL);

            vSaturationCheckMin = _mm_min_epi32(vSaturationCheckMin, vH);
            vSaturationCheckMax = _mm_max_epi32(vSaturationCheckMax, vH);
            vSaturationCheckMax = _mm_max_epi32(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm_max_epi32(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm_max_epi32(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
            vEF_opn = _mm_sub_epi32(vH, vGapO);

            /* Update vE value. */
            vE_ext = _mm_sub_epi32(vE, vGapE);
            vE = _mm_max_epi32(vEF_opn, vE_ext);
            case1 = _mm_cmpgt_epi32(vEF_opn, vE_ext);
            vEM = _mm_blendv_epi8(vEM, vHM, case1);
            vES = _mm_blendv_epi8(vES, vHS, case1);
            vEL = _mm_blendv_epi8(
                    _mm_add_epi32(vEL, vOne),
                    _mm_add_epi32(vHL, vOne),
                    case1);
            _mm_store_si128(pvE + i, vE);
            _mm_store_si128(pvEM + i, vEM);
            _mm_store_si128(pvES + i, vES);
            _mm_store_si128(pvEL + i, vEL);

            /* Update vF value. */
            vF_ext = _mm_sub_epi32(vF, vGapE);
            vF = _mm_max_epi32(vEF_opn, vF_ext);
            case1 = _mm_cmpgt_epi32(vEF_opn, vF_ext);
            vFM = _mm_blendv_epi8(vFM, vHM, case1);
            vFS = _mm_blendv_epi8(vFS, vHS, case1);
            vFL = _mm_blendv_epi8(
                    _mm_add_epi32(vFL, vOne),
                    _mm_add_epi32(vHL, vOne),
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
            vec128i vHp = _mm_load_si128(&pvHLoad[segLen - 1]);
            int64_t tmp = boundary[j+1]-open;
            int32_t tmp2 = tmp < INT32_MIN ? INT32_MIN : tmp;
            vHp = _mm_slli_si128(vHp, 4);
            vF = _mm_slli_si128(vF, 4);
            vFM = _mm_slli_si128(vFM, 4);
            vFS = _mm_slli_si128(vFS, 4);
            vFL = _mm_slli_si128(vFL, 4);
            vHp = _mm_insert_epi32(vHp, boundary[j], 0);
            vF = _mm_insert_epi32(vF, tmp2, 0);
            vFL = _mm_insert_epi32(vFL, 1, 0);
            for (i=0; i<segLen; ++i) {
                vec128i case1;
                vec128i case2;
                vec128i cond;

                vHp = _mm_add_epi32(vHp, _mm_load_si128(vP + i));
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi32(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
                case1 = _mm_cmpeq_epi32(vH, vHp);
                case2 = _mm_cmpeq_epi32(vH, vF);
                cond = _mm_andnot_si128(case1, case2);

                /* calculate vM */
                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_blendv_epi8(vHM, vFM, cond);
                _mm_store_si128(pvHMStore + i, vHM);

                /* calculate vS */
                vHS = _mm_load_si128(pvHSStore + i);
                vHS = _mm_blendv_epi8(vHS, vFS, cond);
                _mm_store_si128(pvHSStore + i, vHS);

                /* calculate vL */
                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_blendv_epi8(vHL, vFL, cond);
                _mm_store_si128(pvHLStore + i, vHL);

                vSaturationCheckMin = _mm_min_epi32(vSaturationCheckMin, vH);
                vSaturationCheckMax = _mm_max_epi32(vSaturationCheckMax, vH);
                vSaturationCheckMax = _mm_max_epi32(vSaturationCheckMax, vHM);
                vSaturationCheckMax = _mm_max_epi32(vSaturationCheckMax, vHS);
                vSaturationCheckMax = _mm_max_epi32(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si128(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si128(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si128(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                /* Update vF value. */
                vEF_opn = _mm_sub_epi32(vH, vGapO);
                vF_ext = _mm_sub_epi32(vF, vGapE);
                if (! _mm_movemask_epi8(
                            _mm_or_si128(
                                _mm_cmpgt_epi32(vF_ext, vEF_opn),
                                _mm_cmpeq_epi32(vF_ext, vEF_opn))))
                    goto end;
                /*vF = _mm_max_epi32(vEF_opn, vF_ext);*/
                vF = vF_ext;
                cond = _mm_cmpgt_epi32(vEF_opn, vF_ext);
                vFM = _mm_blendv_epi8(vFM, vHM, cond);
                vFS = _mm_blendv_epi8(vFS, vHS, cond);
                vFL = _mm_blendv_epi8(
                        _mm_add_epi32(vFL, vOne),
                        _mm_add_epi32(vHL, vOne),
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
                vH = _mm_slli_si128 (vH, 4);
                vHM = _mm_slli_si128 (vHM, 4);
                vHS = _mm_slli_si128 (vHS, 4);
                vHL = _mm_slli_si128 (vHL, 4);
            }
            result->stats->rowcols->score_row[j] = (int32_t) _mm_extract_epi32 (vH, 3);
            result->stats->rowcols->matches_row[j] = (int32_t) _mm_extract_epi32 (vHM, 3);
            result->stats->rowcols->similar_row[j] = (int32_t) _mm_extract_epi32 (vHS, 3);
            result->stats->rowcols->length_row[j] = (int32_t) _mm_extract_epi32 (vHL, 3);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        vec128i vH = _mm_load_si128(pvHStore+i);
        vec128i vHM = _mm_load_si128(pvHMStore+i);
        vec128i vHS = _mm_load_si128(pvHSStore+i);
        vec128i vHL = _mm_load_si128(pvHLStore+i);
        arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
        arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
        arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
        arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        vec128i vH = _mm_load_si128(pvHStore + offset);
        vec128i vHM = _mm_load_si128(pvHMStore + offset);
        vec128i vHS = _mm_load_si128(pvHSStore + offset);
        vec128i vHL = _mm_load_si128(pvHLStore + offset);
        for (k=0; k<position; ++k) {
            vH = _mm_slli_si128 (vH, 4);
            vHM = _mm_slli_si128 (vHM, 4);
            vHS = _mm_slli_si128 (vHS, 4);
            vHL = _mm_slli_si128 (vHL, 4);
        }
        score = (int32_t) _mm_extract_epi32 (vH, 3);
        matches = (int32_t) _mm_extract_epi32 (vHM, 3);
        similar = (int32_t) _mm_extract_epi32 (vHS, 3);
        length = (int32_t) _mm_extract_epi32 (vHL, 3);
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmplt_epi32(vSaturationCheckMin, vNegLimit),
            _mm_cmpgt_epi32(vSaturationCheckMax, vPosLimit)))) {
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
        | PARASAIL_FLAG_BITS_32 | PARASAIL_FLAG_LANES_4;
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


