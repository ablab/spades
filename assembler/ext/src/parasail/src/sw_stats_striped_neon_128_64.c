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

#define FASTSTATS

#define SWAP(A,B) { simde__m128i* tmp = A; A = B; B = tmp; }
#define SWAP3(A,B,C) { simde__m128i* tmp = A; A = B; B = C; C = tmp; }

#define NEG_INF (INT64_MIN/(int64_t)(2))


#ifdef PARASAIL_TABLE
static inline void arr_store(
        int *array,
        simde__m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*(0*seglen+t)*dlen + d] = (int64_t)simde_mm_extract_epi64(vH, 0);
    array[1LL*(1*seglen+t)*dlen + d] = (int64_t)simde_mm_extract_epi64(vH, 1);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        simde__m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int64_t)simde_mm_extract_epi64(vH, 0);
    col[1*seglen+t] = (int64_t)simde_mm_extract_epi64(vH, 1);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_striped_neon_128_64
#define PNAME parasail_sw_stats_table_striped_profile_neon_128_64
#define INAME PNAME
#define STATIC
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_striped_neon_128_64
#define PNAME parasail_sw_stats_rowcol_striped_profile_neon_128_64
#define INAME PNAME
#define STATIC
#else
#define FNAME parasail_sw_stats_striped_neon_128_64
#ifdef FASTSTATS
#define PNAME parasail_sw_stats_striped_profile_neon_128_64_internal
#define INAME parasail_sw_stats_striped_profile_neon_128_64
#define STATIC static
#else
#define PNAME parasail_sw_stats_striped_profile_neon_128_64
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
    parasail_profile_t *profile = parasail_profile_create_stats_neon_128_64(s1, s1Len, matrix);
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
    const int32_t segWidth = 2; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile  = (simde__m128i*)profile->profile64.score;
    simde__m128i* const restrict vProfileM = (simde__m128i*)profile->profile64.matches;
    simde__m128i* const restrict vProfileS = (simde__m128i*)profile->profile64.similar;
    simde__m128i* restrict pvHStore        = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHLoad         = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHMStore       = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHMLoad        = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHSStore       = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHSLoad        = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHLStore       = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHLLoad        = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvE       = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvEM      = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvES      = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvEL      = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHMax          = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHMMax          = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHSMax          = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHLMax          = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i vGapO = simde_mm_set1_epi64x(open);
    simde__m128i vGapE = simde_mm_set1_epi64x(gap);
    simde__m128i vZero = simde_mm_setzero_si128();
    simde__m128i vOne = simde_mm_set1_epi64x(1);
    int64_t score = NEG_INF;
    int64_t matches = NEG_INF;
    int64_t similar = NEG_INF;
    int64_t length = NEG_INF;
    simde__m128i vMaxH = vZero;
    simde__m128i vMaxHUnit = vZero;
    simde__m128i vSaturationCheckMax = vZero;
    simde__m128i vPosLimit = simde_mm_set1_epi64x(INT64_MAX);
    int64_t maxp = INT64_MAX - (int64_t)(matrix->max+1);
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

    parasail_memset_simde__m128i(pvHStore, vZero, segLen);
    parasail_memset_simde__m128i(pvHMStore, vZero, segLen);
    parasail_memset_simde__m128i(pvHSStore, vZero, segLen);
    parasail_memset_simde__m128i(pvHLStore, vZero, segLen);
    parasail_memset_simde__m128i(pvE, simde_mm_set1_epi64x(-open), segLen);
    parasail_memset_simde__m128i(pvEM, vZero, segLen);
    parasail_memset_simde__m128i(pvES, vZero, segLen);
    parasail_memset_simde__m128i(pvEL, vOne, segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vEF_opn;
        simde__m128i vE;
        simde__m128i vE_ext;
        simde__m128i vEM;
        simde__m128i vES;
        simde__m128i vEL;
        simde__m128i vF;
        simde__m128i vF_ext;
        simde__m128i vFM;
        simde__m128i vFS;
        simde__m128i vFL;
        simde__m128i vH;
        simde__m128i vH_dag;
        simde__m128i vHM;
        simde__m128i vHS;
        simde__m128i vHL;
        const simde__m128i* vP = NULL;
        const simde__m128i* vPM = NULL;
        const simde__m128i* vPS = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = vZero;
        vFM = vZero;
        vFS = vZero;
        vFL = vOne;

        /* load final segment of pvHStore and shift left by 8 bytes */
        vH = simde_mm_load_si128(&pvHStore[segLen - 1]);
        vHM = simde_mm_load_si128(&pvHMStore[segLen - 1]);
        vHS = simde_mm_load_si128(&pvHSStore[segLen - 1]);
        vHL = simde_mm_load_si128(&pvHLStore[segLen - 1]);
        vH = simde_mm_slli_si128(vH, 8);
        vHM = simde_mm_slli_si128(vHM, 8);
        vHS = simde_mm_slli_si128(vHS, 8);
        vHL = simde_mm_slli_si128(vHL, 8);

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
            simde__m128i cond_zero;
            simde__m128i case1;
            simde__m128i case2;

            vE = simde_mm_load_si128(pvE+ i);
            vEM = simde_mm_load_si128(pvEM+ i);
            vES = simde_mm_load_si128(pvES+ i);
            vEL = simde_mm_load_si128(pvEL+ i);

            /* Get max from vH, vE and vF. */
            vH_dag = simde_mm_add_epi64(vH, simde_mm_load_si128(vP + i));
            vH_dag = simde_mm_max_epi64(vH_dag, vZero);
            vH = simde_mm_max_epi64(vH_dag, vE);
            vH = simde_mm_max_epi64(vH, vF);
            /* Save vH values. */
            simde_mm_store_si128(pvHStore + i, vH);
            cond_zero = simde_mm_cmpeq_epi64(vH, vZero);

            case1 = simde_mm_cmpeq_epi64(vH, vH_dag);
            case2 = simde_mm_cmpeq_epi64(vH, vF);

            /* calculate vM */
            vHM = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEM, vFM, case2),
                    simde_mm_add_epi64(vHM, simde_mm_load_si128(vPM + i)), case1);
            vHM = simde_mm_andnot_si128(cond_zero, vHM);
            simde_mm_store_si128(pvHMStore + i, vHM);

            /* calculate vS */
            vHS = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vES, vFS, case2),
                    simde_mm_add_epi64(vHS, simde_mm_load_si128(vPS + i)), case1);
            vHS = simde_mm_andnot_si128(cond_zero, vHS);
            simde_mm_store_si128(pvHSStore + i, vHS);

            /* calculate vL */
            vHL = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEL, vFL, case2),
                    simde_mm_add_epi64(vHL, vOne), case1);
            vHL = simde_mm_andnot_si128(cond_zero, vHL);
            simde_mm_store_si128(pvHLStore + i, vHL);

            vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHM);
            vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHS);
            vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
            arr_store(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
            vMaxH = simde_mm_max_epi64(vH, vMaxH);
            vEF_opn = simde_mm_sub_epi64(vH, vGapO);

            /* Update vE value. */
            vE_ext = simde_mm_sub_epi64(vE, vGapE);
            vE = simde_mm_max_epi64(vEF_opn, vE_ext);
            case1 = simde_mm_cmpgt_epi64(vEF_opn, vE_ext);
            vEM = simde_mm_blendv_epi8(vEM, vHM, case1);
            vES = simde_mm_blendv_epi8(vES, vHS, case1);
            vEL = simde_mm_blendv_epi8(
                    simde_mm_add_epi64(vEL, vOne),
                    simde_mm_add_epi64(vHL, vOne),
                    case1);
            simde_mm_store_si128(pvE + i, vE);
            simde_mm_store_si128(pvEM + i, vEM);
            simde_mm_store_si128(pvES + i, vES);
            simde_mm_store_si128(pvEL + i, vEL);

            /* Update vF value. */
            vF_ext = simde_mm_sub_epi64(vF, vGapE);
            vF = simde_mm_max_epi64(vEF_opn, vF_ext);
            case1 = simde_mm_cmpgt_epi64(vEF_opn, vF_ext);
            vFM = simde_mm_blendv_epi8(vFM, vHM, case1);
            vFS = simde_mm_blendv_epi8(vFS, vHS, case1);
            vFL = simde_mm_blendv_epi8(
                    simde_mm_add_epi64(vFL, vOne),
                    simde_mm_add_epi64(vHL, vOne),
                    case1);

            /* Load the next vH. */
            vH = simde_mm_load_si128(pvHLoad + i);
            vHM = simde_mm_load_si128(pvHMLoad + i);
            vHS = simde_mm_load_si128(pvHSLoad + i);
            vHL = simde_mm_load_si128(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            simde__m128i vHp = simde_mm_load_si128(&pvHLoad[segLen - 1]);
            vHp = simde_mm_slli_si128(vHp, 8);
            vF = simde_mm_slli_si128(vF, 8);
            vF = simde_mm_insert_epi64(vF, -open, 0);
            vFM = simde_mm_slli_si128(vFM, 8);
            vFS = simde_mm_slli_si128(vFS, 8);
            vFL = simde_mm_slli_si128(vFL, 8);
            vFL = simde_mm_insert_epi64(vFL, 1, 0);
            for (i=0; i<segLen; ++i) {
                simde__m128i case1;
                simde__m128i case2;
                simde__m128i cond;

                vHp = simde_mm_add_epi64(vHp, simde_mm_load_si128(vP + i));
                vHp = simde_mm_max_epi64(vHp, vZero);
                vH = simde_mm_load_si128(pvHStore + i);
                vH = simde_mm_max_epi64(vH,vF);
                simde_mm_store_si128(pvHStore + i, vH);
                case1 = simde_mm_cmpeq_epi64(vH, vHp);
                case2 = simde_mm_cmpeq_epi64(vH, vF);
                cond = simde_mm_andnot_si128(case1, case2);

                /* calculate vM */
                vHM = simde_mm_load_si128(pvHMStore + i);
                vHM = simde_mm_blendv_epi8(vHM, vFM, cond);
                simde_mm_store_si128(pvHMStore + i, vHM);

                /* calculate vS */
                vHS = simde_mm_load_si128(pvHSStore + i);
                vHS = simde_mm_blendv_epi8(vHS, vFS, cond);
                simde_mm_store_si128(pvHSStore + i, vHS);

                /* calculate vL */
                vHL = simde_mm_load_si128(pvHLStore + i);
                vHL = simde_mm_blendv_epi8(vHL, vFL, cond);
                simde_mm_store_si128(pvHLStore + i, vHL);

                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHM);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHS);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
                arr_store(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
                arr_store(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
                arr_store(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                vMaxH = simde_mm_max_epi64(vH, vMaxH);
                /* Update vF value. */
                vEF_opn = simde_mm_sub_epi64(vH, vGapO);
                vF_ext = simde_mm_sub_epi64(vF, vGapE);
                if (! simde_mm_movemask_epi8(
                            simde_mm_or_si128(
                                simde_mm_cmpgt_epi64(vF_ext, vEF_opn),
                                simde_mm_and_si128(
                                    simde_mm_cmpeq_epi64(vF_ext, vEF_opn),
                                    simde_mm_cmpgt_epi64(vF_ext, vZero)))))
                    goto end;
                /*vF = simde_mm_max_epi64(vEF_opn, vF_ext);*/
                vF = vF_ext;
                cond = simde_mm_cmpgt_epi64(vEF_opn, vF_ext);
                vFM = simde_mm_blendv_epi8(vFM, vHM, cond);
                vFS = simde_mm_blendv_epi8(vFS, vHS, cond);
                vFL = simde_mm_blendv_epi8(
                        simde_mm_add_epi64(vFL, vOne),
                        simde_mm_add_epi64(vHL, vOne),
                        cond);
                vHp = simde_mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = simde_mm_load_si128(pvHStore + offset);
            vHM = simde_mm_load_si128(pvHMStore + offset);
            vHS = simde_mm_load_si128(pvHSStore + offset);
            vHL = simde_mm_load_si128(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = simde_mm_slli_si128(vH, 8);
                vHM = simde_mm_slli_si128(vHM, 8);
                vHS = simde_mm_slli_si128(vHS, 8);
                vHL = simde_mm_slli_si128(vHL, 8);
            }
            result->stats->rowcols->score_row[j] = (int64_t) simde_mm_extract_epi64 (vH, 1);
            result->stats->rowcols->matches_row[j] = (int64_t) simde_mm_extract_epi64 (vHM, 1);
            result->stats->rowcols->similar_row[j] = (int64_t) simde_mm_extract_epi64 (vHS, 1);
            result->stats->rowcols->length_row[j] = (int64_t) simde_mm_extract_epi64 (vHL, 1);
        }
#endif

        {
            simde__m128i vCompare = simde_mm_cmpgt_epi64(vMaxH, vMaxHUnit);
            if (simde_mm_movemask_epi8(vCompare)) {
                score = simde_mm_hmax_epi64(vMaxH);
                /* if score has potential to overflow, abort early */
                if (score > maxp) {
                    result->flag |= PARASAIL_FLAG_SATURATED;
                    break;
                }
                vMaxHUnit = simde_mm_set1_epi64x(score);
                end_ref = j;
            }
        }
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        simde__m128i vH = simde_mm_load_si128(pvHStore+i);
        simde__m128i vHM = simde_mm_load_si128(pvHMStore+i);
        simde__m128i vHS = simde_mm_load_si128(pvHSStore+i);
        simde__m128i vHL = simde_mm_load_si128(pvHLStore+i);
        arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
        arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
        arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
        arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
    }
#endif

    if (score == INT64_MAX
            || simde_mm_movemask_epi8(simde_mm_cmpeq_epi64(vSaturationCheckMax,vPosLimit))) {
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
            int64_t *t = (int64_t*)pvHMax;
            int64_t *m = (int64_t*)pvHMMax;
            int64_t *s = (int64_t*)pvHSMax;
            int64_t *l = (int64_t*)pvHLMax;
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

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->stats->matches = matches;
    result->stats->similar = similar;
    result->stats->length = length;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;
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
    parasail_result_t *result = parasail_sw_striped_profile_neon_128_64(profile, s2, s2Len, open, gap);
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
        profile_new = parasail_profile_create_neon_128_64(
                s1_new, s1Len_new, matrix);
        profile_new->stop = result->score;
        result_new = parasail_sw_striped_profile_neon_128_64(
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
        profile_final = parasail_profile_create_stats_neon_128_64(
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
        profile_final = parasail_profile_create_stats_neon_128_64(
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


