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



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
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
#define FNAME parasail_sg_stats_table_striped_neon_128_64
#define PNAME parasail_sg_stats_table_striped_profile_neon_128_64
#define INAME PNAME
#define STATIC
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_striped_neon_128_64
#define PNAME parasail_sg_stats_rowcol_striped_profile_neon_128_64
#define INAME PNAME
#define STATIC
#else
#define FNAME parasail_sg_stats_striped_neon_128_64
#ifdef FASTSTATS
#define PNAME parasail_sg_stats_striped_profile_neon_128_64_internal
#define INAME parasail_sg_stats_striped_profile_neon_128_64
#define STATIC static
#else
#define PNAME parasail_sg_stats_striped_profile_neon_128_64
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
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
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
    const simde__m128i vGapO = simde_mm_set1_epi64x(open);
    const simde__m128i vGapE = simde_mm_set1_epi64x(gap);
    const int64_t NEG_LIMIT = (-open < matrix->min ?
        INT64_MIN + open : INT64_MIN - matrix->min) + 1;
    const int64_t POS_LIMIT = INT64_MAX - matrix->max - 1;
    const simde__m128i vZero = simde_mm_setzero_si128();
    const simde__m128i vOne = simde_mm_set1_epi64x(1);
    int64_t score = NEG_LIMIT;
    int64_t matches = 0;
    int64_t similar = 0;
    int64_t length = 0;
    simde__m128i vNegLimit = simde_mm_set1_epi64x(NEG_LIMIT);
    simde__m128i vPosLimit = simde_mm_set1_epi64x(POS_LIMIT);
    simde__m128i vSaturationCheckMin = vPosLimit;
    simde__m128i vSaturationCheckMax = vNegLimit;
    simde__m128i vMaxH = vNegLimit;
    simde__m128i vMaxHM = vNegLimit;
    simde__m128i vMaxHS = vNegLimit;
    simde__m128i vMaxHL = vNegLimit;
    simde__m128i vPosMask = simde_mm_cmpeq_epi64(simde_mm_set1_epi64x(position),
            simde_mm_set_epi64x(0,1));
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
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

        /* Initialize F value to neg inf.  Any errors to vH values will
         * be corrected in the Lazy_F loop. */
        vF = vNegLimit;
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

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad,  pvHStore)
        SWAP(pvHMLoad, pvHMStore)
        SWAP(pvHSLoad, pvHSStore)
        SWAP(pvHLLoad, pvHLStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            simde__m128i case1;
            simde__m128i case2;

            vE = simde_mm_load_si128(pvE+ i);
            vEM = simde_mm_load_si128(pvEM+ i);
            vES = simde_mm_load_si128(pvES+ i);
            vEL = simde_mm_load_si128(pvEL+ i);

            /* Get max from vH, vE and vF. */
            vH_dag = simde_mm_add_epi64(vH, simde_mm_load_si128(vP + i));
            vH = simde_mm_max_epi64(vH_dag, vE);
            vH = simde_mm_max_epi64(vH, vF);
            /* Save vH values. */
            simde_mm_store_si128(pvHStore + i, vH);

            case1 = simde_mm_cmpeq_epi64(vH, vH_dag);
            case2 = simde_mm_cmpeq_epi64(vH, vF);

            /* calculate vM */
            vHM = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEM, vFM, case2),
                    simde_mm_add_epi64(vHM, simde_mm_load_si128(vPM + i)),
                    case1);
            simde_mm_store_si128(pvHMStore + i, vHM);

            /* calculate vS */
            vHS = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vES, vFS, case2),
                    simde_mm_add_epi64(vHS, simde_mm_load_si128(vPS + i)),
                    case1);
            simde_mm_store_si128(pvHSStore + i, vHS);

            /* calculate vL */
            vHL = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEL, vFL, case2),
                    simde_mm_add_epi64(vHL, vOne),
                    case1);
            simde_mm_store_si128(pvHLStore + i, vHL);

            vSaturationCheckMin = simde_mm_min_epi64(vSaturationCheckMin, vH);
            vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vH);
            vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHM);
            vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHS);
            vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
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

                vSaturationCheckMin = simde_mm_min_epi64(vSaturationCheckMin, vH);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vH);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHM);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHS);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si128(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si128(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si128(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                /* Update vF value. */
                vEF_opn = simde_mm_sub_epi64(vH, vGapO);
                vF_ext = simde_mm_sub_epi64(vF, vGapE);
                if (! simde_mm_movemask_epi8(
                            simde_mm_or_si128(
                                simde_mm_cmpgt_epi64(vF_ext, vEF_opn),
                                simde_mm_cmpeq_epi64(vF_ext, vEF_opn))))
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

        /* extract vector containing last value from the column */
        {
            simde__m128i cond_max;
            vH = simde_mm_load_si128(pvHStore + offset);
            vHM = simde_mm_load_si128(pvHMStore + offset);
            vHS = simde_mm_load_si128(pvHSStore + offset);
            vHL = simde_mm_load_si128(pvHLStore + offset);
            cond_max = simde_mm_cmpgt_epi64(vH, vMaxH);
            vMaxH = simde_mm_blendv_epi8(vMaxH, vH, cond_max);
            vMaxHM = simde_mm_blendv_epi8(vMaxHM, vHM, cond_max);
            vMaxHS = simde_mm_blendv_epi8(vMaxHS, vHS, cond_max);
            vMaxHL = simde_mm_blendv_epi8(vMaxHL, vHL, cond_max);
            if (simde_mm_movemask_epi8(simde_mm_and_si128(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
            }
        }
#ifdef PARASAIL_ROWCOL
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
#endif
    }

    {
        /* extract last value from the column maximums */
        for (k=0; k<position; ++k) {
            vMaxH  = simde_mm_slli_si128 (vMaxH, 8);
            vMaxHM = simde_mm_slli_si128 (vMaxHM, 8);
            vMaxHS = simde_mm_slli_si128 (vMaxHS, 8);
            vMaxHL = simde_mm_slli_si128 (vMaxHL, 8);
        }
        score = (int64_t) simde_mm_extract_epi64 (vMaxH, 1);
        matches = (int64_t)simde_mm_extract_epi64(vMaxHM, 1);
        similar = (int64_t)simde_mm_extract_epi64(vMaxHS, 1);
        length = (int64_t)simde_mm_extract_epi64(vMaxHL, 1);
    }

    /* max of last column */
    if (INT32_MAX == profile->stop || 0 == profile->stop)
    {
        int64_t score_last;
        vMaxH = vNegLimit;

        if (0 == profile->stop) {
            /* ignore last row contributions */
            score = NEG_LIMIT;
            matches = 0;
            similar = 0;
            length = 0;
            end_query = s1Len;
            end_ref = s2Len - 1;
        }

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            simde__m128i vH = simde_mm_load_si128(pvHStore + i);
#ifdef PARASAIL_ROWCOL
            simde__m128i vHM = simde_mm_load_si128(pvHMStore + i);
            simde__m128i vHS = simde_mm_load_si128(pvHSStore + i);
            simde__m128i vHL = simde_mm_load_si128(pvHLStore + i);
            arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
            arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
            arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
            arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
#endif
            vMaxH = simde_mm_max_epi64(vH, vMaxH);
        }

        /* max in vec */
        score_last = simde_mm_hmax_epi64(vMaxH);
        if (score_last > score || (score_last == score && end_ref == s2Len - 1)) {
            end_query = s1Len;
            end_ref = s2Len - 1;
            /* Trace the alignment ending position on read. */
            {
                int64_t *t = (int64_t*)pvHStore;
                int64_t *m = (int64_t*)pvHMStore;
                int64_t *s = (int64_t*)pvHSStore;
                int64_t *l = (int64_t*)pvHLStore;
                int32_t column_len = segLen * segWidth;
                for (i = 0; i<column_len; ++i, ++t, ++m, ++s, ++l) {
                    int32_t temp = i / segWidth + i % segWidth * segLen;
                    if (temp < s1Len) {
                        if (*t > score || (*t == score && temp < end_query)) {
                            score = *t;
                            end_query = temp;
                            matches = *m;
                            similar = *s;
                            length = *l;
                        }
                    }
                }
            }
        }
    }

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmplt_epi64(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpgt_epi64(vSaturationCheckMax, vPosLimit)))) {
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
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

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
    parasail_result_t *result = parasail_sg_striped_profile_neon_128_64(profile, s2, s2Len, open, gap);
    if (!parasail_result_is_saturated(result)) {
        int s1Len_new = 0;
        int s2Len_new = 0;
        parasail_result_t *result_final = NULL;

        /* using the end loc, call the original stats function */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;

        if (s1Len_new == profile->s1Len) {
            /* special 'stop' value tells stats function not to
             * consider last column results */
            int stop_save = profile->stop;
            ((parasail_profile_t*)profile)->stop = 1;
            result_final = PNAME(
                    profile, s2, s2Len_new, open, gap);
            ((parasail_profile_t*)profile)->stop = stop_save;
        }
        else {
            parasail_profile_t *profile_final = NULL;
            profile_final = parasail_profile_create_stats_neon_128_64(
                    s1, s1Len_new, matrix);
            /* special 'stop' value tells stats function not to
             * consider last row results */
            profile_final->stop = 0;
            result_final = PNAME(
                    profile_final, s2, s2Len_new, open, gap);

            parasail_profile_free(profile_final);
        }

        parasail_result_free(result);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
    }
    else {
        return result;
    }
}
#endif
#endif
#endif


