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



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        simde__m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*( 0*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  0);
    array[1LL*( 1*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  1);
    array[1LL*( 2*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  2);
    array[1LL*( 3*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  3);
    array[1LL*( 4*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  4);
    array[1LL*( 5*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  5);
    array[1LL*( 6*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  6);
    array[1LL*( 7*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  7);
    array[1LL*( 8*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  8);
    array[1LL*( 9*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH,  9);
    array[1LL*(10*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH, 10);
    array[1LL*(11*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH, 11);
    array[1LL*(12*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH, 12);
    array[1LL*(13*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH, 13);
    array[1LL*(14*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH, 14);
    array[1LL*(15*seglen+t)*dlen + d] = (int8_t)simde_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        simde__m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  0);
    col[ 1*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  1);
    col[ 2*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  2);
    col[ 3*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  3);
    col[ 4*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  4);
    col[ 5*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  5);
    col[ 6*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  6);
    col[ 7*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  7);
    col[ 8*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  8);
    col[ 9*seglen+t] = (int8_t)simde_mm_extract_epi8(vH,  9);
    col[10*seglen+t] = (int8_t)simde_mm_extract_epi8(vH, 10);
    col[11*seglen+t] = (int8_t)simde_mm_extract_epi8(vH, 11);
    col[12*seglen+t] = (int8_t)simde_mm_extract_epi8(vH, 12);
    col[13*seglen+t] = (int8_t)simde_mm_extract_epi8(vH, 13);
    col[14*seglen+t] = (int8_t)simde_mm_extract_epi8(vH, 14);
    col[15*seglen+t] = (int8_t)simde_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_table_striped_neon_128_8
#define PNAME parasail_sg_table_striped_profile_neon_128_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_rowcol_striped_neon_128_8
#define PNAME parasail_sg_rowcol_striped_profile_neon_128_8
#else
#define FNAME parasail_sg_striped_neon_128_8
#define PNAME parasail_sg_striped_profile_neon_128_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_neon_128_8(s1, s1Len, matrix);
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
    simde__m128i* const restrict vProfile = (simde__m128i*)profile->profile8.score;
    simde__m128i* restrict pvHStore = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* restrict pvHLoad =  parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvE = parasail_memalign_simde__m128i(16, segLen);
    const simde__m128i vGapO = simde_mm_set1_epi8(open);
    const simde__m128i vGapE = simde_mm_set1_epi8(gap);
    const int8_t NEG_LIMIT = (-open < matrix->min ?
        INT8_MIN + open : INT8_MIN - matrix->min) + 1;
    const int8_t POS_LIMIT = INT8_MAX - matrix->max - 1;
    int8_t score = NEG_LIMIT;
    simde__m128i vNegLimit = simde_mm_set1_epi8(NEG_LIMIT);
    simde__m128i vPosLimit = simde_mm_set1_epi8(POS_LIMIT);
    simde__m128i vSaturationCheckMin = vPosLimit;
    simde__m128i vSaturationCheckMax = vNegLimit;
    simde__m128i vMaxH = vNegLimit;
    simde__m128i vPosMask = simde_mm_cmpeq_epi8(simde_mm_set1_epi8(position),
            simde_mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15));
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    /* initialize H and E */
    parasail_memset_simde__m128i(pvHStore, simde_mm_set1_epi8(0), segLen);
    parasail_memset_simde__m128i(pvE, simde_mm_set1_epi8(-open), segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vE;
        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        simde__m128i vF = vNegLimit;

        /* load final segment of pvHStore and shift left by 2 bytes */
        simde__m128i vH = simde_mm_slli_si128(pvHStore[segLen - 1], 1);

        /* Correct part of the vProfile */
        const simde__m128i* vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        simde__m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = simde_mm_adds_epi8(vH, simde_mm_load_si128(vP + i));
            vE = simde_mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = simde_mm_max_epi8(vH, vE);
            vH = simde_mm_max_epi8(vH, vF);
            /* Save vH values. */
            simde_mm_store_si128(pvHStore + i, vH);
            vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vH);
            vSaturationCheckMax = simde_mm_max_epi8(vSaturationCheckMax, vH);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = simde_mm_subs_epi8(vH, vGapO);
            vE = simde_mm_subs_epi8(vE, vGapE);
            vE = simde_mm_max_epi8(vE, vH);
            simde_mm_store_si128(pvE + i, vE);

            /* Update vF value. */
            vF = simde_mm_subs_epi8(vF, vGapE);
            vF = simde_mm_max_epi8(vF, vH);

            /* Load the next vH. */
            vH = simde_mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            vF = simde_mm_slli_si128(vF, 1);
            vF = simde_mm_insert_epi8(vF, -open, 0);
            for (i=0; i<segLen; ++i) {
                vH = simde_mm_load_si128(pvHStore + i);
                vH = simde_mm_max_epi8(vH,vF);
                simde_mm_store_si128(pvHStore + i, vH);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vH);
                vSaturationCheckMax = simde_mm_max_epi8(vSaturationCheckMax, vH);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = simde_mm_subs_epi8(vH, vGapO);
                vF = simde_mm_subs_epi8(vF, vGapE);
                if (! simde_mm_movemask_epi8(simde_mm_cmpgt_epi8(vF, vH))) goto end;
                /*vF = simde_mm_max_epi8(vF, vH);*/
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            simde__m128i vCompare;
            vH = simde_mm_load_si128(pvHStore + offset);
            vCompare = simde_mm_and_si128(vPosMask, simde_mm_cmpgt_epi8(vH, vMaxH));
            vMaxH = simde_mm_max_epi8(vH, vMaxH);
            if (simde_mm_movemask_epi8(vCompare)) {
                end_ref = j;
                end_query = s1Len - 1;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = simde_mm_slli_si128(vH, 1);
            }
            result->rowcols->score_row[j] = (int8_t) simde_mm_extract_epi8 (vH, 15);
#endif
        }
    }

    /* max last value from all columns */
    {
        for (k=0; k<position; ++k) {
            vMaxH = simde_mm_slli_si128(vMaxH, 1);
        }
        score = (int8_t) simde_mm_extract_epi8(vMaxH, 15);
    }

    /* max of last column */
    {
        int8_t score_last;
        vMaxH = vNegLimit;

        for (i=0; i<segLen; ++i) {
            simde__m128i vH = simde_mm_load_si128(pvHStore + i);
            vMaxH = simde_mm_max_epi8(vH, vMaxH);
#ifdef PARASAIL_ROWCOL
            arr_store_col(result->rowcols->score_col, vH, i, segLen);
#endif
        }

        /* max in vec */
        score_last = simde_mm_hmax_epi8(vMaxH);
        if (score_last > score || (score_last == score && end_ref == s2Len - 1)) {
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                int8_t *t = (int8_t*)pvHStore;
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

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmplt_epi8(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpgt_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

