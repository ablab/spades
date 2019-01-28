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
    array[1LL*(0*seglen+t)*dlen + d] = (int32_t)simde_mm_extract_epi32(vH, 0);
    array[1LL*(1*seglen+t)*dlen + d] = (int32_t)simde_mm_extract_epi32(vH, 1);
    array[1LL*(2*seglen+t)*dlen + d] = (int32_t)simde_mm_extract_epi32(vH, 2);
    array[1LL*(3*seglen+t)*dlen + d] = (int32_t)simde_mm_extract_epi32(vH, 3);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        simde__m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int32_t)simde_mm_extract_epi32(vH, 0);
    col[1*seglen+t] = (int32_t)simde_mm_extract_epi32(vH, 1);
    col[2*seglen+t] = (int32_t)simde_mm_extract_epi32(vH, 2);
    col[3*seglen+t] = (int32_t)simde_mm_extract_epi32(vH, 3);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_stats_table_scan_neon_128_32
#define PNAME parasail_sg_stats_table_scan_profile_neon_128_32
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_scan_neon_128_32
#define PNAME parasail_sg_stats_rowcol_scan_profile_neon_128_32
#else
#define FNAME parasail_sg_stats_scan_neon_128_32
#define PNAME parasail_sg_stats_scan_profile_neon_128_32
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_neon_128_32(s1, s1Len, matrix);
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
    const int32_t segWidth = 4; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    simde__m128i* const restrict pvP  = (simde__m128i*)profile->profile32.score;
    simde__m128i* const restrict pvPm = (simde__m128i*)profile->profile32.matches;
    simde__m128i* const restrict pvPs = (simde__m128i*)profile->profile32.similar;
    simde__m128i* const restrict pvE  = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvEM = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvES = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvEL = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvH  = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHM = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHS = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHL = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHMax  = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHMMax = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHSMax = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvHLMax = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvGapper = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i* const restrict pvGapperL = parasail_memalign_simde__m128i(16, segLen);
    simde__m128i vGapO = simde_mm_set1_epi32(open);
    simde__m128i vGapE = simde_mm_set1_epi32(gap);
    const int32_t NEG_LIMIT = (-open < matrix->min ?
        INT32_MIN + open : INT32_MIN - matrix->min) + 1;
    const int32_t POS_LIMIT = INT32_MAX - matrix->max - 1;
    simde__m128i vZero = simde_mm_setzero_si128();
    simde__m128i vOne = simde_mm_set1_epi32(1);
    int32_t score = NEG_LIMIT;
    int32_t matches = 0;
    int32_t similar = 0;
    int32_t length = 0;
    simde__m128i vNegLimit = simde_mm_set1_epi32(NEG_LIMIT);
    simde__m128i vPosLimit = simde_mm_set1_epi32(POS_LIMIT);
    simde__m128i vSaturationCheckMin = vPosLimit;
    simde__m128i vSaturationCheckMax = vNegLimit;
    simde__m128i vMaxH = vNegLimit;
    simde__m128i vMaxM = vNegLimit;
    simde__m128i vMaxS = vNegLimit;
    simde__m128i vMaxL = vNegLimit;
    simde__m128i vPosMask = simde_mm_cmpeq_epi32(simde_mm_set1_epi32(position),
            simde_mm_set_epi32(0,1,2,3));
    simde__m128i vNegInfFront = vZero;
    simde__m128i vSegLenXgap;
    simde__m128i vSegLen = simde_mm_slli_si128(simde_mm_set1_epi32(segLen), 4);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    vNegInfFront = simde_mm_insert_epi32(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = simde_mm_add_epi32(vNegInfFront,
            simde_mm_slli_si128(simde_mm_set1_epi32(-segLen*gap), 4));

    parasail_memset_simde__m128i(pvH, vZero, segLen);
    parasail_memset_simde__m128i(pvHM, vZero, segLen);
    parasail_memset_simde__m128i(pvHS, vZero, segLen);
    parasail_memset_simde__m128i(pvHL, vZero, segLen);
    parasail_memset_simde__m128i(pvE, vNegLimit, segLen);
    parasail_memset_simde__m128i(pvEM, vZero, segLen);
    parasail_memset_simde__m128i(pvES, vZero, segLen);
    parasail_memset_simde__m128i(pvEL, vZero, segLen);
    {
        simde__m128i vGapper = simde_mm_sub_epi32(vZero,vGapO);
        simde__m128i vGapperL = vOne;
        for (i=segLen-1; i>=0; --i) {
            simde_mm_store_si128(pvGapper+i, vGapper);
            simde_mm_store_si128(pvGapperL+i, vGapperL);
            vGapper = simde_mm_sub_epi32(vGapper, vGapE);
            vGapperL = simde_mm_add_epi32(vGapperL, vOne);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vE;
        simde__m128i vE_ext;
        simde__m128i vE_opn;
        simde__m128i vEM;
        simde__m128i vES;
        simde__m128i vEL;
        simde__m128i vHt;
        simde__m128i vHtM;
        simde__m128i vHtS;
        simde__m128i vHtL;
        simde__m128i vF;
        simde__m128i vF_ext;
        simde__m128i vF_opn;
        simde__m128i vFM;
        simde__m128i vFS;
        simde__m128i vFL;
        simde__m128i vH;
        simde__m128i vHM;
        simde__m128i vHS;
        simde__m128i vHL;
        simde__m128i vHp;
        simde__m128i vHpM;
        simde__m128i vHpS;
        simde__m128i vHpL;
        simde__m128i *pvW;
        simde__m128i vW;
        simde__m128i *pvWM;
        simde__m128i vWM;
        simde__m128i *pvWS;
        simde__m128i vWS;
        simde__m128i case1;
        simde__m128i case2;
        simde__m128i vGapper;
        simde__m128i vGapperL;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = simde_mm_load_si128(pvH+(segLen-1));
        vHpM = simde_mm_load_si128(pvHM+(segLen-1));
        vHpS = simde_mm_load_si128(pvHS+(segLen-1));
        vHpL = simde_mm_load_si128(pvHL+(segLen-1));
        vHp = simde_mm_slli_si128(vHp, 4);
        vHpM = simde_mm_slli_si128(vHpM, 4);
        vHpS = simde_mm_slli_si128(vHpS, 4);
        vHpL = simde_mm_slli_si128(vHpL, 4);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWM= pvPm+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWS= pvPs+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = simde_mm_sub_epi32(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        vFM = vZero;
        vFS = vZero;
        vFL = vZero;
        for (i=0; i<segLen; ++i) {
            vH = simde_mm_load_si128(pvH+i);
            vHM= simde_mm_load_si128(pvHM+i);
            vHS= simde_mm_load_si128(pvHS+i);
            vHL= simde_mm_load_si128(pvHL+i);
            vE = simde_mm_load_si128(pvE+i);
            vEM= simde_mm_load_si128(pvEM+i);
            vES= simde_mm_load_si128(pvES+i);
            vEL= simde_mm_load_si128(pvEL+i);
            vW = simde_mm_load_si128(pvW+i);
            vWM = simde_mm_load_si128(pvWM+i);
            vWS = simde_mm_load_si128(pvWS+i);
            vGapper = simde_mm_load_si128(pvGapper+i);
            vGapperL = simde_mm_load_si128(pvGapperL+i);
            vE_opn = simde_mm_sub_epi32(vH, vGapO);
            vE_ext = simde_mm_sub_epi32(vE, vGapE);
            case1 = simde_mm_cmpgt_epi32(vE_opn, vE_ext);
            vE = simde_mm_max_epi32(vE_opn, vE_ext);
            vEM = simde_mm_blendv_epi8(vEM, vHM, case1);
            vES = simde_mm_blendv_epi8(vES, vHS, case1);
            vEL = simde_mm_blendv_epi8(vEL, vHL, case1);
            vEL = simde_mm_add_epi32(vEL, vOne);
            vGapper = simde_mm_add_epi32(vHt, vGapper);
            case1 = simde_mm_or_si128(
                    simde_mm_cmpgt_epi32(vF, vGapper),
                    simde_mm_cmpeq_epi32(vF, vGapper));
            vF = simde_mm_max_epi32(vF, vGapper);
            vFM = simde_mm_blendv_epi8(vHtM, vFM, case1);
            vFS = simde_mm_blendv_epi8(vHtS, vFS, case1);
            vFL = simde_mm_blendv_epi8(
                    simde_mm_add_epi32(vHtL, vGapperL),
                    vFL, case1);
            vHp = simde_mm_add_epi32(vHp, vW);
            vHpM = simde_mm_add_epi32(vHpM, vWM);
            vHpS = simde_mm_add_epi32(vHpS, vWS);
            vHpL = simde_mm_add_epi32(vHpL, vOne);
            case1 = simde_mm_cmpgt_epi32(vE, vHp);
            vHt = simde_mm_max_epi32(vE, vHp);
            vHtM = simde_mm_blendv_epi8(vHpM, vEM, case1);
            vHtS = simde_mm_blendv_epi8(vHpS, vES, case1);
            vHtL = simde_mm_blendv_epi8(vHpL, vEL, case1);
            simde_mm_store_si128(pvE+i, vE);
            simde_mm_store_si128(pvEM+i, vEM);
            simde_mm_store_si128(pvES+i, vES);
            simde_mm_store_si128(pvEL+i, vEL);
            simde_mm_store_si128(pvH+i, vHp);
            simde_mm_store_si128(pvHM+i, vHpM);
            simde_mm_store_si128(pvHS+i, vHpS);
            simde_mm_store_si128(pvHL+i, vHpL);
            vHp = vH;
            vHpM = vHM;
            vHpS = vHS;
            vHpL = vHL;
        }

        /* pseudo prefix scan on F and H */
        vHt = simde_mm_slli_si128(vHt, 4);
        vHtM = simde_mm_slli_si128(vHtM, 4);
        vHtS = simde_mm_slli_si128(vHtS, 4);
        vHtL = simde_mm_slli_si128(vHtL, 4);
        vGapper = simde_mm_load_si128(pvGapper);
        vGapperL = simde_mm_load_si128(pvGapperL);
        vGapper = simde_mm_add_epi32(vHt, vGapper);
        case1 = simde_mm_or_si128(
                simde_mm_cmpgt_epi32(vGapper, vF),
                simde_mm_cmpeq_epi32(vGapper, vF));
        vF = simde_mm_max_epi32(vF, vGapper);
        vFM = simde_mm_blendv_epi8(vFM, vHtM, case1);
        vFS = simde_mm_blendv_epi8(vFS, vHtS, case1);
        vFL = simde_mm_blendv_epi8(
                vFL,
                simde_mm_add_epi32(vHtL, vGapperL),
                case1);
        for (i=0; i<segWidth-2; ++i) {
            simde__m128i vFt = simde_mm_slli_si128(vF, 4);
            simde__m128i vFtM = simde_mm_slli_si128(vFM, 4);
            simde__m128i vFtS = simde_mm_slli_si128(vFS, 4);
            simde__m128i vFtL = simde_mm_slli_si128(vFL, 4);
            vFt = simde_mm_add_epi32(vFt, vSegLenXgap);
            case1 = simde_mm_or_si128(
                    simde_mm_cmpgt_epi32(vFt, vF),
                    simde_mm_cmpeq_epi32(vFt, vF));
            vF = simde_mm_max_epi32(vF, vFt);
            vFM = simde_mm_blendv_epi8(vFM, vFtM, case1);
            vFS = simde_mm_blendv_epi8(vFS, vFtS, case1);
            vFL = simde_mm_blendv_epi8(
                    vFL,
                    simde_mm_add_epi32(vFtL, vSegLen),
                    case1);
        }

        /* calculate final H */
        vF = simde_mm_slli_si128(vF, 4);
        vFM = simde_mm_slli_si128(vFM, 4);
        vFS = simde_mm_slli_si128(vFS, 4);
        vFL = simde_mm_slli_si128(vFL, 4);
        vF = simde_mm_add_epi32(vF, vNegInfFront);
        case1 = simde_mm_cmpgt_epi32(vF, vHt);
        vH = simde_mm_max_epi32(vF, vHt);
        vHM = simde_mm_blendv_epi8(vHtM, vFM, case1);
        vHS = simde_mm_blendv_epi8(vHtS, vFS, case1);
        vHL = simde_mm_blendv_epi8(vHtL, vFL, case1);
        for (i=0; i<segLen; ++i) {
            vHp = simde_mm_load_si128(pvH+i);
            vHpM = simde_mm_load_si128(pvHM+i);
            vHpS = simde_mm_load_si128(pvHS+i);
            vHpL = simde_mm_load_si128(pvHL+i);
            vE = simde_mm_load_si128(pvE+i);
            vEM = simde_mm_load_si128(pvEM+i);
            vES = simde_mm_load_si128(pvES+i);
            vEL = simde_mm_load_si128(pvEL+i);
            vF_opn = simde_mm_sub_epi32(vH, vGapO);
            vF_ext = simde_mm_sub_epi32(vF, vGapE);
            vF = simde_mm_max_epi32(vF_opn, vF_ext);
            case1 = simde_mm_cmpgt_epi32(vF_opn, vF_ext);
            vFM = simde_mm_blendv_epi8(vFM, vHM, case1);
            vFS = simde_mm_blendv_epi8(vFS, vHS, case1);
            vFL = simde_mm_blendv_epi8(vFL, vHL, case1);
            vFL = simde_mm_add_epi32(vFL, vOne);
            vH = simde_mm_max_epi32(vHp, vE);
            vH = simde_mm_max_epi32(vH, vF);
            case1 = simde_mm_cmpeq_epi32(vH, vHp);
            case2 = simde_mm_cmpeq_epi32(vH, vF);
            vHM = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEM, vFM, case2),
                    vHpM, case1);
            vHS = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vES, vFS, case2),
                    vHpS, case1);
            vHL = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEL, vFL, case2),
                    vHpL, case1);
            simde_mm_store_si128(pvH+i, vH);
            simde_mm_store_si128(pvHM+i, vHM);
            simde_mm_store_si128(pvHS+i, vHS);
            simde_mm_store_si128(pvHL+i, vHL);
            vSaturationCheckMin = simde_mm_min_epi32(vSaturationCheckMin, vH);
            vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vH);
            vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vHM);
            vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vHS);
            vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
#endif
        } 
        /* extract vector containing last value from column */
        {
            simde__m128i cond_max;
            vH = simde_mm_load_si128(pvH + offset);
            vHM = simde_mm_load_si128(pvHM + offset);
            vHS = simde_mm_load_si128(pvHS + offset);
            vHL = simde_mm_load_si128(pvHL + offset);
            cond_max = simde_mm_cmpgt_epi32(vH, vMaxH);
            vMaxH = simde_mm_max_epi32(vMaxH, vH);
            vMaxM = simde_mm_blendv_epi8(vMaxM, vHM, cond_max);
            vMaxS = simde_mm_blendv_epi8(vMaxS, vHS, cond_max);
            vMaxL = simde_mm_blendv_epi8(vMaxL, vHL, cond_max);
            if (simde_mm_movemask_epi8(simde_mm_and_si128(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = simde_mm_slli_si128(vH, 4);
                vHM = simde_mm_slli_si128(vHM, 4);
                vHS = simde_mm_slli_si128(vHS, 4);
                vHL = simde_mm_slli_si128(vHL, 4);
            }
            result->stats->rowcols->score_row[j] = (int32_t) simde_mm_extract_epi32 (vH, 3);
            result->stats->rowcols->matches_row[j] = (int32_t) simde_mm_extract_epi32 (vHM, 3);
            result->stats->rowcols->similar_row[j] = (int32_t) simde_mm_extract_epi32 (vHS, 3);
            result->stats->rowcols->length_row[j] = (int32_t) simde_mm_extract_epi32 (vHL, 3);
#endif
        }
    }

    /* max last value from all columns */
    {
        for (k=0; k<position; ++k) {
            vMaxH = simde_mm_slli_si128(vMaxH, 4);
            vMaxM = simde_mm_slli_si128(vMaxM, 4);
            vMaxS = simde_mm_slli_si128(vMaxS, 4);
            vMaxL = simde_mm_slli_si128(vMaxL, 4);
        }
        score = (int32_t) simde_mm_extract_epi32(vMaxH, 3);
        matches = (int32_t) simde_mm_extract_epi32(vMaxM, 3);
        similar = (int32_t) simde_mm_extract_epi32(vMaxS, 3);
        length = (int32_t) simde_mm_extract_epi32(vMaxL, 3);
    }

    /* max of last column */
    {
        int32_t score_last;
        vMaxH = vNegLimit;

        for (i=0; i<segLen; ++i) {
            simde__m128i vH = simde_mm_load_si128(pvH + i);
#ifdef PARASAIL_ROWCOL
            simde__m128i vHM = simde_mm_load_si128(pvHM + i);
            simde__m128i vHS = simde_mm_load_si128(pvHS + i);
            simde__m128i vHL = simde_mm_load_si128(pvHL + i);
            arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
            arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
            arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
            arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
#endif
            vMaxH = simde_mm_max_epi32(vH, vMaxH);
        }

        /* max in vec */
        score_last = simde_mm_hmax_epi32(vMaxH);
        if (score_last > score || (score_last == score && end_ref == s2Len - 1)) {
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                int32_t *t = (int32_t*)pvH;
                int32_t *m = (int32_t*)pvHM;
                int32_t *s = (int32_t*)pvHS;
                int32_t *l = (int32_t*)pvHL;
                int32_t column_len = segLen * segWidth;
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
    }

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmplt_epi32(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpgt_epi32(vSaturationCheckMax, vPosLimit)))) {
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
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_32 | PARASAIL_FLAG_LANES_4;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(pvGapperL);
    parasail_free(pvGapper);
    parasail_free(pvHLMax);
    parasail_free(pvHSMax);
    parasail_free(pvHMMax);
    parasail_free(pvHMax);
    parasail_free(pvHL);
    parasail_free(pvHS);
    parasail_free(pvHM);
    parasail_free(pvH);
    parasail_free(pvEL);
    parasail_free(pvES);
    parasail_free(pvEM);
    parasail_free(pvE);

    return result;
}


