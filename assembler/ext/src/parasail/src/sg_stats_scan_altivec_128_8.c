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



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        vec128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*( 0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  0);
    array[1LL*( 1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  1);
    array[1LL*( 2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  2);
    array[1LL*( 3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  3);
    array[1LL*( 4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  4);
    array[1LL*( 5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  5);
    array[1LL*( 6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  6);
    array[1LL*( 7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  7);
    array[1LL*( 8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  8);
    array[1LL*( 9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  9);
    array[1LL*(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 10);
    array[1LL*(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 11);
    array[1LL*(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 12);
    array[1LL*(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 13);
    array[1LL*(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 14);
    array[1LL*(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        vec128i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int8_t)_mm_extract_epi8(vH,  0);
    col[ 1*seglen+t] = (int8_t)_mm_extract_epi8(vH,  1);
    col[ 2*seglen+t] = (int8_t)_mm_extract_epi8(vH,  2);
    col[ 3*seglen+t] = (int8_t)_mm_extract_epi8(vH,  3);
    col[ 4*seglen+t] = (int8_t)_mm_extract_epi8(vH,  4);
    col[ 5*seglen+t] = (int8_t)_mm_extract_epi8(vH,  5);
    col[ 6*seglen+t] = (int8_t)_mm_extract_epi8(vH,  6);
    col[ 7*seglen+t] = (int8_t)_mm_extract_epi8(vH,  7);
    col[ 8*seglen+t] = (int8_t)_mm_extract_epi8(vH,  8);
    col[ 9*seglen+t] = (int8_t)_mm_extract_epi8(vH,  9);
    col[10*seglen+t] = (int8_t)_mm_extract_epi8(vH, 10);
    col[11*seglen+t] = (int8_t)_mm_extract_epi8(vH, 11);
    col[12*seglen+t] = (int8_t)_mm_extract_epi8(vH, 12);
    col[13*seglen+t] = (int8_t)_mm_extract_epi8(vH, 13);
    col[14*seglen+t] = (int8_t)_mm_extract_epi8(vH, 14);
    col[15*seglen+t] = (int8_t)_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_stats_table_scan_altivec_128_8
#define PNAME parasail_sg_stats_table_scan_profile_altivec_128_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_scan_altivec_128_8
#define PNAME parasail_sg_stats_rowcol_scan_profile_altivec_128_8
#else
#define FNAME parasail_sg_stats_scan_altivec_128_8
#define PNAME parasail_sg_stats_scan_profile_altivec_128_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_altivec_128_8(s1, s1Len, matrix);
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
    vec128i* const restrict pvP  = (vec128i*)profile->profile8.score;
    vec128i* const restrict pvPm = (vec128i*)profile->profile8.matches;
    vec128i* const restrict pvPs = (vec128i*)profile->profile8.similar;
    vec128i* const restrict pvE  = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvEM = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvES = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvEL = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvH  = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvHM = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvHS = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvHL = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvHMax  = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvHMMax = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvHSMax = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvHLMax = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvGapper = parasail_memalign_vec128i(16, segLen);
    vec128i* const restrict pvGapperL = parasail_memalign_vec128i(16, segLen);
    vec128i vGapO = _mm_set1_epi8(open);
    vec128i vGapE = _mm_set1_epi8(gap);
    const int8_t NEG_LIMIT = (-open < matrix->min ?
        INT8_MIN + open : INT8_MIN - matrix->min) + 1;
    const int8_t POS_LIMIT = INT8_MAX - matrix->max - 1;
    vec128i vZero = _mm_setzero_si128();
    vec128i vOne = _mm_set1_epi8(1);
    int8_t score = NEG_LIMIT;
    int8_t matches = 0;
    int8_t similar = 0;
    int8_t length = 0;
    vec128i vNegLimit = _mm_set1_epi8(NEG_LIMIT);
    vec128i vPosLimit = _mm_set1_epi8(POS_LIMIT);
    vec128i vSaturationCheckMin = vPosLimit;
    vec128i vSaturationCheckMax = vNegLimit;
    vec128i vMaxH = vNegLimit;
    vec128i vMaxM = vNegLimit;
    vec128i vMaxS = vNegLimit;
    vec128i vMaxL = vNegLimit;
    vec128i vPosMask = _mm_cmpeq_epi8(_mm_set1_epi8(position),
            _mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15));
    vec128i vNegInfFront = vZero;
    vec128i vSegLenXgap;
    vec128i vSegLen = _mm_slli_si128(_mm_set1_epi8(segLen), 1);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    vNegInfFront = _mm_insert_epi8(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = _mm_adds_epi8(vNegInfFront,
            _mm_slli_si128(_mm_set1_epi8(-segLen*gap), 1));

    parasail_memset_vec128i(pvH, vZero, segLen);
    parasail_memset_vec128i(pvHM, vZero, segLen);
    parasail_memset_vec128i(pvHS, vZero, segLen);
    parasail_memset_vec128i(pvHL, vZero, segLen);
    parasail_memset_vec128i(pvE, vNegLimit, segLen);
    parasail_memset_vec128i(pvEM, vZero, segLen);
    parasail_memset_vec128i(pvES, vZero, segLen);
    parasail_memset_vec128i(pvEL, vZero, segLen);
    {
        vec128i vGapper = _mm_subs_epi8(vZero,vGapO);
        vec128i vGapperL = vOne;
        for (i=segLen-1; i>=0; --i) {
            _mm_store_si128(pvGapper+i, vGapper);
            _mm_store_si128(pvGapperL+i, vGapperL);
            vGapper = _mm_subs_epi8(vGapper, vGapE);
            vGapperL = _mm_adds_epi8(vGapperL, vOne);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        vec128i vE;
        vec128i vE_ext;
        vec128i vE_opn;
        vec128i vEM;
        vec128i vES;
        vec128i vEL;
        vec128i vHt;
        vec128i vHtM;
        vec128i vHtS;
        vec128i vHtL;
        vec128i vF;
        vec128i vF_ext;
        vec128i vF_opn;
        vec128i vFM;
        vec128i vFS;
        vec128i vFL;
        vec128i vH;
        vec128i vHM;
        vec128i vHS;
        vec128i vHL;
        vec128i vHp;
        vec128i vHpM;
        vec128i vHpS;
        vec128i vHpL;
        vec128i *pvW;
        vec128i vW;
        vec128i *pvWM;
        vec128i vWM;
        vec128i *pvWS;
        vec128i vWS;
        vec128i case1;
        vec128i case2;
        vec128i vGapper;
        vec128i vGapperL;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = _mm_load_si128(pvH+(segLen-1));
        vHpM = _mm_load_si128(pvHM+(segLen-1));
        vHpS = _mm_load_si128(pvHS+(segLen-1));
        vHpL = _mm_load_si128(pvHL+(segLen-1));
        vHp = _mm_slli_si128(vHp, 1);
        vHpM = _mm_slli_si128(vHpM, 1);
        vHpS = _mm_slli_si128(vHpS, 1);
        vHpL = _mm_slli_si128(vHpL, 1);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWM= pvPm+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWS= pvPs+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = _mm_subs_epi8(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        vFM = vZero;
        vFS = vZero;
        vFL = vZero;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vHM= _mm_load_si128(pvHM+i);
            vHS= _mm_load_si128(pvHS+i);
            vHL= _mm_load_si128(pvHL+i);
            vE = _mm_load_si128(pvE+i);
            vEM= _mm_load_si128(pvEM+i);
            vES= _mm_load_si128(pvES+i);
            vEL= _mm_load_si128(pvEL+i);
            vW = _mm_load_si128(pvW+i);
            vWM = _mm_load_si128(pvWM+i);
            vWS = _mm_load_si128(pvWS+i);
            vGapper = _mm_load_si128(pvGapper+i);
            vGapperL = _mm_load_si128(pvGapperL+i);
            vE_opn = _mm_subs_epi8(vH, vGapO);
            vE_ext = _mm_subs_epi8(vE, vGapE);
            case1 = _mm_cmpgt_epi8(vE_opn, vE_ext);
            vE = _mm_max_epi8(vE_opn, vE_ext);
            vEM = _mm_blendv_epi8(vEM, vHM, case1);
            vES = _mm_blendv_epi8(vES, vHS, case1);
            vEL = _mm_blendv_epi8(vEL, vHL, case1);
            vEL = _mm_adds_epi8(vEL, vOne);
            vGapper = _mm_adds_epi8(vHt, vGapper);
            case1 = _mm_or_si128(
                    _mm_cmpgt_epi8(vF, vGapper),
                    _mm_cmpeq_epi8(vF, vGapper));
            vF = _mm_max_epi8(vF, vGapper);
            vFM = _mm_blendv_epi8(vHtM, vFM, case1);
            vFS = _mm_blendv_epi8(vHtS, vFS, case1);
            vFL = _mm_blendv_epi8(
                    _mm_adds_epi8(vHtL, vGapperL),
                    vFL, case1);
            vHp = _mm_adds_epi8(vHp, vW);
            vHpM = _mm_adds_epi8(vHpM, vWM);
            vHpS = _mm_adds_epi8(vHpS, vWS);
            vHpL = _mm_adds_epi8(vHpL, vOne);
            case1 = _mm_cmpgt_epi8(vE, vHp);
            vHt = _mm_max_epi8(vE, vHp);
            vHtM = _mm_blendv_epi8(vHpM, vEM, case1);
            vHtS = _mm_blendv_epi8(vHpS, vES, case1);
            vHtL = _mm_blendv_epi8(vHpL, vEL, case1);
            _mm_store_si128(pvE+i, vE);
            _mm_store_si128(pvEM+i, vEM);
            _mm_store_si128(pvES+i, vES);
            _mm_store_si128(pvEL+i, vEL);
            _mm_store_si128(pvH+i, vHp);
            _mm_store_si128(pvHM+i, vHpM);
            _mm_store_si128(pvHS+i, vHpS);
            _mm_store_si128(pvHL+i, vHpL);
            vHp = vH;
            vHpM = vHM;
            vHpS = vHS;
            vHpL = vHL;
        }

        /* pseudo prefix scan on F and H */
        vHt = _mm_slli_si128(vHt, 1);
        vHtM = _mm_slli_si128(vHtM, 1);
        vHtS = _mm_slli_si128(vHtS, 1);
        vHtL = _mm_slli_si128(vHtL, 1);
        vGapper = _mm_load_si128(pvGapper);
        vGapperL = _mm_load_si128(pvGapperL);
        vGapper = _mm_adds_epi8(vHt, vGapper);
        case1 = _mm_or_si128(
                _mm_cmpgt_epi8(vGapper, vF),
                _mm_cmpeq_epi8(vGapper, vF));
        vF = _mm_max_epi8(vF, vGapper);
        vFM = _mm_blendv_epi8(vFM, vHtM, case1);
        vFS = _mm_blendv_epi8(vFS, vHtS, case1);
        vFL = _mm_blendv_epi8(
                vFL,
                _mm_adds_epi8(vHtL, vGapperL),
                case1);
        for (i=0; i<segWidth-2; ++i) {
            vec128i vFt = _mm_slli_si128(vF, 1);
            vec128i vFtM = _mm_slli_si128(vFM, 1);
            vec128i vFtS = _mm_slli_si128(vFS, 1);
            vec128i vFtL = _mm_slli_si128(vFL, 1);
            vFt = _mm_adds_epi8(vFt, vSegLenXgap);
            case1 = _mm_or_si128(
                    _mm_cmpgt_epi8(vFt, vF),
                    _mm_cmpeq_epi8(vFt, vF));
            vF = _mm_max_epi8(vF, vFt);
            vFM = _mm_blendv_epi8(vFM, vFtM, case1);
            vFS = _mm_blendv_epi8(vFS, vFtS, case1);
            vFL = _mm_blendv_epi8(
                    vFL,
                    _mm_adds_epi8(vFtL, vSegLen),
                    case1);
        }

        /* calculate final H */
        vF = _mm_slli_si128(vF, 1);
        vFM = _mm_slli_si128(vFM, 1);
        vFS = _mm_slli_si128(vFS, 1);
        vFL = _mm_slli_si128(vFL, 1);
        vF = _mm_adds_epi8(vF, vNegInfFront);
        case1 = _mm_cmpgt_epi8(vF, vHt);
        vH = _mm_max_epi8(vF, vHt);
        vHM = _mm_blendv_epi8(vHtM, vFM, case1);
        vHS = _mm_blendv_epi8(vHtS, vFS, case1);
        vHL = _mm_blendv_epi8(vHtL, vFL, case1);
        for (i=0; i<segLen; ++i) {
            vHp = _mm_load_si128(pvH+i);
            vHpM = _mm_load_si128(pvHM+i);
            vHpS = _mm_load_si128(pvHS+i);
            vHpL = _mm_load_si128(pvHL+i);
            vE = _mm_load_si128(pvE+i);
            vEM = _mm_load_si128(pvEM+i);
            vES = _mm_load_si128(pvES+i);
            vEL = _mm_load_si128(pvEL+i);
            vF_opn = _mm_subs_epi8(vH, vGapO);
            vF_ext = _mm_subs_epi8(vF, vGapE);
            vF = _mm_max_epi8(vF_opn, vF_ext);
            case1 = _mm_cmpgt_epi8(vF_opn, vF_ext);
            vFM = _mm_blendv_epi8(vFM, vHM, case1);
            vFS = _mm_blendv_epi8(vFS, vHS, case1);
            vFL = _mm_blendv_epi8(vFL, vHL, case1);
            vFL = _mm_adds_epi8(vFL, vOne);
            vH = _mm_max_epi8(vHp, vE);
            vH = _mm_max_epi8(vH, vF);
            case1 = _mm_cmpeq_epi8(vH, vHp);
            case2 = _mm_cmpeq_epi8(vH, vF);
            vHM = _mm_blendv_epi8(
                    _mm_blendv_epi8(vEM, vFM, case2),
                    vHpM, case1);
            vHS = _mm_blendv_epi8(
                    _mm_blendv_epi8(vES, vFS, case2),
                    vHpS, case1);
            vHL = _mm_blendv_epi8(
                    _mm_blendv_epi8(vEL, vFL, case2),
                    vHpL, case1);
            _mm_store_si128(pvH+i, vH);
            _mm_store_si128(pvHM+i, vHM);
            _mm_store_si128(pvHS+i, vHS);
            _mm_store_si128(pvHL+i, vHL);
            vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vH);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vH);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
#endif
        } 
        /* extract vector containing last value from column */
        {
            vec128i cond_max;
            vH = _mm_load_si128(pvH + offset);
            vHM = _mm_load_si128(pvHM + offset);
            vHS = _mm_load_si128(pvHS + offset);
            vHL = _mm_load_si128(pvHL + offset);
            cond_max = _mm_cmpgt_epi8(vH, vMaxH);
            vMaxH = _mm_max_epi8(vMaxH, vH);
            vMaxM = _mm_blendv_epi8(vMaxM, vHM, cond_max);
            vMaxS = _mm_blendv_epi8(vMaxS, vHS, cond_max);
            vMaxL = _mm_blendv_epi8(vMaxL, vHL, cond_max);
            if (_mm_movemask_epi8(_mm_and_si128(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 1);
                vHM = _mm_slli_si128(vHM, 1);
                vHS = _mm_slli_si128(vHS, 1);
                vHL = _mm_slli_si128(vHL, 1);
            }
            result->stats->rowcols->score_row[j] = (int8_t) _mm_extract_epi8 (vH, 15);
            result->stats->rowcols->matches_row[j] = (int8_t) _mm_extract_epi8 (vHM, 15);
            result->stats->rowcols->similar_row[j] = (int8_t) _mm_extract_epi8 (vHS, 15);
            result->stats->rowcols->length_row[j] = (int8_t) _mm_extract_epi8 (vHL, 15);
#endif
        }
    }

    /* max last value from all columns */
    {
        for (k=0; k<position; ++k) {
            vMaxH = _mm_slli_si128(vMaxH, 1);
            vMaxM = _mm_slli_si128(vMaxM, 1);
            vMaxS = _mm_slli_si128(vMaxS, 1);
            vMaxL = _mm_slli_si128(vMaxL, 1);
        }
        score = (int8_t) _mm_extract_epi8(vMaxH, 15);
        matches = (int8_t) _mm_extract_epi8(vMaxM, 15);
        similar = (int8_t) _mm_extract_epi8(vMaxS, 15);
        length = (int8_t) _mm_extract_epi8(vMaxL, 15);
    }

    /* max of last column */
    {
        int8_t score_last;
        vMaxH = vNegLimit;

        for (i=0; i<segLen; ++i) {
            vec128i vH = _mm_load_si128(pvH + i);
#ifdef PARASAIL_ROWCOL
            vec128i vHM = _mm_load_si128(pvHM + i);
            vec128i vHS = _mm_load_si128(pvHS + i);
            vec128i vHL = _mm_load_si128(pvHL + i);
            arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
            arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
            arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
            arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
#endif
            vMaxH = _mm_max_epi8(vH, vMaxH);
        }

        /* max in vec */
        score_last = _mm_hmax_epi8(vMaxH);
        if (score_last > score || (score_last == score && end_ref == s2Len - 1)) {
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                int8_t *t = (int8_t*)pvH;
                int8_t *m = (int8_t*)pvHM;
                int8_t *s = (int8_t*)pvHS;
                int8_t *l = (int8_t*)pvHL;
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

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmplt_epi8(vSaturationCheckMin, vNegLimit),
            _mm_cmpgt_epi8(vSaturationCheckMax, vPosLimit)))) {
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
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;
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


