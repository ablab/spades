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


static inline __m128i _mm_blendv_epi8_rpl(__m128i a, __m128i b, __m128i mask) {
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
}

static inline __m128i _mm_insert_epi8_rpl(__m128i a, int8_t i, const int imm) {
    __m128i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
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

static inline __m128i _mm_min_epi8_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi8(b, a);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
}


#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*( 0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  0);
    array[1LL*( 1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  1);
    array[1LL*( 2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  2);
    array[1LL*( 3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  3);
    array[1LL*( 4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  4);
    array[1LL*( 5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  5);
    array[1LL*( 6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  6);
    array[1LL*( 7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  7);
    array[1LL*( 8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  8);
    array[1LL*( 9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  9);
    array[1LL*(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 10);
    array[1LL*(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 11);
    array[1LL*(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 12);
    array[1LL*(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 13);
    array[1LL*(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 14);
    array[1LL*(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 15);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  0);
    col[ 1*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  1);
    col[ 2*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  2);
    col[ 3*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  3);
    col[ 4*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  4);
    col[ 5*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  5);
    col[ 6*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  6);
    col[ 7*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  7);
    col[ 8*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  8);
    col[ 9*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  9);
    col[10*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 10);
    col[11*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 11);
    col[12*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 12);
    col[13*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 13);
    col[14*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 14);
    col[15*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_stats_table_scan_sse2_128_8
#define PNAME parasail_nw_stats_table_scan_profile_sse2_128_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_stats_rowcol_scan_sse2_128_8
#define PNAME parasail_nw_stats_rowcol_scan_profile_sse2_128_8
#else
#define FNAME parasail_nw_stats_scan_sse2_128_8
#define PNAME parasail_nw_stats_scan_profile_sse2_128_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_sse_128_8(s1, s1Len, matrix);
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
    __m128i* const restrict pvP  = (__m128i*)profile->profile8.score;
    __m128i* const restrict pvPm = (__m128i*)profile->profile8.matches;
    __m128i* const restrict pvPs = (__m128i*)profile->profile8.similar;
    __m128i* const restrict pvE  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvEM = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvES = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvEL = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvH  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHM = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHS = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHL = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHMax  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHMMax = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHSMax = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHLMax = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvGapper = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvGapperL = parasail_memalign___m128i(16, segLen);
    int8_t* const restrict boundary = parasail_memalign_int8_t(16, s2Len+1);
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    const int8_t NEG_LIMIT = (-open < matrix->min ?
        INT8_MIN + open : INT8_MIN - matrix->min) + 1;
    const int8_t POS_LIMIT = INT8_MAX - matrix->max - 1;
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi8(1);
    int8_t score = NEG_LIMIT;
    int8_t matches = 0;
    int8_t similar = 0;
    int8_t length = 0;
    __m128i vNegLimit = _mm_set1_epi8(NEG_LIMIT);
    __m128i vPosLimit = _mm_set1_epi8(POS_LIMIT);
    __m128i vSaturationCheckMin = vPosLimit;
    __m128i vSaturationCheckMax = vNegLimit;
    __m128i vNegInfFront = vZero;
    __m128i vSegLenXgap;
    __m128i vSegLen = _mm_slli_si128(_mm_set1_epi8(segLen), 1);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    vNegInfFront = _mm_insert_epi8_rpl(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = _mm_adds_epi8(vNegInfFront,
            _mm_slli_si128(_mm_set1_epi8(-segLen*gap), 1));

    parasail_memset___m128i(pvHM, vZero, segLen);
    parasail_memset___m128i(pvHS, vZero, segLen);
    parasail_memset___m128i(pvHL, vZero, segLen);
    parasail_memset___m128i(pvE, vNegLimit, segLen);
    parasail_memset___m128i(pvEM, vZero, segLen);
    parasail_memset___m128i(pvES, vZero, segLen);
    parasail_memset___m128i(pvEL, vZero, segLen);
    {
        __m128i vGapper = _mm_subs_epi8(vZero,vGapO);
        __m128i vGapperL = vOne;
        for (i=segLen-1; i>=0; --i) {
            _mm_store_si128(pvGapper+i, vGapper);
            _mm_store_si128(pvGapperL+i, vGapperL);
            vGapper = _mm_subs_epi8(vGapper, vGapE);
            vGapperL = _mm_adds_epi8(vGapperL, vOne);
        }
    }

    /* initialize H */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            __m128i_8_t h;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
            }
            _mm_store_si128(&pvH[index], h.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT8_MIN ? INT8_MIN : tmp;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vE_ext;
        __m128i vE_opn;
        __m128i vEM;
        __m128i vES;
        __m128i vEL;
        __m128i vHt;
        __m128i vHtM;
        __m128i vHtS;
        __m128i vHtL;
        __m128i vF;
        __m128i vF_ext;
        __m128i vF_opn;
        __m128i vFM;
        __m128i vFS;
        __m128i vFL;
        __m128i vH;
        __m128i vHM;
        __m128i vHS;
        __m128i vHL;
        __m128i vHp;
        __m128i vHpM;
        __m128i vHpS;
        __m128i vHpL;
        __m128i *pvW;
        __m128i vW;
        __m128i *pvWM;
        __m128i vWM;
        __m128i *pvWS;
        __m128i vWS;
        __m128i case1;
        __m128i case2;
        __m128i vGapper;
        __m128i vGapperL;

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
        vHp = _mm_insert_epi8_rpl(vHp, boundary[j], 0);
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
            vE = _mm_max_epi8_rpl(vE_opn, vE_ext);
            vEM = _mm_blendv_epi8_rpl(vEM, vHM, case1);
            vES = _mm_blendv_epi8_rpl(vES, vHS, case1);
            vEL = _mm_blendv_epi8_rpl(vEL, vHL, case1);
            vEL = _mm_adds_epi8(vEL, vOne);
            vGapper = _mm_adds_epi8(vHt, vGapper);
            case1 = _mm_or_si128(
                    _mm_cmpgt_epi8(vF, vGapper),
                    _mm_cmpeq_epi8(vF, vGapper));
            vF = _mm_max_epi8_rpl(vF, vGapper);
            vFM = _mm_blendv_epi8_rpl(vHtM, vFM, case1);
            vFS = _mm_blendv_epi8_rpl(vHtS, vFS, case1);
            vFL = _mm_blendv_epi8_rpl(
                    _mm_adds_epi8(vHtL, vGapperL),
                    vFL, case1);
            vHp = _mm_adds_epi8(vHp, vW);
            vHpM = _mm_adds_epi8(vHpM, vWM);
            vHpS = _mm_adds_epi8(vHpS, vWS);
            vHpL = _mm_adds_epi8(vHpL, vOne);
            case1 = _mm_cmpgt_epi8(vE, vHp);
            vHt = _mm_max_epi8_rpl(vE, vHp);
            vHtM = _mm_blendv_epi8_rpl(vHpM, vEM, case1);
            vHtS = _mm_blendv_epi8_rpl(vHpS, vES, case1);
            vHtL = _mm_blendv_epi8_rpl(vHpL, vEL, case1);
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
        vHt = _mm_insert_epi8_rpl(vHt, boundary[j+1], 0);
        vGapper = _mm_load_si128(pvGapper);
        vGapperL = _mm_load_si128(pvGapperL);
        vGapper = _mm_adds_epi8(vHt, vGapper);
        case1 = _mm_or_si128(
                _mm_cmpgt_epi8(vGapper, vF),
                _mm_cmpeq_epi8(vGapper, vF));
        vF = _mm_max_epi8_rpl(vF, vGapper);
        vFM = _mm_blendv_epi8_rpl(vFM, vHtM, case1);
        vFS = _mm_blendv_epi8_rpl(vFS, vHtS, case1);
        vFL = _mm_blendv_epi8_rpl(
                vFL,
                _mm_adds_epi8(vHtL, vGapperL),
                case1);
        for (i=0; i<segWidth-2; ++i) {
            __m128i vFt = _mm_slli_si128(vF, 1);
            __m128i vFtM = _mm_slli_si128(vFM, 1);
            __m128i vFtS = _mm_slli_si128(vFS, 1);
            __m128i vFtL = _mm_slli_si128(vFL, 1);
            vFt = _mm_adds_epi8(vFt, vSegLenXgap);
            case1 = _mm_or_si128(
                    _mm_cmpgt_epi8(vFt, vF),
                    _mm_cmpeq_epi8(vFt, vF));
            vF = _mm_max_epi8_rpl(vF, vFt);
            vFM = _mm_blendv_epi8_rpl(vFM, vFtM, case1);
            vFS = _mm_blendv_epi8_rpl(vFS, vFtS, case1);
            vFL = _mm_blendv_epi8_rpl(
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
        vH = _mm_max_epi8_rpl(vF, vHt);
        vHM = _mm_blendv_epi8_rpl(vHtM, vFM, case1);
        vHS = _mm_blendv_epi8_rpl(vHtS, vFS, case1);
        vHL = _mm_blendv_epi8_rpl(vHtL, vFL, case1);
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
            vF = _mm_max_epi8_rpl(vF_opn, vF_ext);
            case1 = _mm_cmpgt_epi8(vF_opn, vF_ext);
            vFM = _mm_blendv_epi8_rpl(vFM, vHM, case1);
            vFS = _mm_blendv_epi8_rpl(vFS, vHS, case1);
            vFL = _mm_blendv_epi8_rpl(vFL, vHL, case1);
            vFL = _mm_adds_epi8(vFL, vOne);
            vH = _mm_max_epi8_rpl(vHp, vE);
            vH = _mm_max_epi8_rpl(vH, vF);
            case1 = _mm_cmpeq_epi8(vH, vHp);
            case2 = _mm_cmpeq_epi8(vH, vF);
            vHM = _mm_blendv_epi8_rpl(
                    _mm_blendv_epi8_rpl(vEM, vFM, case2),
                    vHpM, case1);
            vHS = _mm_blendv_epi8_rpl(
                    _mm_blendv_epi8_rpl(vES, vFS, case2),
                    vHpS, case1);
            vHL = _mm_blendv_epi8_rpl(
                    _mm_blendv_epi8_rpl(vEL, vFL, case2),
                    vHpL, case1);
            _mm_store_si128(pvH+i, vH);
            _mm_store_si128(pvHM+i, vHM);
            _mm_store_si128(pvHS+i, vHS);
            _mm_store_si128(pvHL+i, vHL);
            vSaturationCheckMin = _mm_min_epi8_rpl(vSaturationCheckMin, vH);
            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vH);
            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
#endif
        } 

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm_load_si128(pvH + offset);
            vHM = _mm_load_si128(pvHM + offset);
            vHS = _mm_load_si128(pvHS + offset);
            vHL = _mm_load_si128(pvHL + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 1);
                vHM = _mm_slli_si128(vHM, 1);
                vHS = _mm_slli_si128(vHS, 1);
                vHL = _mm_slli_si128(vHL, 1);
            }
            result->stats->rowcols->score_row[j] = (int8_t) _mm_extract_epi8_rpl (vH, 15);
            result->stats->rowcols->matches_row[j] = (int8_t) _mm_extract_epi8_rpl (vHM, 15);
            result->stats->rowcols->similar_row[j] = (int8_t) _mm_extract_epi8_rpl (vHS, 15);
            result->stats->rowcols->length_row[j] = (int8_t) _mm_extract_epi8_rpl (vHL, 15);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m128i vH = _mm_load_si128(pvH+i);
        __m128i vHM = _mm_load_si128(pvHM+i);
        __m128i vHS = _mm_load_si128(pvHS+i);
        __m128i vHL = _mm_load_si128(pvHL+i);
        arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
        arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
        arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
        arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        __m128i vH = _mm_load_si128(pvH + offset);
        __m128i vHM = _mm_load_si128(pvHM + offset);
        __m128i vHS = _mm_load_si128(pvHS + offset);
        __m128i vHL = _mm_load_si128(pvHL + offset);
        for (k=0; k<position; ++k) {
            vH = _mm_slli_si128(vH, 1);
            vHM = _mm_slli_si128(vHM, 1);
            vHS = _mm_slli_si128(vHS, 1);
            vHL = _mm_slli_si128(vHL, 1);
        }
        score = (int8_t) _mm_extract_epi8_rpl (vH, 15);
        matches = (int8_t) _mm_extract_epi8_rpl (vHM, 15);
        similar = (int8_t) _mm_extract_epi8_rpl (vHS, 15);
        length = (int8_t) _mm_extract_epi8_rpl (vHL, 15);
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
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(boundary);
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


