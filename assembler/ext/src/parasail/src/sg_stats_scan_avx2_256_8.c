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


#if HAVE_AVX2_MM256_INSERT_EPI8
#define _mm256_insert_epi8_rpl _mm256_insert_epi8
#else
static inline __m256i _mm256_insert_epi8_rpl(__m256i a, int8_t i, int imm) {
    __m256i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI8
#define _mm256_extract_epi8_rpl _mm256_extract_epi8
#else
static inline int8_t _mm256_extract_epi8_rpl(__m256i a, int imm) {
    __m256i_8_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_cmplt_epi8_rpl(a,b) _mm256_cmpgt_epi8(b,a)

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

static inline int8_t _mm256_hmax_epi8_rpl(__m256i a) {
    a = _mm256_max_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 4));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 2));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 1));
    return _mm256_extract_epi8_rpl(a, 31);
}


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*( 0*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  0);
    array[1LL*( 1*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  1);
    array[1LL*( 2*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  2);
    array[1LL*( 3*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  3);
    array[1LL*( 4*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  4);
    array[1LL*( 5*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  5);
    array[1LL*( 6*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  6);
    array[1LL*( 7*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  7);
    array[1LL*( 8*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  8);
    array[1LL*( 9*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  9);
    array[1LL*(10*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 10);
    array[1LL*(11*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 11);
    array[1LL*(12*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 12);
    array[1LL*(13*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 13);
    array[1LL*(14*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 14);
    array[1LL*(15*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 15);
    array[1LL*(16*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 16);
    array[1LL*(17*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 17);
    array[1LL*(18*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 18);
    array[1LL*(19*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 19);
    array[1LL*(20*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 20);
    array[1LL*(21*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 21);
    array[1LL*(22*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 22);
    array[1LL*(23*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 23);
    array[1LL*(24*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 24);
    array[1LL*(25*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 25);
    array[1LL*(26*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 26);
    array[1LL*(27*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 27);
    array[1LL*(28*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 28);
    array[1LL*(29*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 29);
    array[1LL*(30*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 30);
    array[1LL*(31*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 31);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  0);
    col[ 1*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  1);
    col[ 2*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  2);
    col[ 3*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  3);
    col[ 4*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  4);
    col[ 5*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  5);
    col[ 6*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  6);
    col[ 7*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  7);
    col[ 8*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  8);
    col[ 9*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  9);
    col[10*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 10);
    col[11*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 11);
    col[12*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 12);
    col[13*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 13);
    col[14*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 14);
    col[15*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 15);
    col[16*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 16);
    col[17*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 17);
    col[18*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 18);
    col[19*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 19);
    col[20*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 20);
    col[21*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 21);
    col[22*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 22);
    col[23*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 23);
    col[24*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 24);
    col[25*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 25);
    col[26*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 26);
    col[27*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 27);
    col[28*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 28);
    col[29*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 29);
    col[30*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 30);
    col[31*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 31);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_stats_table_scan_avx2_256_8
#define PNAME parasail_sg_stats_table_scan_profile_avx2_256_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_scan_avx2_256_8
#define PNAME parasail_sg_stats_rowcol_scan_profile_avx2_256_8
#else
#define FNAME parasail_sg_stats_scan_avx2_256_8
#define PNAME parasail_sg_stats_scan_profile_avx2_256_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_avx_256_8(s1, s1Len, matrix);
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
    const int32_t segWidth = 32; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict pvP  = (__m256i*)profile->profile8.score;
    __m256i* const restrict pvPm = (__m256i*)profile->profile8.matches;
    __m256i* const restrict pvPs = (__m256i*)profile->profile8.similar;
    __m256i* const restrict pvE  = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEM = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvES = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEL = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvH  = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHM = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHS = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHL = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHMax  = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHMMax = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHSMax = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHLMax = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvGapper = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvGapperL = parasail_memalign___m256i(32, segLen);
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    const int8_t NEG_LIMIT = (-open < matrix->min ?
        INT8_MIN + open : INT8_MIN - matrix->min) + 1;
    const int8_t POS_LIMIT = INT8_MAX - matrix->max - 1;
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi8(1);
    int8_t score = NEG_LIMIT;
    int8_t matches = 0;
    int8_t similar = 0;
    int8_t length = 0;
    __m256i vNegLimit = _mm256_set1_epi8(NEG_LIMIT);
    __m256i vPosLimit = _mm256_set1_epi8(POS_LIMIT);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
    __m256i vMaxH = vNegLimit;
    __m256i vMaxM = vNegLimit;
    __m256i vMaxS = vNegLimit;
    __m256i vMaxL = vNegLimit;
    __m256i vPosMask = _mm256_cmpeq_epi8(_mm256_set1_epi8(position),
            _mm256_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31));
    __m256i vNegInfFront = vZero;
    __m256i vSegLenXgap;
    __m256i vSegLen = _mm256_slli_si256_rpl(_mm256_set1_epi8(segLen), 1);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    vNegInfFront = _mm256_insert_epi8_rpl(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = _mm256_adds_epi8(vNegInfFront,
            _mm256_slli_si256_rpl(_mm256_set1_epi8(-segLen*gap), 1));

    parasail_memset___m256i(pvH, vZero, segLen);
    parasail_memset___m256i(pvHM, vZero, segLen);
    parasail_memset___m256i(pvHS, vZero, segLen);
    parasail_memset___m256i(pvHL, vZero, segLen);
    parasail_memset___m256i(pvE, vNegLimit, segLen);
    parasail_memset___m256i(pvEM, vZero, segLen);
    parasail_memset___m256i(pvES, vZero, segLen);
    parasail_memset___m256i(pvEL, vZero, segLen);
    {
        __m256i vGapper = _mm256_subs_epi8(vZero,vGapO);
        __m256i vGapperL = vOne;
        for (i=segLen-1; i>=0; --i) {
            _mm256_store_si256(pvGapper+i, vGapper);
            _mm256_store_si256(pvGapperL+i, vGapperL);
            vGapper = _mm256_subs_epi8(vGapper, vGapE);
            vGapperL = _mm256_adds_epi8(vGapperL, vOne);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        __m256i vE_ext;
        __m256i vE_opn;
        __m256i vEM;
        __m256i vES;
        __m256i vEL;
        __m256i vHt;
        __m256i vHtM;
        __m256i vHtS;
        __m256i vHtL;
        __m256i vF;
        __m256i vF_ext;
        __m256i vF_opn;
        __m256i vFM;
        __m256i vFS;
        __m256i vFL;
        __m256i vH;
        __m256i vHM;
        __m256i vHS;
        __m256i vHL;
        __m256i vHp;
        __m256i vHpM;
        __m256i vHpS;
        __m256i vHpL;
        __m256i *pvW;
        __m256i vW;
        __m256i *pvWM;
        __m256i vWM;
        __m256i *pvWS;
        __m256i vWS;
        __m256i case1;
        __m256i case2;
        __m256i vGapper;
        __m256i vGapperL;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = _mm256_load_si256(pvH+(segLen-1));
        vHpM = _mm256_load_si256(pvHM+(segLen-1));
        vHpS = _mm256_load_si256(pvHS+(segLen-1));
        vHpL = _mm256_load_si256(pvHL+(segLen-1));
        vHp = _mm256_slli_si256_rpl(vHp, 1);
        vHpM = _mm256_slli_si256_rpl(vHpM, 1);
        vHpS = _mm256_slli_si256_rpl(vHpS, 1);
        vHpL = _mm256_slli_si256_rpl(vHpL, 1);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWM= pvPm+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWS= pvPs+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = _mm256_subs_epi8(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        vFM = vZero;
        vFS = vZero;
        vFL = vZero;
        for (i=0; i<segLen; ++i) {
            vH = _mm256_load_si256(pvH+i);
            vHM= _mm256_load_si256(pvHM+i);
            vHS= _mm256_load_si256(pvHS+i);
            vHL= _mm256_load_si256(pvHL+i);
            vE = _mm256_load_si256(pvE+i);
            vEM= _mm256_load_si256(pvEM+i);
            vES= _mm256_load_si256(pvES+i);
            vEL= _mm256_load_si256(pvEL+i);
            vW = _mm256_load_si256(pvW+i);
            vWM = _mm256_load_si256(pvWM+i);
            vWS = _mm256_load_si256(pvWS+i);
            vGapper = _mm256_load_si256(pvGapper+i);
            vGapperL = _mm256_load_si256(pvGapperL+i);
            vE_opn = _mm256_subs_epi8(vH, vGapO);
            vE_ext = _mm256_subs_epi8(vE, vGapE);
            case1 = _mm256_cmpgt_epi8(vE_opn, vE_ext);
            vE = _mm256_max_epi8(vE_opn, vE_ext);
            vEM = _mm256_blendv_epi8(vEM, vHM, case1);
            vES = _mm256_blendv_epi8(vES, vHS, case1);
            vEL = _mm256_blendv_epi8(vEL, vHL, case1);
            vEL = _mm256_adds_epi8(vEL, vOne);
            vGapper = _mm256_adds_epi8(vHt, vGapper);
            case1 = _mm256_or_si256(
                    _mm256_cmpgt_epi8(vF, vGapper),
                    _mm256_cmpeq_epi8(vF, vGapper));
            vF = _mm256_max_epi8(vF, vGapper);
            vFM = _mm256_blendv_epi8(vHtM, vFM, case1);
            vFS = _mm256_blendv_epi8(vHtS, vFS, case1);
            vFL = _mm256_blendv_epi8(
                    _mm256_adds_epi8(vHtL, vGapperL),
                    vFL, case1);
            vHp = _mm256_adds_epi8(vHp, vW);
            vHpM = _mm256_adds_epi8(vHpM, vWM);
            vHpS = _mm256_adds_epi8(vHpS, vWS);
            vHpL = _mm256_adds_epi8(vHpL, vOne);
            case1 = _mm256_cmpgt_epi8(vE, vHp);
            vHt = _mm256_max_epi8(vE, vHp);
            vHtM = _mm256_blendv_epi8(vHpM, vEM, case1);
            vHtS = _mm256_blendv_epi8(vHpS, vES, case1);
            vHtL = _mm256_blendv_epi8(vHpL, vEL, case1);
            _mm256_store_si256(pvE+i, vE);
            _mm256_store_si256(pvEM+i, vEM);
            _mm256_store_si256(pvES+i, vES);
            _mm256_store_si256(pvEL+i, vEL);
            _mm256_store_si256(pvH+i, vHp);
            _mm256_store_si256(pvHM+i, vHpM);
            _mm256_store_si256(pvHS+i, vHpS);
            _mm256_store_si256(pvHL+i, vHpL);
            vHp = vH;
            vHpM = vHM;
            vHpS = vHS;
            vHpL = vHL;
        }

        /* pseudo prefix scan on F and H */
        vHt = _mm256_slli_si256_rpl(vHt, 1);
        vHtM = _mm256_slli_si256_rpl(vHtM, 1);
        vHtS = _mm256_slli_si256_rpl(vHtS, 1);
        vHtL = _mm256_slli_si256_rpl(vHtL, 1);
        vGapper = _mm256_load_si256(pvGapper);
        vGapperL = _mm256_load_si256(pvGapperL);
        vGapper = _mm256_adds_epi8(vHt, vGapper);
        case1 = _mm256_or_si256(
                _mm256_cmpgt_epi8(vGapper, vF),
                _mm256_cmpeq_epi8(vGapper, vF));
        vF = _mm256_max_epi8(vF, vGapper);
        vFM = _mm256_blendv_epi8(vFM, vHtM, case1);
        vFS = _mm256_blendv_epi8(vFS, vHtS, case1);
        vFL = _mm256_blendv_epi8(
                vFL,
                _mm256_adds_epi8(vHtL, vGapperL),
                case1);
        for (i=0; i<segWidth-2; ++i) {
            __m256i vFt = _mm256_slli_si256_rpl(vF, 1);
            __m256i vFtM = _mm256_slli_si256_rpl(vFM, 1);
            __m256i vFtS = _mm256_slli_si256_rpl(vFS, 1);
            __m256i vFtL = _mm256_slli_si256_rpl(vFL, 1);
            vFt = _mm256_adds_epi8(vFt, vSegLenXgap);
            case1 = _mm256_or_si256(
                    _mm256_cmpgt_epi8(vFt, vF),
                    _mm256_cmpeq_epi8(vFt, vF));
            vF = _mm256_max_epi8(vF, vFt);
            vFM = _mm256_blendv_epi8(vFM, vFtM, case1);
            vFS = _mm256_blendv_epi8(vFS, vFtS, case1);
            vFL = _mm256_blendv_epi8(
                    vFL,
                    _mm256_adds_epi8(vFtL, vSegLen),
                    case1);
        }

        /* calculate final H */
        vF = _mm256_slli_si256_rpl(vF, 1);
        vFM = _mm256_slli_si256_rpl(vFM, 1);
        vFS = _mm256_slli_si256_rpl(vFS, 1);
        vFL = _mm256_slli_si256_rpl(vFL, 1);
        vF = _mm256_adds_epi8(vF, vNegInfFront);
        case1 = _mm256_cmpgt_epi8(vF, vHt);
        vH = _mm256_max_epi8(vF, vHt);
        vHM = _mm256_blendv_epi8(vHtM, vFM, case1);
        vHS = _mm256_blendv_epi8(vHtS, vFS, case1);
        vHL = _mm256_blendv_epi8(vHtL, vFL, case1);
        for (i=0; i<segLen; ++i) {
            vHp = _mm256_load_si256(pvH+i);
            vHpM = _mm256_load_si256(pvHM+i);
            vHpS = _mm256_load_si256(pvHS+i);
            vHpL = _mm256_load_si256(pvHL+i);
            vE = _mm256_load_si256(pvE+i);
            vEM = _mm256_load_si256(pvEM+i);
            vES = _mm256_load_si256(pvES+i);
            vEL = _mm256_load_si256(pvEL+i);
            vF_opn = _mm256_subs_epi8(vH, vGapO);
            vF_ext = _mm256_subs_epi8(vF, vGapE);
            vF = _mm256_max_epi8(vF_opn, vF_ext);
            case1 = _mm256_cmpgt_epi8(vF_opn, vF_ext);
            vFM = _mm256_blendv_epi8(vFM, vHM, case1);
            vFS = _mm256_blendv_epi8(vFS, vHS, case1);
            vFL = _mm256_blendv_epi8(vFL, vHL, case1);
            vFL = _mm256_adds_epi8(vFL, vOne);
            vH = _mm256_max_epi8(vHp, vE);
            vH = _mm256_max_epi8(vH, vF);
            case1 = _mm256_cmpeq_epi8(vH, vHp);
            case2 = _mm256_cmpeq_epi8(vH, vF);
            vHM = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEM, vFM, case2),
                    vHpM, case1);
            vHS = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vES, vFS, case2),
                    vHpS, case1);
            vHL = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEL, vFL, case2),
                    vHpL, case1);
            _mm256_store_si256(pvH+i, vH);
            _mm256_store_si256(pvHM+i, vHM);
            _mm256_store_si256(pvHS+i, vHS);
            _mm256_store_si256(pvHL+i, vHL);
            vSaturationCheckMin = _mm256_min_epi8(vSaturationCheckMin, vH);
            vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vH);
            vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
            arr_store_si256(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si256(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si256(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
#endif
        } 
        /* extract vector containing last value from column */
        {
            __m256i cond_max;
            vH = _mm256_load_si256(pvH + offset);
            vHM = _mm256_load_si256(pvHM + offset);
            vHS = _mm256_load_si256(pvHS + offset);
            vHL = _mm256_load_si256(pvHL + offset);
            cond_max = _mm256_cmpgt_epi8(vH, vMaxH);
            vMaxH = _mm256_max_epi8(vMaxH, vH);
            vMaxM = _mm256_blendv_epi8(vMaxM, vHM, cond_max);
            vMaxS = _mm256_blendv_epi8(vMaxS, vHS, cond_max);
            vMaxL = _mm256_blendv_epi8(vMaxL, vHL, cond_max);
            if (_mm256_movemask_epi8(_mm256_and_si256(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 1);
                vHM = _mm256_slli_si256_rpl(vHM, 1);
                vHS = _mm256_slli_si256_rpl(vHS, 1);
                vHL = _mm256_slli_si256_rpl(vHL, 1);
            }
            result->stats->rowcols->score_row[j] = (int8_t) _mm256_extract_epi8_rpl (vH, 31);
            result->stats->rowcols->matches_row[j] = (int8_t) _mm256_extract_epi8_rpl (vHM, 31);
            result->stats->rowcols->similar_row[j] = (int8_t) _mm256_extract_epi8_rpl (vHS, 31);
            result->stats->rowcols->length_row[j] = (int8_t) _mm256_extract_epi8_rpl (vHL, 31);
#endif
        }
    }

    /* max last value from all columns */
    {
        for (k=0; k<position; ++k) {
            vMaxH = _mm256_slli_si256_rpl(vMaxH, 1);
            vMaxM = _mm256_slli_si256_rpl(vMaxM, 1);
            vMaxS = _mm256_slli_si256_rpl(vMaxS, 1);
            vMaxL = _mm256_slli_si256_rpl(vMaxL, 1);
        }
        score = (int8_t) _mm256_extract_epi8_rpl(vMaxH, 31);
        matches = (int8_t) _mm256_extract_epi8_rpl(vMaxM, 31);
        similar = (int8_t) _mm256_extract_epi8_rpl(vMaxS, 31);
        length = (int8_t) _mm256_extract_epi8_rpl(vMaxL, 31);
    }

    /* max of last column */
    {
        int8_t score_last;
        vMaxH = vNegLimit;

        for (i=0; i<segLen; ++i) {
            __m256i vH = _mm256_load_si256(pvH + i);
#ifdef PARASAIL_ROWCOL
            __m256i vHM = _mm256_load_si256(pvHM + i);
            __m256i vHS = _mm256_load_si256(pvHS + i);
            __m256i vHL = _mm256_load_si256(pvHL + i);
            arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
            arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
            arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
            arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
#endif
            vMaxH = _mm256_max_epi8(vH, vMaxH);
        }

        /* max in vec */
        score_last = _mm256_hmax_epi8_rpl(vMaxH);
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

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmplt_epi8_rpl(vSaturationCheckMin, vNegLimit),
            _mm256_cmpgt_epi8(vSaturationCheckMax, vPosLimit)))) {
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
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_32;
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


