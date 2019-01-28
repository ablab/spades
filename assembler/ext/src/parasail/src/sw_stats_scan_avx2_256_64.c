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
#include <string.h>

#include <immintrin.h>

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_avx.h"


#if HAVE_AVX2_MM256_INSERT_EPI64
#define _mm256_insert_epi64_rpl _mm256_insert_epi64
#else
static inline __m256i _mm256_insert_epi64_rpl(__m256i a, int64_t i, int imm) {
    __m256i_64_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_SET1_EPI64X
#define _mm256_set1_epi64x_rpl _mm256_set1_epi64x
#else
static inline __m256i _mm256_set1_epi64x_rpl(int64_t i) {
    __m256i_64_t A;
    A.v[0] = i;
    A.v[1] = i;
    A.v[2] = i;
    A.v[3] = i;
    return A.m;
}
#endif

static inline __m256i _mm256_max_epi64_rpl(__m256i a, __m256i b) {
    __m256i_64_t A;
    __m256i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    A.v[2] = (A.v[2]>B.v[2]) ? A.v[2] : B.v[2];
    A.v[3] = (A.v[3]>B.v[3]) ? A.v[3] : B.v[3];
    return A.m;
}

#if HAVE_AVX2_MM256_EXTRACT_EPI64
#define _mm256_extract_epi64_rpl _mm256_extract_epi64
#else
static inline int64_t _mm256_extract_epi64_rpl(__m256i a, int imm) {
    __m256i_64_t A;
    A.m = a;
    return A.v[imm];
}
#endif

static inline __m256i _mm256_min_epi64_rpl(__m256i a, __m256i b) {
    __m256i_64_t A;
    __m256i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]<B.v[1]) ? A.v[1] : B.v[1];
    A.v[2] = (A.v[2]<B.v[2]) ? A.v[2] : B.v[2];
    A.v[3] = (A.v[3]<B.v[3]) ? A.v[3] : B.v[3];
    return A.m;
}

#define _mm256_cmplt_epi64_rpl(a,b) _mm256_cmpgt_epi64(b,a)

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

static inline int64_t _mm256_hmax_epi64_rpl(__m256i a) {
    a = _mm256_max_epi64_rpl(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi64_rpl(a, _mm256_slli_si256(a, 8));
    return _mm256_extract_epi64_rpl(a, 3);
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
    array[1LL*(0*seglen+t)*dlen + d] = (int64_t)_mm256_extract_epi64_rpl(vH, 0);
    array[1LL*(1*seglen+t)*dlen + d] = (int64_t)_mm256_extract_epi64_rpl(vH, 1);
    array[1LL*(2*seglen+t)*dlen + d] = (int64_t)_mm256_extract_epi64_rpl(vH, 2);
    array[1LL*(3*seglen+t)*dlen + d] = (int64_t)_mm256_extract_epi64_rpl(vH, 3);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int64_t)_mm256_extract_epi64_rpl(vH, 0);
    col[1*seglen+t] = (int64_t)_mm256_extract_epi64_rpl(vH, 1);
    col[2*seglen+t] = (int64_t)_mm256_extract_epi64_rpl(vH, 2);
    col[3*seglen+t] = (int64_t)_mm256_extract_epi64_rpl(vH, 3);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_scan_avx2_256_64
#define PNAME parasail_sw_stats_table_scan_profile_avx2_256_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_scan_avx2_256_64
#define PNAME parasail_sw_stats_rowcol_scan_profile_avx2_256_64
#else
#define FNAME parasail_sw_stats_scan_avx2_256_64
#define PNAME parasail_sw_stats_scan_profile_avx2_256_64
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_avx_256_64(s1, s1Len, matrix);
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
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 4; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m256i* const restrict pvP  = (__m256i*)profile->profile64.score;
    __m256i* const restrict pvPm = (__m256i*)profile->profile64.matches;
    __m256i* const restrict pvPs = (__m256i*)profile->profile64.similar;
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
    __m256i vGapO = _mm256_set1_epi64x_rpl(open);
    __m256i vGapE = _mm256_set1_epi64x_rpl(gap);
    const int64_t NEG_LIMIT = (-open < matrix->min ?
        INT64_MIN + open : INT64_MIN - matrix->min) + 1;
    const int64_t POS_LIMIT = INT64_MAX - matrix->max - 1;
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi64x_rpl(1);
    int64_t score = NEG_LIMIT;
    int64_t matches = 0;
    int64_t similar = 0;
    int64_t length = 0;
    __m256i vNegLimit = _mm256_set1_epi64x_rpl(NEG_LIMIT);
    __m256i vPosLimit = _mm256_set1_epi64x_rpl(POS_LIMIT);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
    __m256i vMaxH = vNegLimit;
    __m256i vMaxM = vNegLimit;
    __m256i vMaxS = vNegLimit;
    __m256i vMaxL = vNegLimit;
    __m256i vMaxHUnit = vNegLimit;
    __m256i vNegInfFront = vZero;
    __m256i vSegLenXgap;
    __m256i vSegLen = _mm256_slli_si256_rpl(_mm256_set1_epi64x_rpl(segLen), 8);
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

    vNegInfFront = _mm256_insert_epi64_rpl(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = _mm256_add_epi64(vNegInfFront,
            _mm256_slli_si256_rpl(_mm256_set1_epi64x_rpl(-segLen*gap), 8));

    parasail_memset___m256i(pvH, vZero, segLen);
    parasail_memset___m256i(pvHM, vZero, segLen);
    parasail_memset___m256i(pvHS, vZero, segLen);
    parasail_memset___m256i(pvHL, vZero, segLen);
    parasail_memset___m256i(pvE, vNegLimit, segLen);
    parasail_memset___m256i(pvEM, vZero, segLen);
    parasail_memset___m256i(pvES, vZero, segLen);
    parasail_memset___m256i(pvEL, vZero, segLen);
    {
        __m256i vGapper = _mm256_sub_epi64(vZero,vGapO);
        __m256i vGapperL = vOne;
        for (i=segLen-1; i>=0; --i) {
            _mm256_store_si256(pvGapper+i, vGapper);
            _mm256_store_si256(pvGapperL+i, vGapperL);
            vGapper = _mm256_sub_epi64(vGapper, vGapE);
            vGapperL = _mm256_add_epi64(vGapperL, vOne);
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
        __m256i case0;
        __m256i vGapper;
        __m256i vGapperL;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = _mm256_load_si256(pvH+(segLen-1));
        vHpM = _mm256_load_si256(pvHM+(segLen-1));
        vHpS = _mm256_load_si256(pvHS+(segLen-1));
        vHpL = _mm256_load_si256(pvHL+(segLen-1));
        vHp = _mm256_slli_si256_rpl(vHp, 8);
        vHpM = _mm256_slli_si256_rpl(vHpM, 8);
        vHpS = _mm256_slli_si256_rpl(vHpS, 8);
        vHpL = _mm256_slli_si256_rpl(vHpL, 8);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWM= pvPm+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWS= pvPs+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vZero;
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
            vE_opn = _mm256_sub_epi64(vH, vGapO);
            vE_ext = _mm256_sub_epi64(vE, vGapE);
            case1 = _mm256_cmpgt_epi64(vE_opn, vE_ext);
            vE = _mm256_max_epi64_rpl(vE_opn, vE_ext);
            vEM = _mm256_blendv_epi8(vEM, vHM, case1);
            vES = _mm256_blendv_epi8(vES, vHS, case1);
            vEL = _mm256_blendv_epi8(vEL, vHL, case1);
            vEL = _mm256_add_epi64(vEL, vOne);
            vGapper = _mm256_add_epi64(vHt, vGapper);
            case1 = _mm256_or_si256(
                    _mm256_cmpgt_epi64(vF, vGapper),
                    _mm256_cmpeq_epi64(vF, vGapper));
            vF = _mm256_max_epi64_rpl(vF, vGapper);
            vFM = _mm256_blendv_epi8(vHtM, vFM, case1);
            vFS = _mm256_blendv_epi8(vHtS, vFS, case1);
            vFL = _mm256_blendv_epi8(
                    _mm256_add_epi64(vHtL, vGapperL),
                    vFL, case1);
            vHp = _mm256_add_epi64(vHp, vW);
            vHpM = _mm256_add_epi64(vHpM, vWM);
            vHpS = _mm256_add_epi64(vHpS, vWS);
            vHpL = _mm256_add_epi64(vHpL, vOne);
            case1 = _mm256_cmpgt_epi64(vE, vHp);
            vHt = _mm256_max_epi64_rpl(vE, vHp);
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
        vHt = _mm256_slli_si256_rpl(vHt, 8);
        vHtM = _mm256_slli_si256_rpl(vHtM, 8);
        vHtS = _mm256_slli_si256_rpl(vHtS, 8);
        vHtL = _mm256_slli_si256_rpl(vHtL, 8);
        vGapper = _mm256_load_si256(pvGapper);
        vGapperL = _mm256_load_si256(pvGapperL);
        vGapper = _mm256_add_epi64(vHt, vGapper);
        case1 = _mm256_or_si256(
                _mm256_cmpgt_epi64(vGapper, vF),
                _mm256_cmpeq_epi64(vGapper, vF));
        vF = _mm256_max_epi64_rpl(vF, vGapper);
        vFM = _mm256_blendv_epi8(vFM, vHtM, case1);
        vFS = _mm256_blendv_epi8(vFS, vHtS, case1);
        vFL = _mm256_blendv_epi8(
                vFL,
                _mm256_add_epi64(vHtL, vGapperL),
                case1);
        for (i=0; i<segWidth-2; ++i) {
            __m256i vFt = _mm256_slli_si256_rpl(vF, 8);
            __m256i vFtM = _mm256_slli_si256_rpl(vFM, 8);
            __m256i vFtS = _mm256_slli_si256_rpl(vFS, 8);
            __m256i vFtL = _mm256_slli_si256_rpl(vFL, 8);
            vFt = _mm256_add_epi64(vFt, vSegLenXgap);
            case1 = _mm256_or_si256(
                    _mm256_cmpgt_epi64(vFt, vF),
                    _mm256_cmpeq_epi64(vFt, vF));
            vF = _mm256_max_epi64_rpl(vF, vFt);
            vFM = _mm256_blendv_epi8(vFM, vFtM, case1);
            vFS = _mm256_blendv_epi8(vFS, vFtS, case1);
            vFL = _mm256_blendv_epi8(
                    vFL,
                    _mm256_add_epi64(vFtL, vSegLen),
                    case1);
        }

        /* calculate final H */
        vF = _mm256_slli_si256_rpl(vF, 8);
        vFM = _mm256_slli_si256_rpl(vFM, 8);
        vFS = _mm256_slli_si256_rpl(vFS, 8);
        vFL = _mm256_slli_si256_rpl(vFL, 8);
        vF = _mm256_add_epi64(vF, vNegInfFront);
        case1 = _mm256_cmpgt_epi64(vF, vHt);
        vH = _mm256_max_epi64_rpl(vF, vHt);
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
            vF_opn = _mm256_sub_epi64(vH, vGapO);
            vF_ext = _mm256_sub_epi64(vF, vGapE);
            vF = _mm256_max_epi64_rpl(vF_opn, vF_ext);
            case1 = _mm256_cmpgt_epi64(vF_opn, vF_ext);
            vFM = _mm256_blendv_epi8(vFM, vHM, case1);
            vFS = _mm256_blendv_epi8(vFS, vHS, case1);
            vFL = _mm256_blendv_epi8(vFL, vHL, case1);
            vFL = _mm256_add_epi64(vFL, vOne);
            vH = _mm256_max_epi64_rpl(vHp, vE);
            vH = _mm256_max_epi64_rpl(vH, vF);
            vH = _mm256_max_epi64_rpl(vH, vZero);
            case1 = _mm256_cmpeq_epi64(vH, vHp);
            case2 = _mm256_cmpeq_epi64(vH, vF);
            case0 = _mm256_cmpeq_epi64(vH, vZero);
            vHM = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEM, vFM, case2),
                    vHpM, case1);
            vHM = _mm256_blendv_epi8(vHM, vZero, case0);
            vHS = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vES, vFS, case2),
                    vHpS, case1);
            vHS = _mm256_blendv_epi8(vHS, vZero, case0);
            vHL = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEL, vFL, case2),
                    vHpL, case1);
            vHL = _mm256_blendv_epi8(vHL, vZero, case0);
            _mm256_store_si256(pvH+i, vH);
            _mm256_store_si256(pvHM+i, vHM);
            _mm256_store_si256(pvHS+i, vHS);
            _mm256_store_si256(pvHL+i, vHL);
            vSaturationCheckMin = _mm256_min_epi64_rpl(vSaturationCheckMin, vH);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vH);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
            arr_store_si256(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si256(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si256(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
#endif
            {
                __m256i cond_max = _mm256_cmpgt_epi64(vH, vMaxH);
                vMaxH = _mm256_max_epi64_rpl(vH, vMaxH);
                vMaxM = _mm256_blendv_epi8(vMaxM, vHM, cond_max);
                vMaxS = _mm256_blendv_epi8(vMaxS, vHS, cond_max);
                vMaxL = _mm256_blendv_epi8(vMaxL, vHL, cond_max);
            }
        } 

        {
            __m256i vCompare = _mm256_cmpgt_epi64(vMaxH, vMaxHUnit);
            if (_mm256_movemask_epi8(vCompare)) {
                score = _mm256_hmax_epi64_rpl(vMaxH);
                vMaxHUnit = _mm256_set1_epi64x_rpl(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvH, sizeof(__m256i)*segLen);
                (void)memcpy(pvHMMax, pvHM, sizeof(__m256i)*segLen);
                (void)memcpy(pvHSMax, pvHS, sizeof(__m256i)*segLen);
                (void)memcpy(pvHLMax, pvHL, sizeof(__m256i)*segLen);
            }
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            int32_t k = 0;
            __m256i vH = _mm256_load_si256(pvH + offset);
            __m256i vHM = _mm256_load_si256(pvHM + offset);
            __m256i vHS = _mm256_load_si256(pvHS + offset);
            __m256i vHL = _mm256_load_si256(pvHL + offset);
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 8);
                vHM = _mm256_slli_si256_rpl(vHM, 8);
                vHS = _mm256_slli_si256_rpl(vHS, 8);
                vHL = _mm256_slli_si256_rpl(vHL, 8);
            }
            result->stats->rowcols->score_row[j] = (int64_t) _mm256_extract_epi64_rpl (vH, 3);
            result->stats->rowcols->matches_row[j] = (int64_t) _mm256_extract_epi64_rpl (vHM, 3);
            result->stats->rowcols->similar_row[j] = (int64_t) _mm256_extract_epi64_rpl (vHS, 3);
            result->stats->rowcols->length_row[j] = (int64_t) _mm256_extract_epi64_rpl (vHL, 3);
        }
#endif
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

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m256i vH = _mm256_load_si256(pvH+i);
        __m256i vHM = _mm256_load_si256(pvHM+i);
        __m256i vHS = _mm256_load_si256(pvHS+i);
        __m256i vHL = _mm256_load_si256(pvHL+i);
        arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
        arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
        arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
        arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
    }
#endif

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmplt_epi64_rpl(vSaturationCheckMin, vNegLimit),
            _mm256_cmpgt_epi64(vSaturationCheckMax, vPosLimit)))) {
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
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_4;
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


