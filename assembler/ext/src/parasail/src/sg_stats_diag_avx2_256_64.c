/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

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

#if HAVE_AVX2_MM256_SET_EPI64X
#define _mm256_set_epi64x_rpl _mm256_set_epi64x
#else
static inline __m256i _mm256_set_epi64x_rpl(int64_t e3, int64_t e2, int64_t e1, int64_t e0) {
    __m256i_64_t A;
    A.v[0] = e0;
    A.v[1] = e1;
    A.v[2] = e2;
    A.v[3] = e3;
    return A.m;
}
#endif

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

#define _mm256_srli_si256_rpl(a,imm) _mm256_or_si256(_mm256_slli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)), 16-imm), _mm256_srli_si256(a, imm))


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int64_t)_mm256_extract_epi64_rpl(vWH, 3);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int64_t)_mm256_extract_epi64_rpl(vWH, 2);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int64_t)_mm256_extract_epi64_rpl(vWH, 1);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int64_t)_mm256_extract_epi64_rpl(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        __m256i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int64_t)_mm256_extract_epi64_rpl(vWH, 3);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int64_t)_mm256_extract_epi64_rpl(vWH, 3);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int64_t)_mm256_extract_epi64_rpl(vWH, 2);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int64_t)_mm256_extract_epi64_rpl(vWH, 2);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int64_t)_mm256_extract_epi64_rpl(vWH, 1);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int64_t)_mm256_extract_epi64_rpl(vWH, 1);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int64_t)_mm256_extract_epi64_rpl(vWH, 0);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int64_t)_mm256_extract_epi64_rpl(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_stats_table_diag_avx2_256_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_diag_avx2_256_64
#else
#define FNAME parasail_sg_stats_diag_avx2_256_64
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 4; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int64_t * const restrict s1      = parasail_memalign_int64_t(32, s1Len+PAD);
    int64_t * const restrict s2B     = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _H_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _HM_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _HS_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _HL_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _F_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _FM_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _FS_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _FL_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int64_t * const restrict H_pr = _H_pr+PAD;
    int64_t * const restrict HM_pr = _HM_pr+PAD;
    int64_t * const restrict HS_pr = _HS_pr+PAD;
    int64_t * const restrict HL_pr = _HL_pr+PAD;
    int64_t * const restrict F_pr = _F_pr+PAD;
    int64_t * const restrict FM_pr = _FM_pr+PAD;
    int64_t * const restrict FS_pr = _FS_pr+PAD;
    int64_t * const restrict FL_pr = _FL_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const int64_t NEG_LIMIT = (-open < matrix->min ?
        INT64_MIN + open : INT64_MIN - matrix->min) + 1;
    const int64_t POS_LIMIT = INT64_MAX - matrix->max - 1;
    int64_t score = NEG_LIMIT;
    int64_t matches = NEG_LIMIT;
    int64_t similar = NEG_LIMIT;
    int64_t length = NEG_LIMIT;
    __m256i vNegLimit = _mm256_set1_epi64x_rpl(NEG_LIMIT);
    __m256i vPosLimit = _mm256_set1_epi64x_rpl(POS_LIMIT);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
    __m256i vNegInf = _mm256_set1_epi64x_rpl(NEG_LIMIT);
    __m256i vOpen = _mm256_set1_epi64x_rpl(open);
    __m256i vGap  = _mm256_set1_epi64x_rpl(gap);
    __m256i vZero = _mm256_set1_epi64x_rpl(0);
    __m256i vNegInf0 = _mm256_insert_epi64_rpl(vZero, NEG_LIMIT, 3);
    __m256i vOne = _mm256_set1_epi64x_rpl(1);
    __m256i vN = _mm256_set1_epi64x_rpl(N);
    __m256i vNegOne = _mm256_set1_epi64x_rpl(-1);
    __m256i vI = _mm256_set_epi64x_rpl(0,1,2,3);
    __m256i vJreset = _mm256_set_epi64x_rpl(0,-1,-2,-3);
    __m256i vMaxH = vNegInf;
    __m256i vMaxM = vNegInf;
    __m256i vMaxS = vNegInf;
    __m256i vMaxL = vNegInf;
    __m256i vEndI = vNegInf;
    __m256i vEndJ = vNegInf;
    __m256i vILimit = _mm256_set1_epi64x_rpl(s1Len);
    __m256i vILimit1 = _mm256_sub_epi64(vILimit, vOne);
    __m256i vJLimit = _mm256_set1_epi64x_rpl(s2Len);
    __m256i vJLimit1 = _mm256_sub_epi64(vJLimit, vOne);

    /* convert _s1 from char to int in range 0-23 */
    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    /* pad back of s1 with dummy values */
    for (i=s1Len; i<s1Len_PAD; ++i) {
        s1[i] = 0; /* point to first matrix row because we don't care */
    }

    /* convert _s2 from char to int in range 0-23 */
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }
    /* pad front of s2 with dummy values */
    for (j=-PAD; j<0; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }
    /* pad back of s2 with dummy values */
    for (j=s2Len; j<s2Len_PAD; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }

    /* set initial values for stored row */
    for (j=0; j<s2Len; ++j) {
        H_pr[j] = 0;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = NEG_LIMIT;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        H_pr[j] = 0;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = 0;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = 0;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = 0;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    H_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m256i case1 = vZero;
        __m256i case2 = vZero;
        __m256i vNH = vZero;
        __m256i vNM = vZero;
        __m256i vNS = vZero;
        __m256i vNL = vZero;
        __m256i vWH = vZero;
        __m256i vWM = vZero;
        __m256i vWS = vZero;
        __m256i vWL = vZero;
        __m256i vE = vNegInf0;
        __m256i vE_opn = vNegInf;
        __m256i vE_ext = vNegInf;
        __m256i vEM = vZero;
        __m256i vES = vZero;
        __m256i vEL = vZero;
        __m256i vF = vNegInf0;
        __m256i vF_opn = vNegInf;
        __m256i vF_ext = vNegInf;
        __m256i vFM = vZero;
        __m256i vFS = vZero;
        __m256i vFL = vZero;
        __m256i vJ = vJreset;
        __m256i vs1 = _mm256_set_epi64x_rpl(
                s1[i+0],
                s1[i+1],
                s1[i+2],
                s1[i+3]);
        __m256i vs2 = vNegInf;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        __m256i vIltLimit = _mm256_cmplt_epi64_rpl(vI, vILimit);
        __m256i vIeqLimit1 = _mm256_cmpeq_epi64(vI, vILimit1);
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m256i vMat;
            __m256i vNWH = vNH;
            __m256i vNWM = vNM;
            __m256i vNWS = vNS;
            __m256i vNWL = vNL;
            vNH = _mm256_srli_si256_rpl(vWH, 8);
            vNH = _mm256_insert_epi64_rpl(vNH, H_pr[j], 3);
            vNM = _mm256_srli_si256_rpl(vWM, 8);
            vNM = _mm256_insert_epi64_rpl(vNM, HM_pr[j], 3);
            vNS = _mm256_srli_si256_rpl(vWS, 8);
            vNS = _mm256_insert_epi64_rpl(vNS, HS_pr[j], 3);
            vNL = _mm256_srli_si256_rpl(vWL, 8);
            vNL = _mm256_insert_epi64_rpl(vNL, HL_pr[j], 3);
            vF = _mm256_srli_si256_rpl(vF, 8);
            vF = _mm256_insert_epi64_rpl(vF, F_pr[j], 3);
            vFM = _mm256_srli_si256_rpl(vFM, 8);
            vFM = _mm256_insert_epi64_rpl(vFM, FM_pr[j], 3);
            vFS = _mm256_srli_si256_rpl(vFS, 8);
            vFS = _mm256_insert_epi64_rpl(vFS, FS_pr[j], 3);
            vFL = _mm256_srli_si256_rpl(vFL, 8);
            vFL = _mm256_insert_epi64_rpl(vFL, FL_pr[j], 3);
            vF_opn = _mm256_sub_epi64(vNH, vOpen);
            vF_ext = _mm256_sub_epi64(vF, vGap);
            vF = _mm256_max_epi64_rpl(vF_opn, vF_ext);
            case1 = _mm256_cmpgt_epi64(vF_opn, vF_ext);
            vFM = _mm256_blendv_epi8(vFM, vNM, case1);
            vFS = _mm256_blendv_epi8(vFS, vNS, case1);
            vFL = _mm256_blendv_epi8(vFL, vNL, case1);
            vFL = _mm256_add_epi64(vFL, vOne);
            vE_opn = _mm256_sub_epi64(vWH, vOpen);
            vE_ext = _mm256_sub_epi64(vE, vGap);
            vE = _mm256_max_epi64_rpl(vE_opn, vE_ext);
            case1 = _mm256_cmpgt_epi64(vE_opn, vE_ext);
            vEM = _mm256_blendv_epi8(vEM, vWM, case1);
            vES = _mm256_blendv_epi8(vES, vWS, case1);
            vEL = _mm256_blendv_epi8(vEL, vWL, case1);
            vEL = _mm256_add_epi64(vEL, vOne);
            vs2 = _mm256_srli_si256_rpl(vs2, 8);
            vs2 = _mm256_insert_epi64_rpl(vs2, s2[j], 3);
            vMat = _mm256_set_epi64x_rpl(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWH = _mm256_add_epi64(vNWH, vMat);
            vWH = _mm256_max_epi64_rpl(vNWH, vE);
            vWH = _mm256_max_epi64_rpl(vWH, vF);
            case1 = _mm256_cmpeq_epi64(vWH, vNWH);
            case2 = _mm256_cmpeq_epi64(vWH, vF);
            vWM = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEM, vFM, case2),
                    _mm256_add_epi64(vNWM,
                        _mm256_and_si256(
                            _mm256_cmpeq_epi64(vs1,vs2),
                            vOne)),
                    case1);
            vWS = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vES, vFS, case2),
                    _mm256_add_epi64(vNWS,
                        _mm256_and_si256(
                            _mm256_cmpgt_epi64(vMat,vZero),
                            vOne)),
                    case1);
            vWL = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEL, vFL, case2),
                    _mm256_add_epi64(vNWL, vOne), case1);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m256i cond = _mm256_cmpeq_epi64(vJ,vNegOne);
                vWH = _mm256_andnot_si256(cond, vWH);
                vWM = _mm256_andnot_si256(cond, vWM);
                vWS = _mm256_andnot_si256(cond, vWS);
                vWL = _mm256_andnot_si256(cond, vWL);
                vE = _mm256_blendv_epi8(vE, vNegInf, cond);
                vEM = _mm256_andnot_si256(cond, vEM);
                vES = _mm256_andnot_si256(cond, vES);
                vEL = _mm256_andnot_si256(cond, vEL);
            }
            vSaturationCheckMin = _mm256_min_epi64_rpl(vSaturationCheckMin, vWH);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vWH);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vWM);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vWS);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vWL);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vWL);
            vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vJ);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->stats->tables->score_table, vWH, i, s1Len, j, s2Len);
            arr_store_si256(result->stats->tables->matches_table, vWM, i, s1Len, j, s2Len);
            arr_store_si256(result->stats->tables->similar_table, vWS, i, s1Len, j, s2Len);
            arr_store_si256(result->stats->tables->length_table, vWL, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->stats->rowcols->score_row,   result->stats->rowcols->score_col, vWH, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->matches_row, result->stats->rowcols->matches_col, vWM, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->similar_row, result->stats->rowcols->similar_col, vWS, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->length_row,  result->stats->rowcols->length_col, vWL, i, s1Len, j, s2Len);
#endif
            H_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vWH,0);
            HM_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vWM,0);
            HS_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vWS,0);
            HL_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vWL,0);
            F_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vF,0);
            FM_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vFM,0);
            FS_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vFS,0);
            FL_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vFL,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                __m256i vJeqLimit1 = _mm256_cmpeq_epi64(vJ, vJLimit1);
                __m256i vJgtNegOne = _mm256_cmpgt_epi64(vJ, vNegOne);
                __m256i vJltLimit = _mm256_cmplt_epi64_rpl(vJ, vJLimit);
                __m256i cond_j = _mm256_and_si256(vIltLimit, vJeqLimit1);
                __m256i cond_i = _mm256_and_si256(vIeqLimit1,
                        _mm256_and_si256(vJgtNegOne, vJltLimit));
                __m256i cond_valid_IJ = _mm256_or_si256(cond_i, cond_j);
                __m256i cond_eq = _mm256_cmpeq_epi64(vWH, vMaxH);
                __m256i cond_max = _mm256_cmpgt_epi64(vWH, vMaxH);
                __m256i cond_all = _mm256_and_si256(cond_max, cond_valid_IJ);
                __m256i cond_Jlt = _mm256_cmplt_epi64_rpl(vJ, vEndJ);
                vMaxH = _mm256_blendv_epi8(vMaxH, vWH, cond_all);
                vMaxM = _mm256_blendv_epi8(vMaxM, vWM, cond_all);
                vMaxS = _mm256_blendv_epi8(vMaxS, vWS, cond_all);
                vMaxL = _mm256_blendv_epi8(vMaxL, vWL, cond_all);
                vEndI = _mm256_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = _mm256_blendv_epi8(vEndJ, vJ, cond_all);
                cond_all = _mm256_and_si256(cond_Jlt, cond_eq);
                cond_all = _mm256_and_si256(cond_all, cond_valid_IJ);
                vMaxM = _mm256_blendv_epi8(vMaxM, vWM, cond_all);
                vMaxS = _mm256_blendv_epi8(vMaxS, vWS, cond_all);
                vMaxL = _mm256_blendv_epi8(vMaxL, vWL, cond_all);
                vEndI = _mm256_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = _mm256_blendv_epi8(vEndJ, vJ, cond_all);
            }
            vJ = _mm256_add_epi64(vJ, vOne);
        }
        vI = _mm256_add_epi64(vI, vN);
        vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vI);
    }

    /* alignment ending position */
    {
        int64_t *t = (int64_t*)&vMaxH;
        int64_t *m = (int64_t*)&vMaxM;
        int64_t *s = (int64_t*)&vMaxS;
        int64_t *l = (int64_t*)&vMaxL;
        int64_t *i = (int64_t*)&vEndI;
        int64_t *j = (int64_t*)&vEndJ;
        int32_t k;
        for (k=0; k<N; ++k, ++t, ++m, ++s, ++l, ++i, ++j) {
            if (*t > score) {
                score = *t;
                matches = *m;
                similar = *s;
                length = *l;
                end_query = *i;
                end_ref = *j;
            }
            else if (*t == score) {
                if (*j < end_ref) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *i;
                    end_ref = *j;
                }
                else if (*j == end_ref && *i < end_query) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *i;
                    end_ref = *j;
                }
            }
        }
    }

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
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_4;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(_FL_pr);
    parasail_free(_FS_pr);
    parasail_free(_FM_pr);
    parasail_free(_F_pr);
    parasail_free(_HL_pr);
    parasail_free(_HS_pr);
    parasail_free(_HM_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}


