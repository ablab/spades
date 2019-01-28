/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>



#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_altivec.h"



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        vec128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int8_t)_mm_extract_epi8(vWH, 15);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int8_t)_mm_extract_epi8(vWH, 14);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int8_t)_mm_extract_epi8(vWH, 13);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int8_t)_mm_extract_epi8(vWH, 12);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[1LL*(i+4)*s2Len + (j-4)] = (int8_t)_mm_extract_epi8(vWH, 11);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[1LL*(i+5)*s2Len + (j-5)] = (int8_t)_mm_extract_epi8(vWH, 10);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[1LL*(i+6)*s2Len + (j-6)] = (int8_t)_mm_extract_epi8(vWH, 9);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[1LL*(i+7)*s2Len + (j-7)] = (int8_t)_mm_extract_epi8(vWH, 8);
    }
    if (0 <= i+8 && i+8 < s1Len && 0 <= j-8 && j-8 < s2Len) {
        array[1LL*(i+8)*s2Len + (j-8)] = (int8_t)_mm_extract_epi8(vWH, 7);
    }
    if (0 <= i+9 && i+9 < s1Len && 0 <= j-9 && j-9 < s2Len) {
        array[1LL*(i+9)*s2Len + (j-9)] = (int8_t)_mm_extract_epi8(vWH, 6);
    }
    if (0 <= i+10 && i+10 < s1Len && 0 <= j-10 && j-10 < s2Len) {
        array[1LL*(i+10)*s2Len + (j-10)] = (int8_t)_mm_extract_epi8(vWH, 5);
    }
    if (0 <= i+11 && i+11 < s1Len && 0 <= j-11 && j-11 < s2Len) {
        array[1LL*(i+11)*s2Len + (j-11)] = (int8_t)_mm_extract_epi8(vWH, 4);
    }
    if (0 <= i+12 && i+12 < s1Len && 0 <= j-12 && j-12 < s2Len) {
        array[1LL*(i+12)*s2Len + (j-12)] = (int8_t)_mm_extract_epi8(vWH, 3);
    }
    if (0 <= i+13 && i+13 < s1Len && 0 <= j-13 && j-13 < s2Len) {
        array[1LL*(i+13)*s2Len + (j-13)] = (int8_t)_mm_extract_epi8(vWH, 2);
    }
    if (0 <= i+14 && i+14 < s1Len && 0 <= j-14 && j-14 < s2Len) {
        array[1LL*(i+14)*s2Len + (j-14)] = (int8_t)_mm_extract_epi8(vWH, 1);
    }
    if (0 <= i+15 && i+15 < s1Len && 0 <= j-15 && j-15 < s2Len) {
        array[1LL*(i+15)*s2Len + (j-15)] = (int8_t)_mm_extract_epi8(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        vec128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int8_t)_mm_extract_epi8(vWH, 15);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int8_t)_mm_extract_epi8(vWH, 15);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int8_t)_mm_extract_epi8(vWH, 14);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int8_t)_mm_extract_epi8(vWH, 14);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int8_t)_mm_extract_epi8(vWH, 13);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int8_t)_mm_extract_epi8(vWH, 13);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int8_t)_mm_extract_epi8(vWH, 12);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int8_t)_mm_extract_epi8(vWH, 12);
    }
    if (i+4 == s1Len-1 && 0 <= j-4 && j-4 < s2Len) {
        row[j-4] = (int8_t)_mm_extract_epi8(vWH, 11);
    }
    if (j-4 == s2Len-1 && 0 <= i+4 && i+4 < s1Len) {
        col[(i+4)] = (int8_t)_mm_extract_epi8(vWH, 11);
    }
    if (i+5 == s1Len-1 && 0 <= j-5 && j-5 < s2Len) {
        row[j-5] = (int8_t)_mm_extract_epi8(vWH, 10);
    }
    if (j-5 == s2Len-1 && 0 <= i+5 && i+5 < s1Len) {
        col[(i+5)] = (int8_t)_mm_extract_epi8(vWH, 10);
    }
    if (i+6 == s1Len-1 && 0 <= j-6 && j-6 < s2Len) {
        row[j-6] = (int8_t)_mm_extract_epi8(vWH, 9);
    }
    if (j-6 == s2Len-1 && 0 <= i+6 && i+6 < s1Len) {
        col[(i+6)] = (int8_t)_mm_extract_epi8(vWH, 9);
    }
    if (i+7 == s1Len-1 && 0 <= j-7 && j-7 < s2Len) {
        row[j-7] = (int8_t)_mm_extract_epi8(vWH, 8);
    }
    if (j-7 == s2Len-1 && 0 <= i+7 && i+7 < s1Len) {
        col[(i+7)] = (int8_t)_mm_extract_epi8(vWH, 8);
    }
    if (i+8 == s1Len-1 && 0 <= j-8 && j-8 < s2Len) {
        row[j-8] = (int8_t)_mm_extract_epi8(vWH, 7);
    }
    if (j-8 == s2Len-1 && 0 <= i+8 && i+8 < s1Len) {
        col[(i+8)] = (int8_t)_mm_extract_epi8(vWH, 7);
    }
    if (i+9 == s1Len-1 && 0 <= j-9 && j-9 < s2Len) {
        row[j-9] = (int8_t)_mm_extract_epi8(vWH, 6);
    }
    if (j-9 == s2Len-1 && 0 <= i+9 && i+9 < s1Len) {
        col[(i+9)] = (int8_t)_mm_extract_epi8(vWH, 6);
    }
    if (i+10 == s1Len-1 && 0 <= j-10 && j-10 < s2Len) {
        row[j-10] = (int8_t)_mm_extract_epi8(vWH, 5);
    }
    if (j-10 == s2Len-1 && 0 <= i+10 && i+10 < s1Len) {
        col[(i+10)] = (int8_t)_mm_extract_epi8(vWH, 5);
    }
    if (i+11 == s1Len-1 && 0 <= j-11 && j-11 < s2Len) {
        row[j-11] = (int8_t)_mm_extract_epi8(vWH, 4);
    }
    if (j-11 == s2Len-1 && 0 <= i+11 && i+11 < s1Len) {
        col[(i+11)] = (int8_t)_mm_extract_epi8(vWH, 4);
    }
    if (i+12 == s1Len-1 && 0 <= j-12 && j-12 < s2Len) {
        row[j-12] = (int8_t)_mm_extract_epi8(vWH, 3);
    }
    if (j-12 == s2Len-1 && 0 <= i+12 && i+12 < s1Len) {
        col[(i+12)] = (int8_t)_mm_extract_epi8(vWH, 3);
    }
    if (i+13 == s1Len-1 && 0 <= j-13 && j-13 < s2Len) {
        row[j-13] = (int8_t)_mm_extract_epi8(vWH, 2);
    }
    if (j-13 == s2Len-1 && 0 <= i+13 && i+13 < s1Len) {
        col[(i+13)] = (int8_t)_mm_extract_epi8(vWH, 2);
    }
    if (i+14 == s1Len-1 && 0 <= j-14 && j-14 < s2Len) {
        row[j-14] = (int8_t)_mm_extract_epi8(vWH, 1);
    }
    if (j-14 == s2Len-1 && 0 <= i+14 && i+14 < s1Len) {
        col[(i+14)] = (int8_t)_mm_extract_epi8(vWH, 1);
    }
    if (i+15 == s1Len-1 && 0 <= j-15 && j-15 < s2Len) {
        row[j-15] = (int8_t)_mm_extract_epi8(vWH, 0);
    }
    if (j-15 == s2Len-1 && 0 <= i+15 && i+15 < s1Len) {
        col[(i+15)] = (int8_t)_mm_extract_epi8(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_diag_altivec_128_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_diag_altivec_128_8
#else
#define FNAME parasail_sw_stats_diag_altivec_128_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 16; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int8_t * const restrict s1      = parasail_memalign_int8_t(16, s1Len+PAD);
    int8_t * const restrict s2B     = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _H_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _HM_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _HS_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _HL_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _F_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _FM_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _FS_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _FL_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int8_t * const restrict H_pr = _H_pr+PAD;
    int8_t * const restrict HM_pr = _HM_pr+PAD;
    int8_t * const restrict HS_pr = _HS_pr+PAD;
    int8_t * const restrict HL_pr = _HL_pr+PAD;
    int8_t * const restrict F_pr = _F_pr+PAD;
    int8_t * const restrict FM_pr = _FM_pr+PAD;
    int8_t * const restrict FS_pr = _FS_pr+PAD;
    int8_t * const restrict FL_pr = _FL_pr+PAD;
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
    const int8_t NEG_LIMIT = (-open < matrix->min ?
        INT8_MIN + open : INT8_MIN - matrix->min) + 1;
    const int8_t POS_LIMIT = INT8_MAX - matrix->max - 1;
    int8_t score = NEG_LIMIT;
    int8_t matches = NEG_LIMIT;
    int8_t similar = NEG_LIMIT;
    int8_t length = NEG_LIMIT;
    vec128i vNegLimit = _mm_set1_epi8(NEG_LIMIT);
    vec128i vPosLimit = _mm_set1_epi8(POS_LIMIT);
    vec128i vSaturationCheckMin = vPosLimit;
    vec128i vSaturationCheckMax = vNegLimit;
    vec128i vNegInf = _mm_set1_epi8(NEG_LIMIT);
    vec128i vNegInf0 = _mm_srli_si128(vNegInf, 1); /* shift in a 0 */
    vec128i vOpen = _mm_set1_epi8(open);
    vec128i vGap  = _mm_set1_epi8(gap);
    vec128i vZero = _mm_set1_epi8(0);
    vec128i vOne = _mm_set1_epi8(1);
    vec128i vOne16 = _mm_set1_epi16(1);
    vec128i vNegOne16 = _mm_set1_epi16(-1);
    vec128i vN16 = _mm_set1_epi16(N);
    vec128i vILo16 = _mm_set_epi16(8,9,10,11,12,13,14,15);
    vec128i vIHi16 = _mm_set_epi16(0,1,2,3,4,5,6,7);
    vec128i vJresetLo16 = _mm_set_epi16(-8,-9,-10,-11,-12,-13,-14,-15);
    vec128i vJresetHi16 = _mm_set_epi16(0,-1,-2,-3,-4,-5,-6,-7);
    vec128i vMaxH = vNegInf;
    vec128i vMaxM = vNegInf;
    vec128i vMaxS = vNegInf;
    vec128i vMaxL = vNegInf;
    vec128i vEndILo = vNegInf;
    vec128i vEndIHi = vNegInf;
    vec128i vEndJLo = vNegInf;
    vec128i vEndJHi = vNegInf;
    vec128i vILimit16 = _mm_set1_epi16(s1Len);
    vec128i vJLimit16 = _mm_set1_epi16(s2Len);

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
        H_pr[j] = NEG_LIMIT;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = NEG_LIMIT;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = NEG_LIMIT;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = NEG_LIMIT;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    H_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        vec128i case1 = vZero;
        vec128i case2 = vZero;
        vec128i case0 = vZero;
        vec128i vNH = vNegInf0;
        vec128i vNM = vZero;
        vec128i vNS = vZero;
        vec128i vNL = vZero;
        vec128i vWH = vNegInf0;
        vec128i vWM = vZero;
        vec128i vWS = vZero;
        vec128i vWL = vZero;
        vec128i vE = vNegInf;
        vec128i vE_opn = vNegInf;
        vec128i vE_ext = vNegInf;
        vec128i vEM = vZero;
        vec128i vES = vZero;
        vec128i vEL = vZero;
        vec128i vF = vNegInf;
        vec128i vF_opn = vNegInf;
        vec128i vF_ext = vNegInf;
        vec128i vFM = vZero;
        vec128i vFS = vZero;
        vec128i vFL = vZero;
        vec128i vJLo16 = vJresetLo16;
        vec128i vJHi16 = vJresetHi16;
        vec128i vs1 = _mm_set_epi8(
                s1[i+0],
                s1[i+1],
                s1[i+2],
                s1[i+3],
                s1[i+4],
                s1[i+5],
                s1[i+6],
                s1[i+7],
                s1[i+8],
                s1[i+9],
                s1[i+10],
                s1[i+11],
                s1[i+12],
                s1[i+13],
                s1[i+14],
                s1[i+15]);
        vec128i vs2 = vNegInf;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        const int * const restrict matrow4 = &matrix->matrix[matrix->size*s1[i+4]];
        const int * const restrict matrow5 = &matrix->matrix[matrix->size*s1[i+5]];
        const int * const restrict matrow6 = &matrix->matrix[matrix->size*s1[i+6]];
        const int * const restrict matrow7 = &matrix->matrix[matrix->size*s1[i+7]];
        const int * const restrict matrow8 = &matrix->matrix[matrix->size*s1[i+8]];
        const int * const restrict matrow9 = &matrix->matrix[matrix->size*s1[i+9]];
        const int * const restrict matrow10 = &matrix->matrix[matrix->size*s1[i+10]];
        const int * const restrict matrow11 = &matrix->matrix[matrix->size*s1[i+11]];
        const int * const restrict matrow12 = &matrix->matrix[matrix->size*s1[i+12]];
        const int * const restrict matrow13 = &matrix->matrix[matrix->size*s1[i+13]];
        const int * const restrict matrow14 = &matrix->matrix[matrix->size*s1[i+14]];
        const int * const restrict matrow15 = &matrix->matrix[matrix->size*s1[i+15]];
        vec128i vIltLimit = _mm_packs_epi16(
                    _mm_cmplt_epi16(vILo16, vILimit16),
                    _mm_cmplt_epi16(vIHi16, vILimit16));
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            vec128i vMat;
            vec128i vNWH = vNH;
            vec128i vNWM = vNM;
            vec128i vNWS = vNS;
            vec128i vNWL = vNL;
            vNH = _mm_srli_si128(vWH, 1);
            vNH = _mm_insert_epi8(vNH, H_pr[j], 15);
            vNM = _mm_srli_si128(vWM, 1);
            vNM = _mm_insert_epi8(vNM, HM_pr[j], 15);
            vNS = _mm_srli_si128(vWS, 1);
            vNS = _mm_insert_epi8(vNS, HS_pr[j], 15);
            vNL = _mm_srli_si128(vWL, 1);
            vNL = _mm_insert_epi8(vNL, HL_pr[j], 15);
            vF = _mm_srli_si128(vF, 1);
            vF = _mm_insert_epi8(vF, F_pr[j], 15);
            vFM = _mm_srli_si128(vFM, 1);
            vFM = _mm_insert_epi8(vFM, FM_pr[j], 15);
            vFS = _mm_srli_si128(vFS, 1);
            vFS = _mm_insert_epi8(vFS, FS_pr[j], 15);
            vFL = _mm_srli_si128(vFL, 1);
            vFL = _mm_insert_epi8(vFL, FL_pr[j], 15);
            vF_opn = _mm_subs_epi8(vNH, vOpen);
            vF_ext = _mm_subs_epi8(vF, vGap);
            vF = _mm_max_epi8(vF_opn, vF_ext);
            case1 = _mm_cmpgt_epi8(vF_opn, vF_ext);
            vFM = _mm_blendv_epi8(vFM, vNM, case1);
            vFS = _mm_blendv_epi8(vFS, vNS, case1);
            vFL = _mm_blendv_epi8(vFL, vNL, case1);
            vFL = _mm_adds_epi8(vFL, vOne);
            vE_opn = _mm_subs_epi8(vWH, vOpen);
            vE_ext = _mm_subs_epi8(vE, vGap);
            vE = _mm_max_epi8(vE_opn, vE_ext);
            case1 = _mm_cmpgt_epi8(vE_opn, vE_ext);
            vEM = _mm_blendv_epi8(vEM, vWM, case1);
            vES = _mm_blendv_epi8(vES, vWS, case1);
            vEL = _mm_blendv_epi8(vEL, vWL, case1);
            vEL = _mm_adds_epi8(vEL, vOne);
            vs2 = _mm_srli_si128(vs2, 1);
            vs2 = _mm_insert_epi8(vs2, s2[j], 15);
            vMat = _mm_set_epi8(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]],
                    matrow8[s2[j-8]],
                    matrow9[s2[j-9]],
                    matrow10[s2[j-10]],
                    matrow11[s2[j-11]],
                    matrow12[s2[j-12]],
                    matrow13[s2[j-13]],
                    matrow14[s2[j-14]],
                    matrow15[s2[j-15]]
                    );
            vNWH = _mm_adds_epi8(vNWH, vMat);
            vWH = _mm_max_epi8(vNWH, vE);
            vWH = _mm_max_epi8(vWH, vF);
            vWH = _mm_max_epi8(vWH, vZero);
            case1 = _mm_cmpeq_epi8(vWH, vNWH);
            case2 = _mm_cmpeq_epi8(vWH, vF);
            case0 = _mm_cmpeq_epi8(vWH, vZero);
            vWM = _mm_blendv_epi8(
                    _mm_blendv_epi8(vEM, vFM, case2),
                    _mm_adds_epi8(vNWM,
                        _mm_and_si128(
                            _mm_cmpeq_epi8(vs1,vs2),
                            vOne)),
                    case1);
            vWM = _mm_blendv_epi8(vWM, vZero, case0);
            vWS = _mm_blendv_epi8(
                    _mm_blendv_epi8(vES, vFS, case2),
                    _mm_adds_epi8(vNWS,
                        _mm_and_si128(
                            _mm_cmpgt_epi8(vMat,vZero),
                            vOne)),
                    case1);
            vWS = _mm_blendv_epi8(vWS, vZero, case0);
            vWL = _mm_blendv_epi8(
                    _mm_blendv_epi8(vEL, vFL, case2),
                    _mm_adds_epi8(vNWL, vOne), case1);
            vWL = _mm_blendv_epi8(vWL, vZero, case0);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                vec128i cond = _mm_packs_epi16(
                        _mm_cmpeq_epi16(vJLo16,vNegOne16),
                        _mm_cmpeq_epi16(vJHi16,vNegOne16));
                vWH = _mm_andnot_si128(cond, vWH);
                vWM = _mm_andnot_si128(cond, vWM);
                vWS = _mm_andnot_si128(cond, vWS);
                vWL = _mm_andnot_si128(cond, vWL);
                vE = _mm_blendv_epi8(vE, vNegInf, cond);
                vEM = _mm_andnot_si128(cond, vEM);
                vES = _mm_andnot_si128(cond, vES);
                vEL = _mm_andnot_si128(cond, vEL);
            }
            vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vWH);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vWH);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vWM);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vWS);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vWL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->stats->tables->score_table, vWH, i, s1Len, j, s2Len);
            arr_store_si128(result->stats->tables->matches_table, vWM, i, s1Len, j, s2Len);
            arr_store_si128(result->stats->tables->similar_table, vWS, i, s1Len, j, s2Len);
            arr_store_si128(result->stats->tables->length_table, vWL, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->stats->rowcols->score_row,   result->stats->rowcols->score_col, vWH, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->matches_row, result->stats->rowcols->matches_col, vWM, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->similar_row, result->stats->rowcols->similar_col, vWS, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->length_row,  result->stats->rowcols->length_col, vWL, i, s1Len, j, s2Len);
#endif
            H_pr[j-15] = (int8_t)_mm_extract_epi8(vWH,0);
            HM_pr[j-15] = (int8_t)_mm_extract_epi8(vWM,0);
            HS_pr[j-15] = (int8_t)_mm_extract_epi8(vWS,0);
            HL_pr[j-15] = (int8_t)_mm_extract_epi8(vWL,0);
            F_pr[j-15] = (int8_t)_mm_extract_epi8(vF,0);
            FM_pr[j-15] = (int8_t)_mm_extract_epi8(vFM,0);
            FS_pr[j-15] = (int8_t)_mm_extract_epi8(vFS,0);
            FL_pr[j-15] = (int8_t)_mm_extract_epi8(vFL,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                vec128i cond_valid_J = _mm_and_si128(
                        _mm_packs_epi16(
                            _mm_cmpgt_epi16(vJLo16, vNegOne16),
                            _mm_cmpgt_epi16(vJHi16, vNegOne16)),
                        _mm_packs_epi16(
                            _mm_cmplt_epi16(vJLo16, vJLimit16),
                            _mm_cmplt_epi16(vJHi16, vJLimit16)));
                vec128i cond_valid_IJ = _mm_and_si128(cond_valid_J, vIltLimit);
                vec128i cond_eq = _mm_cmpeq_epi8(vWH, vMaxH);
                vec128i cond_max = _mm_cmpgt_epi8(vWH, vMaxH);
                vec128i cond_all = _mm_and_si128(cond_max, cond_valid_IJ);
                vec128i cond_Jlt = _mm_packs_epi16(
                        _mm_cmplt_epi16(vJLo16, vEndJLo),
                        _mm_cmplt_epi16(vJHi16, vEndJHi));
                vec128i cond_lo = _mm_unpacklo_epi8(cond_all, cond_all);
                vec128i cond_hi = _mm_unpackhi_epi8(cond_all, cond_all);
                vMaxH = _mm_blendv_epi8(vMaxH, vWH, cond_all);
                vMaxM = _mm_blendv_epi8(vMaxM, vWM, cond_all);
                vMaxS = _mm_blendv_epi8(vMaxS, vWS, cond_all);
                vMaxL = _mm_blendv_epi8(vMaxL, vWL, cond_all);
                vEndILo = _mm_blendv_epi8(vEndILo, vILo16, cond_lo);
                vEndIHi = _mm_blendv_epi8(vEndIHi, vIHi16, cond_hi);
                vEndJLo = _mm_blendv_epi8(vEndJLo, vJLo16, cond_lo);
                vEndJHi = _mm_blendv_epi8(vEndJHi, vJHi16, cond_hi);
                cond_all = _mm_and_si128(cond_Jlt, cond_eq);
                cond_all = _mm_and_si128(cond_all, cond_valid_IJ);
                cond_lo = _mm_unpacklo_epi8(cond_all, cond_all);
                cond_hi = _mm_unpackhi_epi8(cond_all, cond_all);
                vMaxM = _mm_blendv_epi8(vMaxM, vWM, cond_all);
                vMaxS = _mm_blendv_epi8(vMaxS, vWS, cond_all);
                vMaxL = _mm_blendv_epi8(vMaxL, vWL, cond_all);
                vEndILo = _mm_blendv_epi8(vEndILo, vILo16, cond_lo);
                vEndIHi = _mm_blendv_epi8(vEndIHi, vIHi16, cond_hi);
                vEndJLo = _mm_blendv_epi8(vEndJLo, vJLo16, cond_lo);
                vEndJHi = _mm_blendv_epi8(vEndJHi, vJHi16, cond_hi);
            }
            vJLo16 = _mm_add_epi16(vJLo16, vOne16);
            vJHi16 = _mm_add_epi16(vJHi16, vOne16);
        }
        vILo16 = _mm_add_epi16(vILo16, vN16);
        vIHi16 = _mm_add_epi16(vIHi16, vN16);
    }

    /* alignment ending position */
    {
        int8_t *t = (int8_t*)&vMaxH;
        int8_t *m = (int8_t*)&vMaxM;
        int8_t *s = (int8_t*)&vMaxS;
        int8_t *l = (int8_t*)&vMaxL;
        int16_t *ilo = (int16_t*)&vEndILo;
        int16_t *jlo = (int16_t*)&vEndJLo;
        int16_t *ihi = (int16_t*)&vEndIHi;
        int16_t *jhi = (int16_t*)&vEndJHi;
        int32_t k;
        for (k=0; k<N/2; ++k, ++t, ++m, ++s, ++l, ++ilo, ++jlo) {
            if (*t > score) {
                score = *t;
                matches = *m;
                similar = *s;
                length = *l;
                end_query = *ilo;
                end_ref = *jlo;
            }
            else if (*t == score) {
                if (*jlo < end_ref) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *ilo;
                    end_ref = *jlo;
                }
                else if (*jlo == end_ref && *ilo < end_query) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *ilo;
                    end_ref = *jlo;
                }
            }
        }
        for (k=N/2; k<N; ++k, ++t, ++m, ++s, ++l, ++ihi, ++jhi) {
            if (*t > score) {
                score = *t;
                matches = *m;
                similar = *s;
                length = *l;
                end_query = *ihi;
                end_ref = *jhi;
            }
            else if (*t == score) {
                if (*jhi < end_ref) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *ihi;
                    end_ref = *jhi;
                }
                else if (*jhi == end_ref && *ihi < end_query) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *ihi;
                    end_ref = *jhi;
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
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;
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


