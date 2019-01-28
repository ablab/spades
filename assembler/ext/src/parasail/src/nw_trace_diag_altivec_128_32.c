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

#define NEG_INF (INT32_MIN/(int32_t)(2))


static inline void arr_store_si128(
        int8_t *array,
        vec128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int8_t)_mm_extract_epi32(vWH, 3);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int8_t)_mm_extract_epi32(vWH, 2);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int8_t)_mm_extract_epi32(vWH, 1);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int8_t)_mm_extract_epi32(vWH, 0);
    }
}

#define FNAME parasail_nw_trace_diag_altivec_128_32

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
    int32_t * const restrict s1 = parasail_memalign_int32_t(16, s1Len+PAD);
    int32_t * const restrict s2B= parasail_memalign_int32_t(16, s2Len+PAD2);
    int32_t * const restrict _H_pr = parasail_memalign_int32_t(16, s2Len+PAD2);
    int32_t * const restrict _F_pr = parasail_memalign_int32_t(16, s2Len+PAD2);
    int32_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int32_t * const restrict H_pr = _H_pr+PAD;
    int32_t * const restrict F_pr = _F_pr+PAD;
    parasail_result_t *result = parasail_result_new_trace(s1Len, s2Len, 16, sizeof(int8_t));
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = s1Len-1;
    int32_t end_ref = s2Len-1;
    int32_t score = NEG_INF;
    vec128i vNegInf = _mm_set1_epi32(NEG_INF);
    vec128i vOpen = _mm_set1_epi32(open);
    vec128i vGap  = _mm_set1_epi32(gap);
    vec128i vOne = _mm_set1_epi32(1);
    vec128i vN = _mm_set1_epi32(N);
    vec128i vGapN = _mm_set1_epi32(gap*N);
    vec128i vNegOne = _mm_set1_epi32(-1);
    vec128i vI = _mm_set_epi32(0,1,2,3);
    vec128i vJreset = _mm_set_epi32(0,-1,-2,-3);
    vec128i vMax = vNegInf;
    vec128i vILimit = _mm_set1_epi32(s1Len);
    vec128i vILimit1 = _mm_sub_epi32(vILimit, vOne);
    vec128i vJLimit = _mm_set1_epi32(s2Len);
    vec128i vJLimit1 = _mm_sub_epi32(vJLimit, vOne);
    vec128i vIBoundary = _mm_set_epi32(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap
            );
    vec128i vTDiag = _mm_set1_epi32(PARASAIL_DIAG);
    vec128i vTIns = _mm_set1_epi32(PARASAIL_INS);
    vec128i vTDel = _mm_set1_epi32(PARASAIL_DEL);
    vec128i vTDiagE = _mm_set1_epi32(PARASAIL_DIAG_E);
    vec128i vTInsE = _mm_set1_epi32(PARASAIL_INS_E);
    vec128i vTDiagF = _mm_set1_epi32(PARASAIL_DIAG_F);
    vec128i vTDelF = _mm_set1_epi32(PARASAIL_DEL_F);
    

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
        H_pr[j] = -open - j*gap;
        F_pr[j] = NEG_INF;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        H_pr[j] = NEG_INF;
        F_pr[j] = NEG_INF;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = NEG_INF;
        F_pr[j] = NEG_INF;
    }
    H_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        vec128i vNH = vNegInf;
        vec128i vWH = vNegInf;
        vec128i vE = vNegInf;
        vec128i vE_opn = vNegInf;
        vec128i vE_ext = vNegInf;
        vec128i vF = vNegInf;
        vec128i vF_opn = vNegInf;
        vec128i vF_ext = vNegInf;
        vec128i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        vNH = _mm_srli_si128(vNH, 4);
        vNH = _mm_insert_epi32(vNH, H_pr[-1], 3);
        vWH = _mm_srli_si128(vWH, 4);
        vWH = _mm_insert_epi32(vWH, -open - i*gap, 3);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            vec128i vMat;
            vec128i vNWH = vNH;
            vNH = _mm_srli_si128(vWH, 4);
            vNH = _mm_insert_epi32(vNH, H_pr[j], 3);
            vF = _mm_srli_si128(vF, 4);
            vF = _mm_insert_epi32(vF, F_pr[j], 3);
            vF_opn = _mm_sub_epi32(vNH, vOpen);
            vF_ext = _mm_sub_epi32(vF, vGap);
            vF = _mm_max_epi32(vF_opn, vF_ext);
            vE_opn = _mm_sub_epi32(vWH, vOpen);
            vE_ext = _mm_sub_epi32(vE, vGap);
            vE = _mm_max_epi32(vE_opn, vE_ext);
            vMat = _mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWH = _mm_add_epi32(vNWH, vMat);
            vWH = _mm_max_epi32(vNWH, vE);
            vWH = _mm_max_epi32(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                vec128i cond = _mm_cmpeq_epi32(vJ,vNegOne);
                vWH = _mm_blendv_epi8(vWH, vIBoundary, cond);
                vF = _mm_blendv_epi8(vF, vNegInf, cond);
                vE = _mm_blendv_epi8(vE, vNegInf, cond);
            }
            
            /* trace table */
            {
                vec128i case1 = _mm_cmpeq_epi32(vWH, vNWH);
                vec128i case2 = _mm_cmpeq_epi32(vWH, vF);
                vec128i vT = _mm_blendv_epi8(
                        _mm_blendv_epi8(vTIns, vTDel, case2),
                        vTDiag,
                        case1);
                vec128i condE = _mm_cmpgt_epi32(vE_opn, vE_ext);
                vec128i condF = _mm_cmpgt_epi32(vF_opn, vF_ext);
                vec128i vET = _mm_blendv_epi8(vTInsE, vTDiagE, condE);
                vec128i vFT = _mm_blendv_epi8(vTDelF, vTDiagF, condF);
                vT = _mm_or_si128(vT, vET);
                vT = _mm_or_si128(vT, vFT);
                arr_store_si128(result->trace->trace_table, vT, i, s1Len, j, s2Len);
            }
            H_pr[j-3] = (int32_t)_mm_extract_epi32(vWH,0);
            F_pr[j-3] = (int32_t)_mm_extract_epi32(vF,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                vec128i cond_valid_I = _mm_cmpeq_epi32(vI, vILimit1);
                vec128i cond_valid_J = _mm_cmpeq_epi32(vJ, vJLimit1);
                vec128i cond_all = _mm_and_si128(cond_valid_I, cond_valid_J);
                vMax = _mm_blendv_epi8(vMax, vWH, cond_all);
            }
            vJ = _mm_add_epi32(vJ, vOne);
        }
        vI = _mm_add_epi32(vI, vN);
        vIBoundary = _mm_sub_epi32(vIBoundary, vGapN);
    }

    /* max in vMax */
    for (i=0; i<N; ++i) {
        int32_t value;
        value = (int32_t) _mm_extract_epi32(vMax, 3);
        if (value > score) {
            score = value;
        }
        vMax = _mm_slli_si128(vMax, 4);
    }

    

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_32 | PARASAIL_FLAG_LANES_4;

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}


