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
#include "parasail/parasail/internal_neon.h"

#define NEG_INF (INT64_MIN/(int64_t)(2))


#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        simde__m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int64_t)simde_mm_extract_epi64(vWH, 1);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int64_t)simde_mm_extract_epi64(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        simde__m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int64_t)simde_mm_extract_epi64(vWH, 1);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int64_t)simde_mm_extract_epi64(vWH, 1);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int64_t)simde_mm_extract_epi64(vWH, 0);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int64_t)simde_mm_extract_epi64(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_table_diag_neon_128_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_rowcol_diag_neon_128_64
#else
#define FNAME parasail_nw_diag_neon_128_64
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 2; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int64_t * const restrict s1 = parasail_memalign_int64_t(16, s1Len+PAD);
    int64_t * const restrict s2B= parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _H_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _F_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int64_t * const restrict H_pr = _H_pr+PAD;
    int64_t * const restrict F_pr = _F_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = s1Len-1;
    int32_t end_ref = s2Len-1;
    int64_t score = NEG_INF;
    simde__m128i vNegInf = simde_mm_set1_epi64x(NEG_INF);
    simde__m128i vOpen = simde_mm_set1_epi64x(open);
    simde__m128i vGap  = simde_mm_set1_epi64x(gap);
    simde__m128i vOne = simde_mm_set1_epi64x(1);
    simde__m128i vN = simde_mm_set1_epi64x(N);
    simde__m128i vGapN = simde_mm_set1_epi64x(gap*N);
    simde__m128i vNegOne = simde_mm_set1_epi64x(-1);
    simde__m128i vI = simde_mm_set_epi64x(0,1);
    simde__m128i vJreset = simde_mm_set_epi64x(0,-1);
    simde__m128i vMax = vNegInf;
    simde__m128i vILimit = simde_mm_set1_epi64x(s1Len);
    simde__m128i vILimit1 = simde_mm_sub_epi64(vILimit, vOne);
    simde__m128i vJLimit = simde_mm_set1_epi64x(s2Len);
    simde__m128i vJLimit1 = simde_mm_sub_epi64(vJLimit, vOne);
    simde__m128i vIBoundary = simde_mm_set_epi64x(
            -open-0*gap,
            -open-1*gap
            );
    

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
        simde__m128i vNH = vNegInf;
        simde__m128i vWH = vNegInf;
        simde__m128i vE = vNegInf;
        simde__m128i vF = vNegInf;
        simde__m128i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        vNH = simde_mm_srli_si128(vNH, 8);
        vNH = simde_mm_insert_epi64(vNH, H_pr[-1], 1);
        vWH = simde_mm_srli_si128(vWH, 8);
        vWH = simde_mm_insert_epi64(vWH, -open - i*gap, 1);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            simde__m128i vMat;
            simde__m128i vNWH = vNH;
            vNH = simde_mm_srli_si128(vWH, 8);
            vNH = simde_mm_insert_epi64(vNH, H_pr[j], 1);
            vF = simde_mm_srli_si128(vF, 8);
            vF = simde_mm_insert_epi64(vF, F_pr[j], 1);
            vF = simde_mm_max_epi64(
                    simde_mm_sub_epi64(vNH, vOpen),
                    simde_mm_sub_epi64(vF, vGap));
            vE = simde_mm_max_epi64(
                    simde_mm_sub_epi64(vWH, vOpen),
                    simde_mm_sub_epi64(vE, vGap));
            vMat = simde_mm_set_epi64x(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]]
                    );
            vNWH = simde_mm_add_epi64(vNWH, vMat);
            vWH = simde_mm_max_epi64(vNWH, vE);
            vWH = simde_mm_max_epi64(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                simde__m128i cond = simde_mm_cmpeq_epi64(vJ,vNegOne);
                vWH = simde_mm_blendv_epi8(vWH, vIBoundary, cond);
                vF = simde_mm_blendv_epi8(vF, vNegInf, cond);
                vE = simde_mm_blendv_epi8(vE, vNegInf, cond);
            }
            
#ifdef PARASAIL_TABLE
            arr_store_si128(result->tables->score_table, vWH, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->rowcols->score_row, result->rowcols->score_col, vWH, i, s1Len, j, s2Len);
#endif
            H_pr[j-1] = (int64_t)simde_mm_extract_epi64(vWH,0);
            F_pr[j-1] = (int64_t)simde_mm_extract_epi64(vF,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                simde__m128i cond_valid_I = simde_mm_cmpeq_epi64(vI, vILimit1);
                simde__m128i cond_valid_J = simde_mm_cmpeq_epi64(vJ, vJLimit1);
                simde__m128i cond_all = simde_mm_and_si128(cond_valid_I, cond_valid_J);
                vMax = simde_mm_blendv_epi8(vMax, vWH, cond_all);
            }
            vJ = simde_mm_add_epi64(vJ, vOne);
        }
        vI = simde_mm_add_epi64(vI, vN);
        vIBoundary = simde_mm_sub_epi64(vIBoundary, vGapN);
    }

    /* max in vMax */
    for (i=0; i<N; ++i) {
        int64_t value;
        value = (int64_t) simde_mm_extract_epi64(vMax, 1);
        if (value > score) {
            score = value;
        }
        vMax = simde_mm_slli_si128(vMax, 8);
    }

    

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}


