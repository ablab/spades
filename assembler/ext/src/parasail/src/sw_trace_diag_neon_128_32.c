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

#define NEG_INF (INT32_MIN/(int32_t)(2))


static inline void arr_store_si128(
        int8_t *array,
        simde__m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int8_t)simde_mm_extract_epi32(vWH, 3);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int8_t)simde_mm_extract_epi32(vWH, 2);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int8_t)simde_mm_extract_epi32(vWH, 1);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int8_t)simde_mm_extract_epi32(vWH, 0);
    }
}

#define FNAME parasail_sw_trace_diag_neon_128_32

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
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int32_t score = NEG_INF;
    simde__m128i vNegInf = simde_mm_set1_epi32(NEG_INF);
    simde__m128i vNegInf0 = simde_mm_srli_si128(vNegInf, 4); /* shift in a 0 */
    simde__m128i vOpen = simde_mm_set1_epi32(open);
    simde__m128i vGap  = simde_mm_set1_epi32(gap);
    simde__m128i vZero = simde_mm_set1_epi32(0);
    simde__m128i vOne = simde_mm_set1_epi32(1);
    simde__m128i vN = simde_mm_set1_epi32(N);
    simde__m128i vNegOne = simde_mm_set1_epi32(-1);
    simde__m128i vI = simde_mm_set_epi32(0,1,2,3);
    simde__m128i vJreset = simde_mm_set_epi32(0,-1,-2,-3);
    simde__m128i vMaxH = vNegInf;
    simde__m128i vEndI = vNegInf;
    simde__m128i vEndJ = vNegInf;
    simde__m128i vILimit = simde_mm_set1_epi32(s1Len);
    simde__m128i vJLimit = simde_mm_set1_epi32(s2Len);
    simde__m128i vTDiag = simde_mm_set1_epi32(PARASAIL_DIAG);
    simde__m128i vTIns = simde_mm_set1_epi32(PARASAIL_INS);
    simde__m128i vTDel = simde_mm_set1_epi32(PARASAIL_DEL);
    simde__m128i vTZero = simde_mm_set1_epi32(PARASAIL_ZERO);
    simde__m128i vTDiagE = simde_mm_set1_epi32(PARASAIL_DIAG_E);
    simde__m128i vTInsE = simde_mm_set1_epi32(PARASAIL_INS_E);
    simde__m128i vTDiagF = simde_mm_set1_epi32(PARASAIL_DIAG_F);
    simde__m128i vTDelF = simde_mm_set1_epi32(PARASAIL_DEL_F);
    

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

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        simde__m128i vNH = vNegInf0;
        simde__m128i vWH = vNegInf0;
        simde__m128i vE = vNegInf;
        simde__m128i vE_opn = vNegInf;
        simde__m128i vE_ext = vNegInf;
        simde__m128i vF = vNegInf;
        simde__m128i vF_opn = vNegInf;
        simde__m128i vF_ext = vNegInf;
        simde__m128i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        simde__m128i vIltLimit = simde_mm_cmplt_epi32(vI, vILimit);
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            simde__m128i vMat;
            simde__m128i vNWH = vNH;
            vNH = simde_mm_srli_si128(vWH, 4);
            vNH = simde_mm_insert_epi32(vNH, H_pr[j], 3);
            vF = simde_mm_srli_si128(vF, 4);
            vF = simde_mm_insert_epi32(vF, F_pr[j], 3);
            vF_opn = simde_mm_sub_epi32(vNH, vOpen);
            vF_ext = simde_mm_sub_epi32(vF, vGap);
            vF = simde_mm_max_epi32(vF_opn, vF_ext);
            vE_opn = simde_mm_sub_epi32(vWH, vOpen);
            vE_ext = simde_mm_sub_epi32(vE, vGap);
            vE = simde_mm_max_epi32(vE_opn, vE_ext);
            vMat = simde_mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWH = simde_mm_add_epi32(vNWH, vMat);
            vNWH = simde_mm_max_epi32(vNWH, vZero);
            vWH = simde_mm_max_epi32(vNWH, vE);
            vWH = simde_mm_max_epi32(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                simde__m128i cond = simde_mm_cmpeq_epi32(vJ,vNegOne);
                vWH = simde_mm_andnot_si128(cond, vWH);
                vF = simde_mm_blendv_epi8(vF, vNegInf, cond);
                vE = simde_mm_blendv_epi8(vE, vNegInf, cond);
            }
            
            /* trace table */
            {
                simde__m128i cond_zero = simde_mm_cmpeq_epi32(vWH, vZero);
                simde__m128i case1 = simde_mm_cmpeq_epi32(vWH, vNWH);
                simde__m128i case2 = simde_mm_cmpeq_epi32(vWH, vF);
                simde__m128i vT = simde_mm_blendv_epi8(
                        simde_mm_blendv_epi8(vTIns, vTDel, case2),
                        simde_mm_blendv_epi8(vTDiag, vTZero, cond_zero),
                        case1);
                simde__m128i condE = simde_mm_cmpgt_epi32(vE_opn, vE_ext);
                simde__m128i condF = simde_mm_cmpgt_epi32(vF_opn, vF_ext);
                simde__m128i vET = simde_mm_blendv_epi8(vTInsE, vTDiagE, condE);
                simde__m128i vFT = simde_mm_blendv_epi8(vTDelF, vTDiagF, condF);
                vT = simde_mm_or_si128(vT, vET);
                vT = simde_mm_or_si128(vT, vFT);
                arr_store_si128(result->trace->trace_table, vT, i, s1Len, j, s2Len);
            }
            H_pr[j-3] = (int32_t)simde_mm_extract_epi32(vWH,0);
            F_pr[j-3] = (int32_t)simde_mm_extract_epi32(vF,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                simde__m128i cond_valid_J = simde_mm_and_si128(
                        simde_mm_cmpgt_epi32(vJ, vNegOne),
                        simde_mm_cmplt_epi32(vJ, vJLimit));
                simde__m128i cond_valid_IJ = simde_mm_and_si128(cond_valid_J, vIltLimit);
                simde__m128i cond_eq = simde_mm_cmpeq_epi32(vWH, vMaxH);
                simde__m128i cond_max = simde_mm_cmpgt_epi32(vWH, vMaxH);
                simde__m128i cond_all = simde_mm_and_si128(cond_max, cond_valid_IJ);
                simde__m128i cond_Jlt = simde_mm_cmplt_epi32(vJ, vEndJ);
                vMaxH = simde_mm_blendv_epi8(vMaxH, vWH, cond_all);
                vEndI = simde_mm_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = simde_mm_blendv_epi8(vEndJ, vJ, cond_all);
                cond_all = simde_mm_and_si128(cond_Jlt, cond_eq);
                cond_all = simde_mm_and_si128(cond_all, cond_valid_IJ);
                vEndI = simde_mm_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = simde_mm_blendv_epi8(vEndJ, vJ, cond_all);
            }
            vJ = simde_mm_add_epi32(vJ, vOne);
        }
        vI = simde_mm_add_epi32(vI, vN);
    }

    /* alignment ending position */
    {
        int32_t *t = (int32_t*)&vMaxH;
        int32_t *i = (int32_t*)&vEndI;
        int32_t *j = (int32_t*)&vEndJ;
        int32_t k;
        for (k=0; k<N; ++k, ++t, ++i, ++j) {
            if (*t > score) {
                score = *t;
                end_query = *i;
                end_ref = *j;
            }
            else if (*t == score) {
                if (*j < end_ref) {
                    end_query = *i;
                    end_ref = *j;
                }
                else if (*j == end_ref && *i < end_query) {
                    end_query = *i;
                    end_ref = *j;
                }
            }
        }
    }

    

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_32 | PARASAIL_FLAG_LANES_4;

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}


