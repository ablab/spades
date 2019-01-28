/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <emmintrin.h>
#include <smmintrin.h>
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_sse.h"

#define NEG_INF INT8_MIN


static inline void arr_store_si128(
        int8_t *array,
        __m128i vWH,
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

#define FNAME parasail_sw_trace_diag_sse41_128_8

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
    int8_t * const restrict s1 = parasail_memalign_int8_t(16, s1Len+PAD);
    int8_t * const restrict s2B= parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _H_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _F_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int8_t * const restrict H_pr = _H_pr+PAD;
    int8_t * const restrict F_pr = _F_pr+PAD;
    parasail_result_t *result = parasail_result_new_trace(s1Len, s2Len, 16, sizeof(int8_t));
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int8_t score = NEG_INF;
    __m128i vNegInf = _mm_set1_epi8(NEG_INF);
    __m128i vNegInf0 = _mm_srli_si128(vNegInf, 1); /* shift in a 0 */
    __m128i vOpen = _mm_set1_epi8(open);
    __m128i vGap  = _mm_set1_epi8(gap);
    __m128i vZero = _mm_set1_epi8(0);
    __m128i vOne16 = _mm_set1_epi16(1);
    __m128i vNegOne16 = _mm_set1_epi16(-1);
    __m128i vN16 = _mm_set1_epi16(N);
    __m128i vILo16 = _mm_set_epi16(8,9,10,11,12,13,14,15);
    __m128i vIHi16 = _mm_set_epi16(0,1,2,3,4,5,6,7);
    __m128i vJresetLo16 = _mm_set_epi16(-8,-9,-10,-11,-12,-13,-14,-15);
    __m128i vJresetHi16 = _mm_set_epi16(0,-1,-2,-3,-4,-5,-6,-7);
    __m128i vMaxH = vNegInf;
    __m128i vEndILo = vNegInf;
    __m128i vEndIHi = vNegInf;
    __m128i vEndJLo = vNegInf;
    __m128i vEndJHi = vNegInf;
    __m128i vILimit16 = _mm_set1_epi16(s1Len);
    __m128i vJLimit16 = _mm_set1_epi16(s2Len);
    __m128i vTDiag = _mm_set1_epi8(PARASAIL_DIAG);
    __m128i vTIns = _mm_set1_epi8(PARASAIL_INS);
    __m128i vTDel = _mm_set1_epi8(PARASAIL_DEL);
    __m128i vTZero = _mm_set1_epi8(PARASAIL_ZERO);
    __m128i vTDiagE = _mm_set1_epi8(PARASAIL_DIAG_E);
    __m128i vTInsE = _mm_set1_epi8(PARASAIL_INS_E);
    __m128i vTDiagF = _mm_set1_epi8(PARASAIL_DIAG_F);
    __m128i vTDelF = _mm_set1_epi8(PARASAIL_DEL_F);
    __m128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
    __m128i vSaturationCheckMin = vPosLimit;
    __m128i vSaturationCheckMax = vNegLimit;

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
        __m128i vNH = vNegInf0;
        __m128i vWH = vNegInf0;
        __m128i vE = vNegInf;
        __m128i vE_opn = vNegInf;
        __m128i vE_ext = vNegInf;
        __m128i vF = vNegInf;
        __m128i vF_opn = vNegInf;
        __m128i vF_ext = vNegInf;
        __m128i vJLo16 = vJresetLo16;
        __m128i vJHi16 = vJresetHi16;
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
        __m128i vIltLimit = _mm_packs_epi16(
                    _mm_cmplt_epi16(vILo16, vILimit16),
                    _mm_cmplt_epi16(vIHi16, vILimit16));
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWH = vNH;
            vNH = _mm_srli_si128(vWH, 1);
            vNH = _mm_insert_epi8(vNH, H_pr[j], 15);
            vF = _mm_srli_si128(vF, 1);
            vF = _mm_insert_epi8(vF, F_pr[j], 15);
            vF_opn = _mm_subs_epi8(vNH, vOpen);
            vF_ext = _mm_subs_epi8(vF, vGap);
            vF = _mm_max_epi8(vF_opn, vF_ext);
            vE_opn = _mm_subs_epi8(vWH, vOpen);
            vE_ext = _mm_subs_epi8(vE, vGap);
            vE = _mm_max_epi8(vE_opn, vE_ext);
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
            vNWH = _mm_max_epi8(vNWH, vZero);
            vWH = _mm_max_epi8(vNWH, vE);
            vWH = _mm_max_epi8(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_packs_epi16(
                        _mm_cmpeq_epi16(vJLo16,vNegOne16),
                        _mm_cmpeq_epi16(vJHi16,vNegOne16));
                vWH = _mm_andnot_si128(cond, vWH);
                vF = _mm_blendv_epi8(vF, vNegInf, cond);
                vE = _mm_blendv_epi8(vE, vNegInf, cond);
            }
            /* check for saturation */
            {
                vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vWH);
                vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vWH);
            }
            /* trace table */
            {
                __m128i cond_zero = _mm_cmpeq_epi8(vWH, vZero);
                __m128i case1 = _mm_cmpeq_epi8(vWH, vNWH);
                __m128i case2 = _mm_cmpeq_epi8(vWH, vF);
                __m128i vT = _mm_blendv_epi8(
                        _mm_blendv_epi8(vTIns, vTDel, case2),
                        _mm_blendv_epi8(vTDiag, vTZero, cond_zero),
                        case1);
                __m128i condE = _mm_cmpgt_epi8(vE_opn, vE_ext);
                __m128i condF = _mm_cmpgt_epi8(vF_opn, vF_ext);
                __m128i vET = _mm_blendv_epi8(vTInsE, vTDiagE, condE);
                __m128i vFT = _mm_blendv_epi8(vTDelF, vTDiagF, condF);
                vT = _mm_or_si128(vT, vET);
                vT = _mm_or_si128(vT, vFT);
                arr_store_si128(result->trace->trace_table, vT, i, s1Len, j, s2Len);
            }
            H_pr[j-15] = (int8_t)_mm_extract_epi8(vWH,0);
            F_pr[j-15] = (int8_t)_mm_extract_epi8(vF,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                __m128i cond_valid_J = _mm_and_si128(
                        _mm_packs_epi16(
                            _mm_cmpgt_epi16(vJLo16, vNegOne16),
                            _mm_cmpgt_epi16(vJHi16, vNegOne16)),
                        _mm_packs_epi16(
                            _mm_cmplt_epi16(vJLo16, vJLimit16),
                            _mm_cmplt_epi16(vJHi16, vJLimit16)));
                __m128i cond_valid_IJ = _mm_and_si128(cond_valid_J, vIltLimit);
                __m128i cond_eq = _mm_cmpeq_epi8(vWH, vMaxH);
                __m128i cond_max = _mm_cmpgt_epi8(vWH, vMaxH);
                __m128i cond_all = _mm_and_si128(cond_max, cond_valid_IJ);
                __m128i cond_Jlt = _mm_packs_epi16(
                        _mm_cmplt_epi16(vJLo16, vEndJLo),
                        _mm_cmplt_epi16(vJHi16, vEndJHi));
                __m128i cond_lo = _mm_unpacklo_epi8(cond_all, cond_all);
                __m128i cond_hi = _mm_unpackhi_epi8(cond_all, cond_all);
                vMaxH = _mm_blendv_epi8(vMaxH, vWH, cond_all);
                vEndILo = _mm_blendv_epi8(vEndILo, vILo16, cond_lo);
                vEndIHi = _mm_blendv_epi8(vEndIHi, vIHi16, cond_hi);
                vEndJLo = _mm_blendv_epi8(vEndJLo, vJLo16, cond_lo);
                vEndJHi = _mm_blendv_epi8(vEndJHi, vJHi16, cond_hi);
                cond_all = _mm_and_si128(cond_Jlt, cond_eq);
                cond_all = _mm_and_si128(cond_all, cond_valid_IJ);
                cond_lo = _mm_unpacklo_epi8(cond_all, cond_all);
                cond_hi = _mm_unpackhi_epi8(cond_all, cond_all);
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
        int16_t *ilo = (int16_t*)&vEndILo;
        int16_t *jlo = (int16_t*)&vEndJLo;
        int16_t *ihi = (int16_t*)&vEndIHi;
        int16_t *jhi = (int16_t*)&vEndJHi;
        int32_t k;
        for (k=0; k<N/2; ++k, ++t, ++ilo, ++jlo) {
            if (*t > score) {
                score = *t;
                end_query = *ilo;
                end_ref = *jlo;
            }
            else if (*t == score) {
                if (*jlo < end_ref) {
                    end_query = *ilo;
                    end_ref = *jlo;
                }
                else if (*jlo == end_ref && *ilo < end_query) {
                    end_query = *ilo;
                    end_ref = *jlo;
                }
            }
        }
        for (k=N/2; k<N; ++k, ++t, ++ihi, ++jhi) {
            if (*t > score) {
                score = *t;
                end_query = *ihi;
                end_ref = *jhi;
            }
            else if (*t == score) {
                if (*jhi < end_ref) {
                    end_query = *ihi;
                    end_ref = *jhi;
                }
                else if (*jhi == end_ref && *ihi < end_query) {
                    end_query = *ihi;
                    end_ref = *jhi;
                }
            }
        }
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            _mm_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}


