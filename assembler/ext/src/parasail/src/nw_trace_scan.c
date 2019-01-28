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

#define NEG_INF_32 (INT32_MIN/2)
#define MAX(a,b) ((a)>(b)?(a):(b))

#define ENAME parasail_nw_trace_scan

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_result_t *result = parasail_result_new_trace(s1Len, s2Len, 16, sizeof(int8_t));
    int * const restrict s1 = parasail_memalign_int(16, s1Len);
    int * const restrict s2 = parasail_memalign_int(16, s2Len);
    int * const restrict HB = parasail_memalign_int(16, s1Len+1);
    int * const restrict H  = HB+1;
    int * const restrict E  = parasail_memalign_int(16, s1Len);
    int * const restrict HtB= parasail_memalign_int(16, s1Len+1);
    int * const restrict Ht = HtB+1;
    int8_t * const restrict HT = (int8_t* const restrict)result->trace->trace_table;
    int * const restrict Ex = parasail_memalign_int(16, s1Len);
    int i = 0;
    int j = 0;

    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }

    /* initialize H */
    H[-1] = 0;
    Ht[-1] = 0;
    for (i=0; i<s1Len; ++i) {
        H[i] = -open - i*gap;
    }

    /* initialize E */
    for (i=0; i<s1Len; ++i) {
        E[i] = NEG_INF_32;
    }

    /* iterate over database */
    for (j=0; j<s2Len; ++j) {
        int Ft = NEG_INF_32;
        const int * const restrict matcol = &matrix->matrix[matrix->size*s2[j]];
        /* calculate E */
        for (i=0; i<s1Len; ++i) {
            int E_opn = H[i]-open;
            int E_ext = E[i]-gap;
            E[i] = MAX(E_ext, E_opn);
            HT[1LL*i*s2Len + j] = (E_opn > E_ext) ? PARASAIL_DIAG_E
                                                  : PARASAIL_INS_E;
        }
        /* calculate Ht */
        for (i=0; i<s1Len; ++i) {
            int H_dag = H[i-1]+matcol[s1[i]];
            Ht[i] = MAX(H_dag, E[i]);
            Ex[i] = (E[i] > H_dag);
        }
        Ht[-1] = -open -j*gap;
        /* calculate H */
        for (i=0; i<s1Len; ++i) {
            int Ft_opn;
            int Ht_pre = Ht[i-1];
            int Ft_ext = Ft-gap;
            if (Ht_pre >= Ft_ext) {
                Ft = Ht_pre;
            }
            else {
                Ft = Ft_ext;
            }
            Ft_opn = Ft-open;
            if (H[i-1] > Ft_ext) {
                HT[1LL*i*s2Len + j] |= PARASAIL_DIAG_F;
            }
            else {
                HT[1LL*i*s2Len + j] |= PARASAIL_DEL_F;
            }
            if (Ht[i] > Ft_opn) {
                H[i] = Ht[i];
                HT[1LL*i*s2Len + j] |= Ex[i] ? PARASAIL_INS : PARASAIL_DIAG;
            }
            else {
                H[i] = Ft_opn;
                if (Ht[i] == Ft_opn) {
                    if (Ex[i]) {
                        HT[1LL*i*s2Len + j] |= PARASAIL_DEL;
                    }
                    else {
                        HT[1LL*i*s2Len + j] |= PARASAIL_DIAG;
                    }
                }
                else {
                    HT[1LL*i*s2Len + j] |= PARASAIL_DEL;
                }
            }
        }
        H[-1] = -open - j*gap;
    }

    result->score = H[s1Len-1];
    result->end_query = s1Len-1;
    result->end_ref = s2Len-1;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_NOVEC_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_INT | PARASAIL_FLAG_LANES_1;

    parasail_free(Ex);
    parasail_free(HtB);
    parasail_free(E);
    parasail_free(HB);
    parasail_free(s2);
    parasail_free(s1);

    return result;
}

