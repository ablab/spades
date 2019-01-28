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

#ifdef PARASAIL_TABLE
#define ENAME parasail_nw_banded_table
#else
#ifdef PARASAIL_ROWCOL
#define ENAME parasail_nw_banded_rowcol
#else
#define ENAME parasail_nw_banded
#endif
#endif

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const int k,
        const parasail_matrix_t *matrix)
{
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    int * const restrict s1 = parasail_memalign_int(16, s1Len);
    int * const restrict s2 = parasail_memalign_int(16, s2Len);
    int * HB = NULL;
    int * H = NULL;
    int * EB = NULL;
    int * E = NULL;
    int i = 0;
    int j = 0;
    int low = 0;
    int colLen = 0;
    int colOff = 0;

    low = s1Len - s2Len - k;
    if (s1Len > s2Len) {
        low = s2Len - s1Len - k;
    }
    colLen = k-low+1;
    colOff = s1Len > s2Len ? -k : low;

    HB = parasail_memalign_int(16, colLen+2);
    H = HB+1;
    EB = parasail_memalign_int(16, colLen+2);
    E = EB+1;
    parasail_memset_int(HB, 0, colLen+2);
    parasail_memset_int(EB, 0, colLen+2);

    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }

#ifdef PARASAIL_TABLE
    /* fill the table with zeros */
    {
        int limit = s1Len*s2Len;
        for (i=0; i<limit; ++i) {
            result->tables->score_table[i] = 0;
        }
    }
#endif

    /* initialize H */
    /* initialize E */
    for (i=-colOff+1,j=0; i<colLen; ++i, ++j) {
        H[i] = -open - j*gap;
        E[i] = NEG_INF_32;
    }
    H[-colOff-1] = -open;
    E[-colOff-1] = NEG_INF_32;
    H[-colOff] = 0;
    E[-colOff] = NEG_INF_32;
    H[-1] = NEG_INF_32;
    E[-1] = NEG_INF_32;
    H[colLen] = NEG_INF_32;
    E[colLen] = NEG_INF_32;

    /* iter over db */
    for (j=0; j<s2Len; ++j) {
        /* substitution matrix row (really we wanted a column) */
        const int * const restrict matcol = &matrix->matrix[matrix->size*s2[j]];
        int F = NEG_INF_32;
        if (colOff < 0) {
            H[-colOff-1] = -open - j*gap;
        }
        /* iter over query */
        for (i=0; i<colLen; ++i) {
            int pos = 0;
            int H_dag = 0;
            int E_opn = 0;
            int E_ext = 0;
            int F_opn = 0;
            int F_ext = 0;
            pos = colOff+i;
            if (pos < 0 || pos >= s1Len) {
                continue;
            }
            H_dag = H[i] + matcol[s1[pos]];
            E_opn = H[i+1] - open;
            E_ext = E[i+1] - gap;
            F_opn = H[i-1] - open;
            F_ext = F - gap;
            E[i] = MAX(E_opn, E_ext);
            F = MAX(F_opn, F_ext);
            H[i] = MAX(E[i], F);
            H[i] = MAX(H[i], H_dag);
#ifdef PARASAIL_TABLE
            result->tables->score_table[pos*s2Len + j] = H[i];
#endif
        }
        colOff += 1;
    }

    result->score = s1Len > s2Len ? H[-low] : H[k];
    result->end_query = s1Len-1;
    result->end_ref = s2Len-1;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_NOVEC
        | PARASAIL_FLAG_BANDED
        | PARASAIL_FLAG_BITS_INT | PARASAIL_FLAG_LANES_1;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(EB);
    parasail_free(HB);
    parasail_free(s2);
    parasail_free(s1);
    
    return result;
}

