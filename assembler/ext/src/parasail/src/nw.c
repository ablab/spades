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
#define ENAME parasail_nw_table
#else
#ifdef PARASAIL_ROWCOL
#define ENAME parasail_nw_rowcol
#else
#define ENAME parasail_nw
#endif
#endif

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
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
    int * const restrict H = parasail_memalign_int(16, s2Len+1);
    int * const restrict F = parasail_memalign_int(16, s2Len+1);
    int i = 0;
    int j = 0;

    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }

    /* upper left corner */
    H[0] = 0;
    F[0] = NEG_INF_32;
    
    /* first row */
    for (j=1; j<=s2Len; ++j) {
        H[j] = -open -(j-1)*gap;
        F[j] = NEG_INF_32;
    }

    /* iter over first sequence */
    for (i=1; i<=s1Len; ++i) {
        const int * const restrict matrow = &matrix->matrix[matrix->size*s1[i-1]];
        /* init first column */
        int NH = H[0];
        int WH = -open - (i-1)*gap;
        int E = NEG_INF_32;
        H[0] = WH;
        for (j=1; j<=s2Len; ++j) {
            int H_dag;
            int H_new;
            int E_opn;
            int E_ext;
            int F_opn;
            int F_ext;
            int NWH = NH;
            NH = H[j];
            F_opn = NH - open;
            F_ext = F[j] - gap;
            F[j] = MAX(F_opn, F_ext);
            E_opn = WH - open;
            E_ext = E    - gap;
            E    = MAX(E_opn, E_ext);
            H_dag = NWH + matrow[s2[j-1]];
            H_new = MAX(H_dag, E);
            H_new = MAX(H_new, F[j]);
            WH = H_new;
            H[j] = WH;
#ifdef PARASAIL_TABLE
            result->tables->score_table[1LL*(i-1)*s2Len + (j-1)] = WH;
#endif
        }
#ifdef PARASAIL_ROWCOL
        result->rowcols->score_col[i-1] = WH;
#endif
    }
#ifdef PARASAIL_ROWCOL
    for (j=1; j<=s2Len; ++j) {
        result->rowcols->score_row[j-1] = H[j];
    }
#endif

    result->score = H[s2Len];
    result->end_query = s1Len-1;
    result->end_ref = s2Len-1;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_NOVEC
        | PARASAIL_FLAG_BITS_INT | PARASAIL_FLAG_LANES_1;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(F);
    parasail_free(H);
    parasail_free(s2);
    parasail_free(s1);
    
    return result;
}

