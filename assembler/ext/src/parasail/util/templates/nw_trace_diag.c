/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

%(HEADER)s

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_%(ISA)s.h"

#define NEG_INF %(NEG_INF)s
%(FIXES)s

static inline void arr_store_si%(BITS)s(
        int8_t *array,
        %(VTYPE)s vWH,
        %(INDEX)s i,
        %(INDEX)s s1Len,
        %(INDEX)s j,
        %(INDEX)s s2Len)
{
%(PRINTER_TRACE)s
}

#define FNAME %(NAME_TRACE)s

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const %(INDEX)s N = %(LANES)s; /* number of values in vector */
    const %(INDEX)s PAD = N-1;
    const %(INDEX)s PAD2 = PAD*2;
    const %(INDEX)s s1Len_PAD = s1Len+PAD;
    const %(INDEX)s s2Len_PAD = s2Len+PAD;
    %(INT)s * const restrict s1 = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s1Len+PAD);
    %(INT)s * const restrict s2B= parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    %(INT)s * const restrict _H_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    %(INT)s * const restrict _F_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    %(INT)s * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    %(INT)s * const restrict H_pr = _H_pr+PAD;
    %(INT)s * const restrict F_pr = _F_pr+PAD;
    parasail_result_t *result = parasail_result_new_trace(s1Len, s2Len, %(ALIGNMENT)s, sizeof(int8_t));
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s end_query = s1Len-1;
    %(INDEX)s end_ref = s2Len-1;
    %(INT)s score = NEG_INF;
    %(VTYPE)s vNegInf = %(VSET1)s(NEG_INF);
    %(VTYPE)s vOpen = %(VSET1)s(open);
    %(VTYPE)s vGap  = %(VSET1)s(gap);
    %(VTYPE)s vOne = %(VSET1)s(1);
    %(VTYPE)s vN = %(VSET1)s(N);
    %(VTYPE)s vGapN = %(VSET1)s(gap*N);
    %(VTYPE)s vNegOne = %(VSET1)s(-1);
    %(VTYPE)s vI = %(VSET)s(%(DIAG_I)s);
    %(VTYPE)s vJreset = %(VSET)s(%(DIAG_J)s);
    %(VTYPE)s vMax = vNegInf;
    %(VTYPE)s vILimit = %(VSET1)s(s1Len);
    %(VTYPE)s vILimit1 = %(VSUB)s(vILimit, vOne);
    %(VTYPE)s vJLimit = %(VSET1)s(s2Len);
    %(VTYPE)s vJLimit1 = %(VSUB)s(vJLimit, vOne);
    %(VTYPE)s vIBoundary = %(VSET)s(
            %(DIAG_IBoundary)s
            );
    %(VTYPE)s vTDiag = %(VSET1)s(PARASAIL_DIAG);
    %(VTYPE)s vTIns = %(VSET1)s(PARASAIL_INS);
    %(VTYPE)s vTDel = %(VSET1)s(PARASAIL_DEL);
    %(VTYPE)s vTDiagE = %(VSET1)s(PARASAIL_DIAG_E);
    %(VTYPE)s vTInsE = %(VSET1)s(PARASAIL_INS_E);
    %(VTYPE)s vTDiagF = %(VSET1)s(PARASAIL_DIAG_F);
    %(VTYPE)s vTDelF = %(VSET1)s(PARASAIL_DEL_F);
    %(SATURATION_CHECK_INIT)s

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
        %(VTYPE)s vNH = vNegInf;
        %(VTYPE)s vWH = vNegInf;
        %(VTYPE)s vE = vNegInf;
        %(VTYPE)s vE_opn = vNegInf;
        %(VTYPE)s vE_ext = vNegInf;
        %(VTYPE)s vF = vNegInf;
        %(VTYPE)s vF_opn = vNegInf;
        %(VTYPE)s vF_ext = vNegInf;
        %(VTYPE)s vJ = vJreset;
        %(DIAG_MATROW_DECL)s
        vNH = %(VRSHIFT)s(vNH, %(BYTES)s);
        vNH = %(VINSERT)s(vNH, H_pr[-1], %(LAST_POS)s);
        vWH = %(VRSHIFT)s(vWH, %(BYTES)s);
        vWH = %(VINSERT)s(vWH, -open - i*gap, %(LAST_POS)s);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            %(VTYPE)s vMat;
            %(VTYPE)s vNWH = vNH;
            vNH = %(VRSHIFT)s(vWH, %(BYTES)s);
            vNH = %(VINSERT)s(vNH, H_pr[j], %(LAST_POS)s);
            vF = %(VRSHIFT)s(vF, %(BYTES)s);
            vF = %(VINSERT)s(vF, F_pr[j], %(LAST_POS)s);
            vF_opn = %(VSUB)s(vNH, vOpen);
            vF_ext = %(VSUB)s(vF, vGap);
            vF = %(VMAX)s(vF_opn, vF_ext);
            vE_opn = %(VSUB)s(vWH, vOpen);
            vE_ext = %(VSUB)s(vE, vGap);
            vE = %(VMAX)s(vE_opn, vE_ext);
            vMat = %(VSET)s(
                    %(DIAG_MATROW_USE)s
                    );
            vNWH = %(VADD)s(vNWH, vMat);
            vWH = %(VMAX)s(vNWH, vE);
            vWH = %(VMAX)s(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                %(VTYPE)s cond = %(VCMPEQ)s(vJ,vNegOne);
                vWH = %(VBLEND)s(vWH, vIBoundary, cond);
                vF = %(VBLEND)s(vF, vNegInf, cond);
                vE = %(VBLEND)s(vE, vNegInf, cond);
            }
            %(SATURATION_CHECK_MID)s
            /* trace table */
            {
                %(VTYPE)s case1 = %(VCMPEQ)s(vWH, vNWH);
                %(VTYPE)s case2 = %(VCMPEQ)s(vWH, vF);
                %(VTYPE)s vT = %(VBLEND)s(
                        %(VBLEND)s(vTIns, vTDel, case2),
                        vTDiag,
                        case1);
                %(VTYPE)s condE = %(VCMPGT)s(vE_opn, vE_ext);
                %(VTYPE)s condF = %(VCMPGT)s(vF_opn, vF_ext);
                %(VTYPE)s vET = %(VBLEND)s(vTInsE, vTDiagE, condE);
                %(VTYPE)s vFT = %(VBLEND)s(vTDelF, vTDiagF, condF);
                vT = %(VOR)s(vT, vET);
                vT = %(VOR)s(vT, vFT);
                arr_store_si%(BITS)s(result->trace->trace_table, vT, i, s1Len, j, s2Len);
            }
            H_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vWH,0);
            F_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vF,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                %(VTYPE)s cond_valid_I = %(VCMPEQ)s(vI, vILimit1);
                %(VTYPE)s cond_valid_J = %(VCMPEQ)s(vJ, vJLimit1);
                %(VTYPE)s cond_all = %(VAND)s(cond_valid_I, cond_valid_J);
                vMax = %(VBLEND)s(vMax, vWH, cond_all);
            }
            vJ = %(VADD)s(vJ, vOne);
        }
        vI = %(VADD)s(vI, vN);
        vIBoundary = %(VSUB)s(vIBoundary, vGapN);
    }

    /* max in vMax */
    for (i=0; i<N; ++i) {
        %(INT)s value;
        value = (%(INT)s) %(VEXTRACT)s(vMax, %(LAST_POS)s);
        if (value > score) {
            score = value;
        }
        vMax = %(VSHIFT)s(vMax, %(BYTES)s);
    }

    %(SATURATION_CHECK_FINAL)s

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

