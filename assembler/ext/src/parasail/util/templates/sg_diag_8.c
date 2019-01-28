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

#ifdef PARASAIL_TABLE
static inline void arr_store_si%(BITS)s(
        int *array,
        %(VTYPE)s vWH,
        %(INDEX)s i,
        %(INDEX)s s1Len,
        %(INDEX)s j,
        %(INDEX)s s2Len)
{
%(PRINTER)s
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        %(VTYPE)s vWH,
        %(INDEX)s i,
        %(INDEX)s s1Len,
        %(INDEX)s j,
        %(INDEX)s s2Len)
{
%(PRINTER_ROWCOL)s
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME %(NAME_TABLE)s
#else
#ifdef PARASAIL_ROWCOL
#define FNAME %(NAME_ROWCOL)s
#else
#define FNAME %(NAME)s
#endif
#endif

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
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    %(INT)s score = NEG_INF;
    %(VTYPE)s vNegInf = %(VSET1)s(NEG_INF);
    %(VTYPE)s vNegInf0 = %(VRSHIFT)s(vNegInf, %(BYTES)s); /* shift in a 0 */
    %(VTYPE)s vOpen = %(VSET1)s(open);
    %(VTYPE)s vGap  = %(VSET1)s(gap);
    %(VTYPE)s vOne16 = %(VSET1x16)s(1);
    %(VTYPE)s vN16 = %(VSET1x16)s(N);
    %(VTYPE)s vNegOne16 = %(VSET1x16)s(-1);
    %(VTYPE)s vILo16 = %(VSETx16)s(%(DIAG_ILO)s);
    %(VTYPE)s vIHi16 = %(VSETx16)s(%(DIAG_IHI)s);
    %(VTYPE)s vJresetLo16 = %(VSETx16)s(%(DIAG_JLO)s);
    %(VTYPE)s vJresetHi16 = %(VSETx16)s(%(DIAG_JHI)s);
    %(VTYPE)s vMaxH = vNegInf;
    %(VTYPE)s vEndILo = vNegInf;
    %(VTYPE)s vEndIHi = vNegInf;
    %(VTYPE)s vEndJLo = vNegInf;
    %(VTYPE)s vEndJHi = vNegInf;
    %(VTYPE)s vILimit16 = %(VSET1x16)s(s1Len);
    %(VTYPE)s vILimit116 = %(VSUBx16)s(vILimit16, vOne16);
    %(VTYPE)s vJLimit16 = %(VSET1x16)s(s2Len);
    %(VTYPE)s vJLimit116 = %(VSUBx16)s(vJLimit16, vOne16);
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
        %(VTYPE)s vNH = vNegInf0;
        %(VTYPE)s vWH = vNegInf0;
        %(VTYPE)s vE = vNegInf;
        %(VTYPE)s vF = vNegInf;
        %(VTYPE)s vJLo16 = vJresetLo16;
        %(VTYPE)s vJHi16 = vJresetHi16;
        %(DIAG_MATROW_DECL)s
        %(VTYPE)s vIltLimit = %(VPACKS)s(
                %(VCMPLTx16)s(vILo16, vILimit16),
                %(VCMPLTx16)s(vIHi16, vILimit16));
        %(VTYPE)s vIeqLimit1 = %(VPACKS)s(
                %(VCMPEQx16)s(vILo16, vILimit116),
                %(VCMPEQx16)s(vIHi16, vILimit116));
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            %(VTYPE)s vMat;
            %(VTYPE)s vNWH = vNH;
            vNH = %(VRSHIFT)s(vWH, %(BYTES)s);
            vNH = %(VINSERT)s(vNH, H_pr[j], %(LAST_POS)s);
            vF = %(VRSHIFT)s(vF, %(BYTES)s);
            vF = %(VINSERT)s(vF, F_pr[j], %(LAST_POS)s);
            vF = %(VMAX)s(
                    %(VSUB)s(vNH, vOpen),
                    %(VSUB)s(vF, vGap));
            vE = %(VMAX)s(
                    %(VSUB)s(vWH, vOpen),
                    %(VSUB)s(vE, vGap));
            vMat = %(VSET)s(
                    %(DIAG_MATROW_USE)s
                    );
            vNWH = %(VADD)s(vNWH, vMat);
            vWH = %(VMAX)s(vNWH, vE);
            vWH = %(VMAX)s(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                %(VTYPE)s cond = %(VPACKS)s(
                        %(VCMPEQx16)s(vJLo16,vNegOne16),
                        %(VCMPEQx16)s(vJHi16,vNegOne16));
                vWH = %(VANDNOT)s(cond, vWH);
                vF = %(VBLEND)s(vF, vNegInf, cond);
                vE = %(VBLEND)s(vE, vNegInf, cond);
            }
            %(SATURATION_CHECK_MID)s
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->tables->score_table, vWH, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->rowcols->score_row, result->rowcols->score_col, vWH, i, s1Len, j, s2Len);
#endif
            H_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vWH,0);
            F_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vF,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                %(VTYPE)s vJeqLimit1 = %(VPACKS)s(
                        %(VCMPEQx16)s(vJLo16, vJLimit116),
                        %(VCMPEQx16)s(vJHi16, vJLimit116));
                %(VTYPE)s vJgtNegOne = %(VPACKS)s(
                        %(VCMPGTx16)s(vJLo16, vNegOne16),
                        %(VCMPGTx16)s(vJHi16, vNegOne16));
                %(VTYPE)s vJltLimit = %(VPACKS)s(
                        %(VCMPLTx16)s(vJLo16, vJLimit16),
                        %(VCMPLTx16)s(vJHi16, vJLimit16));
                %(VTYPE)s cond_j = %(VAND)s(vIltLimit, vJeqLimit1);
                %(VTYPE)s cond_i = %(VAND)s(vIeqLimit1,
                        %(VAND)s(vJgtNegOne, vJltLimit));
                %(VTYPE)s cond_valid_IJ = %(VOR)s(cond_i, cond_j);
                %(VTYPE)s cond_eq = %(VCMPEQ)s(vWH, vMaxH);
                %(VTYPE)s cond_max = %(VCMPGT)s(vWH, vMaxH);
                %(VTYPE)s cond_all = %(VAND)s(cond_max, cond_valid_IJ);
                %(VTYPE)s cond_Jlt = %(VPACKS)s(
                        %(VCMPLTx16)s(vJLo16, vEndJLo),
                        %(VCMPLTx16)s(vJHi16, vEndJHi));
                %(VTYPE)s cond_lo = %(VUNPACKLO)s(cond_all, cond_all);
                %(VTYPE)s cond_hi = %(VUNPACKHI)s(cond_all, cond_all);
                vMaxH = %(VBLEND)s(vMaxH, vWH, cond_all);
                vEndILo = %(VBLEND)s(vEndILo, vILo16, cond_lo);
                vEndIHi = %(VBLEND)s(vEndIHi, vIHi16, cond_hi);
                vEndJLo = %(VBLEND)s(vEndJLo, vJLo16, cond_lo);
                vEndJHi = %(VBLEND)s(vEndJHi, vJHi16, cond_hi);
                cond_all = %(VAND)s(cond_Jlt, cond_eq);
                cond_all = %(VAND)s(cond_all, cond_valid_IJ);
                cond_lo = %(VUNPACKLO)s(cond_all, cond_all);
                cond_hi = %(VUNPACKHI)s(cond_all, cond_all);
                vEndILo = %(VBLEND)s(vEndILo, vILo16, cond_lo);
                vEndIHi = %(VBLEND)s(vEndIHi, vIHi16, cond_hi);
                vEndJLo = %(VBLEND)s(vEndJLo, vJLo16, cond_lo);
                vEndJHi = %(VBLEND)s(vEndJHi, vJHi16, cond_hi);
            }
            vJLo16 = %(VADDx16)s(vJLo16, vOne16);
            vJHi16 = %(VADDx16)s(vJHi16, vOne16);
        }
        vILo16 = %(VADDx16)s(vILo16, vN16);
        vIHi16 = %(VADDx16)s(vIHi16, vN16);
    }

    /* alignment ending position */
    {
        %(INT)s *t = (%(INT)s*)&vMaxH;
        int16_t *ilo = (int16_t*)&vEndILo;
        int16_t *jlo = (int16_t*)&vEndJLo;
        int16_t *ihi = (int16_t*)&vEndIHi;
        int16_t *jhi = (int16_t*)&vEndJHi;
        %(INDEX)s k;
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

    %(SATURATION_CHECK_FINAL)s

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;
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

