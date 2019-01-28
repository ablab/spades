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
#include <string.h>

%(HEADER)s

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_%(ISA)s.h"

%(FIXES)s

static inline void arr_store(
        %(VTYPE)s *array,
        %(VTYPE)s vH,
        %(INDEX)s t,
        %(INDEX)s seglen,
        %(INDEX)s d)
{
    %(VSTORE)s(array + (1LL*d*seglen+t), vH);
}

static inline %(VTYPE)s arr_load(
        %(VTYPE)s *array,
        %(INDEX)s t,
        %(INDEX)s seglen,
        %(INDEX)s d)
{
    return %(VLOAD)s(array + (1LL*d*seglen+t));
}

#define FNAME %(NAME_TRACE)s
#define PNAME %(PNAME_TRACE)s

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_%(ISA)s_%(BITS)s_%(WIDTH)s(s1, s1Len, matrix);
    parasail_result_t *result = PNAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const %(INDEX)s segWidth = %(LANES)s; /* number of values in vector unit */
    const %(INDEX)s segLen = (s1Len + segWidth - 1) / segWidth;
    %(VTYPE)s* const restrict pvP  = (%(VTYPE)s*)profile->profile%(WIDTH)s.score;
    %(VTYPE)s* const restrict pvE  = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHt = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvH  = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHMax  = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvGapper = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s vGapO = %(VSET1)s(open);
    %(VTYPE)s vGapE = %(VSET1)s(gap);
    const %(INT)s NEG_LIMIT = (-open < matrix->min ?
        INT%(WIDTH)s_MIN + open : INT%(WIDTH)s_MIN - matrix->min) + 1;
    const %(INT)s POS_LIMIT = INT%(WIDTH)s_MAX - matrix->max - 1;
    %(VTYPE)s vZero = %(VSET0)s();
    %(INT)s score = NEG_LIMIT;
    %(VTYPE)s vNegLimit = %(VSET1)s(NEG_LIMIT);
    %(VTYPE)s vPosLimit = %(VSET1)s(POS_LIMIT);
    %(VTYPE)s vSaturationCheckMin = vPosLimit;
    %(VTYPE)s vSaturationCheckMax = vNegLimit;
    %(VTYPE)s vMaxH = vNegLimit;
    %(VTYPE)s vMaxHUnit = vNegLimit;
    %(VTYPE)s vNegInfFront = vZero;
    %(VTYPE)s vSegLenXgap;
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, %(ALIGNMENT)s, sizeof(%(VTYPE)s));
    %(VTYPE)s vTZero = %(VSET1)s(PARASAIL_ZERO);
    %(VTYPE)s vTIns  = %(VSET1)s(PARASAIL_INS);
    %(VTYPE)s vTDel  = %(VSET1)s(PARASAIL_DEL);
    %(VTYPE)s vTDiag = %(VSET1)s(PARASAIL_DIAG);
    %(VTYPE)s vTDiagE = %(VSET1)s(PARASAIL_DIAG_E);
    %(VTYPE)s vTInsE = %(VSET1)s(PARASAIL_INS_E);
    %(VTYPE)s vTDiagF = %(VSET1)s(PARASAIL_DIAG_F);
    %(VTYPE)s vTDelF = %(VSET1)s(PARASAIL_DEL_F);

    vNegInfFront = %(VINSERT)s(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = %(VADD)s(vNegInfFront,
            %(VSHIFT)s(%(VSET1)s(-segLen*gap), %(BYTES)s));

    parasail_memset_%(VTYPE)s(pvH, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvE, vNegLimit, segLen);
    {
        %(VTYPE)s vGapper = %(VSUB)s(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            %(VSTORE)s(pvGapper+i, vGapper);
            vGapper = %(VSUB)s(vGapper, vGapE);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        %(VTYPE)s vE;
        %(VTYPE)s vE_ext;
        %(VTYPE)s vE_opn;
        %(VTYPE)s vHt;
        %(VTYPE)s vF;
        %(VTYPE)s vF_ext;
        %(VTYPE)s vF_opn;
        %(VTYPE)s vH;
        %(VTYPE)s vHp;
        %(VTYPE)s *pvW;
        %(VTYPE)s vW;
        %(VTYPE)s case1;
        %(VTYPE)s case2;
        %(VTYPE)s case0;
        %(VTYPE)s vGapper;
        %(VTYPE)s vT;
        %(VTYPE)s vET;
        %(VTYPE)s vFT;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = %(VLOAD)s(pvH+(segLen-1));
        vHp = %(VSHIFT)s(vHp, %(BYTES)s);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vZero;
        vF = vNegLimit;
        for (i=0; i<segLen; ++i) {
            vH = %(VLOAD)s(pvH+i);
            vE = %(VLOAD)s(pvE+i);
            vW = %(VLOAD)s(pvW+i);
            vGapper = %(VLOAD)s(pvGapper+i);
            vE_opn = %(VSUB)s(vH, vGapO);
            vE_ext = %(VSUB)s(vE, vGapE);
            case1 = %(VCMPGT)s(vE_opn, vE_ext);
            vET = %(VBLEND)s(vTInsE, vTDiagE, case1);
            arr_store(result->trace->trace_table, vET, i, segLen, j);
            vE = %(VMAX)s(vE_opn, vE_ext);
            vGapper = %(VADD)s(vHt, vGapper);
            vF = %(VMAX)s(vF, vGapper);
            vHp = %(VADD)s(vHp, vW);
            vHt = %(VMAX)s(vE, vHp);
            %(VSTORE)s(pvE+i, vE);
            %(VSTORE)s(pvHt+i, vHt);
            %(VSTORE)s(pvH+i, vHp);
            vHp = vH;
        }

        /* pseudo prefix scan on F and H */
        vHt = %(VSHIFT)s(vHt, %(BYTES)s);
        vGapper = %(VLOAD)s(pvGapper);
        vGapper = %(VADD)s(vHt, vGapper);
        vF = %(VMAX)s(vF, vGapper);
        for (i=0; i<segWidth-2; ++i) {
            %(VTYPE)s vFt = %(VSHIFT)s(vF, %(BYTES)s);
            vFt = %(VADD)s(vFt, vSegLenXgap);
            vF = %(VMAX)s(vF, vFt);
        }

        /* calculate final H */
        vF = %(VSHIFT)s(vF, %(BYTES)s);
        vF = %(VADD)s(vF, vNegInfFront);
        vH = %(VMAX)s(vF, vHt);
        vH = %(VMAX)s(vH, vZero);
        for (i=0; i<segLen; ++i) {
            vET = arr_load(result->trace->trace_table, i, segLen, j);
            vHp = %(VLOAD)s(pvH+i);
            vHt = %(VLOAD)s(pvHt+i);
            vF_opn = %(VSUB)s(vH, vGapO);
            vF_ext = %(VSUB)s(vF, vGapE);
            vF = %(VMAX)s(vF_opn, vF_ext);
            case1 = %(VCMPGT)s(vF_opn, vF_ext);
            vFT = %(VBLEND)s(vTDelF, vTDiagF, case1);
            vH = %(VMAX)s(vHt, vF);
            vH = %(VMAX)s(vH, vZero);
            case0 = %(VCMPEQ)s(vH, vZero);
            case1 = %(VCMPEQ)s(vH, vHp);
            case2 = %(VCMPEQ)s(vH, vF);
            vT = %(VBLEND)s(
                    %(VBLEND)s(vTIns, vTDel, case2),
                    vTDiag, case1);
            vT = %(VBLEND)s(vT, vTZero, case0);
            vT = %(VOR)s(vT, vET);
            vT = %(VOR)s(vT, vFT);
            arr_store(result->trace->trace_table, vT, i, segLen, j);
            %(VSTORE)s(pvH+i, vH);
            vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
            {
                vMaxH = %(VMAX)s(vH, vMaxH);
            }
        } 

        {
            %(VTYPE)s vCompare = %(VCMPGT)s(vMaxH, vMaxHUnit);
            if (%(VMOVEMASK)s(vCompare)) {
                score = %(VHMAX)s(vMaxH);
                vMaxHUnit = %(VSET1)s(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvH, sizeof(%(VTYPE)s)*segLen);
            }
        }
    }

    /* Trace the alignment ending position on read. */
    {
        %(INT)s *t = (%(INT)s*)pvHMax;
        %(INDEX)s column_len = segLen * segWidth;
        end_query = s1Len;
        for (i = 0; i<column_len; ++i, ++t) {
            if (*t == score) {
                %(INDEX)s temp = i / segWidth + i %% segWidth * segLen;
                if (temp < end_query) {
                    end_query = temp;
                }
            }
        }
    }

    if (%(VMOVEMASK)s(%(VOR)s(
            %(VCMPLT)s(vSaturationCheckMin, vNegLimit),
            %(VCMPGT)s(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;

    parasail_free(pvGapper);
    parasail_free(pvHMax);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}

