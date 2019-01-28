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

%(HEADER)s

#include "parasail/parasail.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/internal_%(ISA)s.h"

#define SWAP(A,B) { %(VTYPE)s* tmp = A; A = B; B = tmp; }

#define NEG_INF %(NEG_INF)s
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
    %(INDEX)s k = 0;
    const int s1Len = profile->s1Len;
    %(INDEX)s end_query = s1Len-1;
    %(INDEX)s end_ref = s2Len-1;
    const parasail_matrix_t *matrix = profile->matrix;
    const %(INDEX)s segWidth = %(LANES)s; /* number of values in vector unit */
    const %(INDEX)s segLen = (s1Len + segWidth - 1) / segWidth;
    const %(INDEX)s offset = (s1Len - 1) %% segLen;
    const %(INDEX)s position = (segWidth - 1) - (s1Len - 1) / segLen;
    %(VTYPE)s* const restrict vProfile = (%(VTYPE)s*)profile->profile%(WIDTH)s.score;
    %(VTYPE)s* restrict pvHStore = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLoad =  parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvE = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvEaStore = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvEaLoad = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHT = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(INT)s* const restrict boundary = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+1);
    %(VTYPE)s vGapO = %(VSET1)s(open);
    %(VTYPE)s vGapE = %(VSET1)s(gap);
    %(VTYPE)s vNegInf = %(VSET1)s(NEG_INF);
    %(INT)s score = NEG_INF;
    %(SATURATION_CHECK_INIT)s
    parasail_result_t *result = parasail_result_new_trace(segLen, s2Len, %(ALIGNMENT)s, sizeof(%(VTYPE)s));
    %(VTYPE)s vTIns  = %(VSET1)s(PARASAIL_INS);
    %(VTYPE)s vTDel  = %(VSET1)s(PARASAIL_DEL);
    %(VTYPE)s vTDiag = %(VSET1)s(PARASAIL_DIAG);
    %(VTYPE)s vTDiagE = %(VSET1)s(PARASAIL_DIAG_E);
    %(VTYPE)s vTInsE = %(VSET1)s(PARASAIL_INS_E);
    %(VTYPE)s vTDiagF = %(VSET1)s(PARASAIL_DIAG_F);
    %(VTYPE)s vTDelF = %(VSET1)s(PARASAIL_DEL_F);
    %(VTYPE)s vTMask = %(VSET1)s(PARASAIL_ZERO_MASK);
    %(VTYPE)s vFTMask = %(VSET1)s(PARASAIL_F_MASK);

%(INIT_H_AND_E)s

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT%(WIDTH)s_MIN ? INT%(WIDTH)s_MIN : tmp;
        }
    }

    for (i=0; i<segLen; ++i) {
        arr_store(result->trace->trace_table, vTDiagE, i, segLen, 0);
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        %(VTYPE)s vEF_opn;
        %(VTYPE)s vE;
        %(VTYPE)s vE_ext;
        %(VTYPE)s vF;
        %(VTYPE)s vF_ext;
        %(VTYPE)s vFa;
        %(VTYPE)s vFa_ext;
        %(VTYPE)s vH;
        %(VTYPE)s vH_dag;
        const %(VTYPE)s* vP = NULL;

        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vF = vNegInf;

        /* load final segment of pvHStore and shift left by %(BYTES)s bytes */
        vH = %(VLOAD)s(&pvHStore[segLen - 1]);
        vH = %(VSHIFT)s(vH, %(BYTES)s);

        /* insert upper boundary condition */
        vH = %(VINSERT)s(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad, pvHStore)
        SWAP(pvEaLoad, pvEaStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = %(VLOAD)s(pvE + i);

            /* Get max from vH, vE and vF. */
            vH_dag = %(VADD)s(vH, %(VLOAD)s(vP + i));
            vH = %(VMAX)s(vH_dag, vE);
            vH = %(VMAX)s(vH, vF);
            /* Save vH values. */
            %(VSTORE)s(pvHStore + i, vH);
            %(SATURATION_CHECK_MID)s

            {
                %(VTYPE)s vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                %(VTYPE)s case1 = %(VCMPEQ)s(vH, vH_dag);
                %(VTYPE)s case2 = %(VCMPEQ)s(vH, vF);
                %(VTYPE)s vT = %(VBLEND)s(
                        %(VBLEND)s(vTIns, vTDel, case2),
                        vTDiag, case1);
                %(VSTORE)s(pvHT + i, vT);
                vT = %(VOR)s(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i, segLen, j);
            }

            vEF_opn = %(VSUB)s(vH, vGapO);

            /* Update vE value. */
            vE_ext = %(VSUB)s(vE, vGapE);
            vE = %(VMAX)s(vEF_opn, vE_ext);
            %(VSTORE)s(pvE + i, vE);
            {
                %(VTYPE)s vEa = %(VLOAD)s(pvEaLoad + i);
                %(VTYPE)s vEa_ext = %(VSUB)s(vEa, vGapE);
                vEa = %(VMAX)s(vEF_opn, vEa_ext);
                %(VSTORE)s(pvEaStore + i, vEa);
                if (j+1<s2Len) {
                    %(VTYPE)s cond = %(VCMPGT)s(vEF_opn, vEa_ext);
                    %(VTYPE)s vT = %(VBLEND)s(vTInsE, vTDiagE, cond);
                    arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                }
            }

            /* Update vF value. */
            vF_ext = %(VSUB)s(vF, vGapE);
            vF = %(VMAX)s(vEF_opn, vF_ext);
            if (i+1<segLen) {
                %(VTYPE)s vTAll = arr_load(result->trace->trace_table, i+1, segLen, j);
                %(VTYPE)s cond = %(VCMPGT)s(vEF_opn, vF_ext);
                %(VTYPE)s vT = %(VBLEND)s(vTDelF, vTDiagF, cond);
                vT = %(VOR)s(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i+1, segLen, j);
            }

            /* Load the next vH. */
            vH = %(VLOAD)s(pvHLoad + i);
        }


        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            int64_t tmp = boundary[j+1]-open;
            %(INT)s tmp2 = tmp < INT%(WIDTH)s_MIN ? INT%(WIDTH)s_MIN : tmp;
            %(VTYPE)s vHp = %(VLOAD)s(&pvHLoad[segLen - 1]);
            vHp = %(VSHIFT)s(vHp, %(BYTES)s);
            vHp = %(VINSERT)s(vHp, boundary[j], 0);
            vEF_opn = %(VSHIFT)s(vEF_opn, %(BYTES)s);
            vEF_opn = %(VINSERT)s(vEF_opn, tmp2, 0);
            vF_ext = %(VSHIFT)s(vF_ext, %(BYTES)s);
            vF_ext = %(VINSERT)s(vF_ext, NEG_INF, 0);
            vF = %(VSHIFT)s(vF, %(BYTES)s);
            vF = %(VINSERT)s(vF, tmp2, 0);
            vFa_ext = %(VSHIFT)s(vFa_ext, %(BYTES)s);
            vFa_ext = %(VINSERT)s(vFa_ext, NEG_INF, 0);
            vFa = %(VSHIFT)s(vFa, %(BYTES)s);
            vFa = %(VINSERT)s(vFa, tmp2, 0);
            for (i=0; i<segLen; ++i) {
                vH = %(VLOAD)s(pvHStore + i);
                vH = %(VMAX)s(vH,vF);
                %(VSTORE)s(pvHStore + i, vH);
                %(SATURATION_CHECK_MID)s
                {
                    %(VTYPE)s vTAll;
                    %(VTYPE)s vT;
                    %(VTYPE)s case1;
                    %(VTYPE)s case2;
                    %(VTYPE)s cond;
                    vHp = %(VADD)s(vHp, %(VLOAD)s(vP + i));
                    case1 = %(VCMPEQ)s(vH, vHp);
                    case2 = %(VCMPEQ)s(vH, vF);
                    cond = %(VANDNOT)s(case1,case2);
                    vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    vT = %(VLOAD)s(pvHT + i);
                    vT = %(VBLEND)s(vT, vTDel, cond);
                    %(VSTORE)s(pvHT + i, vT);
                    vTAll = %(VAND)s(vTAll, vTMask);
                    vTAll = %(VOR)s(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                /* Update vF value. */
                {
                    %(VTYPE)s vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    %(VTYPE)s cond = %(VCMPGT)s(vEF_opn, vFa_ext);
                    %(VTYPE)s vT = %(VBLEND)s(vTDelF, vTDiagF, cond);
                    vTAll = %(VAND)s(vTAll, vFTMask);
                    vTAll = %(VOR)s(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                vEF_opn = %(VSUB)s(vH, vGapO);
                vF_ext = %(VSUB)s(vF, vGapE);
                {
                    %(VTYPE)s vEa = %(VLOAD)s(pvEaLoad + i);
                    %(VTYPE)s vEa_ext = %(VSUB)s(vEa, vGapE);
                    vEa = %(VMAX)s(vEF_opn, vEa_ext);
                    %(VSTORE)s(pvEaStore + i, vEa);
                    if (j+1<s2Len) {
                        %(VTYPE)s cond = %(VCMPGT)s(vEF_opn, vEa_ext);
                        %(VTYPE)s vT = %(VBLEND)s(vTInsE, vTDiagE, cond);
                        arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                    }
                }
                if (! %(VMOVEMASK)s(
                            %(VOR)s(
                                %(VCMPGT)s(vF_ext, vEF_opn),
                                %(VCMPEQ)s(vF_ext, vEF_opn))))
                    goto end;
                /*vF = %(VMAX)s(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = %(VSUB)s(vFa, vGapE);
                vFa = %(VMAX)s(vEF_opn, vFa_ext);
                vHp = %(VLOAD)s(pvHLoad + i);
            }
        }
end:
        {
        }
    }

    /* extract last value from the last column */
    {
        %(VTYPE)s vH = %(VLOAD)s(pvHStore + offset);
        for (k=0; k<position; ++k) {
            vH = %(VSHIFT)s (vH, %(BYTES)s);
        }
        score = (%(INT)s) %(VEXTRACT)s (vH, %(LAST_POS)s);
    }

    %(SATURATION_CHECK_FINAL)s

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;

    parasail_free(boundary);
    parasail_free(pvHT);
    parasail_free(pvEaLoad);
    parasail_free(pvEaStore);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

