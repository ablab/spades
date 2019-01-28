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

#define FASTSTATS

#define SWAP(A,B) { %(VTYPE)s* tmp = A; A = B; B = tmp; }
#define SWAP3(A,B,C) { %(VTYPE)s* tmp = A; A = B; B = C; C = tmp; }

#define NEG_INF %(NEG_INF)s
%(FIXES)s

#ifdef PARASAIL_TABLE
static inline void arr_store(
        int *array,
        %(VTYPE)s vH,
        %(INDEX)s t,
        %(INDEX)s seglen,
        %(INDEX)s d,
        %(INDEX)s dlen,
        %(INDEX)s bias)
{
%(PRINTER_BIAS)s
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        %(VTYPE)s vH,
        %(INDEX)s t,
        %(INDEX)s seglen,
        %(INDEX)s bias)
{
%(PRINTER_BIAS_ROWCOL)s
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME %(NAME_TABLE)s
#define PNAME %(PNAME_TABLE)s
#define INAME PNAME
#define STATIC
#else
#ifdef PARASAIL_ROWCOL
#define FNAME %(NAME_ROWCOL)s
#define PNAME %(PNAME_ROWCOL)s
#define INAME PNAME
#define STATIC
#else
#define FNAME %(NAME)s
#ifdef FASTSTATS
#define PNAME %(PNAME)s_internal
#define INAME %(PNAME)s
#define STATIC static
#else
#define PNAME %(PNAME)s
#define INAME PNAME
#define STATIC
#endif
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_%(ISA)s_%(BITS)s_%(WIDTH)s(s1, s1Len, matrix);
    parasail_result_t *result = INAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

STATIC parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s k = 0;
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const %(INDEX)s segWidth = %(LANES)s; /* number of values in vector unit */
    const %(INDEX)s segLen = (s1Len + segWidth - 1) / segWidth;
    %(VTYPE)s* const restrict vProfile  = (%(VTYPE)s*)profile->profile%(WIDTH)s.score;
    %(VTYPE)s* const restrict vProfileM = (%(VTYPE)s*)profile->profile%(WIDTH)s.matches;
    %(VTYPE)s* const restrict vProfileS = (%(VTYPE)s*)profile->profile%(WIDTH)s.similar;
    %(VTYPE)s* restrict pvHStore        = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLoad         = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHMStore       = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHMLoad        = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHSStore       = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHSLoad        = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLStore       = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLLoad        = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvE       = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvEM      = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvES      = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvEL      = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHMax          = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHMMax         = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHSMax         = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLMax         = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s vGapO = %(VSET1)s(open);
    %(VTYPE)s vGapE = %(VSET1)s(gap);
    %(VTYPE)s vZero = %(VSET0)s();
    %(VTYPE)s vOne = %(VSET1)s(1);
    %(INT)s bias = INT%(WIDTH)s_MIN;
    %(INT)s score = bias;
    %(INT)s matches = bias;
    %(INT)s similar = bias;
    %(INT)s length = bias;
    %(VTYPE)s vBias = %(VSET1)s(bias);
    %(VTYPE)s vBias1 = %(VADD)s(vBias,vOne);
    %(VTYPE)s vMaxH = vBias;
    %(VTYPE)s vMaxHUnit = vBias;
    %(VTYPE)s insert_mask = %(VCMPGT)s(
            %(VSET)s(%(STRIPED_INSERT_MASK)s),
            vZero);
    %(VTYPE)s vSaturationCheckMax = vBias;
    %(VTYPE)s vPosLimit = %(VSET1)s(INT%(WIDTH)s_MAX);
    %(INT)s maxp = INT%(WIDTH)s_MAX - (%(INT)s)(matrix->max+1);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
    const %(INDEX)s offset = (s1Len - 1) %% segLen;
    const %(INDEX)s position = (segWidth - 1) - (s1Len - 1) / segLen;
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    parasail_memset_%(VTYPE)s(pvHStore, vBias, segLen);
    parasail_memset_%(VTYPE)s(pvHMStore, vBias, segLen);
    parasail_memset_%(VTYPE)s(pvHSStore, vBias, segLen);
    parasail_memset_%(VTYPE)s(pvHLStore, vBias, segLen);
    parasail_memset_%(VTYPE)s(pvE, vBias, segLen);
    parasail_memset_%(VTYPE)s(pvEM, vBias, segLen);
    parasail_memset_%(VTYPE)s(pvES, vBias, segLen);
    parasail_memset_%(VTYPE)s(pvEL, vBias1, segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        %(VTYPE)s vEF_opn;
        %(VTYPE)s vE;
        %(VTYPE)s vE_ext;
        %(VTYPE)s vEM;
        %(VTYPE)s vES;
        %(VTYPE)s vEL;
        %(VTYPE)s vF;
        %(VTYPE)s vF_ext;
        %(VTYPE)s vFM;
        %(VTYPE)s vFS;
        %(VTYPE)s vFL;
        %(VTYPE)s vH;
        %(VTYPE)s vH_dag;
        %(VTYPE)s vHM;
        %(VTYPE)s vHS;
        %(VTYPE)s vHL;
        const %(VTYPE)s* vP = NULL;
        const %(VTYPE)s* vPM = NULL;
        const %(VTYPE)s* vPS = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = vBias;
        vFM = vBias;
        vFS = vBias;
        vFL = vBias1;

        /* load final segment of pvHStore and shift left by %(BYTES)s bytes */
        vH = %(VLOAD)s(&pvHStore[segLen - 1]);
        vHM = %(VLOAD)s(&pvHMStore[segLen - 1]);
        vHS = %(VLOAD)s(&pvHSStore[segLen - 1]);
        vHL = %(VLOAD)s(&pvHLStore[segLen - 1]);
        vH = %(VSHIFT)s(vH, %(BYTES)s);
        vHM = %(VSHIFT)s(vHM, %(BYTES)s);
        vHS = %(VSHIFT)s(vHS, %(BYTES)s);
        vHL = %(VSHIFT)s(vHL, %(BYTES)s);
        vH = %(VBLEND)s(vH, vBias, insert_mask);
        vHM = %(VBLEND)s(vHM, vBias, insert_mask);
        vHS = %(VBLEND)s(vHS, vBias, insert_mask);
        vHL = %(VBLEND)s(vHL, vBias, insert_mask);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPM = vProfileM + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + matrix->mapper[(unsigned char)s2[j]] * segLen;

        if (end_ref == j-2) {
            /* Swap in the max buffer. */
            SWAP3(pvHMax,  pvHLoad,  pvHStore)
            SWAP3(pvHMMax, pvHMLoad, pvHMStore)
            SWAP3(pvHSMax, pvHSLoad, pvHSStore)
            SWAP3(pvHLMax, pvHLLoad, pvHLStore)
        }
        else {
            /* Swap the 2 H buffers. */
            SWAP(pvHLoad,  pvHStore)
            SWAP(pvHMLoad, pvHMStore)
            SWAP(pvHSLoad, pvHSStore)
            SWAP(pvHLLoad, pvHLStore)
        }

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            %(VTYPE)s cond_zero;
            %(VTYPE)s case1;
            %(VTYPE)s case2;

            vE = %(VLOAD)s(pvE+ i);
            vEM = %(VLOAD)s(pvEM+ i);
            vES = %(VLOAD)s(pvES+ i);
            vEL = %(VLOAD)s(pvEL+ i);

            /* Get max from vH, vE and vF. */
            vH_dag = %(VADD)s(vH, %(VLOAD)s(vP + i));
            /*vH = %(VMAX)s(vH_dag, vBias);*/
            vH = %(VMAX)s(vH_dag, vE);
            vH = %(VMAX)s(vH, vF);
            /* Save vH values. */
            %(VSTORE)s(pvHStore + i, vH);
            cond_zero = %(VCMPEQ)s(vH, vBias);

            case1 = %(VCMPEQ)s(vH, vH_dag);
            case2 = %(VCMPEQ)s(vH, vF);

            /* calculate vM */
            vHM = %(VBLEND)s(
                    %(VBLEND)s(vEM, vFM, case2),
                    %(VADD)s(vHM, %(VLOAD)s(vPM + i)), case1);
            vHM = %(VBLEND)s(vHM, vBias, cond_zero);
            %(VSTORE)s(pvHMStore + i, vHM);

            /* calculate vS */
            vHS = %(VBLEND)s(
                    %(VBLEND)s(vES, vFS, case2),
                    %(VADD)s(vHS, %(VLOAD)s(vPS + i)), case1);
            vHS = %(VBLEND)s(vHS, vBias, cond_zero);
            %(VSTORE)s(pvHSStore + i, vHS);

            /* calculate vL */
            vHL = %(VBLEND)s(
                    %(VBLEND)s(vEL, vFL, case2),
                    %(VADD)s(vHL, vOne), case1);
            vHL = %(VBLEND)s(vHL, vBias, cond_zero);
            %(VSTORE)s(pvHLStore + i, vHL);

            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHM);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHS);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len, bias);
            arr_store(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len, bias);
            arr_store(result->stats->tables->length_table, vHL, i, segLen, j, s2Len, bias);
            arr_store(result->stats->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
            vMaxH = %(VMAX)s(vH, vMaxH);
            vEF_opn = %(VSUB)s(vH, vGapO);

            /* Update vE value. */
            vE_ext = %(VSUB)s(vE, vGapE);
            vE = %(VMAX)s(vEF_opn, vE_ext);
            case1 = %(VCMPGT)s(vEF_opn, vE_ext);
            vEM = %(VBLEND)s(vEM, vHM, case1);
            vES = %(VBLEND)s(vES, vHS, case1);
            vEL = %(VBLEND)s(
                    %(VADD)s(vEL, vOne),
                    %(VADD)s(vHL, vOne),
                    case1);
            %(VSTORE)s(pvE + i, vE);
            %(VSTORE)s(pvEM + i, vEM);
            %(VSTORE)s(pvES + i, vES);
            %(VSTORE)s(pvEL + i, vEL);

            /* Update vF value. */
            vF_ext = %(VSUB)s(vF, vGapE);
            vF = %(VMAX)s(vEF_opn, vF_ext);
            case1 = %(VCMPGT)s(vEF_opn, vF_ext);
            vFM = %(VBLEND)s(vFM, vHM, case1);
            vFS = %(VBLEND)s(vFS, vHS, case1);
            vFL = %(VBLEND)s(
                    %(VADD)s(vFL, vOne),
                    %(VADD)s(vHL, vOne),
                    case1);

            /* Load the next vH. */
            vH = %(VLOAD)s(pvHLoad + i);
            vHM = %(VLOAD)s(pvHMLoad + i);
            vHS = %(VLOAD)s(pvHSLoad + i);
            vHL = %(VLOAD)s(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            %(VTYPE)s vHp = %(VLOAD)s(&pvHLoad[segLen - 1]);
            vHp = %(VSHIFT)s(vHp, %(BYTES)s);
            vF = %(VSHIFT)s(vF, %(BYTES)s);
            vFM = %(VSHIFT)s(vFM, %(BYTES)s);
            vFS = %(VSHIFT)s(vFS, %(BYTES)s);
            vFL = %(VSHIFT)s(vFL, %(BYTES)s);
            vHp = %(VBLEND)s(vHp, vBias, insert_mask);
            vF = %(VBLEND)s(vF, vBias, insert_mask);
            vFM = %(VBLEND)s(vFM, vBias, insert_mask);
            vFS = %(VBLEND)s(vFS, vBias, insert_mask);
            vFL = %(VBLEND)s(vFL, vBias1, insert_mask);
            for (i=0; i<segLen; ++i) {
                %(VTYPE)s case1;
                %(VTYPE)s case2;
                %(VTYPE)s cond;

                vHp = %(VADD)s(vHp, %(VLOAD)s(vP + i));
                vH = %(VLOAD)s(pvHStore + i);
                vH = %(VMAX)s(vH, vF);
                %(VSTORE)s(pvHStore + i, vH);
                case1 = %(VCMPEQ)s(vH, vHp);
                case2 = %(VCMPEQ)s(vH, vF);
                cond = %(VANDNOT)s(case1, case2);

                /* calculate vM */
                vHM = %(VLOAD)s(pvHMStore + i);
                vHM = %(VBLEND)s(vHM, vFM, cond);
                %(VSTORE)s(pvHMStore + i, vHM);

                /* calculate vS */
                vHS = %(VLOAD)s(pvHSStore + i);
                vHS = %(VBLEND)s(vHS, vFS, cond);
                %(VSTORE)s(pvHSStore + i, vHS);

                /* calculate vL */
                vHL = %(VLOAD)s(pvHLStore + i);
                vHL = %(VBLEND)s(vHL, vFL, cond);
                %(VSTORE)s(pvHLStore + i, vHL);

                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHM);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHS);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len, bias);
                arr_store(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len, bias);
                arr_store(result->stats->tables->length_table, vHL, i, segLen, j, s2Len, bias);
                arr_store(result->stats->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
                vMaxH = %(VMAX)s(vH, vMaxH);
                /* Update vF value. */
                vEF_opn = %(VSUB)s(vH, vGapO);
                vF_ext = %(VSUB)s(vF, vGapE);
                if (! %(VMOVEMASK)s(
                            %(VOR)s(
                                %(VCMPGT)s(vF_ext, vEF_opn),
                                %(VAND)s(
                                    %(VCMPEQ)s(vF_ext, vEF_opn),
                                    %(VCMPGT)s(vF_ext, vBias)))))
                    goto end;
                /*vF = %(VMAX)s(vEF_opn, vF_ext);*/
                vF = vF_ext;
                cond = %(VCMPGT)s(vEF_opn, vF_ext);
                vFM = %(VBLEND)s(vFM, vHM, cond);
                vFS = %(VBLEND)s(vFS, vHS, cond);
                vFL = %(VBLEND)s(
                        %(VADD)s(vFL, vOne),
                        %(VADD)s(vHL, vOne),
                        cond);
                vHp = %(VLOAD)s(pvHLoad + i);
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = %(VLOAD)s(pvHStore + offset);
            vHM = %(VLOAD)s(pvHMStore + offset);
            vHS = %(VLOAD)s(pvHSStore + offset);
            vHL = %(VLOAD)s(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = %(VSHIFT)s(vH, %(BYTES)s);
                vHM = %(VSHIFT)s(vHM, %(BYTES)s);
                vHS = %(VSHIFT)s(vHS, %(BYTES)s);
                vHL = %(VSHIFT)s(vHL, %(BYTES)s);
            }
            result->stats->rowcols->score_row[j] = (%(INT)s) %(VEXTRACT)s (vH, %(LAST_POS)s) - bias;
            result->stats->rowcols->matches_row[j] = (%(INT)s) %(VEXTRACT)s (vHM, %(LAST_POS)s) - bias;
            result->stats->rowcols->similar_row[j] = (%(INT)s) %(VEXTRACT)s (vHS, %(LAST_POS)s) - bias;
            result->stats->rowcols->length_row[j] = (%(INT)s) %(VEXTRACT)s (vHL, %(LAST_POS)s) - bias;
        }
#endif

        {
            %(VTYPE)s vCompare = %(VCMPGT)s(vMaxH, vMaxHUnit);
            if (%(VMOVEMASK)s(vCompare)) {
                score = %(VHMAX)s(vMaxH);
                /* if score has potential to overflow, abort early */
                if (score > maxp) {
                    result->flag |= PARASAIL_FLAG_SATURATED;
                    break;
                }
                vMaxHUnit = %(VSET1)s(score);
                end_ref = j;
            }
        }
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        %(VTYPE)s vH = %(VLOAD)s(pvHStore+i);
        %(VTYPE)s vHM = %(VLOAD)s(pvHMStore+i);
        %(VTYPE)s vHS = %(VLOAD)s(pvHSStore+i);
        %(VTYPE)s vHL = %(VLOAD)s(pvHLStore+i);
        arr_store_col(result->stats->rowcols->score_col, vH, i, segLen, bias);
        arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen, bias);
        arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen, bias);
        arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen, bias);
    }
#endif

    if (score == INT%(WIDTH)s_MAX
            || %(VMOVEMASK)s(%(VCMPEQ)s(vSaturationCheckMax,vPosLimit))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
    }

    if (parasail_result_is_saturated(result)) {
        score = 0;
        end_query = 0;
        end_ref = 0;
        matches = 0;
        similar = 0;
        length = 0;
    }
    else {
        if (end_ref == j-1) {
            /* end_ref was the last store column */
            SWAP(pvHMax,  pvHStore)
            SWAP(pvHMMax, pvHMStore)
            SWAP(pvHSMax, pvHSStore)
            SWAP(pvHLMax, pvHLStore)
        }
        else if (end_ref == j-2) {
            /* end_ref was the last load column */
            SWAP(pvHMax,  pvHLoad)
            SWAP(pvHMMax, pvHMLoad)
            SWAP(pvHSMax, pvHSLoad)
            SWAP(pvHLMax, pvHLLoad)
        }
        /* Trace the alignment ending position on read. */
        {
            %(INT)s *t = (%(INT)s*)pvHMax;
            %(INT)s *m = (%(INT)s*)pvHMMax;
            %(INT)s *s = (%(INT)s*)pvHSMax;
            %(INT)s *l = (%(INT)s*)pvHLMax;
            %(INDEX)s column_len = segLen * segWidth;
            end_query = s1Len;
            for (i = 0; i<column_len; ++i, ++t, ++m, ++s, ++l) {
                if (*t == score) {
                    %(INDEX)s temp = i / segWidth + i %% segWidth * segLen;
                    if (temp < end_query) {
                        end_query = temp;
                        matches = *m;
                        similar = *s;
                        length = *l;
                    }
                }
            }
        }
    }

    result->score = score - bias;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->stats->matches = matches - bias;
    result->stats->similar = similar - bias;
    result->stats->length = length - bias;
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(pvHLMax);
    parasail_free(pvHSMax);
    parasail_free(pvHMMax);
    parasail_free(pvHMax);
    parasail_free(pvEL);
    parasail_free(pvES);
    parasail_free(pvEM);
    parasail_free(pvE);
    parasail_free(pvHLLoad);
    parasail_free(pvHLStore);
    parasail_free(pvHSLoad);
    parasail_free(pvHSStore);
    parasail_free(pvHMLoad);
    parasail_free(pvHMStore);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

#ifdef FASTSTATS
#ifdef PARASAIL_TABLE
#else
#ifdef PARASAIL_ROWCOL
#else
#include <assert.h>
parasail_result_t* INAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    const char *s1 = profile->s1;
    const parasail_matrix_t *matrix = profile->matrix;

    /* find the end loc first with the faster implementation */
    parasail_result_t *result = %(PNAME_BASE)s(profile, s2, s2Len, open, gap);
    if (!parasail_result_is_saturated(result)) {
#if 0
        int s1Len_new = 0;
        int s2Len_new = 0;
        char *s1_new = NULL;
        char *s2_new = NULL;
        parasail_profile_t *profile_new = NULL;
        parasail_result_t *result_new = NULL;
        int s1_begin = 0;
        int s2_begin = 0;
        int s1Len_final = 0;
        int s2Len_final = 0;
        parasail_profile_t *profile_final = NULL;
        parasail_result_t *result_final = NULL;

        /* using the end loc and the non-stats version of the function,
         * reverse the inputs and find the beg loc */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;
        s1_new = parasail_reverse(s1, s1Len_new);
        s2_new = parasail_reverse(s2, s2Len_new);
        profile_new = parasail_profile_create_%(ISA)s_%(BITS)s_%(WIDTH)s(
                s1_new, s1Len_new, matrix);
        profile_new->stop = result->score;
        result_new = %(PNAME_BASE)s(
                profile_new, s2_new, s2Len_new, open, gap);

        /* using both the beg and end loc, call the original stats func */
        s1_begin = s1Len_new - result_new->end_query - 1;
        s2_begin = s2Len_new - result_new->end_ref - 1;
        s1Len_final = s1Len_new - s1_begin;
        s2Len_final = s2Len_new - s2_begin;
        assert(s1_begin >= 0);
        assert(s2_begin >= 0);
        assert(s1Len_new > s1_begin);
        assert(s2Len_new > s2_begin);
        profile_final = parasail_profile_create_stats_%(ISA)s_%(BITS)s_%(WIDTH)s(
                &s1[s1_begin], s1Len_final, matrix);
        result_final = PNAME(
                profile_final, &s2[s2_begin], s2Len_final, open, gap);

        /* clean up all the temporary profiles, sequences, and results */
        free(s1_new);
        free(s2_new);
        parasail_profile_free(profile_new);
        parasail_profile_free(profile_final);
        parasail_result_free(result);
        parasail_result_free(result_new);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
#else
        int s1Len_new = 0;
        int s2Len_new = 0;
        parasail_profile_t *profile_final = NULL;
        parasail_result_t *result_final = NULL;

        /* using the end loc, call the original stats function */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;
        profile_final = parasail_profile_create_stats_%(ISA)s_%(BITS)s_%(WIDTH)s(
                s1, s1Len_new, matrix);
        result_final = PNAME(
                profile_final, s2, s2Len_new, open, gap);

        /* clean up all the temporary profiles, sequences, and results */
        parasail_profile_free(profile_final);
        parasail_result_free(result);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
#endif
    }
    else {
        return result;
    }
}
#endif
#endif
#endif

