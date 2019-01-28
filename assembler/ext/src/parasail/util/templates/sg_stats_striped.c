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

%(FIXES)s

#ifdef PARASAIL_TABLE
static inline void arr_store_si%(BITS)s(
        int *array,
        %(VTYPE)s vH,
        %(INDEX)s t,
        %(INDEX)s seglen,
        %(INDEX)s d,
        %(INDEX)s dlen)
{
%(PRINTER)s
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        %(VTYPE)s vH,
        %(INDEX)s t,
        %(INDEX)s seglen)
{
%(PRINTER_ROWCOL)s
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
    const %(INDEX)s offset = (s1Len - 1) %% segLen;
    const %(INDEX)s position = (segWidth - 1) - (s1Len - 1) / segLen;
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
    const %(VTYPE)s vGapO = %(VSET1)s(open);
    const %(VTYPE)s vGapE = %(VSET1)s(gap);
    const %(INT)s NEG_LIMIT = (-open < matrix->min ?
        INT%(WIDTH)s_MIN + open : INT%(WIDTH)s_MIN - matrix->min) + 1;
    const %(INT)s POS_LIMIT = INT%(WIDTH)s_MAX - matrix->max - 1;
    const %(VTYPE)s vZero = %(VSET0)s();
    const %(VTYPE)s vOne = %(VSET1)s(1);
    %(INT)s score = NEG_LIMIT;
    %(INT)s matches = 0;
    %(INT)s similar = 0;
    %(INT)s length = 0;
    %(VTYPE)s vNegLimit = %(VSET1)s(NEG_LIMIT);
    %(VTYPE)s vPosLimit = %(VSET1)s(POS_LIMIT);
    %(VTYPE)s vSaturationCheckMin = vPosLimit;
    %(VTYPE)s vSaturationCheckMax = vNegLimit;
    %(VTYPE)s vMaxH = vNegLimit;
    %(VTYPE)s vMaxHM = vNegLimit;
    %(VTYPE)s vMaxHS = vNegLimit;
    %(VTYPE)s vMaxHL = vNegLimit;
    %(VTYPE)s vPosMask = %(VCMPEQ)s(%(VSET1)s(position),
            %(VSET)s(%(POSITION_MASK)s));
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    parasail_memset_%(VTYPE)s(pvHStore, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvHMStore, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvHSStore, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvHLStore, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvE, %(VSET1)s(-open), segLen);
    parasail_memset_%(VTYPE)s(pvEM, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvES, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvEL, vOne, segLen);

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

        /* Initialize F value to neg inf.  Any errors to vH values will
         * be corrected in the Lazy_F loop. */
        vF = vNegLimit;
        vFM = vZero;
        vFS = vZero;
        vFL = vOne;

        /* load final segment of pvHStore and shift left by %(BYTES)s bytes */
        vH = %(VLOAD)s(&pvHStore[segLen - 1]);
        vHM = %(VLOAD)s(&pvHMStore[segLen - 1]);
        vHS = %(VLOAD)s(&pvHSStore[segLen - 1]);
        vHL = %(VLOAD)s(&pvHLStore[segLen - 1]);
        vH = %(VSHIFT)s(vH, %(BYTES)s);
        vHM = %(VSHIFT)s(vHM, %(BYTES)s);
        vHS = %(VSHIFT)s(vHS, %(BYTES)s);
        vHL = %(VSHIFT)s(vHL, %(BYTES)s);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPM = vProfileM + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad,  pvHStore)
        SWAP(pvHMLoad, pvHMStore)
        SWAP(pvHSLoad, pvHSStore)
        SWAP(pvHLLoad, pvHLStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            %(VTYPE)s case1;
            %(VTYPE)s case2;

            vE = %(VLOAD)s(pvE+ i);
            vEM = %(VLOAD)s(pvEM+ i);
            vES = %(VLOAD)s(pvES+ i);
            vEL = %(VLOAD)s(pvEL+ i);

            /* Get max from vH, vE and vF. */
            vH_dag = %(VADD)s(vH, %(VLOAD)s(vP + i));
            vH = %(VMAX)s(vH_dag, vE);
            vH = %(VMAX)s(vH, vF);
            /* Save vH values. */
            %(VSTORE)s(pvHStore + i, vH);

            case1 = %(VCMPEQ)s(vH, vH_dag);
            case2 = %(VCMPEQ)s(vH, vF);

            /* calculate vM */
            vHM = %(VBLEND)s(
                    %(VBLEND)s(vEM, vFM, case2),
                    %(VADD)s(vHM, %(VLOAD)s(vPM + i)),
                    case1);
            %(VSTORE)s(pvHMStore + i, vHM);

            /* calculate vS */
            vHS = %(VBLEND)s(
                    %(VBLEND)s(vES, vFS, case2),
                    %(VADD)s(vHS, %(VLOAD)s(vPS + i)),
                    case1);
            %(VSTORE)s(pvHSStore + i, vHS);

            /* calculate vL */
            vHL = %(VBLEND)s(
                    %(VBLEND)s(vEL, vFL, case2),
                    %(VADD)s(vHL, vOne),
                    case1);
            %(VSTORE)s(pvHLStore + i, vHL);

            vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHM);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHS);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si%(BITS)s(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si%(BITS)s(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si%(BITS)s(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
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
            vF = %(VINSERT)s(vF, -open, 0);
            vFM = %(VSHIFT)s(vFM, %(BYTES)s);
            vFS = %(VSHIFT)s(vFS, %(BYTES)s);
            vFL = %(VSHIFT)s(vFL, %(BYTES)s);
            vFL = %(VINSERT)s(vFL, 1, 0);
            for (i=0; i<segLen; ++i) {
                %(VTYPE)s case1;
                %(VTYPE)s case2;
                %(VTYPE)s cond;

                vHp = %(VADD)s(vHp, %(VLOAD)s(vP + i));
                vH = %(VLOAD)s(pvHStore + i);
                vH = %(VMAX)s(vH,vF);
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

                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHM);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHS);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store_si%(BITS)s(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si%(BITS)s(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si%(BITS)s(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si%(BITS)s(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                /* Update vF value. */
                vEF_opn = %(VSUB)s(vH, vGapO);
                vF_ext = %(VSUB)s(vF, vGapE);
                if (! %(VMOVEMASK)s(
                            %(VOR)s(
                                %(VCMPGT)s(vF_ext, vEF_opn),
                                %(VCMPEQ)s(vF_ext, vEF_opn))))
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

        /* extract vector containing last value from the column */
        {
            %(VTYPE)s cond_max;
            vH = %(VLOAD)s(pvHStore + offset);
            vHM = %(VLOAD)s(pvHMStore + offset);
            vHS = %(VLOAD)s(pvHSStore + offset);
            vHL = %(VLOAD)s(pvHLStore + offset);
            cond_max = %(VCMPGT)s(vH, vMaxH);
            vMaxH = %(VBLEND)s(vMaxH, vH, cond_max);
            vMaxHM = %(VBLEND)s(vMaxHM, vHM, cond_max);
            vMaxHS = %(VBLEND)s(vMaxHS, vHS, cond_max);
            vMaxHL = %(VBLEND)s(vMaxHL, vHL, cond_max);
            if (%(VMOVEMASK)s(%(VAND)s(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
            }
        }
#ifdef PARASAIL_ROWCOL
        for (k=0; k<position; ++k) {
            vH = %(VSHIFT)s(vH, %(BYTES)s);
            vHM = %(VSHIFT)s(vHM, %(BYTES)s);
            vHS = %(VSHIFT)s(vHS, %(BYTES)s);
            vHL = %(VSHIFT)s(vHL, %(BYTES)s);
        }
        result->stats->rowcols->score_row[j] = (%(INT)s) %(VEXTRACT)s (vH, %(LAST_POS)s);
        result->stats->rowcols->matches_row[j] = (%(INT)s) %(VEXTRACT)s (vHM, %(LAST_POS)s);
        result->stats->rowcols->similar_row[j] = (%(INT)s) %(VEXTRACT)s (vHS, %(LAST_POS)s);
        result->stats->rowcols->length_row[j] = (%(INT)s) %(VEXTRACT)s (vHL, %(LAST_POS)s);
#endif
    }

    {
        /* extract last value from the column maximums */
        for (k=0; k<position; ++k) {
            vMaxH  = %(VSHIFT)s (vMaxH, %(BYTES)s);
            vMaxHM = %(VSHIFT)s (vMaxHM, %(BYTES)s);
            vMaxHS = %(VSHIFT)s (vMaxHS, %(BYTES)s);
            vMaxHL = %(VSHIFT)s (vMaxHL, %(BYTES)s);
        }
        score = (%(INT)s) %(VEXTRACT)s (vMaxH, %(LAST_POS)s);
        matches = (%(INT)s)%(VEXTRACT)s(vMaxHM, %(LAST_POS)s);
        similar = (%(INT)s)%(VEXTRACT)s(vMaxHS, %(LAST_POS)s);
        length = (%(INT)s)%(VEXTRACT)s(vMaxHL, %(LAST_POS)s);
    }

    /* max of last column */
    if (INT32_MAX == profile->stop || 0 == profile->stop)
    {
        %(INT)s score_last;
        vMaxH = vNegLimit;

        if (0 == profile->stop) {
            /* ignore last row contributions */
            score = NEG_LIMIT;
            matches = 0;
            similar = 0;
            length = 0;
            end_query = s1Len;
            end_ref = s2Len - 1;
        }

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            %(VTYPE)s vH = %(VLOAD)s(pvHStore + i);
#ifdef PARASAIL_ROWCOL
            %(VTYPE)s vHM = %(VLOAD)s(pvHMStore + i);
            %(VTYPE)s vHS = %(VLOAD)s(pvHSStore + i);
            %(VTYPE)s vHL = %(VLOAD)s(pvHLStore + i);
            arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
            arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
            arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
            arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
#endif
            vMaxH = %(VMAX)s(vH, vMaxH);
        }

        /* max in vec */
        score_last = %(VHMAX)s(vMaxH);
        if (score_last > score || (score_last == score && end_ref == s2Len - 1)) {
            end_query = s1Len;
            end_ref = s2Len - 1;
            /* Trace the alignment ending position on read. */
            {
                %(INT)s *t = (%(INT)s*)pvHStore;
                %(INT)s *m = (%(INT)s*)pvHMStore;
                %(INT)s *s = (%(INT)s*)pvHSStore;
                %(INT)s *l = (%(INT)s*)pvHLStore;
                %(INDEX)s column_len = segLen * segWidth;
                for (i = 0; i<column_len; ++i, ++t, ++m, ++s, ++l) {
                    %(INDEX)s temp = i / segWidth + i %% segWidth * segLen;
                    if (temp < s1Len) {
                        if (*t > score || (*t == score && temp < end_query)) {
                            score = *t;
                            end_query = temp;
                            matches = *m;
                            similar = *s;
                            length = *l;
                        }
                    }
                }
            }
        }
    }

    if (%(VMOVEMASK)s(%(VOR)s(
            %(VCMPLT)s(vSaturationCheckMin, vNegLimit),
            %(VCMPGT)s(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        matches = 0;
        similar = 0;
        length = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->stats->matches = matches;
    result->stats->similar = similar;
    result->stats->length = length;
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

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
        int s1Len_new = 0;
        int s2Len_new = 0;
        parasail_result_t *result_final = NULL;

        /* using the end loc, call the original stats function */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;

        if (s1Len_new == profile->s1Len) {
            /* special 'stop' value tells stats function not to
             * consider last column results */
            int stop_save = profile->stop;
            ((parasail_profile_t*)profile)->stop = 1;
            result_final = PNAME(
                    profile, s2, s2Len_new, open, gap);
            ((parasail_profile_t*)profile)->stop = stop_save;
        }
        else {
            parasail_profile_t *profile_final = NULL;
            profile_final = parasail_profile_create_stats_%(ISA)s_%(BITS)s_%(WIDTH)s(
                    s1, s1Len_new, matrix);
            /* special 'stop' value tells stats function not to
             * consider last row results */
            profile_final->stop = 0;
            result_final = PNAME(
                    profile_final, s2, s2Len_new, open, gap);

            parasail_profile_free(profile_final);
        }

        parasail_result_free(result);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
    }
    else {
        return result;
    }
}
#endif
#endif
#endif

