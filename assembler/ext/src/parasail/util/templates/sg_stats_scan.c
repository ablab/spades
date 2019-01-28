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
#else
#ifdef PARASAIL_ROWCOL
#define FNAME %(NAME_ROWCOL)s
#define PNAME %(PNAME_ROWCOL)s
#else
#define FNAME %(NAME)s
#define PNAME %(PNAME)s
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_%(ISA)s_%(BITS)s_%(WIDTH)s(s1, s1Len, matrix);
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
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const %(INDEX)s segWidth = %(LANES)s; /* number of values in vector unit */
    const %(INDEX)s segLen = (s1Len + segWidth - 1) / segWidth;
    const %(INDEX)s offset = (s1Len - 1) %% segLen;
    const %(INDEX)s position = (segWidth - 1) - (s1Len - 1) / segLen;
    %(VTYPE)s* const restrict pvP  = (%(VTYPE)s*)profile->profile%(WIDTH)s.score;
    %(VTYPE)s* const restrict pvPm = (%(VTYPE)s*)profile->profile%(WIDTH)s.matches;
    %(VTYPE)s* const restrict pvPs = (%(VTYPE)s*)profile->profile%(WIDTH)s.similar;
    %(VTYPE)s* const restrict pvE  = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvEM = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvES = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvEL = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvH  = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHM = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHS = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHL = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHMax  = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHMMax = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHSMax = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHLMax = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvGapper = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvGapperL = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s vGapO = %(VSET1)s(open);
    %(VTYPE)s vGapE = %(VSET1)s(gap);
    const %(INT)s NEG_LIMIT = (-open < matrix->min ?
        INT%(WIDTH)s_MIN + open : INT%(WIDTH)s_MIN - matrix->min) + 1;
    const %(INT)s POS_LIMIT = INT%(WIDTH)s_MAX - matrix->max - 1;
    %(VTYPE)s vZero = %(VSET0)s();
    %(VTYPE)s vOne = %(VSET1)s(1);
    %(INT)s score = NEG_LIMIT;
    %(INT)s matches = 0;
    %(INT)s similar = 0;
    %(INT)s length = 0;
    %(VTYPE)s vNegLimit = %(VSET1)s(NEG_LIMIT);
    %(VTYPE)s vPosLimit = %(VSET1)s(POS_LIMIT);
    %(VTYPE)s vSaturationCheckMin = vPosLimit;
    %(VTYPE)s vSaturationCheckMax = vNegLimit;
    %(VTYPE)s vMaxH = vNegLimit;
    %(VTYPE)s vMaxM = vNegLimit;
    %(VTYPE)s vMaxS = vNegLimit;
    %(VTYPE)s vMaxL = vNegLimit;
    %(VTYPE)s vPosMask = %(VCMPEQ)s(%(VSET1)s(position),
            %(VSET)s(%(POSITION_MASK)s));
    %(VTYPE)s vNegInfFront = vZero;
    %(VTYPE)s vSegLenXgap;
    %(VTYPE)s vSegLen = %(VSHIFT)s(%(VSET1)s(segLen), %(BYTES)s);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif

    vNegInfFront = %(VINSERT)s(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = %(VADD)s(vNegInfFront,
            %(VSHIFT)s(%(VSET1)s(-segLen*gap), %(BYTES)s));

    parasail_memset_%(VTYPE)s(pvH, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvHM, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvHS, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvHL, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvE, vNegLimit, segLen);
    parasail_memset_%(VTYPE)s(pvEM, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvES, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvEL, vZero, segLen);
    {
        %(VTYPE)s vGapper = %(VSUB)s(vZero,vGapO);
        %(VTYPE)s vGapperL = vOne;
        for (i=segLen-1; i>=0; --i) {
            %(VSTORE)s(pvGapper+i, vGapper);
            %(VSTORE)s(pvGapperL+i, vGapperL);
            vGapper = %(VSUB)s(vGapper, vGapE);
            vGapperL = %(VADD)s(vGapperL, vOne);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        %(VTYPE)s vE;
        %(VTYPE)s vE_ext;
        %(VTYPE)s vE_opn;
        %(VTYPE)s vEM;
        %(VTYPE)s vES;
        %(VTYPE)s vEL;
        %(VTYPE)s vHt;
        %(VTYPE)s vHtM;
        %(VTYPE)s vHtS;
        %(VTYPE)s vHtL;
        %(VTYPE)s vF;
        %(VTYPE)s vF_ext;
        %(VTYPE)s vF_opn;
        %(VTYPE)s vFM;
        %(VTYPE)s vFS;
        %(VTYPE)s vFL;
        %(VTYPE)s vH;
        %(VTYPE)s vHM;
        %(VTYPE)s vHS;
        %(VTYPE)s vHL;
        %(VTYPE)s vHp;
        %(VTYPE)s vHpM;
        %(VTYPE)s vHpS;
        %(VTYPE)s vHpL;
        %(VTYPE)s *pvW;
        %(VTYPE)s vW;
        %(VTYPE)s *pvWM;
        %(VTYPE)s vWM;
        %(VTYPE)s *pvWS;
        %(VTYPE)s vWS;
        %(VTYPE)s case1;
        %(VTYPE)s case2;
        %(VTYPE)s vGapper;
        %(VTYPE)s vGapperL;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = %(VLOAD)s(pvH+(segLen-1));
        vHpM = %(VLOAD)s(pvHM+(segLen-1));
        vHpS = %(VLOAD)s(pvHS+(segLen-1));
        vHpL = %(VLOAD)s(pvHL+(segLen-1));
        vHp = %(VSHIFT)s(vHp, %(BYTES)s);
        vHpM = %(VSHIFT)s(vHpM, %(BYTES)s);
        vHpS = %(VSHIFT)s(vHpS, %(BYTES)s);
        vHpL = %(VSHIFT)s(vHpL, %(BYTES)s);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWM= pvPm+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvWS= pvPs+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = %(VSUB)s(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        vFM = vZero;
        vFS = vZero;
        vFL = vZero;
        for (i=0; i<segLen; ++i) {
            vH = %(VLOAD)s(pvH+i);
            vHM= %(VLOAD)s(pvHM+i);
            vHS= %(VLOAD)s(pvHS+i);
            vHL= %(VLOAD)s(pvHL+i);
            vE = %(VLOAD)s(pvE+i);
            vEM= %(VLOAD)s(pvEM+i);
            vES= %(VLOAD)s(pvES+i);
            vEL= %(VLOAD)s(pvEL+i);
            vW = %(VLOAD)s(pvW+i);
            vWM = %(VLOAD)s(pvWM+i);
            vWS = %(VLOAD)s(pvWS+i);
            vGapper = %(VLOAD)s(pvGapper+i);
            vGapperL = %(VLOAD)s(pvGapperL+i);
            vE_opn = %(VSUB)s(vH, vGapO);
            vE_ext = %(VSUB)s(vE, vGapE);
            case1 = %(VCMPGT)s(vE_opn, vE_ext);
            vE = %(VMAX)s(vE_opn, vE_ext);
            vEM = %(VBLEND)s(vEM, vHM, case1);
            vES = %(VBLEND)s(vES, vHS, case1);
            vEL = %(VBLEND)s(vEL, vHL, case1);
            vEL = %(VADD)s(vEL, vOne);
            vGapper = %(VADD)s(vHt, vGapper);
            case1 = %(VOR)s(
                    %(VCMPGT)s(vF, vGapper),
                    %(VCMPEQ)s(vF, vGapper));
            vF = %(VMAX)s(vF, vGapper);
            vFM = %(VBLEND)s(vHtM, vFM, case1);
            vFS = %(VBLEND)s(vHtS, vFS, case1);
            vFL = %(VBLEND)s(
                    %(VADD)s(vHtL, vGapperL),
                    vFL, case1);
            vHp = %(VADD)s(vHp, vW);
            vHpM = %(VADD)s(vHpM, vWM);
            vHpS = %(VADD)s(vHpS, vWS);
            vHpL = %(VADD)s(vHpL, vOne);
            case1 = %(VCMPGT)s(vE, vHp);
            vHt = %(VMAX)s(vE, vHp);
            vHtM = %(VBLEND)s(vHpM, vEM, case1);
            vHtS = %(VBLEND)s(vHpS, vES, case1);
            vHtL = %(VBLEND)s(vHpL, vEL, case1);
            %(VSTORE)s(pvE+i, vE);
            %(VSTORE)s(pvEM+i, vEM);
            %(VSTORE)s(pvES+i, vES);
            %(VSTORE)s(pvEL+i, vEL);
            %(VSTORE)s(pvH+i, vHp);
            %(VSTORE)s(pvHM+i, vHpM);
            %(VSTORE)s(pvHS+i, vHpS);
            %(VSTORE)s(pvHL+i, vHpL);
            vHp = vH;
            vHpM = vHM;
            vHpS = vHS;
            vHpL = vHL;
        }

        /* pseudo prefix scan on F and H */
        vHt = %(VSHIFT)s(vHt, %(BYTES)s);
        vHtM = %(VSHIFT)s(vHtM, %(BYTES)s);
        vHtS = %(VSHIFT)s(vHtS, %(BYTES)s);
        vHtL = %(VSHIFT)s(vHtL, %(BYTES)s);
        vGapper = %(VLOAD)s(pvGapper);
        vGapperL = %(VLOAD)s(pvGapperL);
        vGapper = %(VADD)s(vHt, vGapper);
        case1 = %(VOR)s(
                %(VCMPGT)s(vGapper, vF),
                %(VCMPEQ)s(vGapper, vF));
        vF = %(VMAX)s(vF, vGapper);
        vFM = %(VBLEND)s(vFM, vHtM, case1);
        vFS = %(VBLEND)s(vFS, vHtS, case1);
        vFL = %(VBLEND)s(
                vFL,
                %(VADD)s(vHtL, vGapperL),
                case1);
        for (i=0; i<segWidth-2; ++i) {
            %(VTYPE)s vFt = %(VSHIFT)s(vF, %(BYTES)s);
            %(VTYPE)s vFtM = %(VSHIFT)s(vFM, %(BYTES)s);
            %(VTYPE)s vFtS = %(VSHIFT)s(vFS, %(BYTES)s);
            %(VTYPE)s vFtL = %(VSHIFT)s(vFL, %(BYTES)s);
            vFt = %(VADD)s(vFt, vSegLenXgap);
            case1 = %(VOR)s(
                    %(VCMPGT)s(vFt, vF),
                    %(VCMPEQ)s(vFt, vF));
            vF = %(VMAX)s(vF, vFt);
            vFM = %(VBLEND)s(vFM, vFtM, case1);
            vFS = %(VBLEND)s(vFS, vFtS, case1);
            vFL = %(VBLEND)s(
                    vFL,
                    %(VADD)s(vFtL, vSegLen),
                    case1);
        }

        /* calculate final H */
        vF = %(VSHIFT)s(vF, %(BYTES)s);
        vFM = %(VSHIFT)s(vFM, %(BYTES)s);
        vFS = %(VSHIFT)s(vFS, %(BYTES)s);
        vFL = %(VSHIFT)s(vFL, %(BYTES)s);
        vF = %(VADD)s(vF, vNegInfFront);
        case1 = %(VCMPGT)s(vF, vHt);
        vH = %(VMAX)s(vF, vHt);
        vHM = %(VBLEND)s(vHtM, vFM, case1);
        vHS = %(VBLEND)s(vHtS, vFS, case1);
        vHL = %(VBLEND)s(vHtL, vFL, case1);
        for (i=0; i<segLen; ++i) {
            vHp = %(VLOAD)s(pvH+i);
            vHpM = %(VLOAD)s(pvHM+i);
            vHpS = %(VLOAD)s(pvHS+i);
            vHpL = %(VLOAD)s(pvHL+i);
            vE = %(VLOAD)s(pvE+i);
            vEM = %(VLOAD)s(pvEM+i);
            vES = %(VLOAD)s(pvES+i);
            vEL = %(VLOAD)s(pvEL+i);
            vF_opn = %(VSUB)s(vH, vGapO);
            vF_ext = %(VSUB)s(vF, vGapE);
            vF = %(VMAX)s(vF_opn, vF_ext);
            case1 = %(VCMPGT)s(vF_opn, vF_ext);
            vFM = %(VBLEND)s(vFM, vHM, case1);
            vFS = %(VBLEND)s(vFS, vHS, case1);
            vFL = %(VBLEND)s(vFL, vHL, case1);
            vFL = %(VADD)s(vFL, vOne);
            vH = %(VMAX)s(vHp, vE);
            vH = %(VMAX)s(vH, vF);
            case1 = %(VCMPEQ)s(vH, vHp);
            case2 = %(VCMPEQ)s(vH, vF);
            vHM = %(VBLEND)s(
                    %(VBLEND)s(vEM, vFM, case2),
                    vHpM, case1);
            vHS = %(VBLEND)s(
                    %(VBLEND)s(vES, vFS, case2),
                    vHpS, case1);
            vHL = %(VBLEND)s(
                    %(VBLEND)s(vEL, vFL, case2),
                    vHpL, case1);
            %(VSTORE)s(pvH+i, vH);
            %(VSTORE)s(pvHM+i, vHM);
            %(VSTORE)s(pvHS+i, vHS);
            %(VSTORE)s(pvHL+i, vHL);
            vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHM);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHS);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
            arr_store_si%(BITS)s(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si%(BITS)s(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si%(BITS)s(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
#endif
        } 
        /* extract vector containing last value from column */
        {
            %(VTYPE)s cond_max;
            vH = %(VLOAD)s(pvH + offset);
            vHM = %(VLOAD)s(pvHM + offset);
            vHS = %(VLOAD)s(pvHS + offset);
            vHL = %(VLOAD)s(pvHL + offset);
            cond_max = %(VCMPGT)s(vH, vMaxH);
            vMaxH = %(VMAX)s(vMaxH, vH);
            vMaxM = %(VBLEND)s(vMaxM, vHM, cond_max);
            vMaxS = %(VBLEND)s(vMaxS, vHS, cond_max);
            vMaxL = %(VBLEND)s(vMaxL, vHL, cond_max);
            if (%(VMOVEMASK)s(%(VAND)s(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
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
    }

    /* max last value from all columns */
    {
        for (k=0; k<position; ++k) {
            vMaxH = %(VSHIFT)s(vMaxH, %(BYTES)s);
            vMaxM = %(VSHIFT)s(vMaxM, %(BYTES)s);
            vMaxS = %(VSHIFT)s(vMaxS, %(BYTES)s);
            vMaxL = %(VSHIFT)s(vMaxL, %(BYTES)s);
        }
        score = (%(INT)s) %(VEXTRACT)s(vMaxH, %(LAST_POS)s);
        matches = (%(INT)s) %(VEXTRACT)s(vMaxM, %(LAST_POS)s);
        similar = (%(INT)s) %(VEXTRACT)s(vMaxS, %(LAST_POS)s);
        length = (%(INT)s) %(VEXTRACT)s(vMaxL, %(LAST_POS)s);
    }

    /* max of last column */
    {
        %(INT)s score_last;
        vMaxH = vNegLimit;

        for (i=0; i<segLen; ++i) {
            %(VTYPE)s vH = %(VLOAD)s(pvH + i);
#ifdef PARASAIL_ROWCOL
            %(VTYPE)s vHM = %(VLOAD)s(pvHM + i);
            %(VTYPE)s vHS = %(VLOAD)s(pvHS + i);
            %(VTYPE)s vHL = %(VLOAD)s(pvHL + i);
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
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                %(INT)s *t = (%(INT)s*)pvH;
                %(INT)s *m = (%(INT)s*)pvHM;
                %(INT)s *s = (%(INT)s*)pvHS;
                %(INT)s *l = (%(INT)s*)pvHL;
                %(INDEX)s column_len = segLen * segWidth;
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
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(pvGapperL);
    parasail_free(pvGapper);
    parasail_free(pvHLMax);
    parasail_free(pvHSMax);
    parasail_free(pvHMMax);
    parasail_free(pvHMax);
    parasail_free(pvHL);
    parasail_free(pvHS);
    parasail_free(pvHM);
    parasail_free(pvH);
    parasail_free(pvEL);
    parasail_free(pvES);
    parasail_free(pvEM);
    parasail_free(pvE);

    return result;
}

