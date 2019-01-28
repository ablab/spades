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
#define ENAME parasail_sg_stats_table_scan
#else
#ifdef PARASAIL_ROWCOL
#define ENAME parasail_sg_stats_rowcol_scan
#else
#define ENAME parasail_sg_stats_scan
#endif
#endif

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif
    int * const restrict s1 = parasail_memalign_int(16, s1Len);
    int * const restrict s2 = parasail_memalign_int(16, s2Len);
    int * const restrict HB = parasail_memalign_int(16, s1Len+1);
    int * const restrict H  = HB+1;
    int * const restrict HMB = parasail_memalign_int(16, s1Len+1);
    int * const restrict HM  = HMB+1;
    int * const restrict HSB = parasail_memalign_int(16, s1Len+1);
    int * const restrict HS  = HSB+1;
    int * const restrict HLB = parasail_memalign_int(16, s1Len+1);
    int * const restrict HL  = HLB+1;
    int * const restrict E  = parasail_memalign_int(16, s1Len);
    int * const restrict EM = parasail_memalign_int(16, s1Len);
    int * const restrict ES = parasail_memalign_int(16, s1Len);
    int * const restrict EL = parasail_memalign_int(16, s1Len);
    int * const restrict HtB= parasail_memalign_int(16, s1Len+1);
    int * const restrict Ht = HtB+1;
    int * const restrict HtMB= parasail_memalign_int(16, s1Len+1);
    int * const restrict HtM = HtMB+1;
    int * const restrict HtSB= parasail_memalign_int(16, s1Len+1);
    int * const restrict HtS = HtSB+1;
    int * const restrict HtLB= parasail_memalign_int(16, s1Len+1);
    int * const restrict HtL = HtLB+1;
    int * const restrict Ex = parasail_memalign_int(16, s1Len);
    int i = 0;
    int j = 0;
    int score = NEG_INF_32;
    int matches = 0;
    int similar = 0;
    int length = 0;
    int end_query = s1Len;
    int end_ref = s2Len;

    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }

    /* initialize H,HM,HS,HL */
    H[-1] = 0;
    HM[-1] = 0;
    HS[-1] = 0;
    HL[-1] = 0;
    Ht[-1] = 0;
    for (i=0; i<s1Len; ++i) {
        H[i] = 0;
        HM[i] = 0;
        HS[i] = 0;
        HL[i] = 0;
    }

    /* initialize E */
    for (i=0; i<s1Len; ++i) {
        E[i] = NEG_INF_32;
        EM[i] = 0;
        ES[i] = 0;
        EL[i] = 0;
    }

    /* iterate over database */
    for (j=0; j<s2Len; ++j) {
        int Ft = NEG_INF_32;
        int FtM = 0;
        int FtS = 0;
        int FtL = 0;
        const int * const restrict matcol = &matrix->matrix[matrix->size*s2[j]];
        /* calculate E */
        for (i=0; i<s1Len; ++i) {
            int E_opn = H[i]-open;
            int E_ext = E[i]-gap;
            if (E_opn > E_ext) {
                E[i] = E_opn;
                EM[i] = HM[i];
                ES[i] = HS[i];
                EL[i] = HL[i];
            }
            else {
                E[i] = E_ext;
            }
            EL[i] += 1;
        }
        /* calculate Ht */
        Ht[-1] = 0;
        for (i=0; i<s1Len; ++i) {
            int H_dag = H[i-1]+matcol[s1[i]];
            Ex[i] = (E[i] > H_dag);
            if (H_dag >= E[i]) {
                Ht[i] = H_dag;
                HtM[i] = HM[i-1] + (s1[i]==s2[j]);
                HtS[i] = HS[i-1] + (matcol[s1[i]] > 0);
                HtL[i] = HL[i-1] + 1;
            }
            else {
                Ht[i] = E[i];
                HtM[i] = EM[i];
                HtS[i] = ES[i];
                HtL[i] = EL[i];
            }
        }
        /* calculate H,HM,HS,HL */
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
                FtM = HM[i-1];
                FtS = HS[i-1];
                FtL = HL[i-1];
            }
            FtL += 1;
            if (Ht[i] > Ft_opn) {
                H[i] = Ht[i];
                HM[i] = HtM[i];
                HS[i] = HtS[i];
                HL[i] = HtL[i];
            }
            else {
                H[i] = Ft_opn;
                if (Ht[i] == Ft_opn) {
                    if (Ex[i]) {
                        /* we favor F/up/del when F and E scores tie */
                        HM[i] = FtM;
                        HS[i] = FtS;
                        HL[i] = FtL;
                    }
                    else {
                        HM[i] = HtM[i];
                        HS[i] = HtS[i];
                        HL[i] = HtL[i];
                    }
                }
                else {
                    HM[i] = FtM;
                    HS[i] = FtS;
                    HL[i] = FtL;
                }
            }
#ifdef PARASAIL_TABLE
            result->stats->tables->score_table[1LL*i*s2Len + j] = H[i];
            result->stats->tables->matches_table[1LL*i*s2Len + j] = HM[i];
            result->stats->tables->similar_table[1LL*i*s2Len + j] = HS[i];
            result->stats->tables->length_table[1LL*i*s2Len + j] = HL[i];
#endif
        }
        H[-1] = 0;
#ifdef PARASAIL_ROWCOL
        if (j == s2Len-1) {
            for (i=0; i<s1Len; ++i) {
                result->stats->rowcols->score_col[i] = H[i];
                result->stats->rowcols->matches_col[i] = HM[i];
                result->stats->rowcols->similar_col[i] = HS[i];
                result->stats->rowcols->length_col[i] = HL[i];
            }
        }
        result->stats->rowcols->score_row[j] = H[s1Len-1];
        result->stats->rowcols->matches_row[j] = HM[s1Len-1];
        result->stats->rowcols->similar_row[j] = HS[s1Len-1];
        result->stats->rowcols->length_row[j] = HL[s1Len-1];
#endif
        /* last value from column */
        if (H[s1Len-1] > score) {
            score = H[s1Len-1];
            matches = HM[s1Len-1];
            similar = HS[s1Len-1];
            length = HL[s1Len-1];
            end_query = s1Len-1;
            end_ref = j;
        }
    }
    /* max of last column */
    for (i=0; i<s1Len; ++i) {
        if (H[i] > score) {
            score = H[i];
            matches = HM[i];
            similar = HS[i];
            length = HL[i];
            end_query = i;
            end_ref = s2Len-1;
        }
        else if (H[i] == score && end_ref == s2Len-1 && i < end_query) {
            matches = HM[i];
            similar = HS[i];
            length = HL[i];
            end_query = i;
            end_ref = s2Len-1;
        }
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->stats->matches = matches;
    result->stats->similar = similar;
    result->stats->length = length;
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_NOVEC_SCAN
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_INT | PARASAIL_FLAG_LANES_1;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(Ex);
    parasail_free(HtLB);
    parasail_free(HtSB);
    parasail_free(HtMB);
    parasail_free(HtB);
    parasail_free(EL);
    parasail_free(ES);
    parasail_free(EM);
    parasail_free(E);
    parasail_free(HLB);
    parasail_free(HSB);
    parasail_free(HMB);
    parasail_free(HB);
    parasail_free(s2);
    parasail_free(s1);

    return result;
}

