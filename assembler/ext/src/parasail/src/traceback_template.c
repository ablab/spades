/* this file is only included by traceback.c, multiple times */

#define D CONCAT3(int, T, _t)

#if defined(STRIPED)
#define NAME parasail_result_get_traceback_striped_
#define LOC LOC_STRIPED
#else
#define NAME parasail_result_get_traceback_
#define LOC LOC_NOVEC
#endif

static inline parasail_traceback_t* CONCAT(NAME, T) (
        parasail_result_t *result,
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix,
        char match, char pos, char neg)
{
    parasail_traceback_t *traceback = NULL;
    char *q = malloc(sizeof(char)*(lena+lenb));
    char *d = malloc(sizeof(char)*(lena+lenb));
    char *a = malloc(sizeof(char)*(lena+lenb));
    char *qc = q;
    char *dc = d;
    char *ac = a;
    int64_t i = result->end_query;
    int64_t j = result->end_ref;
    int where = PARASAIL_DIAG;
    D *HT = (D*)result->trace->trace_table;
#if defined(STRIPED)
    int64_t segWidth = 0;
    int64_t segLen = 0;
    if (result->flag & PARASAIL_FLAG_LANES_1) {
        segWidth = 1;
    }
    if (result->flag & PARASAIL_FLAG_LANES_2) {
        segWidth = 2;
    }
    if (result->flag & PARASAIL_FLAG_LANES_4) {
        segWidth = 4;
    }
    if (result->flag & PARASAIL_FLAG_LANES_8) {
        segWidth = 8;
    }
    if (result->flag & PARASAIL_FLAG_LANES_16) {
        segWidth = 16;
    }
    if (result->flag & PARASAIL_FLAG_LANES_32) {
        segWidth = 32;
    }
    if (result->flag & PARASAIL_FLAG_LANES_64) {
        segWidth = 64;
    }
    segLen = (lena + segWidth - 1) / segWidth;
#endif
    /* semi-global alignment includes the end gaps */
    if (result->flag & PARASAIL_FLAG_SG) {
        int k;
        if (result->end_query+1 == lena) {
            k = lenb-1;
            while (k > j) {
                *(qc++) = '-';
                *(dc++) = seqB[k];
                *(ac++) = ' ';
                --k;
            }
        }
        else if (result->end_ref+1 == lenb) {
            k = lena-1;
            while (k > i) {
                *(qc++) = seqA[k];
                *(dc++) = '-';
                *(ac++) = ' ';
                --k;
            }
        }
        else {
            assert(0);
        }
    }
    while (i >= 0 || j >= 0) {
        LOC
        if (i < 0) {
            if (!(result->flag & PARASAIL_FLAG_SW)) {
                while (j >= 0) {
                    *(qc++) = '-';
                    *(dc++) = seqB[j];
                    *(ac++) = ' ';
                    --j;
                }
            }
            break;
        }
        if (j < 0) {
            if (!(result->flag & PARASAIL_FLAG_SW)) {
                while (i >= 0) {
                    *(qc++) = seqA[i];
                    *(dc++) = '-';
                    *(ac++) = ' ';
                    --i;
                }
            }
            break;
        }
        if (PARASAIL_DIAG == where) {
            if (HT[loc] & PARASAIL_DIAG) {
                *(qc++) = seqA[i];
                *(dc++) = seqB[j];
                *(ac++) = match_char(seqA[i], seqB[j], matrix, match, pos, neg);
                --i;
                --j;
            }
            else if (HT[loc] & PARASAIL_INS) {
                where = PARASAIL_INS;
            }
            else if (HT[loc] & PARASAIL_DEL) {
                where = PARASAIL_DEL;
            }
            else {
                break;
            }
        }
        else if (PARASAIL_INS == where) {
            *(qc++) = '-';
            *(dc++) = seqB[j];
            *(ac++) = ' ';
            --j;
            if (HT[loc] & PARASAIL_DIAG_E) {
                where = PARASAIL_DIAG;
            }
            else if (HT[loc] & PARASAIL_INS_E) {
                where = PARASAIL_INS;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_DEL == where) {
            *(qc++) = seqA[i];
            *(dc++) = '-';
            *(ac++) = ' ';
            --i;
            if (HT[loc] & PARASAIL_DIAG_F) {
                where = PARASAIL_DIAG;
            }
            else if (HT[loc] & PARASAIL_DEL_F) {
                where = PARASAIL_DEL;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_ZERO == where) {
            break;
        }
        else {
            assert(0);
        }
    }
    *(qc++) = '\0';
    *(dc++) = '\0';
    *(ac++) = '\0';

    traceback = malloc(sizeof(parasail_traceback_t));
    traceback->query = parasail_reverse(q, strlen(q));
    traceback->comp = parasail_reverse(a, strlen(a));
    traceback->ref = parasail_reverse(d, strlen(d));

    free(q);
    free(d);
    free(a);

    return traceback;
}

#undef NAME
#undef LOC

#if defined(STRIPED)
#define NAME parasail_traceback_striped_
#define LOC LOC_STRIPED
#else
#define NAME parasail_traceback_
#define LOC LOC_NOVEC
#endif

static inline void CONCAT(NAME, T) (
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const char *nameA,
        const char *nameB,
        const parasail_matrix_t *matrix,
        parasail_result_t *result,
        char match, char pos, char neg,
        int width,
        int name_width,
        int use_stats,
        int int_width,
        FILE *stream)
{
    char *q = malloc(sizeof(char)*(lena+lenb));
    char *d = malloc(sizeof(char)*(lena+lenb));
    char *a = malloc(sizeof(char)*(lena+lenb));
    char *qc = q;
    char *dc = d;
    char *ac = a;
    int64_t i = result->end_query;
    int64_t j = result->end_ref;
    int where = PARASAIL_DIAG;
    int64_t c_ins = 0;
    int64_t c_del = 0;
    D *HT = (D*)result->trace->trace_table;
    int64_t namelenA = (NULL == nameA) ? 0 : (int64_t)strlen(nameA);
    int64_t namelenB = (NULL == nameB) ? 0 : (int64_t)strlen(nameB);
    char tmp[32];
    int _int_width = 0;
#if defined(STRIPED)
    int64_t segWidth = 0;
    int64_t segLen = 0;
    if (result->flag & PARASAIL_FLAG_LANES_1) {
        segWidth = 1;
    }
    if (result->flag & PARASAIL_FLAG_LANES_2) {
        segWidth = 2;
    }
    if (result->flag & PARASAIL_FLAG_LANES_4) {
        segWidth = 4;
    }
    if (result->flag & PARASAIL_FLAG_LANES_8) {
        segWidth = 8;
    }
    if (result->flag & PARASAIL_FLAG_LANES_16) {
        segWidth = 16;
    }
    if (result->flag & PARASAIL_FLAG_LANES_32) {
        segWidth = 32;
    }
    if (result->flag & PARASAIL_FLAG_LANES_64) {
        segWidth = 64;
    }
    segLen = (lena + segWidth - 1) / segWidth;
#endif
    /* how wide does our index label need to be? */
    _int_width = snprintf(tmp, 32, "%llu", (unsigned long long)lena+(unsigned long long)lenb);
    /* if requested int width is too small, override it */
    if (int_width < _int_width) {
        int_width = _int_width;
    }
    /* semi-global alignment includes the end gaps */
    if (result->flag & PARASAIL_FLAG_SG) {
        int k;
        if (result->end_query+1 == lena) {
            k = lenb-1;
            while (k > j) {
                ++c_ins;
                *(qc++) = '-';
                *(dc++) = seqB[k];
                *(ac++) = ' ';
                --k;
            }
        }
        else if (result->end_ref+1 == lenb) {
            k = lena-1;
            while (k > i) {
                ++c_del;
                *(qc++) = seqA[k];
                *(dc++) = '-';
                *(ac++) = ' ';
                --k;
            }
        }
        else {
            assert(0);
        }
    }
    while (i >= 0 || j >= 0) {
        LOC
        if (i < 0) {
            if (!(result->flag & PARASAIL_FLAG_SW)) {
                while (j >= 0) {
                    ++c_ins;
                    *(qc++) = '-';
                    *(dc++) = seqB[j];
                    *(ac++) = ' ';
                    --j;
                }
            }
            break;
        }
        if (j < 0) {
            if (!(result->flag & PARASAIL_FLAG_SW)) {
                while (i >= 0) {
                    ++c_del;
                    *(qc++) = seqA[i];
                    *(dc++) = '-';
                    *(ac++) = ' ';
                    --i;
                }
            }
            break;
        }
        if (PARASAIL_DIAG == where) {
            if (HT[loc] & PARASAIL_DIAG) {
                *(qc++) = seqA[i];
                *(dc++) = seqB[j];
                *(ac++) = match_char(seqA[i], seqB[j], matrix, match, pos, neg);
                --i;
                --j;
            }
            else if (HT[loc] & PARASAIL_INS) {
                where = PARASAIL_INS;
            }
            else if (HT[loc] & PARASAIL_DEL) {
                where = PARASAIL_DEL;
            }
            else {
                break;
            }
        }
        else if (PARASAIL_INS == where) {
            ++c_ins;
            *(qc++) = '-';
            *(dc++) = seqB[j];
            *(ac++) = ' ';
            --j;
            if (HT[loc] & PARASAIL_DIAG_E) {
                where = PARASAIL_DIAG;
            }
            else if (HT[loc] & PARASAIL_INS_E) {
                where = PARASAIL_INS;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_DEL == where) {
            ++c_del;
            *(qc++) = seqA[i];
            *(dc++) = '-';
            *(ac++) = ' ';
            --i;
            if (HT[loc] & PARASAIL_DIAG_F) {
                where = PARASAIL_DIAG;
            }
            else if (HT[loc] & PARASAIL_DEL_F) {
                where = PARASAIL_DEL;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_ZERO == where) {
            break;
        }
        else {
            assert(0);
        }
    }
    *(qc++) = '\0';
    *(dc++) = '\0';
    *(ac++) = '\0';

    if (1) {
        char *qr = NULL;
        char *ar = NULL;
        char *dr = NULL;
        int64_t mch = 0;
        int64_t sim = 0;
        int64_t gap = 0;
        int64_t len = strlen(a);
        int64_t q_pindex = 0;
        int64_t d_pindex = 0;
        int64_t qi = 0;
        int64_t ai = 0;
        int64_t di = 0;

        if (result->flag & PARASAIL_FLAG_SW) {
            q_pindex = result->end_query + 1 - len + c_ins;
            d_pindex = result->end_ref + 1 - len + c_del;
        }
        else {
            q_pindex = 0;
            d_pindex = 0;
        }
        qr = parasail_reverse(q, strlen(q));
        ar = parasail_reverse(a, strlen(a));
        dr = parasail_reverse(d, strlen(d));
        for (i=0; i<len; i+=width) {
            fprintf(stream, "\n");
            for (j=0; j<name_width; ++j) {
                if (j >= namelenB) break;
                fprintf(stream, "%c", nameB[j]);
            }
            for (; j<name_width; ++j) {
                fprintf(stream, " ");
            }
            fprintf(stream, " %*lld ", int_width, (long long)d_pindex+1);
            for (j=0; j<len&&j<width&&di<len; ++j) {
                if (dr[di] != '-') ++d_pindex;
                fprintf(stream, "%c", dr[di]);
                ++di;
            }
            fprintf(stream, " %*lld\n", int_width, (long long)d_pindex);
            for (j=0; j<name_width+1+int_width+1; ++j) {
                fprintf(stream, " ");
            }
            for (j=0; j<len&&j<width&&ai<len; ++j) {
                if (ar[ai] == match) { ++mch; ++sim; }
                else if (ar[ai] == pos) ++sim;
                else if (ar[ai] == neg) ;
                else if (ar[ai] == ' ') ++gap;
                else {
                    fprintf(stderr, "bad char in traceback '%c'\n", ar[ai]);
                    assert(0);
                }
                fprintf(stream, "%c", ar[ai]);
                ++ai;
            }
            fprintf(stream, "\n");
            for (j=0; j<name_width; ++j) {
                if (j >= namelenA) break;
                fprintf(stream, "%c", nameA[j]);
            }
            for (; j<name_width; ++j) {
                fprintf(stream, " ");
            }
            fprintf(stream, " %*lld ", int_width, (long long)q_pindex+1);
            for (j=0; j<len&&j<width&&qi<len; ++j) {
                if (qr[qi] != '-') ++q_pindex;
                fprintf(stream, "%c", qr[qi]);
                ++qi;
            }
            fprintf(stream, " %*lld\n", int_width, (long long)q_pindex);
        }
        if (use_stats) {
            fprintf(stream, "\n");
            fprintf(stream, "Length: %lld\n", (long long)len);
            fprintf(stream, "Identity:   %*lld/%lld (%4.1f%%)\n", int_width, (long long)mch, (long long)len, 100.0*mch/len);
            fprintf(stream, "Similarity: %*lld/%lld (%4.1f%%)\n", int_width, (long long)sim, (long long)len, 100.0*sim/len);
            fprintf(stream, "Gaps:       %*lld/%lld (%4.1f%%)\n", int_width, (long long)gap, (long long)len, 100.0*gap/len);
            fprintf(stream, "Score: %d\n", result->score);
        }
        free(qr);
        free(ar);
        free(dr);
    }
    else {
        fprintf(stream, "%s\n", q);
        fprintf(stream, "%s\n", a);
        fprintf(stream, "%s\n", d);
    }

    free(q);
    free(d);
    free(a);
}

#undef NAME
#undef LOC
#undef D

