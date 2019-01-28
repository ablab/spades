#include "config.h"

/* getopt needs _POSIX_C_SOURCE 2 */
#define _POSIX_C_SOURCE 2

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#if defined(_MSC_VER)
#include "wingetopt/src/getopt.h"
#else
#include <unistd.h>
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/cpuid.h"
#include "parasail/parasail/io.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/matrix_lookup.h"

#include "func_verify_rowcols.h"

static int verbose = 0;

typedef struct gap_score {
    int open;
    int extend;
} gap_score_t;

gap_score_t gap_scores[] = {
    {10,1},
    {10,2},
    {14,2},
    {40,2},
    {INT_MIN,INT_MIN}
};

static inline unsigned long binomial_coefficient(
        unsigned long n,
        unsigned long k)
{
    /* from http://blog.plover.com/math/choose.html */
    unsigned long r = 1;
    unsigned long d;
    if (k > n) {
        return 0;
    }
    for (d = 1; d <= k; d++) {
        r *= n--;
        r /= d;
    }
    return r;
}

static inline void k_combination2(
        unsigned long pos,
        unsigned long *a,
        unsigned long *b)
{
    double s;
    double i = floor(sqrt(2.0 * pos)) - 1.0;
    if (i <= 1.0) {
        i = 1.0;
    }
    s = i * (i - 1.0) / 2.0;
    while (pos - s >= i) {
        s += i;
        i += 1;
    }
    *a = (unsigned long)(pos - s);
    *b = (unsigned long)(i);
}

static inline int diff_array(
        unsigned long size,
        int *a,
        int *b)
{
    unsigned long i = 0;
    for (i=0; i<size; ++i) {
        if (a[i] != b[i]) return 1;
    }
    return 0;
}

static void check_functions(
        parasail_function_group_t f,
        parasail_sequences_t *sequences,
        unsigned long pair_limit_,
        const parasail_matrix_t *matrix_,
        gap_score_t gap)
{
    const parasail_function_info_t *functions = f.fs;
    unsigned long matrix_index = 0;
    unsigned long gap_index = 0;
    unsigned long function_index = 0;
    long long pair_index = 0;
    long long pair_limit = (long long)pair_limit_;
    parasail_function_t *reference_function = NULL;
    const parasail_matrix_t ** matrices = parasail_matrices;
    const parasail_matrix_t * single_matrix[] = {
        matrix_,
        NULL
    };

    if (NULL != matrix_) {
        matrices = single_matrix;
    }

    printf("checking %s functions\n", f.name);
    for (matrix_index=0; NULL!=matrices[matrix_index]; ++matrix_index) {
        const parasail_matrix_t *matrix = matrices[matrix_index];
        const char *matrixname = matrix->name;
        if (verbose) printf("\t%s\n", matrixname);
        for (gap_index=0; INT_MIN!=gap_scores[gap_index].open; ++gap_index) {
            int open = gap_scores[gap_index].open;
            int extend = gap_scores[gap_index].extend;
            if (gap.open != INT_MIN && gap.extend != INT_MIN) {
                open = gap.open;
                extend = gap.extend;
            }
            if (verbose) printf("\t\topen=%d extend=%d\n", open, extend);
            reference_function = functions[0].pointer;
            for (function_index=1;
                    NULL!=functions[function_index].pointer;
                    ++function_index) {
                unsigned long saturated = 0;
                if (verbose) printf("\t\t\t%s\n", functions[function_index].name);
#pragma omp parallel for
                for (pair_index=0; pair_index<pair_limit; ++pair_index) {
                    parasail_result_t *reference_result = NULL;
                    parasail_result_t *result = NULL;
                    unsigned long a = 0;
                    unsigned long b = 1;
                    int *ref_score_row = NULL;
                    int *ref_score_col = NULL;
                    int *score_row = NULL;
                    int *score_col = NULL;
                    size_t size_a = 0;
                    size_t size_b = 0;
                    k_combination2(pair_index, &a, &b);
                    size_a = sequences->seqs[a].seq.l;
                    size_b = sequences->seqs[b].seq.l;
                    /*printf("\t\t\t\tpair=%lld (%lu,%lu)\n", pair_index, a, b);*/
                    reference_result = reference_function(
                            sequences->seqs[a].seq.s, size_a,
                            sequences->seqs[b].seq.s, size_b,
                            open, extend,
                            matrix);
                    result = functions[function_index].pointer(
                            sequences->seqs[a].seq.s, size_a,
                            sequences->seqs[b].seq.s, size_b,
                            open, extend,
                            matrix);
                    if (parasail_result_is_saturated(result)) {
                        /* no point in comparing a result that saturated */
                        parasail_result_free(reference_result);
                        parasail_result_free(result);
#pragma omp atomic
                        saturated += 1;
                        continue;
                    }
                    ref_score_row = parasail_result_get_score_row(reference_result);
                    ref_score_col = parasail_result_get_score_col(reference_result);
                    score_row = parasail_result_get_score_row(result);
                    score_col = parasail_result_get_score_col(result);
                    if (reference_result->score != result->score) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) wrong score (%d!=%d)\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname,
                                    reference_result->score, result->score);
                        }
                    }
                    if (diff_array(size_b, ref_score_row, score_row)) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) bad score row\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname);
                        }
                    }
                    if (diff_array(size_a, ref_score_col, score_col)) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) bad score col\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname);
                        }
                    }
                    if (parasail_result_is_stats(result)) {
                        int *ref_matches_row = parasail_result_get_matches_row(reference_result);
                        int *ref_matches_col = parasail_result_get_matches_col(reference_result);
                        int *ref_similar_row = parasail_result_get_similar_row(reference_result);
                        int *ref_similar_col = parasail_result_get_similar_col(reference_result);
                        int *ref_length_row = parasail_result_get_length_row(reference_result);
                        int *ref_length_col = parasail_result_get_length_col(reference_result);
                        int *matches_row = parasail_result_get_matches_row(result);
                        int *matches_col = parasail_result_get_matches_col(result);
                        int *similar_row = parasail_result_get_similar_row(result);
                        int *similar_col = parasail_result_get_similar_col(result);
                        int *length_row = parasail_result_get_length_row(result);
                        int *length_col = parasail_result_get_length_col(result);
                        if (diff_array(size_b, ref_matches_row, matches_row)) {
#pragma omp critical(printer)
                            {
                                printf("%s(%lu,%lu,%d,%d,%s) bad matches row\n",
                                        functions[function_index].name,
                                        a, b, open, extend,
                                        matrixname);
                            }
                        }
                        if (diff_array(size_a, ref_matches_col, matches_col)) {
#pragma omp critical(printer)
                            {
                                printf("%s(%lu,%lu,%d,%d,%s) bad matches col\n",
                                        functions[function_index].name,
                                        a, b, open, extend,
                                        matrixname);
                            }
                        }
                        if (diff_array(size_b, ref_similar_row, similar_row)) {
#pragma omp critical(printer)
                            {
                                printf("%s(%lu,%lu,%d,%d,%s) bad similar row\n",
                                        functions[function_index].name,
                                        a, b, open, extend,
                                        matrixname);
                            }
                        }
                        if (diff_array(size_a, ref_similar_col, similar_col)) {
#pragma omp critical(printer)
                            {
                                printf("%s(%lu,%lu,%d,%d,%s) bad similar col\n",
                                        functions[function_index].name,
                                        a, b, open, extend,
                                        matrixname);
                            }
                        }
                        if (diff_array(size_b, ref_length_row, length_row)) {
#pragma omp critical(printer)
                            {
                                printf("%s(%lu,%lu,%d,%d,%s) bad length row\n",
                                        functions[function_index].name,
                                        a, b, open, extend,
                                        matrixname);
                            }
                        }
                        if (diff_array(size_a, ref_length_col, length_col)) {
#pragma omp critical(printer)
                            {
                                printf("%s(%lu,%lu,%d,%d,%s) bad length col\n",
                                        functions[function_index].name,
                                        a, b, open, extend,
                                        matrixname);
                            }
                        }
                    }
                    parasail_result_free(reference_result);
                    parasail_result_free(result);
                }
                if (verbose && saturated) {
                    printf("%s %d %d %s saturated %lu times\n",
                            functions[function_index].name,
                            open, extend,
                            matrixname,
                            saturated);
                }
            }
            if (gap.open != INT_MIN && gap.extend != INT_MIN) {
                /* user-specified gap, don't loop */
                break;
            }
        }
    }
}

int main(int argc, char **argv)
{
    unsigned long seq_count = 0;
    unsigned long limit = 0;
    parasail_sequences_t *sequences = NULL;
    char *endptr = NULL;
    char *filename = NULL;
    int c = 0;
    int test_scores = 1;
    int test_stats = 0;
    char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    gap_score_t gap = {INT_MIN,INT_MIN};

    while ((c = getopt(argc, argv, "f:m:n:o:e:vsS")) != -1) {
        switch (c) {
            case 'f':
                filename = optarg;
                break;
            case 'm':
                matrixname = optarg;
                break;
            case 'n':
                errno = 0;
                seq_count = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'o':
                errno = 0;
                gap.open = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol gap.open");
                    exit(1);
                }
                break;
            case 'e':
                errno = 0;
                gap.extend = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol gap.extend");
                    exit(1);
                }
                break;
            case 'v':
                verbose = 1;
                break;
            case 's':
                test_stats = 1;
                break;
            case 'S':
                test_scores = 0;
                break;
            case '?':
                if (optopt == 'f' || optopt == 'n') {
                    fprintf(stderr,
                            "Option -%c requires an argument.\n",
                            optopt);
                }
                else if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option `-%c'.\n",
                            optopt);
                }
                else {
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                }
                exit(1);
            default:
                fprintf(stderr, "default case in getopt\n");
                exit(1);
        }
    }

    if (filename) {
        sequences = parasail_sequences_from_file(filename);
        if (0 == seq_count) {
            seq_count = sequences->l;
        }
    }
    else {
        fprintf(stderr, "no filename specified\n");
        exit(1);
    }

    /* select the matrix */
    if (matrixname) {
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            fprintf(stderr, "Specified substitution matrix not found.\n");
            exit(1);
        }
    }

    limit = binomial_coefficient(seq_count, 2);
    printf("%lu choose 2 is %lu\n", seq_count, limit);


#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        if (test_scores) {
            check_functions(parasail_nw_rowcol_sse2, sequences, limit, matrix, gap);
            check_functions(parasail_sg_rowcol_sse2, sequences, limit, matrix, gap);
            check_functions(parasail_sw_rowcol_sse2, sequences, limit, matrix, gap);
        }
        if (test_stats) {
            check_functions(parasail_nw_stats_rowcol_sse2, sequences, limit, matrix, gap);
            check_functions(parasail_sg_stats_rowcol_sse2, sequences, limit, matrix, gap);
            check_functions(parasail_sw_stats_rowcol_sse2, sequences, limit, matrix, gap);
        }
    }
#endif

#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        if (test_scores) {
            check_functions(parasail_nw_rowcol_sse41, sequences, limit, matrix, gap);
            check_functions(parasail_sg_rowcol_sse41, sequences, limit, matrix, gap);
            check_functions(parasail_sw_rowcol_sse41, sequences, limit, matrix, gap);
        }
        if (test_stats) {
            check_functions(parasail_nw_stats_rowcol_sse41, sequences, limit, matrix, gap);
            check_functions(parasail_sg_stats_rowcol_sse41, sequences, limit, matrix, gap);
            check_functions(parasail_sw_stats_rowcol_sse41, sequences, limit, matrix, gap);
        }
    }
#endif

#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        if (test_scores) {
            check_functions(parasail_nw_rowcol_avx2, sequences, limit, matrix, gap);
            check_functions(parasail_sg_rowcol_avx2, sequences, limit, matrix, gap);
            check_functions(parasail_sw_rowcol_avx2, sequences, limit, matrix, gap);
        }
        if (test_stats) {
            check_functions(parasail_nw_stats_rowcol_avx2, sequences, limit, matrix, gap);
            check_functions(parasail_sg_stats_rowcol_avx2, sequences, limit, matrix, gap);
            check_functions(parasail_sw_stats_rowcol_avx2, sequences, limit, matrix, gap);
        }
    }
#endif

#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        if (test_scores) {
            check_functions(parasail_nw_rowcol_altivec, sequences, limit, matrix, gap);
            check_functions(parasail_sg_rowcol_altivec, sequences, limit, matrix, gap);
            check_functions(parasail_sw_rowcol_altivec, sequences, limit, matrix, gap);
        }
        if (test_stats) {
            check_functions(parasail_nw_stats_rowcol_altivec, sequences, limit, matrix, gap);
            check_functions(parasail_sg_stats_rowcol_altivec, sequences, limit, matrix, gap);
            check_functions(parasail_sw_stats_rowcol_altivec, sequences, limit, matrix, gap);
        }
    }
#endif

#if HAVE_NEON
    if (parasail_can_use_neon()) {
        if (test_scores) {
            check_functions(parasail_nw_rowcol_neon, sequences, limit, matrix, gap);
            check_functions(parasail_sg_rowcol_neon, sequences, limit, matrix, gap);
            check_functions(parasail_sw_rowcol_neon, sequences, limit, matrix, gap);
        }
        if (test_stats) {
            check_functions(parasail_nw_stats_rowcol_neon, sequences, limit, matrix, gap);
            check_functions(parasail_sg_stats_rowcol_neon, sequences, limit, matrix, gap);
            check_functions(parasail_sw_stats_rowcol_neon, sequences, limit, matrix, gap);
        }
    }
#endif

    if (test_scores) {
        check_functions(parasail_nw_rowcol_disp, sequences, limit, matrix, gap);
        check_functions(parasail_sg_rowcol_disp, sequences, limit, matrix, gap);
        check_functions(parasail_sw_rowcol_disp, sequences, limit, matrix, gap);
    }
    if (test_stats) {
        check_functions(parasail_nw_stats_rowcol_disp, sequences, limit, matrix, gap);
        check_functions(parasail_sg_stats_rowcol_disp, sequences, limit, matrix, gap);
        check_functions(parasail_sw_stats_rowcol_disp, sequences, limit, matrix, gap);
    }
    
    parasail_sequences_free(sequences);

    return 0;
}

