#include "config.h"

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

#include "func_verify_traces.h"

static int verbose = 0;

typedef struct gap_score {
    int open;
    int extend;
} gap_score_t;

gap_score_t gap_scores[] = {
    {9,1},
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
        unsigned long s1Len,
        unsigned long s2Len,
        int8_t *a,
        int8_t *b_,
        parasail_result_t *result)
{
    unsigned long i = 0;
    unsigned long size = s1Len * s2Len;
    if ((result->flag & PARASAIL_FLAG_TRACE)
            && ((result->flag & PARASAIL_FLAG_STRIPED)
            ||  (result->flag & PARASAIL_FLAG_SCAN))) {
        int *b = parasail_striped_unwind(s1Len, s2Len, result, b_);
        for (i=0; i<size; ++i) {
            if (a[i] != b[i]) return 1;
        }
        free(b);
    }
    else {
        for (i=0; i<size; ++i) {
            if (a[i] != b_[i]) return 1;
        }
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
            if (verbose) printf("\t\topen=%d extend=%d ref=%s\n", open, extend, functions[0].name);
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
                    int8_t *ref_trace_table = NULL;
                    int *trace_table = NULL;
                    size_t size_a = 0;
                    size_t size_b = 0;
                    k_combination2(pair_index, &a, &b);
                    size_a = sequences->seqs[a].seq.l;
                    size_b = sequences->seqs[b].seq.l;
                    if (verbose) printf("\t\t\t\tpair=%lld (%lu,%lu)\n", pair_index, a, b);
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
                    ref_trace_table = parasail_result_get_trace_table(reference_result);
                    trace_table = parasail_result_get_trace_table(result);
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
                    if (diff_array(size_a, size_b, ref_trace_table, trace_table, result)) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) bad trace table\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname);
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
    char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    gap_score_t gap = {INT_MIN,INT_MIN};
    int do_sse2 = 1;
    int do_sse41 = 1;
    int do_avx2 = 1;
    int do_altivec = 1;
    int do_neon = 1;
    int do_disp = 1;
    int do_nw = 1;
    int do_sg = 1;
    int do_sw = 1;

    while ((c = getopt(argc, argv, "f:m:n:o:e:vSi:")) != -1) {
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
            case 'S':
                test_scores = 0;
                break;
            case 'i':
                do_sse2 = (NULL == strstr(optarg, "sse2"));
                do_sse41 = (NULL == strstr(optarg, "sse41"));
                do_avx2 = (NULL == strstr(optarg, "avx2"));
                do_altivec = (NULL == strstr(optarg, "altivec"));
                do_neon = (NULL == strstr(optarg, "neon"));
                do_disp = (NULL == strstr(optarg, "disp"));
                do_nw = (NULL == strstr(optarg, "nw"));
                do_sg = (NULL == strstr(optarg, "sg"));
                do_sw = (NULL == strstr(optarg, "sw"));
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
    if (do_sse2 && parasail_can_use_sse2()) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_trace_sse2, sequences, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_trace_sse2, sequences, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_trace_sse2, sequences, limit, matrix, gap);
        }
    }
#endif

#if HAVE_SSE41
    if (do_sse41 && parasail_can_use_sse41()) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_trace_sse41, sequences, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_trace_sse41, sequences, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_trace_sse41, sequences, limit, matrix, gap);
        }
    }
#endif

#if HAVE_AVX2
    if (do_avx2 && parasail_can_use_avx2()) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_trace_avx2, sequences, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_trace_avx2, sequences, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_trace_avx2, sequences, limit, matrix, gap);
        }
    }
#endif

#if HAVE_ALTIVEC
    if (do_altivec && parasail_can_use_altivec()) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_trace_altivec, sequences, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_trace_altivec, sequences, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_trace_altivec, sequences, limit, matrix, gap);
        }
    }
#endif

#if HAVE_NEON
    if (do_neon && parasail_can_use_neon()) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_trace_neon, sequences, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_trace_neon, sequences, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_trace_neon, sequences, limit, matrix, gap);
        }
    }
#endif

    if (do_disp) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_trace_disp, sequences, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_trace_disp, sequences, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_trace_disp, sequences, limit, matrix, gap);
        }
    }
    
    parasail_sequences_free(sequences);

    return 0;
}

