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

#if defined(_OPENMP)
#include <omp.h>
#endif

#if HAVE_SSE2
#include "ssw.h"
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/io.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/stats.h"
/*#include "timer.h"*/
#include "timer_real.h"

#if HAVE_SSE2
static parasail_result_t* ssw_(
        const char * const restrict s1, const int s1_len,
        const char * const restrict s2, const int s2_len,
        const int open, const int gap, const parasail_matrix_t * pmatrix,
        int score_size)
{
    parasail_result_t *result = parasail_result_new();
    s_profile *profile = NULL;
    int8_t *s1_num = (int8_t*)malloc(sizeof(int8_t) * s1_len);
    int8_t *s2_num = (int8_t*)malloc(sizeof(int8_t) * s2_len);
    int8_t *matrix = (int8_t*)malloc(sizeof(int8_t) * 24 * 24);
    s_align *ssw_result = NULL;
    int m = 0;

    /* This table is used to transform amino acid letters into numbers. */
    static const int8_t table[128] = {
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
    };

    /* initialize score matrix */
    for (m = 0; m < s1_len; ++m) s1_num[m] = table[(int)s1[m]];
    for (m = 0; m < s2_len; ++m) s2_num[m] = table[(int)s2[m]];
    for (m = 0; m < 24*24; ++m) matrix[m] = pmatrix->matrix[m];
    profile = ssw_init(s1_num, s1_len, matrix, 24, score_size);
    ssw_result = ssw_align(profile, s2_num, s2_len, -open, -gap, 2, 0, 0, s1_len/2);
    result->score = ssw_result->score1;
    if (ssw_result->saturated) {
        result->flag |= PARASAIL_FLAG_SATURATED;
    }
    align_destroy(ssw_result);
    init_destroy(profile);
    free(s1_num);
    free(s2_num);
    free(matrix);

    return result;
}

static parasail_result_t* ssw(
        const char * const restrict s1, const int s1_len,
        const char * const restrict s2, const int s2_len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    return ssw_(s1, s1_len, s2, s2_len, open, gap, matrix, 2);
}

static parasail_result_t* ssw_16(
        const char * const restrict s1, const int s1_len,
        const char * const restrict s2, const int s2_len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    return ssw_(s1, s1_len, s2, s2_len, open, gap, matrix, 1);
}
#endif

static inline unsigned long binomial_coefficient(unsigned long n, unsigned long k)
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

static inline void k_combination2(unsigned long pos, unsigned long *a, unsigned long *b)
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

int main(int argc, char **argv)
{
    double timer_clock = 0.0;
    unsigned long i = 0;
    size_t limit = 0;
    char *filename_database = NULL;
    parasail_sequences_t *sequences_database = NULL;
    size_t seq_count_database = 0;
    char *filename_queries = NULL;
    parasail_sequences_t *sequences_queries = NULL;
    size_t seq_count_queries = 0;
    char *endptr = NULL;
    char *funcname1 = NULL;
    char *funcname2 = NULL;
    parasail_function_t *function1 = NULL;
    parasail_function_t *function2 = NULL;
    int banded = 0;
    int kbandsize = 3;
    int c = 0;
    char *matrixname = "blosum62";
    const parasail_matrix_t *matrix = NULL;
    int gap_open = 10;
    int gap_extend = 1;
    int N = 1;
    int saturated = 0;
    int smallest_first = 0;
    int biggest_first = 0;
    int truncate = 0;
    int iterations = 1;
    int func_cutoff = 0;
    int iter = 0;
    stats_t stats_time;

    stats_clear(&stats_time);

    while ((c = getopt(argc, argv, "a:A:c:b:f:k:q:o:e:slt:i:")) != -1) {
        switch (c) {
            case 'a':
                funcname1 = optarg;
                break;
            case 'A':
                funcname2 = optarg;
                break;
            case 'b':
                matrixname = optarg;
                break;
            case 'c':
                errno = 0;
                func_cutoff = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'f':
                filename_database = optarg;
                break;
            case 'q':
                filename_queries = optarg;
                break;
            case 'k':
                errno = 0;
                kbandsize = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'i':
                errno = 0;
                iterations = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'o':
                errno = 0;
                gap_open = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'e':
                errno = 0;
                gap_extend = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 's':
                smallest_first = 1;
                break;
            case 'l':
                biggest_first = 1;
                break;
            case 't':
                errno = 0;
                truncate = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'b'
                        || optopt == 'e'
                        || optopt == 'f'
                        || optopt == 'i'
                        || optopt == 'n'
                        || optopt == 'o'
                        || optopt == 't')
                {
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

    if (smallest_first && biggest_first) {
        fprintf(stderr, "cannot choose both smallest and biggest first\n");
        exit(1);
    }

    /* select the function */
    if (funcname1) {
        function1 = parasail_lookup_function(funcname1);
#if HAVE_SSE2
        if (NULL == function1) {
            if (0 == strcmp(funcname1, "ssw_16")) {
                function1 = ssw_16;
            }
            else if (0 == strcmp(funcname1, "ssw_8")) {
                function1 = ssw;
            }
        }
#endif
        if (NULL == function1 && NULL != strstr(funcname1, "nw_banded")) {
            banded = 1;
        }
        if (NULL == function1 && 0 == banded) {
            fprintf(stderr, "Specified function1 not found.\n");
            exit(1);
        }
    }
    else {
        fprintf(stderr, "No alignment function1 specified.\n");
        exit(1);
    }

    if (funcname2) {
        function2 = parasail_lookup_function(funcname2);
        if (NULL == function2) {
            fprintf(stderr, "Specified function2 not found.\n");
            exit(1);
        }
    }

    /* select the substitution matrix */
    if (matrixname) {
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            fprintf(stderr, "Specified substitution matrix not found.\n");
            exit(1);
        }
    }

    if (filename_database) {
        sequences_database = parasail_sequences_from_file(filename_database);
        seq_count_database = sequences_database->l;
    }
    else {
        fprintf(stderr, "missing database filename\n");
        exit(1);
    }

    limit = binomial_coefficient(seq_count_database, 2);
    /*printf("%lu choose 2 is %lu\n", seq_count_database, limit);*/

#if defined(_OPENMP)
#pragma omp parallel
    {
#pragma omp single
        {
            N = omp_get_max_threads();
            /*printf("omp_get_max_threads()=%d\n", N);*/
        }
    }
#endif

    if (filename_queries) {
        double total_timer = 0.0;
        sequences_queries = parasail_sequences_from_file(filename_queries);
        seq_count_queries = sequences_queries->l;
        for (i=0; i<seq_count_queries; ++i) {
            int saturated_query = 0;
            double local_timer = 0.0;
            parasail_function_t *function = function1;

            if (func_cutoff > 0) {
                if (sequences_queries->seqs[i].seq.l > (unsigned long)func_cutoff) {
                    function = function2;
                }
            }

            local_timer = timer_real();
#pragma omp parallel
            {
				long long int seq_count_database_signed = seq_count_database;
				long long int j_signed = 0;
#pragma omp for schedule(guided)
                for (j_signed=0; j_signed<seq_count_database_signed; ++j_signed) {
                    parasail_result_t *result = function(
                            sequences_queries->seqs[i].seq.s,
                            sequences_queries->seqs[i].seq.l,
                            sequences_database->seqs[j_signed].seq.s,
                            sequences_database->seqs[j_signed].seq.l,
                            gap_open, gap_extend, matrix);
#pragma omp atomic
                    saturated_query += parasail_result_is_saturated(result);
                    parasail_result_free(result);
                }
            }
            local_timer = timer_real() - local_timer;
            total_timer += local_timer;
            printf("%lu\t %lu\t %d\t %f\n",
                    i, (unsigned long)sequences_queries->seqs[i].seq.l,
                    saturated_query, local_timer);
            fflush(stdout);
        }
        printf("total_time=%f\n", total_timer);
        fflush(stdout);
    }
    else {
        for (iter=0; iter<iterations; ++iter) {
            timer_clock = timer_real();
#pragma omp parallel
            {
                unsigned long a=0;
                unsigned long b=1;
                unsigned long swap=0;
				long long i_signed = 0;
				long long limit_signed = limit;
#pragma omp for schedule(guided)
                for (i_signed=0; i_signed<limit_signed; ++i_signed) {
                    parasail_function_t *function = function1;
                    parasail_result_t *result = NULL;
                    unsigned long query_size;
                    k_combination2(i_signed, &a, &b);
                    if (smallest_first) {
                        if (sequences_database->seqs[a].seq.l > sequences_database->seqs[b].seq.l) {
                            swap = a;
                            a = b;
                            b = swap;
                        }
                    }
                    else if (biggest_first) {
                        if (sequences_database->seqs[a].seq.l < sequences_database->seqs[b].seq.l) {
                            swap = a;
                            a = b;
                            b = swap;
                        }
                    }
                    query_size = sequences_database->seqs[a].seq.l;
                    if (truncate > 0) {
                        if (query_size > (unsigned long)truncate) {
                            query_size = truncate;
                        }
                    }
                    if (func_cutoff > 0) {
                        if (query_size > (unsigned long)func_cutoff) {
                            function = function2;
                        }
                    }
                    if (NULL != function1) {
                        result = function(
                                sequences_database->seqs[a].seq.s, query_size,
                                sequences_database->seqs[b].seq.s, sequences_database->seqs[b].seq.l,
                                gap_open, gap_extend, matrix);
                    }
                    else if (0 != banded) {
                        result = parasail_nw_banded(
                                sequences_database->seqs[a].seq.s, query_size,
                                sequences_database->seqs[b].seq.s, sequences_database->seqs[b].seq.l,
                                gap_open, gap_extend, kbandsize, matrix);
                    }
#pragma omp atomic
                    saturated += parasail_result_is_saturated(result);
                    if (parasail_result_is_saturated(result)) {
#pragma omp critical(printer)
                        {
                            printf("saturated -- (%lu,%lu)\n", a, b);
                        }
                    }
                    parasail_result_free(result);
                }
            }
            timer_clock = timer_real() - timer_clock;
            stats_sample_value(&stats_time, timer_clock);
        }
        printf("%s\t %s\t %d\t %d\t %d\t %d\t %f\t %f\t %f\t %f\n",
                funcname1, matrixname, gap_open, gap_extend, N,
                saturated,
                stats_time._mean, stats_stddev(&stats_time),
                stats_time._min, stats_time._max);
        fflush(stdout);
    }

    return 0;
}

