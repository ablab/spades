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
#include "parasail/parasail/function_lookup.h"
#include "parasail/parasail/io.h"
/*#include "timer.h"*/
#include "timer_real.h"

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

static inline char* rand_string(size_t size)
{
    char *str = NULL;
    const char charset[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    if (size) {
        size_t n;
        --size;
        str = malloc(size + 1);
        for (n = 0; n < size; n++) {
            int key = rand() % (int) (sizeof charset - 1);
            str[n] = charset[key];
        }
        str[size] = '\0';
    }
    return str;
}

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
    unsigned long shortest = INT_MAX;
    unsigned long longest = 0;
    double timer_clock = 0.0;
    unsigned long i = 0;
    size_t seq_count = 10;
    size_t limit = 0;
    parasail_sequences_t *sequences = NULL;
    char *endptr = NULL;
    char *funcname = NULL;
    parasail_function_t *function = NULL;
    int lanes = 1;
    char *filename = NULL;
    int c = 0;
    const char *matrixname = "blosum62";
    const parasail_matrix_t *matrix = NULL;
    int gap_open = 10;
    int gap_extend = 1;
    int saturated = 0;

    while ((c = getopt(argc, argv, "a:b:f:n:o:e:")) != -1) {
        switch (c) {
            case 'a':
                funcname = optarg;
                break;
            case 'b':
                matrixname = optarg;
                break;
            case 'f':
                filename = optarg;
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

    /* select the function */
    if (funcname) {
        int index = 0;
        parasail_function_info_t f;
        f = functions[index++];
        while (f.pointer) {
            if (0 == strcmp(funcname, f.name)) {
                function = f.pointer;
                lanes = f.lanes;
                break;
            }
            f = functions[index++];
        }
        if (NULL == function) {
            fprintf(stderr, "Specified function not found.\n");
            exit(1);
        }
    }
    else {
        fprintf(stderr, "No alignment function specified.\n");
        exit(1);
    }

    /* select the substitution matrix */
    if (matrixname) {
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            fprintf(stderr, "Specified substitution matrix not found.\n");
            exit(1);
        }
    }

    if (filename) {
        sequences = parasail_sequences_from_file(filename);
        seq_count = sequences->l;
    }
    else {
        /* generate 'seq_count' number of random strings */
        sequences = (parasail_sequences_t*)malloc(sizeof(parasail_sequences_t));
        sequences->l = seq_count;
        sequences->seqs = (parasail_sequence_t*)malloc(sizeof(parasail_sequence_t)*seq_count);
        for (i=0; i<seq_count; ++i) {
            size_t size = (rand()%32767)+10;
            sequences->seqs[i].seq.l = size;
            shortest = size < shortest ? size : shortest;
            longest = size > longest ? size : longest;
            sequences->seqs[i].seq.s = rand_string(size);
        }
    }

    limit = binomial_coefficient(seq_count, 2);

    printf("size_A,segLen,size_B,score,matches,similar,length,corrections,cells,time,");
    printf("A_a,R_a,N_a,D_a,C_a,Q_a,E_a,G_a,H_a,I_a,L_a,K_a,M_a,F_a,P_a,S_a,T_a,W_a,Y_a,V_a,B_a,Z_a,X_a,NA_a,");
    printf("A_b,R_b,N_b,D_b,C_b,Q_b,E_b,G_b,H_b,I_b,L_b,K_b,M_b,F_b,P_b,S_b,T_b,W_b,Y_b,V_b,B_b,Z_b,X_b,NA_b,");
    printf("CUPS\n");

    timer_clock = timer_real();
#pragma omp parallel
    {
        unsigned long a=0;
        unsigned long b=1;
        double timer_local = 0.0;
        unsigned long a_counts[24];
        unsigned long b_counts[24];
        unsigned long j;
#pragma omp for schedule(dynamic)
        for (i=0; i<limit; ++i) {
            parasail_result_t *result = NULL;
            char *seq_a = sequences->seqs[a].seq.s;
            char *seq_b = sequences->seqs[b].seq.s;
            size_t size_a = 0;
            size_t size_b = 0;
            k_combination2(i, &a, &b);
            size_a = sequences->seqs[a].seq.l;
            size_b = sequences->seqs[b].seq.l;
            timer_local = timer_real();
            result = function(seq_a, size_a, seq_b, size_b,
                    gap_open, gap_extend, matrix);
            timer_local = timer_real() - timer_local;
            for (j=0; j<24; ++j) {
                a_counts[j] = 0;
                b_counts[j] = 0;
            }
            for (j=0; j<size_a; ++j) {
                a_counts[table[(unsigned)seq_a[j]]] += 1;
            }
            for (j=0; j<size_b; ++j) {
                b_counts[table[(unsigned)seq_b[j]]] += 1;
            }
#pragma omp critical
            printf("%lu,%lu,%lu,%d,%d,%d,%d,%lu,%f",
                    (unsigned long)size_a,
                    (unsigned long)(size_a+lanes-1)/lanes,
                    (unsigned long)size_b,
                    parasail_result_get_score(result),
                    parasail_result_get_matches(result),
                    parasail_result_get_similar(result),
                    parasail_result_get_length(result),
                    (unsigned long)(size_a*size_b),
                    timer_local);
            for (j=0; j<24; ++j) {
                /*printf(",%lu", a_counts[j]);*/
                printf(",%f", (double)(a_counts[j])/size_a);
            }
            for (j=0; j<24; ++j) {
                /*printf(",%lu", b_counts[j]);*/
                printf(",%f", (double)(b_counts[j])/size_b);
            }
            printf(",%f\n", size_a*size_b/timer_local);
#pragma omp atomic
            saturated += parasail_result_is_saturated(result);
            parasail_result_free(result);
        }
    }
    timer_clock = timer_real() - timer_clock;

    return 0;
}

