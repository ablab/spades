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
#include <unistd.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/io.h"
#include "parasail/parasail/stats.h"
/*#include "timer.h"*/
#include "timer_real.h"


int main(int argc, char **argv)
{
    unsigned long i = 0;
    unsigned long j = 0;
    char *filename_database = NULL;
    parasail_sequences_t *sequences_database = NULL;
    size_t seq_count_database = 0;
    char *filename_queries = NULL;
    parasail_sequences_t *sequences_queries = NULL;
    size_t seq_count_queries = 0;
    char *endptr = NULL;
    char *funcname = NULL;
    parasail_function_t *function = NULL;
    int c = 0;
    const char *matrixname = "blosum62";
    const parasail_matrix_t *matrix = NULL;
    int gap_open = 10;
    int gap_extend = 1;
    int truncate = 0;
    int exact_length = 0;
    stats_t stats_time;
    size_t biggest = 0;
    char *seen = NULL;

    stats_clear(&stats_time);

    while ((c = getopt(argc, argv, "a:b:f:q:o:e:t:x:")) != -1) {
        switch (c) {
            case 'a':
                funcname = optarg;
                break;
            case 'b':
                matrixname = optarg;
                break;
            case 'f':
                filename_database = optarg;
                break;
            case 'q':
                filename_queries = optarg;
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
            case 't':
                errno = 0;
                truncate = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'x':
                errno = 0;
                exact_length = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'b'
                        || optopt == 'f'
                        || optopt == 'q'
                        || optopt == 'o'
                        || optopt == 'e'
                        || optopt == 't'
                        || optopt == 'x')
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

    /* select the function */
    if (funcname) {
        function = parasail_lookup_function(funcname);
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

    if (filename_database) {
        sequences_database = parasail_sequences_from_file(filename_database);
        seq_count_database = sequences_database->l;
    }
    else {
        fprintf(stderr, "missing database filename\n");
        exit(1);
    }

    if (filename_queries) {
        sequences_queries = parasail_sequences_from_file(filename_queries);
        seq_count_queries = sequences_queries->l;
        biggest = 0;
        for (i=0; i<seq_count_queries; ++i) {
            if (sequences_queries->seqs[i].seq.l > biggest) {
                biggest = sequences_queries->seqs[i].seq.l;
            }
        }
    }
    else {
        fprintf(stderr, "missing query filename\n");
        exit(1);
    }

    seen = malloc(sizeof(char)*biggest);
    for (i=0; i<biggest; ++i) {
        seen[i] = 0;
    }

    for (i=0; i<seq_count_queries; ++i) {
        int saturated_query = 0;
        double local_timer = timer_real();
        if (truncate > 0 && sequences_queries->seqs[i].seq.l > (unsigned long)truncate) {
            continue;
        }
        if (exact_length > 0 && sequences_queries->seqs[i].seq.l != (unsigned long)exact_length) {
            continue;
        }
        if (seen[sequences_queries->seqs[i].seq.l]) {
            /*printf("skipping %d\n", i);*/
            continue;
        }
        else {
            seen[sequences_queries->seqs[i].seq.l] = 1;
        }
        for (j=0; j<seq_count_database; ++j) {
            parasail_result_t *result = function(
                    sequences_queries->seqs[i].seq.s,
                    sequences_queries->seqs[i].seq.l,
                    sequences_database->seqs[j].seq.s,
                    sequences_database->seqs[j].seq.l,
                    gap_open, gap_extend, matrix);
            saturated_query += parasail_result_is_saturated(result);
            parasail_result_free(result);
        }
        local_timer = timer_real() - local_timer;
        printf("%lu\t %lu\t %d\t %f\n",
                i, (unsigned long)sequences_queries->seqs[i].seq.l,
                saturated_query,
                local_timer);
        if (exact_length != 0) {
            /* if we got this far, we found our query, so break */
            break;
        }
    }

    return 0;
}

