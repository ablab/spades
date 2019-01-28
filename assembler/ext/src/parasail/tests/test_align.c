#include "config.h"

/* getopt needs _POSIX_C_SOURCE 2 */
#define _POSIX_C_SOURCE 2

#include <ctype.h>
#include <limits.h>
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
#include "parasail/parasail/function_lookup.h"
#include "parasail/parasail/io.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/stats.h"
#include "timer.h"
#include "timer_real.h"

static double pctf(double orig, double new)
{
    return orig / new;
}

static void print_array(
        const char * filename_,
        int * array_,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        parasail_result_t *result)
{
    int * array = NULL;
    int i;
    int j;
    FILE *f = NULL;
    const char *filename = filename_;

    if ((result->flag & PARASAIL_FLAG_TRACE)
            && ((result->flag & PARASAIL_FLAG_STRIPED)
                || (result->flag & PARASAIL_FLAG_SCAN))) {
        array = parasail_striped_unwind(s1Len, s2Len, result, array_);
    }
    else {
        array = array_;
    }
    f = fopen(filename, "w");
    if (NULL == f) {
        printf("fopen(\"%s\") error: %s\n", filename, strerror(errno));
        exit(-1);
    }
    fprintf(f, " ");
    for (j=0; j<s2Len; ++j) {
        fprintf(f, "%4c", s2[j]);
    }
    fprintf(f, "\n");
    for (i=0; i<s1Len; ++i) {
        fprintf(f, "%c", s1[i]);
        for (j=0; j<s2Len; ++j) {
            fprintf(f, "%4d", array[i*s2Len + j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    if ((result->flag & PARASAIL_FLAG_TRACE)
            && ((result->flag & PARASAIL_FLAG_STRIPED)
                || (result->flag & PARASAIL_FLAG_SCAN))) {
        free(array);
    }
}

static void print_rowcol(
        const char * filename_,
        const int * const restrict row,
        const int * const restrict col,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len)
{
    int i;
    int j;
    FILE *f = NULL;
    const char *filename = filename_;

    f = fopen(filename, "w");
    if (NULL == f) {
        printf("fopen(\"%s\") error: %s\n", filename, strerror(errno));
        exit(-1);
    }
    fprintf(f, "%c", s1[s1Len-1]);
    if (NULL == row) {
        for (j=0; j<s2Len; ++j) {
            fprintf(f, "%4c", '!');
        }
    }
    else {
        for (j=0; j<s2Len; ++j) {
            fprintf(f, "%4d", row[j]);
        }
    }
    fprintf(f, "\n");
    fprintf(f, "%c", s2[s2Len-1]);
    if (NULL == col) {
        for (i=0; i<s1Len; ++i) {
            fprintf(f, "%4c", '!');
        }
    }
    else {
        for (i=0; i<s1Len; ++i) {
            fprintf(f, "%4d", col[i]);
        }
    }
    fprintf(f, "\n");
    fclose(f);
}


int main(int argc, char **argv)
{
    long seqA_index = LONG_MAX;
    long seqB_index = LONG_MAX;
    const char *seqA = NULL;
    const char *seqB = NULL;
    int lena = 0;
    int lenb = 0;
    int score = 0;
    int matches = 0;
    int similar = 0;
    int length = 0;
    int end_query = 0;
    int end_ref = 0;
    unsigned long long timer_rdtsc = 0;
    unsigned long long timer_rdtsc_single = 0;
    double timer_nsecs = 0.0;
    double timer_nsecs_single = 0.0;
    double timer_nsecs_ref_mean = 0.0;
    double timer_rdtsc_ref_mean = 0.0;
    int limit = 2;
    int i = 0;
    int index = 0;
    parasail_function_info_t f;
    parasail_result_t *result = NULL;
    stats_t stats_rdtsc;
    stats_t stats_nsecs;
    int c = 0;
    char *filename = NULL;
    parasail_sequences_t *sequences = NULL;
    char *endptr = NULL;
    char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    int open = 10;
    int extend = 1;
    int match = 1;
    int mismatch = 0;
    int do_normal = 1;
    int do_stats = 1;
    int do_nonstats = 1;
    int do_table = 1;
    int do_rowcol = 1;
    int do_trace = 1;
    int use_rdtsc = 0;
    int do_sse2 = 1;
    int do_sse41 = 1;
    int do_avx2 = 1;
    int do_disp = 1;
    int do_nw = 1;
    int do_sg = 1;
    int do_sw = 1;
    int use_dna = 0;

    while ((c = getopt(argc, argv, "a:b:df:i:m:M:n:o:e:rRTtNSsBX:")) != -1) {
        switch (c) {
            case 'a':
                errno = 0;
                seqA_index = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol seqA_index");
                    fprintf(stderr, "invalid seqA index\n");
                    exit(1);
                }
                break;
            case 'b':
                errno = 0;
                seqB_index = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol seqB_index");
                    fprintf(stderr, "invalid seqB index\n");
                    exit(1);
                }
                break;
            case 'd':
                use_dna = 1;
                break;
            case 'f':
                filename = optarg;
                break;
            case 'm':
                matrixname = optarg;
                break;
            case 'n':
                errno = 0;
                limit = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol limit");
                    exit(1);
                }
                break;
            case 'o':
                errno = 0;
                open = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol open");
                    exit(1);
                }
                break;
            case 'e':
                errno = 0;
                extend = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol extend");
                    exit(1);
                }
                break;
            case 'M':
                match = atoi(optarg);
                if (match < 0) {
                    fprintf(stderr, "match (-M) must be >= 0");
                    exit(1);
                }
                break;
            case 'X':
                mismatch = atoi(optarg);
                if (mismatch < 0) {
                    fprintf(stderr, "mismatch (-X) must be >= 0");
                    exit(1);
                }
                break;
            case 'r':
                use_rdtsc = 1;
                break;
            case 'R':
                do_rowcol = 0;
                break;
            case 'T':
                do_table = 0;
                break;
            case 't':
                do_trace = 0;
                break;
            case 'N':
                do_normal = 0;
                break;
            case 'S':
                do_stats = 0;
                break;
            case 's':
                do_nonstats = 0;
                break;
            case 'i':
                do_sse2 = (NULL == strstr(optarg, "sse2"));
                do_sse41 = (NULL == strstr(optarg, "sse41"));
                do_avx2 = (NULL == strstr(optarg, "avx2"));
                do_disp = (NULL == strstr(optarg, "disp"));
                do_nw = (NULL == strstr(optarg, "nw"));
                do_sg = (NULL == strstr(optarg, "sg"));
                do_sw = (NULL == strstr(optarg, "sw"));
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'b'
                        || optopt == 'f'
                        || optopt == 'n') {
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
    }
    else {
        fprintf(stderr, "no filename specified\n");
        exit(1);
    }

    /* select the substitution matrix */
    if (NULL == matrixname && use_dna) {
        matrixname = "ACGT";
        matrix = parasail_matrix_create("ACGT", match, -mismatch);
    }
    else {
        if (NULL == matrixname) {
            matrixname = "blosum62";
        }
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            /* try as a filename */
            matrix = parasail_matrix_from_file(matrixname);
        }
        if (NULL == matrix) {
            fprintf(stderr, "Specified substitution matrix not found.\n");
            exit(EXIT_FAILURE);
        }
    }

    if (seqA_index == LONG_MAX) {
        fprintf(stderr, "seqA index not specified\n");
        exit(1);
    }

    if (seqB_index == LONG_MAX) {
        fprintf(stderr, "seqB index not specified\n");
        exit(1);
    }

    if ((unsigned long)seqA_index >= sequences->l) {
        fprintf(stderr, "seqA index out of bounds\n");
        exit(1);
    }

    if ((unsigned long)seqB_index >= sequences->l) {
        fprintf(stderr, "seqB index out of bounds\n");
        exit(1);
    }

    seqA = sequences->seqs[seqA_index].seq.s;
    seqB = sequences->seqs[seqB_index].seq.s;
    lena = sequences->seqs[seqA_index].seq.l;
    lenb = sequences->seqs[seqB_index].seq.l;

    printf("file: %s\n", filename);
    printf("matrix: %s\n", matrixname);
    printf("gap open: %d\n", open);
    printf("gap extend: %d\n", extend);
    printf("seq pair %lu,%lu\n", seqA_index, seqB_index);

    printf("%-15s %8s %6s %4s %5s %5s %4s "
           "%8s %8s %8s %8s %9s %8s "
           "%8s %5s %8s %8s %8s\n",
            "name", "type", "isa", "bits", "width", "elem", "sat",
            "score", "matches", "similar", "length", "end_query", "end_ref",
            "avg", "imp", "stddev", "min", "max");

    stats_clear(&stats_rdtsc);
    stats_clear(&stats_nsecs);
    index = 0;
    f = functions[index++];
    while (f.pointer) {
        char name[16] = {'\0'};
        int new_limit = f.is_table || f.is_rowcol || f.is_trace ? 1 : limit;
        int saturated = 0;
        if ((0 == strncmp(f.isa, "sse2",  4) && 0 == parasail_can_use_sse2()) 
                || (0 == strncmp(f.isa, "sse41", 5) && 0 == parasail_can_use_sse41())
                || (0 == strncmp(f.isa, "avx2",  4) && 0 == parasail_can_use_avx2())) {
            f = functions[index++];
            continue;
        }
        if (!do_sse2 && strstr(f.name, "sse2")) {
            f = functions[index++];
            continue;
        }
        if (!do_sse41 && strstr(f.name, "sse41")) {
            f = functions[index++];
            continue;
        }
        if (!do_avx2 && strstr(f.name, "avx2")) {
            f = functions[index++];
            continue;
        }
        if (!do_disp && strstr(f.name, "disp")) {
            f = functions[index++];
            continue;
        }
        if (!do_nw && strstr(f.name, "nw")) {
            f = functions[index++];
            continue;
        }
        if (!do_sg && strstr(f.name, "sg")) {
            f = functions[index++];
            continue;
        }
        if (!do_sw && strstr(f.name, "sw")) {
            f = functions[index++];
            continue;
        }
        if (f.is_stats) {
            if (!do_stats) {
                f = functions[index++];
                continue;
            }
        }
        else {
            if (!do_nonstats) {
                f = functions[index++];
                continue;
            }
        }
        if (f.is_table) {
            if (!do_table) {
                f = functions[index++];
                continue;
            }
        }
        else if (f.is_rowcol) {
            if (!do_rowcol) {
                f = functions[index++];
                continue;
            }
        }
        else if (f.is_trace) {
            if (!do_trace) {
                f = functions[index++];
                continue;
            }
        }
        else {
            if (!do_normal) {
                f = functions[index++];
                continue;
            }
        }
        fflush(stdout);
        stats_clear(&stats_rdtsc);
        stats_clear(&stats_nsecs);
        timer_rdtsc = timer_start();
        timer_nsecs = timer_real();
        for (i=0; i<new_limit; ++i) {
            timer_rdtsc_single = timer_start();
            timer_nsecs_single = timer_real();
            result = f.pointer(seqA, lena, seqB, lenb, open, extend, matrix);
            timer_rdtsc_single = timer_start()-(timer_rdtsc_single);
            timer_nsecs_single = timer_real() - timer_nsecs_single;
            stats_sample_value(&stats_rdtsc, timer_rdtsc_single);
            stats_sample_value(&stats_nsecs, timer_nsecs_single);
            score = parasail_result_get_score(result);
            if (f.is_stats) {
                similar = parasail_result_get_similar(result);
                matches = parasail_result_get_matches(result);
                length = parasail_result_get_length(result);
            }
            else {
                similar = 0;
                matches = 0;
                length = 0;
            }
            end_query = parasail_result_get_end_query(result);
            end_ref = parasail_result_get_end_ref(result);
            saturated = parasail_result_is_saturated(result);
            parasail_result_free(result);
        }
        timer_rdtsc = timer_start()-(timer_rdtsc);
        timer_nsecs = timer_real() - timer_nsecs;
        if (f.is_ref) {
            timer_nsecs_ref_mean = stats_nsecs._mean;
            timer_rdtsc_ref_mean = stats_rdtsc._mean;
        }
        strcpy(name, f.alg);
        if (f.is_table) {
            strcat(name, "_table");
        }
        else if (f.is_rowcol) {
            strcat(name, "_rowcol");
        }
        else if (f.is_trace) {
            strcat(name, "_trace");
        }
        printf("%-15s %8s %6s %4s %5s %5d ",
                name, f.type, f.isa, f.bits, f.width, f.lanes);
        /* xeon phi was unable to perform I/O running natively */
        if (f.is_table) {
            char suffix[256] = {0};
            if (strlen(f.type)) {
                strcat(suffix, "_");
                strcat(suffix, f.type);
            }
            if (strlen(f.isa)) {
                strcat(suffix, "_");
                strcat(suffix, f.isa);
            }
            if (strlen(f.bits)) {
                strcat(suffix, "_");
                strcat(suffix, f.bits);
            }
            if (strlen(f.width)) {
                strcat(suffix, "_");
                strcat(suffix, f.width);
            }
            strcat(suffix, ".txt");
            result = f.pointer(seqA, lena, seqB, lenb, open, extend, matrix);
            {
                int *table = parasail_result_get_score_table(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_scr");
                strcat(filename, suffix);
                print_array(filename, table, seqA, lena, seqB, lenb, result);
            }
            if (f.is_stats) {
                int *table = parasail_result_get_matches_table(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_mch");
                strcat(filename, suffix);
                print_array(filename, table, seqA, lena, seqB, lenb, result);
            }
            if (f.is_stats) {
                int *table = parasail_result_get_similar_table(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_sim");
                strcat(filename, suffix);
                print_array(filename, table, seqA, lena, seqB, lenb, result);
            }
            if (f.is_stats) {
                int *table = parasail_result_get_length_table(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_len");
                strcat(filename, suffix);
                print_array(filename, table, seqA, lena, seqB, lenb, result);
            }
            parasail_result_free(result);
        }
        else if (f.is_rowcol) {
            char suffix[256] = {0};
            if (strlen(f.type)) {
                strcat(suffix, "_");
                strcat(suffix, f.type);
            }
            if (strlen(f.isa)) {
                strcat(suffix, "_");
                strcat(suffix, f.isa);
            }
            if (strlen(f.bits)) {
                strcat(suffix, "_");
                strcat(suffix, f.bits);
            }
            if (strlen(f.width)) {
                strcat(suffix, "_");
                strcat(suffix, f.width);
            }
            strcat(suffix, ".txt");
            result = f.pointer(seqA, lena, seqB, lenb, open, extend, matrix);
            {
                int *row = parasail_result_get_score_row(result);
                int *col = parasail_result_get_score_col(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_rowcol_scr");
                strcat(filename, suffix);
                print_rowcol(filename, row, col, seqA, lena, seqB, lenb);
            }
            if (f.is_stats) {
                int *row = parasail_result_get_matches_row(result);
                int *col = parasail_result_get_matches_col(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_rowcol_mch");
                strcat(filename, suffix);
                print_rowcol(filename, row, col, seqA, lena, seqB, lenb);
            }
            if (f.is_stats) {
                int *row = parasail_result_get_similar_row(result);
                int *col = parasail_result_get_similar_col(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_rowcol_sim");
                strcat(filename, suffix);
                print_rowcol(filename, row, col, seqA, lena, seqB, lenb);
            }
            if (f.is_stats) {
                int *row = parasail_result_get_length_row(result);
                int *col = parasail_result_get_length_col(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_rowcol_len");
                strcat(filename, suffix);
                print_rowcol(filename, row, col, seqA, lena, seqB, lenb);
            }
            parasail_result_free(result);
        }
        else if (f.is_trace) {
            char suffix[256] = {0};
            if (strlen(f.type)) {
                strcat(suffix, "_");
                strcat(suffix, f.type);
            }
            if (strlen(f.isa)) {
                strcat(suffix, "_");
                strcat(suffix, f.isa);
            }
            if (strlen(f.bits)) {
                strcat(suffix, "_");
                strcat(suffix, f.bits);
            }
            if (strlen(f.width)) {
                strcat(suffix, "_");
                strcat(suffix, f.width);
            }
            strcat(suffix, ".txt");
            result = f.pointer(seqA, lena, seqB, lenb, open, extend, matrix);
            {
                int *table = parasail_result_get_trace_table(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_dag");
                strcat(filename, suffix);
                print_array(filename, table, seqA, lena, seqB, lenb, result);
            }
            {
                int *table = parasail_result_get_trace_ins_table(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_ins");
                strcat(filename, suffix);
                print_array(filename, table, seqA, lena, seqB, lenb, result);
            }
            {
                int *table = parasail_result_get_trace_del_table(result);
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_del");
                strcat(filename, suffix);
                print_array(filename, table, seqA, lena, seqB, lenb, result);
            }
            parasail_result_free(result);
        }
        if (use_rdtsc) {
            printf(
                "%4d "
                "%8d %8d %8d %8d "
                "%9d %8d "
                "%8.1f %5.1f %8.1f %8.0f %8.0f\n",
                saturated,
                score, matches, similar, length,
                end_query, end_ref,
                saturated ? 0 : stats_rdtsc._mean,
                saturated ? 0 : pctf(timer_rdtsc_ref_mean, stats_rdtsc._mean),
                saturated ? 0 : stats_stddev(&stats_rdtsc),
                saturated ? 0 : stats_rdtsc._min,
                saturated ? 0 : stats_rdtsc._max);
        }
        else {
            printf(
                "%4d "
                "%8d %8d %8d %8d "
                "%9d %8d "
                "%8.3f %5.2f %8.3f %8.3f %8.3f\n",
                saturated,
                score, matches, similar, length,
                end_query, end_ref,
                saturated ? 0 : stats_nsecs._mean,
                saturated ? 0 : pctf(timer_nsecs_ref_mean, stats_nsecs._mean),
                saturated ? 0 : stats_stddev(&stats_nsecs),
                saturated ? 0 : stats_nsecs._min,
                saturated ? 0 : stats_nsecs._max);
        }
        f = functions[index++];
    }
    /* banded test */
    if (do_nw) {
        int saturated = 0;
        stats_clear(&stats_rdtsc);
        stats_clear(&stats_nsecs);
        timer_rdtsc = timer_start();
        timer_nsecs = timer_real();
        for (i=0; i<limit; ++i) {
            timer_rdtsc_single = timer_start();
            timer_nsecs_single = timer_real();
            result = parasail_nw_banded(seqA, lena, seqB, lenb, open, extend, 3, matrix);
            timer_rdtsc_single = timer_start()-(timer_rdtsc_single);
            timer_nsecs_single = timer_real() - timer_nsecs_single;
            stats_sample_value(&stats_rdtsc, timer_rdtsc_single);
            stats_sample_value(&stats_nsecs, timer_nsecs_single);
            score = parasail_result_get_score(result);
            similar = 0;
            matches = 0;
            length = 0;
            end_query = parasail_result_get_end_query(result);
            end_ref = parasail_result_get_end_ref(result);
            saturated = parasail_result_is_saturated(result);
            parasail_result_free(result);
        }
        timer_rdtsc = timer_start()-(timer_rdtsc);
        timer_nsecs = timer_real() - timer_nsecs;
        if (use_rdtsc) {
            printf(
                    "%-15s %8s %6s %4s %5s %5d %4d "
                    "%8d %8d %8d %8d "
                    "%9d %8d "
                    "%8.1f %5.1f %8.1f %8.0f %8.0f\n",
                    "nw", "banded", "NA", "32", "32", 1, saturated,
                    score, matches, similar, length,
                    end_query, end_ref,
                    saturated ? 0 : stats_rdtsc._mean,
                    saturated ? 0 : pctf(timer_rdtsc_ref_mean, stats_rdtsc._mean),
                    saturated ? 0 : stats_stddev(&stats_rdtsc),
                    saturated ? 0 : stats_rdtsc._min,
                    saturated ? 0 : stats_rdtsc._max);
        }
        else {
            printf(
                    "%-15s %8s %6s %4s %5s %5d %4d "
                    "%8d %8d %8d %8d "
                    "%9d %8d "
                    "%8.3f %5.2f %8.3f %8.3f %8.3f\n",
                    "nw", "banded", "NA", "32", "32", 1, saturated,
                    score, matches, similar, length,
                    end_query, end_ref,
                    saturated ? 0 : stats_nsecs._mean,
                    saturated ? 0 : pctf(timer_nsecs_ref_mean, stats_nsecs._mean),
                    saturated ? 0 : stats_stddev(&stats_nsecs),
                    saturated ? 0 : stats_nsecs._min,
                    saturated ? 0 : stats_nsecs._max);
        }
    }
    /* ssw test */
    if (do_sw) {
        int saturated = 0;
        stats_clear(&stats_rdtsc);
        stats_clear(&stats_nsecs);
        timer_rdtsc = timer_start();
        timer_nsecs = timer_real();
        for (i=0; i<limit; ++i) {
            parasail_result_ssw_t *ssw_result = NULL;
            timer_rdtsc_single = timer_start();
            timer_nsecs_single = timer_real();
            ssw_result = parasail_ssw(seqA, lena, seqB, lenb, open, extend, matrix);
            timer_rdtsc_single = timer_start()-(timer_rdtsc_single);
            timer_nsecs_single = timer_real() - timer_nsecs_single;
            stats_sample_value(&stats_rdtsc, timer_rdtsc_single);
            stats_sample_value(&stats_nsecs, timer_nsecs_single);
            score = ssw_result->score1;
            similar = 0;
            matches = 0;
            length = 0;
            end_query = ssw_result->read_end1;
            end_ref = ssw_result->ref_end1;
            parasail_result_ssw_free(ssw_result);
        }
        timer_rdtsc = timer_start()-(timer_rdtsc);
        timer_nsecs = timer_real() - timer_nsecs;
        if (use_rdtsc) {
            printf(
                    "%-15s %8s %6s %4s %5s %5d %4d "
                    "%8d %8d %8d %8d "
                    "%9d %8d "
                    "%8.1f %5.1f %8.1f %8.0f %8.0f\n",
                    "ssw", "striped", "NA", "16", "16", 1, saturated,
                    score, matches, similar, length,
                    end_query, end_ref,
                    saturated ? 0 : stats_rdtsc._mean,
                    saturated ? 0 : pctf(timer_rdtsc_ref_mean, stats_rdtsc._mean),
                    saturated ? 0 : stats_stddev(&stats_rdtsc),
                    saturated ? 0 : stats_rdtsc._min,
                    saturated ? 0 : stats_rdtsc._max);
        }
        else {
            printf(
                    "%-15s %8s %6s %4s %5s %5d %4d "
                    "%8d %8d %8d %8d "
                    "%9d %8d "
                    "%8.3f %5.2f %8.3f %8.3f %8.3f\n",
                    "ssw", "striped", "NA", "16", "16", 1, saturated,
                    score, matches, similar, length,
                    end_query, end_ref,
                    saturated ? 0 : stats_nsecs._mean,
                    saturated ? 0 : pctf(timer_nsecs_ref_mean, stats_nsecs._mean),
                    saturated ? 0 : stats_stddev(&stats_nsecs),
                    saturated ? 0 : stats_nsecs._min,
                    saturated ? 0 : stats_nsecs._max);
        }
    }

    if (filename) {
        parasail_sequences_free(sequences);
    }

    return 0;
}
