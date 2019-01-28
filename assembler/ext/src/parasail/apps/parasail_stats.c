/**
 * @file parasail_stats
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2015 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Reads fasta/fastq file of database sequences and computes statistics.
 */
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/io.h"

static void print_help(const char *progname, int status) {
    fprintf(stderr, "\nusage: %s "
            "file "
            "\n\n",
            progname);
    exit(status);
}

int main(int argc, char **argv) {
    parasail_sequences_t *sequences = NULL;
    const char *progname = "parasail_stats";

    /* Check arguments. */
    if (argc > 2) {
        fprintf(stderr, "Too many arguments.\n");
        print_help(progname, EXIT_FAILURE);
    }
    else if (argc < 2) {
        fprintf(stderr, "Missing input file.\n");
        print_help(progname, EXIT_FAILURE);
    }

    /* open file */
    sequences = parasail_sequences_from_file(argv[1]);

    /* print the stats */
    fprintf(stdout,
            "%25s: %lu\n"
            "%25s: %lu\n"
            "%25s: %lu\n"
            "%25s: %lu\n"
            "%25s: %f\n"
            "%25s: %f\n",
            "sequence count", (long unsigned)sequences->l,
            "character count", (long unsigned)sequences->characters,
            "shortest sequence", (long unsigned)sequences->shortest,
            "longest sequence", (long unsigned)sequences->longest,
            "sequence length mean", sequences->mean,
            "sequence length stddev", sequences->stddev
            );

    parasail_sequences_free(sequences);

    return 0;
}

