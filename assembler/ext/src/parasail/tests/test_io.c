#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "parasail/parasail/io.h"

int main(int argc, char **argv)
{
    parasail_sequences_t *sequences = NULL;
    char *packed = NULL;
    unsigned long packed_size = 0;

    if (2 != argc) {
        fprintf(stderr, "only takes a single argument, a fasta/fastq file\n");
        exit(EXIT_FAILURE);
    }

    sequences = parasail_sequences_from_file(argv[1]);
    printf("--sequences--\n");
    printf("count = %lu\n", (long unsigned)sequences->l);
    printf("characters = %lu\n", (long unsigned)sequences->characters);
    printf("shortest = %lu\n", (long unsigned)sequences->shortest);
    printf("longest = %lu\n", (long unsigned)sequences->longest);
    printf("mean = %lu\n", (long unsigned)sequences->mean);
    printf("stddev = %lu\n", (long unsigned)sequences->stddev);

    packed = parasail_sequences_pack(sequences, &packed_size);
    printf("\n--packed--\n");
    printf("size = %lu\n", packed_size);
    printf("%s\n", packed);

    free(packed);
    parasail_sequences_free(sequences);

    return 0;
}

