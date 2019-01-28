/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER)
#include <io.h>
#define READ_FUNCTION _read
#else
#include <unistd.h>
#define READ_FUNCTION read
#endif

#include <kseq.h>

#if HAVE_ZLIB
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)
#else
KSEQ_INIT(int, READ_FUNCTION)
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/io.h"
#include "parasail/parasail/stats.h"

static char* strdup_rpl(const char *str)
{
    size_t l = strlen(str);
    char *s = (char*)calloc(l + 1, sizeof(char));
    if(!s) {
        return NULL;
    }
    (void)memcpy(s, str, l);
    return s;
}

parasail_sequences_t* parasail_sequences_from_file(const char *filename)
{
#if HAVE_ZLIB
    gzFile fp;
#else
    FILE* fp;
    int is_stdin = 0;
#endif
    kseq_t *seq = NULL;
    int l = 0;
    parasail_sequences_t *retval = NULL;
    parasail_sequence_t *sequences = NULL;
    unsigned long count = 0;
    unsigned long capacity = 1000;
    size_t chars = 0;
    stats_t stats;

    stats_clear(&stats);

    retval = (parasail_sequences_t*)malloc(sizeof(parasail_sequences_t));
    if (NULL == retval) {
        perror("malloc");
        exit(1);
    }
    sequences = (parasail_sequence_t*)malloc(sizeof(parasail_sequence_t) * capacity);
    if (NULL == sequences) {
        perror("malloc");
        exit(1);
    }

    /* check for stdin instead of a normal filename */
    if (0 == strncmp("stdin", filename, 5)) {
#if HAVE_ZLIB
        fp = gzdopen(fileno(stdin), "r");
        if (fp == Z_NULL) {
            perror("gzdopen");
            exit(1);
        }
#else
        fp = stdin;
        is_stdin = 1;
#endif
    }
    else {
        /* open the file */
        errno = 0;
#if HAVE_ZLIB
        fp = gzopen(filename, "r");
        if (fp == Z_NULL) {
            perror("gzopen");
            exit(1);
        }
#else
        fp = fopen(filename, "r");
        if (fp == NULL) {
            perror("fopen");
            exit(1);
        }
#endif
    }
    
    /* initialize kseq structures */
#if HAVE_ZLIB
    seq = kseq_init(fp);
#else
    seq = kseq_init(fileno(fp));
#endif

    /* parse file */
    while ((l = kseq_read(seq)) >= 0) {
        errno = 0;
        sequences[count].name.l = seq->name.l;
        sequences[count].comment.l = seq->comment.l;
        sequences[count].seq.l = seq->seq.l;
        sequences[count].qual.l = seq->qual.l;
        sequences[count].name.s = NULL;
        sequences[count].comment.s = NULL;
        sequences[count].seq.s = NULL;
        sequences[count].qual.s = NULL;

        stats_sample_value(&stats, seq->seq.l);
        chars += seq->seq.l;

        if (sequences[count].name.l) {
            sequences[count].name.s = strdup_rpl(seq->name.s);
            if (NULL == sequences[count].name.s) {
                perror("strdup name");
                exit(1);
            }
        }
        if (sequences[count].comment.l) {
            sequences[count].comment.s = strdup_rpl(seq->comment.s);
            if (NULL == sequences[count].comment.s) {
                perror("strdup comment");
                exit(1);
            }
        }
        if (sequences[count].seq.l) {
            sequences[count].seq.s = strdup_rpl(seq->seq.s);
            if (NULL == sequences[count].seq.s) {
                perror("strdup seq");
                exit(1);
            }
        }
        if (sequences[count].qual.l) {
            sequences[count].qual.s = strdup_rpl(seq->qual.s);
            if (NULL == sequences[count].qual.s) {
                perror("strdup qual");
                exit(1);
            }
        }
        ++count;

        /* allocate more space for sequences if we ran out */
        if (count >= capacity) {
            parasail_sequence_t *new_sequences = NULL;
            capacity *= 2;
            errno = 0;
            new_sequences = (parasail_sequence_t*)realloc(sequences, sizeof(parasail_sequence_t) * capacity);
            if (NULL == new_sequences) {
                perror("realloc");
                exit(1);
            }
            sequences = new_sequences;
            errno = 0;
        }
    }
    kseq_destroy(seq);
#if HAVE_ZLIB
    gzclose(fp);
#else
    if (!is_stdin) {
        fclose(fp);
    }
#endif

    retval->seqs = sequences;
    retval->l = count;
    retval->characters = chars;
    retval->shortest = stats._min;
    retval->longest = stats._max;
    retval->mean = stats._mean;
    retval->stddev = stats_stddev(&stats);

    return retval;
}

void parasail_sequences_free(parasail_sequences_t *sequences)
{
    size_t i;
    for (i=0; i<sequences->l; ++i) {
        if (sequences->seqs[i].name.s)    free(sequences->seqs[i].name.s);
        if (sequences->seqs[i].comment.s) free(sequences->seqs[i].comment.s);
        if (sequences->seqs[i].seq.s)     free(sequences->seqs[i].seq.s);
        if (sequences->seqs[i].qual.s)    free(sequences->seqs[i].qual.s);
    }
    free(sequences->seqs);
    free(sequences);
}

char* parasail_sequences_pack(const parasail_sequences_t *sequences, size_t *size)
{
    size_t i = 0;
    size_t offset = 0;
    char *packed = NULL;

    packed = (char*)malloc(sizeof(char) * (sequences->characters+sequences->l+1));

    for (i=0; i<sequences->l; ++i) {
        memcpy(&packed[offset], sequences->seqs[i].seq.s, sequences->seqs[i].seq.l);
        offset += sequences->seqs[i].seq.l;
        packed[offset++] = '$';
    }
    packed[offset] = '\0';

    *size = offset;

    return packed;
}

