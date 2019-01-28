/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "parasail/parasail.h"
#include "parasail/parasail/io.h"
#include "parasail/parasail/memory.h"
#include "parasail/parasail/stats.h"

#ifdef UNUSED
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#elif defined(__cplusplus)
# define UNUSED(x)
#else
# define UNUSED(x) x
#endif

/* special functions may exist for windows, msys, mingw */
#if defined(HAVE_STRUCT___STAT64) && defined(HAVE__STAT64) && defined(HAVE__FSTAT64)
#define STATBUF struct __stat64
#define STATFUNC _stat64
#define FSTATFUNC _fstat64
#else
#define STATBUF struct stat
#define STATFUNC stat
#define FSTATFUNC fstat
#endif

parasail_file_t* parasail_open(const char *fname)
{
    parasail_file_t *pf = NULL;
    char *buf = NULL;

#if defined(HAVE_SYS_MMAN_H)
    int fd = -1;
    STATBUF fs;

    fd = open(fname, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Cannot open input file `%s': ", fname);
        perror("open");
        exit(EXIT_FAILURE);
    }

    if (-1 == FSTATFUNC(fd, &fs)) {
        fprintf(stderr, "Cannont stat input file `%s': ", fname);
        perror("fstat");
        exit(EXIT_FAILURE);
    }

    buf = (char*)mmap(NULL, fs.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (MAP_FAILED == buf) {
        fprintf(stderr, "Cannont mmap input file `%s': ", fname);
        perror("mmap");
        exit(EXIT_FAILURE);
    }
#else
    FILE *fd = NULL;
    STATBUF fs;

    fd = fopen(fname, "rb");
    if (NULL == fd) {
        fprintf(stderr, "Cannot open input file `%s': ", fname);
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    if (0 != STATFUNC(fname, &fs)) {
        fprintf(stderr, "Cannont stat input file `%s': ", fname);
        perror("_stat");
        exit(EXIT_FAILURE);
    }

    /* Allocate a buffer to hold the whole file */
    buf = (char*)malloc(fs.st_size + 1);
    if (NULL == buf) {
        fprintf(stderr, "Cannont malloc buffer for input file `%s': ", fname);
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    /* Slurp file into buffer */
    if (fs.st_size != fread(buf, 1, fs.st_size, fd)) {
        fprintf(stderr, "Cannont read input file `%s': ", fname);
        free(buf);
        perror("fread");
        exit(EXIT_FAILURE);
    }
    /* Close the file early */
    fclose(fd);
#endif

    pf = (parasail_file_t*)malloc(sizeof(parasail_file_t));
    if (NULL == pf) {
        fprintf(stderr, "Cannont allocate parasail_file_t");
        perror("malloc");
        exit(EXIT_FAILURE);
    }

#if defined(HAVE_SYS_MMAN_H)
    pf->fd = fd;
#else
    pf->fd = 0;
#endif
    pf->size = fs.st_size;
    pf->buf = buf;
    return pf;
}

void parasail_close(parasail_file_t *pf)
{
#if defined(HAVE_SYS_MMAN_H)
    munmap((void*)pf->buf, pf->size);
    close(pf->fd);
#else
    free((void*)pf->buf);
    /* file was already closed */
#endif
    free(pf);
}

int parasail_is_fasta(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail/parasail_is_fasta given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_is_fasta_buffer(pf->buf, pf->size);
}

int parasail_is_fasta_buffer(const char *buf, off_t UNUSED(size))
{
    if (NULL == buf) {
        fprintf(stderr, "parasail/parasail_is_fasta_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return '>' == buf[0];
}

int parasail_is_fastq(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail/parasail_is_fastq given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_is_fastq_buffer(pf->buf, pf->size);
}

int parasail_is_fastq_buffer(const char *buf, off_t UNUSED(size))
{
    if (NULL == buf) {
        fprintf(stderr, "parasail/parasail_is_fastq_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return '@' == buf[0];
}

parasail_file_stat_t* parasail_stat(const parasail_file_t *pf)
{
    parasail_file_stat_t *stat = NULL;

    if (NULL == pf) {
        fprintf(stderr, "parasail/parasail_stat given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (parasail_is_fasta(pf)) {
        stat = parasail_stat_fasta(pf);
    }
    else if (parasail_is_fastq(pf)) {
        stat = parasail_stat_fastq(pf);
    }
    else {
        fprintf(stderr, "parasail/parasail_stat: cannot determine file format\n");
        exit(EXIT_FAILURE);
    }

    return stat;
}

parasail_file_stat_t* parasail_stat_buffer(const char *buf, off_t size)
{
    parasail_file_stat_t *stat = NULL;

    if (NULL == buf) {
        fprintf(stderr, "parasail/parasail_stat_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (parasail_is_fasta_buffer(buf, size)) {
        stat = parasail_stat_fasta_buffer(buf, size);
    }
    else if (parasail_is_fastq_buffer(buf, size)) {
        stat = parasail_stat_fastq_buffer(buf, size);
    }
    else {
        fprintf(stderr, "parasail/parasail_stat: cannot determine file format\n");
        exit(EXIT_FAILURE);
    }

    return stat;
}

/* increments i until T[i] points to final newline character, returns index */
inline static off_t skip_line(const char *T, off_t i)
{
    while (T[i] != '\n' && T[i] != '\r') {
        ++i;
    }

    /* for the case of "\r\n" or "\n\r" */
    if (T[i+1] == '\n' || T[i+1] == '\r') {
        ++i;
    }

    return i;
}

parasail_file_stat_t* parasail_stat_fasta(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail/parasail_stat_fasta given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_stat_fasta_buffer(pf->buf, pf->size);
}

parasail_file_stat_t* parasail_stat_fasta_buffer(const char *T, off_t size)
{
    off_t i = 0;
    unsigned long seq = 0;
    unsigned long c = 0;
    unsigned long c_tot = 0;
    stats_t stats;
    parasail_file_stat_t *pfs = NULL;

    stats_clear(&stats);

    if (NULL == T) {
        fprintf(stderr, "parasail/parasail_stat_fasta_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    /* first line is always first sequence ID */
    if (T[i] != '>') {
        fprintf(stderr, "poorly formatted FASTA file\n");
        exit(EXIT_FAILURE);
    }

    i = skip_line(T, i);
    ++i;

    /* count that first sequence */
    ++seq;

    /* read rest of file */
    while (i<size) {
        if (T[i] == '>') {
            /* encountered a new sequence */
            ++seq;
            stats_sample_value(&stats, c);
            c = 0;
            i = skip_line(T, i);
        }
        else if (isalpha(T[i])) {
            ++c;
            ++c_tot;
        }
        else if (T[i] == '\n' || T[i] == '\r') {
            /* ignore newline */
            /* for the case of "\r\n" or "\n\r" */
            if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
                ++i;
            }
        }
        else if (isprint(T[i])) {
            fprintf(stderr, "error: non-alpha character ('%c')\n", T[i]);
            exit(EXIT_FAILURE);
        }
        else {
            fprintf(stderr, "error: non-printing character ('%d')\n", T[i]);
            exit(EXIT_FAILURE);
        }
        ++i;
    }

    /* still should have one sequence in the pipe */
    if (0 == c) {
        fprintf(stderr, "error: empty sequence at end of input\n");
        exit(EXIT_FAILURE);
    }
    stats_sample_value(&stats, c);

    pfs = (parasail_file_stat_t*)malloc(sizeof(parasail_file_stat_t));
    if (NULL == pfs) {
        fprintf(stderr, "Cannont allocate parasail_file_stat_t");
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    pfs->sequences = seq;
    pfs->characters = c_tot;
    pfs->shortest = (unsigned long)stats._min;
    pfs->longest = (unsigned long)stats._max;
    pfs->mean = (float)stats._mean;
    pfs->stddev = (float)stats_stddev(&stats);

    return pfs;
}

/*
 * Line 1 begins with a '@' character and is followed by a sequence
 * identifier and an optional description (like a FASTA title line).
 *
 * Line 2 is the raw sequence letters.
 *
 * Line 3 begins with a '+' character and is optionally followed by the
 * same sequence identifier (and any description) again.
 *
 * Line 4 encodes the quality values for the sequence in Line 2, and
 * must contain the same number of symbols as letters in the sequence.
 */
parasail_file_stat_t* parasail_stat_fastq(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail/parasail_stat_fastq given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_stat_fastq_buffer(pf->buf, pf->size);
}

parasail_file_stat_t* parasail_stat_fastq_buffer(const char *T, off_t size)
{
    int first = 1;
    off_t i = 0;
    unsigned long seq = 0;
    unsigned long c = 0;
    unsigned long c_tot = 0;
    unsigned long line = 0;
    stats_t stats;
    parasail_file_stat_t *pfs = NULL;

    stats_clear(&stats);

    if (NULL == T) {
        fprintf(stderr, "parasail/parasail_stat_fastq_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    /* read file */
    while (i<size) {
        if (T[i] != '@') {
            fprintf(stderr, "poorly formatted FASTQ file\n");
            fprintf(stderr, "line %lu\n", line);
            exit(EXIT_FAILURE);
        }

        /* encountered a new sequence */
        ++seq;
        if (first) {
            first = 0;
        }
        else {
            stats_sample_value(&stats, c);
        }
        c = 0;

        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;

        /* rest of next line is the sequence */
        while (T[i] != '\n' && T[i] != '\r') {
            ++c;
            ++i;
        }

        /* for the case of "\r\n" or "\n\r" */
        if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
            ++i;
        }

        stats_sample_value(&stats, c);

        /* go to next line */
        ++i;
        ++line;

        if (T[i] != '+') {
            fprintf(stderr, "poorly formatted FASTQ file\n");
            fprintf(stderr, "line %lu\n", line);
            exit(EXIT_FAILURE);
        }

        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;

        /* rest of next line are the quality control values */
        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;
    }

    pfs = (parasail_file_stat_t*)malloc(sizeof(parasail_file_stat_t));
    if (NULL == pfs) {
        fprintf(stderr, "Cannont allocate parasail_file_stat_t");
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    pfs->sequences = seq;
    pfs->characters = c_tot;
    pfs->shortest = (unsigned long)stats._min;
    pfs->longest = (unsigned long)stats._max;
    pfs->mean = (float)stats._mean;
    pfs->stddev = (float)stats_stddev(&stats);

    return pfs;
}

char * parasail_read(const parasail_file_t *pf, long * size)
{
    char * buffer = (char*)malloc(sizeof(char) * (pf->size+1));
    if (NULL == buffer) {
        fprintf(stderr, "Cannont malloc buffer for input file");
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    (void)memcpy(buffer, pf->buf, pf->size);
    buffer[pf->size] = '\0';
    *size = pf->size;
    return buffer;
}

char * parasail_pack(const parasail_file_t *pf, long * size)
{
    char *packed = NULL;

    if (NULL == pf) {
        fprintf(stderr, "parasail/parasail_pack given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (parasail_is_fasta(pf)) {
        packed = parasail_pack_fasta(pf, size);
    }
    else if (parasail_is_fastq(pf)) {
        packed = parasail_pack_fastq(pf, size);
    }
    else {
        fprintf(stderr, "parasail/parasail_pack: cannot determine file format\n");
        exit(EXIT_FAILURE);
    }

    return packed;
}

char * parasail_pack_buffer(const char *buf, off_t size, long * packed_size)
{
    char *packed = NULL;

    if (NULL == buf) {
        fprintf(stderr, "parasail/parasail_pack_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (parasail_is_fasta_buffer(buf, size)) {
        packed = parasail_pack_fasta_buffer(buf, size, packed_size);
    }
    else if (parasail_is_fastq_buffer(buf, size)) {
        packed = parasail_pack_fastq_buffer(buf, size, packed_size);
    }
    else {
        fprintf(stderr, "parasail/parasail_pack: cannot determine file format\n");
        exit(EXIT_FAILURE);
    }

    return packed;
}

char * parasail_pack_fasta(const parasail_file_t *pf, long * packed_size)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail/parasail_pack_fasta given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail/parasail_pack_fasta given NULL size pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_pack_fasta_buffer(pf->buf, pf->size, packed_size);
}

char * parasail_pack_fasta_buffer(const char *T, off_t size, long * packed_size)
{
    parasail_file_stat_t *pfs = NULL;
    off_t i = 0;
    off_t w = 0;
    char *P = NULL;

    if (NULL == T) {
        fprintf(stderr, "parasail/parasail_pack_fasta_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail/parasail_pack_fasta_buffer given NULL size pointer\n");
        exit(EXIT_FAILURE);
    }

    pfs = parasail_stat_fasta_buffer(T, size);

    P = (char*)malloc(sizeof(char) * (pfs->characters+pfs->sequences+1));

    /* first line is always first sequence ID */
    if (T[i] != '>') {
        fprintf(stderr, "poorly formatted FASTA file\n");
        exit(EXIT_FAILURE);
    }

    i = skip_line(T, i);
    ++i;

    /* read rest of file */
    while (i<size) {
        if (T[i] == '>') {
            /* encountered a new sequence */
            P[w++] = '$';
            i = skip_line(T, i);
        }
        else if (isalpha(T[i])) {
            P[w++] = T[i];
        }
        else if (T[i] == '\n' || T[i] == '\r') {
            /* ignore newline */
            /* for the case of "\r\n" or "\n\r" */
            if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
                ++i;
            }
        }
        else if (isprint(T[i])) {
            fprintf(stderr, "error: non-alpha character ('%c')\n", T[i]);
            exit(EXIT_FAILURE);
        }
        else {
            fprintf(stderr, "error: non-printing character ('%d')\n", T[i]);
            exit(EXIT_FAILURE);
        }
        ++i;
    }

    free(pfs);

    P[w++] = '$';
    P[w] = '\0';
    *packed_size = w;
    return P;
}

char * parasail_pack_fastq(const parasail_file_t *pf, long * size)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail/parasail_pack_fastq given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (NULL == size) {
        fprintf(stderr, "parasail/parasail_pack_fastq given NULL size pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_pack_fastq_buffer(pf->buf, pf->size, size);
}

char * parasail_pack_fastq_buffer(const char *T, off_t size, long * packed_size)
{
    char *P = NULL;
    int first = 1;
    off_t i = 0;
    off_t w = 0;
    unsigned long line = 0;
    parasail_file_stat_t *pfs = NULL;

    if (NULL == T) {
        fprintf(stderr, "parasail/parasail_pack_fastq_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail/parasail_pack_fastq_buffer given NULL size pointer\n");
        exit(EXIT_FAILURE);
    }

    pfs = parasail_stat_fastq_buffer(T, size);

    P = (char*)malloc(sizeof(char) * (pfs->characters+pfs->sequences+1));

    /* read file */
    while (i<size) {
        if (T[i] != '@') {
            fprintf(stderr, "poorly formatted FASTQ file\n");
            fprintf(stderr, "line %lu\n", line);
            exit(EXIT_FAILURE);
        }

        /* encountered a new sequence */
        if (first) {
            first = 0;
        }
        else {
            P[w++] = '$';
        }

        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;

        /* rest of next line is the sequence */
        while (T[i] != '\n' && T[i] != '\r') {
            P[w++] = T[i];
            ++i;
        }

        /* for the case of "\r\n" or "\n\r" */
        if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
            ++i;
        }

        /* go to next line */
        ++i;
        ++line;

        if (T[i] != '+') {
            fprintf(stderr, "poorly formatted FASTQ file\n");
            fprintf(stderr, "line %lu\n", line);
            exit(EXIT_FAILURE);
        }

        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;

        /* rest of next line are the quality control values */
        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;
    }

    free(pfs);

    P[w++] = '$';
    P[w] = '\0';
    *packed_size = w;
    return P;
}

/* increments i until T[i] points non-number, returns number */
#define TOKEN_MAX 10
inline static off_t get_num(const char *T, off_t i, int *result)
{
    int retval = 0;
    int p = 0;
    char token[TOKEN_MAX];

    /* skip any whitespace */
    while (T[i] == ' ' || T[i] == '\t') {
        ++i;
    }

    if (isdigit(T[i]) || T[i] == '-') {
        token[0] = T[i];
    }
    else {
        fprintf(stderr, "poorly formed matrix file\n");
        exit(EXIT_FAILURE);
    }

    ++i;
    for (p=1; p<TOKEN_MAX; ++p) {
        if (isdigit(T[i])) {
            token[p] = T[i];
            ++i;
        }
        else {
            break;
        }
    }
    if (TOKEN_MAX == p) {
        fprintf(stderr, "poorly formed matrix file\n");
        exit(EXIT_FAILURE);
    }
    token[p] = '\0';

    retval = sscanf(token, "%d", result);
    if (1 != retval) {
        fprintf(stderr, "poorly formed matrix file\n");
        exit(EXIT_FAILURE);
    }

    return i;
}

#define ALPHABET_MAX 256
inline static char*  get_alphabet(const char *T, off_t i, off_t size)
{
    char *alphabet = NULL;
    off_t _i = i;
    size_t count = 0;

    /* count number of letters first */
    while (i<size) {
        if (T[i] == '\n' || T[i] == '\r' || T[i] == '#') {
            break;
        }
        else if (T[i] == ' ' || T[i] == '\t') {
            /* ignore newline */
        }
        else if (isalpha(T[i]) || T[i] == '*') {
            ++count;
        }
        else {
            fprintf(stderr, "error: poorly formed matrix file alphabet\n");
            exit(EXIT_FAILURE);
        }
        ++i;
    }

    if (0 == count) {
        fprintf(stderr, "error: poorly formed matrix file alphabet\n");
        exit(EXIT_FAILURE);
    }

    alphabet = (char*)malloc(sizeof(char)*(count+1));
    if (NULL == alphabet) {
        fprintf(stderr, "Cannont malloc buffer for matrix alphabet\n");
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    i = _i;
    count = 0;
    while (i<size) {
        if (T[i] == '\n' || T[i] == '\r' || T[i] == '#') {
            break;
        }
        else if (T[i] == ' ' || T[i] == '\t') {
            /* ignore newline */
        }
        else if (isalpha(T[i]) || T[i] == '*') {
            alphabet[count] = T[i];
            ++count;
        }
        else {
            fprintf(stderr, "error: poorly formed matrix file alphabet\n");
            exit(EXIT_FAILURE);
        }
        ++i;
    }

    alphabet[count] = '\0';
    return alphabet;
}

/* expects a matrix in the form
 * # FREQS A 0.325 C 0.175 G 0.175 T 0.325
 *     A   R   G   C   Y   T   K   M   S   W   N   X
 * A   8   3  -7 -17 -19 -21 -14  -4 -12  -6  -1 -30
 * R   0   3   2 -16 -18 -19  -8  -8  -7 -10  -1 -30
 * G -10  12  12 -16 -17 -18  -2 -13  -1 -14  -1 -30
 * C -18 -17 -16  12  12 -10 -13  -2  -1 -14  -1 -30
 * Y -19 -18 -16   2   0   0  -8  -8  -7 -10  -1 -30
 * T -21 -19 -17  -7   3   8  -4 -14 -12  -6  -1 -30
 * K -15  -9  -2 -11  -8  -4  -3 -13  -7 -10  -1 -30
 * M  -4  -8 -11  -2  -9 -15 -13  -3  -7 -10  -1 -30
 * S -14  -8  -1  -1  -8 -14  -8  -8  -1 -14  -1 -30
 * W  -6  -9 -12 -12  -9  -6  -9  -9 -12  -6  -1 -30
 * N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
 * X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
 */
parasail_matrix_t* parasail_matrix_from_file(const char *filename)
{
    parasail_matrix_t *retval = NULL;
    int *matrix = NULL;
    int *mapper = NULL;
    char *alphabet = NULL;
    parasail_file_t *pf = NULL;
    const char *T = NULL;
    off_t i = 0;
    off_t size = 0;
    int first_alpha = 1;
    size_t count = 0;
    size_t asize = 0;
    int max = INT_MIN;
    int min = INT_MAX;
    size_t c = 0;

    mapper = (int*)malloc(sizeof(int)*256);

    pf = parasail_open(filename);
    T = pf->buf;
    size = pf->size;

    while (i<size) {
        if (T[i] == '#') {
            /* ignore comments */
            i = skip_line(T, i);
        }
        else if (isalpha(T[i]) || T[i] == '*') {
            if (first_alpha) {
                first_alpha = 0;
                alphabet = get_alphabet(T, i, size);
                asize = strlen(alphabet);
                matrix = (int*)malloc(sizeof(int)*asize*asize);
                if (NULL == matrix) {
                    fprintf(stderr, "Cannont malloc buffer for matrix");
                    perror("malloc");
                    exit(EXIT_FAILURE);
                }
                i = skip_line(T, i);
            }
            else {
                size_t j=0;
                /* make sure it is in same order as first line */
                if (T[i] != alphabet[count]) {
                    fprintf(stderr, "error: matrix header out of order\n");
                    exit(EXIT_FAILURE);
                }
                ++i; /* skip over the letter */
                ++count;
                for (j=0; j<asize; ++j) {
                    int val = 0;
                    i = get_num(T, i, &val);
                    matrix[c++] = val;
                    max = val > max ? val : max;
                    min = val < min ? val : min;
                }
            }
        }
        else if (T[i] == '\n' || T[i] == '\r') {
            /* ignore newline */
            /* for the case of "\r\n" or "\n\r" */
            if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
                ++i;
            }
        }
        else if (T[i] == ' ' || T[i] == '\t') {
            /* ignore spaces and tabs */
        }
        else if (isprint(T[i])) {
            fprintf(stderr, "error: non-alpha character in matrix file ('%c')\n", T[i]);
            exit(EXIT_FAILURE);
        }
        else {
            fprintf(stderr, "error: non-printing character in matrix file ('%d')\n", T[i]);
            exit(EXIT_FAILURE);
        }
        ++i;
    }

    parasail_close(pf);

    if (c != asize*asize) {
        fprintf(stderr, "matrix is missing values");
        exit(EXIT_FAILURE);
    }
    if (count != asize) {
        fprintf(stderr, "matrix is missing rows");
        exit(EXIT_FAILURE);
    }

    mapper = (int*)malloc(sizeof(int)*256);
    if (NULL == mapper) {
        fprintf(stderr, "Cannont malloc buffer for matrix mapper for file `%s': ", filename);
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    parasail_memset_int(mapper, size-1, 256);
    for (c=0; c<asize; ++c) {
        mapper[toupper((unsigned char)alphabet[c])] = (int)c;
        mapper[tolower((unsigned char)alphabet[c])] = (int)c;
    }

    free(alphabet);

    retval = (parasail_matrix_t*)malloc(sizeof(parasail_matrix_t));
    if (NULL == retval) {
        fprintf(stderr, "Cannont malloc buffer for matrix file `%s': ", filename);
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    retval->name = filename;
    retval->matrix = matrix;
    retval->mapper = mapper;
    retval->size = (int)asize;
    retval->max = max;
    retval->min = min;
    retval->user_matrix = matrix;
    return retval;
}

