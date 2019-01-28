/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_IO_H_
#define _PARASAIL_IO_H_

#include <sys/stat.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct parasail_file {
    int fd;
    off_t size;
    const char *buf;
} parasail_file_t;

typedef struct parasail_file_stat {
    unsigned long sequences;
    unsigned long characters;
    unsigned long shortest;
    unsigned long longest;
    float mean;
    float stddev;
} parasail_file_stat_t;

typedef struct parasail_string {
    size_t l;
    char *s;
} parasail_string_t;

typedef struct parasail_sequence {
    parasail_string_t name;
    parasail_string_t comment;
    parasail_string_t seq;
    parasail_string_t qual;
} parasail_sequence_t;

typedef struct parasail_sequences {
    parasail_sequence_t *seqs;
    size_t l;
    size_t characters;
    size_t shortest;
    size_t longest;
    float mean;
    float stddev;
} parasail_sequences_t;

parasail_sequences_t* parasail_sequences_from_file(const char *fname);

char* parasail_sequences_pack(const parasail_sequences_t *sequences, size_t *size);

void parasail_sequences_free(parasail_sequences_t *sequences);

parasail_file_t* parasail_open(const char *fname);

/** Closes file and frees file parameter. */
void parasail_close(parasail_file_t *file);

int parasail_is_fasta(const parasail_file_t *pf);

int parasail_is_fastq(const parasail_file_t *pf);

parasail_file_stat_t* parasail_stat(const parasail_file_t *pf);

parasail_file_stat_t* parasail_stat_fasta(const parasail_file_t *pf);

parasail_file_stat_t* parasail_stat_fastq(const parasail_file_t *pf);

char * parasail_read(const parasail_file_t *pf, long * size);

char * parasail_pack(const parasail_file_t *pf, long * size);

char * parasail_pack_fasta(const parasail_file_t *pf, long * size);

char * parasail_pack_fastq(const parasail_file_t *pf, long * size);


/* char buffer versions of io functions */

int parasail_is_fasta_buffer(const char *, off_t size);

int parasail_is_fastq_buffer(const char *pf, off_t size);

parasail_file_stat_t* parasail_stat_buffer(const char *buffer, off_t size);

parasail_file_stat_t* parasail_stat_fasta_buffer(const char *buffer, off_t size);

parasail_file_stat_t* parasail_stat_fastq_buffer(const char *buffer, off_t size);

char * parasail_pack_buffer(const char *buffer, off_t size, long * packed_size);

char * parasail_pack_fasta_buffer(const char *buffer, off_t size, long * packed_size);

char * parasail_pack_fastq_buffer(const char *buffer, off_t size, long * packed_size);

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_IO_H_ */
