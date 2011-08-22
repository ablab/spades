/*
 * Copyright (c) Medical Research Council 1994. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by James Bonfield, Simon Dear, Rodger Staden,
 * as part of the Staden Package at the MRC Laboratory of Molecular
 * Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 */

/*
 * File: scf.h
 * Version: 3.00
 *
 * Description: file structure definitions for SCF file
 *
 * Created: 19 November 1992
 *
 */

#ifndef _SCF_H_
#define _SCF_H_

#include <stdio.h>
#include <sys/types.h>

#include "io_lib/mFILE.h"
#include "io_lib/os.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 *-----------------------------------------------------------------------------
 * Macros
 *-----------------------------------------------------------------------------
 */

/* The SCF magic number */
#define SCF_MAGIC (((((uint_4)'.'<<8)+(uint_4)'s'<<8)+(uint_4)'c'<<8)+(uint_4)'f')

/* prior to this was a different format */
#define SCF_VERSION_OLDEST 2.00
#define SCF_VERSION_OLD 2.02

/* The current SCF format level */
#define SCF_VERSION 3.00

/* Uncertainty code sets supported */
#define CSET_DEFAULT 0  /* {A,C,G,T,-} */
#define CSET_STADEN  1 
#define CSET_NC_IUB  2  /* Pharmacia A.L.F. */
#define CSET_ALF     3  /* extended NC_IUB */
#define CSET_ABI     4  /* {A,C,G,T,N} */
#define CSET_IBI     5  /* IBI/Pustell */
#define CSET_DNASTAR 6  /* DNA* */
#define CSET_DNASIS  7
#define CSET_PCGENE  8  /* IG/PC-Gene */
#define CSET_GENIE   9  /* MicroGenie */

/* define samples to delta_delta values */
#define DELTA_IT 1

/* What components to read */
#define READ_BASES	(1<<0)
#define READ_SAMPLES	(1<<1)
#define READ_COMMENTS	(1<<2)
#define READ_ALL	(READ_BASES | READ_SAMPLES | READ_COMMENTS)

/*
 *-----------------------------------------------------------------------------
 * Structures and typedefs
 *-----------------------------------------------------------------------------
 */

/*
 * Type definition for the Header structure
 */
typedef struct {
    uint_4 magic_number;       /* SCF_MAGIC */
    uint_4 samples;            /* Number of elements in Samples matrix */
    uint_4 samples_offset;     /* Byte offset from start of file */
    uint_4 bases;              /* Number of bases in Bases matrix */
    uint_4 bases_left_clip;    /* OBSOLETE: No. bases in left clip (vector) */
    uint_4 bases_right_clip;   /* OBSOLETE: No. bases in right clip (qual) */
    uint_4 bases_offset;       /* Byte offset from start of file */
    uint_4 comments_size;      /* Number of bytes in Comment section */
    uint_4 comments_offset;    /* Byte offset from start of file */
    char   version[4];	       /* "version.revision" */
    uint_4 sample_size;	       /* precision of samples (in bytes) */
    uint_4 code_set;	       /* uncertainty codes used */
    uint_4 private_size;       /* size of private data, 0 if none */
    uint_4 private_offset;     /* Byte offset from start of file */
    uint_4 spare[18];          /* Unused */
} Header;

/*
 * Header.sample_size == 1.
 */
typedef struct {
    uint_1 sample_A;			/* Sample for A trace */
    uint_1 sample_C;			/* Sample for C trace */
    uint_1 sample_G;			/* Sample for G trace */
    uint_1 sample_T;			/* Sample for T trace */
} Samples1;

/*
 * Header.sample_size == 2.
 */
typedef struct {
    uint_2 sample_A;			/* Sample for A trace */
    uint_2 sample_C;			/* Sample for C trace */
    uint_2 sample_G;			/* Sample for G trace */
    uint_2 sample_T;			/* Sample for T trace */
} Samples2;

/*
 * Type definition for the sequence data
 */
typedef struct {
    uint_4 peak_index;        /* Index into Samples matrix for base position */
    uint_1 prob_A;            /* Probability of it being an A */
    uint_1 prob_C;            /* Probability of it being an C */
    uint_1 prob_G;            /* Probability of it being an G */
    uint_1 prob_T;            /* Probability of it being an T */
    char base;		      /* Base called */
    uint_1 spare[3];          /* Spare */
} Bases;


/*
 * Type definition for the comments
 */
typedef char Comments;      /* Zero terminated list of \n separated entries */


/*
 * All of the above structs in a single scf format.
 */
typedef struct {
    Header header;
    union Samples {
	Samples1 *samples1;
	Samples2 *samples2;
    } samples;
    Bases *bases;
    Comments *comments;
    char *private_data;
} Scf;

/*
 *-----------------------------------------------------------------------------
 * Function prototypes
 *-----------------------------------------------------------------------------
 */

/*
 * Reading SCF Files
 * -----------------
 */

/*
 * Read the Header struct.
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_header(mFILE *fp, Header *h);

/*
 * Read a single 8bit sample
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_sample1(mFILE *fp, Samples1 *s);

/*
 * Read several 8bit samples
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_samples1(mFILE *fp, Samples1 *s, size_t num_samples);

/*
 * Read several 8bit samples in delta_delta format
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_samples31(mFILE *fp, Samples1 *s, size_t num_samples);

/*
 * Read a single 16bit sample
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_sample2(mFILE *fp, Samples2 *s);

/*
 * Read several 16bit samples
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_samples2(mFILE *fp, Samples2 *s, size_t num_samples);

/*
 * Read several 16bit samples in delta_delta format
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_samples32(mFILE *fp, Samples2 *s, size_t num_samples);

/*
 * Read a single Bases structure
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_base(mFILE *fp, Bases *b);

/*
 * Read several Bases structures consecutively
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_bases(mFILE *fp, Bases *b, size_t num_bases);

/*
 * Read Bases, peak_indexes and probs
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_bases3(mFILE *fp, Bases *b, size_t num_bases);

/*
 * Read the SCF Comments.
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int read_scf_comment(mFILE *fp, Comments *c, size_t l);

/*
 * Reads a whole SCF file into a Scf structure. This memory for this
 * structure is allocated by this routine. To free this memory use
 * scf_deallocate().
 * Returns:
 *    Scf *	- Success, the Scf structure read.
 *    NULL	- Failure.
 * On failure NULL is returned, otherwise the Scf struct.
 */
Scf *read_scf(char *fn);
Scf *fread_scf(FILE *fp);
Scf *mfread_scf(mFILE *fp);


/*
 * Writing SCF Files
 * -----------------
 */

/*
 * Write the Header struct.
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_header(mFILE *fp, Header *h);

/*
 * Write a single 8bit sample
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_sample1(mFILE *fp, Samples1 *s);

/*
 * Write several 8bit samples
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_samples1(mFILE *fp, Samples1 *s, size_t num_samples);

/*
 * Write several 8bit samples in delta_delta format
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_samples31(mFILE *fp, Samples1 *s, size_t num_samples);

/*
 * Write 16bit samples
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_sample2(mFILE *fp, Samples2 *s);

/*
 * Write several 16bit samples
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_samples2(mFILE *fp, Samples2 *s, size_t num_samples);

/*
 * Write several 16bit samples in delta_delta format
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_samples32(mFILE *fp, Samples2 *s, size_t num_samples);

/*
 * Write the Bases structure
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_base(mFILE *fp, Bases *b);

/*
 * Write the several Bases structures consecutively
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_bases(mFILE *fp, Bases *b, size_t num_bases);

/*
 * Write the bases, then peak indexes, then probs
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_bases3(mFILE *fp, Bases *b, size_t num_bases);

/*
 * Write the SCF Comments.
 * Returns:
 *    0 - success
 *   -1 - failure
 */
int write_scf_comment(mFILE *fp, Comments *c, size_t l);


/*
 * Writes a whole Scf structure to filename "fn".
 * This initialises several fields in the Header struct for you. These are:
 *     samples_offset
 *     bases_offset
 *     comments_offset
 *     magic_number
 *
 * All other fields are assumed to be correctly set.
 *
 * Returns:
 *     0 for success
 *    -1 for failure
 */
int write_scf(Scf *scf, char *fn);
int fwrite_scf(Scf *scf, FILE *fp);
int mfwrite_scf(Scf *scf, mFILE *fp);

/*
 * Request which (major) version of scf to use when writing.
 * Defaults to the latest. Currently suitable fields are
 * 2 and 3.
 *
 * Returns 0 for success, -1 for failure.
 */
int set_scf_version(int version);


/*
 * Miscellaneous SCF utilities
 * ---------------------------
 */

/*
 * Converts an SCF version string (eg "2.00") to a float
 */
float scf_version_str2float(char version[]);

/*
 * Converts an SCF version float (eg 2.00) to a string
 * Returns:
 *    A statically allocated 5 character string.
 */
char *scf_version_float2str(float f);

/*
 * Allocates memory for the scf elements based upon arguments passed.
 * Returns;
 *    Scf *	- Success. The scf structure and it's samples, bases,
 *                and comments fields have been allocated.
 *    NULL	- Failure.
 */
Scf *scf_allocate(int num_samples, int sample_size, int num_bases,
		  int comment_size, int private_size);

/*
 * Frees memory allocated by scf_allocate.
 */
void scf_deallocate(Scf *scf);

/*
 * Checks to see if the file with name "fn" is in SCF format.
 * Returns:
 *   1  - is in SCF format
 *   0  - is not in SCF format
 *  -1  - failure
 */
int is_scf(char *fn);

/*
 * Change sample points to delta_delta values for uint1
 */
void scf_delta_samples1 ( int1 samples[], int num_samples, int job);

/*
 * Change sample points to delta_delta values for uint2
 */
void scf_delta_samples2 ( uint2 samples[], int num_samples, int job);

#ifdef __cplusplus
}
#endif

#endif /*_SCF_H_*/

