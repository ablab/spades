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

#ifndef _Read_h_
#define _Read_h_

/* 
 * Title:	Read
 *
 * File: 	Read.h
 * Purpose:	Read data type
 * Last update:	June  14 1994
 */

/*
 * This module encodes the `Read' sequence data structure.
 *
 * A `Read' contains information about bases and traces which are laid
 * out along a single dimension of points. The number of points in a
 * paricular sequence is given by `getNPoints', and these are numbered
 * 0..getNPoints-1. At each point there are four trace readings, one
 * for each base.
 *
 * The number of bases is `getNBases' which are numbered 0..N-1. 
 * Bases are represented by `char's. Every base is located at a 
 * particular point.
 *
 * The behaviour of these routines is undefined if given NULLRead or
 * an undefined sequence.
 */

#include "io_lib/os.h"
#include "io_lib/scf.h"
#include "io_lib/mFILE.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 *-----------------------------------------------------------------------------
 * Macros
 *-----------------------------------------------------------------------------
 */

#define NULLRead     ((Read *)NULL)

/* Trace file formats */
#define TT_ERR -1
#define TT_UNK 0
#define TT_SCF 1
#define TT_ABI 2
#define TT_ALF 3
#define TT_PLN 4
#define TT_EXP 5
#define TT_CTF 6
#define TT_ZTR 7
#define TT_ZTR1 8
#define TT_ZTR2 9
#define TT_ZTR3 10
#define TT_BIO 11
#define TT_SFF 12
#define TT_ANY TT_UNK
/* ANYTR is specifically any *trace* type and not EXP or PLN format */
#define TT_ANYTR 13

#define READ_BASES	(1<<0)
#define READ_SAMPLES	(1<<1)
#define READ_COMMENTS	(1<<2)
#define READ_ALL	(READ_BASES | READ_SAMPLES | READ_COMMENTS)

/*
 *-----------------------------------------------------------------------------
 * Structures and typedefs
 *-----------------------------------------------------------------------------
 */

typedef uint_2 TRACE;        /* for trace heights */

typedef struct
{
    int		format;	     /* Trace file format */
    char       *trace_name;  /* Trace file name */

    int         NPoints;     /* No. of points of data */
    int         NBases;      /* No. of bases */

    /* Traces */
    TRACE      *traceA;      /* Array of length `NPoints' */
    TRACE      *traceC;      /* Array of length `NPoints' */
    TRACE      *traceG;      /* Array of length `NPoints' */
    TRACE      *traceT;      /* Array of length `NPoints' */
    TRACE       maxTraceVal; /* The maximal value in any trace */
    int         baseline;    /* The zero offset for TRACE values */ 

    /* Bases */
    char       *base;        /* Array of length `NBases' */
    uint_2     *basePos;     /* Array of length `NBases' */

    /* Cutoffs */
    int         leftCutoff;  /* Number of unwanted bases */
    int         rightCutoff; /* First unwanted base at right end */

    /* Miscellaneous Sequence Information */
    char       *info;        /* misc seq info, eg comments */

    /* Probability information */
    char       *prob_A;      /* Array of length 'NBases' */
    char       *prob_C;      /* Array of length 'NBases' */
    char       *prob_G;      /* Array of length 'NBases' */
    char       *prob_T;      /* Array of length 'NBases' */

    /* The original input format data, or NULL if inapplicable */
    int orig_trace_format;
    void (*orig_trace_free)(void *ptr);
    void *orig_trace;

    char       *ident;	     /* Seq id, NULL for unknown. Malloced data.
				Owned and freed by io_lib. */

    /* Pyrosequencing "peaks" (more like spikes). NULL if not used */
    int          nflows;     /* Number of "flows" */
    char        *flow_order; /* Bases flowed across */
    float       *flow;       /* Processed to be 1 base unit oriented */
    unsigned int*flow_raw;   /* Unprocessed data */

    void *private_data;	     /* The 'private data' block and size from SCF, */
    int private_size;        /*         NULL & 0 if not present.            */
} Read;


/*
 *-----------------------------------------------------------------------------
 * Function prototypes
 *-----------------------------------------------------------------------------
 */


/* ----- Main I/O routines ----- */

/*
 * Read a sequence from a file "fn" of format "format". If "format" is 0
 * (TT_ANY), we automatically determine the correct format.
 *
 * Returns:
 *   Read *   for success
 *   NULLRead for failure
 */
Read *read_reading(char *fn, int format);
Read *fread_reading(FILE *fp, char *fn, int format);
Read *mfread_reading(mFILE *fp, char *fn, int format);


/*
 * Write a sequence to a file "fn" of format "format". If "format" is 0,
 * we choose our favourite - SCF.
 *
 * Returns:
 *   0 for success
 *  -1 for failure
 */
int write_reading(char *fn, Read *read, int format);
int fwrite_reading(FILE *fp, Read *read, int format);
int mfwrite_reading(mFILE *fp, Read *read, int format);


/* ----- Utility routines ----- */

/*
 * Allocate a new sequence, with the given sizes.
 * Returns:
 *   "Read *" for success
 *   "NULLRead" for failure
 */
Read *read_allocate(int num_points, int num_bases);


/*
 * Duplicates the read structure and optionally gives it a new filename.
 * The following fields are not duplicated:
 *    
 *  int  orig_trace_format;
 *  void (*orig_trace_free)(void *ptr);
 *  void *orig_trace;
 *  char *ident;
 *
 * Returns:
 *   "Read *" for success
 *   "NULLRead" for failure
 */
Read* read_dup( Read* src, const char* new_name );


/*
 * Free memory allocated to a sequence by read_allocate().
 */
void read_deallocate(Read *read);

/* unix specific file deletion routine */

int remove_file(char *fn);

Read *read_abi(char *fn);
Read *fread_abi(FILE *fp);
Read *mfread_abi(mFILE *fp);
int write_abi(char *fn, Read *read);
int fwrite_abi(FILE *fp, Read *read);
int mfwrite_abi(mFILE *fp, Read *read);

int write_alf(char *fn, Read *read);
int fwrite_alf(FILE *fp, Read *read);
int mfwrite_alf(mFILE *fp, Read *read);
Read *read_alf(char *fn);
Read *fread_alf(FILE *fp);
Read *mfread_alf(mFILE *fp);

int write_pln(char *fn, Read *read);
int fwrite_pln(FILE *fp, Read *read);
int mfwrite_pln(mFILE *fp, Read *read);
Read *read_pln(char *fn);
Read *fread_pln(FILE *fp);
Read *mfread_pln(mFILE *fp);

Read *read_ctf(char *fn);
Read *fread_ctf(FILE *fp);
Read *mfread_ctf(mFILE *fp);
int write_ctf(char *fn, Read *read);
int fwrite_ctf(FILE *fp, Read *read);
int mfwrite_ctf(mFILE *fp, Read *read);

int read_sections(int sec);

#include "io_lib/translate.h"
#include "io_lib/compress.h"

#ifdef __cplusplus
}
#endif

#endif /* _Read_h_ */
