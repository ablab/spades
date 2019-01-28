/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_H_
#define _PARASAIL_H_

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Version macros for compile-time API version detection */
#define PARASAIL_VERSION_MAJOR 2
#define PARASAIL_VERSION_MINOR 3
#define PARASAIL_VERSION_PATCH 0

#define PARASAIL_MAKE_VERSION(major, minor, patch) \
    ((major) * 10000 + (minor) * 100 + (patch))
#define PARASAIL_VERSION \
    PARASAIL_MAKE_VERSION(PARASAIL_VERSION_MAJOR, \
                          PARASAIL_VERSION_MINOR, \
                          PARASAIL_VERSION_PATCH)

/* for traceback */
#define PARASAIL_ZERO_MASK 120 /* all bits set except the first three */
#define PARASAIL_E_MASK 103 /* all bits set except the E bits */
#define PARASAIL_F_MASK 31 /* all bits set except the F bits */
#define PARASAIL_ZERO   0
#define PARASAIL_INS    1
#define PARASAIL_DEL    2
#define PARASAIL_DIAG   4
#define PARASAIL_DIAG_E 8
#define PARASAIL_INS_E  16
#define PARASAIL_DIAG_F 32
#define PARASAIL_DEL_F  64

/*                                            3         2         1          */
/*                                           10987654321098765432109876543210*/
#define PARASAIL_FLAG_NW          (1 << 0) /*00000000000000000000000000000001*/
#define PARASAIL_FLAG_SG          (1 << 1) /*00000000000000000000000000000010*/
#define PARASAIL_FLAG_SW          (1 << 2) /*00000000000000000000000000000100*/
#define PARASAIL_FLAG_SATURATED   (1 << 6) /*00000000000000000000000001000000*/
#define PARASAIL_FLAG_BANDED      (1 << 7) /*00000000000000000000000010000000*/
#define PARASAIL_FLAG_NOVEC       (1 << 8) /*00000000000000000000000100000000*/
#define PARASAIL_FLAG_NOVEC_SCAN  (1 << 9) /*00000000000000000000001000000000*/
#define PARASAIL_FLAG_SCAN        (1 <<10) /*00000000000000000000010000000000*/
#define PARASAIL_FLAG_STRIPED     (1 <<11) /*00000000000000000000100000000000*/
#define PARASAIL_FLAG_DIAG        (1 <<12) /*00000000000000000001000000000000*/
#define PARASAIL_FLAG_BLOCKED     (1 <<13) /*00000000000000000010000000000000*/
#define PARASAIL_FLAG_STATS       (1 <<16) /*00000000000000010000000000000000*/
#define PARASAIL_FLAG_TABLE       (1 <<17) /*00000000000000100000000000000000*/
#define PARASAIL_FLAG_ROWCOL      (1 <<18) /*00000000000001000000000000000000*/
#define PARASAIL_FLAG_TRACE       (1 <<19) /*00000000000010000000000000000000*/
#define PARASAIL_FLAG_BITS_8      (1 <<20) /*00000000000100000000000000000000*/
#define PARASAIL_FLAG_BITS_16     (1 <<21) /*00000000001000000000000000000000*/
#define PARASAIL_FLAG_BITS_32     (1 <<22) /*00000000010000000000000000000000*/
#define PARASAIL_FLAG_BITS_64     (1 <<23) /*00000000100000000000000000000000*/
#define PARASAIL_FLAG_LANES_1     (1 <<24) /*00000001000000000000000000000000*/
#define PARASAIL_FLAG_LANES_2     (1 <<25) /*00000010000000000000000000000000*/
#define PARASAIL_FLAG_LANES_4     (1 <<26) /*00000100000000000000000000000000*/
#define PARASAIL_FLAG_LANES_8     (1 <<27) /*00001000000000000000000000000000*/
#define PARASAIL_FLAG_LANES_16    (1 <<28) /*00010000000000000000000000000000*/
#define PARASAIL_FLAG_LANES_32    (1 <<29) /*00100000000000000000000000000000*/
#define PARASAIL_FLAG_LANES_64    (1 <<30) /*01000000000000000000000000000000*/
#define PARASAIL_FLAG_INVALID  0x8000C038  /*10000000000000001100000000111000*/

/*
 * This helps users not familiar with the restrict keyword.
 */
#if !defined(restrict) && (defined(__cplusplus) || __STDC_VERSION__ < 199901L)
#define restrict
#define PARASAIL_RESTRICT_REMOVED
#endif

/* all the possible combinations of results so we can reduce the memory
 * footprint of the result */

typedef struct parasail_result_extra_stats_tables {
    int * restrict score_table;     /* DP table of scores */
    int * restrict matches_table;   /* DP table of exact match counts */
    int * restrict similar_table;   /* DP table of similar substitution counts */
    int * restrict length_table;    /* DP table of lengths */
} parasail_result_extra_stats_tables_t;

typedef struct parasail_result_extra_stats_rowcols {
    int * restrict score_row;       /* last row of DP table of scores */
    int * restrict matches_row;     /* last row of DP table of exact match counts */
    int * restrict similar_row;     /* last row of DP table of similar substitution counts */
    int * restrict length_row;      /* last row of DP table of lengths */
    int * restrict score_col;       /* last col of DP table of scores */
    int * restrict matches_col;     /* last col of DP table of exact match counts */
    int * restrict similar_col;     /* last col of DP table of similar substitution counts */
    int * restrict length_col;      /* last col of DP table of lengths */
} parasail_result_extra_stats_rowcols_t;

typedef struct parasail_result_extra_stats {
    int matches;    /* number of exactly matching characters in the alignment */
    int similar;    /* number of similar characters (positive substitutions) in the alignment */
    int length;     /* length of the alignment */
    union {
        void *extra;
        parasail_result_extra_stats_tables_t *tables;
        parasail_result_extra_stats_rowcols_t *rowcols;
    };
} parasail_result_extra_stats_t;

typedef struct parasail_result_extra_tables {
    int * restrict score_table;     /* DP table of scores */
} parasail_result_extra_tables_t;

typedef struct parasail_result_extra_rowcols {
    int * restrict score_row;       /* last row of DP table of scores */
    int * restrict score_col;       /* last col of DP table of scores */
} parasail_result_extra_rowcols_t;

typedef struct parasail_result_extra_trace {
    void * restrict trace_table;    /* DP table of traceback */
    void * restrict trace_ins_table;/* DP table of insertions traceback */
    void * restrict trace_del_table;/* DP table of deletions traceback */
} parasail_result_extra_trace_t;

typedef struct parasail_result {
    int score;      /* alignment score */
    int end_query;  /* end position of query sequence */
    int end_ref;    /* end position of reference sequence */
    int flag;       /* bit field for various flags */
    /* union of pointers to extra result data based on the flag */
    union {
        void *extra;
        parasail_result_extra_stats_t *stats;
        parasail_result_extra_tables_t *tables;
        parasail_result_extra_rowcols_t *rowcols;
        parasail_result_extra_trace_t *trace;
    };
} parasail_result_t;

typedef struct parasail_matrix {
    const char * name;
    const int *matrix;
    const int *mapper;
    int size;
    int max;
    int min;
    int *user_matrix;
} parasail_matrix_t;

typedef struct parasail_profile_data {
    void * score;
    void * matches;
    void * similar;
} parasail_profile_data_t;

typedef struct parasail_profile {
    const char *s1;
    int s1Len;
    const parasail_matrix_t *matrix;
    struct parasail_profile_data profile8;
    struct parasail_profile_data profile16;
    struct parasail_profile_data profile32;
    struct parasail_profile_data profile64;
    void (*free)(void * profile);
    int stop;
} parasail_profile_t;

extern void parasail_profile_free(parasail_profile_t *profile);

typedef parasail_result_t* parasail_function_t(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix);

typedef struct parasail_function_info {
    parasail_function_t * pointer;
    const char * name;
    const char * alg;
    const char * type;
    const char * isa;
    const char * bits;
    const char * width;
    int lanes;
    char is_table;
    char is_rowcol;
    char is_trace;
    char is_stats;
    char is_ref;
} parasail_function_info_t;

typedef parasail_result_t* parasail_pfunction_t(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

typedef parasail_profile_t* parasail_pcreator_t(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix);

typedef struct parasail_pfunction_info {
    parasail_pfunction_t * pointer;
    parasail_pcreator_t * creator;
    const char * name;
    const char * alg;
    const char * type;
    const char * isa;
    const char * bits;
    const char * width;
    int lanes;
    char is_table;
    char is_rowcol;
    char is_trace;
    char is_stats;
    char is_ref;
} parasail_pfunction_info_t;

/* Run-time API version detection */
extern void parasail_version(int *major, int *minor, int *patch);

/** Deallocate result. */
extern void parasail_result_free(parasail_result_t *result);

/** Lookup function by name. */
extern parasail_function_t * parasail_lookup_function(const char *funcname);

/** Lookup pfunction by name. */
extern parasail_pfunction_t * parasail_lookup_pfunction(const char *funcname);

/** Lookup pcreator by name. */
extern parasail_pcreator_t * parasail_lookup_pcreator(const char *funcname);

/** Lookup function info by name. */
extern const parasail_function_info_t * parasail_lookup_function_info(const char *funcname);

/** Lookup function info by name. */
extern const parasail_pfunction_info_t * parasail_lookup_pfunction_info(const char *funcname);

/** Current time in seconds with nanosecond resolution. */
extern double parasail_time(void);

/** Lookup substitution matrix by name. */
extern const parasail_matrix_t* parasail_matrix_lookup(const char *matrixname);

/** Create simple substitution matrix. */
extern parasail_matrix_t* parasail_matrix_create(
        const char *alphabet, const int match, const int mismatch);

/** Create substitution matrix from a filename. */
extern parasail_matrix_t* parasail_matrix_from_file(const char *filename);

/** Deallocate substitution matrix. */
extern void parasail_matrix_free(parasail_matrix_t *matrix);

/** Copy any matrix into writeable matrix. */
extern parasail_matrix_t* parasail_matrix_copy(const parasail_matrix_t *original);

/** Modify a user matrix. */
extern void parasail_matrix_set_value(parasail_matrix_t *matrix, int row, int col, int value);

extern parasail_result_t* parasail_nw_banded(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int k,
        const parasail_matrix_t* matrix);

typedef struct parasail_traceback_{
    char *query;
    char *comp;
    char *ref;
} parasail_traceback_t;

extern parasail_traceback_t* parasail_result_get_traceback(
        parasail_result_t *result,
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix,
        char match, char pos, char neg);

extern void parasail_traceback_free(parasail_traceback_t *traceback);

extern void parasail_traceback_generic(
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const char *nameA,
        const char *nameB,
        const parasail_matrix_t *matrix,
        parasail_result_t *result,
        char match, char pos, char neg,
        int width,
        int name_width,
        int use_stats);

extern void parasail_traceback_generic_extra(
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const char *nameA,
        const char *nameB,
        const parasail_matrix_t *matrix,
        parasail_result_t *result,
        char match, char pos, char neg,
        int width,
        int name_width,
        int use_stats,
        int int_width,
        FILE *stream);

extern const uint8_t parasail_cigar_encoded_ops[];

typedef struct parasail_cigar_ {
    uint32_t *seq;
    int len;
    int beg_query;
    int beg_ref;
} parasail_cigar_t;

/**
 * Produce CIGAR 32-bit unsigned integer from CIGAR operation and CIGAR length.
 *
 * @param[in] length    length of CIGAR
 * @param[in] op_letter CIGAR operation character ('M', 'I', etc)
 * @return              encoded CIGAR operation and length
 */
extern uint32_t parasail_cigar_encode(uint32_t length, char op_letter);

/**
 * Produce CIGAR 32-bit unsigned integer array from CIGAR string.
 *
 * @param[in] cigar CIGAR string, e.g., '3=2I2=1X4D14='
 * @return          encoded CIGAR
 */
extern parasail_cigar_t* parasail_cigar_encode_string(const char *cigar);

/**
 * Extract CIGAR operation character from CIGAR 32-bit unsigned integer.
 *
 * @param[in] cigar_int   32-bit unsigned integer, representing encoded
 *                        CIGAR operation and length
 * @return                CIGAR operation character ('M', 'I', etc)
 */
extern char parasail_cigar_decode_op(uint32_t cigar_int);

/**
 * Extract length of a CIGAR operation from CIGAR 32-bit unsigned integer.
 *
 * @param[in] cigar_int   32-bit unsigned integer, representing encoded
 *                        CIGAR operation and length
 * @return                length of CIGAR operation
 */
extern uint32_t parasail_cigar_decode_len(uint32_t cigar_int);

/**
 * Convert CIGAR array into a character array.
 *
 * @param[in] cigar   32-bit unsigned integers, representing encoded
 *                    CIGAR operations and lengths
 * @return            CIGAR string, e.g., '3=2I2=1X4D14='
 */
extern char* parasail_cigar_decode(parasail_cigar_t *cigar);

/* allocate and return the cigar for the given alignment */
extern parasail_cigar_t* parasail_result_get_cigar(
        parasail_result_t *result,
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix);

/* free the cigar structure */
extern void parasail_cigar_free(parasail_cigar_t *cigar);

typedef struct parasail_result_ssw {
    uint16_t score1;
    int32_t ref_begin1;
    int32_t ref_end1;
    int32_t read_begin1;
    int32_t read_end1;
    uint32_t *cigar;
    int32_t cigarLen;
} parasail_result_ssw_t;

extern parasail_result_ssw_t* parasail_ssw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_ssw_t* parasail_ssw_profile(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_profile_t* parasail_ssw_init(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix, const int8_t score_size);

extern void parasail_result_ssw_free(parasail_result_ssw_t *result);

/* The following functions help access result attributes. */

extern int parasail_result_is_nw(const parasail_result_t * const restrict result);
extern int parasail_result_is_sg(const parasail_result_t * const restrict result);
extern int parasail_result_is_sw(const parasail_result_t * const restrict result);
extern int parasail_result_is_saturated(const parasail_result_t * const restrict result);
extern int parasail_result_is_banded(const parasail_result_t * const restrict result);
extern int parasail_result_is_scan(const parasail_result_t * const restrict result);
extern int parasail_result_is_striped(const parasail_result_t * const restrict result);
extern int parasail_result_is_diag(const parasail_result_t * const restrict result);
extern int parasail_result_is_blocked(const parasail_result_t * const restrict result);
extern int parasail_result_is_stats(const parasail_result_t * const restrict result);
extern int parasail_result_is_stats_table(const parasail_result_t * const restrict result);
extern int parasail_result_is_stats_rowcol(const parasail_result_t * const restrict result);
extern int parasail_result_is_table(const parasail_result_t * const restrict result);
extern int parasail_result_is_rowcol(const parasail_result_t * const restrict result);
extern int parasail_result_is_trace(const parasail_result_t * const restrict result);

extern int parasail_result_get_score(const parasail_result_t * const restrict result);
extern int parasail_result_get_end_query(const parasail_result_t * const restrict result);
extern int parasail_result_get_end_ref(const parasail_result_t * const restrict result);

extern int parasail_result_get_matches(const parasail_result_t * const restrict result);
extern int parasail_result_get_similar(const parasail_result_t * const restrict result);
extern int parasail_result_get_length(const parasail_result_t * const restrict result);

extern int* parasail_result_get_score_table(const parasail_result_t * const restrict result);
extern int* parasail_result_get_matches_table(const parasail_result_t * const restrict result);
extern int* parasail_result_get_similar_table(const parasail_result_t * const restrict result);
extern int* parasail_result_get_length_table(const parasail_result_t * const restrict result);
extern int* parasail_result_get_score_row(const parasail_result_t * const restrict result);
extern int* parasail_result_get_matches_row(const parasail_result_t * const restrict result);
extern int* parasail_result_get_similar_row(const parasail_result_t * const restrict result);
extern int* parasail_result_get_length_row(const parasail_result_t * const restrict result);
extern int* parasail_result_get_score_col(const parasail_result_t * const restrict result);
extern int* parasail_result_get_matches_col(const parasail_result_t * const restrict result);
extern int* parasail_result_get_similar_col(const parasail_result_t * const restrict result);
extern int* parasail_result_get_length_col(const parasail_result_t * const restrict result);
extern int* parasail_result_get_trace_table(const parasail_result_t * const restrict result);
extern int* parasail_result_get_trace_ins_table(const parasail_result_t * const restrict result);
extern int* parasail_result_get_trace_del_table(const parasail_result_t * const restrict result);

/* The following function signatures were generated by the 'names.py'
 * script located in the 'util' directory of the main distribution. */

/* BEGIN GENERATED NAMES */

extern parasail_result_t* parasail_nw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse2_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse2_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse2_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse2_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse2_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse41_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse41_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sse41_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_avx2_256_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_avx2_256_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_avx2_256_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_avx2_256_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_avx2_256_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse2_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse2_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse2_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse2_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse2_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse41_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse41_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse41_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse41_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sse41_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_avx2_256_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_avx2_256_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_avx2_256_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_avx2_256_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_avx2_256_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_altivec_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_altivec_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_altivec_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_altivec_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_altivec_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_neon_128_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_neon_128_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_neon_128_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_neon_128_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_neon_128_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_blocked_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_blocked_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_blocked_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_blocked_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_blocked_sse41_128_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_blocked_sse41_128_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_64(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_trace_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_trace_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_trace_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_table_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_64(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_32(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_16(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_8(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_profile_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap);

extern parasail_profile_t* parasail_profile_create_sse_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_sse_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_sse_128_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_sse_128_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_sse_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_avx_256_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_avx_256_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_avx_256_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_avx_256_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_avx_256_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_sse_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_sse_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_sse_128_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_sse_128_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_sse_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_avx_256_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_avx_256_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_avx_256_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_avx_256_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_avx_256_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_altivec_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_altivec_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_altivec_128_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_altivec_128_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_altivec_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_neon_128_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_neon_128_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_16(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_8(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_profile_t* parasail_profile_create_stats_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_trace_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_nw_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_trace_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sg_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_trace_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_table_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_scan_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_striped_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

extern parasail_result_t* parasail_sw_stats_rowcol_diag_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);

/* END GENERATED NAMES */

#ifdef __cplusplus
}
#endif

#ifdef PARASAIL_RESTRICT_REMOVED
#undef restrict
#endif

#endif /* _PARASAIL_H_ */
