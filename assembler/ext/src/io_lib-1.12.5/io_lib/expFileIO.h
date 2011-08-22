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
 * File: expFileIO.h
 * Version:
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */

#ifndef _EXPFILEIO_H_
#define _EXPFILEIO_H_

#include <stdio.h>

#include "io_lib/mFILE.h"
#include "io_lib/array.h"
#include "io_lib/os.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Definitions
 */
#define MAXIMUM_EFLT_LENGTH     4
#define MAXIMUM_EFLTS          60
#define EFLT_FILE_LINE_LENGTH 128
#define EXP_FILE_LINE_LENGTH  128

typedef Array Exp_entries;

typedef struct {
    Array entries[MAXIMUM_EFLTS]; /* array of array of entries */
    int Nentries[MAXIMUM_EFLTS];
    mFILE *fp;
} Exp_info;

#define NULL_Exp_info ( (Exp_info *) NULL )



#define exp_Nentries(E,I) ((E)->Nentries[I]) /* get last entry for line I */

#define exp_get_entry(E,I) (arr(char *,(E)->entries[I],(E)->Nentries[I] - 1)) /* get last entry for line I */

/*
 * Allocate an set a new experiment file entry
 */
extern char *exp_set_entry(Exp_info *e, int eflt, char *str);



/*************************************************************
 * Experiment file line types
 *************************************************************/ 

extern char eflt_feature_ids[MAXIMUM_EFLTS][MAXIMUM_EFLT_LENGTH+1];

#define EFLT_CF  0
#define EFLT_CN  1
#define EFLT_CS  2
#define EFLT_CV  3
#define EFLT_DR  4
#define EFLT_DT  5
#define EFLT_EN  6
#define EFLT_EX  7
#define EFLT_FM  8
#define EFLT_LN  9
#define EFLT_LT 10
#define EFLT_MC 11
#define EFLT_MN 12
#define EFLT_MT 13
#define EFLT_OP 14
#define EFLT_PN 15
#define EFLT_QR 16
#define EFLT_SC 17
#define EFLT_SF 18
#define EFLT_SI 19
#define EFLT_SL 20
#define EFLT_SP 21
#define EFLT_SQ 22
#define EFLT_SR 23
#define EFLT_ST 24
#define EFLT_SV 25
#define EFLT_TN 26
#define EFLT_QL 27
#define EFLT_PS 28
#define EFLT_CC 29
#define EFLT_SS 30
#define EFLT_TG 31
#define EFLT_ID 32
#define EFLT_AQ 33
#define EFLT_PR 34
#define EFLT_LI 35
#define EFLT_LE 36
#define EFLT_TC 37
#define EFLT_AC 38
#define EFLT_BC 39
#define EFLT_ON 40
#define EFLT_AV 41
#define EFLT_PC 42
#define EFLT_SE 43
#define EFLT_CL 44
#define EFLT_CR 45
#define EFLT_AP 46
#define EFLT_CH 47
#define EFLT_PD 48
#define EFLT_WT 49
#define EFLT_NT 50
#define EFLT_GD 51
#define EFLT_WL 52
#define EFLT_WR 53
#define EFLT_FT 54
#define EFLT_LG 55


/*************************************************************************************/


/*
 * Creates a string of 'range format' from the start and end points.
 * The string (of form start..end) is also returned.
 */
extern char *exp_create_range(char *str, int start, int end);

/*
 * Extracts the start and end points from a range string.
 * Returns 0 for success and -1 for failure.
 */
extern int exp_extract_range(char *str, int *start, int *end);

/*
 * Output an experiment file line
 */
extern int exp_print_line(mFILE *fp, Exp_info *e, int eflt, int i);

/*
 * Output an experiment file multi-line
 */
extern int exp_print_mline(mFILE *fp, Exp_info *e, int eflt, int i);

extern int exp_print_seq(mFILE *fp, Exp_info *e, int eflt, int i);
/*
 * Output an experiment file multi line
 */



extern int exp_get_feature_index(char *e);

extern void exp_destroy_info(Exp_info *e);
/*
 * Destroy experiment file information
 */



extern Exp_info *exp_create_info(void);
/*
 * Allocate space for new experiment file information
 */







extern Exp_info *exp_fread_info(FILE *fp);
extern Exp_info *exp_mfread_info(mFILE *fp);
extern Exp_info *exp_read_info(char *file);
/*
 * Read in an experiment file and return handle
 */


char *opos2str(int2 *opos, int len, char *buf);
int   str2opos(int2 *opos, int len, char *buf);
char *conf2str(int1 *conf, int len, char *buf);
int   str2conf(int1 *conf, int len, char *buf);

extern int exp_get_int(Exp_info *e, int id, int *val);
/*
 * Get the integer for entry id
 * returns:
 *    0 - success
 *    1 - no entry
 */


extern int exp_get_rng(Exp_info *e, int id, int *from, int *to);
/*
 * Get the integer pair for entry id
 * returns:
 *    0 - success
 *    1 - no entry
 */


extern int exp_get_str(Exp_info *e, int id, char *s, f_implicit s_l);
/*
 * Get the string for entry id
 * returns:
 *    0 - success
 *    1 - no entry
 */


extern int exp_put_int(Exp_info *e, int id, int *val);
/*
 * Append the integer for entry id to the experiment file
 * returns:
 *    0 - success
 *    1 - no update
 */


extern int exp_put_rng(Exp_info *e, int id, int *from, int *to);
/*
 * Append the integer pair for entry id to the experiment file
 * returns:
 *    0 - success
 *    1 - no update
 */



extern int exp_put_str(Exp_info *e, int id, char *s, f_implicit s_l);
/*
 * Append the string for entry id to the experiment file
 * returns:
 *    0 - success
 *    1 - no update
 */


extern void exp_close(Exp_info *e);
/*
 * Closes an experiment file (if open), but does not free it.
 */

/*
 * FORTRAN INTERFACE
 */



extern f_int expopn_(char *fn, f_implicit fn_l);
/*
 * FORTRAN interface to exp_open_file()
 */

extern f_proc_ret expkil_(f_int *le);
/*
 * FORTRAN interface to exp_destroy_info
 */

extern f_int expri_(f_int *le, f_int *id, f_int *val);
/*
 * FORTRAN interface to exp_get_int
 */


extern f_int exprr_(f_int *le, f_int *id, f_int *from, f_int *to);
/*
 * FORTRAN interface to exp_get_rng
 */


extern f_int exprsa_(f_int *le, f_int *id, char *s, f_int *max_len, f_implicit s_l);
/*
 * FORTRAN interface to exp_get_str workalike
 * NOTE: for use with FORTRAN CHARACTER arrays instead CHARACTER strings
 */

extern f_int exprs_(f_int *le, f_int *id, char *s, f_implicit s_l);
/*
 * FORTRAN interface to exp_get_str workalike
 * NOTE: for use with FORTRAN CHARACTER strings instead CHARACTER arrays
 */

extern f_int expwi_(f_int *le, f_int *id, f_int *val);
/*
 * FORTRAN interface to exp_put_int
 */


extern f_int expwr_(f_int *le, f_int *id, f_int *from, f_int *to);
/*
 * FORTRAN interface to exp_put_rng
 */



extern f_int expwsa_(f_int *le, f_int *id, char *s, f_int *max_len, f_implicit s_l);
/*
 * FORTRAN interface to exp_put_str workalike
 * NOTE: for use with FORTRAN CHARACTER arrays instead CHARACTER strings
 */



extern f_int expws_(f_int *le, f_int *id, char *s, f_implicit s_l);
/*
 * FORTRAN interface to exp_put_str workalike
 * NOTE: for use with FORTRAN CHARACTER strings instead CHARACTER arrays
 */


extern void exp_print_file(FILE *fp, Exp_info *e);
extern void exp_print_mfile(mFILE *fp, Exp_info *e);

/*
 * FORTRAN interface to exp_create_range()
 */
extern void expcr_(char *str, f_int *start, f_int *end, f_implicit str_l);

/*
 * FORTRAN interface to exp_extract_range()
 */
extern f_int exper_(char *str, f_int *start, f_int *end, f_implicit str_l);

#ifdef __cplusplus
}
#endif

#endif /* _EXPFILEIO_H_ */

