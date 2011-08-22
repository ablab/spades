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

#ifndef _traceType_h
#define _traceType_h

/*
 * Title:       traceType.h
 * 
 * File:        traceType.h
 * Purpose:     determining traceType of traces
 * Last update: Tue Jan 15 1991
 *
 * Change log :-
 */

/* ---- Imports ---- */

#include <stdio.h>      /* IMPORT: fopen, fclose, fseek, ftell, fgetc */
#include <ctype.h>      
#include <string.h>     /* IMPORT: isprint*/

#include "io_lib/Read.h"	/* IMPORT: TT_xxx defines */
#include "io_lib/mFILE.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Determine the trace type for file 'fn'.
 *
 * Returns:
 *     TT_SCF, TT_ABI, TT_ALF, or TT_PLN for success.
 *     TT_UNK for unknown type.
 *     TT_ERR for error.
 */
extern int determine_trace_type(char *fn);

/*
 * Determine the trace type for FILE * 'fp'.
 *
 * Returns:
 *     TT_SCF, TT_ABI, TT_ALF, or TT_PLN for success.
 *     TT_UNK for unknown type.
 *     TT_ERR for error.
 */
extern int fdetermine_trace_type(mFILE *fp);

/*
 * Returns a statically declared string containing a 3 character
 * identifier for this trace type.
 * "ERR" represents error, and "UNK" for unknown.
 * Successful values are "SCF", "ABI", "ALF" and "PLN".
 */
extern char *trace_type_str(char *traceName);

/*
 * Converts a trace type string to an integer.
 */
extern int trace_type_str2int(char *str);

/*
 * Converts a trace type integer to a string.
 */
char *trace_type_int2str(int type);


#ifdef __cplusplus
}
#endif

#endif /*_traceType_h*/
