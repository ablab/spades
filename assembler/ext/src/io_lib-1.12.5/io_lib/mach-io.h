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

#ifndef _mach_io_h
#define _mach_io_h
/*
 * Machine independant io
 * For reading and writing to big-endian and little-endian files
 *
 * Routines available:
 *     be_write_int_1()
 *     be_write_int_2()
 *     be_write_int_4()
 *     be_read_int_1()
 *     be_read_int_2()
 *     be_read_int_4()
 *     le_write_int_1()
 *     le_write_int_2()
 *     le_write_int_4()
 *     le_read_int_1()
 *     le_read_int_2()
 *     le_read_int_4()
 *
 * All routine return:
 *    0 - an error has occurred during io operation
 *    1 - value successfully read or written
 */

#include <stdio.h>
#include "io_lib/os.h"
#include "io_lib/mFILE.h"

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/
/* IO for big-endian files                                            */
/**********************************************************************/

/*
 * Write a big-endian int1
 */
extern int be_write_int_1(mFILE *fp, uint1 *i1);

/*
 * Write a big-endian int2
 */
extern int be_write_int_2(mFILE *fp, uint2 *i2);

/*
 * Write a big-endian int4
 */
extern int be_write_int_4(mFILE *fp, uint4 *i4);

/*
 * Write a big-endian int8
 */
extern int be_write_int_8(mFILE *fp, uint8 *i8);

/*
 * Read a big-endian int1
 */
extern int be_read_int_1(mFILE *fp, uint1 *i1);

/*
 * Read a big-endian int2
 */
extern int be_read_int_2(mFILE *fp, uint2 *i2);

/*
 * Read a big-endian int4
 */
extern int be_read_int_4(mFILE *fp, uint4 *i4);

/*
 * Read a big-endian int8
 */
extern int be_read_int_8(mFILE *fp, uint8 *i8);

/**********************************************************************/
/* IO for little-endian files                                         */
/**********************************************************************/

/*
 * Write a little-endian int1
 */
extern int le_write_int_1(mFILE *fp, uint1 *i1);

/*
 * Write a little-endian int2
 */
extern int le_write_int_2(mFILE *fp, uint2 *i2);

/*
 * Write a little-endian int4
 */
extern int le_write_int_4(mFILE *fp, uint4 *i4);

/*
 * Write a little-endian int8
 */
extern int le_write_int_8(mFILE *fp, uint8 *i8);

/*
 * Read a little-endian int1
 */
extern int le_read_int_1(mFILE *fp, uint1 *i1);

/*
 * Read a little-endian int2
 */
extern int le_read_int_2(mFILE *fp, uint2 *i2);

/*
 * Read a little-endian int4
 */
extern int le_read_int_4(mFILE *fp, uint4 *i4);

/*
 * Read a little-endian int8
 */
extern int le_read_int_8(mFILE *fp, uint8 *i8);

#ifdef __cplusplus
}
#endif

#endif /* _mach_io_h */
