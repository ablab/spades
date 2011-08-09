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

#ifndef _COMPRESS_H
#define _COMPRESS_H

#include <stdio.h>

#include "io_lib/mFILE.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Compress a file using the method set in the compression_used value
 * (set by set_compression_method and fopen_compressed).
 *
 * If compression succeeds, we rename the file back its original name.
 *
 * When compression_used is 0 no compression is done.
 */
int compress_file(char *file);
int fcompress_file(mFILE *fp);

/*
 * Returns a file pointer of an uncompressed copy of 'file'.
 * 'file' need not exist if 'file'.ext (eg file.gz)
 * exists and can be uncompressed.
 *
 * If ofp is non NULL then the original file pointer will also be returned
 * (opened for update) to allow writing back to the original file. In cases
 * of uncompressed data this is the same as the returned file pointer.
 */
mFILE *fopen_compressed(char *file, mFILE **ofp);

/*
 * Returns a file pointer of an uncompressed copy of 'fp'.
 *
 * If ofp is non NULL then the original file pointer will also be returned
 * (opened for update) to allow writing back to the original file. In cases
 * of uncompressed data this is the same as the returned file pointer.
 */
mFILE *freopen_compressed(mFILE *fp, mFILE **ofp);

/*
 * Sets the desired compression method. The below macros relate to entries
 * in the compression magic numbers table.
 */
void set_compression_method(int method);
int get_compression_method(void);

/*
 * Converts compress mode strings (eg "gzip") to numbers.
 */
int compress_str2int(char *mode);

/*
 * Converts compress mode numbers to strings (eg "gzip").
 */
char *compress_int2str(int mode);

#define COMP_METHOD_NONE	0
#define COMP_METHOD_BZIP	1
#define COMP_METHOD_GZIP	2
#define COMP_METHOD_COMPRESS	3
#define COMP_METHOD_BZIP2	5
#define COMP_METHOD_SZIP	6

/*
 * In memory gzip and gunzip using zlib. In each case data is the input
 * data (with 'size') and the returned value is the output data along with
 * [cu]data_size being the returned data size.
 *
 * Both return NULL on failure
 */
char *memgzip(char *data, size_t size, size_t *cdata_size);
char *memgunzip(char *data, size_t size, size_t *udata_size);

/*
 * Given a filename remove a known compression extension
 *
 * Returns: None
 */
void remove_extension(char *file);

#ifdef __cplusplus
}
#endif

#endif /* _COMPRESS_H */
