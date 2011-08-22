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

#ifndef _MY_MALLOC_H_
#define _MY_MALLOC_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void *xmalloc(size_t size);
extern void *xrealloc(void *ptr, size_t size);
extern void *xcalloc(size_t num, size_t size);
extern void xfree(void *ptr);

#ifdef __cplusplus
}
#endif

#endif
