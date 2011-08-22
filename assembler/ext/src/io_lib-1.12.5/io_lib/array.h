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
 * File: array.h
 * Version:
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */

#ifndef _ARRAY_H_
#define _ARRAY_H_

/* 12/1/99 johnt - use stddef.h not sys/types.h for size_t */
#include <stddef.h>		/* IMPORT: size_t */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t size;		/* element size */
    size_t dim;			/* allocated number of elements */
    size_t max;			/* elements accessed */
    void *base;			/* base address of array */
} ArrayStruct, *Array;



extern Array ArrayCreate(size_t size, size_t dim);

extern int ArrayExtend(Array a, size_t dim);

extern void *ArrayRef(Array a, size_t i);

extern int ArrayDestroy(Array a);

#define ArrayMax(a) ( (a)->max )

#define ArrayBase(t,a) ( (t *)((a)->base) )

/*
#define arr(t,a,n) \
    (*(t*)((a)->base + (a)->size*(n)))

#define arrp(t,a,n) \
    ((t*)((a)->base + (a)->size*(n)))
*/




#define arr(t,a,n) \
    ((t*)((a)->base))[n]

#define ARR(t,a,n) \
    (*((t*)ArrayRef((a),(n))))

#define arrp(t,a,n) \
    &((t*)((a)->base))[n]
#define ARRP(t,a,n) \
    ((t*)ArrayRef(a,n))

#define ARRAY_NO_ERROR			 0
#define ARRAY_FULL			-1
#define ARRAY_INVALID_ARGUMENTS		-2
#define ARRAY_OUT_OF_MEMORY		-3

extern int ArrayError;

extern char *ArrayErrorString(int error);


#ifdef __cplusplus
}
#endif

#endif /*_ARRAY_H_*/
