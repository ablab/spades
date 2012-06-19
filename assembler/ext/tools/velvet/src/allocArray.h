/*
Copyright 2009 Sylvain Foret (sylvain.foret@anu.edu.au) 

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#ifndef _ALLOC_ARRAY_H_
#define _ALLOC_ARRAY_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include "globals.h"

typedef struct AllocArray_st AllocArray;
typedef struct AllocArrayFreeElement_st AllocArrayFreeElement;

struct AllocArray_st
{
	void **blocks;
	AllocArrayFreeElement *freeElements;
	size_t elementSize;
	size_t blockSize;
	size_t maxBlocks;
	size_t currentBlocks;
	size_t maxElements;
	size_t currentElements;
#ifdef DEBUG
	char *name;
	size_t elementsRecycled;
	size_t elementsAllocated;
#endif
#ifdef _OPENMP
	int nbThreads;
#endif
};

AllocArray* newAllocArray (size_t elementSize, char *name);
void destroyAllocArray (AllocArray *array);
ArrayIdx allocArrayAllocate (AllocArray *array);
void allocArrayFree (AllocArray *array, ArrayIdx idx);

#define DECLARE_FAST_ACCESSORS(name, type, array) \
/* Fast version, without null pointer checks */ \
static inline type* name##_FI2P(ArrayIdx idx) \
{ \
	const ArrayIdx i = idx - 1; \
	const ArrayIdx blockIdx = i / array->maxElements; \
	const ArrayIdx elementIdx = i % array->maxElements; \
	return &((type*)(array->blocks[blockIdx]))[elementIdx]; \
} \
/* Slower version, with null pointer checks */ \
static inline type* name##_I2P(ArrayIdx idx) \
{ \
	if (idx != NULL_IDX) \
		return name##_FI2P(idx); \
	return NULL; \
} 

#ifdef _OPENMP
// For multithreading: thread-specific alloc arrays 
AllocArray *newAllocArrayArray(unsigned int n,
			       size_t elementSize,
			       char * name);
void destroyAllocArrayArray(AllocArray * allocArray);
ArrayIdx allocArrayArrayAllocate (AllocArray *array);
void allocArrayArrayFree (AllocArray *array, ArrayIdx idx);
#endif

#endif /* _ALLOC_ARRAY_H_ */
