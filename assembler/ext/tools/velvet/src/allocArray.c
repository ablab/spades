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

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "allocArray.h"
#include "utility.h"

#ifdef DEBUG
#define NB_PAGES_ALLOC    1
#define BLOCKS_ALLOC_SIZE 1
#else
#define NB_PAGES_ALLOC    128
#define BLOCKS_ALLOC_SIZE 128
#endif

#ifndef VBIGASSEMBLY
#define INDEX_LENGTH 32
#else
#define INDEX_LENGTH 64
#endif

struct AllocArrayFreeElement_st
{
	AllocArrayFreeElement *next;
	ArrayIdx idx;
};

static void initAllocArray(AllocArray * array, size_t elementSize, char * name) 
{
	array->freeElements = NULL;
	array->elementSize = elementSize;
	array->blockSize = sysconf (_SC_PAGESIZE) * NB_PAGES_ALLOC;
	array->maxElements = array->blockSize / array->elementSize;
	array->maxBlocks = BLOCKS_ALLOC_SIZE;
	array->blocks = mallocOrExit (array->maxBlocks, void*);
	array->blocks[0] = mallocOrExit (array->blockSize, char);
	array->currentBlocks = 1;
	array->currentElements = 0;
#ifdef DEBUG
	array->elementsRecycled = 0;
	array->elementsAllocated = 0;
	array->name = name;
#endif
#ifdef _OPENMP
	array->nbThreads = 0;
#endif
}

AllocArray*
newAllocArray (size_t elementSize, char *name)
{
	AllocArray *array;

	if (elementSize < sizeof(AllocArrayFreeElement)) {
		velvetLog("Elements too small to create an AllocArray!\n");
		exit(-1);
	}
	array = mallocOrExit (1, AllocArray);
	initAllocArray(array, elementSize, name);
	return array;
}

void
destroyAllocArrayChunks (AllocArray *array)
{
	size_t i;
	for (i = 0; i < array->currentBlocks; i++)
		free (array->blocks[i]);
	free (array->blocks);
}

void
destroyAllocArray (AllocArray *array)
{
	if (array)
	{
		destroyAllocArrayChunks(array);
#ifdef DEBUG
		velvetLog(">>> Allocation summary for %s\n", array->name);
		velvetLog(">>> Alloc'ed %ld bytes\n", array->blockSize * array->currentBlocks);
		velvetLog(">>> Alloc'ed %ld elements\n", array->elementsAllocated);
		velvetLog(">>> Recycled %ld elements\n", array->elementsRecycled);
#endif
		free (array);
	}
}

ArrayIdx
allocArrayAllocate (AllocArray *array)
{
	if (array->freeElements != NULL)
	{
		AllocArrayFreeElement *element;

		element = array->freeElements;
		array->freeElements = element->next;
#ifdef DEBUG
		array->elementsRecycled++;
#endif
		return element->idx;
	}
	if (array->currentElements >= array->maxElements)
	{
		if (array->currentBlocks == array->maxBlocks)
		{
			array->maxBlocks += BLOCKS_ALLOC_SIZE;
			array->blocks = reallocOrExit (array->blocks, array->maxBlocks, void*);
		}
		array->blocks[array->currentBlocks] = mallocOrExit (array->blockSize, char);
		array->currentBlocks++;
		array->currentElements = 0;
	}
	array->currentElements++;
#ifndef VBIGASSEMBLY
	if (array->maxElements * (array->currentBlocks - 1) + array->currentElements == UINT32_MAX)
	{
#ifdef DEBUG
		velvetLog (">>> Reached maximum `%s' addressable with %i bits\n", array->name, INDEX_LENGTH);
#else
		velvetLog (">>> Reached maximum elements addressable with %i bits\n", INDEX_LENGTH);
#endif
		abort();
	}
#endif

#ifdef DEBUG
	array->elementsAllocated++;
#endif

	return array->maxElements * (array->currentBlocks - 1) + array->currentElements;
}

static void*
allocArrayGetElement (AllocArray *array, ArrayIdx idx)
{
	if (idx != NULL_IDX)
	{
		const ArrayIdx i = idx - 1;
		const ArrayIdx blockIdx = i / array->maxElements;
		const ArrayIdx elementIdx = i % array->maxElements;
		return ((char*)(array->blocks[blockIdx])) + elementIdx * array->elementSize;
	}
	return NULL;
}

void
allocArrayFree (AllocArray *array, ArrayIdx idx)
{
	if (idx != NULL_IDX)
	{
		AllocArrayFreeElement *freeElem;

		freeElem = allocArrayGetElement (array, idx);
		freeElem->idx = idx;
		freeElem->next = array->freeElements;
		array->freeElements = freeElem;
	}
}

#ifdef _OPENMP

#define BLOCKS_ALLOC_SHIFT 16

static void initAllocArrayArray(AllocArray *array,
				void *blocks,
				size_t nbBlocks,
				size_t elementSize,
				char *name,
				int nbThreads,
				int thisThread)
{
	array->freeElements = NULL;
	array->elementSize = elementSize;
	array->blockSize = (((size_t) 1) << BLOCKS_ALLOC_SHIFT) * elementSize;
	array->maxElements = (((size_t) 1) << BLOCKS_ALLOC_SHIFT);
	array->maxBlocks = nbBlocks;
	array->blocks = blocks;
	array->blocks[thisThread] = callocOrExit(array->blockSize, char);
	array->currentBlocks = thisThread;
	array->currentElements = 0;
	array->nbThreads = nbThreads;
#ifdef DEBUG
	array->elementsRecycled = 0;
	array->elementsAllocated = 0;
	array->name = name;
#endif
}

AllocArray *newAllocArrayArray(unsigned int n,
			       size_t elementSize,
			       char * name)
{
	AllocArray *allocArray;
	void **blocks;
	size_t nbBlocks;
	int i;

	allocArray = callocOrExit (n + 1, AllocArray);
	nbBlocks = (((size_t) 1) << (INDEX_LENGTH - BLOCKS_ALLOC_SHIFT));
	blocks = callocOrExit(nbBlocks, void*);
	for (i = 0; i < n; i++)
		initAllocArrayArray(allocArray + i,
				    blocks,
				    nbBlocks,
				    elementSize,
				    name,
				    n, i);
	/* Last element marker */
	allocArray[n].currentBlocks = n;

	return allocArray;
}

void destroyAllocArrayArray(AllocArray *allocArray)
{
	int i;

	for (i = 0; i < allocArray[0].maxBlocks; i++)
		if (allocArray[0].blocks[i] != NULL)
			free(allocArray[0].blocks[i]);

	free(allocArray[0].blocks);
	free(allocArray);
}

ArrayIdx allocArrayArrayAllocate(AllocArray *array)
{
	int thread = omp_get_thread_num();

	AllocArray * lastArray = array + array->nbThreads;
	array += thread;
	if (array->freeElements != NULL)
	{
		AllocArrayFreeElement *element;

		element = array->freeElements;
		array->freeElements = element->next;
#ifdef DEBUG
		array->elementsRecycled++;
#endif
		return element->idx;
	}
	if (array->currentElements >= array->maxElements)
	{
		#pragma omp critical 
		array->currentBlocks = lastArray->currentBlocks++;

		if (array->currentBlocks >= array->maxBlocks)
		{
#ifdef DEBUG
			velvetLog(">>> Reached maximum `%s' addressable with %i bits\n",
				  array->name, INDEX_LENGTH);
#else
			velvetLog(">>> Reached maximum elements addressable with %i bits\n", INDEX_LENGTH);
#endif
			abort();
		}
		array->blocks[array->currentBlocks] = callocOrExit(array->blockSize, char);
		array->currentElements = 0;
	}
	array->currentElements++;
#ifdef DEBUG
	array->elementsAllocated++;
#endif
	return array->maxElements * array->currentBlocks + array->currentElements;
}

void
allocArrayArrayFree(AllocArray *array, ArrayIdx idx)
{
	if (idx != NULL_IDX)
	{
		AllocArrayFreeElement *freeElem;

		freeElem = allocArrayGetElement(array, idx);
		freeElem->idx = idx;
		freeElem->next = array->freeElements;
		array->freeElements = freeElem;
	}
}
#endif
