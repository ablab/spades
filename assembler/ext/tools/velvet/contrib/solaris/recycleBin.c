/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

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
/****************************************************************\
*                                                                *
*  Efficient Memory Allocation Routines                          *
*                                                                *
*  Guy St.C. Slater..   mailto:guy@ebi.ac.uk                     *
*  Copyright (C) 2000-2005.  All Rights Reserved.                *
*                                                                *
*  This source code is distributed under the terms of the        *
*  GNU Lesser General Public License. See the file COPYING       *
*  or http://www.fsf.org/copyleft/lesser.html for details        *
*                                                                *
*  If you use this code, please keep this notice intact.         *
*                                                                *
\****************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "recycleBin.h"
#include "utility.h"

typedef struct RecycleBin_Node {
	struct RecycleBin_Node *next;
} RecycleBin_Node;

typedef struct chunk_st {
	struct chunk_st *next;
} Chunk;

struct recycleBin_st {
	Chunk **chunk_list;
	size_t chunk_list_size;
	int chunk_index;
	int count;
	RecycleBin_Node *recycle;
	size_t node_size;
	int chunk_pos;
	int nodes_per_chunk;
};

RecycleBin *newRecycleBin(size_t node_size, int nodes_per_chunk)
{
	register RecycleBin *recycle_bin = malloc(sizeof(RecycleBin));
	size_t chunckSize, allocSize;

	if (recycle_bin == NULL) {
		velvetLog("Allocation failed!\n");
#ifdef DEBUG
		abort();
#endif
		exit(-1);
	}

	if (node_size < sizeof(RecycleBin_Node)) {
		velvetLog("Too small elements to create a recycle bin!\n");
#ifdef DEBUG
		abort();
#endif
		exit(-1);
	}
	chunckSize = sizeof(Chunk) + nodes_per_chunk * node_size;
	allocSize = 1;
	/* Get nearest power of 2 */
	while (allocSize < chunckSize)
		allocSize <<= 1;
	nodes_per_chunk = (allocSize - sizeof(Chunk)) / node_size;
	recycle_bin->chunk_list = malloc(8 * sizeof(void *));

	if (recycle_bin->chunk_list == NULL) {
		puts("Could not create recycle bin! Exiting");
		exit(-1);
	}
	recycle_bin->chunk_list_size = 8;
	recycle_bin->chunk_index = 0;
	recycle_bin->chunk_pos = nodes_per_chunk;
	recycle_bin->nodes_per_chunk = nodes_per_chunk;
	recycle_bin->node_size = node_size;
	recycle_bin->recycle = NULL;
	return recycle_bin;
}

void destroyRecycleBin(RecycleBin * recycle_bin)
{
	register int i;

	if (recycle_bin == NULL)
		return;

	for (i = 0; i < recycle_bin->chunk_index; i++)
		free(recycle_bin->chunk_list[i]);
	free(recycle_bin->chunk_list);
	free(recycle_bin);
	return;
}

void *allocatePointer(RecycleBin * recycle_bin)
{
	register RecycleBin_Node *node;
	register Chunk *chunk;
	void * ptr;

	if (recycle_bin == NULL) {
		velvetLog("Null recycle bin!\n");
#ifdef DEBUG
		abort();
#endif 
		exit(-1);
	}
	if (recycle_bin->recycle) {
		ptr = (void *) recycle_bin->recycle;
		recycle_bin->recycle = recycle_bin->recycle->next;
	} else {
		if (recycle_bin->chunk_pos == recycle_bin->nodes_per_chunk) {
			chunk = malloc(recycle_bin->nodes_per_chunk
				 * recycle_bin->node_size);

			if (chunk == NULL) {
				puts("Could not allocate memory, exiting");
				exit(-1);
			}
			if (recycle_bin->chunk_index >= recycle_bin->chunk_list_size) {
				recycle_bin->chunk_list_size *= 2;
				recycle_bin->chunk_list = realloc(recycle_bin->chunk_list, recycle_bin->chunk_list_size * sizeof(void *));
			}
			recycle_bin->chunk_list[recycle_bin->chunk_index++] = (void *) chunk;

			recycle_bin->chunk_pos = 1;
			ptr = (void *) chunk;
		} else {
			chunk = (char *) recycle_bin->chunk_list
				[recycle_bin->chunk_index - 1];
			ptr = (void *) (chunk + (recycle_bin->node_size
					 * recycle_bin->chunk_pos++));
		}
	}

	return ptr;
}

void deallocatePointer(RecycleBin * recycle_bin, void *data)
{
	register RecycleBin_Node *node = data;

	node->next = recycle_bin->recycle;
	recycle_bin->recycle = node;

	return;
}

size_t RecycleBin_memory_usage(RecycleBin * recycle_bin)
{
	return recycle_bin->node_size
	    * recycle_bin->chunk_index * recycle_bin->nodes_per_chunk;
}

size_t recycleBinFreeSpace(RecycleBin * bin)
{
	RecycleBin_Node *freeNode = bin->recycle;
	size_t result = 0;
	while (freeNode != NULL) {
		freeNode = freeNode->next;
		result++;
	}

	return bin->node_size * (result +
				 (bin->nodes_per_chunk - bin->chunk_pos));
}

size_t recycleBinAvailablePointers(RecycleBin * bin)
{
	return bin->chunk_index * bin->nodes_per_chunk;
}
