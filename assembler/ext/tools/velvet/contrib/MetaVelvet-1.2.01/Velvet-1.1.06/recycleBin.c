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
	Chunk *chunk_list;
	RecycleBin_Node *recycle;
	size_t node_size;
	int chunk_pos;
	int nodes_per_chunk;
};

static void initRecycleBin(RecycleBin *recycleBin,
			   size_t node_size, int nodes_per_chunk)
{
	size_t chunckSize, allocSize;

	chunckSize = sizeof(Chunk) + nodes_per_chunk * node_size;
	allocSize = 1;
	/* Get nearest power of 2 */
	while (allocSize < chunckSize)
		allocSize <<= 1;
	nodes_per_chunk = (allocSize - sizeof(Chunk)) / node_size;
	recycleBin->chunk_list = NULL;
	recycleBin->chunk_pos = nodes_per_chunk;
	recycleBin->nodes_per_chunk = nodes_per_chunk;
	recycleBin->node_size = node_size;
	recycleBin->recycle = NULL;
}

RecycleBin *newRecycleBin(size_t node_size, int nodes_per_chunk)
{
	RecycleBin *recycleBin;

	if (node_size < sizeof(RecycleBin_Node)) {
		velvetLog("Too small elements to create a recycle bin!\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(-1);
	}
	recycleBin = mallocOrExit(1, RecycleBin);
	initRecycleBin (recycleBin, node_size, nodes_per_chunk);

	return recycleBin;
}

static void destroyRecycleBinChunks(RecycleBin * recycleBin)
{
	while (recycleBin->chunk_list != NULL)
	{
		Chunk *chunk;

		chunk = recycleBin->chunk_list;
		recycleBin->chunk_list = recycleBin->chunk_list->next;
		free(chunk);
	}
}

void destroyRecycleBin(RecycleBin * recycleBin)
{
	if (recycleBin == NULL)
		return;

	destroyRecycleBinChunks(recycleBin);
	free(recycleBin);
}

void *allocatePointer(RecycleBin * recycle_bin)
{
	RecycleBin_Node *node;
	Chunk *chunk;

	if (recycle_bin == NULL) {
		velvetLog("Null recycle bin!\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(-1);
	}

	if (recycle_bin->recycle != NULL) {
		node = recycle_bin->recycle;
		recycle_bin->recycle = node->next;
		return node;
	}

	if (recycle_bin->chunk_pos == recycle_bin->nodes_per_chunk) {
	  chunk = (Chunk*)malloc(sizeof(Chunk) + recycle_bin->nodes_per_chunk
				 * recycle_bin->node_size);
		if (chunk == NULL) {
			velvetLog("No more memory for memory chunk!\n");
#ifdef DEBUG 
			abort();
#endif 
			exit(-1);
		}
		chunk->next = recycle_bin->chunk_list;
		recycle_bin->chunk_list = chunk;
		recycle_bin->chunk_pos = 1;
		return (RecycleBin_Node *) ((size_t) (void *) chunk +
					    sizeof(Chunk));
	}

	chunk = recycle_bin->chunk_list;
	return (RecycleBin_Node *) ((size_t) (void *) chunk + sizeof(Chunk)
				    +
				    (recycle_bin->
				     node_size *
				     recycle_bin->chunk_pos++));
}

void deallocatePointer(RecycleBin * recycle_bin, void *data)
{
  RecycleBin_Node *node = (RecycleBin_Node*)data;

  node->next = recycle_bin->recycle;
  recycle_bin->recycle = node;
  
  return;
}

#ifdef _OPENMP
RecycleBin *newRecycleBinArray(unsigned int n,
			       size_t node_size, int nodes_per_chunk)
{
	RecycleBin *recycleBin;
	int i;

	recycleBin = mallocOrExit (n + 1, RecycleBin);
	for (i = 0; i < n; i++)
		initRecycleBin(recycleBin + i, node_size, nodes_per_chunk);
	/* Last element marker */
	recycleBin[n].node_size = 0;

	return recycleBin;
}

void destroyRecycleBinArray(RecycleBin * recycleBin)
{
	int i;

	for (i = 0; recycleBin[i].node_size != 0; i++)
		destroyRecycleBinChunks(recycleBin + i);
	free(recycleBin);
}

RecycleBin *getRecycleBinInArray(RecycleBin *recycleBin,
				 int	     position)
{
	return recycleBin + position;
}
#endif
