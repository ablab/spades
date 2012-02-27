/***************************************************************************
 * Title:          free_graph.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

void free_graph(NODES **vertex, int num_vertex);
READPOSITION *free_readposition(READPOSITION *readposition);

void free_graph(NODES **vertex, int num_vertex)
{
	int	i, j, k, l, m;

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			if(vertex[i] -> nextedge[j] -> seq)	{
				free((void *) vertex[i] -> nextedge[j] -> seq);
			}
			if(vertex[i] -> nextedge[j] -> multip > 0)
				free((void *) vertex[i] -> nextedge[j] -> readinterval);
			free((void *) vertex[i] -> nextedge[j]);
		}
		if(vertex[i] -> num_lastedge > 0)	free((void **) vertex[i] -> lastedge);
		if(vertex[i] -> num_nextedge > 0)	free((void **) vertex[i] -> nextedge);
		while(vertex[i] -> readposition)	{
			vertex[i] -> readposition = free_readposition(vertex[i] -> readposition);
		}
		free((void *) vertex[i]);
	}
}

void free_graph_ver(VERTEX **vertex, int num_vertex)
{
	int	i, j, k, l, m;

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			if(vertex[i] -> nextedge[j] -> seq)	{
				free((void *) vertex[i] -> nextedge[j] -> seq);
			}
			free((void *) vertex[i] -> nextedge[j] -> readinterval);
			free((void *) vertex[i] -> nextedge[j]);
		}
		if(vertex[i] -> num_lastedge > 0)	free((void **) vertex[i] -> lastedge);
		if(vertex[i] -> num_nextedge > 0)	free((void **) vertex[i] -> nextedge);
		free((void *) vertex[i] -> path_index);
		free((void *) vertex[i] -> path_pos);
		free((void *) vertex[i]);
	}
}

READPOSITION *free_readposition(READPOSITION *readposition)
{
	READPOSITION *pos;

	pos = readposition -> next;
	free((void *) readposition);
	return(pos);
}
