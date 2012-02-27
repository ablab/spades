/***************************************************************************
 * Title:          splitbeg.c
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

int splitbeg(NODES ***vertexListPtr, int num_vertex, PATH *path, int num_path);
NODES *new_vertex(NODES *vertex);
void move_path_next(EDGE *nextedge, NODES *vertex, NODES *vertex0, PATH *path);
void move_path_last(EDGE *lastedge, NODES *vertex, NODES *vertex0, PATH *path);

int splitbeg(NODES ***vertexListPtr, int num_vertex, PATH *path, int num_path)
{
	int	i, j, k, l, m, n, n1, n2, k1, k2, p, q;
	NODES	*begin, *end, *vertex0, *b0, *e0;
	EDGE	**tmpedge, *edge;
	NODES **vertex = *vertexListPtr;
	tmpedge = (EDGE **) ckalloc(MAX_BRA * sizeof(EDGE *));

	i = 0;
	k1 = k2 = 0;
	n1 = n2 = 0;
	while(i < num_vertex)	{
		if(vertex[i] -> num_nextedge > 1 && vertex[i] -> num_lastedge == 0)	{
			k1 ++;
			n1 += vertex[i] -> num_nextedge;
			for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
				vertex0 = new_vertex(vertex[i]);
				move_path_next(vertex[i] -> nextedge[j], vertex[i], vertex0, path);
				*vertexListPtr = (NODES**) realloc(*vertexListPtr, (++num_vertex)*sizeof(NODES*));
				vertex = *vertexListPtr;
				vertex[num_vertex-1] = vertex0;
			}
			free((void **) vertex[i] -> lastedge);
			free((void **) vertex[i] -> nextedge);
			free((void *) vertex[i] -> path_index);
			free((void *) vertex[i] -> path_pos);
			free((void *) vertex[i]);
			for(j = i; j < num_vertex - 1; j ++)	{
				vertex[j] = vertex[j + 1];
			}
			num_vertex --;
		} else if(vertex[i] -> num_lastedge > 1 && vertex[i] -> num_nextedge == 0)     {
			k2 ++;
			n2 += vertex[i] -> num_lastedge;
			for(j = 0; j < vertex[i] -> num_lastedge; j ++)	{
				vertex0 = new_vertex(vertex[i]);
				move_path_last(vertex[i] -> lastedge[j], vertex[i], vertex0, path);
				*vertexListPtr = (NODES**) realloc(*vertexListPtr, (num_vertex+1)*sizeof(NODES*));
				vertex = *vertexListPtr;
				vertex[num_vertex ++] = vertex0;
			}
			free((void **) vertex[i] -> lastedge);
			free((void **) vertex[i] -> nextedge);
			free((void *) vertex[i] -> path_index);
			free((void *) vertex[i] -> path_pos);
			free((void *) vertex[i]);
			for(j = i; j < num_vertex - 1; j ++)	{
				vertex[j] = vertex[j + 1];
			}
			num_vertex --;
		} else	{
			i ++;
		}
	}

	free((void **) tmpedge);

	return(num_vertex);
}


NODES *new_vertex(NODES *vertex)
{
	int	i, j, k, l;
	NODES	*vertex0;

	vertex0 = (NODES *) ckalloc(1 * sizeof(NODES));
	*vertex0 = *vertex;
	vertex0 -> path_index = (int *) ckalloc(vertex0 -> num_path * sizeof(int));
	vertex0 -> path_pos = (int *) ckalloc(vertex0 -> num_path * sizeof(int));
	return(vertex0);
}

void move_path_next(EDGE *nextedge, NODES *vertex, NODES *vertex0, PATH *path)
{
	int	i, j, k, l, n;

	vertex0 -> num_lastedge = 0;
	vertex0 -> num_nextedge = 1;
	vertex0 -> lastedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
	vertex0 -> nextedge = (EDGE **) ckalloc(vertex0 -> num_nextedge * sizeof(EDGE *));
	vertex0 -> nextedge[0] = nextedge;
	nextedge -> begin = vertex0;
	vertex0 -> num_path = 0;
	for(k = vertex -> num_path - 1; k >= 0; k --)	{
		n = vertex -> path_index[k];
		if(path[n].len_path <= 0)	{
			for(j = k; j < vertex -> num_path - 1; j ++)	{
				vertex -> path_index[j] = vertex -> path_index[j + 1];
				vertex -> path_pos[j] = vertex -> path_pos[j + 1];
			}
			vertex -> num_path --;
			continue;
		}
		if(vertex -> path_pos[k] == 0 && path[n].edge[0] == nextedge)	{
			add_path(vertex0, n, 0);
			for(j = k; j < vertex -> num_path - 1; j ++)	{
				vertex -> path_index[j] = vertex -> path_index[j + 1];
				vertex -> path_pos[j] = vertex -> path_pos[j + 1];
			}
			vertex -> num_path --;
		}
	}
}

void move_path_last(EDGE *lastedge, NODES *vertex, NODES *vertex0, PATH *path)
{
	int	i, j, k, l, n;

	vertex0 -> num_lastedge = 1;
	vertex0 -> num_nextedge = 0;
	vertex0 -> lastedge = (EDGE **) ckalloc(vertex0 -> num_lastedge * sizeof(EDGE *));
	vertex0 -> nextedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
	vertex0 -> lastedge[0] = lastedge;
	lastedge -> end = vertex0;
	vertex0 -> num_path = 0;
	for(k = vertex -> num_path - 1; k >= 0; k --)	{
		n = vertex -> path_index[k];
		if(path[n].len_path <= 0)	{
			for(j = k; j < vertex -> num_path - 1; j ++)	{
				vertex -> path_index[j] = vertex -> path_index[j + 1];
				vertex -> path_pos[j] = vertex -> path_pos[j + 1];
			}
			vertex -> num_path --;
			continue;
		}
		if(vertex -> path_pos[k] == path[n].len_path && path[n].edge[path[n].len_path - 1] == lastedge)	{
			add_path(vertex0, n, path[n].len_path);
			for(j = k; j < vertex -> num_path - 1; j ++)	{
				vertex -> path_index[j] = vertex -> path_index[j + 1];
				vertex -> path_pos[j] = vertex -> path_pos[j + 1];
			}
			vertex -> num_path --;
		}
	}
}
