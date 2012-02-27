/***************************************************************************
 * Title:          newvertex.c
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

int newvertex(NODES ***vertexListPtr, int num_vertex, EDGE *edge, PATH *path);

int newvertex(NODES ***vertexListPtr, int num_vertex, EDGE *edge, PATH *path)
{
	int	i, j, k, l, m, n, n1, n2, c, len;
	int	*label;
	double	r;
	NODES	*begin, *end, *vertex_new;
	EDGE	*edge1, *edge2;
	NODES **vertex = *vertexListPtr;
	begin = edge -> begin;
	end = edge -> end;
	n1 = 0;
	n2 = 0;
	for(i = 0; i < begin -> num_path; i ++)	{
		j = begin -> path_index[i];
		k = begin -> path_pos[i];
		if(k < path[j].len_path && path[j].edge[k] == edge)	{
			if(k > 0 || begin -> num_lastedge == 0)	{
				n1 ++;
			}
			if(k < path[j].len_path - 1 || end -> num_nextedge == 0)	{
				n2 ++;
			}
		}
	}

	if(n1 == 0 && begin -> num_lastedge > 0)	{
		vertex_new = (NODES *) ckalloc(1 * sizeof(NODES));
		*vertex_new = *begin;
		edge -> begin = vertex_new;
		vertex_new -> num_nextedge = 1;
		vertex_new -> num_lastedge = 0;
		vertex_new -> lastedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		vertex_new -> nextedge = (EDGE **) ckalloc(vertex_new -> num_nextedge * sizeof(EDGE *));
		vertex_new -> nextedge[0] = edge;
		vertex_new -> num_path = 0;
		label = (int *) ckalloc(begin -> num_path * sizeof(int));
		vertex_new -> path_index = (int *) ckalloc(begin -> num_path * sizeof(int));
		vertex_new -> path_pos = (int *) ckalloc(begin -> num_path * sizeof(int));
		for(i = 0; i < begin -> num_path; i ++)	{
			j = begin -> path_index[i];
			k = begin -> path_pos[i];
			if(k < path[j].len_path && path[j].edge[k] == edge)	{
				vertex_new -> path_index[vertex_new -> num_path] = j;
				vertex_new -> path_pos[vertex_new -> num_path ++] = k;
				label[i] = 1;
			}
		}
		for(i = begin -> num_path - 1; i >= 0; i --)	{
			if(label[i] == 1)	{
				for(j = i; j < begin -> num_path - 1; j ++)	{
					begin -> path_index[j] = begin -> path_index[j + 1];
					begin -> path_pos[j] = begin -> path_pos[j + 1];
				}
				begin -> num_path --;
			}
		}
		*vertexListPtr = (NODES**) realloc(*vertexListPtr, num_vertex+1);
		vertex = *vertexListPtr;
		vertex[num_vertex ++] = vertex_new;
		n = searcherase(begin -> nextedge, edge, begin -> num_nextedge);
		erasenext(begin, n);
		free((void *) label);
	}

	if(n2 == 0 && end -> num_nextedge > 0)	{
		vertex_new = (NODES *) ckalloc(1 * sizeof(NODES));
		*vertex_new = *end;
		edge -> end = vertex_new;
		vertex_new -> num_nextedge = 0;
		vertex_new -> num_lastedge = 1;
		vertex_new -> lastedge = (EDGE **) ckalloc(vertex_new -> num_lastedge * sizeof(EDGE *));
		vertex_new -> nextedge = (EDGE **) ckalloc(1 * sizeof(EDGE *));
		vertex_new -> lastedge[0] = edge;
		vertex_new -> num_path = 0;
		label = (int *) ckalloc(end -> num_path * sizeof(int));
		vertex_new -> path_index = (int *) ckalloc(end -> num_path * sizeof(int));
		vertex_new -> path_pos = (int *) ckalloc(end -> num_path * sizeof(int));
		for(i = 0; i < end -> num_path; i ++)	{
			j = end -> path_index[i];
			k = end -> path_pos[i];
			if(k > 0 && path[j].edge[k - 1] == edge)	{
				vertex_new -> path_index[vertex_new -> num_path] = j;
				vertex_new -> path_pos[vertex_new -> num_path ++] = k;
				label[i] = 1;
			}
		}
		for(i = end -> num_path - 1; i >= 0; i --)	{
			if(label[i] == 1)	{
				for(j = i; j < end -> num_path - 1; j ++)	{
					end -> path_index[j] = end -> path_index[j + 1];
					end -> path_pos[j] = end -> path_pos[j + 1];
				}
				end -> num_path --;
			}
		}
		
		vertex[num_vertex ++] = vertex_new;
		n = searcherase(end -> lastedge, edge, end -> num_lastedge);
		eraselast(end, n);
		free((void *) label);
	}
	*vertexListPtr = (NODES**) realloc(*vertexListPtr, num_vertex+1);
	vertex = *vertexListPtr;
	return(num_vertex);
}
