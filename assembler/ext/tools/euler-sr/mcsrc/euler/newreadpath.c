/***************************************************************************
 * Title:          newreadpath.c
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

int newreadpath(NODES **vertex, int num_vertex, PATH *path, int *len_seq, int num_seq);

int newreadpath(NODES **vertex, int num_vertex, PATH *path, int *len_seq, int num_seq)
{
	int	i, j, k, l, m, n;
	int	num_path;
	int	nch;
	int	reads;
	char	*label;
	EDGE	*edge, *bal_edge;
	NODES	*v;

	label = (char *) ckalloc(2 * num_seq * sizeof(char));

/*	Define read paths	*/

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			findpath(path, edge, len_seq, num_seq, label);
		}
	}

/*	Reorder paths by removing single read edges paths with length 0 (reads in a single edge)	*/

	num_path = 0;
	for(i = 0; i < num_seq * 2; i ++)	{
		if(path[i].len_path > 1)	{
			path[num_path ++] = path[i];
		}
	}

/*	Set up path index & positions for each vertices	*/

	for(i = 0; i < num_path; i ++)	{
		for(j = 0; j < path[i].len_path; j ++)	{
			v = path[i].edge[j] -> begin;
			v -> num_path ++;
		}
		v = path[i].edge[j - 1] -> end;
		v -> num_path ++;
	}

	for(i = 0; i < num_vertex; i ++)	{
		vertex[i] -> path_pos = (int *) ckalloc(vertex[i] -> num_path * sizeof(int));
		vertex[i] -> path_index = (int *) ckalloc(vertex[i] -> num_path * sizeof(int));
		vertex[i] -> num_path = 0;
	}

	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path <= 1)	continue;
		for(j = 0; j < path[i].len_path; j ++)	{
			v = path[i].edge[j] -> begin;
			v -> path_index[v -> num_path] = i;
			v -> path_pos[v -> num_path ++] = j;
		}
		v = path[i].edge[j - 1] -> end;
		v -> path_index[v -> num_path] = i;
		v -> path_pos[v -> num_path ++] = j;
	}
	free((void *) label);
	return(num_path);
}
