/***************************************************************************
 * Title:          merge_vertex_path.c
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

EDGE *merge_vertex_path(NODES *vertex, PATH *path, int num_path);

EDGE *merge_vertex_path(NODES *vertex, PATH *path, int num_path)
{
	int	i, j, k, l, m, n, k1, k2, n1, n2, c, q, p, num;
	NODES	*begin, *end, *b0, *e0;
	int	*del_path;
	int	*beginlist, *endlist;
	EDGE	*lastedge, *nextedge, *newedge;
	double	v, r;
	int bne, ble, ene, ele;
	int bne2, ble2, ene2, ele2;
	begin = vertex -> lastedge[0] -> begin;
	end = vertex -> nextedge[0] -> end;
	bne = begin->num_nextedge;
	ble = begin->num_lastedge;
	ene = end->num_nextedge; 
	ele = end->num_lastedge;
	lastedge = vertex -> lastedge[0];
	nextedge = vertex -> nextedge[0];
	if(lastedge == nextedge)	{
		return(NULL);
	}

	del_path = (int *) ckalloc((MAX_BRA + vertex -> num_path) * sizeof(int));
	/* Create a new edge out of lastedge, nextedge. */
	newedge = new_edge(vertex, begin, end, lastedge, nextedge, beginlist, endlist, 0, 0);

	bne2 = begin->num_nextedge;
	ble2 = begin->num_lastedge;
	ene2 = end->num_nextedge; 
	ele2 = end->num_lastedge;
	
	num = vertex -> num_path;
	/* Vertex was bypassed, so we need to do something with the paths
		 that went through this.
	*/
	for(i = 0; i < num; i ++)	{
		j = vertex -> path_index[i];
		k = vertex -> path_pos[i];

		if(path[j].len_path == 0)	continue;
		if(k == 0)	{
			/*
				Newedge is at the beginning of the path. 
				No need to update this.
			*/
			int obl = path[j].begin_length;
			int obr = path[j].edge[k]->length - path[j].begin_length;
			path[j].edge[k] = newedge;
			add_path(newedge -> begin, j, k);
			path[j].begin_length = lastedge->length + path[j].begin_length - VERTEX_SIZE;
			if(vertex != newedge -> begin)	{
				del_path[i] = 1;
			}
		} else if(k == path[j].len_path)	{
			/*
				Newedge is at the end of the path.
			*/
			path[j].edge[k - 1] = newedge;
			add_path(newedge -> end, j, k);
			if(vertex != newedge -> end)	{
				del_path[i] = 1;
			}
		} else	{
			/* 
				 Newedge is in the middle of a path, update this.
			*/
			if(path[j].edge[k - 1] == lastedge && path[j].edge[k] == nextedge)	{
				del_path[i] = 1;
				path[j].edge[k - 1] = newedge;
				remove_edge(&path[j], j, k);
			}
		}
	}

	if(vertex -> num_path != num)	{
		printf("lastedge %d %d %d %d nextedge %d %d %d %d\n", lastedge,
			lastedge -> length, lastedge -> begin, lastedge -> end,
			nextedge, nextedge -> length, nextedge -> begin,
			nextedge -> end);
		printf("num %d %d\n", vertex -> num_path, num);
	}

	n = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		if(del_path[i] == 0)	{
			vertex -> path_index[n] = vertex -> path_index[i];
			vertex -> path_pos[n] = vertex -> path_pos[i];
			n ++;
		}
	}
	vertex -> num_path = n;

	free((void *) del_path);

	n1 = searcherase(begin -> nextedge, lastedge, begin -> num_nextedge);
	n2 = searcherase(end -> lastedge, nextedge, end -> num_lastedge);
	erasenext(begin, n1);
	eraselast(end, n2);

	/*	printf("merging %d prev: %d %d %d %d,  cur: %d %d %d %d, finally: %d %d %d %d\n",
				 vertex->index, bne, ble, ene, ele, 
				 bne2, ble2, ene2, ele2,
				 begin->num_nextedge,
				 begin->num_lastedge,
				 end->num_nextedge ,
				 end->num_lastedge);
	*/

	free((void *) lastedge);
	free((void *) nextedge);
	vertex -> num_nextedge = vertex -> num_lastedge = 0;
/*	printf("vertex %d degree: in: %d out: %d\n", vertex->index, vertex->num_lastedge, vertex->num_nextedge);*/

	return(newedge);
}
