/***************************************************************************
 * Title:          merge_graph_path.c
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

int merge_graph_path(NODES **vertex, int num_vertex, PATH *path, int num_path, int small_cover);
int merge_graph_path(NODES **vertex, int num_vertex, PATH *path, int num_path, int small_cover)
{
	int	i, j, k, l, m, n, k1, k2, n1, n2, q, p, len, label, c;
	int	num_vertex_old;
	NODES	*begin, *end, *vertex0;
	EDGE	*lastedge, *nextedge, *edge, *edge1, *edge2, *edge01, *edge02;
	double	v, r;
	
	c = 0;
	do	{
		num_vertex_old = num_vertex;
		i = 0;
		while(i < num_vertex)	{
		  /* Remove 0-degree vertices */
			if(vertex[i] -> num_nextedge == 0 && vertex[i] -> num_lastedge == 0)	{
				free((void **) vertex[i] -> lastedge);
				free((void **) vertex[i] -> nextedge);
				free((void *) vertex[i] -> path_index);
				free((void *) vertex[i] -> path_pos);
				free((void *) vertex[i]);
				for(j = i; j < num_vertex - 1; j ++)	{
					vertex[j] = vertex[j + 1];
				}
				num_vertex --;
			} 
			/*
			  vertex[i] is along a simple path.  Remove it and merge paths from lastedge to nextedge
			*/
			else if(vertex[i] -> num_nextedge == 1 && vertex[i] -> num_lastedge == 1)	{
				vertex0 = vertex[i] -> lastedge[0] -> bal_edge -> begin;
				lastedge = vertex[i] -> lastedge[0];
				nextedge = vertex[i] -> nextedge[0];
				if(lastedge == nextedge)	{
					i ++;
				} else	{
					edge01 = lastedge -> bal_edge;
					edge02 = nextedge -> bal_edge;
					edge1 = merge_vertex_path(vertex[i], path, num_path);
					if(vertex0 == vertex[i])	{
						edge1 -> bal_edge = edge1;
					} else if(vertex0 -> lastedge[0] != vertex0 -> nextedge[0])	{
						edge2 = merge_vertex_path(vertex0, path, num_path);
						if(edge01 == lastedge || edge02 == nextedge)	{
							edge2 -> bal_edge = edge2;
						} else	{
							edge1 -> bal_edge = edge2;
							edge2 -> bal_edge = edge1;
						}
					} else	{
						edge1 -> bal_edge = edge1;
					}
					i ++;
				}
			} else	{
				c = rm_edge(vertex[i], path, num_path);
				if(!c)	{
					i ++;
				}
			}
		}
	} while(num_vertex < num_vertex_old);

	return(num_vertex);
}
