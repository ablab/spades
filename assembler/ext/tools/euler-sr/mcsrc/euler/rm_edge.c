/***************************************************************************
 * Title:          rm_edge.c
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

int rm_edge(NODES *vertex, PATH *path, int num_path);
void replace1edge(PATH *path, int num_path, EDGE *edge1, EDGE *edge2);

int rm_edge(NODES *vertex, PATH *path, int num_path)
{

  /* 
     vertex: a vertex to check to see if it has parallel edges from it.
     path: a list of all num_path paths.

     Result: if 'vertex' is the dest of parallel edges, or the source of 
     parallel edges, the edge of lower multiplicity is removed, and all paths
     are routed through 'vertex'.
  */
       

	int	i, j, k, l, m, n, k1, k2, n1, n2, c, q, p, label;
	int	num_vertex_old;
	NODES	*begin, *end, *v0;
	EDGE	*lastedge, *nextedge, *edge1, *edge2;
	double	v, r;

	label = 0;

	j = 0;
	while(j < vertex -> num_lastedge - 1)	{
		k = j + 1;
		while(k < vertex -> num_lastedge)	{
		  /* 
		     If lastedge[j] has the same source as lastedge[k], and the edges are about the same length,
		     we want to merge them together.  Remove the lower multiplicity edge.
		  */
			if(vertex -> lastedge[j] -> begin == vertex -> lastedge[k] -> begin)	{
				if(abs(vertex -> lastedge[j] -> length - vertex -> lastedge[k] -> length) < MIN_INT)	{
					if(vertex -> lastedge[j] -> multip < vertex -> lastedge[k] -> multip)	{
						lastedge = vertex -> lastedge[j];
						vertex -> lastedge[j] = vertex -> lastedge[k];
						vertex -> lastedge[k] = lastedge;
					}
					lastedge = vertex -> lastedge[j];
					nextedge = vertex -> lastedge[k];
					edge1 = nextedge -> bal_edge;
					edge2 = lastedge -> bal_edge;
					if(edge1 != lastedge && edge2 != lastedge && edge1 != nextedge && edge2 != nextedge)	{
						/*						replace1edge(path, num_path, nextedge, lastedge);
						replace1edge(path, num_path, edge1, edge2);
						*/
						k++;
						//						label = 1;
						continue;
					}
				}
			}
			k ++;
		}
		j ++;
	}

	/* 
	   The same thing as the first part, except remove the out edges that have the same source. 
	   I'm not sure why this is necessary, since if two edge have the same source and same dest, they 
	   will eventually be removed by the first part.
	*/
	j = 0;
	while(j < vertex -> num_nextedge - 1)	{
		k = j + 1;
		while(k < vertex -> num_nextedge)	{
			if(vertex -> nextedge[j] -> end == vertex -> nextedge[k] -> end)	{
				if(abs(vertex -> nextedge[j] -> length - vertex -> nextedge[k] -> length) < MIN_INT)	{
					if(vertex -> nextedge[j] -> multip < vertex -> nextedge[k] -> multip)	{
						nextedge = vertex -> nextedge[j];
						vertex -> nextedge[j] = vertex -> nextedge[k];
						vertex -> nextedge[k] = nextedge;
					}
					lastedge = vertex -> nextedge[j];
					nextedge = vertex -> nextedge[k];
					edge1 = nextedge -> bal_edge;
					edge2 = lastedge -> bal_edge;
					if(edge1 != lastedge && edge2 != lastedge && edge1 != nextedge && edge2 != nextedge)	{
						/*
						replace1edge(path, num_path, nextedge, lastedge);
						replace1edge(path, num_path, edge1, edge2);
						*/
						k++;
						//						label = 1;
						continue;
					}
				}
			}
			k ++;
		}
		j ++;
	}

	return(label);
}

void replace1edge(PATH *path, int num_path, EDGE *edge1, EDGE *edge2)
{
	int	n;
	NODES	*vertex;

	if(edge1 -> begin == edge2 -> begin)	{
		vertex = edge1 -> begin;
	} else if(edge1 -> end == edge2 -> end)	{
		vertex = edge1 -> end;
	}
	replacepath(path, num_path, vertex, edge1, edge2);
	movereadinterval(edge2, edge1);
	n = searcherase(edge1 -> begin -> nextedge, edge1, edge1 -> begin -> num_nextedge);
	erasenext(edge1 -> begin, n);
	n = searcherase(edge1 -> end -> lastedge, edge1, edge1 -> end -> num_lastedge);
	eraselast(edge1 -> end, n);
	free((void *) edge1);
}
