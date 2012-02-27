/***************************************************************************
 * Title:          locpath.c
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
#include <set>
#define MAX_COPY 3

int locpath(EDGE *edge1, EDGE *edge2, int min_leg, int max_leg, int len, PATH *path,
						std::set<EDGE*> &traversedEdges);

int locpath(EDGE *edge1, EDGE *edge2, int min_leg, int max_leg, int len, PATH *path,
						std::set<EDGE*> &traversedEdges)
{
  // count the number of paths from edge1 to edge2 with a dfs
	int	i, j, k, l, num_path, n;
	NODES	*vertex;
	EDGE	*edge;
	PATH	*path_temp;

	if(len > max_leg)	{
	  return(0);
	}

	vertex = edge1 -> end;

	path_temp = (PATH *) ckalloc(vertex -> num_nextedge * sizeof(PATH));
	for(i = 0; i < vertex -> num_nextedge; i ++)	{
		path_temp[i].edge = (EDGE **) ckalloc(max_leg * sizeof(EDGE *));
		path_temp[i].readindex = 0;
	}
	num_path = 0;
	int validPathEdgeIndex;
	for(i = 0; i < vertex -> num_nextedge; i ++)	{
	  
		edge = vertex -> nextedge[i];
		traversedEdges.insert(edge);
		// If this edge has been traversed too many times, continue
		if(edge -> start_cover > MAX_COPY)	continue;

		// mark this edge as having been traversed in this search.
		edge -> start_cover ++;
		
		// path_temp represents paths starting from the end of edge1
		path_temp[i].edge[0] = edge;
		path_temp[i].len_path = 1;
		
		// if this edge represents a valid path to edge2, count it as a valid path
		// and reference this path from edge1 as 'k'
		if(edge == edge2 && len >= min_leg)	{
			num_path ++;
			validPathEdgeIndex = i;
		}
		// keep looking from this edge
		n = locpath(edge, edge2, min_leg, 
								max_leg, len + edge -> length - VERTEX_SIZE, 
								&path_temp[i], traversedEdges);
		if(n == 1)	{
			validPathEdgeIndex = i;
		}
		num_path += n;
		if(num_path > 1)	{
			break;
		}
	}

	if(num_path == 1)	{
		for(i = 0; i < path_temp[validPathEdgeIndex].len_path; i ++)	{
			path -> edge[path -> len_path + i] = path_temp[validPathEdgeIndex].edge[i];
		}
		path -> len_path += path_temp[validPathEdgeIndex].len_path;
	}

	for(i = 0; i < vertex -> num_nextedge; i ++)	{
		free((void **) path_temp[i].edge);
	}
	free((void *) path_temp);
	return(num_path);
}


void ResetTraversedEdges(std::set<EDGE*> &traversedEdges) {
	std::set<EDGE*>::iterator setIt;
	for (setIt = traversedEdges.begin(); setIt != traversedEdges.end(); ++setIt) {
		(*setIt)->start_cover = 0;
	}
}
