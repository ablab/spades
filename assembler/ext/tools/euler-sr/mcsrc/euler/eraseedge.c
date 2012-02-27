/***************************************************************************
 * Title:          eraseedge.c
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
#include <assert.h>

void erasedge(EDGE *edge);
int crossedge(EDGE *edge, PATH *path);


void RemovePath(PATH *paths, int numPaths, int path) {
	if (path >= numPaths) {
		printf("ERROR! Referencing a path that does not exist\n");
		assert(0);
	}
	// Remove the path from the vertices that reference it
	int pos;
	for (pos = 0; pos < paths[path].len_path; pos++) {
		RemovePathFromVertex(paths[path].edge[pos]->begin, path);
		remove_readinterval(paths[path].edge[pos], path);
	}
	RemovePathFromVertex(paths[path].edge[pos-1]->end, path);
	paths[path].len_path = 0;
}

void MarkEdgeForRemoval(EDGE *edge) {
	edge->deleted = 1;
	edge->bal_edge->deleted = 1;
}

void RemovePathsWithRemovedEdges(PATH *paths, int numPaths ) {
	int p, i;
	for (p = 0; p < numPaths; p++) {
		for (i = 0; i < paths[p].len_path; i++) {
			if (paths[p].edge[i]->deleted == 1) {
				RemovePath(paths, numPaths, p); 
			}
		}
	}
}

void EraseMarkedEdgesAndPaths(NODES **vertices, int numVertices,
														 PATH *paths, int numPaths) {
	RemovePathsWithRemovedEdges(paths, numPaths);
	EraseEdgesMarkedForDeletion(vertices, numVertices);
}

void EraseEdgesMarkedForDeletion(NODES **vertices, int numVertices) {
	int v, e;
	/* First do a sanity check to make sure all edges 
		 have balanced edges that are also to be erased.*/
	for (v = 0; v < numVertices; v++ ){ 
		e = 0;
		for (e = 0; e < vertices[v]->num_nextedge; e++) {
			if (vertices[v]->nextedge[e]->deleted == 1 and
					vertices[v]->nextedge[e]->bal_edge->deleted != 1) {
				printf("ERROR, edge: %d is marked for deletion, but not: %d\n",
							 vertices[v]->nextedge[e]->index,
							 vertices[v]->nextedge[e]->bal_edge->index);
				exit(1);
			}
		}
	}
		

	for (v = 0; v < numVertices; v++ ){ 
		e = 0;
		while (e < vertices[v]->num_nextedge) {
			if (vertices[v]->nextedge[e]->deleted == 1) {
				/*				printf("erasing edge: %d\n", vertices[v]->nextedge[e]->index);*/
				erasedge(vertices[v]->nextedge[e]);
			}
			else {
				e++;
			}
		}
	}
}



void EraseEdgeAndPassingPaths(EDGE *edge, PATH *paths, int numPaths) {
	int prevMultip;
	while (edge->multip > 0) {
		prevMultip = edge->multip;
		RemovePath(paths, numPaths, edge->readinterval[0].eq_read);
		assert(prevMultip > edge->multip);
	}
	erasedge(edge);
}

void erasedge(EDGE *edge)
{
	int	i, j, k, l, n, m;
	NODES	*begin, *end;

	begin = edge -> begin;
	end = edge -> end;
	k = searcherase(begin -> nextedge, edge, begin -> num_nextedge);
	erasenext(begin, k);
	k = searcherase(end -> lastedge, edge, end -> num_lastedge);
	eraselast(end, k);

	/*free((void *) edge -> readinterval);
		free((void *) edge);*/
	edge->deleted = 1;
}

int crossedge(EDGE *edge, PATH *path)
{
	int	i, j, k, l, n1, n2;
	NODES	*begin, *end;

	n1 = 0;
	n2 = 0;
	begin = edge -> begin;
	for(i = 0; i < begin -> num_path; i ++)	{
		j = begin -> path_index[i];
		k = begin -> path_pos[i];
		if(k < path[j].len_path && path[j].edge[k] == edge)	{
			if(k > 0 || begin -> num_lastedge == 0)	{
				n1 ++;
			}
			if(k < path[j].len_path - 1 || edge -> end -> num_nextedge == 0)	{
				n2 ++;
			}
		}
	}

	if(n1 == 0 || n2 == 0)	{
		return(0);
	} else	{
		return(1);
	}
}
