/***************************************************************************
 * Title:          path.c
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
void replacepath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *newedge);
int count_path(PATH *path, int num_path);
int cutbegpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *nextedge);
int cutendpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *lastedge);
void remove_edge(PATH *path, int path_index, int path_pos);
int chk_path(int p, int *vp, int nvp);

int cutbegpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *nextedge)
{
	int	i, j, k, n, m, l, label, *pos;

	pos = (int *) ckalloc(num_path * sizeof(int));
	for(i = 0; i < vertex -> num_path; i ++)	{
		pos[i] = MAX_BRA;
	}

	m = 0;
	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path <= 0)	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
			continue;
		}
		if(j < path[k].len_path && path[k].edge[j] == nextedge)	{
			if(j > m)	{
				m = j;
			}
			for(l = 0; l < m; l ++)	{
				remove_edge(&path[k], k, 0);
			}
			pos[k] = j;
		}
	}
	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j < pos[k])	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
		}
	}

	free((void *) pos);

	return(m);
}

int cutendpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *lastedge)
{
	int	i, j, k, n, m, label, *pos;

	pos = (int *) ckalloc(num_path * sizeof(int));
	for(i = 0; i < vertex -> num_path; i ++)	{
		pos[i] = MAX_BRA;
	}

	m = 0;
	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path <= 0)	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
			continue;
		}
		if(j > 0 && path[k].edge[j - 1] == lastedge)	{
			if(path[k].len_path - j > m)	{
				m = path[k].len_path - j;
			}
			path[k].len_path = j;
			pos[k] = j;
		}
	}
	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j > pos[k])	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
		}
	}

	free((void *) pos);

	return(m);
}


int chk_path(int p, int *vp, int nvp)
{
	int	i;

	for(i = 0; i < nvp; i ++)	{
		if(vp[i] == p)	{
			break;
		}
	}

	if(i == nvp)	{
		return(0);
	} else	{
		return(1);
	}
}

void replacepath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *newedge)
{
	int	i, j, k, l, n, m;
	NODES	*begin;

	for(i = vertex -> num_path - 1; i >= 0; i --)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path <= 0)	{
			for(n = i; n < vertex -> num_path - 1; n ++)	{
				vertex -> path_index[n] = vertex -> path_index[n + 1];
				vertex -> path_pos[n] = vertex -> path_pos[n + 1];
			}
			vertex -> num_path --;
			continue;
		}
		if(j < path[k].len_path)	{
			if(path[k].edge[j] == edge)	{
				path[k].edge[j] = newedge;
			}
		}
	}
}

void remove_edge(PATH *path, int path_index, int path_pos)
{
	int	i, j, k, l, n;
	int	**pos;
	NODES	*vertex, **vertex_temp;

	/* 
	   remove the edge at path_pos from path.
	*/
	int *lengths;
	vertex_temp = (NODES **) ckalloc((path -> len_path - 1) * sizeof(NODES *));
	lengths = (int*) ckalloc(path->len_path * sizeof(int));
	pos = (int **) ckalloc((path -> len_path - path_pos - 1) * sizeof(int *));
	/*	printf("allocating path: ");*/
	for(i = path_pos; i < path -> len_path -1; i ++)	{ 
	  if (path->edge[i]->begin->num_path > 0) {
	    /*	    printf(" %d", path->edge[i]->begin->num_path); */
	    pos[i-path_pos] = (int *) ckalloc(path->edge[i+1]->begin->num_path * sizeof(int));
	    lengths[i-path_pos] = path->edge[i+1]->begin->num_path;
	  }
	}
	/*
	  printf("\n");
	*/

	/*
	  vertex temp will hold the 'begin' vertex of the 
	  edges after path_pos.
	*/
	for(i = path_pos; i < path -> len_path - 1; i ++)	{	
	        /* move all paths back by 1*/
		path -> edge[i] = path -> edge[i + 1];
		vertex = vertex_temp[i - path_pos] = path -> edge[i] -> begin;
		/* for each path that goes through this vertex */
		for(j = 0; j < vertex -> num_path; j ++)	{
		  /* If the path 'j' is the one being deleted AND
		   */
		  assert(j < lengths[i - path_pos]);
			if(vertex -> path_index[j] == path_index && vertex -> path_pos[j] == i + 1)	{
				pos[i - path_pos][j] = i;
			} else	{
				pos[i - path_pos][j] = -1;
			}
		}
	}

	vertex = path -> edge[path -> len_path - 1] -> end;
	for(j = vertex -> num_path - 1; j >= 0; j --)	{
		if(vertex -> path_index[j] == path_index && vertex -> path_pos[j] == path -> len_path)	{
			vertex -> path_pos[j] = path -> len_path - 1;
			break;
		}
	}

	path -> len_path -= 1;
	
	for(i = path_pos; i < path -> len_path; i ++)	{
		vertex = vertex_temp[i - path_pos];
		for(j = 0; j < vertex -> num_path; j ++)	{
		  assert(j < lengths[i - path_pos]);
			if(pos[i - path_pos][j] >= 0)	{
				vertex -> path_pos[j] = pos[i - path_pos][j];
			}
		}
	}

	for(i = 0; i < path -> len_path - path_pos ; i ++)	{
	  if (lengths[i] > 0)
	    free( pos[i]);
	}
	if (path->len_path > 0) {
	  free( pos);
	  free(lengths);
	  free((void **) vertex_temp);
	}
	/*
	  printf("removed a path\n");*/

}



void RemoveUninformativePaths(PATH *path, int num_path) {
	int p, i, r, e;
	int numRemoved = 0;
	NODES *source, *dest;
	for (p = 0; p < num_path; p++) {
		if (path[p].len_path == 2) {
			/*			printf("removing path: %d\n", p);*/
			// remove this path from the edges 
			// and vertices that reference it
			int numDel;
			int eq_read = path[p].readindex - 1;
			for (e = 0; e < path[p].len_path; e++) {
				numDel = 0;
				source = path[p].edge[e]->begin;
				dest   = path[p].edge[e]->end;
				RemovePathIntervalFromVertex(source, p, e);
				source->num_path -= numDel;
				remove_readinterval(path[p].edge[e], eq_read);
			}
			RemovePathIntervalFromVertex(dest, p, path[p].len_path);
			path[p].len_path = 0;
			numRemoved++;
		}
	}
	printf("removed %d uninformative paths\n", numRemoved);
}

void RemovePath(PATH *paths, int pathIndex) {
	int i;
	for (i = 0; i < paths[pathIndex].len_path; i++) {
		RemovePathFromVertex(paths[pathIndex].edge[i]->begin, pathIndex);
	}
	RemovePathFromVertex(paths[pathIndex].edge[i-1]->end, pathIndex);
	paths[pathIndex].len_path = 0;
}

void RemovePathsFromVertex(NODES *vertex, int *pathIndices, int *pathPos, int numPaths) {
	int i;
	// The vertex paths aren't necessarily in order, so do n^2 work
	// for now.
	for (i = 0; i < numPaths; i++ ) {
		RemovePathIntervalFromVertex(vertex, pathIndices[i], pathPos[i]);
	}
}

void RemovePathFromVertex(NODES *vertex, int path) {
	// This differs from the next function because
	// all intervals of that path are removed from this vertex
	int r;
	int numDel = 0;
	for (r = 0; r < vertex->num_path; r++) {
		if (vertex->path_index[r] == path) {
			numDel++;
		}
		else {
			vertex->path_index[r-numDel] = vertex->path_index[r];
			vertex->path_pos[r-numDel]   = vertex->path_pos[r];
		}
	}
	vertex->num_path -= numDel;
}

int RemovePathIntervalFromVertex(NODES *vertex, int path, int pos) {
	int r;
	int numDel = 0;
	for (r = 0; r < vertex->num_path; r++) {
		if (vertex->path_index[r] == path and
				vertex->path_pos[r] == pos) {
			numDel++;
		}
		else {
			vertex->path_index[r-numDel] = vertex->path_index[r];
			vertex->path_pos[r-numDel]   = vertex->path_pos[r];
		}
	}
	vertex->num_path -= numDel;
	return numDel;
}

void UpdatePathIndices(NODES *vertex, int path, int pos, int offset) {
	int p;
	for (p = 0; p < vertex->num_path; p++) { 
		if (vertex->path_index[p] == path and
				vertex->path_pos[p] == pos) {
			vertex->path_pos[p] += offset;
			break;
		}
	}
}
void UpdatePathIndices(NODES *vertex, int path, int offset) {
	int p;
	for (p = 0; p < vertex->num_path; p++) { 
		if (vertex->path_index[p] == path)
			vertex->path_pos[p] += offset;
	}
}

int GetMaxPathLength(PATH *paths, int numPaths) {
	int p, maxLength;
	maxLength = 0;
	for (p = 0; p < numPaths; p++) {
		if (paths[p].len_path > maxLength) maxLength = paths[p].len_path;
	}
	return maxLength;
}


void TrimShortEnds(PATH *paths, int numPaths, int minPathLength, int maxEdgeLength) {
	int p;
	int toRemove[2];
	int numToRemove = 0;
	int numTrimmed = 0;
	for (p = 0; p < numPaths; p++) {
		if (paths[p].len_path >= minPathLength ) {
			/*			printf("path %d beginedge: %d end: %d length: %d\n",
						 p, paths[p].edge[0]->length,
						 paths[p].edge[paths[p].len_path-1]->length,
						 paths[p].len_path);*/
			if (paths[p].edge[paths[p].len_path-1]->length < maxEdgeLength) {
				/* The path ends in a short edge, truncate that.*/
				/*				printf("trimming %d end\n", p);*/
				RemovePathIntervalFromVertex(paths[p].edge[paths[p].len_path-1]->end, p, paths[p].len_path);
				remove_readinterval(paths[p].edge[paths[p].len_path-1], paths[p].readindex-1);
				paths[p].len_path--;
				paths[p].end_length = 0;
				++numTrimmed;
			}
			if (paths[p].edge[0]->length < maxEdgeLength) {
				/*				printf("trimming %d begin\n", p);*/
				/* The path starts in a short edge, truncate that interval.*/
				RemovePathIntervalFromVertex(paths[p].edge[0]->begin, p, 0);
				remove_readinterval(paths[p].edge[0], paths[p].readindex-1);
				int r;
				for (r = 1; r < paths[p].len_path; r++) {
					paths[p].edge[r-1] = paths[p].edge[r];
					UpdatePathIndices(paths[p].edge[r-1]->begin, p, -1);
				}
				UpdatePathIndices(paths[p].edge[r-1]->end, p, -1);
				++numTrimmed;
				paths[p].len_path--;
				paths[p].begin_length = 0;
			}
		}
	}
	/*	printf("trimmed: %d paths.\n", numTrimmed);*/
}

int count_path(PATH *path, int num_path)
{
	int	i, n, m, q, j, k;
	NODES	*vertex;

	n = 0;
	m = 0;
	q = 0;
	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path > 1)	{
			q += path[i].len_path;
			n ++;
			if(path[i].len_path == 1)	{
				m ++;
			}
		} else if(path[i].len_path < 0)	{
			printf("Path length error!\n");
			exit(-1);
		}
	}

	return(q);
}

