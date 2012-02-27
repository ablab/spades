/***************************************************************************
* Title:          readpath.c
* Author:         Haixu Tang
* Created:        Jun. 2002
* Last modified:  May. 2004
*
* Copyright (c) 2001-2004 The Regents of the University of California
* All Rights Reserved
* See file LICENSE for details.
***************************************************************************/
#include <algorithm>
#include <vector>
#include <assert.h>
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>
void remove_readinterval(EDGE *edge, int index);
int filter_path_read(PATH *path, int num, int num_path, int filterMates);
int filter_path_path(PATH *path, int num, int num_path);
void allocpath(PATH *path, int len);
void findpath(PATH *path, EDGE *edge, int *len_seq, int num_seq, char *label);
void backpath(PATH *path, NODES *vertex, int reads, int end, int *len_seq);
void forpath(PATH *path, NODES *vertex, int reads, int begin, int *len_seq);
int readpath(NODES **vertex, int *num, PATH *path, int *len_seq, int num_seq, int *chim, int *num_chim);
int comppath(PATH *path1, PATH *path2);


int EdgeConnectsPaths(EDGE *edge) {
	if (edge->begin->num_nextedge == 1 or
			edge->end->num_lastedge == 1) 
		return 1;
	else {
		int e;
		/* Make sure the edge is not a short directed cycle.*/
		for (e = 0; e < edge->end->num_nextedge; e++) {
			if (edge->end->nextedge[e]->end == edge->begin) {
				/*				printf("edge: %d creates a short cycle!\n", edge->index);*/
				return 1;
			}
		}
		/* The edge does not create a short directed cycle, get rid of it.*/
		return 0;
	}
}

int RemoveLowCoverageEdges(NODES **vertices, int numVertices,
														PATH *paths,  int numPaths,
														int totalLength, int totalReadStarts,
														float maxStddev, int minCoverage) {
	int v, e;
	double expNumReads, stddevNumReads;
	double readsPerNucleotide = totalReadStarts / (1.0*totalLength);
	double normstddev;
	EDGE *edge;
	int removedAnEdge;
	do {
		removedAnEdge=  0;
		for (v = 0; v < numVertices; v++) {
			for (e = 0; e < vertices[v]->num_nextedge;e++) {
				expNumReads = vertices[v]->nextedge[e]->length * readsPerNucleotide;
				stddevNumReads = sqrt(expNumReads * readsPerNucleotide);
				edge = vertices[v]->nextedge[e];
				normstddev = (expNumReads - edge->multip)/ stddevNumReads;
				if (edge->multip < minCoverage and 
						( normstddev > maxStddev) and
						edge->length < 120 and
						! EdgeConnectsPaths(edge)) {
					/*
					printf("removing edge %d multip: %d len: %d  stddev: %f srcdeg: %d destdeg: %d\n",
								 edge->index, edge->multip, edge->length,
								 normstddev, edge->begin->num_nextedge, edge->end->num_lastedge);
					*/
					removedAnEdge = 1;
					MarkEdgeForRemoval(edge);
				}
			}
		}
		RemovePathsWithRemovedEdges(paths, numPaths);
		EraseEdgesMarkedForDeletion(vertices, numVertices);
	} while (removedAnEdge == 1);
	numVertices = merge_graph_path(vertices, numVertices, paths, numPaths, 0);
	return numVertices;
}

void CheckEdgeMultiplicities(NODES **vertices, int numVertices, int totalLength, int totalReadStarts) {
	double expNumReads, stddevNumReads;
	double readsPerNucleotide = totalReadStarts / (1.0*totalLength);
	int v, e;
	printf("rpn: %f\n", readsPerNucleotide);
	for (v = 0; v < numVertices; v++ ) {
		for (e = 0; e < vertices[v]->num_nextedge; e++ ) {
			expNumReads = vertices[v]->nextedge[e]->length * readsPerNucleotide;
			stddevNumReads = sqrt(expNumReads * (1- readsPerNucleotide));
			if ( expNumReads - vertices[v]->nextedge[e]->multip> 2*stddevNumReads) {
				printf("rpn: %f, edge: %d  multip: %d  length: %d expected: %f deviation: %f\n",
							 readsPerNucleotide,
							 vertices[v]->nextedge[e]->index,
							 vertices[v]->nextedge[e]->multip,
							 vertices[v]->nextedge[e]->length,
							 expNumReads, 
							 fabs(expNumReads -  vertices[v]->nextedge[e]->multip) / stddevNumReads);
			}
		}
	}
}

void CountCoverage(NODES **vertices, int numVertices, int *totalLengthParam, int *totalReadsParam) {
	int totalLength, totalReads;
	
	int v, e, ri;
	EDGE *edge;
	totalReads = 0; totalLength = 0;
	for (v = 0; v < numVertices; v++ ) {
		for (e = 0; e < vertices[v]->num_nextedge; e++ ) {
			edge = vertices[v]->nextedge[e];
			totalLength += edge->length;
			for (ri = 0; ri < edge->multip; ri++) {
				if (edge->readinterval[ri].begin == 0) {
					++totalReads;
				}
			}
		}
	}
	*totalLengthParam = totalLength;
	*totalReadsParam  = totalReads;
}


int readpath(NODES **vertex, int *num, PATH *path, int *len_seq, int num_seq, int *chim, int *num_chim)
{
	int	i, j, k, l, m, n;
	int	num_path;
	int	nch;
	int	num_vertex;
	int	reads;
	char	*label;
	EDGE	*edge, *bal_edge;

	num_vertex = *num;
	label = (char *) ckalloc(2 * num_seq * sizeof(char) + 1);

/*	Define read paths	*/
	
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
		  edge = vertex[i] -> nextedge[j];
			findpath(path, edge, len_seq, num_seq, label);
		}
	}

	for(i = 0; i < num_seq; i ++)	{
		if(path[i].len_path != path[i + num_seq].len_path)	{
		  printf("i %d %d len_path %d %d\n", i, i + num_seq, path[i].len_path, path[i + num_seq].len_path);
			printf("PATH %d ", i);
			for(j = 0; j < path[i].len_path; j ++)	{
				printf("%d(%d) ", path[i].edge[j] -> length, path[i].edge[j] -> multip);
			}
			printf("\n");
			printf("PATH %d ", i + num_seq);
			for(j = 0; j < path[i + num_seq].len_path; j ++)	{
				printf("%d(%d) ", path[i + num_seq].edge[j] -> length, path[i + num_seq].edge[j] -> multip);
			}
			printf("\n");
			printf("len_seq %d %d\n", len_seq[i], len_seq[i + num_seq]);
			for(m = 0; m < path[i].len_path; m ++)	{
				remove_readinterval(path[i].edge[m], i);
			}
			for(m = 0; m < path[i + num_seq].len_path; m ++)	{
				remove_readinterval(path[i + num_seq].edge[m], i + num_seq);
			}
			path[i].len_path = path[i + num_seq].len_path = 0;
/*
			getchar();
*/
		}
	}

	*num = num_vertex;
	free((void *) label);
	return(2 * num_seq);
}

int  LookupReadintervalIndex(EDGE *edge, int index, int pos) {
	int i;
	while (i < edge->multip) {
		if (edge->readinterval[i].eq_read == index and 
				((pos == 0 and edge->readinterval[i].begin == 0) or
				 (pos != 0 and edge->readinterval[i].begin != 0))) {
			return i;
		}
	}
	return -1;
}

int MarkIndexedReadIntervalForRemoval(EDGE *edge, int pos){ 
	/* Make this read interval unrecognizeable.*/
	edge->readinterval[pos].begin = -1;
	edge->readinterval[pos].eq_read = -1;
}

int RemoveMarkedReadIntervals(EDGE *edge) {
	int numRemoved = 0;
	int intv;
	for (intv = 0; intv <edge->multip; intv++) {
		if (edge->readinterval[intv].begin >= 0) 
			edge->readinterval[intv-numRemoved] = edge->readinterval[intv];
		else
			/* this interval has been marked as removed,don't copy it.*/
			++numRemoved;
	}
	edge->multip -= numRemoved;

}

int RemoveIndexedReadinterval(EDGE *edge, int pos) {
	int i;
	i = pos;
	assert(i < edge->multip);
	i++;
	while (i < edge->multip) {
		edge->readinterval[i-1] = edge->readinterval[i];
		i++;
	}
	edge->multip--;
}
	

void remove_readinterval(EDGE *edge, int index)
{
	int	i, j, k, l;

	i = 0;
	while(i < edge -> multip)	{
		if(edge -> readinterval[i].eq_read == index)	{
			edge -> readinterval[i] = edge -> readinterval[edge -> multip - 1];
			edge -> multip --;
		} else	{
			i ++;
		}
	}
}

void findpath(PATH *path, EDGE *edge, int *len_seq, int num_seq, char *label)
{
  /*
    Input:
    path* An array of paths, one path per read.
    edge* a list of all edges in the graph
    len_seq* the list of all sequence lengths in the read list
    num_seq  the number of sequences
    label   Boolean-array that marks whether or not a read has been processed.

    Output:
    Store the list of edges that each read traverses in 'path'.

  */

	int i, j, k, l, reads, m;
	int	i1, i2;
	READINTERVAL *readinterval;
	PATH path1, path2;
	
	for(i = 0; i < edge -> multip; i ++)	{
		readinterval = &(edge -> readinterval[i]);
		reads = readinterval -> eq_read;
		/*
		  If this read has already been processed (an edge that
		  includes this read has been traversed), continue.
		*/
		if(label[reads] == 1)	{
			continue;
		}
		
		label[reads] = 1;
		path1.len_path = path2.len_path = 0;

		/* The length of the path is at most the length of a read.  The actual length
		   of the path will be stored later.
		 */
		allocpath(&path1, len_seq[reads]);
		allocpath(&path2, len_seq[reads]);

		path1.begin_length = readinterval -> offset;
		path2.end_length   = edge -> length - (readinterval -> offset + readinterval -> length);

		if(readinterval -> begin > 0)	{
			backpath(&path1, edge -> begin, reads, readinterval -> begin, len_seq);
		}
		forpath(&path2, edge -> end, reads, readinterval -> begin + readinterval -> length - VERTEX_SIZE, len_seq);
		path[reads].len_path = path1.len_path + path2.len_path + 1;
		allocpath(&path[reads], path[reads].len_path);

		k = 0;
		for(m = 0; m < path1.len_path; m ++)	{
			path[reads].edge[k ++] = path1.edge[path1.len_path - 1 - m];
		}
		path[reads].edge[k ++] = edge;
		for(m = 0; m < path2.len_path; m ++)	{
			path[reads].edge[k ++] = path2.edge[m];
		}
		path[reads].len_path = k;
		path[reads].begin_length = path1.begin_length;
		path[reads].end_length = path2.end_length;
		/*		path[reads].index = reads;*/
		path[reads].seqLength = 0;
		if (path[reads].len_path > 0) {
		  if (path[reads].len_path == 1) {
		     path[reads].seqLength = path[reads].edge[0]->length - (path[reads].end_length + path[reads].begin_length);
 		  }
	  	  else {
		     path[reads].seqLength += path[reads].edge[0]->length - path[reads].begin_length;;
		     path[reads].seqLength += path[reads].edge[path[reads].len_path-1]->length - path[reads].end_length;
		     path[reads].seqLength -= VERTEX_SIZE;
  	             for (m = 1; m  < path[reads].len_path - 1; m++ ) {
	                path[reads].seqLength += path[reads].edge[m]->length - VERTEX_SIZE;
 		     }
                  }
		  /*
 		     printf("path: %d %d edges %d seqlen bl: %d ebl: %d el: %d eel %d\n", 
			 reads, path[reads].len_path, path[reads].seqLength,
			 path[reads].begin_length, path[reads].edge[0]->length,
			 path[reads].end_length, path[reads].edge[path[reads].len_path-1]->length);
                  */
 		}
	
		if(path[reads].len_path == 0)	free((void **) path[reads].edge);
		free((void **) path1.edge);
		free((void **) path2.edge);
	}
}

void backpath(PATH *path, NODES *vertex, int reads, int end, int *len_seq)
{
	int	i, j, k, l;
	READINTERVAL *readinterval;

	for(i = 0; i < vertex -> num_lastedge; i ++)	{
		for(j = 0; j < vertex -> lastedge[i] -> multip; j ++)	{
			readinterval = &(vertex -> lastedge[i] -> readinterval[j]);
			if(readinterval -> eq_read == reads &&
			   readinterval -> begin + readinterval -> length - VERTEX_SIZE - end == 0)	{
				path -> edge[path -> len_path ++] = vertex -> lastedge[i];
				path -> begin_length = readinterval -> offset;
		//		path -> begin_length = vertex -> lastedge[i] -> length - readinterval -> offset + 1;
				if(readinterval -> begin > 0)	{
					backpath(path, vertex -> lastedge[i] -> begin, readinterval -> eq_read,
					 	 readinterval -> begin, len_seq);
				}
				return;
			}
		}
	}
}

void forpath(PATH *path, NODES *vertex, int reads, int begin, int *len_seq)
{
	int	i, j, k, l;
	int	cut;
	READINTERVAL *readinterval;

	for(i = 0; i < vertex -> num_nextedge; i ++)	{
		for(j = 0; j < vertex -> nextedge[i] -> multip; j ++)	{
			readinterval = &(vertex -> nextedge[i] -> readinterval[j]);
			if(readinterval -> eq_read == reads && 
			   readinterval -> begin - begin == 0)	{
				path -> edge[path -> len_path ++] = vertex -> nextedge[i];
				path -> end_length   = vertex -> nextedge[i] -> length - (readinterval -> offset + readinterval -> length);
//				path -> end_length = readinterval -> offset + readinterval -> length + 1;
				forpath(path, vertex -> nextedge[i] -> end, readinterval -> eq_read,
				 	 readinterval -> begin + readinterval -> length - VERTEX_SIZE, len_seq);
				return;
			}
		}
	}
}

void	allocpath(PATH *path, int len)
{
	path -> edge = (EDGE **) ckalloc(len * sizeof(EDGE *));
}

int EqualPaths(PATH *path1, PATH *path2) {

  if (path1->len_path != path2->len_path) 
    return 0;

  int i;
  for (i = 0; i < path1->len_path; i++ ) {
    if (path1->edge[i] != path2->edge[i]) 
      return 0;
  }
  return 1;
}


int comppath(PATH *path1, PATH *path2)
{
	int	i, j, k, l, m;

	l = min(path1 -> len_path, path2 -> len_path);
	if(path1 -> len_path == path2 -> len_path)	{
		for(k = 0; k < l; k ++)	{
			if(path1 -> edge[k] != path2 -> edge[k])	{
				break;
			}
		}
		if(k == l)	{
			return(0);
		}
	} else if(path1 -> len_path > path2 -> len_path)	{
		for(i = 0; i <= path1 -> len_path - path2 -> len_path; i ++)	{
			for(k = 0; k < l; k ++)	{
				if(path1 -> edge[i + k] != path2 -> edge[k])	{
					break;
				}
			}
			if(k == l)	{
				return(2);
			}
		}
	} else	{
		for(i = 0; i <= path2 -> len_path - path1 -> len_path; i ++)	{
			for(k = 0; k < l; k ++)	{
				if(path2 -> edge[i + k] != path1 -> edge[k])	{
					break;
				}
			}
			if(k == l)	{
				return(1);
			}
		}
	}
	return(-1);
}


int ComparePaths(const PATH& path1, const PATH& path2) {
  if (path1.len_path == 0) {
    return  0;
  }
  if (path2.len_path == 0) {
    return 0;
  }
  int minLength;
  minLength = path1.len_path < path2.len_path ? path1.len_path : path2.len_path;
  int i;
  for (i = 0; i < minLength; i++ ){
    if (path1.edge[i] != path2.edge[i]) {
      return path1.edge[i] < path2.edge[i];
    }
  }
  // all else equal, the lesser path is the shorter one
  return path1.len_path < path2.len_path;
}

class ComparePathToList { 
 public:
  PATH *paths;
  int operator()(int &p1, const PATH &p2) {
    return ComparePaths(paths[p1], p2);
  }
};

class SortPaths {
 public:
  PATH *paths;
  int operator()(int p1, int p2) {
    return ComparePaths(paths[p1], paths[p2]);
  }
};

int filter_path_read(PATH *path, int num_unpaired_path, int num_path, int filterMates)
{
	int	i, j, k, m, n, l;
	int	*label;

/*	Reorder paths by removing single read edges paths with length 0 (reads in a single edge)	*/

	n = 0;
	int nseq = num_unpaired_path/2;
	int numReadPath = 0;
	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path > 1)	{
			if (i < nseq) {
				assert(path[i].len_path == path[i+nseq].len_path);
			}
			path[n] = path[i];
			if(i < num_unpaired_path)	{
				numReadPath++;
				path[n].readindex = i + 1;
			}
			n ++;
		} else if(path[i].len_path > 0)	{
			free((void **) path[i].edge);
		}
	}
	num_path = n;

/*	Filter low coverage paths	*/

	label = (int *) ckalloc(num_path * sizeof(int));
	int *sortedPaths;
	sortedPaths = (int*) ckalloc(num_path * sizeof(int));
	
	for(i = 0; i < num_path; i ++)	{
		label[i] = 1;
		sortedPaths[i] = i;
	}
	printf("sorting paths\n");
	SortPaths pathSort;
	pathSort.paths = path;
	std::sort(sortedPaths, sortedPaths + num_path, pathSort);
	printf("done\n");
	int curSortedPath;
	int subPath1, subPath2;
	PATH tempPath;
	int* pathIndexIt;
	int pathIndex;
	ComparePathToList pathToList;
	pathToList.paths = path;
	int sortedPathIndex;
	int pathMultiplicity = 1;
	int nextSortedPath;
	curSortedPath = 0;
	while (curSortedPath < num_path) {
	  pathIndex = sortedPaths[curSortedPath];
	  // add a little blurb here to advance while the paths are the same
	  nextSortedPath = curSortedPath + 1;
	  while (nextSortedPath < num_path &&
		 EqualPaths(&path[sortedPaths[curSortedPath]],
								&path[sortedPaths[nextSortedPath]])) {
	    nextSortedPath++;
	  }
	  pathMultiplicity = nextSortedPath - curSortedPath;
		 
	  // now search all subpaths of this path, and increment by this path's multiplicity
	  int sortedSearchPath;
	  for (subPath1 = 0; subPath1 < path[sortedPaths[curSortedPath]].len_path - 1; subPath1++ ) {
	    // only search paths of length 2
	    for (subPath2 = subPath1 + 1; subPath2 < path[sortedPaths[curSortedPath]].len_path; subPath2++ ) {
	      tempPath.edge = &(path[sortedPaths[curSortedPath]].edge[subPath1]);
	      tempPath.len_path = subPath2 - subPath1 + 1;
	      
	      // search for this path in the list of sorted paths
	      pathIndexIt = std::lower_bound(sortedPaths, sortedPaths + num_path, tempPath, pathToList);
	      // compute the start of the path index
	      sortedSearchPath = pathIndexIt - sortedPaths;

	      // increment all paths that are the same as this substring (from subpath1 ... subpath2)
	      while (sortedSearchPath < num_path && EqualPaths(&path[sortedPaths[sortedSearchPath]], &tempPath)) {
					label[sortedPaths[sortedSearchPath]]+= pathMultiplicity;
					sortedSearchPath++;
	      }
	    }
	  }
	  // advance to the next path, skipping duplicates
	  curSortedPath += pathMultiplicity;
	}

	k = 0;
	int numReadPathsDeleted = 0;

	// This loop removes low coverage paths.
	// If the path 'i' is a read path, path[i].readindex > 0
	// 
	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path >= 2 && label[i] <= LOW_COV_PATH && path[i].seqLength < 200)	{
			if (!filterMates  || path[i].readindex == 0) {
								printf("removing path %d, coverage: %d\n", i, label[i]);
				for(m = 0; m < path[i].len_path; m ++)	{
					remove_readinterval(path[i].edge[m], path[i].readindex - 1);
				}
				free((void **) path[i].edge);
				path[i].len_path = 0;
				
				
				if (!filterMates)
					numReadPathsDeleted++;
				k ++;
			}
		}
	}
	free((void *) label);

/*	Again, order paths by removing single read edges paths with length 1 (reads in a single edge)	*/

	n = 0;
	int numReads = numReadPath/2;
	nseq = num_path/2;
	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path > 1)	{
			if (i < numReads) 
				assert(path[i].len_path == path[i+numReads].len_path);
			path[n] = path[i];
			// The read index was set to read+1 so that the 0 could be3 
			// used as a sentinal for deletion.  Fix it back to the
			// correct read index here.
			// path[n].readindex--;
			n ++;
		} else if(path[i].len_path > 0)	{
			free((void **) path[i].edge);
		}
	}
	num_path = n;
	return(num_path);
}

int filter_path_mate(PATH *path, int num, int num_path)
{
	int	i, j, k, m, n, l;
	int	*label;

/*	Reorder paths by removing single read edges paths with length 0 (reads in a single edge)	*/

	n = 0;
	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path > 1)	{
			path[n] = path[i];
			if(i < num)	{
				path[n].readindex = i + 1;
			}
			n ++;
		} else if(path[i].len_path > 0)	{
			free((void **) path[i].edge);
		}
	}
	num_path = n;

/*	Filter low coverage paths	*/

	label = (int *) ckalloc(num_path * sizeof(int));
	for(i = 0; i < num_path; i ++)	{
		label[i] = 1;
	}
	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path == 0)	continue;
		for(j = i + 1; j < num_path; j ++)	{
			if(path[j].len_path == 0)	continue;
			k = comppath(&path[i], &path[j]);
			if(k == 0)	{
				label[i] ++;
				label[j] ++;
			} else if(k == 1)	{
				label[i] ++;
			} else if(k == 2)	{
				label[j] ++;
			}
		}
	}

	k = 0;
	for(i = 0; i < num_path; i ++)	{
/*	Only mate-paths should be filtered here	(therefore path[i].readinex == 0 */
		if(path[i].len_path >= 2 && label[i] <= LOW_COV_PATH)	{
			if(path[i].readindex == 0)	{
				for(m = 0; m < path[i].len_path; m ++)	{
					remove_readinterval(path[i].edge[m], path[i].readindex - 1);
				}
				free((void **) path[i].edge);
				path[i].len_path = 0;
				k ++;
			}
		}
	}
	free((void *) label);

/*	Again, order paths by removing single read edges paths with length 0 (reads in a single edge)	*/

	n = 0;
	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path > 1)	{
			path[n] = path[i];
			n ++;
		} else if(path[i].len_path > 0)	{
			free((void **) path[i].edge);
		}
	}
	num_path = n;
	return(num_path);
}


void free_path(int num_path, PATH *path)
{
  int i;

  for (i = 0; i < num_path; i ++) {
    if(path[i].edge && path[i].edge > 0)	{
    	free((void **) path[i].edge);
    }
  }
  free((void *) path);
}
