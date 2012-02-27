/***************************************************************************
 * Title:          countedge.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <assert.h>
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

typedef struct TPathAddress {
	int path_index;
	int path_pos;
} PathAddress;

extern int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;

int count_edge(NODES **vertex, int num_vertex, int **num_pa);
void count_multip(NODES **vertex, int num_vertex);
int count_edge_simp(NODES **vertex, int num_vertex, int **num_pa);
int forwardtanglelen(EDGE *edge, int *ave, int *multip);
int backtanglelen(EDGE *edge, int *ave, int *multip);
int count_tangle(NODES **vertex, int num_vertex, int disttangle[][7]);
int count_bal(NODES *vertex);
void getmaxedge(EDGE *edge, EDGE **maxedge);
int LookupPathIntervalIndicesForward(PATH *paths, 
																		 int pathIndex, 
																		 int *intervalIndices, 
																		 int *lenTracedPath);

int LookupPathIntervalIndicesReverse(PATH *paths, int pathIndex, 
																		 EDGE *startEdge, int startInterval,
																		 int *intervalIndices, int *lenTracedPath);

int LookupPathIntervalIndices(PATH *paths, int pathIndex, int *intervalIndices);

void AppendFirstReadIntervals(READINTERVAL **intvList, int &intvListLen,
															READINTERVAL **appList, int appListLen, int offset);
void GrowFirstReadInterval(READINTERVAL **intvList,
													 int *intvListLen,
													 int numPath);

void GrowReadInterval(READINTERVAL *srcInterval,
											READINTERVAL *interval);

int StorePathsThroughVertex(PATH *paths,
														EDGE *enterEdge, NODES *vertex, EDGE *exitEdge, 
														PathAddress **pathAddresses);

int CountPathsThroughVertex(PATH *paths, EDGE *enterEdge, NODES *vertex, EDGE *exitEdge,
														PathAddress *pathAddresses);

int StorePathsInEdge(PATH *paths, NODES *src, EDGE *edge, PathAddress **pathIndices);
int CountPathsInEdge(PATH *paths, NODES *src, EDGE *edge, PathAddress *pathIndices);

void AppendReadIntervalsToEdge(EDGE *srcEdge, int *srcEdgeIntervalIndices, int numIntv,
															 EDGE *destEdge, int offset);

void RemovePathsBetweenEdges(PATH *paths, PathAddress *pathAddrs, int numPathIndices, 
														 EDGE *enterEdge, EDGE *exitEdge,
														 int **enterEdgeIntervalIndices);

void AppendPathAddresses(NODES *vertex, PathAddress *pathAddrs, int numPathAddrs, int posOffset);
void AppendVertexPaths(NODES *enterEnd, NODES *exitBegin);

int RemoveDuplicatePathIndices(PathAddress **pathIndices, int numIndices);
void RepointPathEdges(PATH *paths,
											EDGE *newEdge, int pos, 
											PathAddress *pathAddrs, int numPathAddrs, int newOffset );

void StoreEdgePathIntervalIndices(PATH *paths, PathAddress *pathAddrs, int numPathAddrs, 
																	int **edgeIntervalInidices );

int StraightenShortDirectedCycle(NODES ***vertices, int numVertices,
																 PATH *paths,
																 EDGE *enterEdge, EDGE *repeat, EDGE *link, EDGE *exitEdge,
																 EDGE **repeatCopyP) {

	/*	printf("starting ssdc with %d vertices, lengths: enter: %d (%d)  repeat: %d (%d) link: %d (%d) end: %d(%d)\n", 
				 numVertices, 
				 enterEdge->length, enterEdge->multip, 
				 repeat->length, repeat->multip,
				 link->length, link->multip, 
				 exitEdge->length, exitEdge->multip);
	
	printf("Cycle vertices: %d (%d) %d (%d)\n",
				 enterEdge->end, enterEdge->end->index,
				 exitEdge->begin, exitEdge->begin->index);

	printf("removing directed cycle: %d(%d)%d(%d) %d(%d) %d(%d)\n",
				 enterEdge->index, enterEdge->length,
				 repeat->index, repeat->length,
				 link->index, link->length,
				 exitEdge->index, exitEdge->length);
	*/

	/*
		Currently there is a cycle A->B->C->B->D (where B is the same
		vertex).  We need to replace this with a nonbranchign path
		A->B->C->B'->D.  This requires one new edge, B', and two new
		vertices, dest(B), and dest(C).
		
		I don't think there shoudl be a path that goes all the way through
		C, because the path should have been transformed.  If it does,
		we'll just have to deal with it.
		
		Here we call 'A' the enter edge and 'D' the exit edge (since they
		enter and exit the cycle).
	*/

	/*
		Allocate the new structures.  Will need a new edge, and two
		new vertices.
	*/

	EDGE *repeatCopy;
	NODES *repeatSrc, *repeatDest, *linkDest, *repeatCopyDest;
	*repeatCopyP = NULL;
	repeatDest = repeat->end;
	repeatSrc  = repeat->begin;

	int enterLength, repeatLength, linkLength;
	enterLength  = enterEdge->length;
	repeatLength = repeat->length;
	linkLength   = link->length;
	/*
		Ok, fixing the toplogy is done.
		Now the paths. 
	*/

	/*
		Collect paths that go through the whirl.
		 
		Now the paths may be added back to the graph appropriately:
		P1 may pass through the entire cycle.
		P2 may go from repeat->link->repeatCopy->exit.
		P3 may go from linker->repeat2->exit.
	*/

	PathAddress *enterRepeatAddrs,*repeatLinkAddrs, *linkRepeatAddrs, *repeatExitAddrs,
		*repeatAddrs, *linkAddrs;
	enterRepeatAddrs = repeatLinkAddrs = linkRepeatAddrs = repeatExitAddrs = NULL;

	int numEnterRepeat, numRepeatLink, numLinkRepeat, numRepeatExit;
	int numRepeat, numLink;

	int *enterRepeatReadIntervals, *repeatLinkReadIntervals, 
		*linkRepeatReadIntervals, *repeatExitReadIntervals;
	
	enterRepeatReadIntervals = repeatLinkReadIntervals =
		linkRepeatReadIntervals = repeatExitReadIntervals = NULL;
	

	/* Look for paths that go through the repeat more than once,
		 or that pass through the repeat without using the loop.*/
	int vX, pX;
	for (vX = 0; vX < numVertices; vX++) {
		for (pX = 0; pX < (*vertices)[vX]->num_path; pX++) {
			assert((*vertices)[vX]->path_pos[pX] <= paths[(*vertices)[vX]->path_index[pX]].len_path);
		}
	}

	int repPath, repPathIndex;
	int repCount, linkCount, enterCount, exitCount;
	int repPathPos;
	int *multiLoopPaths = (int*) ckalloc(sizeof(int) * repeatSrc->num_path);
	int numMultiLoops = 0;
	int enterIndex, exitIndex;
	/*
		printf("repeat source has %d paths crossing it\n", repeatSrc->num_path);
	*/
	for (repPath = 0; repPath < repeatSrc->num_path; repPath++ ) {
		repCount = 0;
		linkCount = 0;
		enterCount = 0; exitCount = 0;
		enterIndex = -1; exitIndex = -1;
		repPathIndex = repeatSrc->path_index[repPath];
		for (repPathPos = 0; repPathPos < paths[repPathIndex].len_path; repPathPos++) {
			if (paths[repPathIndex].edge[repPathPos] == repeat)
				++repCount;
			else if (paths[repPathIndex].edge[repPathPos] == link)
				++linkCount;
			else if (paths[repPathIndex].edge[repPathPos] == enterEdge) {
				++enterCount; enterIndex = repPathPos;
			}
			else if (paths[repPathIndex].edge[repPathPos] == exitEdge) {
				++exitCount; exitIndex = repPathPos;
			}
		}	
		/*		printf("path: %d %d %d %d %d\n", repPathIndex, enterCount, 
					repCount, linkCount, exitCount);*/
		if (repCount > 2 || linkCount > 1) {
			/* 			printf("path: %d is a multi-loop path\n", repPathIndex);*/
			RemovePath(paths, repPathIndex);
			++numMultiLoops;
			repPath--;
		}
		else if (enterCount >= 1 and exitCount >= 1 and repCount >= 1 and linkCount == 0) {
			/*			printf("path: %d goes straight through the loopg, truncating it\n",
						 repPathIndex);
			*/
			RemovePath(paths, repPathIndex);
			repPath--;
			/*
			for (repPathPos = 0; repPathPos < paths[repPathIndex].len_path; repPathPos++) {
				if (paths[repPathIndex].edge[repPathPos] == exitEdge)
					break;
			}
			int truncPathLength = repPathPos;
			for (; repPathPos < paths[repPathIndex].len_path; repPathPos++) {
				RemovePathFromVertex(paths[repPathIndex].edge[repPathPos]->end,
														 repPathIndex);
				remove_readinterval(paths[repPathIndex].edge[repPathPos],
														paths[repPathIndex].readindex-1);
			}
			paths[repPathIndex].len_path = truncPathLength;
			paths[repPathIndex].end_length = 0;
			*/
		}
		else if (enterCount > 1 or exitCount > 1) {
			/*
				printf("path: %d has complicated loops, removing it\n", repPathIndex);
			*/
			RemovePath(paths, repPathIndex);
			repPath--;
		}
		else if (enterIndex != -1 and exitIndex != -1 and enterIndex > exitIndex) {
			/*
				printf("path: %d has return loops, removing it\n", repPathIndex);
			*/
			RemovePath(paths, repPathIndex);
			repPath--;
		}
	}


	for (vX = 0; vX < numVertices; vX++) {
		for (pX = 0; pX < (*vertices)[vX]->num_path; pX++) {
			assert((*vertices)[vX]->path_pos[pX] <= paths[(*vertices)[vX]->path_index[pX]].len_path);
		}
	}


	/*	printf("removed: %d multi loops\n", numMultiLoops);*/
	free(multiLoopPaths);
	
	/* Lookup paths enter->repeat*/
	numEnterRepeat = 
		StorePathsThroughVertex(paths, enterEdge, repeatSrc, repeat, &enterRepeatAddrs);
	numEnterRepeat = 
		RemoveDuplicatePathIndices( &enterRepeatAddrs, numEnterRepeat);

	/* Remove them from edges and vertices.*/

	RemovePathsBetweenEdges(paths, enterRepeatAddrs, numEnterRepeat,
													enterEdge, exitEdge, &enterRepeatReadIntervals);


	AppendPathAddresses(enterEdge->begin, enterRepeatAddrs, numEnterRepeat, 0);
	AppendPathAddresses(enterEdge->end, enterRepeatAddrs, numEnterRepeat, 1);

	/* Find paths that are in only the repeat.*/
	numRepeat = 
		StorePathsInEdge(paths, repeatSrc, repeat, &repeatAddrs);

	int *repeatIntervalIndices;
	StoreEdgePathIntervalIndices(paths, repeatAddrs, numRepeat, &repeatIntervalIndices);
	AppendReadIntervalsToEdge(repeat, repeatIntervalIndices, numRepeat,
														enterEdge, enterLength);
	
	RepointPathEdges(paths, enterEdge, 0, repeatAddrs, numRepeat, 
									 enterLength - VERTEX_SIZE);
	
	for (vX = 0; vX < numVertices; vX++) {
		for (pX = 0; pX < (*vertices)[vX]->num_path; pX++) {
			assert((*vertices)[vX]->path_pos[pX] <= paths[(*vertices)[vX]->path_index[pX]].len_path);
		}
	}

	
	/* Lookup paths that start in the repeat and continue to the link.*/
	numRepeatLink = 
		StorePathsThroughVertex(paths, repeat, repeatDest, link, &repeatLinkAddrs);

	numRepeatLink = 
		RemoveDuplicatePathIndices(&repeatLinkAddrs, numRepeatLink);


	RemovePathsBetweenEdges(paths, repeatLinkAddrs, numRepeatLink,
													repeat, exitEdge, &repeatLinkReadIntervals);

	AppendReadIntervalsToEdge(repeat, repeatLinkReadIntervals, numRepeatLink, 
														enterEdge, enterEdge->length);


	enterEdge->length += repeat->length - VERTEX_SIZE;
	for (vX = 0; vX < numVertices; vX++) {
		for (pX = 0; pX < (*vertices)[vX]->num_path; pX++) {
			assert((*vertices)[vX]->path_pos[pX] <= paths[(*vertices)[vX]->path_index[pX]].len_path);
		}
	}


	/* Move edges that are only in the link edge.*/
	numLink =
		StorePathsInEdge(paths, repeatDest, link, &linkAddrs);

	int *linkIntervalIndices;
	StoreEdgePathIntervalIndices(paths, linkAddrs, numLink, &linkIntervalIndices);
	AppendReadIntervalsToEdge(link, linkIntervalIndices, numLink,
														enterEdge, enterLength + repeatLength - VERTEX_SIZE);
	
	RepointPathEdges(paths, enterEdge, 0, linkAddrs, numLink, 
									 enterLength + repeatLength + - (2*VERTEX_SIZE));

	for (vX = 0; vX < numVertices; vX++) {
		for (pX = 0; pX < (*vertices)[vX]->num_path; pX++) {
			assert((*vertices)[vX]->path_pos[pX] <= paths[(*vertices)[vX]->path_index[pX]].len_path);
		}
	}

	/* Lookup paths that start in the link, and go to the repeat. */


	numLink = StorePathsInEdge(paths, repeatDest, link, &linkAddrs);


	numLinkRepeat = 
		StorePathsThroughVertex(paths, link, repeatSrc, repeat, &linkRepeatAddrs);

	numLinkRepeat = 
		RemoveDuplicatePathIndices( &linkRepeatAddrs, numLinkRepeat);

	/* Remove them from edges and vertices.*/

	RemovePathsBetweenEdges(paths, linkRepeatAddrs, numLinkRepeat,
													link, exitEdge, &linkRepeatReadIntervals);
	
	/* Lookup paths that start in the repeat, and continue to the exit.*/
	for (vX = 0; vX < numVertices; vX++) {
		for (pX = 0; pX < (*vertices)[vX]->num_path; pX++) {
			assert((*vertices)[vX]->path_pos[pX] <= paths[(*vertices)[vX]->path_index[pX]].len_path);
		}
	}

	AppendReadIntervalsToEdge(link, linkRepeatReadIntervals, numLinkRepeat, 
														enterEdge, enterEdge->length);

	enterEdge->length += link->length - VERTEX_SIZE;

	/* Lookup paths that start in the repeat, and continue to the exit.*/
	for (vX = 0; vX < numVertices; vX++) {
		for (pX = 0; pX < (*vertices)[vX]->num_path; pX++) {
			assert((*vertices)[vX]->path_pos[pX] <= paths[(*vertices)[vX]->path_index[pX]].len_path);
		}
	}

	numRepeatExit = 
		StorePathsThroughVertex(paths, repeat, repeatDest, exitEdge, &repeatExitAddrs);
	numRepeatExit = 
		RemoveDuplicatePathIndices( &repeatExitAddrs, numRepeatExit);

	/* Remove them from edges and vertices.*/

	RemovePathsBetweenEdges(paths, repeatExitAddrs, numRepeatExit,
													repeat, exitEdge, &repeatExitReadIntervals);

	AppendReadIntervalsToEdge(repeat, repeatExitReadIntervals, numRepeatExit, 
														enterEdge, enterEdge->length);

	for (vX = 0; vX < numVertices; vX++) {
		for (pX = 0; pX < (*vertices)[vX]->num_path; pX++) {
			assert((*vertices)[vX]->path_pos[pX] <= paths[(*vertices)[vX]->path_index[pX]].len_path);
		}
	}

	AppendPathAddresses(enterEdge->begin, repeatLinkAddrs, numRepeatLink, 0);
	AppendPathAddresses(enterEdge->end, repeatLinkAddrs, numRepeatLink, 1);

	RepointPathEdges(paths, enterEdge, 0, repeatLinkAddrs, numRepeatLink, 
									 enterLength - VERTEX_SIZE);

	AppendPathAddresses(enterEdge->begin, linkRepeatAddrs, numLinkRepeat, 0);
	AppendPathAddresses(enterEdge->end, linkRepeatAddrs, numLinkRepeat, 1);

 
	RepointPathEdges(paths, enterEdge, 0, linkRepeatAddrs, numLinkRepeat,
									 enterLength + repeatLength - 2*VERTEX_SIZE);
	

	AppendPathAddresses(enterEdge->begin, repeatExitAddrs, numRepeatExit, 0);
	AppendPathAddresses(enterEdge->end, repeatExitAddrs, numRepeatExit, 1);
											
	RepointPathEdges(paths, enterEdge, 0, repeatExitAddrs, numRepeatExit,
									 enterLength + repeatLength + linkLength - 3*VERTEX_SIZE);


	AppendVertexPaths(enterEdge->end, exitEdge->begin);

	enterEdge->length  += repeat->length - VERTEX_SIZE;
	

	/* Fix the connectivity of the graph. */
	/* Remove the loop from the graph.*/
	enterEdge->end->num_lastedge = 1;
	enterEdge->end->lastedge[0]  = enterEdge;
	enterEdge->end->num_nextedge = 1;
	enterEdge->end->nextedge[0]  = exitEdge;

	/* Fix the  exit edge to skip the loop.*/
	exitEdge->begin = enterEdge->end;

	/* Fix the paths at the vertex between enter and exit.*/


	/* Remove the repeatDest from the graph.*/
	int v;
	int numRemoved = 0;
	for (v = 0; v < numVertices; v++) {
		if ((*vertices)[v] == repeatDest) 
			numRemoved++;
		else 
			(*vertices)[v - numRemoved] = (*vertices)[v];
	}
	numVertices--;
	for (vX = 0; vX < numVertices; vX++) {
		for (pX = 0; pX < (*vertices)[vX]->num_path; pX++) {
			assert((*vertices)[vX]->path_pos[pX] <= 
						 paths[(*vertices)[vX]->path_index[pX]].len_path);
		}
	}


	free(repeat);
	free(link);
	/*	*vertices = (NODES**) realloc(*vertices, sizeof(NODES*) * numVertices);*/
	
	printf("ending with: %d vertices  enterlen: %d exit: %d\n", numVertices, 
				 enterEdge->length, exitEdge->length);
	printf("enter has:%d multp, exit: %d\n", enterEdge->multip, exitEdge->multip);
	printf("------------------------------------------------------------\n");
	return numVertices;
}


void GrowFirstReadInterval(READINTERVAL **intvList,
													 int *intvListLen,
													 int numPath) {
	int path, intv;
	for (path = 0; path < numPath; path++ ) {
		for (intv = 1; intv < intvListLen[path]; intv++) {
			intvList[path][0].length += intvList[path][intv].length - VERTEX_SIZE;
		}
	}
}

int CompareIndices(const void *a, const void *b) {
	return (((PathAddress*)a)->path_index < ((PathAddress*) b)->path_index);
}


int FindIndexBoundByEdge(PATH *paths, int pathIndex, int pathPos,
												 EDGE *exitEdge) {

	int exitEdgeIndex  = -1;

	/* Find what part of the path passes through 'exit'.*/
	for (pathPos++; pathPos < paths[pathIndex].len_path; pathPos++ ) {
		if (paths[pathIndex].edge[pathPos] == exitEdge) {
			return pathPos;
			break;
		}
	}
	return paths[pathIndex].len_path;
}

int FindReadPos(PATH *paths, int pathIndex, int pathPos) {
	int e;
	int readPos = 0;
	int intv;
	int numIntv, eq_read;
	READINTERVAL *intvs;
	for (e = 0; e < pathPos; e++ ) {
		intvs   = paths[pathIndex].edge[e]->readinterval;
		numIntv = paths[pathIndex].edge[e]->multip;
		eq_read = paths[pathIndex].readindex - 1;
		sortreadinterval_index(intvs, numIntv);
		for (intv = 0; intv < numIntv; intv++) {
			if (intvs[intv].eq_read == eq_read and
					intvs[intv].begin == readPos) {
				readPos += intvs[intv].length - VERTEX_SIZE;
				break;
			}
		}
		/* put the read intervals back as they were.*/
		sortreadinterval(intvs, numIntv);
	}
	return readPos;
}

int LookupPathIntervalIndices(PATH *paths, int pathIndex, int *intervalIndices) {

	int p;

	if (paths[pathIndex].len_path == 0) return 1;

	int *forIndices, *revIndices, numForIndices, numRevIndices;
	forIndices = (int*) ckalloc(sizeof(int) * paths[pathIndex].len_path);
	revIndices = (int*) ckalloc(sizeof(int) * paths[pathIndex].len_path);
	
	LookupPathIntervalIndicesForward(paths, pathIndex, forIndices, &numForIndices);
	if (numForIndices == 0) return 0;

	LookupPathIntervalIndicesReverse(paths, pathIndex, paths[pathIndex].edge[0],
																	 forIndices[0], revIndices, &numRevIndices);

	if (numRevIndices + numForIndices != paths[pathIndex].len_path) {
		printf("ERROR, could not trace full path %d %d %d\n", 
					 pathIndex, paths[pathIndex].len_path, numRevIndices + numForIndices);
		return 0;
	}
	int i;
	for (i = numRevIndices-1; i >= 0; i-- ) 
		intervalIndices[numRevIndices - i - 1] = revIndices[i];
	for (i = 0; i < numForIndices; i++) 
		intervalIndices[i+numRevIndices] = forIndices[i];



	free(forIndices);
	free(revIndices);
	return 1;
}


int LookupPathIntervalIndicesForward(PATH *paths, int pathIndex, 
																		 int *intervalIndices, int *lenTracedPath) {

	/* Store the indices of the read intervals startign at startEdge, 
		 moving forwards.  Just go through the first interval found.
	*/
	int intv, eq_read;
	NODES *dest;
	EDGE *curEdge, *nextEdge, *startEdge;
	eq_read = paths[pathIndex].readindex-1;

	curEdge = paths[pathIndex].edge[0];
	
	/* 
		 Locate an instance of this path in curEdge. It may not be
		 correspond to paths[pathIndex].edges[0], since
		 some paths[pathIndex].edges[i] i > 0 may be in the list.
	*/
	for (intv = 0; intv < curEdge->multip; intv++) {
		if (curEdge->readinterval[intv].eq_read == eq_read)
			break;
	}
	int curIntv = 0;
	int readPos, intvLength;
	int nextIntvIndex;
	if (intv < curEdge->multip) {
		intervalIndices[curIntv++] = intv;
		readPos = curEdge->readinterval[intv].begin;
		intvLength = curEdge->readinterval[intv].length;

		do {
			int e;
			dest = curEdge->end;
			/* advance past cur edge.*/
			readPos += intvLength - VERTEX_SIZE;
			nextIntvIndex = -1;
			for (e = 0; e < dest->num_nextedge; e++ ) {
				nextEdge = dest->nextedge[e];
				/* look for the continuation of this path in nextEdge.*/
				nextIntvIndex = FindIntervalBeginningAtPos(nextEdge->readinterval, nextEdge->multip, eq_read, readPos);
				if (nextIntvIndex >= 0) {
					curEdge = nextEdge;
					intervalIndices[curIntv++] = nextIntvIndex;
					intvLength = nextEdge->readinterval[nextIntvIndex].length;
					break;
				}
			}
		} while (nextIntvIndex >= 0 and dest->num_nextedge > 0);
	}
	return *lenTracedPath = curIntv;
}

int LookupPathIntervalIndicesReverse(PATH *paths, int pathIndex, 
																		 EDGE *startEdge, int startInterval,
																		 int *intervalIndices, int *lenTracedPath) {

	/* 
		 Store the indices of the readintervals for the path.  This looks 
		 in reverse, given a specific starting interval.  It's assumed
		 that when looking forward the first interval was stored, so we don't 
		 store that here.
	*/
	int intv, eq_read;
	NODES *src;
	EDGE  *curEdge, *prevEdge;
	eq_read = paths[pathIndex].readindex-1;

	int curIntv = 0;
	int readPos, intvLength;
	int prevIntvIndex;

	curEdge = startEdge;
	readPos = startEdge->readinterval[startInterval].begin;
	intvLength = startEdge->readinterval[startInterval].length;
	do {
		int e;
		src  = curEdge->begin;

		/* 
			 Look for an in-edge that has an interval that ends at 'readPos'.
		*/
		prevIntvIndex = -1;
		for (e = 0; e < src->num_lastedge; e++ ) {
			prevEdge = src->lastedge[e];
			/* look for the continuation of this path in nextEdge.*/
			prevIntvIndex = FindIntervalEndingAtPos(prevEdge->readinterval, prevEdge->multip, eq_read, readPos);
			if (prevIntvIndex >= 0) {
				curEdge = prevEdge;
				readPos = prevEdge->readinterval[prevIntvIndex].begin;
				intervalIndices[curIntv++] = prevIntvIndex;
				break;
			}
		} 
	} while (prevIntvIndex >= 0 and src->num_lastedge > 0);
	return *lenTracedPath = curIntv;
}

	
int FindIntervalIndex(EDGE *edge, int readIndex, int readPos) {
	int intv;
	for (intv = 0; intv < edge->multip; intv++) {
		if (readIndex == edge->readinterval[intv].eq_read) {
			if (edge->readinterval[intv].begin == readPos) {
				return intv;
			}
		}
	}
	printf("Interval not found for edge: %d multip: %d read: %d readPos %d\n",
				 edge->index, edge->multip, readIndex, readPos);
	assert(0);
	return -1;
}

int LookupVertexPathPos(PATH *paths, NODES *vertex, int readindex) {
	int p;
	for (p = 0; p < vertex->num_path; p++) {
		if (paths[vertex->path_index[p]].readindex == readindex) 
			return p;
	}
	return -1;
}

void AppendVertexPaths(NODES *enterEnd, NODES *exitBegin) {
	/* Add the paths in exitBegin to enterEnd.*/

	/* Count the number of paths that go from enter->exit*/
	enterEnd->path_index = (int*) realloc(enterEnd->path_index, sizeof(int) * (enterEnd->num_path + exitBegin->num_path));
	enterEnd->path_pos   = (int*) realloc(enterEnd->path_pos, sizeof(int) * (enterEnd->num_path + exitBegin->num_path));

	int curPath = enterEnd->num_path;
	int p;
	for (p = 0; p < exitBegin->num_path; p++) {
		enterEnd->path_index[curPath] = exitBegin->path_index[p];
		enterEnd->path_pos[curPath]   = exitBegin->path_pos[p];
		curPath++;
	}
	enterEnd->num_path = curPath;
}

void RepointPathEdges(PATH *paths,
											EDGE *newEdge, int pos, PathAddress *pathAddrs, int numPathAddrs, int beginOffset) {
	int p;
	for (p = 0; p < numPathAddrs; p++) {
		if (pathAddrs[p].path_index  >= 0) {
			paths[pathAddrs[p].path_index].edge[pos] = newEdge;
			if (pos == 0) 
				paths[pathAddrs[p].path_index].begin_length += beginOffset;
		}
	}
}


void AppendPathAddresses(NODES *vertex, PathAddress *pathAddrs, int numPathAddrs, int posOffset) {

	int p;
	int numOkPathAddrs = 0;
	if (numPathAddrs == 0)
		return;

	for (p = 0; p < numPathAddrs; p++) 
		if (pathAddrs[p].path_index >= 0) ++numOkPathAddrs;

	vertex->path_index = (int*) realloc(vertex->path_index, sizeof(int) * (vertex->num_path + numOkPathAddrs));
	vertex->path_pos   = (int*) realloc(vertex->path_pos, sizeof(int) * (vertex->num_path + numOkPathAddrs));

	int curPath = vertex->num_path;
	for (p = 0; p < numPathAddrs; p++) { 
		if (pathAddrs[p].path_index >= 0) {
			vertex->path_index[curPath] = pathAddrs[p].path_index;
			vertex->path_pos[curPath] = pathAddrs[p].path_pos + posOffset;
			curPath++;
		}
	}
	vertex->num_path = curPath;
}

void DeletePathRange(PATH *paths, int pathIndex, int pathDelStart, int pathDelEnd) {
	int p;
	int delLen = pathDelEnd - pathDelStart;
	for (p = pathDelEnd; p < paths[pathIndex].len_path; p++) {
		paths[pathIndex].edge[p-delLen] = paths[pathIndex].edge[p];
	}
	paths[pathIndex].len_path -= delLen;
	EDGE **delEdges = (EDGE**) ckalloc(sizeof(EDGE*) * paths[pathIndex].len_path);
	int numDelEdges = 0;
	int edgeAlreadyDel = 0;
	int de;
	for (p = pathDelStart; p < paths[pathIndex].len_path; p++ ) {
		edgeAlreadyDel = 0;
		for (de = 0; de < numDelEdges; de++ )
			if (delEdges[de] == paths[pathIndex].edge[p]) edgeAlreadyDel = 1;
		if (!edgeAlreadyDel) {
			delEdges[numDelEdges++] = paths[pathIndex].edge[p];
			UpdatePathIndices(paths[pathIndex].edge[p]->end, pathIndex, -delLen);
		}
	}
	free(delEdges);
}


void StoreEdgePathIntervalIndices(PATH *paths, PathAddress *pathAddrs, int numPathAddrs, 
																	int **edgeIntervalIndices ){

	(*edgeIntervalIndices) = (int*) ckalloc(sizeof(int) * numPathAddrs);
	int p;
	int edgeIndex[1];
	for (p = 0; p < numPathAddrs; p++) {
		assert(paths[pathAddrs[p].path_index].len_path == 1);
		LookupPathIntervalIndices(paths, pathAddrs[p].path_index, edgeIndex);
		(*edgeIntervalIndices)[p] = edgeIndex[0];
	}
}

void RemovePathsBetweenEdges(PATH *paths, PathAddress *pathAddrs, int numPathIndices, 
														 EDGE *enterEdge, EDGE *exitEdge,
														 int **enterEdgeIntervalIndices) {
	/*
		Input: pahts - all paths in the graph.
		pathIndices, numPathIndices  - a set of path indices to care
		about.
		enter, exit- the edges that enter, or exit a loop, respectively.
	*/
	int p;
	int pathIndex, pathPos;
	int edgeIndex, enterEdgeIndex, exitEdgeIndex;

	int numRemovedIntervals; 
	int *readIntervalIndices; 

	/* Will store the indices of the read intervals in enterEdge. */
	*enterEdgeIntervalIndices = (int*) ckalloc(sizeof(int) * numPathIndices);

	for (p = 0; p < numPathIndices; p++) {
		pathIndex = pathAddrs[p].path_index;
		pathPos   = pathAddrs[p].path_pos;
		readIntervalIndices = (int*) ckalloc(paths[pathIndex].len_path* sizeof(int));
		if (LookupPathIntervalIndices(paths, pathIndex, readIntervalIndices) == 0) {
			printf("TRACING path %d failed\n", pathIndex);
			//			assert(0);
			paths[pathIndex].len_path = 0;
			(*enterEdgeIntervalIndices)[p] = -1;
			pathAddrs[p].path_index = -1;
			continue;
		}

		/* We are what part of the path we start on as a parameter.*/
		enterEdgeIndex = pathPos;

		/* Store the index of the read interval in enterEdge. */
		(*enterEdgeIntervalIndices)[p] = readIntervalIndices[enterEdgeIndex];

		/* Find the index of the edge that is after the one we want to remove.*/
		exitEdgeIndex = FindIndexBoundByEdge(paths, pathIndex, pathPos, exitEdge);
		/*		
				printf("removing path %d pos %d to %d (len: %d) enterindex: %d exitIndex: %d\n",
				pathIndex, pathPos, exitEdgeIndex, paths[pathIndex].len_path, enterEdge->index,
					exitEdge->index);
		*/
		/* Leave the first interval intact, remove all intervals after that.*/		
		numRemovedIntervals = exitEdgeIndex - (enterEdgeIndex + 1);
		/*
		if (numRemovedIntervals > 0) {
			(*removedIntervals)[p] = (READINTERVAL*) ckalloc(sizeof(READINTERVAL) * numRemovedIntervals);
			(*numPathIntervals)[p] = numRemovedIntervals;
		}
		*/
		
		int readIntervalIndex;
		/*		readPos += enterEdge->length - VERTEX_SIZE;*/
		int numRemovedEdges = 0;
		int ri;
		int removedEdgeIndex;
		EDGE *removedEdge = paths[pathIndex].edge[enterEdgeIndex];
		// Copy the read intervals that will be removed.
		RemovePathIntervalFromVertex(removedEdge->begin, pathIndex, enterEdgeIndex);
		for (removedEdgeIndex = enterEdgeIndex + 1; 
				 removedEdgeIndex < exitEdgeIndex; ++removedEdgeIndex) {
			/* easy access to the removed edge.*/
			removedEdge = paths[pathIndex].edge[removedEdgeIndex];

			/*
				Find the read interval that corresponds to the path
				at positoin 'pathIndex' that passes through 'edge'.
			*/
			readIntervalIndex = readIntervalIndices[removedEdgeIndex];

			/* 
				 Grow the interval stored in 'enterindex' by the length
				 of the removed interval.  This is done to extend
				 the 'enterEdge' into the tandem repeat edges.
			*/
			GrowReadInterval(&(enterEdge->readinterval[(*enterEdgeIntervalIndices)[p]]),
											 &(removedEdge->readinterval[readIntervalIndex]));
			
			/*
				Store that read interval.  When removing intervals that
				are from paths that start in the whirl, they will later be 
				added to the 'enterEdge'.

			(*removedIntervals)[p][numRemovedEdges] = 
				removedEdge->readinterval[readIntervalIndex];
			*/
			++numRemovedEdges;
			/*
				Remove the interval. 
			*/
			MarkIndexedReadIntervalForRemoval(removedEdge, readIntervalIndex);

			/*
				Remove this path from the dest vertex.
			*/
			if (!RemovePathIntervalFromVertex(removedEdge->begin, pathIndex, removedEdgeIndex))
				printf("could not remove %d,%d from vertex: %d\n", pathIndex, removedEdgeIndex, removedEdge->begin->index);
		}

		/* 
			 Remove this path from the last vertex of the removed path.
		*/

		RemovePathIntervalFromVertex(removedEdge->end, pathIndex, removedEdgeIndex);

		if (numRemovedIntervals > 0) {
			DeletePathRange(paths, pathIndex, enterEdgeIndex+1, exitEdgeIndex);
		}

		free(readIntervalIndices);

		/* Fix the remainder of the path, if it continues past the end.*/
		if (numRemovedIntervals  > 0 and exitEdgeIndex < paths[pathIndex].len_path) {
			int remainingEdgeIndex;
			for (remainingEdgeIndex = exitEdgeIndex; 
					 remainingEdgeIndex < paths[pathIndex].len_path;
					 remainingEdgeIndex++) {
				UpdatePathIndices(paths[pathIndex].edge[remainingEdgeIndex]->end, pathIndex, -numRemovedIntervals);
			}
		}
		/* remove the edges marked for deletion.
			 But not for nwo.
			 for (removedEdgeIndex = enterEdgeIndex ; 
			 removedEdgeIndex < exitEdgeIndex; ++removedEdgeIndex) {
			 removedEdge = paths[pathIndex].edge[removedEdgeIndex];  
			 RemoveMarkedReadIntervals(removedEdge);
		*/
	}
}
 

void GrowReadInterval(READINTERVAL *srcInterval, READINTERVAL *interval) {
	srcInterval->length += interval->length - VERTEX_SIZE;
}


void AppendReadIntervalsToEdge(EDGE *srcEdge, int *srcEdgeIntervalIndices, int numIntv,
															 EDGE *destEdge, int offset) {

	/*
		First count the number of intervals that are good.
	*/
	int numGoodIntervals, intv;
	numGoodIntervals = 0;
	for (intv = 0; intv < numIntv; intv++) {
		if (srcEdgeIntervalIndices[intv] >= 0)  ++numGoodIntervals;
	}

	destEdge->readinterval = (READINTERVAL*) realloc(destEdge->readinterval,
																									 sizeof(READINTERVAL) * (destEdge->multip  + numGoodIntervals));

	int curInterval = 0;
	int intvIndex;
	for (intv = 0; intv < numIntv; intv++) {
		intvIndex = srcEdgeIntervalIndices[intv];
		if (intvIndex >= 0) {
			srcEdge->readinterval[intvIndex].offset += offset;
			destEdge->readinterval[destEdge->multip++] = srcEdge->readinterval[intvIndex];
		}
	}
}

int RemoveDuplicatePathIndices(PathAddress **pathIndices, int numIndices) {
	qsort(*pathIndices, numIndices, sizeof(PathAddress), CompareIndices);

	int numRemoved = 0;
	int p;
	for (p = 1; p < numIndices; p++ ) {
		if ((*pathIndices)[p-1].path_index == (*pathIndices)[p].path_index) {
			numRemoved++;
		}
		else {
			(*pathIndices)[p-numRemoved] = (*pathIndices)[p];
		}
	}
	
	if (numRemoved > 0) {
		numIndices -= numRemoved;
		*pathIndices = (PathAddress*) realloc(*pathIndices, sizeof(PathAddress) * (numIndices));
	}
	return numIndices;
}


int CountPathsStartingAtEdge(PATH *paths, EDGE *firstEdge, int *pathIndices) {
	int p;
	int pathIndex, pathPos;
	int numStartingPaths = 0;
	NODES *secondVertex = firstEdge->end;
	for (p = 0; p < secondVertex->num_path; p++ ){ 
		pathIndex = secondVertex->path_index[p];
		pathPos   = secondVertex->path_pos[p];
		/* If this vertex is not an end vertex on the path.*/
		if (pathPos == 0 and paths[pathIndex].edge[0] == firstEdge)
			if (pathIndices != NULL ) {
				pathIndices[numStartingPaths] = pathIndex;
				++numStartingPaths;
			}
	}
	return numStartingPaths;
}


int StorePathsStartingAtEdge(PATH *paths, EDGE *firstEdge, int **pathIndices) {
	int numStartingPaths = 0;
	numStartingPaths = CountPathsStartingAtEdge(paths, firstEdge, NULL);

	if (numStartingPaths > 0) {
		*pathIndices = (int*) ckalloc(sizeof(int) * numStartingPaths);
		numStartingPaths = CountPathsStartingAtEdge(paths, firstEdge, *pathIndices);
	}
	else {
		*pathIndices = NULL;
	}
	return numStartingPaths;
}

int CountPathsInEdge(PATH *paths, NODES *src, EDGE *edge, PathAddress *pathIndices) {
	int p;
	int numPaths = 0;
	for (p = 0; p < src->num_path; p++ ) {
		if (paths[src->path_index[p]].len_path == 1 and
				paths[src->path_index[p]].edge[0] == edge) {
			if (pathIndices != NULL) {
				pathIndices[numPaths].path_index = src->path_index[p];
				pathIndices[numPaths].path_pos   = 0;
			}
			numPaths++;
		}
	}
	return numPaths;
}

int StorePathsInEdge(PATH *paths, NODES *src, EDGE *edge, PathAddress **pathIndices) {
	int numPaths = 0;
	numPaths = CountPathsInEdge(paths, src, edge, NULL);
	*pathIndices = (PathAddress*) ckalloc(sizeof(PathAddress) * numPaths);
	CountPathsInEdge(paths, src, edge, *pathIndices);
	printf("edge: %d contains %d single paths\n", edge->index, numPaths);
	return numPaths;
}

int CountPathsThroughVertex(PATH *paths, EDGE *enterEdge, NODES *vertex, EDGE *exitEdge,
														PathAddress *pathAddresses) {
	int p;
	int pathIndex, pathPos;
	int numPassingPaths = 0;
	for (p = 0; p < vertex->num_path; p++ ){ 
		pathIndex = vertex->path_index[p];
		pathPos   = vertex->path_pos[p];
		/* If this vertex is not an end vertex on the path.*/
		if (pathPos > 0 and
				pathPos < paths[pathIndex].len_path) {
			if (paths[pathIndex].edge[pathPos-1] == enterEdge and
					paths[pathIndex].edge[pathPos] == exitEdge) {
				if (pathAddresses != NULL) {
					pathAddresses[numPassingPaths].path_index = pathIndex;
					pathAddresses[numPassingPaths].path_pos   = pathPos-1;
				}
				numPassingPaths++;
			}
		}
	}
	return numPassingPaths;
}
					

int StorePathsThroughVertex(PATH *paths,
														EDGE *enterEdge, NODES *vertex, EDGE *exitEdge, 
														PathAddress **pathAddresses) {
	printf("Storing paths %d-> %d(%d) %d\n",
				 enterEdge->index, vertex->index, vertex->num_path, exitEdge->index);

	int numPassingPaths;
	numPassingPaths = CountPathsThroughVertex(paths, enterEdge, vertex, exitEdge, NULL);
	printf("counted %d passing paths\n", numPassingPaths);
	if (numPassingPaths == 0) {
		*pathAddresses = NULL;
		return 0;
	}

	*pathAddresses = (PathAddress*) ckalloc(numPassingPaths * sizeof(PathAddress));
	numPassingPaths = CountPathsThroughVertex(paths, enterEdge, vertex, exitEdge, *pathAddresses);
	return numPassingPaths;
}


int EdgeIsShortSelfCycle(EDGE *edge, int maxLength) {
	return (edge->length < maxLength and edge->begin == edge->end);
}


int StraightenShortTandemRepeats(NODES ***vertices, int numVertices, 
																 PATH *paths, int cycleLength) {
	
	int v, e;
	EDGE *enterEdge, *repeatEdge, *exitEdge, *linkEdge;
	EDGE *balEnterEdge, *balRepeatEdge, *balExitEdge, *balLinkEdge;
	EDGE *repeatCopy, *balRepeatCopy;
	for (v = 0; v < numVertices; v++) {
		for (e = 0; e < (*vertices)[v]->num_nextedge; e++) {
			linkEdge = (*vertices)[v]->nextedge[e];
			if (EdgeConnectsShortTandemRepeats(linkEdge, cycleLength,
																				 &enterEdge, &repeatEdge, &exitEdge)) {

				if (repeatEdge != repeatEdge->bal_edge and
						enterEdge->bal_edge != exitEdge and
						exitEdge->bal_edge != enterEdge) {
					balLinkEdge = (*vertices)[v]->nextedge[e]->bal_edge;
					balEnterEdge = enterEdge->bal_edge;
					balExitEdge  = exitEdge->bal_edge;
					balRepeatEdge = repeatEdge->bal_edge;

					/*					if (enterEdge->index == 1811 or
							balExitEdge->index == 1811)
							continue;*/
					numVertices = StraightenShortDirectedCycle(vertices, numVertices, 
																										 paths,
																										 enterEdge, repeatEdge, linkEdge,
																										 exitEdge, &repeatCopy);
					
					int v2, p;
					for (v2 = 0; v2 < numVertices; v2++) {
						for (p = 0; p < (*vertices)[v2]->num_path; p++) {
							assert((*vertices)[v2]->path_pos[p] <= paths[(*vertices)[v2]->path_index[p]].len_path);
							if ((*vertices)[v2]->path_pos[p] < paths[(*vertices)[v2]->path_index[p]].len_path) {
								assert(paths[(*vertices)[v2]->path_index[p]].edge[(*vertices)[v2]->path_pos[p]] != repeatEdge);
								assert(paths[(*vertices)[v2]->path_index[p]].edge[(*vertices)[v2]->path_pos[p]] != linkEdge);
							}
						}
					}

					numVertices = StraightenShortDirectedCycle(vertices, numVertices,
																										 paths,
																										 balExitEdge, balRepeatEdge, balLinkEdge, 
																										 balEnterEdge, &balRepeatCopy);
				}
			}
		}
	}
}


int FindShortTandemRepeats(NODES **vertices, int numVertices,
													 PATH *paths, int cycleLength) {

	int v, e;
	EDGE *enterEdge, *repeatEdge, *exitEdge;
	for (v = 0; v < numVertices; v++) {
		for (e = 0; e < vertices[v]->num_nextedge; e++) {
			if (EdgeConnectsShortTandemRepeats(vertices[v]->nextedge[e], cycleLength,
																				 &enterEdge, &repeatEdge, &exitEdge)) {
				printf("Found whirl: enter: %d (%d) repeat %d (%d) link: %d (%d)  exit: %d (%d)\n",
							 enterEdge->index, enterEdge->length,
							 repeatEdge->index, repeatEdge->length,
							 vertices[v]->nextedge[e]->index, vertices[v]->nextedge[e]->length,
							 exitEdge->index, exitEdge->length);
			}
		}
	}
}
								
int sortedge(const void*a, const void*b) {
	return *((int*)a) < *((int*)b);
}

void CountN50(NODES **vertices, int numVertices) {
	/* Count the number of edges.*/
	int v, e;
	int numEdges = 0;
	int totalLength = 0;
	for (v = 0; v < numVertices; v++) {
		numEdges += vertices[v]->num_nextedge;
	}
	int *edgeLengths = (int*) ckalloc(numEdges *sizeof(int));
	int edgeIndex =0;
	for (v = 0; v < numVertices; v++) {
		for (e = 0; e < vertices[v]->num_nextedge; e++) {
			edgeLengths[edgeIndex++] = vertices[v]->nextedge[e]->length;
			totalLength += vertices[v]->nextedge[e]->length;
		}
	}
	qsort(edgeLengths, numEdges, sizeof(int), sortedge);
	
	int top10 = 10;
	if (numEdges > 20) {
		printf(" the top 10 edges are: \n");
		for (e = 0; e < 20; e+=2) {
			printf("%d ", edgeLengths[e]);
		}
		printf("\n");
	}
	int n50Len = 0;
	e = 0;
	while (n50Len < totalLength/2 and e < numEdges) {
		n50Len += edgeLengths[e];
		e++;
	}
	printf("the n50 length is: %d, out of %d edges\n", edgeLengths[e], numEdges);
	free(edgeLengths);
}

int EdgeConnectsShortTandemRepeats(EDGE *edge, int maxLength, 
																	 EDGE **enterEdge, EDGE **repeatEdge, EDGE **exitEdge) {
	NODES *dest, *src;
	dest = edge->end;
	src  = edge->begin;
	int e;
	if (dest->num_lastedge == 2 and dest->num_nextedge == 1 and
			src->num_lastedge == 1 and src->num_nextedge == 2 and
			dest->nextedge[0]->end == edge->begin and
			edge->length + dest->nextedge[0]->length - VERTEX_SIZE < maxLength) {
		for (e = 0; e < dest->num_lastedge; e++) {
			if (dest->lastedge[e] != edge) {
				*enterEdge = dest->lastedge[e];
				break;
			}
		}
		*repeatEdge = dest->nextedge[0];
		for (e = 0; e < src->num_nextedge; e++) {
			if (src->nextedge[e] != edge) {
				*exitEdge = src->nextedge[e];
				break;
			}
		}
		return 1;
	}
	else {
		return 0;
	}
}

int CountShortDirectedCycles(NODES **vertices, int numVertices, int maxCycleSize) {
	int v, e;
	for (v = 0; v < numVertices; v++ ) {
		for (e = 0; e < vertices[v]->num_nextedge; e++ ){
			if (vertices[v]->nextedge[e]->end ==
					vertices[v]->nextedge[e]->begin and
					vertices[v]->nextedge[e]->length < maxCycleSize) {
				printf("Edge: %d is a short cycle of length : %d\n",
							 vertices[v]->nextedge[e]->index,
							 vertices[v]->nextedge[e]->length);
			}
			else {
				NODES *dest;
				int destEdge;
				dest = vertices[v]->nextedge[e]->end;
				for (destEdge = 0; destEdge < dest->num_nextedge; destEdge++ ){
					if (dest->nextedge[destEdge]->end == vertices[v] and
							dest->nextedge[destEdge]->length +
							vertices[v]->nextedge[e]->length - VERTEX_SIZE < maxCycleSize) {
						printf("Edge: %d is a short cycle of length : %d\n",
									 vertices[v]->nextedge[e]->index,
									 vertices[v]->nextedge[e]->length);
					}
				}
			}
		}
	}
}


int count_edge_simp(NODES **vertex, int num_vertex, int **num_pa)
{
	/* Count the number of edges in the graph.
	 */
	int	i, j, k, l, m;
	int	l1, l2;

	for(i = 0; i < MAX_BRA; i ++)	{
		for(j = 0; j < MAX_BRA; j ++)	{
			num_pa[i][j] = 0;
		}
	}

	l1 = l2 = 0;
	for(i = 0; i < num_vertex; i ++)	{
		if(vertex[i] -> num_lastedge < 6 && vertex[i] -> num_nextedge < 6)	{
			num_pa[vertex[i] -> num_lastedge][vertex[i] -> num_nextedge] ++;
		} else if(vertex[i] -> num_lastedge < 6)	{
			num_pa[vertex[i] -> num_lastedge][6] ++;
		} else if(vertex[i] -> num_nextedge < 6)	{
			num_pa[6][vertex[i] -> num_nextedge] ++;
		} else	{
			num_pa[6][6] ++;
		}
		/*		printf("vertex: %d (%d) next: %d nlast: %d \n",
					 i, vertex[i]->index, 
					 vertex[i]->num_nextedge, vertex[i]->num_lastedge);*/
		l1 += vertex[i] -> num_nextedge;
		l2 += vertex[i] -> num_lastedge;
		//		printf("vertex: %d num next: %d num last %d\n", i, vertex[i]->num_nextedge, vertex[i]->num_lastedge);
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			vertex[i] -> nextedge[j] -> visit = 0;
		}
	}
	if(l1 != l2)	{
		printf("edge not balanced  %d %d\n", l1, l2);
		//		exit(-1);
	}
	return(l1);
}

int count_edge(NODES **vertex, int num_vertex, int **num_pa)
{
	int	i, j, k, l, m;
	int	l1, l2;
	int	ave, multip;
	EDGE	**maxedge;

	nsuper = 0;

	for(i = 0; i < MAX_BRA; i ++)	{
		for(j = 0; j < MAX_BRA; j ++)	{
			num_pa[i][j] = 0;
		}
	}

	l1 = l2 = 0;
	for(i = 0; i < num_vertex; i ++)	{
		if(vertex[i] -> num_lastedge < 6 && vertex[i] -> num_nextedge < 6)	{
			num_pa[vertex[i] -> num_lastedge][vertex[i] -> num_nextedge] ++;
		} else if(vertex[i] -> num_lastedge < 6)	{
			num_pa[vertex[i] -> num_lastedge][6] ++;
		} else if(vertex[i] -> num_nextedge < 6)	{
			num_pa[6][vertex[i] -> num_nextedge] ++;
		} else	{
			num_pa[6][6] ++;
		}
		l1 += vertex[i] -> num_nextedge;
		l2 += vertex[i] -> num_lastedge;
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			vertex[i] -> nextedge[j] -> visit = 0;
		}
	}
	if(l1 != l2)	{
		printf("edge not balanced %d %d\n", l1, l2);
		//		exit(-1);
	}

	maxedge = (EDGE **) ckalloc(l1 * sizeof(EDGE *));
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			if(vertex[i] -> nextedge[j] -> start_cover > 1)	{
				if(vertex[i] -> nextedge[j] -> visit == 0)	{
					getmaxedge(vertex[i] -> nextedge[j], maxedge);
				}
			}
		}
		/*		if(numtangle[nsuper] > 0)	nsuper ++;*/
	}

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			vertex[i] -> nextedge[j] -> visit = 0;
		}
	}

	for(i = 0; i < nsuper; i ++)	{
		ave = maxedge[i] -> length * maxedge[i] -> start_cover;
		multip = maxedge[i] -> start_cover;
		m = backtanglelen(maxedge[i], &ave, &multip);
		m += maxedge[i] -> length;
		m += forwardtanglelen(maxedge[i], &ave, &multip);
		maxlength[i] = m;
		avelength[i] = ave / multip;
	}

	free((void **) maxedge);
	return(l1);
}

int backtanglelen(EDGE *edge, int *ave, int *multip)
{
	int	i, j, k, l;
	NODES	*vertex;

	vertex = edge -> begin;
	l = 0;
	for(j = 0; j < vertex -> num_lastedge; j ++)	{
		if(vertex -> lastedge[j] -> visit == 0)	{
			vertex -> lastedge[j] -> visit = 1;
		} else	{
			continue;
		}
		if(vertex -> lastedge[j] -> start_cover > 1)	{
			*ave += (vertex -> lastedge[j] -> length - 1) * vertex -> lastedge[j] -> start_cover;
			*multip += vertex -> lastedge[j] -> start_cover;
			k = vertex -> lastedge[j] -> length + backtanglelen(vertex -> lastedge[j], ave, multip) - 1;
			if(k > l)	 l = k;
		}
	}
	return(l);
}

int forwardtanglelen(EDGE *edge, int *ave, int *multip)
{
	int	i, j, k, l;
	NODES	*vertex;

	vertex = edge -> end;
	l = 0;
	for(j = 0; j < vertex -> num_nextedge; j ++)	{
		if(vertex -> nextedge[j] -> visit == 0)	{
			vertex -> nextedge[j] -> visit = 1;
		} else	{
			continue;
		}
		if(vertex -> nextedge[j] -> start_cover > 1)	{
			*ave += (vertex -> nextedge[j] -> length - 1) * vertex -> nextedge[j] -> start_cover;
			*multip += vertex -> nextedge[j] -> start_cover;
			k = vertex -> nextedge[j] -> length + forwardtanglelen(vertex -> nextedge[j], ave, multip) - 1;
			if(k > l)	 l = k;
		}
	}
	return(l);
}

void count_multip(NODES **vertex, int num_vertex)
{
	int	i, j, k, l, n1, n2;
	EDGE	*edge, *edge0;
	int	max_red, l_red;
	NODES	*v;
	int	num_unb, num_unb_copy;

	num_unb_copy = 2 * num_vertex + 1;
	num_unb = num_vertex * 2;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			edge -> start_cover = 1;
		}
	}
	for(i = 0; i < num_vertex; i ++)	{
		if(vertex[i] -> num_lastedge == 0)	{
			edge = vertex[i] -> nextedge[0];
			v = edge -> end;
			if(v -> num_nextedge == 1)	{
				v -> nextedge[0] -> start_cover = v -> num_lastedge;
			}
		} else if(vertex[i] -> num_nextedge == 0)	{
			edge = vertex[i] -> lastedge[0];
			v = edge -> begin;
			if(v -> num_lastedge == 1)	{
				v -> lastedge[0] -> start_cover = v -> num_nextedge;
			}
		}
	}
	while(num_unb_copy > num_unb)	{
		num_unb_copy = num_unb;
		num_unb = 0;
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_nextedge != 0 && vertex[i] -> num_lastedge != 0)	{
				num_unb += abs(count_bal(vertex[i]));
			}
		}
		max_red = 0;
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_nextedge != 0 && vertex[i] -> num_lastedge != 0)	{
				for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
					edge = vertex[i] -> nextedge[j];
					n1 = count_bal(edge -> begin);
					n2 = count_bal(edge -> end);
					if(n1 > 0 && n2 < 0)	{
						if(max_red < 2)		{
							edge0 = edge;
							l_red = 1;
							max_red = 2;
						}
					} else if(n1 < 0 && n2 > 0 && edge -> start_cover > 1)	{
						if(max_red < 2)		{
							edge0 = edge;
							l_red = 0;
							max_red = 2;
						}
					}
				}
			}
			if(max_red == 2)	break;
		}
		if(max_red == 2)	{
			if(l_red == 1)	{
				edge0 -> start_cover ++;
			} else	{
				edge0 -> start_cover --;
			}
		}
	}
}

int count_tangle(NODES **vertex, int num_vertex, int disttangle[][7])
{
	int	i, j, k, l, n1, n2;
	EDGE	*edge, *edge0;
	int	num_tangle;

	for(i = 0; i < 8; i ++)	{
		for(j = 0; j < 7; j ++)	{
			disttangle[i][j] = 0;
		}
	}

	num_tangle = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			if(edge -> start_cover <= 6)	{
				n1 = edge -> start_cover - 1;
			} else	{
				n1 = 6;
			}
			if(edge -> length < 500)	{
				n2 = 0;
			} else if(edge -> length < 1000)	{
				n2 = 1;
			} else if(edge -> length < 5000)	{
				n2 = 2;
			} else if(edge -> length < 10000)	{
				n2 = 3;
			} else	{
				n2 = 4;
			}
			disttangle[n1][n2] ++;
			disttangle[7][n2] ++;
			disttangle[n1][5] ++;
			disttangle[7][5] ++;
			if(edge -> start_cover > 1)	num_tangle ++;
		}
	}
	return(num_tangle);
}

int count_bal(NODES *vertex)
{
	int	i, j, k, l, n1, n2;

	n1 = n2 = 0;
	for(j = 0; j < vertex -> num_lastedge; j ++)	{
		n1 += vertex -> lastedge[j] -> start_cover;
	}
	for(j = 0; j < vertex -> num_nextedge; j ++)	{
		n2 += vertex -> nextedge[j] -> start_cover;
	}

	return(n1 - n2);
}

void getmaxedge(EDGE *edge, EDGE **maxedge)
{
	int	i, j, k, l;

	edge -> visit = 1;
	numtangle[nsuper] ++;
	if(edge -> start_cover > maxmultip[nsuper])	{
		maxmultip[nsuper] = edge -> start_cover;
		mlength[nsuper] = edge -> length;
		maxedge[nsuper] = edge;
	}
	for(j = 0; j < edge -> end -> num_nextedge; j ++)	{
		if(edge -> end -> nextedge[j] -> start_cover > 1)	{
			if(edge -> end -> nextedge[j] -> visit == 0)	{
				getmaxedge(edge -> end -> nextedge[j], maxedge);
			}
		}
	}
	for(j = 0; j < edge -> begin -> num_lastedge; j ++)	{
		if(edge -> begin -> lastedge[j] -> start_cover > 1)	{
			if(edge -> begin -> lastedge[j] -> visit == 0)	{
				getmaxedge(edge -> begin -> lastedge[j], maxedge);
			}
		}
	}
}
