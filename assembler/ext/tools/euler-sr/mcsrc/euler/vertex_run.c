#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>
#include <assert.h>

int vertex_run(NODES **vertex, int num_vertex, PATH *path, int num_path, int small_cover, int num_seq);
int run1vertex(NODES *curVertex, NODES **vertex, int num_vertex, PATH *path, int num_path, EDGE *inedge, int num_seq, 
							 EDGE **newedge);
int createvertex(NODES ***vertexListPtr, int num_vertex, NODES *curVertex, PATH *path, int num_path);
int chkcross(NODES *curVertex, PATH *path, int num_path);
int chkforbid(EDGE *edge1, EDGE *edge2, PATH *path);


void CopyVertex(NODES *src, NODES **dest) {
	if (*dest == NULL) 
		*dest = (NODES*) ckalloc(sizeof(NODES));
	
	int i;

	/* Copy edge structure.*/
	(*dest)->num_nextedge = src->num_nextedge;
	(*dest)->nextedge = (EDGE**) ckalloc(sizeof(EDGE*)*(*dest)->num_nextedge);
	for (i = 0; i < src->num_nextedge; i++ ) {
		(*dest)->nextedge[i] = src->nextedge[i];
		(*dest)->nextedge[i]->deleted = 0;
	}
	(*dest)->num_lastedge = src->num_lastedge;
	(*dest)->lastedge     = (EDGE**) ckalloc(sizeof(EDGE)*(*dest)->num_lastedge);
	for (i = 0; i < src->num_lastedge; i++ ) {
		(*dest)->lastedge[i] = src->lastedge[i];
		(*dest)->lastedge[i]->deleted = 0;
	}

	/* Copy paths. */
	(*dest)->num_path = src->num_path;
	(*dest)->path_index = (int*) ckalloc(sizeof(int) * src->num_path);
	(*dest)->path_pos   = (int*) ckalloc(sizeof(int) * src->num_path);
	for (i = 0; i < src->num_path; i++) {
		(*dest)->path_index[i] = src->path_index[i];
		(*dest)->path_pos[i]   = src->path_pos[i];
	}
	
	/* Copy ancillary information that I don't think will have any effect here.*/
	(*dest)->npos = src->npos;
	(*dest)->nlinks = src->nlinks;
	(*dest)->index = src->index;
	(*dest)->comp = src->comp;
}
	
void FreeVertex(NODES *vertex) {
	//	free((void*) vertex->nextedge);
	free((void*) vertex->lastedge);
	free((void*) vertex->path_index);
	free((void*) vertex->path_pos);
	free((void*) vertex);
}

void CopyEdge(EDGE *src, EDGE**dest) {
	if (*dest == NULL) {
		*dest = (EDGE*) ckalloc(sizeof(EDGE));
		(*dest)->deleted = 0;
	}
	else {
		assert((*dest)->deleted == 0);
	}
	(*dest)->index = src->index;
	(*dest)->begin = src->begin;
	(*dest)->end   = src->end;
	(*dest)->start_cover = src->start_cover;
	(*dest)->bal_edge = src->bal_edge;
	
	/* Copy sequence for hte edge.*/
	(*dest)->length = src->length;
	if (src->seq == NULL)
		(*dest)->seq = NULL;
	else {
		(*dest)->seq = (char*) ckalloc(sizeof(char)*src->length);
		memcpy((*dest)->seq, src->seq, src->length);
	}
	
	/* Copy the read intervals */
	int i;
	(*dest)->multip = src->multip;
	(*dest)->readinterval = (READINTERVAL*) ckalloc(sizeof(READINTERVAL)*src->multip);
	for (i = 0; i< src->multip; i++) {
		(*dest)->readinterval[i].begin = src->readinterval[i].begin;
		(*dest)->readinterval[i].eq_read = src->readinterval[i].eq_read;
		(*dest)->readinterval[i].length = src->readinterval[i].length;
		(*dest)->readinterval[i].offset = src->readinterval[i].offset;
		(*dest)->readinterval[i].cov = src->readinterval[i].cov;
		(*dest)->readinterval[i].next = src->readinterval[i].next;
	}
}

void FreeEdge(EDGE* edge) {
	/*
	free((void*)edge->seq);
	free((void*)edge->readinterval);
	free((void*)edge);
	*/
	/*	printf("(FreeEdge) deleting edge: %d\n", edge);*/
	edge->deleted = 1;
}


int LookupPathIndex(NODES *vertex, int pathIndex) {
	int p;
	for (p = 0; p < vertex->num_path; p++) {
		if (vertex->path_index[p] == pathIndex)
			break;
	}
	return p;
}

void InsertEdge(PATH *paths, int pathIndex, int pathEdge, EDGE *newEdge) {
	/* 
		 Insert 'edge' after 'pathEdge'.
	*/
	EDGE **tmpPathEdges;
	/* Make a new edge list that has space for the array. */
	tmpPathEdges = (EDGE**) ckalloc(sizeof(EDGE*)*paths[pathIndex].len_path + 1);
	int e;
	int destVertexPathIndex, srcVertexPathIndex;

	/* Copy over the old path with a gap for the end. */
	if (paths[pathIndex].len_path > 0) 
		for (e = 0; e <= pathEdge; e++) 
			tmpPathEdges[e] = paths[pathIndex].edge[e];

	if (pathEdge + 1 < paths[pathIndex].len_path) {
		for (e = pathEdge + 1; e < paths[pathIndex].len_path; e++ ) {
			/* Insert the shifted edges. */
			tmpPathEdges[e+1] = paths[pathIndex].edge[e];
			
			/* Now the indices of the vertices along this path are increased
				 by 1. Move along the path and update the vertex indices.
			*/
			destVertexPathIndex = LookupPathIndex(paths[pathIndex].edge[e]->end, pathIndex);
			assert(destVertexPathIndex >= 0);
			paths[pathIndex].edge[e]->end->path_pos[destVertexPathIndex]++;
		}
	}
	/* add the new edge */
	assert(pathEdge+1 < paths[pathIndex].len_path+ 1);
	tmpPathEdges[pathEdge+1] = newEdge;

	/* fix the index of the destination of the new edge*/
	destVertexPathIndex = LookupPathIndex(newEdge->end, pathIndex);
	assert(destVertexPathIndex >= 0);
	newEdge->end->path_pos[destVertexPathIndex]++;

	/* Replace the path. */
	free((void**) paths[pathIndex].edge);
	paths[pathIndex].edge = tmpPathEdges;

	/* 
	
	 The dest vertex of the new edge does not have the interval
		 corresponding to this edge, add it.  The location of the end
		 vertex along the length of the path is one after the current edge.
	
		 AppendPathInterval(path[pathIndex].edge[pathIndex]->end, pathIndex, pathEdge+1);
	*/
	/* 
		 Finally update the length of the path to hold the new edge.
	*/
	paths[pathIndex].len_path++;
	
}

void PrintPaths(NODES *vertex, PATH *paths) {
	int p;
	int pathIndex, pathEdge;
	int pathPos;
	printf("vertex: %d contains paths: \n", vertex->index);
	for (p = 0; p < vertex->num_path; p++ ) {
		pathIndex = vertex->path_index[p];
		printf("path: %d : ", pathIndex);
		for (pathPos = 0; pathPos < paths[pathIndex].len_path; pathPos++) {
			printf("src %d edge %d (%d) > ", paths[pathIndex].edge[pathPos]->begin->index,
						 paths[pathIndex].edge[pathPos],
						 paths[pathIndex].edge[pathPos]->index);

		}
		printf("dest %d\n", paths[pathIndex].edge[pathPos-1]->end->index);
	}
}

void RestorePaths(EDGE *newEdge, EDGE *inEdge, EDGE *outEdge, NODES *origVertex, PATH *paths, int vertexSize) {
	/* Make the paths that were re-routed to newedge and put them back in to 'inEdge' and 
		 'outEdge'.  
	*/
	NODES *srcVertex, *destVertex;
	srcVertex  = newEdge->begin;
	destVertex = newEdge->end;

	int p;
	int pathIndex, pathEdge;
	int indexInOrig;
	for (p = 0; p < srcVertex->num_path; p++) {
		pathIndex = srcVertex->path_index[p];
		pathEdge  = srcVertex->path_pos[p];

		if (paths[pathIndex].edge[pathEdge] == newEdge) {
			
			/* 
				 The path goes throug the new edge, so we need to update i
				 to use the old 'inEdge->origVertex->outEdge'.
			*/
			if (pathEdge == 0) {
				/* 
					 The path begins in new edge.  There are many possible cases here.
					 The path either begins and ends in newEdge, if it does, then it either:
					    A: Began and ended in 'inEdge'
							B: Began and ended in 'outEdge'
							C: Began in 'inEdge', and passed through 'outEdge'
							D: Began in 'inEdge', and ended in 'outEdge'.
							
							If it was case A, C, or D, then origVertex is the second vertex along the path (path_pos == 1).
							If it was case A, B, the length of the new and original paths are 1.
							If it was case D, the length of the new path is 1, and the original path is 2.
							If it was case B, then origVertex is the first vertex in the path (path_pos == 0).
				*/
				indexInOrig = LookupPathIndex(origVertex, pathIndex);
				assert(indexInOrig < origVertex->num_path);

				if (origVertex->path_pos[indexInOrig] == 0) {
					/* 
						 The path starts at the original vertex, so skips inEdge and starts in outEdge.
					*/
					paths[pathIndex].edge[pathEdge] = outEdge;

					/* 
						 Add this interval to the edges corresponding to the out
						 edge.
					*/
					
					/* 
						 Since we're subtracting inEdge from newEdge, this changes the distance to the beginning
						 of the edge.
					*/
					paths[pathIndex].begin_length = paths[pathIndex].begin_length - inEdge->length + origVertex->size;
				}
				else {
					/* The paths starts before this vertex, so it should start in inedge*/
					assert(origVertex->path_pos[indexInOrig] == 1);
					
					paths[pathIndex].edge[pathEdge] = inEdge;


					/* Either the path starts in 'inEdge' and passes through 'outEdge', or
						 'inEdge' was long enough to contain the entire path.  If 'inEdge'
						 contains the entire path, the 'end_length' field of path[pathIndex]
						 needs to be updated to reflect the shorter length of 'inEdge', and 
						 no new edge is needed on the path.
						 Otherwise, even if the path ends in 'outEdge', the length of the
						 'end_length' does not need to be changed.
					*/
					if (paths[pathIndex].len_path == 1 and /* there is only 1 edge in the new path*/
							paths[pathIndex].end_length > outEdge->length - vertexSize /* the path ends before the inEdge ends */
							) 
						/* 
							 The edge ends in 'inEdge', update the end_length. No need
							 to add 'outEdge' to the path.
						 */
						paths[pathIndex].end_length = paths[pathIndex].end_length - outEdge->length + vertexSize;
					else
						/* The edge passes through 'inEdge', so we need to add an outEdge to the path.*/
						InsertEdge(paths, pathIndex, pathEdge, outEdge);
				}
			}
			else if (pathEdge == paths[pathIndex].len_path - 1) {
				/* 
					 The path ends in 'newEdge'. That means the path may have:
					 A: Began and ended in 'inEdge'
					 B: Began and ended in 'outEdge'
					 C: Ended in 'inEdge',
					 D: Passed through 'inEdge', and ended in 'outEdge'.
					 
					 The first two cases were already handled, since it's also the case that 
					 paths[pathIndex].len_path == 0, so all that needs to be taken care of are
					 when the path started before 'inEdge'.

					 If it ended in inEge, the length of the path does not
					 change, so we can simply update the path.  If it ended in
					 outedge, we need replace inEdge as the end of the path, and add outEdge to the path.
				*/
				indexInOrig = LookupPathIndex(origVertex, pathIndex);
				if (origVertex->path_pos[indexInOrig] == paths[pathIndex].len_path) {
					/* 
						 The path ends before orig vertex, so it ends in inEdge.
					*/
					paths[pathIndex].edge[pathEdge] = inEdge;
					paths[pathIndex].end_length = paths[pathIndex].end_length - outEdge->length + origVertex->size;
				}
				else {
					paths[pathIndex].edge[pathEdge] = inEdge;
					/* 
						 If the path begins after inEdge.
					*/
					InsertEdge(paths, pathIndex, pathEdge, outEdge);
				}
			}
			else {
				/*
					The path passes through 'inEdge' and 'outEdge'.  
					replace newedge with 'inEdge', and add 'outEdge' to the path.
				*/
				paths[pathIndex].edge[pathEdge] = inEdge;
				
				InsertEdge(paths, pathIndex, pathEdge, outEdge);
			}
		}
	}
}


void UndetachEdges(NODES *begin, NODES*end, EDGE*newEdge,
									 PATH *paths,
									 int deletedEdges[2],
									 EDGE* inEdgeCopy, NODES *vertexCopy, EDGE *outEdgeCopy,
									 EDGE* inEdge, NODES*vertex, EDGE*outEdge ) {
	
	/* 
		 Steps to undetach a previously detached edge.
		 there used to be:
		 inEdge.src -inEdge-> vertex -outEdge-> outEdge.dest
		 
		 this is replaced by :
		 inEdge.src -- newEdge --> outEdge.dest

		 So to move back we need to:
		   Replace connectivity:
			 Remove newedge from begin and end.

		   Add inEdge to inEdge.src's out edges
			 Add inEdge to vertex's inEdges
			 Add outEdge to vertex's outEdges
			 Add outEdge to outEdge.dest's in edges.
			 
			 Fix the paths:
			 Modify the paths that were re-routed to newEdge to go back through
			 inEdge, vertex, and outEdge.
			 Add the paths back to vertex.
	*/

	int newEdgeInIndex, newEdgeOutIndex;
	newEdgeOutIndex = searcherase(begin->nextedge, newEdge, begin->num_nextedge);
	assert(newEdgeOutIndex < begin->num_nextedge);
	erasenext(begin, newEdgeOutIndex);
	
	newEdgeInIndex = searcherase(end->lastedge, newEdge, end->num_lastedge);
	assert(newEdgeInIndex < end->num_lastedge);
	eraselast(end, newEdgeInIndex);

	
	EDGE *origInEdge, *origOutEdge;
	origInEdge  = inEdge;
	origOutEdge = outEdge;
	if (deletedEdges[0]) {
		inEdge = (EDGE*) ckalloc(sizeof(EDGE));
		/* The in edge was erased. In order to undetach the edges, we need to 
			 make the edge again.*/
		add_lastedge(vertex, inEdge);
		add_nextedge(begin, inEdge);
	}
	if (deletedEdges[1]) {
		outEdge = (EDGE*) ckalloc(sizeof(EDGE));
		add_nextedge(vertex, outEdge);
		add_lastedge(end, outEdge);
	}
	/*	printf("copying %d into %d\n", inEdgeCopy, inEdge);*/
	
	CopyEdge(inEdgeCopy, &inEdge);
	/*	printf("copying %d into %d\n", outEdgeCopy, outEdge);*/
	CopyEdge(outEdgeCopy, &outEdge);

	/* Fix balanced edges if need be. */
	if (inEdgeCopy->bal_edge == origInEdge and
			deletedEdges[0])
		inEdge->bal_edge = inEdge;
	else if (inEdgeCopy->bal_edge == origOutEdge and
					 deletedEdges[1])
		inEdgeCopy->bal_edge = outEdge;

	if (outEdgeCopy->bal_edge == origOutEdge and
			deletedEdges[1])
		outEdge->bal_edge = outEdge;
	else if (outEdgeCopy->bal_edge == origInEdge and
					 deletedEdges[0]) 
		outEdge->bal_edge = inEdge;

	


	/* Copy back the paths from the vertex copy.
		 The entire vertex cannot be copied back because
		 some inEdge/outEdge may have been deleted.*/

	/* Copy paths. */
	if (vertex->num_path > 0) {
		free((void**) vertex->path_index);
		free((void**) vertex->path_pos);
	}

	vertex->num_path = vertexCopy->num_path;
	vertex->path_index = (int*) ckalloc(sizeof(int) * vertexCopy->num_path);
	vertex->path_pos   = (int*) ckalloc(sizeof(int) * vertexCopy->num_path);
	int i;
	for (i = 0; i < vertexCopy->num_path; i++) {
		vertex->path_index[i] = vertexCopy->path_index[i];
		vertex->path_pos[i]   = vertexCopy->path_pos[i];
	}
	
	/* Copy ancillary information that I don't think will have any effect here.*/
	vertex->npos = vertexCopy->npos;
	vertex->nlinks = vertexCopy->nlinks;
	vertex->index = vertexCopy->index;
	vertex->comp = vertexCopy->comp;

	/* Connectivity.  Put back old edges.*/
	/*	a bug: add_nextedge(begin, inEdge);*/


	/*
		Replace the modified 'vertex' with its original, 'vertexCopy'.
	*/

	int in, out;
	for (in = 0; in < vertex->num_lastedge; in++ ) {
		/*		printf("writing to in: %d &(in->end): %d, in->end %d\n", 
					 vertex->lastedge[in], &(vertex->lastedge[in]->end), 
					 vertex->lastedge[in]->end);*/
		vertex->lastedge[in]->end = vertex;
	}
	for (out = 0; out < vertex->num_nextedge; out++) {
		vertex->nextedge[out]->begin = vertex;
	}
	int vs = VERTEX_SIZE;

	RestorePaths(newEdge, inEdge, outEdge, vertex, paths, vs);

	FreeVertex(vertexCopy);
	assert(inEdge->bal_edge != NULL);
	assert(outEdge->bal_edge != NULL);
	int tin, tout;
	for (tin = 0; tin < vertex->num_lastedge; tin++ )
		assert(vertex->lastedge[tin]->bal_edge != NULL);
	
	for (tout = 0; tout < vertex->num_nextedge; tout++) 
		assert(vertex->nextedge[tout]->bal_edge != NULL);
}

int ShortXCut(NODES *vertex, EDGE* edge,
							PATH *path, int numPaths,
							int **savedEndPathIndices, int **savedEndPathPos,
							int **savedEndIntervalGap,
							READINTERVAL **savedEndPathIntervals) {
	int vertexPath;
	int pathIndex, pathPos;
	int numEndPaths = 0;
	// Remove path ends that start in 'vertex', and pass through 'edge'.
	int *endPathIndices, *endPathPos, *endIntervalGap;
	READINTERVAL *endPathIntervals;
	endPathIndices = (int*) ckalloc(sizeof(int) * numEndPaths);
	endPathPos     = (int*) ckalloc(sizeof(int) * numEndPaths);
	endIntervalGap = (int*) ckalloc(sizeof(int) * numEndPaths);
	endPathIntervals = (READINTERVAL*) ckalloc(sizeof(READINTERVAL)*numEndPaths);
	for (vertexPath = 0; vertexPath < vertex->num_path; vertexPath++) { 
		pathIndex = vertex->path_index[vertexPath];
		pathPos   = vertex->path_pos[vertexPath];
		if ((pathPos < path[pathIndex].len_path and 
				 path[pathIndex].edge[pathPos] == edge) and
				(pathPos == 0 or pathPos == path[pathIndex].len_path-1)) {
			// This path either begins or ends in 'edge'
			// Look to see if 'edge' is short.
			endPathIndices[numEndPaths] = pathIndex;
			endPathPos[numEndPaths] = pathPos;
			++numEndPaths;
		}
		else {
			if (numEndPaths > 0) { // if not needed, but saves computation
				vertex->path_index[vertexPath - numEndPaths] = vertex->path_index[vertexPath];
				vertex->path_pos[vertexPath - numEndPaths] = vertex->path_pos[vertexPath];
			}
		}
	}

	// Now remove all read intervals corresponding to removed path segments.

	int ri;
	int readIntvIndex;
	int ep;
	for (ep = 0; ep < numEndPaths; ep++ ) {
		readIntvIndex = LookupReadintervalIndex(edge, path[endPathIndices[ep]].readindex-1, 
																						endPathPos[endPathIndices[ep]]);
		assert(readIntvIndex >= 0);
		endPathIntervals[ep] = edge->readinterval[readIntvIndex];
		RemoveIndexedReadinterval(edge, readIntvIndex);
	}
	
	// Adjust the length of each path.
	for (ep = 0; ep < numEndPaths; ep++) {
		path[endPathIndices[ep]].len_path--;
	}

	// Now update the indices of the remainder of the path for each 
	// path that had a segment removed.
	int intv;
	PATH *endPath;
	for (ep = 0; ep < numEndPaths; ep++) {
		if (endPathPos[ep] == 0) {
			// the segment was removed from the beginnign of the
			// path.  We need to update the path positions for 
			// the reaminder of the path.
			endPath = &path[endPathIndices[ep]];
			for (intv = 0; intv < path->len_path; intv++) {
				UpdatePathIndices(path->edge[intv]->begin, endPathIndices[ep], -1);
			}
			UpdatePathIndices(path->edge[intv]->end, endPathIndices[ep], -1);

			// This path now begins exactly at the beginning of the edge.
			// Record the old one in case it needs to be restored.
			endIntervalGap[ep] = path->begin_length;
			path->begin_length = 0;
		}
		else {
			// This path now ends exactly at the
			// end of an edge. Record that
			endIntervalGap[ep] = path->end_length;
			path->end_length = 0;
		}
					
	}

	// All of the intervals that have been removed are saved, in case
	// they need to be restored if the vertex_run doesn't work.
	*savedEndPathIndices = endPathIndices;
	*savedEndPathPos     = endPathPos;
	*savedEndPathIntervals = endPathIntervals;
	return numEndPaths;
}

void RestoreXCut(EDGE* edge, 
								 PATH *path, int numPaths,
								 int *endPathIndices, int *endPathPos, int *endIntervalGap, 
								 READINTERVAL *endPathIntervals, 
								 int numEndPaths) {

	// Given a list of intervals that have been removed from 
	// the graph, add them back.

	int ep;
	int intv;
	int pathIndex;
	int numVertexIntervals;
	// First add back all the paths at the beginning	
	for (ep = 0; ep <  numEndPaths; ep++) {
		pathIndex = endPathIndices[ep];
		if (endPathPos[ep] == 0) {
			// Add 'edge' back to the beginning of the path.
			// Fortunately, since no memory was freed, we can 
			// just shift everything back into place.
			for (intv = path[pathIndex].len_path; intv > 0; intv--) {
				path[pathIndex].edge[intv] = path[pathIndex].edge[intv-1];
				UpdatePathIndices(path[pathIndex].edge[intv]->end, pathIndex, 1);
			}
			path[pathIndex].len_path++;
			path[pathIndex].edge[0] = edge;
			edge->begin->path_index[edge->begin->num_path] = pathIndex;
			edge->begin->path_pos[edge->begin->num_path] = 0;
			edge->begin->num_path++;
			path[pathIndex].begin_length = endIntervalGap[ep];
		}
	}
	// Now add back all paths at the end.  This is easy because
	// no path indices need to be modified.

	for (ep = 0; ep < numEndPaths; ep++ ){
		pathIndex = endPathIndices[ep];
		if (endPathPos[ep] != 0) {
			// Simply append this path interval to the 'end' 
			// of 'edge'.
			edge->end->path_index[edge->end->num_path] = pathIndex;
			edge->end->path_pos[edge->end->num_path] = endPathPos[ep];
			edge->end->num_path++;
		}
		path[pathIndex].end_length = endIntervalGap[ep]; 
	}

	// Finally add the removed read intervals back to 'edge'
	// since they were deleted.  Since they are appended
	// to the back, and order of the readintervals is important,
	// the read intervals need to be resorted.

	for (ep = 0; ep < numEndPaths; ep++) {
		edge->readinterval[++edge->multip] = endPathIntervals[ep];
	}
	sortreadinterval(edge->readinterval, edge->multip);
}

int ShortEdgeCutTransform(NODES **vertices, int numVertices,
													EDGE *edge, EDGE *balEdge,
													PATH *path, int numPaths, int *posttransNumVertices,
													int numSeq) {
	// Look to see if an X-cut will allow an additional equivalent transformation.
	// If it does, keep it.  If it doesn't restore the paths.
	int *edgeEndPathIndices, *edgeEndPathPos, 
		*balEndPathIndices, *balEndPathPos;
	READINTERVAL *edgeEndIntervals, *balEndIntervals;
	int *edgeEndPathGap, *balEndPathGap;
	edgeEndPathIndices = edgeEndPathPos = balEndPathIndices	= balEndPathPos = NULL;
	edgeEndIntervals = balEndIntervals = NULL;

	int numEdgeEndPaths, numBalEndPaths;

	numEdgeEndPaths = ShortXCut(edge->begin, edge, 
															path, numPaths,
															&edgeEndPathIndices, &edgeEndPathPos, &edgeEndPathGap,
															&edgeEndIntervals);


	if (edge != balEdge) {
		numBalEndPaths  = ShortXCut(balEdge->begin, balEdge,
																path, numPaths,
																&balEndPathIndices, &balEndPathPos, &balEndPathGap,
																&balEndIntervals);

		if (numBalEndPaths != numEdgeEndPaths) {
			printf("WARNING! When cutting end path intervals from an edge\n"
						 "and its balance, the number of cut intervals should be\n"
						 "the same.  Here it is %d (edge) and %d (bal)\n",
						 numEdgeEndPaths, numBalEndPaths);
		}
	}

	int numLastedge, numNextedge;

	EDGE *newEdge;
	newEdge = edge;
	numLastedge = -1;
	numNextedge = -1;
	int origNumVertices;
	int eqtransWorked = 0;
	while (edge != NULL and edge->end->num_nextedge > 0 and
				 (numLastedge != edge->end->num_lastedge or
					numNextedge != edge->end->num_nextedge)) {
		// The counts will change if an equivalent transformation
		// was successful.  If it wasn't, don't try any further
		// transformations.
		numLastedge = edge->end->num_lastedge;
		numNextedge = edge->end->num_nextedge;
		origNumVertices = numVertices;
		numVertices = run1vertex(edge->end, vertices, numVertices, 
														 path, numPaths, edge, numSeq, &newEdge);
		if (numVertices != origNumVertices){ 
			eqtransWorked = 1;
		}
		edge = newEdge;
	}		

	*posttransNumVertices = numVertices;
	if (eqtransWorked) {
		return 1;
	}
	
	else {
		// No equivalent transformation was possible.  It is possible
		// that other edges may be transformed using the intervals
		// that were deleted here, so we should put them back.

		// Put them back in the reverse order that they were removed 
		// in.  This way, if the balanced edge was part of the 
		// path of the original edge, it will be restored in the 
		// context that it was removed (with the interval for the 
		// original edge still removed).  Next restore the original
		// edge intervals, in the context of an existing balanced edge.
		if (edge != balEdge ) {
			RestoreXCut(balEdge, path, numPaths,
									balEndPathIndices, balEndPathPos, balEndPathGap,
									balEndIntervals, numBalEndPaths);
		}
		// 11:13
		RestoreXCut(edge, path, numPaths,
								edgeEndPathIndices, edgeEndPathPos, edgeEndPathGap,
								edgeEndIntervals, numEdgeEndPaths);
		return 0;
	}
}


int TryXCuts(NODES **vertex, int numVertices, PATH *path, int numPaths, int numSeq, 
						 int maxEdgeLength,
						 int *transformedNumVertices) {

	int foundHelpingCut = 0;
	int v, e;

	// Find the maximum number of nextedges
	int *tmpNextedge;
	int eqtransWorked = 1;
	int postTransNumVertices;	
	do {
		foundHelpingCut = 0;
		for (v = 0; v < numVertices; v++) {
			// Copy the nextedge so that no edge gets reprocessed.
			eqtransWorked = 1;
			while (eqtransWorked) {
				eqtransWorked = 0;
				while (e < vertex[v]->num_nextedge and eqtransWorked == 0) {
					if (vertex[v]->nextedge[e]->length < maxEdgeLength) {
						if (ShortEdgeCutTransform(vertex, numVertices, 
																			vertex[v]->nextedge[e],
																			vertex[v]->nextedge[e]->bal_edge,
																			path, numPaths, &postTransNumVertices, numSeq)) {
							eqtransWorked = 1;
							foundHelpingCut = 1;
							numVertices = postTransNumVertices;
						}
					}
					e++;
				}
			}
		}
	} while (foundHelpingCut == 1);
	if (foundHelpingCut) {
		*transformedNumVertices = postTransNumVertices;
		return 1;
	}
	else {
		*transformedNumVertices = numVertices;
		return 0;
	}
}

int vertex_run(NODES **vertex, int num_vertex, PATH *path, int num_path, int small_cover, int num_seq)
{
	int	i, j, k, l, m, n, q;
	int	tot_edge, tot_edge_old;
	int	**num_pa;
	int	n1, n2;
	NODES	*begin, **start_vertex;
	EDGE	*edge, *edge1, *edge2;
	/* some debugging code
	 */

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}

	start_vertex = (NODES **) ckalloc(num_vertex * sizeof(NODES *));
	do	{
		tot_edge_old = tot_edge;
		n = 0;
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_lastedge == 0)	{
				start_vertex[n ++] = vertex[i];
			}
		}
		/*
		  Process source vertices first.
		*/
		EDGE *newedge;
		for(i = 0; i < n; i ++)	{
			for(j = 0; j < start_vertex[i] -> num_nextedge; j ++)	{
				edge = start_vertex[i] -> nextedge[j];
				begin = edge -> end;
				/* don't do anything to this if it is circular*/
				if(begin == start_vertex[i])	continue;
				/* 
				   Attempt to do equivalent transformation on this vertex 
				   until there are no paths going through the vertex 
				   or the vertex has not been changed by transformations.
				*/
				while(begin -> num_nextedge > 0 and edge != NULL)	{
					n1 = begin -> num_lastedge;
					n2 = begin -> num_nextedge;
					/*
						printf("runvertex, running source: %d %d\n", begin->index, edge->index);
					*/
					num_vertex = run1vertex(begin, vertex, num_vertex, path, num_path, edge, num_seq, &newedge);
					if(n1 == begin -> num_lastedge || n2 == begin -> num_nextedge)	break;
					edge = newedge;
				}
			}
		}
		/* 
		   Now transform all other vertices.
		*/
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_lastedge == 0 || vertex[i] -> num_nextedge == 0)	continue;
			for(j = 0; j < vertex[i] -> num_lastedge; j ++)	{
			  /* try to transform this vertex by considering the paths that run through it.*/
				num_vertex = run1vertex(vertex[i], vertex, num_vertex, path, num_path,
																vertex[i] -> lastedge[j], num_seq, &newedge);
			}
		}
		/* 
		   Now clean up the graph.  If there are vertices that have been separated from the
		   graph, remove them.  Also remove parallel edges.
		*/
		num_vertex = merge_graph_path(vertex, num_vertex, path, num_path, small_cover);
		
		tot_edge = count_edge_simp(vertex, num_vertex, num_pa);
	} while(tot_edge_old != tot_edge);

	for(i = 0; i < MAX_BRA; i ++)	{
		free((void *) num_pa[i]);
	}
	free((void **) num_pa);
	free((void **) start_vertex);

	return(num_vertex);
}

int run1vertex(NODES *curVertex, NODES **vertex, int num_vertex, PATH *path, int num_path, EDGE *inedge, int num_seq,
							 EDGE **ranEdge)
{
  /* Input:
     The graph structure is:
       -- inedge --> curVertex 

       vertex :all vertices
       path   :all paths
  */


	int	i, j, k, l, m, n, n1, n2, c, q, p, t, f_edge[2];
	int freedDetachedEdge[2];
	int	num_nextedge, num_lastedge;
	int	*match;
	EDGE	*edge, *detachedEdge, *detachedBalEdge, *lastedge, *nextedge, *balEdgeNext, *balEdge;
	NODES	*vertex0;

	NODES *curVertexCopy;
	EDGE  *edgeCopy;
	EDGE  *nextEdgeCopy;
	/*	printf("running vertex %d -> %d\n", inedge->index, curVertex->index);*/
	if(curVertex -> num_lastedge == 0 || curVertex -> num_nextedge == 0)	return(num_vertex);
	l = rm_edge(curVertex, path, num_path);

	n = chkcross(curVertex, path, num_path);
	if(n > 0)	{
		return(num_vertex);
	}

	match = (int *) ckalloc(MAX_BRA * sizeof(int));

	num_lastedge = curVertex -> num_lastedge;
	num_nextedge = curVertex -> num_nextedge;
	/*
		printf("running vertex: %d edge: %d last: %d next %d\n", curVertex->index, inedge->index, num_lastedge, num_nextedge);
	*/
	int detachBalEdgeWorked = 0;
	if(num_lastedge == 1 && num_nextedge == 1)	{
	  /* 
	     This is a simple path. Remove 'curVertex' from the graph
	     by merging all paths that were in nextedge into lastedge.
	     
	     There are some caveats:
	       - If this vertex is the source 
	         of it's own balanced edge, there is no need to 
		 do merging on the balanced edge.
	       - If this vertex is a self-loop, no need to do merging 
	         on the balanced edge.
	  */
		lastedge = curVertex -> lastedge[0];
		vertex0 = lastedge -> bal_edge -> begin;
		nextedge = curVertex -> nextedge[0];
		balEdgeNext = lastedge -> bal_edge;
		balEdge = nextedge -> bal_edge;
		detachedEdge = merge_vertex_path(curVertex, path, num_path);
		if(!detachedEdge)	{
			return(num_vertex);
		}
		if(vertex0 == curVertex)	{
			detachedEdge -> bal_edge = detachedEdge;
		} else if(vertex0 -> lastedge[0] != vertex0 -> nextedge[0])	{
			detachedBalEdge = merge_vertex_path(vertex0, path, num_path);
			/*
			  If previous balance is the same as the previous edge
			  or the next balance is the same as the next edge
			  the merged balanced edge is the same as the balanced edge.
			*/
			if(balEdgeNext == lastedge || balEdge == nextedge)	{
				detachedBalEdge -> bal_edge = detachedBalEdge;
			} else	{
			  /* 
			     Otherwise, the new balanced edge may be assigned as normal.
			  */
				detachedEdge -> bal_edge = detachedBalEdge;
				detachedBalEdge -> bal_edge = detachedEdge;
			}
		} else	{
			detachedEdge -> bal_edge = detachedEdge;
		}
	} else	{
	  /* 
	     This vertex is degree(in) > 1 or degree(out) > 1.  First find out how many paths
	   */
		edge = inedge;
		//		assert(edge->deleted == 0);
		/* 
		   Find out how many paths leave 'curVertex' that came from 'edge'.
		 */
		n = searchnext(curVertex, edge, path, num_path, match);
		for(k = 0; k < n; k ++)	{
		  /* Attempt to do equivalent transformations on 'vertex' for 
		     for paths  'edge'->curVertex->'nextedge'
		  */
			nextedge = curVertex -> nextedge[match[k]];
			//			assert(nextedge->deleted == 0);
			m = chkforbid(edge, nextedge, path);
			/* 
			   Equivalent transformation is ok only if there are edges that pass through 
			   edge*,edge,nextedge,edge** where edge* and edge** are other edges in the graph.
			   Alternatively, if no such path exists, if all paths that start in 'edge' continue
			   to the same next edge, or all paths that end in 'nextedge' originate from the
			   same edge, it is ok to do transformation.
			*/
			if(m)	{
				/*
					printf("edge: %d, %d contains forbidden paths\n", edge, nextedge);
				*/
				continue;
			}

			/* 
			   It is ok to do transformation.  We'll want to do the transformation both
			   on the edge, nad it's balanced edges.
			*/
			balEdgeNext = nextedge -> bal_edge;
			balEdge = edge -> bal_edge;
			/*
			  printf("detaching edge (%d) and nextedge (%d)\n", edge->multip, nextedge->multip);
			*/
			curVertexCopy = NULL;
			edgeCopy      = NULL;
			nextEdgeCopy  = NULL;

			CopyVertex(curVertex, &curVertexCopy);
			CopyEdge(edge, &edgeCopy);
			CopyEdge(nextedge, &nextEdgeCopy);

			/*
				struct DetachedConfiguration detachedConfig;
				detachedConfig.srcNumPaths = edge->begin->num_path;
				detachedConfig.destNumPaths = edge->end->num_path;
			*/
			detachedEdge = detach_bal(edge, nextedge, path, num_path, f_edge, num_seq);

			if (detachedEdge && detachedEdge->begin && detachedEdge->end)
			detachedEdge->deleted = 0;

			/* 
			   If the detachment was successful -- it isn't if the edges are not adjacent.  Maybe this should 
			   be changed to an assert?
			*/
			if(detachedEdge)	{
			  /* 
			     If the edge is its' own complement, don't detach the balanced edge.
			     The correct assignment of the balanced edge is the edge itself.
			  */
				if(balEdge == nextedge && balEdgeNext == edge)	{
					detachedEdge -> bal_edge = detachedEdge;
				} else	{
				  /*
				    Standard case: the two balanced edges are distinct from the edges that were detached.
				  */
					if(balEdgeNext != nextedge && balEdge != edge)	{
					  /* 
					     
					   */
						vertex0 = balEdgeNext -> end;
						if(vertex0 == curVertex && 
						   vertex0 -> num_lastedge == 1 &&
						   vertex0 -> num_nextedge == 1)	{
						  /*
						    The edge that was detached was routed through the end of the balanced edge.
						    If that moved the degree from 4 down to 2, then the vertex does not need to 
						    be detached to be balanced.
						    This really doesn't need to be here, since it's possible that the edges 
						    may be detached in the first part of this.
						   */
							detachedBalEdge = merge_vertex_path(vertex0, path, num_path);
						} else	{
							detachedBalEdge = detach_bal(balEdgeNext, balEdge, path, num_path, 
																					 freedDetachedEdge, num_seq);
						}
						if (detachedBalEdge) {
							detachBalEdgeWorked = 1;
							detachedBalEdge->deleted = 0;
							detachedEdge -> bal_edge = detachedBalEdge;
							detachedBalEdge -> bal_edge = detachedEdge;
						}
					} else if(balEdgeNext == nextedge)	{
					  /* Strange case 1: 
					     'nextedge' is palindromic, so it's balanced edge is itself
					  */
						if(f_edge[1])	{
							vertex0 = detachedEdge -> end;
							if(vertex0 -> num_lastedge == 1 &&
							   vertex0 -> num_nextedge == 1)	{
								detachedBalEdge = merge_vertex_path(vertex0, path, num_path);
							} else	{
								detachedBalEdge = detach_bal(detachedEdge, balEdge, path, num_path, 
																						 freedDetachedEdge, num_seq);
							}
							if (detachedBalEdge) 
								detachBalEdgeWorked = 1;

							if (detachedBalEdge) {
								detachedBalEdge -> bal_edge = detachedBalEdge;
								detachedBalEdge->deleted = 0;
							}
						} else	{
							n1 = countmatch(detachedEdge, balEdge, path, num_path);
							n2 = countmatch(nextedge, balEdge, path, num_path);
							if(n1 > 0)	{
							   detachedBalEdge = detach_bal(detachedEdge, balEdge, path, num_path, freedDetachedEdge, num_seq);
								 if (detachedBalEdge) {
									 detachedBalEdge -> bal_edge = detachedBalEdge;
									 detachedBalEdge->deleted = 0;
								 }

								if (detachedBalEdge) 
									detachBalEdgeWorked = 1;
							}
							if(n2 > 0)	{
								detachedBalEdge = detach_bal(nextedge, balEdge, path, num_path, freedDetachedEdge, num_seq);
								if (detachedBalEdge) {
									detachedEdge -> bal_edge = detachedBalEdge;
									detachedBalEdge -> bal_edge = detachedEdge;
									detachBalEdgeWorked = 1;
								}
							}
						}
					} /* End dealing with balEdge == nextedge */
					else if(balEdge == edge)	{
						if(f_edge[0])	{
							vertex0 = balEdgeNext -> end;
							if(vertex0 -> num_lastedge == 1 &&
							   vertex0 -> num_nextedge == 1)	{
								detachedBalEdge = merge_vertex_path(vertex0, path, num_path);
							} else	{
								detachedBalEdge = detach_bal(balEdgeNext, detachedEdge, path, num_path, freedDetachedEdge, num_seq);
							}
							if (detachedBalEdge) {
								detachedBalEdge -> bal_edge = detachedBalEdge;
								detachBalEdgeWorked = 1;
							}
						} else	{
							n1 = countmatch(balEdgeNext, detachedEdge, path, num_path);
							n2 = countmatch(balEdgeNext, edge, path, num_path);
							if(n1 > 0)	{
							   detachedBalEdge = detach_bal(balEdgeNext, detachedEdge, path, num_path, freedDetachedEdge, num_seq);
								 if (detachedBalEdge) {
									 detachedBalEdge -> bal_edge = detachedBalEdge;
									 detachBalEdgeWorked = 1;
								 }
							}
							if(n2 > 0)	{
								detachedBalEdge = detach_bal(balEdgeNext, edge, path, num_path, freedDetachedEdge, num_seq);
								if (detachedBalEdge and detachedEdge) {
									detachedEdge -> bal_edge = detachedBalEdge;
									detachedBalEdge -> bal_edge = detachedEdge;
									detachBalEdgeWorked = 1;
								}
							}
						}
					}
				} /* end dealing with balEdge == edge */
				if (!detachBalEdgeWorked) {
					printf("detaching balanced edge FAILED!!!\n");

					printf("trying to put back edges: src: %d - %d , %d , %d -> %d \n",
								 detachedEdge->begin->index, edgeCopy->index, curVertex->index, nextEdgeCopy->index,
								 detachedEdge->end->index);
					UndetachEdges(detachedEdge->begin, detachedEdge->end, detachedEdge,
												path,	f_edge,
												edgeCopy, curVertexCopy, nextEdgeCopy,
												edge, curVertex, nextedge);

					/*
						printf("the paths after restoring are: \n");
						PrintPaths(curVertex, path);
					*/
				}
				break;
			} // end if (detachedEdge) ... if the first detach bal worked
		}
	}
	if (detachBalEdgeWorked) {
		if (edge != balEdge and nextedge != balEdgeNext)
			*ranEdge = detachedEdge;
		else
			*ranEdge = detachedBalEdge;

		if (curVertexCopy != NULL)
			FreeVertex(curVertexCopy);
		if (edgeCopy != NULL)
			FreeEdge(edgeCopy);
		if (nextEdgeCopy != NULL)
			FreeEdge(nextEdgeCopy);
	}
	else {
		*ranEdge = NULL;
	}

	free((void *) match);
	return(num_vertex);
}

int createvertex(NODES ***vertexListPtr, int num_vertex, NODES *curVertex, PATH *path, int num_path)
{
	printf("THIS CODE DOES NOTHING\n");
	assert(0);
	int	i, j, k, l, m, n, q;
	EDGE	*edge, *edge0;
	NODES **vertex = *vertexListPtr;
	for(i = 0; i < curVertex -> num_nextedge; i ++)	{
		k = crossedge(curVertex -> nextedge[i], path);
		if(k == 0)	{
			edge = curVertex -> nextedge[i];
			edge0 = edge -> bal_edge;
			num_vertex = newvertex(vertexListPtr, num_vertex, edge, path);
			if(edge != edge0)	{
				num_vertex = newvertex(vertexListPtr, num_vertex, edge0, path);
			}
		}
	}
	for(i = 0; i < curVertex -> num_lastedge; i ++)	{
		k = crossedge(curVertex -> lastedge[i], path);
		if(k == 0)	{
			edge = curVertex -> lastedge[i];
			edge0 = edge -> bal_edge;
			num_vertex = newvertex(vertexListPtr, num_vertex, edge, path);
			if(edge != edge0)	{
				num_vertex = newvertex(vertexListPtr, num_vertex, edge0, path);
			}
		}
	}

	return(num_vertex);
}

int chkcross(NODES *curVertex, PATH *path, int num_path)
{
  /* 
     Input: a vertex, and all paths in the graph
     Output: the number of in-edges and out-edges that do not have any paths going through them.
  */

	int	i, j, k, l, n, n1, n2;
	int	**match, *match1, *match2;
	int *saneLastedge, *saneNextedge;
	/* 
		 I have no idea what this next statment is for. 
	*/
	if((curVertex -> num_lastedge > 1  || curVertex -> lastedge[0] -> length < 500) &&
		 (curVertex -> num_nextedge > 1 || curVertex -> nextedge[0] -> length < 500))	return(0);

	match1 = (int *) ckalloc(curVertex -> num_lastedge * sizeof(int));
	match2 = (int *) ckalloc(curVertex -> num_nextedge * sizeof(int));
	saneLastedge  = (int *) ckalloc(curVertex -> num_lastedge * sizeof(int));
	saneNextedge  = (int *) ckalloc(curVertex -> num_nextedge * sizeof(int));

	for (i = 0; i < curVertex->num_lastedge; i++) {match1[i] = 0; saneLastedge[i] = 0;};
	for (i = 0; i < curVertex->num_nextedge; i++) {match2[i] = 0; saneNextedge[i] = 0;};

	match = (int **) ckalloc(curVertex -> num_lastedge * sizeof(int *));
	for(j = 0; j < curVertex -> num_lastedge; j ++)	{
		match[j] = (int *) ckalloc(curVertex -> num_nextedge * sizeof(int));
	}

	/* Count the routing of paths through this vertex. For example, how many
	   edges enter through the first in edge and exit through the first out edge,
	   enter the first in edge and out the second, and so on.
	*/
	   
	for(i = 0; i < curVertex -> num_path; i ++)	{
		j = curVertex -> path_index[i];
		k = curVertex -> path_pos[i];
		if(k > 0 && k < path[j].len_path)	{
		  // get the 'lastedge' index where path[j] enters this vertex
			n1 = searcherase(curVertex -> lastedge, path[j].edge[k - 1], curVertex -> num_lastedge);
			//  get the 'nextedge' index where path[j] leaves this vertex
			n2 = searcherase(curVertex -> nextedge, path[j].edge[k], curVertex -> num_nextedge);
			// record that
			match[n1][n2] ++;
			saneLastedge[n1]++;
			saneNextedge[n2]++;
		}
	}

	for(j = 0; j < curVertex -> num_lastedge; j ++)	{
		for(i = 0; i < curVertex -> num_nextedge; i ++)	{
		  // match1 stores the number of paths that have entered the vertex through [j]
			match1[j] += match[j][i];
			// match2 stores the number of paths that leave the vertex through i
			match2[i] += match[j][i];
		}
	}
	
	// count the number of 'blank' edges, edges of zero multiplicity going into
	// vertex 'curVertex'.

	/* 
	   This checks the marginal sums of the paths entering each edge. I have no idea why this isn't 
	   simply summing the number of paths that enter each edge, and I might add a sanity
	   check just to show that's the case.
	*/
	n = 0;
	for(i = 0; i < curVertex -> num_lastedge; i ++)	{
		if(match1[i] == 0)	{
		  assert(saneLastedge[i]==0);
			n ++;
		}
	}
	for(i = 0; i < curVertex -> num_nextedge; i ++)	{
		if(match2[i] == 0)	{
		  assert(saneNextedge[i] ==0);
			n ++;
		}
	}

	for(j = 0; j < curVertex -> num_lastedge; j ++)	{
		free((void *) match[j]);
	}
	free((void **) match);
	free((void *) match1);
	free((void *) match2);

	free((void*) saneLastedge);
	free((void*) saneNextedge);

	return(n);
}

int chkforbid(EDGE *edge1, EDGE *edge2, PATH *path)
{
	int	i, j, k, l, n1, n2;
	NODES	*vertex;
	EDGE	*edge;
	/*
	  Input:
	  Part of a graph of the format: --edge1-->vertex--edge2--> .  So edge1->end must be 
	  edge2->begin.
	  
	  It is ok to do equivalent transformation on edge1,edge2 if there is a path that
	  goes through edge1,edge2 and continues on: (edge1,edge2,edge*), where edge* is any 
	  other edge. 
	*/
	  
	if(edge1 -> end != edge2 -> begin)	{
		printf("Edges are not linked!. %d(%d-%d) %d(%d-%d)\n", edge1, edge1 -> begin,
			edge1 -> end, edge2, edge2 -> begin, edge2 -> end);
		exit(0);
	}
	if(edge1 -> begin == edge1 -> end || edge2 -> begin == edge2 -> end)	return(0);
	if(edge1 -> length > SMALL_EDGE)	{
		vertex = edge1 -> end;
		edge = (EDGE *) NULL;
		n1 = 0;
		for(i = 0; i < vertex -> num_path; i ++)	{
		  /* 
		     Don't consider the path through this vertex if :
		       - it starts at this vertex
		       - the previous edge isn't the one we're looking for, 'edge1'
		  */
			if(vertex -> path_pos[i] == 0 ||
			   path[vertex -> path_index[i]].edge[vertex -> path_pos[i] - 1] != edge1)	{
				continue;
			}
			/* 
			   If the path going through this vertex ('i') is not at the end of the path
			   and the index of the vertex along the path is not the first, the first edge is ok.
			*/
			if(vertex -> path_pos[i] < path[vertex -> path_index[i]].len_path
			   && vertex -> path_pos[i] > 1)	{
				n1 = 0;
				break;
			} 
			/* If this is the second vertex on this path, and there are more than 2 vertices on this 
			   path:
			*/
			else if(vertex -> path_pos[i] == 1 && 
							vertex -> path_pos[i] < path[vertex -> path_index[i]].len_path)	{
			  /*
			    The vertex is on the first edge of path 'i'. We would like to make sure that
			    all paths that have vertex as the dest of the first edge continue to the same edge.
			  */
				if(!edge)	{
				  /* 
				     If this is the first path we've seen that has 'vertex' as the first dest, save 
				     the edge.
				  */
				  edge = path[vertex -> path_index[i]].edge[vertex -> path_pos[i]];
				} else if(edge != path[vertex -> path_index[i]].edge[vertex -> path_pos[i]])	{
				  /*
				    This is not the first path that we've seen that's reached this
				    vertex on the first edge. We want to make sure all the paths that go 
				    through this vertex come from the same edge, but this one doesn't.  Record that
				    by setting n1 to 1 (retval > 0 means forbidden)
				  */
					n1 = 1;
				}
			}
		}
	} else {
		//		printf("edge1length: %d\n", edge1->length);
		n1 = 0;
	}
	if(edge2 -> length > SMALL_EDGE)	{
		vertex = edge2 -> begin;
		edge = (EDGE *) NULL;
		n2 = 0;
		for(i = 0; i < vertex -> num_path; i ++)	{
			/* 
			   Don't consider paths that are not ending after this vertex.
			*/
			if(vertex -> path_pos[i] == path[vertex -> path_index[i]].len_path ||
			   path[vertex -> path_index[i]].edge[vertex -> path_pos[i]] != edge2)	{
				continue;
			}
			/* 
			   There is an edge that passes through this vertex, and continues on. 
			   It is not forbidden.
			 */
			if(vertex -> path_pos[i] < path[vertex -> path_index[i]].len_path - 1
			   && vertex -> path_pos[i] > 0)	{
				n2 = 0;
				break;
			} 
			/* 
			   A path ends at the edge leaving this vertex. If all paths end after
			   this vertex, they need have came from the same edge.  If they came from
			   different edges, equivalent transformation on this edge is forbidden.
			*/
			else if(vertex -> path_pos[i] > 0 && 
			   vertex -> path_pos[i] == path[vertex -> path_index[i]].len_path - 1)	{
				if(!edge)	{
					edge = path[vertex -> path_index[i]].edge[vertex -> path_pos[i] - 1];
				} else if(edge != path[vertex -> path_index[i]].edge[vertex -> path_pos[i] - 1]) {
					n2 = 1;
				}
			}
		}
	} else {
		n2 = 0;
	}

	return(n1 || n2);
}
