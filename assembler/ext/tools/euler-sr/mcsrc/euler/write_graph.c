/***************************************************************************
 * Title:          write_graph.c
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

void write_graph(NODES **vertex, int num_vertex, FILE *fp, FILE *fp1);

void write_paths(NODES **vertex, int numVertices, 
								 PATH *path, int numPaths, FILE *fp) {
	fprintf(fp, "%d\n", numPaths);

	int v, e, i;

	int curReadPos;
	int read, ri;

	for (read = 0; read < numPaths; read++) { 
		for (e = 0; e < path[read].len_path; e++) {
			for (ri = 0; ri < path[read].edge[e]->multip; ri++) {
				path[read].edge[e]->readinterval[ri] = path[read].edge[e]->readinterval[ri];
			} 
		}
	}

	
	// Sort the read intervals on edges in order of their
	// index, then by the order of the position in the read.
	for (v =0 ; v < numVertices; v++ ) {
		for (e = 0; e < vertex[v]->num_nextedge; e++) {
			// record current order before sorting
			for (i = 0; i < vertex[v]->nextedge[e]->multip; i++) {
				vertex[v]->nextedge[e]->readinterval[i].pos = i;
			}
			sortreadinterval_index(vertex[v]->nextedge[e]->readinterval,
														 vertex[v]->nextedge[e]->multip);
			//		  printf("ri: %d\n", path[10].edge[2]->readinterval);
		}
	}
	
	int pathStart;
	int readIndex;

	for (read = 0; read < numPaths; read++) { 
		for (e = 0; e < path[read].len_path; e++) {
			for (ri = 0; ri < path[read].edge[e]->multip; ri++) {
				path[read].edge[e]->readinterval[ri] = path[read].edge[e]->readinterval[ri];
			} 
		}
	}


	for (read = 0; read < numPaths; read++) { 
		fprintf(fp, "%d %d ", path[read].readindex-1, path[read].len_path);
		curReadPos = -1;
		// for each edge in a path, find where 
		// the corresponding read interval is stored
		// for the path.
		
		for (e = 0; e < path[read].len_path; e++) {
			// Locate the read in the list of read intervals
			// sorted by read index.
			readIndex = path[read].readindex-1;
			for (ri = 0; ri < path[read].edge[e]->multip; ri++) {
				path[read].edge[e]->readinterval[ri] = path[read].edge[e]->readinterval[ri];
			} 
			if (FindRead(path[read].edge[e]->readinterval, path[read].edge[e]->multip,
									 readIndex, &pathStart)) {
				// The read intervals are re-sorted according to the position
				// in the read. If one edge has multiple intervals of the same
				// read, the position of the correct read interval is the next read
				// interval that starts after curReadPos.
				// The FindRead function always references the beginning of the 
				// run of read intervals from the same read.
				while (pathStart < path[read].edge[e]->multip and 
							 path[read].edge[e]->readinterval[pathStart].eq_read == readIndex and
							 path[read].edge[e]->readinterval[pathStart].begin < curReadPos) {
					pathStart++;
				}
				if (path[read].edge[e]->readinterval[pathStart].eq_read == readIndex ) {
					fprintf(fp, " %d %d", path[read].edge[e]->index, 
									path[read].edge[e]->readinterval[pathStart].pos);
					curReadPos = path[read].edge[e]->readinterval[pathStart].begin;
				}
			}
			else {
				/*
					printf("ERROR! A path references an edge that does not\n" 
					"have a read interval for that path.");
				*/
			}
			for (ri = 0; ri < path[read].edge[e]->multip; ri++) {
				path[read].edge[e]->readinterval[ri] = path[read].edge[e]->readinterval[ri];
			} 
		}
		fprintf(fp, "\n");
	}

	// Put the read intervals back in the original order
	for (v =0 ; v < numVertices; v++ ) {
		for (e = 0; e < vertex[v]->num_nextedge; e++) {
			sortreadinterval(vertex[v]->nextedge[e]->readinterval,
											 vertex[v]->nextedge[e]->multip);
		}
	}
}

void write_branching_graph(NODES **vertex, int num_vertex, FILE *fp) {
	fprintf(fp, "Number_of_vertices %d\n", num_vertex);

	int numEdges = 0;
	int v, e;
	int edgeIndex = 0;
	for (v = 0; v < num_vertex; v++) {
		vertex[v]->index = v;
		numEdges += vertex[v]->num_nextedge;
		for (e = 0; e < vertex[v]->num_nextedge; e++) {
			vertex[v]->nextedge[e]->index = edgeIndex;
			edgeIndex++;
		}
	}
	
	EDGE **edgeList = (EDGE**) ckalloc(sizeof(EDGE*) * numEdges);
	edgeIndex = 0;
	for (v = 0; v < num_vertex; v++) {
		for (e = 0; e < vertex[v]->num_nextedge; e++) {
			edgeList[edgeIndex] = vertex[v]->nextedge[e];
			edgeIndex++;
		}
	}
	int edgeOut[4];
	int edgeIn[4];


	for (v = 0; v < num_vertex; v++ ) {
		fprintf(fp, "%d %d ", v, vertex[v]->index);
		fprintf(fp, "out %d", vertex[v]->num_nextedge);
		for (e = 0; e < vertex[v]->num_nextedge; e++) {
			edgeIndex = vertex[v]->nextedge[e]->index;
			fprintf(fp, " %d", edgeIndex);
		}
		fprintf(fp, " in %d", vertex[v]->num_lastedge);
		for (e = 0 ; e < vertex[v]->num_lastedge; e++) {
				edgeIndex = vertex[v]->lastedge[e]->index;
				fprintf(fp, " %d", edgeIndex);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "Number_of_edges %d\n", numEdges);	
	for (e = 0; e < numEdges; e++) {
		fprintf(fp, "%d %d %d %d %d %d %d\n", e, 
						edgeList[e]->begin->index, edgeList[e]->end->index, edgeList[e]->multip,
						edgeList[e]->index, edgeList[e]->bal_edge->index,
						edgeList[e]->length);
	}
}

void write_graph(NODES **vertex, int num_vertex, FILE *fp, FILE *fp1)
{
	int	i, j, k, l, n;
	char	temp[100];
	EDGE	*edge;

	n = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			sprintf(temp, "edge%d", n + 1);
			writeseq(fp, edge -> seq, temp, edge -> length);
			edge -> start_cover = n;
			n ++;
		}
	}

	fprintf(fp1, "Number_of_Vertex %d\n", num_vertex);
	for(i = 0; i < num_vertex; i ++)	{
		fprintf(fp1, "Vertex %d %d %d\n", i, vertex[i] -> num_nextedge, vertex[i] -> num_lastedge);
		fprintf(fp1, "Last_edge ");
		for(j = 0; j < vertex[i] -> num_lastedge; j ++)	{
			fprintf(fp1, "%d %d ", vertex[i] -> lastedge[j] -> start_cover, vertex[i] -> lastedge[j] -> bal_edge -> start_cover);
		}
		fprintf(fp1, "\n");
		fprintf(fp1, "Next_edge ");
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			fprintf(fp1, "%d %d ", vertex[i] -> nextedge[j] -> start_cover, vertex[i] -> nextedge[j] -> bal_edge -> start_cover);
		}
		fprintf(fp1, "\n");
	}
}
