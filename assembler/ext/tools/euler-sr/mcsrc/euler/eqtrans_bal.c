/***************************************************************************
 * Title:          eqtrans_bal.c
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
extern int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;


int eqtrans_bal(NODES ***vertexPtr, int num_vertex, PATH *path, int num_path, int num_seq, int cut);

int eqtrans_bal(NODES ***vertexPtr, int num_vertex, PATH *path, int num_path, int num_seq, int cut)
{
	int	i, j, k, l, n, m;
	int	tot_edge, tot_edge_old, tot_path, tot_path_old;
	int	**num_pa;
	int	votethresh;
	NODES **vertex = *vertexPtr;
	EDGE	*edge, **all_edge, *edge1, *edge2, *edge0;
	int v, e;
	votethresh = LOW_COV;

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}

	tot_edge = count_edge_simp(vertex, num_vertex, num_pa);
	int totalLength, totalReadStarts;
	do	{
		tot_edge_old = tot_edge;
		tot_path_old = tot_path;
		if (cut) {
			RemoveUninformativePaths(path, num_path);
		}
		num_vertex = vertex_run(vertex, num_vertex, path, num_path, votethresh, num_seq);
		CountCoverage(vertex, num_vertex, &totalLength, &totalReadStarts);
		if (!cut) {
			/* We only remove low coverage if the uninformative paths are not cut.*/

			num_vertex = RemoveLowCoverageEdges(vertex, num_vertex, path, num_path,
																					totalLength, totalReadStarts,
																				3.5, 5);
			
		}

		tot_edge = count_edge_simp(vertex, num_vertex, num_pa);
		tot_path = count_path(path, num_path);
	} while(tot_edge_old > tot_edge || tot_path_old > tot_path);

	/*
		Some debugging information that I'll have to clear up soon.
	for (v= 0; v < num_vertex; v++) {
		for (e= 0 ; e < vertex[v]->num_nextedge; e++) {
			if (vertex[v]->nextedge[e]->multip !=
					vertex[v]->nextedge[e]->bal_edge->multip) {
								printf("in eqtrans_bal, edge: %d, multip: %d %d bal: %d multip: %d\n", vertex[v]->nextedge[e]->index,
							 vertex[v]->nextedge[e]->multip,
							 vertex[v]->nextedge[e]->bal_edge->index,
							 vertex[v]->nextedge[e]->bal_edge->multip,
							 vertex[v]->nextedge[e]->bal_edge->length);
			}
		}
	}
	*/

	num_vertex = splitbeg(vertexPtr, num_vertex, path, num_path);
	vertex = *vertexPtr;
	num_vertex = merge_graph_path(vertex, num_vertex, path, num_path, votethresh);
	tot_edge = count_edge_simp(vertex, num_vertex, num_pa);

	for(i = 0; i < MAX_BRA; i ++)	{
		free((void *) num_pa[i]);
	}
	free((void **) num_pa);
	return(num_vertex);
}
