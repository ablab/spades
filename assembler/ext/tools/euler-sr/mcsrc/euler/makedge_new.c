/***************************************************************************
 * Title:          makedge_new.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <branch_gviz.h>
#include <clean_graph.h>
#include <extvab.h>
#include <extfunc.h>

int makedge_new(NODES *node, EDGE **edge, int num_edge, EDGE *edge1, NODES **allnodes);
EDGE *newedge_new(EDGE **midedge, int num_midedge, NODES *begin, NODES *end);

int makedge_new(NODES *node, EDGE **edge, int num_edge, EDGE *edge1, NODES **allnodes)
{
	int	i, j, k, l;
	int	num_midedge;
	EDGE	**midedge, *nedge;
	NODES	*node1, *end;

	num_midedge = 0;
	midedge = (EDGE **) ckalloc(MAX_TMP_LEG * sizeof(EDGE *));
	midedge[num_midedge ++] = edge1;
	node1 = edge1 -> end;
	while((nedge = find_unique_oedge(node1)) != (EDGE *) NULL)	{
		midedge[num_midedge ++] = nedge;
		node1 = nedge -> end;
	}
	nedge = newedge_new(midedge, num_midedge, node, allnodes[node1 -> nlinks]);
	if(nedge)	{
		edge[num_edge ++] = nedge;
	} else	{
		printf("No returned edge\n");
		exit(0);
	}
	free((void **) midedge);
	return(num_edge);
}

EDGE *newedge_new(EDGE **midedge, int num_midedge, NODES *begin, NODES *end)
{
	int	i, j, k, l, n;
	INSERT	*insert;
	READINTERVAL	*readcov;
	NODES	*node, *node1;
	EDGE	*tmpedge;
	READPOSITION *readposition;

	readposition = (READPOSITION *) ckalloc(1 * sizeof(READPOSITION));
/*	Sort the readpositions every nodes cover by the index of the reads	*/

	for(i = 0; i < num_midedge; i ++)	{
		sortreadinterval_index(midedge[i] -> readinterval, midedge[i] -> multip);
	}

/*	Define read coverage of the edge	*/

	readcov = (READINTERVAL *) NULL;
	k = 0;
	for(i = 0; i < num_midedge; i ++)	{
		for(j = 0; j < midedge[i] -> multip; j ++)	{
			readcov = insert_readcov(readcov, midedge[i] -> readinterval[j].eq_read, midedge[i] -> readinterval[j].begin,
					 midedge[i] -> readinterval[j].length, k);
		}
		k += midedge[i] -> length - 1;
	}
	free((void *) readposition);

/*	Setup the beginning and ending vertices of new edge	*/

	tmpedge = (EDGE *) ckalloc(1 * sizeof(EDGE));
	tmpedge -> length = midedge[0] -> length;
	tmpedge -> begin = begin;
	tmpedge -> end = end;
	add_lastedge(tmpedge -> end, tmpedge);

/*	remove single read covers, copy and sort read covers	*/
	n = size_readinterval(readcov) + 1;
	tmpedge -> readinterval = (READINTERVAL *) ckalloc(n * sizeof(READINTERVAL));
	tmpedge -> multip = copyreadinterval(tmpedge -> readinterval, readcov);
	if(tmpedge -> multip == 0 || tmpedge -> multip != n - 1)	{
		printf("Warning: No read cover: tmpedge %d length %d.\n", tmpedge, num_midedge);
		readposition = tmpedge -> begin -> readposition;
		while(readposition)	{
			printf("node %d readposition %d %d\n", tmpedge -> begin,
				readposition -> readindex, readposition -> position);
			readposition = readposition -> next;
		}
		readposition = tmpedge -> end -> readposition;
		while(readposition)	{
			printf("node %d readposition %d %d\n", tmpedge -> end,
				readposition -> readindex, readposition -> position);
			readposition = readposition -> next;
		}
	}
	while(readcov)	{
		readcov = free_readinterval(readcov);
	}
	if(tmpedge -> multip > 0)	{
		sortreadinterval(tmpedge -> readinterval, tmpedge -> multip);
	}

/*	Set up the length of the edge	*/

	for(i = 1; i < num_midedge; i ++)	{
		tmpedge -> length += midedge[i] -> length - 1;
	}
	return(tmpedge);
}
