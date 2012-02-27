/***************************************************************************
 * Title:          branch_graph.c
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

int branch_graph(int num_seq, int *len_seq, NODES **nodes, EDGE **edge, int num_nodes);

int branch_graph(int num_seq, int *len_seq, NODES **nodes, EDGE **edge, int num_nodes)
{
	int	i, j, k, l, m, n;
	int	num_edge;
	EDGE	*edge1;

	num_edge = 0;
	for(i = 0; i < num_nodes; i ++)	{
		for(j = 0; j < nodes[i] -> num_nextedge; j ++)	{
			edge1 = nodes[i] -> nextedge[j];
			num_edge = makedge_new(nodes[i], edge, num_edge, edge1, nodes);
			nodes[i] -> nextedge[j] = edge[num_edge - 1];
		}
	}

/*	Remove edges that are not covered	*/

	i = 0;
	while(i < num_edge)	{
		if(edge[i] -> multip == 0)	{
			erasedge(edge[i]);
			edge[i] = edge[num_edge - 1];
			num_edge --;
		} else	{
			i ++;
		}
	}

/*	Set up balancing node	*/
	for(i = 0; i < num_nodes; i ++)	{
		n = nodes[i] -> bal_node -> nlinks;
		nodes[i] -> bal_node = nodes[n];
		if(nodes[i] -> bal_node -> num_nextedge != nodes[i] -> num_lastedge ||
		   nodes[i] -> bal_node -> num_lastedge != nodes[i] -> num_nextedge)	{
			printf("nlinks %d %d\n", i, n);
			printf("i %d nodes %d %d %d bal_nodes %d %d %d\n", i, nodes[i],
				nodes[i] -> num_lastedge, nodes[i] -> num_nextedge,
				nodes[i] -> bal_node, nodes[i] -> bal_node -> num_lastedge,
				nodes[i] -> bal_node -> num_nextedge);
			getchar();
		}
	}
	printf("Edge made.\n");
	return(num_edge);
}
