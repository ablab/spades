/***************************************************************************
 * Title:          graph.c
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

int graph(int num_seq, char **src_name, int *len_seq, READLIST **readlist, EDGE **edge);
void movereadinterval(EDGE *edge1, EDGE *edge);

int graph(int num_seq, char **src_name, int *len_seq, READLIST **readlist, EDGE **edge)
{
	int	i, j, k, l, m, n;
	char	c, c1;
	int	num_edge;
	int	pos1, pos2;
	READPOSITION	*readposition;
	INSERT	*insert, *insert1;
	EDGE	*edge1, *edge2;
	NODES	*node, *node_next, *node1, *node2, *node0;

	for(i = 0; i < 2 * num_seq; i ++)	{
		pos1 = 0;
		node = readlist[i][0].node;
		for(j = 1; j < len_seq[i]; j ++)	{
			pos2 = j;
			node_next = readlist[i][j].node;
			edge1 = insert_edge(node, node_next, i, pos1, pos2);
			node = node_next;
			pos1 = pos2;
		}
	}
	printf("Edges assigned.\n");

	for(i = 0; i < 2 * num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			node = readlist[i][j].node;
			readposition = node -> readposition;
			while(readposition)	{
				if(readposition -> position >= len_seq[readposition -> readindex] ||
				   readlist[readposition -> readindex][readposition -> position].node != node)	{
					printf("readposition not correct.\n");
					printf("readposition %d %d\n", readposition -> readindex, readposition -> position);
					exit(0);
				}
				readposition = readposition -> next;
			}
			if(node)	{
				node_next = node -> bal_node;
				if(node -> num_lastedge != node_next -> num_nextedge ||
				   node_next -> num_lastedge != node -> num_nextedge)	{
					printf("pos2 i %d j %d %d(%d-%d) %d(%d-%d)\n", i, j, node, node -> num_lastedge,
						node -> num_nextedge, node_next, node_next -> num_lastedge,
						node_next -> num_nextedge);
					printf("multip %d %d\n", node -> num_path, node_next -> num_path);
					getchar();
				}
			}
		}
	}

/*	Determine long edges and remove 1-in-1-out nodes	*/

	num_edge = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			node = readlist[i][j].node;
			if(node && node -> visit == 0)	{
				if(node -> num_lastedge == 1 && node -> num_nextedge == 1)	continue;
				num_edge = makedge(node, edge, num_edge, readlist);
			}
		}
	}
	printf("Edge made: %d.\n", num_edge);
	cleannode(readlist, len_seq, num_seq);

	for(i = 0; i < 2 * num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			node = readlist[i][j].node;
			if(node)	{
				node_next = node -> bal_node;
				if((node -> num_lastedge == 1 && node -> num_nextedge == 1) ||
				   (node -> num_lastedge != node_next -> num_nextedge ||
				   node_next -> num_lastedge != node -> num_nextedge))	{
					printf("pos3 i %d j %d %d(%d-%d) %d(%d-%d)\n", i, j, node, node -> num_lastedge,
						node -> num_nextedge, node_next, node_next -> num_lastedge,
						node_next -> num_nextedge);
					printf("multip %d %d\n", node -> npos, node_next -> npos);
					getchar();
				}
			}
		}
	}
	return(num_edge);
}

void movereadinterval(EDGE *edge1, EDGE *edge)
{
	int	i, j, k, l, m;
	READINTERVAL	*tmpreadinterval;

	tmpreadinterval = (READINTERVAL *) ckalloc(edge1 -> multip * sizeof(READINTERVAL));
	for(m = 0; m < edge1 -> multip; m ++)	{
		tmpreadinterval[m] = edge1 -> readinterval[m];
	}
	free((void *) edge1 -> readinterval);
	l = edge1 -> multip + edge -> multip;
	edge1 -> readinterval = (READINTERVAL *) ckalloc(l * sizeof(READINTERVAL));
	for(m = 0; m < edge1 -> multip; m ++)	{
		edge1 -> readinterval[m] = tmpreadinterval[m];
	}
	for(m = 0; m < edge -> multip; m ++)	{
		edge1 -> readinterval[edge1 -> multip ++] = edge -> readinterval[m];
	}
	free((void *) tmpreadinterval);
}
