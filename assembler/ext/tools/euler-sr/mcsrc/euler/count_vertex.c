/***************************************************************************
 * Title:          count_vertex.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <stdinc.h>
#include <assert.h>

int count_vertex(EDGE **edge, int num_edge, NODES **vertex);
int collect_vertex(NODES *v, NODES **vertex, int num_vertex);
EDGE *find_bal_edge(EDGE *edge, int *len_seq, int num_seq, int index);
char chk_readinterval(READINTERVAL *readinterval1, READINTERVAL *readinterval2, int n1, int n2, int *len_seq, int num_seq);

int count_vertex(EDGE **edge, int num_edge, NODES **vertex)
{
	printf("THIS CODE IS NOT CALLED\n");
	assert(0);
	int	i, j, k, l;
	int	num_vertex;

	for(i = 0; i < num_edge; i ++)	{
		edge[i] -> begin -> visit = edge[i] -> end -> visit = 0;
	}

	num_vertex = 0;
	for(i = 0; i < num_edge; i ++)	{
		num_vertex = collect_vertex(edge[i] -> begin, vertex, num_vertex);
		num_vertex = collect_vertex(edge[i] -> end, vertex, num_vertex);
	}
	return(num_vertex);
}

int collect_vertex(NODES *v, NODES **vertex, int num_vertex)
{
	if(v -> visit == 1)	return(num_vertex);
	vertex[num_vertex ++] = v;
	v -> visit = 1;
	return(num_vertex);
}

EDGE *find_bal_edge(EDGE *edge, int *len_seq, int num_seq, int index)
{
	int	i, j, k, l;
	char	c;
	NODES	*node, *bal_node;
	EDGE	*edge1;

	node = edge -> begin;
	bal_node = node -> bal_node;

	for(i = 0; i < bal_node -> num_lastedge; i ++)	{
		if(bal_node -> lastedge[i] -> length != edge -> length ||
		   bal_node -> lastedge[i] -> begin != edge -> end -> bal_node ||
		   bal_node -> lastedge[i] -> multip != edge -> multip)
			continue;
		c = chk_readinterval(edge -> readinterval, bal_node -> lastedge[i] -> readinterval, edge -> multip,
			      bal_node -> lastedge[i] -> multip, len_seq, num_seq);
		if(c)	{
			return(bal_node -> lastedge[i]);
		}
	}
	printf("%d Bal_edge not found. %d %d %d %d(%d-%d)-%d(%d-%d)\n", index,
		edge, edge -> multip, edge -> length, edge -> begin, edge -> begin -> num_lastedge,
		edge -> begin -> num_nextedge, edge -> end, edge -> end -> num_lastedge, edge -> end -> num_nextedge);
	for(j = 0; j < edge -> multip; j ++)	{
		printf("read %d %d position %d %d %d %d\n", edge -> readinterval[j].eq_read,
			reverse_read(edge -> readinterval[j].eq_read, num_seq),
			len_seq[edge -> readinterval[j].eq_read],
			edge -> readinterval[j].begin, edge -> readinterval[j].length, edge -> readinterval[j].offset);
	}
	for(i = 0; i < bal_node -> num_lastedge; i ++)	{
		edge1 = bal_node -> lastedge[i];
		printf("Bal_edge %d %d %d %d(%d-%d)-%d(%d-%d)\n",
			edge1, edge1 -> multip, edge1 -> length, edge1 -> begin, edge1 -> begin -> num_lastedge,
			edge1 -> begin -> num_nextedge, edge1 -> end, edge1 -> end -> num_lastedge,
			edge1 -> end -> num_nextedge);
		for(j = 0; j < edge1 -> multip; j ++)	{
			printf("read %d %d position %d %d %d %d\n", edge1 -> readinterval[j].eq_read,
				reverse_read(edge1 -> readinterval[j].eq_read, num_seq),
				len_seq[edge1 -> readinterval[j].eq_read],
				edge1 -> readinterval[j].begin, edge1 -> readinterval[j].length, edge1 -> readinterval[j].offset);
		}
	}
	printf("end %d(%d) %d %d\n", edge -> end, edge -> end -> visit, edge -> end -> num_lastedge,
		edge -> end -> num_nextedge);
	printf("end Bal %d(%d) %d %d\n", edge -> end -> bal_node, edge -> end -> bal_node -> visit,
		edge -> end -> bal_node -> num_lastedge, edge -> end -> bal_node -> num_nextedge);
	exit(-1);
}

char chk_readinterval(READINTERVAL *readinterval1, READINTERVAL *readinterval2, int n1, int n2, int *len_seq, int num_seq)
{
	int	i, j, k, i1, i2;

	i1 = readinterval1[0].eq_read;
	i2 = readinterval1[0].begin;
	i1 = reverse_read(i1, num_seq);
	i2 = len_seq[i1] - i2 - readinterval1[0].length;
	for(i = n2 - 1; i >= 0; i --)	{
		if(readinterval2[i].eq_read == i1 && readinterval2[i].begin == i2 && readinterval1[0].length ==
		   readinterval2[i].length)
			return(1);
	}
	return(0);
}
