/***************************************************************************
 * Title:          initialize_edge.c
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

void initial_edge(NODES **vertex, int num_vertex, char **src_seq, int *len_seq, int num_seq);
void initedge(EDGE *edge, int *len_seq, char **src_seq);

void initial_edge(NODES **vertex, int num_vertex, char **src_seq, int *len_seq, int num_seq)
{
	int	i, j, k, l;
	EDGE	*edge, *bal_edge;

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			bal_edge = edge -> bal_edge;
			initedge(edge, len_seq, src_seq);
			if(bal_edge != edge)	{
				bal_edge -> length = edge -> length;
				bal_edge -> seq = (char *) ckalloc(bal_edge -> length * sizeof(char));
				for(k = 0; k < edge -> length; k ++)	{
					bal_edge -> seq[edge -> length - k - 1] = rev(edge -> seq[k]);
				}
			}
		}
	}
}

void initedge(EDGE *edge, int *len_seq, char **src_seq)
{
	int	i, j, k, l, n, len;
	int	reads, pos, min, max;
	char	*seq;

	seq = (char *) ckalloc(edge -> length * sizeof(char));
	/*	for (i = 0; i < edge->length; i++) 
		seq[i] = 5;
	*/

	sortreadinterval(edge -> readinterval, edge -> multip);
	for(i = 0; i < edge -> multip; i ++)	{
		if (edge->readinterval[i].eq_read >= 0) {
			reads = edge -> readinterval[i].eq_read;
			pos = edge -> readinterval[i].offset;
			l = edge -> readinterval[i].length;
			
			for(j = 0; j < l; j ++)	{
				/*
					
				*/
				if(edge -> readinterval[i].begin + j < len_seq[reads] &&
			   pos + j < edge -> length && pos + j >= 0)
					seq[pos + j] = src_seq[reads][edge -> readinterval[i].begin + j] + 1;
			}
		}
	}
	edge -> seq = (char *) ckalloc((edge -> length + 1) * sizeof(char));
	
	min = max = 0;
	for(i = 0; i < edge -> length; i ++)	{
		if(seq[i] > 0)	{
			min = i;
			break;
		}
	}
	for(i = edge -> length - 1; i >= 0; i --)	{
		if(seq[i] > 0)	{
			max = i;
			break;
		}
	}
	n = 0;
	for(i = min; i <= max; i ++)	{
		if(seq[i] > 0)	{
			edge -> seq[n ++] = seq[i] - 1;
		}
	}
	if(n > 0)	{
		edge -> length = n;
	}
	free((void *) seq);
}
