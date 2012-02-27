/***************************************************************************
 * Title:          consensus.c
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

#define CORRSTEPS 3

void consensus(NODES **vertex, int num_vertex, char **src_seq, int num_seq, int *len_seq);
int revise_consensus(EDGE *edge, char **src_seq, int num_seq);
char profileedge(char **edgeprof, char *seq, int len, int offset, char *src_seq, int length);

void consensus(NODES **vertex, int num_vertex, char **src_seq, int num_seq, int *len_seq)
{
	int	i, j, k, l, m, n, num_corr;
	EDGE	*edge, *bal_edge;

	initial_edge(vertex, num_vertex, src_seq, len_seq, num_seq);
	n = 0;
	do	{
		n ++;
		num_corr = 0;
		for(i = 0; i < num_vertex; i ++)	{
			for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
				edge = vertex[i] -> nextedge[j];
				edge -> visit = 0;
			}
		}
		for(i = 0; i < num_vertex; i ++)	{
			for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
				edge = vertex[i] -> nextedge[j];
				if(edge -> visit == 1)	continue;
				k = revise_consensus(edge, src_seq, num_seq);
				edge -> visit = 1;
				num_corr += k;
				bal_edge = edge -> bal_edge;
				if(bal_edge != edge)	{
					free((void *) bal_edge -> seq);
					bal_edge -> length = edge -> length;
					bal_edge -> seq = (char *) ckalloc(bal_edge -> length * sizeof(char));
					for(m = 0; m < edge -> length; m ++)	{
						bal_edge -> seq[m] = rev(edge -> seq[edge -> length - m - 1]);
					}
					bal_edge -> visit = 1;
				}
			}
		}
		printf("Step: %d; Corrected errors by consensus: %d.\n", n, num_corr);
	} while(num_corr > 0 && n < CORRSTEPS);
}

int revise_consensus(EDGE *edge, char **src_seq, int num_seq)
{
	int	i, j, k, l, m, n;
	char	**edgeprof;
	READINTERVAL	*readinterval;

	l = edge -> length;
	edgeprof = (char **) ckalloc(edge -> length * sizeof(char *));
	for(i = 0; i < l; i ++)
		edgeprof[i] = (char *) ckalloc(10 * sizeof(char));
	for(i = 0; i < edge -> multip; i ++)	{
		readinterval = &(edge -> readinterval[i]);
		profileedge(edgeprof, edge -> seq, edge -> length, readinterval -> offset, 
			&src_seq[readinterval -> eq_read][readinterval -> begin], readinterval -> length);
	}
	n = updateedge(edgeprof, edge);
	for(i = 0; i < l; i ++)
		free((void *) edgeprof[i]);
	free((void **) edgeprof);
	return(n);
}

char profileedge(char **edgeprof, char *seq, int len, int offset, char *src_seq, int length)
{
	int	i, j, k, l, m, n, t, l0;
	int	pos1[2], pos2[2], l1, l2;
	int	*sapp;

	offset = max(offset - band, 0);
	sapp = (int *) ckalloc(length * 2 * sizeof(int));
	l0 = min(length + band * 2, len - offset);
	LOCAL_ALIGN0(&seq[offset - 1], src_seq - 1, l0, length, -band * 4, band * 4, W, g, h, &pos1[0], &pos2[0], 
		&pos1[1], &pos2[1], length * 2);
	pos1[0] --;
	pos1[1] --;
	pos2[0] --;
	pos2[1] --;
	offset += pos1[0];
	l1 = pos1[1] - pos1[0] + 1;
	l2 = pos2[1] - pos2[0] + 1;
	ALIGN0(&seq[offset - 1], &src_seq[pos2[0] - 1], l1, l2, -band * 2, band * 2, W, g, h, sapp, length * 2, length * 2);
	k = t = 0;
	i = offset;
	j = pos2[0];
	while(i < offset + l1 || j < pos2[1])	{
		if(sapp[k] == 0)	{
			edgeprof[i][src_seq[j]] ++;
			i ++;
			j ++;
			t ++;
		} else if(sapp[k] < 0)	{
			if(j > pos2[0])	{
				edgeprof[i][4] ++;
			}
			i -= sapp[k];
		} else	{
			if(i > offset)	{
				edgeprof[i][5 + src_seq[j]] ++;
			}
			j += sapp[k];
		}
		k ++;
	}
	free((void *) sapp);
}

int updateedge(char **edgeprof, EDGE *edge)
{
	int	i, j, k, l, n, max_k, n1, n2, n3;
	int	*posindex;
	char	*newseq;

	l = n = 0;
	newseq = (char *) ckalloc(2 * edge -> length * sizeof(char));
	posindex = (int *) ckalloc(edge -> length * sizeof(int));
	for(i = 0; i < edge -> length; i ++)	{
		k = max_k = 0;
		for(j = 0; j < 9; j ++)	{
			if(edgeprof[i][j] > k)	{
				k = edgeprof[i][j];
				max_k = j;
			}
		}
		if(max_k < 4)	{
			if(k == 0)	max_k = edge -> seq[i];
			if(max_k != edge -> seq[i])	{
				n ++;
			}
			posindex[i] = l;
			newseq[l ++] = max_k;
		} else if(max_k == 4)	{
			posindex[i] = l;
			if(k <= 1)	{
				newseq[l ++] = 1;
				n ++;
			}
		} else {
			posindex[i] = l + 1;
			if(k > 1)	 {
				newseq[l ++] = max_k - 5;
				n ++;
			}
			newseq[l ++] = edge -> seq[i];
		}
	}
	free((void *) edge -> seq);
	edge -> seq = (char *) ckalloc(l * sizeof(char));
	for(i = 0; i < l; i ++)	{
		edge -> seq[i] = newseq[i];
	}
	edge -> length = l;
	for(i = 0; i < edge -> multip; i ++)	{
		edge -> readinterval[i].offset = posindex[edge -> readinterval[i].offset];
	}
	free((void *) newseq);
	free((void *) posindex);
	return(n);
}
