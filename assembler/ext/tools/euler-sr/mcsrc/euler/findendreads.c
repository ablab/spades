/***************************************************************************
 * Title:          findendreads.c
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

#define SHORT_LEN 1000

static char	*label;

void findendreads(NODES **vertex, int num_vertex, int *endreads, int *beginreads, int *num_endreads, int num_seq);
void findmatereads(NODES **vertex, int num_vertex, int *endreads, int *beginreads, int *num_endreads, int num_seq, MATEPAIRTABLE *MP);
int findmate(int k, MATEPAIRTABLE *MP);

void findendreads(NODES **vertex, int num_vertex, int *endreads, int *beginreads, int *num_endreads, int num_seq)
{
	int	i, j, k, l, m;
	EDGE 	*edge;

	label = (char *) ckalloc(2 * num_seq * sizeof(char));
	num_endreads[0] = num_endreads[1] = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(m = 0; m < vertex[i] -> num_nextedge; m ++)	{
			edge = vertex[i] -> nextedge[m];
			if(edge -> multip > 1)	{
/*
				if(edge -> begin -> num_lastedge == 0 || edge -> length < SHORT_LEN)	{
*/
				if(edge -> begin -> num_lastedge == 0)	{
					for(j = 0; j < edge -> multip; j ++)	{
						if(edge -> readinterval[j].offset < MIN_END2)	{
							k = edge -> readinterval[j].eq_read;	
							if(!label[k])	{
								beginreads[num_endreads[0] ++] = k;
								l = reverse_read(k, num_seq);
								endreads[num_endreads[1] ++] = l;
								label[k] = 1;
							}
						}
					}
				}
			}
		}
	}
}

void findmatereads(NODES **vertex, int num_vertex, int *endreads, int *beginreads, int *num_endreads, int num_seq, MATEPAIRTABLE *MP)
{
	int	i, j, k, l, m, n, c;
	EDGE 	*edge;

	for(i = 0; i < num_vertex; i ++)	{
		for(m = 0; m < vertex[i] -> num_nextedge; m ++)	{
			edge = vertex[i] -> nextedge[m];
			for(j = 0; j < edge -> multip; j ++)	{
				k = edge -> readinterval[j].eq_read;	
				label[k] = 1;
			}
		}
	}
	for(i = 0; i < num_vertex; i ++)	{
		for(m = 0; m < vertex[i] -> num_nextedge; m ++)	{
			edge = vertex[i] -> nextedge[m];
			if(edge -> multip > 1 && (vertex[i] -> num_lastedge == 0 || edge -> length < SHORT_LEN))	{
				for(j = 0; j < edge -> multip; j ++)	{
					k = edge -> readinterval[j].eq_read;	
					if(k < num_seq)	{
						l = findmate(k, MP);
						n = l + num_seq;
					} else	{
						n = findmate(k - num_seq, MP);
						l = n + num_seq;
					}
					if(n >= 0 && !label[l])	{
						endreads[num_endreads[1] ++] = n;
						beginreads[num_endreads[0] ++] = l;
						label[l] = 1;
					}
				}
			}
		}
	}
}

int findmate(int k, MATEPAIRTABLE *MP)
{
	int	i, j;

	for(i = 0; i < MP -> num_pair; i ++)	{
		if(MP -> pair1[i] == k)	{
			return(MP -> pair2[i]);
		} else if(MP -> pair2[i] == k)	{
			return(MP -> pair1[i]);
		}
	}
	return(-1);
}

void free_label()
{
	free((void *) label);
}
