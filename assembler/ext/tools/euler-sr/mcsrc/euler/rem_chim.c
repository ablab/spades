/***************************************************************************
 * Title:          rem_chim.c
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

void rem_chim(NODES **vertex, int num_vertex, int *chim, int num_chim, int num_seq);
int chk_chimi2(int *chim, int num_chim, int index, int num_seq);

void rem_chim(NODES **vertex, int num_vertex, int *chim, int num_chim, int num_seq)
{
	int	i, j, k, l, m, n, c;
	EDGE	*edge;

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			m = 0;
			while(m < edge -> multip)	{
				c = chk_chim2(chim, num_chim, edge -> readinterval[m].eq_read, num_seq);
				if(c)	{
					edge -> multip --;
					edge -> readinterval[m] = edge -> readinterval[edge -> multip];
				} else	{
					m ++;
				}
			}
		}
	}
}

int chk_chim2(int *chim, int num_chim, int index, int num_seq)
{
	int	i;

	for(i = 0; i < num_chim; i ++)	{
		if(chim[i] == index || reverse_read(chim[i], num_seq) == index)	return(1);
	}
	return(0);
}
