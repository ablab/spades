/***************************************************************************
 * Title:          erasenext.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <assert.h>
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

void erasenext(NODES *vertex, int n);
void eraselast(NODES *vertex, int n);
int searcherase(EDGE **edge, EDGE *e, int num);

void erasenext(NODES *vertex, int n)
{
	int	i;

	for(i = n; i < vertex -> num_nextedge - 1; i ++)	{
		vertex -> nextedge[i] = vertex -> nextedge[i + 1];
	}
	vertex -> num_nextedge --;
}

void eraselast(NODES *vertex, int n)
{
	int	i;

	for(i = n; i < vertex -> num_lastedge - 1; i ++)	{
		vertex -> lastedge[i] = vertex -> lastedge[i + 1];
	}
	vertex -> num_lastedge --;
}

int searcherase(EDGE **edge, EDGE *e, int num)
{
	int	i, n;
	/* 
	   Look for edge e within the list of edges 'edge'.
	   Return the index of that edge.
	*/
	n = -1;
	for(i = 0; i < num; i ++)	{
		if(edge[i] == e)	{
			n = i;
			break;
		}
	}

	if(i == num)	{
		printf("e %d %d %d\n", e, e -> begin, e -> end);
		for(i = 0; i < num; i ++)	{
			printf("%d %d\n", i, edge[i]);
		}
		printf("Not found\n");
		assert(0);
	}

	return(n);
}
