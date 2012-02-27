/***************************************************************************
 * Title:          cleangraph.c
 * Author:         Haixu Tang
 * Created:        Jun. 2003
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

void cleangraph(NODES **nodes, int num_nodes);

void cleangraph(NODES **nodes, int num_nodes)
{
	int	i, j, k;

	for(i = 0; i < num_nodes; i ++)	{
		nodes[i] -> subg_flag = SUBG_IN;
		for(j = 0; j < nodes[i] -> num_nextedge; j ++)	{
			nodes[i] -> nextedge[j] -> subg_flag = SUBG_IN;
		}
	}
}
