/***************************************************************************
 * Title:          countdegree.c
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

int countdegree(EDGE **edges, int num_edges);

int countdegree(EDGE **edges, int num_edges)
{
	int	i, n;

	n = 0;
	for(i = 0; i < num_edges; i ++)	{
		if(edges[i] -> subg_flag)	n ++;
	}
	return(n);
}
