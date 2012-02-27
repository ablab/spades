/***************************************************************************
 * Title:          visit_comp.c
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

void visit_comp(NODES **node, NODES *begin_node, int *comp, int num_comp, int num_seq);

void visit_comp(NODES **node, NODES *begin_node, int *comp, int num_comp, int num_seq)
{
	int	i, j, k;
	NODES	*nextnode;

	k = begin_node -> npos;
	if(k >= num_seq)	k -= num_seq; 
	comp[k] = num_comp;
	for(i = 0; i < begin_node -> num_nextedge; i ++)	{
		nextnode = begin_node -> nextedge[i] -> end;
		k = nextnode -> npos;
		if(k >= num_seq)	k -= num_seq; 
		if(comp[k] == 0)	{
			visit_comp(node, nextnode, comp, num_comp, num_seq);
		}
	}
	for(i = 0; i < begin_node -> num_lastedge; i ++)	{
		nextnode = begin_node -> lastedge[i] -> begin;
		k = nextnode -> npos;
		if(k >= num_seq)	k -= num_seq; 
		if(comp[k] == 0)	{
			visit_comp(node, nextnode, comp, num_comp, num_seq);
		}
	}
}
