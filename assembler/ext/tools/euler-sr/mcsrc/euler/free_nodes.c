/***************************************************************************
 * Title:          free_nodes.c
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

NODES *free_nodes(NODES *node);

NODES *free_nodes(NODES *node)
{
	READPOSITION	*readposition;

	while(node -> readposition)	{
		readposition = node -> readposition -> next;
		free((void *) node -> readposition);
		node -> readposition = readposition;
	}
	if(node -> path_pos)	free((void **) node -> path_pos);
	if(node -> path_index)	free((void **) node -> path_index);
	if(node -> lastedge)	free((void **) node -> lastedge);
	if(node -> nextedge)	free((void **) node -> nextedge);
	free((void *) node);
	return((NODES *) NULL);
}
