/***************************************************************************
 * Title:          index.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <stdinc.h>
#include <extfunc.h>


INDEX *insert_index(INDEX *index, int p);
INDEX *free_index(INDEX *index);

INDEX *insert_index(INDEX *index, int p)
{
	INDEX	*index0;

	index0 = (INDEX *) ckalloc(1 * sizeof(INDEX));
	index0 -> index = p;
	index0 -> next = index;
	return(index0);
}

INDEX *free_index(INDEX *index)
{
	INDEX	*index0;

	index0 = index -> next;
	free((void *) index);
	return(index0);
}
