/***************************************************************************
 * Title:          free_class.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>

READOVERLAP *free_readoverlap(READOVERLAP *readoverlap);
READINTERVAL *free_readinterval(READINTERVAL *readinterval);
int size_readinterval(READINTERVAL *readinterval);

READOVERLAP *free_readoverlap(READOVERLAP *readoverlap)
{
	READOVERLAP	*cl;

	cl = readoverlap -> next;
	free((void *) readoverlap);
	return(cl);
}

READINTERVAL *free_readinterval(READINTERVAL *readinterval)
{
	READINTERVAL	*cl;

	cl = readinterval -> next;
	free((void *) readinterval);
	return(cl);
}

int size_readinterval(READINTERVAL *readinterval)
{
	int	n;
	READINTERVAL	*cl;

	n = 0;
	cl = readinterval;
	while(cl)	{
		n ++;
		cl = cl -> next;
	}
	return(n);
}
