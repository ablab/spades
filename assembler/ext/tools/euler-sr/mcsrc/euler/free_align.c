/***************************************************************************
 * Title:          free_align.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <stdinc.h>

ALIGN *free_align(ALIGN *align);
int size_align(ALIGN *align);

ALIGN *free_align(ALIGN *align)
{
	int	i;
	ALIGN	*cl;

	cl = align -> next;
	for(i = 0; i < 2; i ++)
		free((void *) align -> pos[i]);
	if(align -> prev)	{
		align -> prev -> next = align -> next;
	}
	if(align -> next)	{
		align -> next -> prev = align -> prev;
	}
	free((void *) align);
	return(cl);
}

int size_align(ALIGN *align)
{
	int	n;
	ALIGN	*cl;

	n = 0;
	cl = align;
	while(cl)	{
		n ++;
		cl = cl -> next;
	}
	return(n);
}
