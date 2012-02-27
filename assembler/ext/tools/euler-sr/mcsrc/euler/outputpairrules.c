/***************************************************************************
 * Title:          outputpairrules.c
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

void outputpairrules();

void outputpairrules()
{
	int	i;

	printf("EULER mate-pair INPUT:\n\n");
	for(i = 0; i < ntypepair; i ++)	{
		printf("%c-%c(plate number from %d to %d) %d %d\n", pairrule[i][0], pairrule[i][1], 
			platerule[i][0], platerule[i][1], pairrange[i][0], pairrange[i][1]);
	}
}
