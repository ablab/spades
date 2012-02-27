/***************************************************************************
 * Title:          pairrule.c
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

char pairrules(char ext_name, int plate, char *library, int *range);

char pairrules(char ext_name, int plate, char *library, int *range)
{
	int	i, j, k, l;

	for(i = 0; i < ntypepair; i ++)	{
		if(pairrule[i][0] == ext_name && plate >= platerule[i][0] && plate <= platerule[i][1] &&
		   (!strcmp(librule[i], "ALL") || !strcmp(library, librule[i])))	{
			range[0] = i;
			return(pairrule[i][1]);
		}
	}

	return(ext_name);
}
