/***************************************************************************
 * Title:          readoverlap.c
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

int readoverlap(READINTERVAL **readinterval, FILE *fp);

int readoverlap(READINTERVAL **readinterval, FILE *fp)
{
	int	i, j, k, l, n;
	int	n1, n2, n3, n4;
	int	offset;
	READINTERVAL	*nreadinterval;
	char	str[100];

	n = 0;
	while(fgets(str, 90, fp))	{
		sscanf(str, "%d%d%d%d", &n1, &n2, &n3, &n4);
		offset = (n3 + n4) / 2;
		k = max(n1, n2);
		j = min(n1, n2);
		nreadinterval = (READINTERVAL *) ckalloc(1 * sizeof(READINTERVAL));
		nreadinterval -> eq_read = j;
		if(k == n1)	{
			nreadinterval -> offset = offset;
		} else	{
			nreadinterval -> offset = -offset;
		}
		nreadinterval -> next = readinterval[k];
		readinterval[k] = nreadinterval;
		n ++;
	}
	return(n);
}
