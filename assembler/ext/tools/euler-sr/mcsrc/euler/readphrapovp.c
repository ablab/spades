/***************************************************************************
 * Title:          readphrapovp.c
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

int readphrapovp(READOVERLAP **readoverlap, int threshold, FILE *fp);

int readphrapovp(READOVERLAP **readoverlap, int threshold, FILE *fp)
{
	int	i, j, k, l, n;
	int	n1, n2, n3, n4, n5, n6;
	int	offset;
	READOVERLAP	*nreadoverlap;
	char	str[100];

	n = 0;
	while(fgets(str, 90, fp))	{
		sscanf(str, "%d%d%d%d%d%d%d", &n1, &n2, &n3, &n4, &n5, &n6, &k);
		if(k < threshold)	continue;
		nreadoverlap = (READOVERLAP *) ckalloc(1 * sizeof(READOVERLAP));
		nreadoverlap -> read1 = n1;
		nreadoverlap -> read2 = n2;
		nreadoverlap -> pos1[0] = n3;
		nreadoverlap -> pos1[1] = n4;
		nreadoverlap -> pos2[0] = n5;
		nreadoverlap -> pos2[1] = n6;
		nreadoverlap -> LLR_score = k;
		nreadoverlap -> next = readoverlap[n1];
		readoverlap[n1] = nreadoverlap;
		n ++;
	}
	return(n);
}
