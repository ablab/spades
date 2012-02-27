/***************************************************************************
 * Title:          trim_align.c
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

#define MIN_ID 0.90

void trim_align(ALIGN *align, char **src_seq, int *len_seq);

void trim_align(ALIGN *align, char **src_seq, int *len_seq)
{
	int	i, j, k, l;
	int	k1, k2;
	int	read1, read2, pos1, pos2;
	double	*identity;

	read1 = align -> reads[0];
	read2 = align -> reads[1];
	k = l = 0;
	identity = (double *) ckalloc(align -> length * sizeof(double));
	for(i = 0; i < align -> length - 1; i ++)	{
		pos1 = align -> pos[0][i];
		pos2 = align -> pos[1][i];
		while(pos1 < align -> pos[0][i + 1] && pos2 < align -> pos[1][i + 1])	{
			if(src_seq[read1][pos1] == src_seq[read2][pos2])	{
				k ++;
			} else	{
				l ++;
			}
			pos1 ++;
			pos2 ++;
		}
		l += (align -> pos[0][i + 1] - pos1 + align -> pos[1][i + 1] - pos2);
		identity[i] = ((double) l) / (l + k);
	}
	for(i = 0; i < align -> length - 1; i ++)	{
		if(identity[i] > MIN_ID)	break;
	}
	k1 = i;
	for(i = align -> length - 2; i >= 0; i --)	{
		if(identity[i] > MIN_ID)	break;
	}
	k2 = i;
	if(k1 < k2)	{
		for(i = k1; i <= k2 + 1; i ++)	{
			align -> pos[0][i - k1] = align -> pos[0][i];
			align -> pos[1][i - k1] = align -> pos[1][i];
		}
		align -> length = k2 - k1 + 2;
	}
	free((void *) identity);
}
