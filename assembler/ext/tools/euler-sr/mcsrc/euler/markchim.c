/***************************************************************************
 * Title:          markchim.c
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

int chk_chim(int *reads, int *chim, int num_chim);
void markchim(ALIGN **eq_class, int num_seq, int *chim, int num_chim);

void markchim(ALIGN **eq_class, int num_seq, int *chim, int num_chim)
{
	int	i, j, k, l;
	ALIGN	*align, *align0, *align1;

	if (chim == (int *) 0) return;

	for(i = 0; i < num_seq * 2; i ++)	{
		align = eq_class[i];
		align0 = (ALIGN *) NULL;
		while(align)	{
			k = chk_chim(align -> reads, chim, num_chim);
			if(k)	{
				align1 = align -> next;
				if(!align0)	{
					eq_class[i] = align1;
				} else	{
					align0 -> next = align1;
				}
				free_align(align);
				align = align1;
			} else	{
				align0 = align;
				align = align -> next;
			}
		}
	}
}

int chk_chim(int *reads, int *chim, int num_chim)
{
	int	i;

	for(i = 0; i < num_chim; i ++)	{
		if(chim[i] == reads[0] || chim[i] == reads[1])	return(1);
	}
	return(0);
}
