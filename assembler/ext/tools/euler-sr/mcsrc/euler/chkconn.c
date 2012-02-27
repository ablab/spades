/***************************************************************************
 * Title:          chkconn.c
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

int alignlength(ALIGN *align);
ALIGN *remove_trans(ALIGN *align1, ALIGN *align2, ALIGN *align3);
ALIGN *chkconn(int i1, int i2, ALIGN **eq_class, int num_seq);

ALIGN *chkconn(int i1, int i2, ALIGN **eq_class, int num_seq)
{
	int	i, j, k, l;
	int	c1, c2;
	ALIGN	*align;

	align = eq_class[i1];
	while(align)	{
		if(align -> reads[1] == i2)	{
			return(align);
		}
		align = align -> next;
	}
	align = eq_class[i2];
	while(align)	{
		if(align -> reads[1] == i1)	{
			return(align);
		}
		align = align -> next;
	}
	c1 = reverse_read(i1, num_seq);
	c2 = reverse_read(i2, num_seq);
	align = eq_class[c1];
	while(align)	{
		if(align -> reads[1] == c2)	{
			return(align);
		}
		align = align -> next;
	}
	align = eq_class[c2];
	while(align)	{
		if(align -> reads[1] == c1)	{
			return(align);
		}
		align = align -> next;
	}
	return((ALIGN *) NULL);
}

ALIGN *remove_trans(ALIGN *align1, ALIGN *align2, ALIGN *align3)
{
	int	i1, i2, i3, n;

	i1 = alignlength(align1);
	i2 = alignlength(align2);
	i3 = alignlength(align3);
	if(i1 <= i2 && i1 <= i3)	{
		return(align1);
	} else if(i2 <= i3)	{
		return(align2);
	} else	{
		return(align3);
	}
}

int alignlength(ALIGN *align)
{
	int	n;

	n = min(align -> pos[0][align -> length - 1] - align -> pos[0][0] + 1,
	        align -> pos[1][align -> length - 1] - align -> pos[1][0] + 1);
	return(n);
}

