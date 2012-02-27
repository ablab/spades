/***************************************************************************
 * Title:          readclass.c
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

int readclass(ALIGN **eq_class, int num_seq, FILE *fp);

int readclass(ALIGN **eq_class, int num_seq, FILE *fp)
{
	int	i, j, k, l, n;
	int	num_readinterval;
	ALIGN	*readinterval;

	l = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		fread(&n, sizeof(int), 1, fp);
		l += n;
		for(j = 0; j < n; j ++)	{
			readinterval = (ALIGN *) ckalloc(1 * sizeof(ALIGN));
			readinterval -> reads[0] = i;
			fread(&(readinterval -> reads[1]), sizeof(int), 1, fp);
			fread(&(readinterval -> mis_match), sizeof(int), 1, fp);
			fread(&(readinterval -> length), sizeof(int), 1, fp);
			for(k = 0; k < 2; k ++)	{
				readinterval -> pos[k] = (int *) ckalloc(readinterval -> length * sizeof(int));
			}
			fread(readinterval -> pos[0], sizeof(int), readinterval -> length, fp);
			fread(readinterval -> pos[1], sizeof(int), readinterval -> length, fp);
			readinterval -> next = eq_class[i];
			readinterval -> prev = NULL;
			if(eq_class[i])	{
				eq_class[i] -> prev = readinterval;
			}
			eq_class[i] = readinterval;
		}
	}
	return(l);
}
