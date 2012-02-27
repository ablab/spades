/***************************************************************************
 * Title:          trimpath.c
 * Author:         Haixu Tang
 * Created:        Jun. 2003
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

#define MIN_INTV 10

void trimpath(PATH *path, PATH *path_rev, int readindex, int readindex_rev);
int findreadinterval(EDGE *edge, int readindex);
void remove_readinterval_len(EDGE *edge, int index);

void trimpath(PATH *path, PATH *path_rev, int readindex, int readindex_rev)
{
	int	i, j, k, l, m, n, k1, k2;
	int	len1, len2;

	l = path -> len_path;
	k1 = l;
	k2 = 0;
	/* 
		 Find the first edge on path 'readindex' that has length > MIN_INTV.
	*/
	for(i = 0; i < l; i ++)	{
		j = findreadinterval(path -> edge[i], readindex);
		len1 = path -> edge[i] -> readinterval[j].length;
		if(len1 >= MIN_INTV || len1 == path -> edge[i] -> length)	{
			k1 = i;
			break;
		}
	}
	/* 
		 Find the last edge on path 'readindex' that has length > MIN_INTV.
	 */
	j = findreadinterval(path -> edge[l - 1], readindex);
	for(i = l - 1; i >= 0; i --)	{
		j = findreadinterval(path -> edge[i], readindex);
		len1 = path -> edge[i] -> readinterval[j].length;
		if(len1 >= MIN_INTV || len1 == path -> edge[i] -> length)	{
			k2 = i;
			break;
		}
	}
	/*
		If all intervals are short, just ignore this path. 
	*/
	if(k1 > k2)	{
		path -> len_path = path_rev -> len_path = 0;
		return;
	}
	/*
		Remove all short intervals at the beginnign of a path.
	*/
	for(i = 0; i < k1; i ++)	{
		remove_readinterval_len(path -> edge[i], readindex);
		remove_readinterval_len(path_rev -> edge[l - i - 1], readindex_rev);
	}
	/* 
		 Remove all short intervals at the end of a path.
	*/
	for(i = k2 + 1; i < l; i ++)	{
		remove_readinterval_len(path -> edge[i], readindex);
		remove_readinterval_len(path_rev -> edge[l - i - 1], readindex_rev);
	}
	/* 
		 Pack the path to the beginning of the path array.
	*/
	for(i = k1; i <= k2; i ++)	{
		path -> edge[i - k1] = path -> edge[i];
		path_rev -> edge[k2 - i] = path -> edge[i] -> bal_edge;
	}
	path_rev -> len_path = path -> len_path = k2 - k1 + 1;
	if(k1 > 0)	{
		path -> begin_length = 0;
		path_rev -> end_length = 0;
	}
	if(k2 < l - 1)	{
		path -> end_length = 0;
		path_rev -> begin_length = 0;
	}
}

int findreadinterval(EDGE *edge, int readindex)
{
	int	i, j;

	for(i = 0; i < edge -> multip; i ++)	{
		if(edge -> readinterval[i].eq_read == readindex)	{
			return(i);
		}
	}
	printf("Read interval not found.\n");
	exit(0);
}

void remove_readinterval_len(EDGE *edge, int index)
{
	int	i, j, k, l;

	i = 0;
	while(i < edge -> multip)	{
		if(edge -> readinterval[i].eq_read == index && edge -> readinterval[i].length < edge -> length)	{
			edge -> readinterval[i] = edge -> readinterval[edge -> multip - 1];
			edge -> multip --;
		} else	{
			i ++;
		}
	}
}

