/***************************************************************************
 * Title:          width.c
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

int countwidth(NODES *node);
int countthickness(NODES *node);
void countallnode(READLIST **readlist, int *len_seq, int num_seq);

int countwidth(NODES *node)
{
	int	i, j, k, l, n, width;
	int	read, readold;
	int	num_read, *readall, *widthall;
	READPOSITION	*readposition;

	sort_nodepos(node);
	readall = (int *) ckalloc(node -> npos * sizeof(int));
	widthall = (int *) ckalloc(node -> npos * sizeof(int));
	readposition = node -> readposition;
	num_read = 0;
	k = 1;
	readold = -1;
	while(readposition)	{
		read = readposition -> readindex;
		if(read != readold)	{
			readall[num_read] = read;
			widthall[num_read] = k;
			num_read ++;
			if(readposition -> position >= 0)	{
				k = 1;
			} else	{
				k = 0;
			}
			readold = read;
		} else if(readposition -> position >= 0)	{
			k ++;
		}
		readposition = readposition -> next;
	}
	width = 0;
	for(i = 0; i < num_read; i ++)	{
		if(widthall[i] > width)	{
			width = widthall[i];
		}
	}
	free((void *) readall);
	free((void *) widthall);
	return(width);
}

int countthickness(NODES *node)
{
	int	i, j, k, l, n, thickness;
	int	read, readold;
	int	mini, maxi;
	READPOSITION	*readposition, *pos1;

	sort_nodepos(node);
	readposition = node -> readposition;
	mini = MAX_TMP_LEG;
	while(readposition)	{
		read = readposition -> readindex;
		pos1 = readposition -> next;
		while(pos1 && pos1 -> readindex == read)	{
			if(mini > abs(readposition -> position - pos1 -> position) + 1)	{
				mini = abs(readposition -> position - pos1 -> position) + 1;
			}
			pos1 = pos1 -> next;
		}
		readposition = readposition -> next;
	}
	if(mini >= WhirlLength)	{
		thickness = 1;
	} else	{
		thickness = mini;
	}
/*
	if(mini == MAX_TMP_LEG)	thickness = 1;
	else			thickness = mini;
*/
	return(thickness);
}

void countallnode(READLIST **readlist, int *len_seq, int num_seq)
{
	int	i, j, k, l, n;

	for(i = 0; i < num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			if(readlist[i][j].node && readlist[i][j].node -> visit == 0)	{
				readlist[i][j].node -> visit = 1;
				n = countwidth(readlist[i][j].node);
				readlist[i][j].node -> num_path = n;
			}
		}
	}
	cleannode(readlist, len_seq, num_seq/2);
}
