/***************************************************************************
 * Title:          makedge.c
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
#include <assert.h>
int makedge(NODES *node, EDGE **edge, int num_edge, READLIST **readlist);
EDGE *newedge(EDGE **midedge, int num_midedge, READLIST **readlist);
int copyreadinterval(READINTERVAL *readinterval, READINTERVAL *readcov);
void sortreadinterval(READINTERVAL *readinterval, int multip);
int readintervalcompar(const void *a, const void *b);
void sortreadinterval_index(READINTERVAL *readinterval, int multip);
int readintervalcompar_index(const void *a, const void *b);
READINTERVAL *insert_readcov(READINTERVAL *readcov, int readindex, int begin, int length, int offset);
void sort_nodepos(NODES *node);
int poscompar(const void *a, const void *b);

int FindRead(READINTERVAL *intervals, int numIntervals, int readIndex, int *intervalPos) {
	// Assume intervals are sorted by read index, 
	// do binay search until readindex is found.

	int begin, end;
	begin = 0; end = numIntervals;
	int cur;
	cur = (begin + end) / 2;
	while (begin < end and intervals[cur].eq_read != readIndex) {
		if (intervals[cur].eq_read > readIndex) 
			end = cur;
		else if (intervals[cur].eq_read < readIndex) {
			begin = cur+1;
		}
		cur = (end + begin) / 2;
		assert(cur < numIntervals);
	}

	if (cur >= numIntervals or intervals[cur].eq_read != readIndex) {
		*intervalPos = -1;
		return 0;
	}
	else {
		while (cur > 0 and intervals[cur-1].eq_read ==  readIndex)
			cur--;
	}
	assert(cur >= 0);
	*intervalPos = cur;
	return 1;
}


int makedge(NODES *node, EDGE **edge, int num_edge, READLIST **readlist)
{
	int	i, j, k, l;
	int	num_midedge;
	EDGE	**midedge, *nedge;
	NODES	*node1;

	node -> visit = 1;
	for(j = 0; j < node -> num_nextedge; j ++)	{
		num_midedge = 0;
		midedge = (EDGE **) ckalloc(MAX_TMP_LEG * sizeof(EDGE *));
		midedge[num_midedge ++] = node -> nextedge[j];
		node1 = midedge[num_midedge - 1] -> end;
		while(node1 -> num_lastedge == 1 && node1 -> num_nextedge == 1)	{
			midedge[num_midedge ++] = node1 -> nextedge[0];
			node1 = node1 -> nextedge[0] -> end;
		}
		if(num_midedge > 1)	{
			nedge = newedge(midedge, num_midedge, readlist);
			if(nedge)	{
				nedge -> visit = 1;
				edge[num_edge ++] = nedge;
			} else	{
				printf("No returned edge\n");
				exit(0);
			}
		} else if(midedge[0] -> visit == 0)	{
			nedge = midedge[0];
			nedge -> visit = 1;
			edge[num_edge ++] = nedge;
		}
		free((void **) midedge);
		if(node1 -> visit == 0)	{
			num_edge = makedge(node1, edge, num_edge, readlist);
		}
	}
	return(num_edge);
}

EDGE *newedge(EDGE **midedge, int num_midedge, READLIST **readlist)
{
	int	i, j, k, l, n;
	int	readindex;
	INSERT	*insert;
	READINTERVAL	*readcov, *readcov0;
	NODES	*node, *node1;
	EDGE	*tmpedge;
	READPOSITION *readposition;


/*	Sort the readpositions every nodes cover by the index of the reads	*/

	for(i = 0; i < num_midedge; i ++)	{
		sortreadinterval_index(midedge[i] -> readinterval, midedge[i] -> multip);
	}

/*	Define read coverage of the edge	*/

	readcov = (READINTERVAL *) NULL;
	k = 0;
	for(i = 0; i < num_midedge; i ++)	{
		for(j = 0; j < midedge[i] -> multip; j ++)	{
			readcov = insert_readcov(readcov, midedge[i] -> readinterval[j].eq_read, midedge[i] -> readinterval[j].begin,
					 midedge[i] -> readinterval[j].length, k);
		}
		k += midedge[i] -> length - 1;
	}

/*	Setup the beginning and ending vertices of new edge	*/

	tmpedge = midedge[0];
	if(num_midedge > 1)	{
		tmpedge -> end = midedge[num_midedge - 1] -> end;
		add_lastedge(tmpedge -> end, tmpedge);
	}

/*	remove single read covers, copy and sort read covers	*/
	n = size_readinterval(readcov);
	tmpedge -> readinterval = (READINTERVAL *) ckalloc(n * sizeof(READINTERVAL));
	tmpedge -> multip = copyreadinterval(tmpedge -> readinterval, readcov);
	if(tmpedge -> multip == 0)	{
		printf("No read cover: tmpedge %d length %d.\n", tmpedge, num_midedge);
		readposition = tmpedge -> begin -> readposition;
		while(readposition)	{
			printf("node %d readposition %d %d\n", tmpedge -> begin,
				readposition -> readindex, readposition -> position);
			readposition = readposition -> next;
		}
		readposition = tmpedge -> end -> readposition;
		while(readposition)	{
			printf("node %d readposition %d %d\n", tmpedge -> end,
				readposition -> readindex, readposition -> position);
			readposition = readposition -> next;
		}
		readcov0 = readcov;
		while(readcov0)	{
			printf("readcov %d %d %d %d\n", readcov0 -> eq_read, readcov0 -> length,
				readcov0 -> offset, readcov0 -> cov);
			readcov0 = readcov0 -> next;
		}
		exit(-1);
	}
	while(readcov)	{
		readcov = free_readinterval(readcov);
	}
	sortreadinterval(tmpedge -> readinterval, tmpedge -> multip);

/*	Set up the length of the edge	*/

	for(i = 1; i < num_midedge; i ++)	{
		tmpedge -> length += midedge[i] -> length - 1;
	}

/*	Delete nodes on the path	*/

	for(i = 1; i < num_midedge; i ++)	{
		node = midedge[i] -> begin;
		readposition = node -> readposition;
		while(readposition)	{
			k = readposition -> readindex;
			l = readposition -> position;
			if(l >= 0)	{
				readlist[k][l].node = (NODES *) NULL;
			}
			readposition = readposition -> next;
		}
	}
	if(i > 1)	{
		k = searcherase(midedge[i - 1] -> end -> lastedge, midedge[i - 1],
				 midedge[i - 1] -> end -> num_lastedge);
		eraselast(midedge[i - 1] -> end, k);
	}
	for(i = 1; i < num_midedge; i ++)	{
		node = midedge[i] -> begin;
		if(node)	 {
			free_nodes(node);
		}
		free((void *) midedge[i]);
	}
	return(tmpedge);
}

READINTERVAL *insert_readcov(READINTERVAL *readcov, int readindex, int begin, int length, int offset)
{
	int	i, j, k, l;
	READINTERVAL	*readcov0, *readcov_new, *readcov_last;
	int	pos;

	readcov0 = readcov;
	readcov_last = (READINTERVAL *) NULL;
	while(readcov0)	{
/*	Check if some readcov overlapped with the new one	*/
		if(readcov0 -> eq_read == readindex && begin > readcov0 -> begin &&
		   begin <= readcov0 -> begin + readcov0 -> length - 1)	{
			readcov0 -> length = (begin + length - readcov0 -> begin);
			return(readcov);
		} else if(readcov0 -> eq_read > readindex)	{
			break;
		}
		readcov_last = readcov0;
		readcov0 = readcov0 -> next;
	}
	readcov_new = (READINTERVAL *) ckalloc(1 * sizeof(READINTERVAL));
	readcov_new -> eq_read = readindex;
	readcov_new -> begin = begin;
	readcov_new -> length = length;
	readcov_new -> offset = offset;
	readcov_new -> next = readcov0;
	readcov_new -> cov = 1;
	if(readcov_last)	{
		readcov_last -> next = readcov_new;
		return(readcov);
	} else	{
		return(readcov_new);
	}
}

int copyreadinterval(READINTERVAL *readinterval, READINTERVAL *readcov)
{
	int	i, n;

	i = 0;
	while(readcov)	{
		if(readcov -> cov == 1)	{
			readinterval[i] = *readcov;
			readinterval[i ++].next = (READINTERVAL *) NULL;
		}
		readcov = readcov -> next;
	}
	return(i);
}

void sortreadinterval(READINTERVAL *readinterval, int multip)
{
	int	n;

	qsort((void *) readinterval, multip, sizeof(READINTERVAL),  readintervalcompar);
}

int readintervalcompar(const void *va, const void *vb)
{
  READINTERVAL *a, *b;
  a = (READINTERVAL*) va;
  b = (READINTERVAL*) vb;
	if(a -> offset > b -> offset)	return(1);
	if(a -> offset == b -> offset && a -> length > b -> length)	return(1);
	if(a -> offset < b -> offset)	return(-1);
	return(0);
}

void sortreadinterval_index(READINTERVAL *readinterval, int multip)
{
	int	n;

	qsort((void *) readinterval, multip, sizeof(READINTERVAL), readintervalcompar_index);
}

int readintervalcompar_index(const void *a, const void *b)
{
	if(((READINTERVAL*)a) -> eq_read > ((READINTERVAL*)b) -> eq_read ||
	   ((READINTERVAL*)a) -> eq_read == ((READINTERVAL*)b) -> eq_read && ((READINTERVAL*)a) -> begin > ((READINTERVAL*)b) -> begin)	return(1);
	if(((READINTERVAL*)a) -> eq_read < ((READINTERVAL*)b) -> eq_read ||
	   ((READINTERVAL*)a) -> eq_read == ((READINTERVAL*)b) -> eq_read && ((READINTERVAL*)a) -> begin < ((READINTERVAL*)b) -> begin)	return(-1);
	return(0);
}

void sort_nodepos(NODES *node)
{
	int	i, j, k;
	READPOSITION *readposition;
	POSPOINT *tmppos;

	if (node->npos == 0) return;

	tmppos = (POSPOINT *) ckalloc(node -> npos * sizeof(POSPOINT));
	i = 0;
	readposition = node -> readposition;
	while(readposition)	{
		tmppos[i ++].readposition = readposition;
		readposition = readposition -> next;
	}
	if(i != node -> npos)	{
		printf("npos not equal %d %d\n", node -> npos, i);
		exit(0);
	}
	qsort((void *) tmppos, node -> npos, sizeof(POSPOINT), poscompar);
	node -> readposition = tmppos[0].readposition;
	for(i = 0; i < node -> npos - 1; i ++)	{
		tmppos[i].readposition -> next = tmppos[i + 1].readposition;
	}
	tmppos[i].readposition -> next = (READPOSITION *) NULL;
	free((void *) tmppos);
}

int poscompar(const void *a, const void *b)
{
	if(((POSPOINT*)a) -> readposition -> readindex > ((POSPOINT*)b) -> readposition -> readindex)	 return(1);
	if(((POSPOINT*)a) -> readposition -> readindex < ((POSPOINT*)b) -> readposition -> readindex)	 return(-1);
	return(0);
}
