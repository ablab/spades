/***************************************************************************
 * Title:          initialize.c
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

void rem_position(NODES *node, int reads, int pos);
void initialize(READLIST **readlist, int *len_seq, int num_seq);
void makenode(READLIST **readlist, int i, int j, int *len_seq, int num_seq);
void separatenode(NODES *node, READLIST **readlist, int *len_seq, int num_seq);
void separatenode0(NODES *node, READLIST **readlist, int *len_seq, int num_seq, int mthresh);
int add_prenode(NODES **prenode, int nnode, int *mnode, NODES *node);
void 	splitnode(NODES *node, NODES *prenode, READLIST **readlist, int *len_seq, int num_seq);

void initialize(READLIST **readlist, int *len_seq, int num_seq)
{
	int	i, j, k, l;
	NODES	*vertex1, *vertex2;

	for(i = 0; i < num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			makenode(readlist, i, j, len_seq, num_seq);
		}
	}
}

void makenode(READLIST **readlist, int i, int j, int *len_seq, int num_seq) 
{
	int	k;
	NODES	*vertex1, *vertex2;

	vertex1 = (NODES *) ckalloc(1 * sizeof(NODES));
	vertex1 -> npos = 1;
	vertex1 -> readposition = (READPOSITION *) ckalloc(1 * sizeof(READPOSITION));
	vertex1 -> readposition[0].readindex = i;
	vertex1 -> readposition[0].position = j;
	readlist[i][j].node = vertex1;
	vertex2 = (NODES *) ckalloc(1 * sizeof(NODES));
	vertex2 -> npos = 1;
	vertex2 -> readposition = (READPOSITION *) ckalloc(1 * sizeof(READPOSITION));
	k = reverse_read(i, num_seq);
	vertex2 -> readposition[0].readindex = k;
	vertex2 -> readposition[0].position = len_seq[i] - j - 1;
	readlist[k][len_seq[i] - j - 1].node = vertex2;
	vertex1 -> bal_node = vertex2;
	vertex2 -> bal_node = vertex1;
}

void separatenode(NODES *node, READLIST **readlist, int *len_seq, int num_seq)
{
	int	i, j, k, l;
	READPOSITION	*readposition;

	readposition = node -> readposition;
	while(readposition)	{
		if(readlist[readposition -> readindex][readposition -> position].node == node)	{
			if(readposition -> readindex < num_seq)	{
				makenode(readlist, readposition -> readindex, readposition -> position, len_seq, num_seq);
			} else	{
				makenode(readlist, readposition -> readindex - num_seq, 
					 len_seq[readposition -> readindex] - 1 - readposition -> position, len_seq, num_seq);
			}
		}
		readposition = readposition -> next;
	}
	node -> npos = 0;
}

void separatenode0(NODES *node, READLIST **readlist, int *len_seq, int num_seq, int mthresh)
{
	int	i, j, k, l;
	READPOSITION	*readposition;
	int	nnode, *mnode;
	NODES	**prenode, *bal_node;

	nnode = 0;
	prenode = (NODES **) ckalloc(node -> npos * sizeof(NODES *));
	mnode = (int *) ckalloc(node -> npos * sizeof(int));
	readposition = node -> readposition;
	while(readposition)	{
		if(readposition -> position > 0 && readposition -> position < len_seq[readposition -> readindex] - 1)	{
			if(readlist[readposition -> readindex][readposition -> position - 1].node -> num_path <= 1) {
				nnode = add_prenode(prenode, nnode, mnode, readlist[readposition -> readindex][readposition -> position - 1].node);
			}
		}
		readposition = readposition -> next;
	}
	for(i = 0; i < nnode; i ++)	{
		if(node -> npos > 0 && mnode[i] >= mthresh)	{
			splitnode(node, prenode[i], readlist, len_seq, num_seq);
		}
	}
	free((void **) prenode);
	free((void *) mnode);
}

void separatenode1(NODES *node, READLIST **readlist, int *len_seq, int num_seq, int mthresh)
{
	int	i, j, k, l;
	READPOSITION	*readposition;
	int	nnode, *mnode;
	NODES	**prenode, *bal_node;

	nnode = 0;
	prenode = (NODES **) ckalloc(node -> npos * sizeof(NODES *));
	mnode = (int *) ckalloc(node -> npos * sizeof(int));
	readposition = node -> readposition;
	while(readposition)	{
		if(readposition -> position > 0 && readposition -> position < len_seq[readposition -> readindex] - 1)	{
			if(readlist[readposition -> readindex][readposition -> position - 1].node -> num_path <= 1) {
				nnode = add_prenode(prenode, nnode, mnode, readlist[readposition -> readindex][readposition -> position - 1].node);
			}
		}
		readposition = readposition -> next;
	}
	for(i = 0; i < nnode; i ++)	{
		if(node -> npos > 0 && mnode[i] >= mthresh)	{
			splitnode(node, prenode[i], readlist, len_seq, num_seq);
		}
	}
	free((void **) prenode);
	free((void *) mnode);
}

void 	splitnode(NODES *node, NODES *prenode, READLIST **readlist, int *len_seq, int num_seq)
{
	int	i, j, k, l, m, n;
	NODES	*newnode, *bal_node;
	READPOSITION	*readposition, *pos_last, *pos_next;

	bal_node = node -> bal_node;
	newnode = (NODES *) NULL;
	pos_last = (READPOSITION *) NULL;
	readposition = node -> readposition;
	while(readposition)	{
		pos_next = readposition -> next;
		m = reverse_read(readposition -> readindex, num_seq);
		l = len_seq[readposition -> readindex] - 1 - readposition -> position;
		if(readposition -> position > 0 && l > 0)	{
			if(readlist[readposition -> readindex][readposition -> position - 1].node == prenode)	{
				if(!newnode)	{
					makenode(readlist, readposition -> readindex, readposition -> position, len_seq, num_seq);
					newnode = readlist[readposition -> readindex][readposition -> position].node;
				} else	{
					insert_position(newnode, readposition -> readindex, readposition -> position);
					readlist[readposition -> readindex][readposition -> position].node = newnode;
					insert_position(newnode -> bal_node, m, l);
					readlist[m][l].node = newnode -> bal_node;
				}
				if(pos_last)	{
					pos_last -> next = pos_next;
				} else	{
					node -> readposition = pos_next;
				}
				if(bal_node == node)	{
					if(pos_next && pos_next -> readindex == m && pos_next -> position == l)	{
						pos_next = pos_next -> next;
					}
				}
				rem_position(bal_node, m, l);
				node -> npos --;
				free((void *) readposition);
			} else	{
				pos_last = readposition;
			}
		} else	{
			pos_last = readposition;
		}
		readposition = pos_next;
	}
}

void rem_position(NODES *node, int reads, int pos)
{
	READPOSITION *readposition, *pos_last, *pos_next;

	pos_last = (READPOSITION *) NULL;
	readposition = node -> readposition;
	while(readposition)	{
		pos_next = readposition -> next;
		if(readposition -> readindex == reads && readposition -> position == pos)	{
			if(pos_last)	{
				pos_last -> next = pos_next;
			} else	{
				node -> readposition = pos_next;
			}
			node -> npos --;
			free((void *) readposition);
			return;
		} else	{
			pos_last = readposition;
		}
		readposition = pos_next;
	}
	printf("Position not found: node %d read %d readposition %d.\n", node, reads, pos);
	printf("node %d %d\n", node, node -> npos);
	readposition = node -> readposition;
	while(readposition)	{
		printf("readposition %d %d\n", readposition -> readindex, readposition -> position);
		pos_next = readposition -> next;
		readposition = pos_next;
	}
	exit(-1);
}

int add_prenode(NODES **prenode, int nnode, int *mnode, NODES *node)
{
	int	i;

	for(i = 0; i < nnode; i ++)	{
		if(prenode[i] == node)	{
			mnode[i] ++;
			return(nnode);;
		}
	}
	mnode[nnode] = 1;
	prenode[nnode ++] = node;
	return(nnode);
}
