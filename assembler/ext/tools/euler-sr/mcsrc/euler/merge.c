/***************************************************************************
 * Title:          merge.c
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

char findnode(INSERT *insert, NODES **node0, int num, int read, int pos);
void ins_node(READLIST **readlist, int read1, int read2, int pos1, int pos2, int endpos);
INSERT *insert_nodes(INSERT *insert, NODES **nodes, int num, int read, int pos);
INSERT *free_insert(INSERT *insert);
NODES *chk_merge_node(READLIST **readlist, int read1, int read2, int pos1, int pos2);
void insert_position(NODES *node, int read, int pos);
void merge(int num_seq, int *len_seq, ALIGN **eq_class, int num_class, READLIST **readlist);
NODES *combine_nodes(NODES *node1, NODES *node2);
NODES *free_nodes(NODES *node);
void update_link(NODES *node, int r1, int s1, int r2, int s2);

void merge(int num_seq, int *len_seq, ALIGN **eq_class, int num_class, READLIST **readlist)
{
  /****************************************************
    !!!!!!!  THIS CODE IS NOT USED !!!!!!!!!
  ***************************************************/
  printf("ERROR, running deprecated code\n");
  assert(0);

	int 	i, j, k, l, m, n;
	char	c1, c2;
	ALIGN	*align;
	int	read1, read2, read_rev1, read_rev2, pos1, pos2, pos_rev1, pos_rev2;
	NODES	*node, *node_rev, *node_next, *node1, *node2;
	INSERT	*insert;

	n = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			read_rev1 = reverse_read(read1, num_seq);
			read_rev2 = reverse_read(read2, num_seq);
			for(j = 0; j < align -> length - 1; j ++)	{
				pos1 = align -> pos[0][j];
				pos2 = align -> pos[1][j];
				pos_rev1 = len_seq[read1] - pos1 - 1;
				pos_rev2 = len_seq[read2] - pos2 - 1;
				while(pos1 < align -> pos[0][j + 1] && pos2 < align -> pos[1][j + 1])	{
					node = chk_merge_node(readlist, read1, read2, pos1, pos2);
					node_rev = chk_merge_node(readlist, read_rev1, read_rev2, pos_rev1, pos_rev2);
					node = readlist[read1][pos1].node;
					node -> bal_node = node_rev;
					if(node != node_rev)	{
						node_rev -> bal_node = node;
					}
					pos1 ++;
					pos2 ++;
					pos_rev1 --;
					pos_rev2 --;
				}
			}
			n ++;
			align = align -> next;
		}
		if(i % 500 == 1)	printf("..");
	}
	printf("# merged overlaps: %d\n", n);
}

NODES *chk_merge_node(READLIST **readlist, int read1, int read2, int pos1, int pos2)
{
	int	i, j, k, l;
	char	c;
	READPOSITION *readposition;
	INSERT	*insert, *insert0, *insert1;
	NODES	*node1, *node2, *node3, *new_node;
	/* 
	   Merge two nodes.  This requires merging the readpositions that are stored at each 
	   node, and updating the read positions to point to the new merged node.
	*/

	if(pos1 >= 0 && pos2 >= 0)	{
	  /* 
	     Each readlist .
	  */
		node1 = readlist[read1][pos1].node;
		node2 = readlist[read2][pos2].node;
		if(node1 == node2)	{
		  /*
		    The two node-pointers already point to the same node, but
		    they got there with different alignments. Increment the number of times
		    the node is linked to by 1.
		  */
			update_link(node1, read1, pos1, read2, pos2);
			return(node1);
		} else	{
/*	new_node == node1	*/
			new_node = combine_nodes(node1, node2);
			update_link(new_node, read1, pos1, read2, pos2);
			readposition = node2 -> readposition;
			/* Each read position (nucleotide in a read) corresponds to a node in the graph.
			   Before it correponded to 'node2'.  Now it needs to correspond to 'new_node'.
			   Update the list.
			*/
			while(readposition)	{
				read2 = readposition -> readindex;
				pos2 = readposition -> position;
/*	Link readlists to new node	*/
				if(pos2 >= 0)	{
					readlist[read2][pos2].node = new_node;
				}
				readposition = readposition -> next;
			}
			free((void *) node2);
			return(new_node);
		}
	} else	{
		printf("Negative readpositions: %d %d\n", pos1, pos2);
		exit(-1);
	}
}

void ins_node(READLIST **readlist, int read1, int read2, int pos1, int pos2, int endpos)
{
	int	i, j, k, l, m, n;
	char	c;
	NODES	**nodes;

	if(pos2 < 0)	{
		l = endpos - pos1 + 1;
		nodes = (NODES **) ckalloc(l * sizeof(NODES *));
		for(i = pos1; i <= endpos; i ++)	{
			nodes[i - pos1] = readlist[read1][i].node;
		}
		c = findnode(readlist[read2][-pos2-1].insert, nodes, l, read2, pos2);
		if(c)	{
			readlist[read2][-pos2-1].insert = insert_nodes(readlist[read2][-pos2-1].insert, nodes, l,
				read2, pos2);
		}
		for(i = 0; i < l; i ++)	{
			insert_position(nodes[i], read2, pos2);
			update_link(nodes[i], read1, pos1 + i, read2, pos2);
		}
		free((void **) nodes);
	} else if(pos1 < 0)	{
		l = endpos - pos2 + 1;
		nodes = (NODES **) ckalloc(l * sizeof(NODES *));
		for(i = pos2; i <= endpos; i ++)	{
			nodes[i - pos2] = readlist[read2][i].node;
		}
		c = findnode(readlist[read1][-pos1-1].insert, nodes, l, read1, pos1);
		if(c)	{
			readlist[read1][-pos1-1].insert = insert_nodes(readlist[read1][-pos1-1].insert, nodes, l,
				read1, pos1);
		}
		for(i = 0; i < l; i ++)	{
			insert_position(nodes[i], read1, pos1);
			update_link(nodes[i], read1, pos1, read2, pos2 + i);
		}
		free((void **) nodes);
	} else	{
		printf("Both positive readpositions: %d %d\n", pos1, pos2);
		exit(-1);
	}
}

INSERT *insert_nodes(INSERT *insert, NODES **nodes, int num, int read, int pos)
{
	int	i;
	INSERT	*insert0;

	insert0 = (INSERT *) ckalloc(1 * sizeof(INSERT));
	insert0 -> node = (NODES **) ckalloc(num * sizeof(NODES *));
	insert0 -> num_nodes = num;
	for(i = 0; i < num; i ++)	{
		insert0 -> node[i] = nodes[i];
	}
	insert0 -> next = insert;
	return(insert0);
}

char findnode(INSERT *insert, NODES **node0, int num, int read, int pos)
{
	int	i, k;
	INSERT	*insert1;

	insert1 = insert;
	while(insert1)	{
		if(insert1 -> num_nodes == num)	{
			for(i = 0; i < insert -> num_nodes; i ++)	{
				if(node0[i] != insert1 -> node[i])	{
					break;
				}
			}
			if(i == insert -> num_nodes)	{
				return(0);
			}
		}
		insert1 = insert1 -> next;
	}
	return(1);
}

void insert_position(NODES *node, int read, int pos)
{
	READPOSITION *readposition, *readposition1;

	readposition1 = node -> readposition;
	while(readposition1)	{
		if(readposition1 -> readindex == read && readposition1 -> position == pos)	{
			return;
		}
		readposition1 = readposition1 -> next;
	}
	readposition = (READPOSITION *) ckalloc(1 * sizeof(READPOSITION));
	readposition -> readindex = read;
	readposition -> position = pos;
	readposition -> next = node -> readposition;
	node -> readposition = readposition;
	node -> npos ++;
}

void update_link(NODES *node, int r1, int s1, int r2, int s2)
{
	int	i, j, k, l, t1, t2;
	READPOSITION	*readposition;
	node -> nlinks ++;
}

NODES *free_nodes(NODES *node)
{
	READPOSITION	*readposition;

	while(node -> readposition)	{
		readposition = node -> readposition -> next;
		free((void *) node -> readposition);
		node -> readposition = readposition;
	}
	if(node -> path_pos)	free((void **) node -> path_pos);
	if(node -> path_index)	free((void **) node -> path_index);
	if(node -> lastedge)	free((void **) node -> lastedge);
	if(node -> nextedge)	free((void **) node -> nextedge);
	free((void *) node);
	return((NODES *) NULL);
}

NODES *combine_nodes(NODES *node1, NODES *node2)
{
	int	i, k, j, l, n;
	READPOSITION *pos, *pos1, *pos2, *pos3;
	INSERT	*insert;
	NODES	*new_node;
	/* 
	   Merge two nodes together.  The new node has a number of links
	   equal to the sum of the links in node1 and node2 (nlinks are the number of different 
	   vertices joined together, the 'size' of a vertex), and the readpositions
	   of the two nodes are joined.
	*/

	new_node = node1;
	new_node -> npos += node2 -> npos;
	new_node -> nlinks += node2 -> nlinks;
	pos = new_node -> readposition;
	/* find the end of the list  of read positions.*/
	while(pos -> next)	{
		pos = pos -> next;
	}
	pos -> next = node2 -> readposition;
	return(new_node);
}

INSERT *free_insert(INSERT *insert)
{
	INSERT	*insert1;

	insert1 = insert -> next;
	free((void **) insert -> node);
	free((void *) insert);
	return(insert1);
}
