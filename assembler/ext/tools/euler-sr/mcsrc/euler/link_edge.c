/***************************************************************************
 * Title:          link_edge.c
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

int link_edge(int num_seq, int *len_seq, READLIST **readlist, NODES **nodes);

int link_edge(int num_seq, int *len_seq, READLIST **readlist, NODES **nodes)
{
  /* 
     Use adjacencies in the readlist to define edges between nodes.
   */

	int	i, j, k, l, m, n;
	EDGE	*edge1, *edge2;
	NODES	*node, *node_next, *node1, *node2;
	int len_seq_i;
	
	n = m = 0;
	for(i = 0; i < num_seq; i ++)	{
	        len_seq_i = len_seq[i];
		if (len_seq_i == 0)
		  continue;
		node = readlist[i][0].node;
		for(j = 0; j < len_seq_i - 1; j ++)	{
			node_next = readlist[i][j + 1].node;
			/* 
			   node, and node_next are adjacent in read, readlist[i], connect
			   them with an edge.  Add the position of the read to the list
			   of intervals of the read.
			*/
			edge1 = insert_edge(node, node_next, i, j, j + 1);
			/* 
			   Link the vertices corersponding to the reverse complement.
			*/
			node1 = readlist[i + num_seq][len_seq_i - j - 2].node;
			node2 = readlist[i + num_seq][len_seq_i - j - 1].node;
			
			edge2 = insert_edge(node1, node2, i + num_seq,
					    len_seq_i - j - 2, len_seq_i - j - 1);

			/* 
			   If both edge1 and edge2 are new, then the balance should be assigned now.
			   Otherwise, edge1 is not the balanced edge of edge2, neither edge1 nor 
			   edge2 should be the same.  On the other hand, if edge1 was created,
			   but edge2 was not, that mean that when trying to create edge2 it was found to already
			   exist, therefore edge1 is the balance of itself.
			*/
			if(edge1 && edge2)	{
				edge1 -> bal_edge = edge2;
				edge2 -> bal_edge = edge1;
				m += 2;
			} else if(edge1)	{
				edge1 -> bal_edge = edge1;
				m ++;
			}

			/* Nodes stores a list of unique nodes in the graph. Assign 'nodes[n] to this 
			   node if it is new and it hasn't been visited (by any other aligned  position).
			*/
			if(node -> visit == 0)	{
				nodes[n ++] = node;
				node -> visit = 1;
			}
			/* 
			   Do the same with the balanced node.
			*/
			if(node -> bal_node -> visit == 0)	{
				nodes[n ++] = node -> bal_node;
				node -> bal_node -> visit = 1;
			}
			node = node_next;
		}
		/* 
		   Handle the last node in the list.
		*/
		if(node -> visit == 0)	{
			nodes[n ++] = node;
			node -> visit = 1;
		}
		if(node -> bal_node -> visit == 0)	{
			nodes[n ++] = node -> bal_node;
			node -> bal_node -> visit = 1;
		}
	}
	printf("# nodes: %d, # edges: %d.\n", n, m);
	return(m);
}
