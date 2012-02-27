/***************************************************************************
 * Title:          merge_gap.c
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

void merge_gap(char **src_seq, int num_seq, int *len_seq, ALIGN **eq_class, READLIST **readlist, int gap_k);

void merge_gap(char **src_seq, int num_seq, int *len_seq, ALIGN **eq_class, READLIST **readlist, int gap_k)
{
  /*
    src_seq - a table of reads (the actual ACTG sequence, translated to lower case.
    len_seq - list of lengths of arrays.
    
    eq_class - a list of lists of alignments, one list per read.
    readlist - a list of node-arrays, one node-array per read.
    gap_k    - not sure.  it's 0 in default runs of over-repeat-new

    Result:
    The alignments stored in eq_class are used to link nodes in 'readlist' together.
  */


	int 	i, j, k, l, m, n;
	char	c1, c2;
	ALIGN	*align;
	int	read1, read2, read_rev1, read_rev2, pos1, pos2, pos_rev1, pos_rev2;
	int	pnext1, pnext2;
	NODES	*node, *node_rev, *node_next, *node1, *node2;
	INSERT	*insert;

	n = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			/*
			  Get the indices of the reverse complements of a read.
			*/
			read_rev1 = reverse_read(read1, num_seq);
			read_rev2 = reverse_read(read2, num_seq);
			
			for(j = 0; j < align -> length - 1; j ++)	{
			  /* END_MERGE is a flag for whether or not the ends of reads are trusted.
			     in over-repeat-new, END_MERGE=1, so the ends are always trusted.
			  */
				if(j == 0 && END_MERGE)	{
				  /* pos1 is in read1 is aligned with pos2 in read 2 */
					pos1 = align -> pos[0][j];
					pos2 = align -> pos[1][j];
				} else	{
					pos1 = align -> pos[0][j] + gap_k;
					pos2 = align -> pos[1][j] + gap_k;
				}
				pos_rev1 = len_seq[read1] - pos1 - 1;
				pos_rev2 = len_seq[read2] - pos2 - 1;
				if(j == align -> length - 2 && END_MERGE)	{
					pnext1 = align -> pos[0][j + 1];
					pnext2 = align -> pos[1][j + 1];
				} else	{
					pnext1 = align -> pos[0][j + 1] - gap_k;
					pnext2 = align -> pos[1][j + 1] - gap_k;
				}
				/* 
				   read1 is aligned with an ugapped alignment from pos1 ... pnext1 
				   to read2 pos2 .. pnext2
				   
				   Glue together all nodes that corresponded to these positions.
				*/
				while(pos1 < pnext1 && pos2 < pnext2)	{
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
		if(i % 500 == 0)	{
			printf("..");
			fflush(stdout);
		}
	}
	printf("# merged overlaps: %d\n", n);
}
