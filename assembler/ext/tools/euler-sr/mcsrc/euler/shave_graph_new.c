/***************************************************************************
 * Title:          shave_graph_new.c
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

void insert_chim(READTABLE *RT, int index);
int shave_graph_new(NODES **vertex, int num_vertex, READTABLE *RT, int MIN_LENGTH, int MIN_MULTIP);

int shave_graph_new(NODES **vertex, int num_vertex, READTABLE *RT, int MIN_LENGTH, int MIN_MULTIP)
{
	int	i, j, k, l, n, m, n1, reads;
	int	tot_edge;
	int	nsh, nbul, nch;
	int	maxmlt, maxl, maxk, multip;
	int	true_multip;
	NODES	*begin, *end, *bal_node;
	EDGE	*edge, *edge1, *edge2, *bal_edge, *bal_edge1, *bal_edge2;

	nch = nsh = 0;
	do	{
		n1 = 0;

/*	Remove single end edges */

		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_lastedge != 0)	continue;
			/*
			  vertex[i] is a source.
			*/
			j = 0;
			while(j < vertex[i] -> num_nextedge)	{
			  edge = vertex[i] -> nextedge[j];
				if(edge -> end -> num_lastedge != 1 && edge -> multip <= LOW_COV)	{
				  /* Edge: vertex[i]->nextedge[j] is short, remove it*/
					bal_edge = edge -> bal_edge;
					erasedge(edge);
					if(edge != bal_edge)	{
						erasedge(bal_edge);
						n1 ++;
						nsh ++;
					}
					n1 ++;
					nsh ++;
				} else	{
					j ++;
				}
			}
		}
		num_vertex = merge_graph(vertex, num_vertex, NULL, 0);

/*	Erase chimeric edges	*/

		for(i = 0; i < num_vertex; i ++)	{
			j = 0;
			while(j < vertex[i] -> num_nextedge)	{
				edge = vertex[i] -> nextedge[j];
				multip = edge -> multip;
				if(edge -> length < ChimericTerm && vertex[i] -> num_nextedge > 1 &&
				   edge -> end -> num_lastedge > 1 && multip <= MIN_MULTIP)	{
					bal_edge = edge -> bal_edge;
					for(m = 0; m < edge -> multip; m ++)	{
						reads = edge -> readinterval[m].eq_read;
						insert_chim(RT, reads);
						nch ++;
						n1 ++;
					}
					erasedge(edge);
					if(edge != bal_edge)	{
						for(m = 0; m < bal_edge -> multip; m ++)	{
							reads = bal_edge -> readinterval[m].eq_read;
							insert_chim(RT, reads);
							nch ++;
							n1 ++;
						}
						erasedge(bal_edge);
					}
				} else	{
					j ++;
				}
			}
		}
		num_vertex = merge_graph(vertex, num_vertex, NULL, 0);

/*	Remove end edges */

		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_lastedge != 0)	continue;
			j = 0;
			while(j < vertex[i] -> num_nextedge)	{
				edge = vertex[i] -> nextedge[j];
				multip = realmultip1(edge, EndLength);
				if(edge -> end -> num_lastedge != 1 && multip <= LOW_COV)	{
					bal_edge = edge -> bal_edge;
					erasedge(edge);
					if(edge != bal_edge)	{
						erasedge(bal_edge);
						n1 ++;
						nsh ++;
					}
					n1 ++;
					nsh ++;
				} else	{
					j ++;
				}
			}
		}
		num_vertex = merge_graph(vertex, num_vertex, NULL, 0);
	} while(n1 > 0);

	tot_edge = 0;
	for(i = 0; i < num_vertex; i ++)	{
		tot_edge += vertex[i] -> num_nextedge;
	}

	return(num_vertex);
}

void insert_chim(READTABLE *RT, int index)
{
	int	i, *chim;

	if(!(RT -> chim))	{
		RT -> chim = (int *) ckalloc(100 * sizeof(int));
		RT -> num_chim_alloc = 100;
		RT -> num_chim = 0;
	} else if(RT -> num_chim >= RT -> num_chim_alloc)	{
		chim = (int *) ckalloc(RT -> num_chim_alloc * sizeof(int));
		for(i = 0; i < RT -> num_chim; i ++)	{
			chim[i] = RT -> chim[i];
		}
		free((void *) RT -> chim);
		RT -> num_chim_alloc += 100;
		RT -> chim = (int *) ckalloc(RT -> num_chim_alloc * sizeof(int));
		for(i = 0; i < RT -> num_chim; i ++)	{
			RT -> chim[i] = chim[i];
		}
		free((void *) chim);
	}
	RT -> chim[RT -> num_chim ++] = index;
}
