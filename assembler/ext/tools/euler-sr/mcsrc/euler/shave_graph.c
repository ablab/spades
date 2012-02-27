/***************************************************************************
 * Title:          shave_graph.c
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

#define MIN_LENGTH 1000
#define MIN_MULTIP 10


int realmultip1(EDGE *edge, int l);
int realmultip2(EDGE *edge, int l);
int merge_graph(NODES **vertex, int num_vertex, PATH *paths, int numPaths);
int shave_graph(NODES **vertex, int num_vertex, int *chim, int *num_chim, PATH *paths, int numPaths);
EDGE *merge_vertex(NODES *vertex, PATH *paths, int numPaths);
int chklist(int *list, int n, int index);
EDGE *new_edge(NODES *vertex, NODES *begin, NODES *end, EDGE *edge1, EDGE *edge2, 
							 int *beginlist, int *endlist,
							 int nbegin, int nend);
void combine_readinterval(EDGE *edge1, EDGE *edge2, EDGE *newedge, int *beginlist, int *endlist, int nbegin, int nend, PATH *paths, int numPaths);
int FindIntervalEndingAtPos(READINTERVAL *readinterval, int num, int readindex, int readposition);
int FindIntervalBeginningAtPos(READINTERVAL *readinterval, int num, int readindex, int readposition);

int shave_graph(NODES **vertex, int num_vertex, int *chim, int *num_chim, PATH *paths, int numPaths)
{
	int	i, j, k, l, n, m, n1, reads;
	int	tot_edge;
	int	nsh, nbul, nch;
	int	maxmlt, maxl, maxk, multip;
	int	true_multip;
	NODES	*begin, *end, *bal_node;
	EDGE	*edge, *edge1, *edge2, *bal_edge, *bal_edge1, *bal_edge2;

/*	Remove bulges and shave edges linking sources & sinks	*/

	nsh = nbul = nch = 0;
	do	{
		m = 0;
		for(i = 0; i < num_vertex; i ++)	{
			for(k = 0; k < vertex[i] -> num_nextedge; k ++)	{
				edge1 = vertex[i] -> nextedge[k];
				bal_edge1 = edge1 -> bal_edge;
				j = k + 1;
				while(j < vertex[i] -> num_nextedge)	{
					edge2 = vertex[i] -> nextedge[j];
					bal_edge2 = edge2 -> bal_edge;
					if(edge2 -> end == edge1 -> end && edge1 != bal_edge2 &&
					   abs(edge1 -> length - edge2 -> length) < MIN_INT)	{
						movereadinterval(edge1, edge2);
						erasedge(edge2);
						if(bal_edge2 != edge2)	{
							movereadinterval(bal_edge1, bal_edge2);
							erasedge(bal_edge2);
							m ++;
						}
						m ++;
					} else	{
						j ++;
					}
				}
			}
		}
		nbul += m;

/*	Remove end edges	*/

		n = 0;
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_lastedge != 0)	continue;
			j = 0;
			while(j < vertex[i] -> num_nextedge)	{
				edge = vertex[i] -> nextedge[j];
				multip = realmultip1(edge, 100);
				if(edge -> end -> num_lastedge > 1 && multip <= LOW_COV)	{
					bal_edge = edge -> bal_edge;
					erasedge(edge);
					if(edge != bal_edge)	{
						erasedge(bal_edge);
						n ++;
					}
					n ++;
				} else	{
					j ++;
				}
			}
		}
		nsh += n;
		num_vertex = merge_graph(vertex, num_vertex, paths, numPaths);
	} while(n > 0 || m > 0);

/*	Erase chimeric edges	*/

	do	{
		n1 = 0;
		for(i = 0; i < num_vertex; i ++)	{
			j = 0;
			while(j < vertex[i] -> num_nextedge)	{
				edge = vertex[i] -> nextedge[j];
				if(vertex[i] -> num_nextedge > 1 &&
				   edge -> end -> num_lastedge > 1 && edge -> length < MIN_LENGTH &&
				   edge -> multip < MIN_MULTIP)	{
					multip = realmultip2(edge, overlaplen);
					if(multip <= MID_COV)	{
						bal_edge = edge -> bal_edge;
						for(m = 0; m < edge -> multip; m ++)	{
							reads = edge -> readinterval[m].eq_read;
							chim[nch ++] = reads;
							n1 ++;
						}
						erasedge(edge);
						if(edge != bal_edge)	{
							for(m = 0; m < bal_edge -> multip; m ++)	{
								reads = bal_edge -> readinterval[m].eq_read;
								chim[nch ++] = reads;
								n1 ++;
							}
							erasedge(bal_edge);
						}
					} else	{
						j ++;
					}
				} else	{
					j ++;
				}
			}
		}
		num_vertex = merge_graph(vertex, num_vertex, paths, numPaths);

		m = 0;
		for(i = 0; i < num_vertex; i ++)	{
			for(k = 0; k < vertex[i] -> num_nextedge; k ++)	{
				edge1 = vertex[i] -> nextedge[k];
				bal_edge1 = edge1 -> bal_edge;
				j = k + 1;
				while(j < vertex[i] -> num_nextedge)	{
					edge2 = vertex[i] -> nextedge[j];
					bal_edge2 = edge2 -> bal_edge;
					if(edge2 -> end == edge1 -> end && edge1 != bal_edge2 &&
					   abs(edge1 -> length - edge2 -> length) < MIN_INT)	{
						movereadinterval(edge1, edge2);
						erasedge(edge2);
						if(bal_edge2 != edge2)	{
							movereadinterval(bal_edge1, bal_edge2);
							erasedge(bal_edge2);
							m ++;
						}
						m ++;
					} else	{
						j ++;
					}
				}
			}
		}
		nbul += m;

/*	Remove end edges	*/

		n = 0;
		for(i = 0; i < num_vertex; i ++)	{
			if(vertex[i] -> num_lastedge != 0)	continue;
			j = 0;
			while(j < vertex[i] -> num_nextedge)	{
				edge = vertex[i] -> nextedge[j];
				multip = realmultip1(edge, END_LEG);
				if(edge -> end -> num_lastedge > 1 && multip <= LOW_COV)	{
					bal_edge = edge -> bal_edge;
					erasedge(edge);
					if(edge != bal_edge)	{
						erasedge(bal_edge);
						n ++;
					}
					n ++;
				} else	{
					j ++;
				}
			}
		}
		nsh += n;
		num_vertex = merge_graph(vertex, num_vertex, paths, numPaths);

		for(i = 0; i < num_vertex; i ++)	{
			bal_node = vertex[i] -> bal_node;
			if(vertex[i] -> num_lastedge != bal_node -> num_nextedge ||
			   vertex[i] -> num_nextedge != bal_node -> num_lastedge)	{
				/*				printf("pos2 vertex %d(%d-%d) bal_node %d(%d-%d)\n",
					vertex[i], vertex[i] -> num_lastedge, vertex[i] -> num_nextedge,
					bal_node, bal_node -> num_lastedge, bal_node -> num_nextedge);*/
				getchar();
			}
		}
	} while(n1 > 0);

	tot_edge = 0;
	for(i = 0; i < num_vertex; i ++)	{
		tot_edge += vertex[i] -> num_nextedge;
	}
	printf("%d short end edges, %d bulges and %d chimeric edges shaved, %d vertices %d edges left.\n",
		nsh, nbul, nch, num_vertex, tot_edge);

	*num_chim = nch;
	return(num_vertex);
}

int realmultip1(EDGE *edge, int l)
{
	int	i, j, k, k1, k2;

	if(edge -> length < l)	return(0);
	k = 0;
	for(i = 0; i < edge -> multip; i ++)	{
		if(edge -> readinterval[i].length >= l)	{
			k ++;
		}
	}
	return(k);
}

int realmultip2(EDGE *edge, int l)
{
	int	i, j, k1, k2, k;

	k = 0;
	for(i = 0; i < edge -> multip; i ++)	{
		if(edge -> readinterval[i].length >= min(l, edge -> length))	{
			k ++;
		}
	}
	return(k);
}


int merge_graph(NODES **vertex, int num_vertex, PATH *paths, int numPaths)
{
	int	i, j, k, l;
	int	num_vertex_old;
	NODES	*bal_node;
	READPOSITION	*readposition;
	EDGE	*edge1, *edge2;

	do	{
		num_vertex_old = num_vertex;
		i = 0;
		while(i < num_vertex)	{
			if(vertex[i] -> num_nextedge == 0 && vertex[i] -> num_lastedge == 0)	{
				free_nodes(vertex[i]);
				for(j = i; j < num_vertex - 1; j ++)	{
					vertex[j] = vertex[j + 1];
				}
				num_vertex --;
			} else if(vertex[i] -> num_nextedge == 1 && vertex[i] -> num_lastedge == 1 &&
				  vertex[i] -> nextedge[0] != vertex[i] -> lastedge[0]) {
				if (vertex[i]->num_nextedge > 0)
					bal_node = vertex[i]->nextedge[0]->bal_edge->end;
				else if (vertex[i]->num_lastedge > 0)
					bal_node = vertex[i]->lastedge[0]->bal_edge->begin;
						/*				bal_node = vertex[i] -> bal_node;*/
				edge1 = merge_vertex(vertex[i], paths, numPaths);
				if(bal_node != vertex[i])	{
					if(bal_node -> num_lastedge == 1 && bal_node -> num_nextedge == 1)	{
						if(bal_node == edge1 -> end || bal_node == edge1 -> begin)	{
							edge2 = merge_vertex(bal_node, paths, numPaths);
							edge2 -> bal_edge = edge2;
						} else	{
							edge2 = merge_vertex(bal_node, paths, numPaths);
							edge1 -> bal_edge = edge2;
							edge2 -> bal_edge = edge1;
						}
					} else	{
						printf("pos3 vertex %d(%d-%d) bal_node %d(%d-%d)\n",
							vertex[i], vertex[i] -> num_lastedge, vertex[i] -> num_nextedge,
							bal_node, bal_node -> num_lastedge, bal_node -> num_nextedge);
						exit(0);
					}
				} else	{
					edge1 -> bal_edge = edge1;
				}
			} else	{
				i ++;
			}
		}
	} while(num_vertex_old > num_vertex);
	return(num_vertex);
}

EDGE *merge_vertex(NODES *vertex, PATH *paths, int numPaths)
{
	int	i, j, k, l, m, n, k1, k2, n1, n2, c, q, p, num;
	NODES	*begin, *end, *b0, *e0;
	int	*beginlist, *endlist;
	int	*del_path;
	EDGE	*lastedge, *nextedge, *newedge, *new1edge;
	double	v, r;

	begin = vertex -> lastedge[0] -> begin;
	end = vertex -> nextedge[0] -> end;
	lastedge = vertex -> lastedge[0];
	nextedge = vertex -> nextedge[0];
	if(lastedge == nextedge)	{
		return(lastedge);
	}
	newedge = new_edge(vertex, begin, end, lastedge, nextedge, beginlist, endlist, 0, 0);
	printf("merged : %d (%d) %d (%d) via %d into %d %d \n", 
				 lastedge, lastedge->index, nextedge, nextedge->index,
				 vertex, newedge, newedge->index);
	erasedge(lastedge);
	erasedge(nextedge);
	return(newedge);
}

EDGE *new_edge(NODES *vertex, NODES *begin, NODES *end, 
							 EDGE *edge1, EDGE *edge2, 
							 int *beginlist, int *endlist,
							 int nbegin, int nend)
{
	
	/* Create a new edge, 'edge', and add it to the nextlist from 'begin'
		 and the last list from 'end'.
	*/
	int	i, j, k, l, m;
	EDGE	**tmpedge, *edge;

	edge = (EDGE *) ckalloc(1 * sizeof(EDGE));
	edge->deleted = 0;
	edge -> begin = begin;
	edge -> end = end;
	edge -> length = edge1 -> length + edge2 -> length - VERTEX_SIZE;
	/* Make read intervals that pass from edge1 to edge2 have one less segment.*/
	combine_readinterval(edge1, edge2, edge, beginlist, endlist, nbegin, nend);
	edge->index = edge2->index;
	/*	printf("creating new edge with index: %d\n", edge->index);*/
	if (begin->num_nextedge > 0) 
	  tmpedge = (EDGE **) ckalloc(begin->num_nextedge * sizeof(EDGE *));
	else
	  tmpedge = NULL;

	/* add  newedge to begin*/
	for(i = 0; i < begin -> num_nextedge; i ++)	{
	  tmpedge[i] = begin -> nextedge[i];
	}
	if (begin->nextedge != NULL)
	  free((void **) begin -> nextedge);
	  
	begin -> nextedge = (EDGE **) ckalloc((begin -> num_nextedge + 1) * sizeof(EDGE *));
	for(i = 0; i < begin -> num_nextedge; i ++)	{
	  begin -> nextedge[i] = tmpedge[i];
	}
	begin -> nextedge[i] = edge;
	begin -> num_nextedge ++;
	
	if (tmpedge != NULL) free(tmpedge);
	/* add newedge to end*/
	if (end->num_lastedge > 0) 
	  tmpedge = (EDGE**) ckalloc(end->num_lastedge * sizeof(EDGE*));
	else
	  tmpedge = NULL;

	for(i = 0; i < end -> num_lastedge; i ++)	{
		tmpedge[i] = end -> lastedge[i];
	}
	if (end->lastedge != NULL)
	  free((void **) end -> lastedge);
	
	end -> lastedge = (EDGE **) ckalloc((end -> num_lastedge + 1) * sizeof(EDGE *));
	for(i = 0; i < end -> num_lastedge; i ++)	{
		end -> lastedge[i] = tmpedge[i];
	}
	end -> lastedge[i] = edge;
	end -> num_lastedge ++;
	
	if (tmpedge != NULL)
	  free((void **) tmpedge);

	return(edge);
}

void combine_readinterval(EDGE *edge1, EDGE *edge2, EDGE *newedge, 
													int *beginlist, int *endlist, int nbegin, int nend)
{
	int	i, j, k, l, m, n, c, num1, t, q;
	int	reads;
	char	*mark1, *mark2;
	double	rm;
	int readIndex, pathLength;
	NODES	*vertex0;
	READINTERVAL *readinterval;
	EDGE	*edge01, *edge02, *edge3;

	vertex0 = edge1 -> end;

	l = edge1 -> length - VERTEX_SIZE;
	k = edge1 -> multip + edge2 -> multip;
	mark1 = (char *) ckalloc(k * sizeof(char));
	mark2 = (char *) ckalloc(k * sizeof(char));
	//	readinterval = NULL; //(READINTERVAL *) ckalloc(k * sizeof(READINTERVAL));
	newedge -> readinterval = (READINTERVAL *) ckalloc(k * sizeof(READINTERVAL));
	vertex0 = edge1 -> end;
	n = 0;
	for(i = 0; i < edge2 -> multip; i ++)	{
	  /*
	    Find the interval that starts in edge1 that continues to interval[i] in edge2.
	   */
		c = FindIntervalEndingAtPos(edge1 -> readinterval, edge1 -> multip, edge2 -> readinterval[i].eq_read,
																edge2 -> readinterval[i].begin);
		/*
		  interval 'i' in edge2 started in edge1
		*/

		if(c >= 0)	{
			newedge -> readinterval[n].eq_read = edge1 -> readinterval[c].eq_read;
			newedge -> readinterval[n].offset  = edge1 -> readinterval[c].offset;
			newedge -> readinterval[n].begin   = edge1 -> readinterval[c].begin;
			newedge -> readinterval[n].length  = edge1 -> readinterval[c].length +
				edge2 -> readinterval[i].length - VERTEX_SIZE;
			/*
			  Mark this interval in edge1 ('c') as going to newedge, 
			  and this interval in edge2 ('i') as going to newedge.
			*/
			mark1[c] = 1;
			mark2[i] = 1;

			n ++;
		} else	{
		  /*
		    Interval 'i' just starts in edge2.
		  */
			c = -1;
			/* 
			   Check to see if the interval starts in a different edge.
			*/
			for(m = 0; m < vertex0 -> num_lastedge; m ++)	{
				if(vertex0 -> lastedge[m] == edge1)	continue;
				/* edge3 is an edge otehr than edge1 and edge2*/
				edge3 = vertex0 -> lastedge[m];
				/* Look to see if interval 'i' starts in edge 3*/
				c = FindIntervalEndingAtPos(edge3 -> readinterval, edge3 -> multip,
					 edge2 -> readinterval[i].eq_read,
					 edge2 -> readinterval[i].begin);
				if(c >= 0)	break;
			}
			/*  interval 'i' starts in edge 2.
			 */
			k = chklist(beginlist, nbegin, edge2 -> readinterval[i].eq_read);
			/*
			  If there are no other in edges than edge1 going into edge2, 
			  interval 'i' belongs in newedge.
			*/
			if(vertex0 -> num_lastedge == 1 || c < 0 && k < 0)	{
				newedge -> readinterval[n].eq_read = edge2 -> readinterval[i].eq_read;
				newedge -> readinterval[n].offset = edge2  -> readinterval[i].offset + l;
				newedge -> readinterval[n].begin  = edge2  -> readinterval[i].begin;
				newedge -> readinterval[n].length = edge2  -> readinterval[i].length;

				// Update the path.  It now starts
				// in newedge instead of edge2
				/*
				readIndex = edge2->readinterval[i].eq_read;
				if (paths != NULL ) {
					assert(paths[readIndex].edge[0] == edge2);
					paths[readIndex].edge[0] = newedge;
					paths[readIndex].begin_length = newedge->readinterval[n].offset;
				}
				*/
				if(edge2 -> readinterval[i].offset + edge2 -> readinterval[i].length >= edge2 -> length) {
					mark2[i] = 1;
				}
				n ++;
			}
		}
	}
	

	/*
	  Now check paths in edge1 to see if they go into edge2.
	*/
	for(i = 0; i < edge1 -> multip; i ++)	{
	  /* 
	     Already processed this edge.
	  */
		if(mark1[i] == 1)	continue;
		c = -1;
		/* 
		   Look to see if this read interval goes out into a different edge than edge2.
		 */
		for(m = 0; m < vertex0 -> num_nextedge; m ++)	{
			if(vertex0 -> nextedge[m] == edge2)	continue;
			edge3 = vertex0 -> nextedge[m];
			c = FindIntervalBeginningAtPos(edge3 -> readinterval, edge3 -> multip, edge1 -> readinterval[i].eq_read,
				 edge1 -> readinterval[i].begin + edge1 -> readinterval[i].length - VERTEX_SIZE);
			if(c >= 0)	break;
		}
		k = chklist(endlist, nend, edge1 -> readinterval[i].eq_read);
		// This path ends in edge1. 
		if(vertex0 -> num_nextedge == 1 || c < 0 && k < 0)	{
			newedge -> readinterval[n].eq_read = edge1 -> readinterval[i].eq_read;
			newedge -> readinterval[n].offset  = edge1 -> readinterval[i].offset;
			newedge -> readinterval[n].begin   = edge1 -> readinterval[i].begin;
			newedge -> readinterval[n].length  = edge1 -> readinterval[i].length;
			if(edge1 -> readinterval[i].offset == 0)	{
				mark1[i] = 1;
			}
			n ++;
		}
	}


	newedge -> multip = n;
	k = 0;
	/* Delete the re-routed read intervals.*/
	READINTERVAL *tmpIntvList;
	int newMult = 0;
	for(i = 0; i < edge1 -> multip; i ++)	{
	  if(!mark1[i])	{
	    edge1->readinterval[k++] = edge1->readinterval[i];
	  }
	}
	if (k > 0)
	  edge1->readinterval = (READINTERVAL*) realloc(edge1->readinterval, k * sizeof(READINTERVAL));
	else
	  edge1->readinterval = NULL;
	edge1->multip = k;
				      
	k = 0;
	for(i = 0; i < edge2 -> multip; i ++)	{
	  if(!mark2[i]) {
	    edge2->readinterval[k++] = edge2->readinterval[i];
	  }
	}
	if (k > 0)
	  edge2->readinterval = (READINTERVAL*) realloc(edge2->readinterval, k * sizeof(READINTERVAL));
	else
	  edge2->readinterval = NULL;
	edge2->multip = k;

	if (n > 0);
	//	  newedge->readinterval = (READINTERVAL*) realloc(newedge->readinterval,n * sizeof(READINTERVAL));
	else
	  newedge->readinterval = NULL;

	//	free((void *) readinterval);
	free((void *) mark1);
	free((void *) mark2);
}

int chklist(int *list, int n, int index)
{
	int	i, j, k, l;
	
	for(i = 0; i < n; i ++)	{
		if(list[i] == index)	return(i);
	}
	return(-1);
}

int FindIntervalEndingAtPos(READINTERVAL *readinterval, int num, int readindex, int readposition)
{
	int	i, j, k, l;
	/* Find the index of a read interval that ends before readposition*/
	for(i = 0; i < num; i ++)	{
		if(readinterval[i].eq_read == readindex &&
		   readinterval[i].begin + readinterval[i].length - VERTEX_SIZE == readposition)	{
			return(i);
		}
	}
	return(-1);
}

int FindIntervalBeginningAtPos(READINTERVAL *readinterval, int num, int readindex, int readposition)
{
  int	i, j, k, l;
  /*
    Find the read that begins at the same pos as readpos.
  */
	for(i = 0; i < num; i ++)	{
		if(readinterval[i].eq_read == readindex && readinterval[i].begin == readposition)	{
			return(i);
		}
	}
	return(-1);
}
