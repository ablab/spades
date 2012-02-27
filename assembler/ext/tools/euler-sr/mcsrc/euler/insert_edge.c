/***************************************************************************
 * Title:          insert_edge.c
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

EDGE *insert_edge(NODES *node, NODES *node_next, int read, int pos1, int pos2);
void add_nextedge(NODES *node, EDGE *edge);
void add_lastedge(NODES *node, EDGE *edge);
void insert_interval(EDGE *edge, int read, int pos1, int pos2);

EDGE *insert_edge(NODES *node, NODES *node_next, int read, int pos1, int pos2)
{
  /* 
		 Insert an edge into a read path.
     Link node to node_next, and add an interval representing read, pos1, and pos2.
     If the edge already exists, only add the interval.
  */
        int     i, j;
        EDGE *edge;

        for(i = 0; i < node -> num_nextedge; i ++)      {
                if(node -> nextedge[i] -> end == node_next)     break;
        }
        if(i < node -> num_nextedge)    {
					insert_interval(node -> nextedge[i], read, pos1, pos2);
					return((EDGE *) NULL);
        }

        edge = (EDGE *) ckalloc(1 * sizeof(EDGE));
        edge -> begin = node;
        edge -> end = node_next;
				edge->deleted = 0;
				insert_interval(edge, read, pos1, pos2);
	/*
	  Incorporate this edge into the graph.
	*/
        add_nextedge(node, edge);
        add_lastedge(node_next, edge);
        return(edge);
}

void insert_interval(EDGE *edge, int read, int pos1, int pos2)
{
	int	i, j, k, l;
	READINTERVAL	*readinterval0, readinterval;

	readinterval.eq_read = read;
	readinterval.offset = 0;
	readinterval.begin = pos1;
	readinterval.length = pos2 - pos1 + 1;
	readinterval.next = (READINTERVAL *) NULL;

	edge->readinterval = (READINTERVAL*) realloc(edge->readinterval, (1 + edge->multip)* sizeof(READINTERVAL));
        edge -> multip ++;
	edge->readinterval[edge->multip-1] = readinterval;
	  /*
	    readinterval0 = (READINTERVAL *) ckalloc((1 + edge -> multip) * sizeof(READINTERVAL));
	    for(j = 0; j < edge -> multip; j ++)	{
	    readinterval0[j] = edge -> readinterval[j];
	    }
	    if(edge -> multip > 0)		 free((void *) edge -> readinterval);
	    edge -> readinterval = (READINTERVAL *) ckalloc(edge -> multip * sizeof(READINTERVAL));
	    for(j = 0; j < edge -> multip - 1; j ++)	{
	    edge -> readinterval[j] = readinterval0[j];
	    }
	    edge -> readinterval[j] = readinterval;
	    free((void *) readinterval0);
	  */
	if(readinterval.length > edge -> length)	{
		edge -> length = readinterval.length;
	}
}

void add_nextedge(NODES *node, EDGE *edge)
{
  /* 
     node:  a node in the graph to add an edge to. 
     edge:  an edge. The edge should be fully initialized.

     Result:
     If 'edge' is already in the list of edges for 'node', nothing.
     Otherwise, append it to the list of nextedges for node.  I'm not
     sure why a temporary array is created, perhaps just a realloc 
     will be fine.
  */

        int     i, j, k, l;
        EDGE    **tmpedge;

        tmpedge = (EDGE **) ckalloc((node -> num_nextedge + 1) * sizeof(EDGE *));
        for(i = 0; i < node -> num_nextedge; i ++)      {
                tmpedge[i] = node -> nextedge[i];
                if(edge == tmpedge[i])  {
                        break;
                }
        }
	/*
	  If 'edge' has not been found in node->nextedge, 
	  adde 'edge' to the list of nextedge edges.
	*/
        if(i == node -> num_nextedge)   {
                if(node -> num_nextedge > 0)    free((void **) node -> nextedge);
                node -> num_nextedge ++;
                node -> nextedge = (EDGE **) ckalloc(node -> num_nextedge * sizeof(EDGE *));
                for(i = 0; i < node -> num_nextedge - 1; i ++)  {
                        node -> nextedge[i] = tmpedge[i];
                }
                node -> nextedge[i] = edge;
        }
        free((void **) tmpedge);
}

void add_lastedge(NODES *node, EDGE *edge)
{
  /* 
     node:  a node in the graph to add a last (in) edge to.
     edge:  an edge. The edge shoudl be fully initialized.

     Result:
     If 'edge' is already in the list of lastedge for 'node', nothing.
     Otherwise, append it to the list of nextedges for node.  
     Same comment as before... not sure why a temporary array is created.
  */
        int     i, j, k, l;
        EDGE    **tmpedge;

        tmpedge = (EDGE **) ckalloc((node -> num_lastedge + 1) * sizeof(EDGE *));
        for(i = 0; i < node -> num_lastedge; i ++)      {
                tmpedge[i] = node -> lastedge[i];
                if(edge == tmpedge[i])  {
                        break;
                }
        }
        if(i == node -> num_lastedge)   {
                if(node -> num_lastedge > 0)    free((void **) node -> lastedge);
                node -> num_lastedge ++;
                node -> lastedge = (EDGE **) ckalloc(node -> num_lastedge * sizeof(EDGE *));
                for(i = 0; i < node -> num_lastedge - 1; i ++)  {
                        node -> lastedge[i] = tmpedge[i];
                }
                node -> lastedge[i] = edge;
        }
        free((void **) tmpedge);
}
