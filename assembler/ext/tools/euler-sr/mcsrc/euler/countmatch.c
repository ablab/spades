/***************************************************************************
 * Title:          countmatch.c
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

int searchlast(NODES *vertex, EDGE *edge, PATH *path, int num_path, int *match);
int searchnext(NODES *vertex, EDGE *edge, PATH *path, int num_path, int *match);
int countmatch(EDGE *edge1, EDGE *edge2, PATH *path, int num_path);
int countstartmatch(EDGE *edge, NODES *vertex, PATH *path, int num_path);
int countendmatch(EDGE *edge, NODES *vertex, PATH *path, int num_path);

int countstartmatch(EDGE *edge, NODES *vertex, PATH *path, int num_path)
{
  /*
		Count the number of edges leaving 'vertex' for which there is a path that
		enters through edge and passes through 'vertex' (vertex is not at the end of the path).

  */
	int	i, j, k, l, n, m;
	int	*match;
	match = (int *) ckalloc(vertex -> num_nextedge * sizeof(int));
	/*
	  For each path that passes through this vertex
	*/
	for (i = 0; i < vertex->num_nextedge; i++ ) {
	  match[i] = 0;
	}
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		/* 
		   Fisrt two conditions:
		      j != 0 - the path does not start in 'edge'.
  		      j != path[i].len_path - the path does not end in 'edge'
		      path[k].edge[j-1] == edge - the is path follows the route we're looking for
		         (through edge).
		*/
		   
		if(j != 0 && j != path[k].len_path && path[k].edge[j - 1] == edge)	{
			/* 
				 There is a path that enters 'vertex' through 'edge'.  
			   The next edge on the path is path[k].edge[j].  Count the number of 
			   paths that leave edge 'n' through 'vertex'
			*/
			n = searcherase(vertex -> nextedge, path[k].edge[j], vertex -> num_nextedge);
			match[n] ++;
		}
	}
	/* Count the number of paths that leave 'vertex' that came from 'edge' */

	n = 0;
	for(i = 0; i < vertex -> num_nextedge; i ++)	{
		if(match[i] > 0)	{
			n ++;
		}
	}
	free((void *) match);
	return(n);
}

int countmatch(EDGE *edge1, EDGE *edge2, PATH *path, int num_path)
{
  /* 
     Return the number of times that paths tha contain 'edge1->vertex->edge2'.
  */

	int	i, j, k, l, n, m;
	NODES	*vertex;

	if(edge1 -> end != edge2 -> begin)	{
		printf("edge1 %d %d edge2 %d %d\n", edge1, edge1 -> end, edge2, edge2 -> begin);
		printf("Edge not connected in countmatch\n");
		exit(-1);
	}

	vertex = edge1 -> end;
	n = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j != 0 && j != path[k].len_path && path[k].edge[j - 1] == edge1 &&
		   path[k].edge[j] == edge2)	{
			n ++;
		}
	}

	return(n);
}

int searchnext(NODES *vertex, EDGE *edge, PATH *path, int num_path, int *match)
{
  /* input: 
     vertex: a vertex in the graph
     edge : 
     path: all paths in the graph
     match: 
  */
	int	i, j, k, l, n, m, *next;
	
	next = (int *) ckalloc(vertex -> num_nextedge * sizeof(int));
	for (i = 0; i < vertex->num_nextedge;i++) 
	  next[i] = 0;

	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		/* if previous edge on this path is 'edge', and 
		   the next edge on the path isn't 'edge' (the path
		   is not a simple cycle)
		*/
		if(j > 0 && j < path[k].len_path && path[k].edge[j - 1] == edge &&
		   path[k].edge[j] != edge)	{
		   /* 
		     m is the index of path[k].edge[j] in the list of edges vertex->nextedge
		   */
		  m = searcherase(vertex -> nextedge, path[k].edge[j], vertex -> num_nextedge);
		  next[m] = 1;
		}
	}

	n = 0;
	for(i = 0; i < vertex -> num_nextedge; i ++)	{
	  /* If there is a path that leaves this vertex through next [i]
	     and that edge doesn't end at this vertex (a cycle), record the 
	     index of the next edge in 'match'
	  */
		if(next[i] == 1 && vertex -> nextedge[i] -> end != edge -> begin)	{
			match[n ++] = i;
		}
	}

	free((void *) next);

	return(n);
}

int searchlast(NODES *vertex, EDGE *edge, PATH *path, int num_path, int *match)
{
	int	i, j, k, l, n, m, *last;
	/* 
	   Count the number of times paths leave 'vertex' through 'edge' for each 
	   edge into 'vertex'.
	*/
	last = (int *) ckalloc(vertex -> num_lastedge * sizeof(int));
	for (i = 0; i < vertex->num_lastedge; i++ ){
	  last[i] = 0;
	}
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j > 0 && j < path[k].len_path && path[k].edge[j] == edge &&
		   path[k].edge[j - 1] != edge)	{
			m = searcherase(vertex -> lastedge, path[k].edge[j], vertex -> num_lastedge);
			last[m] = 1;
		}
	}

	n = 0;
	for(i = 0; i < vertex -> num_lastedge; i ++)	{
		if(last[i] == 1 && vertex -> lastedge[i] -> begin != edge -> end)	{
			match[n ++] = i;
		}
	}

	free((void *) last);

	return(n);
}

int countendmatch(EDGE *edge, NODES *vertex, PATH *path, int num_path)
{
	int	i, j, k, l, n, m;
	int	*match;

	/* Count the number of edges entering vertex for which there exists
		 paths that go through vertex.
	 */
	match = (int *) ckalloc(vertex -> num_lastedge * sizeof(int));
	for (i = 0; i < vertex->num_lastedge; i++ )
	  match[i] = 0;

	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j != 0 && j != path[k].len_path && path[k].edge[j] == edge)	{
			n = searcherase(vertex -> lastedge, path[k].edge[j - 1], vertex -> num_lastedge);
			match[n] ++;
		}
	}

	n = 0;
	for(i = 0; i < vertex -> num_lastedge; i ++)	{
		if(match[i] > 0)	{
			n ++;
		}
	}
	free((void *) match);
	return(n);
}
