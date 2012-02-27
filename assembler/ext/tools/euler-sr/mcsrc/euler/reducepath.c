/***************************************************************************
 * Title:          reducepath.c
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


void derivelist(PATH *path, int num_path, NODES *vertex, EDGE *edge1, EDGE *edge2,
		int *edgematch1, int *edgematch2, int *beginlist, int *endlist, int *nbegin, int *nend);
void reducepath(PATH *path, int num_path, NODES *vertex, EDGE *edge1, EDGE *edge2, EDGE *newedge, 
		int *edgematch1, int *edgematch2);

void reducepath(PATH *path, int num_path, NODES *vertex, EDGE *edge1, EDGE *edge2, EDGE *newedge, 
		int *edgematch1, int *edgematch2)
{
  /* INPUT:
     path, num_path: all paths in the graph
     
     Transform the set edge1->vertex->edge2 into newedge.

  */
     

	int	i, j, k, l, num, m, c, n1, n2, n;
	int	*del_path;
	NODES	*begin, *end;

	del_path = (int *) ckalloc((MAX_BRA + vertex -> num_path) * sizeof(int));
	for (i = 0; i < MAX_BRA + vertex->num_path; i++) del_path[i] = 0;

	n1 = n2 = 0;
	num = vertex -> num_path;
	for(i = 0; i < num; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		/* If there are edges to delete on the path
			 and the path possibly ends in 'edge1' (path[k].len_path == j)
		*/
		if(path[k].len_path > 0 && j == path[k].len_path)	{
			/* This path indeed comes from edge1*/
			if(path[k].edge[j - 1] == edge1)	{
				/* This path ends in edge1, so it should be changed to end in newedge.*/
				if(edgematch1[n1])	{
					path[k].edge[j - 1] = newedge;
					if(vertex != newedge -> end)	{
						del_path[i] = 1;
					}
					/* append the path components to the dest of the new edge*/
					add_path(newedge -> end, k, j);
				}
				n1 ++;
			}
		} else if(path[k].len_path > 0 && j == 0)	{
			/* the path begins in this vertex, make sure
				 it leaves through edge2
			*/
			if(path[k].edge[j] == edge2)	{
				/* This path starts in edge2, change it to start in newedge*/
				if(edgematch2[n2])	{
					/* the path is consistent with all paths that pass through edge1,edge2*/
					int ol;
					ol = path[k].edge[0]->length - path[k].begin_length;
					path[k].edge[j] = newedge;
					if(vertex != newedge -> begin)	{
						del_path[i] = 1;
					}
					/* add this to the source of the new edge */
					add_path(newedge -> begin, k, j);
					path[k].begin_length = edge1->length + path[k].begin_length - VERTEX_SIZE;
				}
				n2 ++;
			}
		} else if(path[k].len_path > 0 && path[k].edge[j - 1] == edge1 && path[k].edge[j] == edge2)	{
			/* This path passes through edge1,edge2.
				 Update the paths accordingly.
			*/
			del_path[i] = 1;
			path[k].edge[j - 1] = newedge;
			remove_edge(&path[k], k, j);
		}
	}

	/* Shrink all the removed paths from this vertex.*/
	n = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		if(del_path[i] == 0)	{
			vertex -> path_index[n] = vertex -> path_index[i];
			vertex -> path_pos[n] = vertex -> path_pos[i];
			n ++;
		}
	}
	vertex -> num_path = n;

	free((void *) del_path);
}

void derivelist(PATH *path, int num_path, NODES *vertex, EDGE *edge1, EDGE *edge2,
		int *edgematch1, int *edgematch2, int *beginlist, int *endlist, int *nbegin, int *nend)
{
  /* 
     INPUT:
     path, num_path - all paths in the graph
     
     edge1->vertex->edge2
		 Find the edges that end in edge1 that are consistent with paths that leave into edge2.

  */
  /* 
     OUTPUT;
     
  */
	int	i, j, k, l, m, n1, n2, num;
	
	n1 = n2 = 0;
	num = vertex -> num_path;
	/* check each path that passes through this vertex*/
	for(i = 0; i < num; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		/* if the path does not begin/end at this vertex */
		if(path[k].len_path > 0 && j == path[k].len_path)	{
		  /* if the path 'k' passes through edge 1*/
			if(path[k].edge[j - 1] == edge1)	{
			  /* not sure what the readindex is talking about*/
			  
				if(path[k].readindex > 0 && !edgematch1[n1])	{
					endlist[(*nend) ++] = path[k].readindex - 1;
				} else if(path[k].pairindex[0] > 0 && !edgematch1[n1]) 	{
					endlist[(*nend) ++] = path[k].pairindex[1] - 1;
				}
				n1 ++;
			}
		} else if(path[k].len_path > 0 && j == 0)	{
			if(path[k].edge[j] == edge2)	{
				if(path[k].readindex > 0 && !edgematch2[n2])	{
					beginlist[(*nbegin) ++] = path[k].readindex - 1;
				} else if(path[k].pairindex[0] > 0 && !edgematch2[n2]) 	{
					beginlist[(*nbegin) ++] = path[k].pairindex[0] - 1;
				}
				n2 ++;
			}
		}
	}
}
