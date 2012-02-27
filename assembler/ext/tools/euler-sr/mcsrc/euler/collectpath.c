/***************************************************************************
 * Title:          collectpath.c
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

int collect2forpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midforpath,
		     PATH *path, int num);
int collect2aftpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midaftpath,
		     PATH *path, int num);
int collectstartpaths(NODES *vertex, EDGE *edge, PATH *startpath, PATH *path, int num);
int collectendpaths(NODES *vertex, EDGE *edge, PATH *endpath, PATH *path, int num);

int collect2forpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midforpath,
		     PATH *path, int num)
{
  /* 
     vertex, edge1, and edge2: part of the graph of the format --edge1-->vertex--edge2-->
     midforpath: storage for paths that are before edge1
     path: all paths
     num: NOT USED

     Result: Stores in midforpath all unique, maximal paths that contain (edge1,edge2). So if there
     are two paths:  (c,b,a,edge1,edge2), and (d,c,b,a,edge1,edge2), then midforpath will
     contain [(d,c,b,a,edge1,edge2)] only.

     Return: the number of paths stored in midforpath.
  */
     
        int	i, j, k, l, m, n, c, q;
	int	num_path;
	PATH	tmppath;
	
	/*	tmppath.edge = (EDGE **) ckalloc(MAX_TMP_LEG * sizeof(EDGE *));
	 */
	
	num_path = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		/* If this vertex is not the beginning vertex not end vertex of a path
		   AND the source edge for this vertex is 'edge1' and dest edge on this 
		   path is 'edge2' 
		*/
		if(j != 0 && j != path[k].len_path && path[k].edge[j - 1] == edge1 && path[k].edge[j] == edge2)	{
			tmppath.len_path = j - 1;
			if(tmppath.len_path == 0)	continue;
		  tmppath.edge = (EDGE**) ckalloc(tmppath.len_path * sizeof(EDGE*));
			/* 
			   Store the edges in the path before 'edge1' in 'tmppath'.
			*/
			for(m = 0; m < j - 1; m ++)	{
				tmppath.edge[m] = path[k].edge[j - 2 - m];
			}
			/* 
			   Find the path in 'midforpath' is contained in tmppath.

			*/
			c = chk_consist(&tmppath, midforpath, num_path, &q);
			if(c > 0)	{
			  /*
			    'tmppath' contains a path in midforpath.  Extend the length of that path
			    to include the entire new tmppath.
			  */
				midforpath[c - 1].len_path = tmppath.len_path;
				free((void **) midforpath[c - 1].edge);
				midforpath[c - 1].edge = (EDGE **) ckalloc((midforpath[c - 1].len_path + 1) * sizeof(EDGE *));
				for(m = 0; m < tmppath.len_path; m ++)	{
					midforpath[c - 1].edge[m] = tmppath.edge[m];
				}
			} else if(c < 0)	{
			  /*
			    'tmppath' is not contained in midforpath.  Append it to the list of paths.
			  */
				midforpath[num_path].len_path = tmppath.len_path;
				midforpath[num_path].edge = (EDGE **) ckalloc((midforpath[num_path].len_path + 1) * sizeof(EDGE *));
				for(m = 0; m < tmppath.len_path; m ++)	{
					midforpath[num_path].edge[m] = tmppath.edge[m];
				}
				num_path ++;
			}
			free(tmppath.edge);
		}
	}

/*
  free((void **) tmppath.edge);
*/
	return(num_path);
}

int collect2aftpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midaftpath,
		     PATH *path, int num)
{
	int	i, j, k, l, m, n, c, q;
	int	num_path;
	PATH	tmppath;

	/*
	  tmppath.edge = (EDGE **) ckalloc(MAX_TMP_LEG * sizeof(EDGE *));
	*/

	num_path = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(j != 0 && j != path[k].len_path && path[k].edge[j - 1] == edge1 && path[k].edge[j] == edge2)	{
			tmppath.len_path = path[k].len_path - j - 1;
			if(tmppath.len_path == 0)	continue;
			tmppath.edge = (EDGE**) ckalloc(tmppath.len_path * sizeof(EDGE*));
			for(m = j + 1; m < path[k].len_path; m ++)	{
				tmppath.edge[m - j - 1] = path[k].edge[m];
			}
			c = chk_consist(&tmppath, midaftpath, num_path, &q);
			if(c > 0)	{
				midaftpath[c - 1].len_path = tmppath.len_path;
				free((void **) midaftpath[c - 1].edge);
				midaftpath[c - 1].edge = (EDGE **) ckalloc((midaftpath[c - 1].len_path + 1) * sizeof(EDGE *));
				for(m = 0; m < tmppath.len_path; m ++)	{
					midaftpath[c - 1].edge[m] = tmppath.edge[m];
				}
			} else if(c < 0)	{
				midaftpath[num_path].len_path = tmppath.len_path;
				midaftpath[num_path].edge = (EDGE **) ckalloc((midaftpath[num_path].len_path + 1) * sizeof(EDGE *));
				for(m = 0; m < tmppath.len_path; m ++)	{
					midaftpath[num_path].edge[m] = tmppath.edge[m];
				}
				num_path ++;
			}
			free((void **) tmppath.edge);
		}
	}


	return(num_path);
}

int collectendpaths(NODES *vertex, EDGE *edge, PATH *endpath, PATH *path, int num)
{
	int	i, j, k, l, m, n, c;
	int	num_path;
	/*
		Store the paths that end in 'edge'.  Edge is an in-edge to 'vertex'.
	  In terms of the paper, P-->vertex
	*/
	num_path = 0;
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path > 0 && j == path[k].len_path && path[k].edge[j - 1] == edge)	{
			endpath[num_path].len_path = j - 1;
			endpath[num_path].edge = (EDGE **) ckalloc((endpath[num_path].len_path + 1) * sizeof(EDGE *));
			for(m = 0; m < j - 1; m ++)	{
				endpath[num_path].edge[m] = path[k].edge[j - 2 - m];
			}
			endpath[num_path].readindex = path[k].readindex;
			num_path ++;
		}
	}

	return(num_path);
}

int collectstartpaths(NODES *vertex, EDGE *edge, PATH *startpath, PATH *path, int num)
{
	int	i, j, k, l, m, n, c;
	int	num_path;

	num_path = 0;
	/* 
	   Store paths that start in 'edge' in 'startpath'. 
	   In terms of the paper, this is Pvertex-->
	*/
	for(i = 0; i < vertex -> num_path; i ++)	{
		k = vertex -> path_index[i];
		j = vertex -> path_pos[i];
		if(path[k].len_path > 0 && j == 0 && path[k].edge[j] == edge)	{
			startpath[num_path].len_path = path[k].len_path - 1;
			startpath[num_path].edge = (EDGE **) ckalloc((startpath[num_path].len_path + 1) * sizeof(EDGE *));
			for(m = 0; m < startpath[num_path].len_path; m ++)	{
				startpath[num_path].edge[m] = path[k].edge[m + 1];
			}
			startpath[num_path].readindex = path[k].readindex;
			num_path ++;
		}
	}

	return(num_path);
}
