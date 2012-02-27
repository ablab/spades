/***************************************************************************
 * Title:          chk_consist.c
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

int chk_consist(PATH *startpath, PATH *midpath, int num_midpath, int *n);

int chk_consist(PATH *startpath, PATH *midpath, int num_midpath, int *n)
{
	int	i, j, k, l, c;
	/* 
	   startpath: a path from an edge stored in reverse order
	   midpath: an array of paths in reverse order from the same edge
	   n: will store the last index of a path in midpath that contains
	   path as a subpath.
	   
	   Returnval: 
	     Returns -1 if 'startpath' is not contained in any midpath.
	     Returns 0 if 'startpath' is contained in a midpath.
	     Return > 0 if 'startpath' contains a path in midpath.
	*/
	c = -1;
	for(i = 0; i < num_midpath; i ++)	{
	  /* If startpath may be contained in midpath[i], look to see
	     if it is indeed contained in it.
	  */
		if(startpath -> len_path <= midpath[i].len_path)	{
			for(j = 0; j < startpath -> len_path; j ++)	{
				if(startpath -> edge[j] != midpath[i].edge[j])	{
					break;
				}
			}
			/* 
			   Startpath is contained in midpath.  
			*/
			if(j == startpath -> len_path)	{
				c = 0;
				*n = i;
				break;
			}
		} else	{
		  /* 
		     'path' is longer than midpath[i].  If 'path' contains 'midpath[i]'
		     as a subpath, store the length of the longest path in c.
		  */
			for(j = 0; j < midpath[i].len_path; j ++)	{
				if(startpath -> edge[j] != midpath[i].edge[j])	{
					break;
				}
			}
			if(j == midpath[i].len_path)	{
				c = i + 1;
				break;
			}
		}
	}

	return(c);
}
