/***************************************************************************
 * Title:          mapread.c
 * Author:         Haixu Tang
 * Created:        Jun. 2003
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

#define BACK_GAP 10

int searchnextpath(NODES *node1, NODES *node2, READLIST *readlist, int max_len, int len, char *spath);
void destroyreadlist(READLIST *newreadlist, int index, int len);
int makebothreadlist(READLIST *newreadlist, int index, int len, int num_seq, int len_or);
int fillreadlist(READLIST *readlist, NODES *node_last, NODES *node, int pos_last, int pos, IGRAPH *G);
void mapread(IGRAPH *G, READTABLE *RT);

void mapread(IGRAPH *G, READTABLE *RT)
{
  int	i, j, k, l, m, n, t, p1, p2;
  int	last, last_n, gapl;
  char	*clist;
  READPOSITION	*position;
  EDGE	*edge;
  NODES	*node, *node_next, *node_last, *node1, *node2;
  READLIST	*newreadlist, *newreadlist_copy;
  IVEDGE	*itver;
  IEDGE	*itever;

  /*	Remode all read intevals in each edge	*/
  itever = it_e_new(G, E_G);
  while((edge = it_e_next(itever)) != (EDGE *) NULL)	{
    edge -> multip = 0;
    free((void *) edge -> readinterval);
  }
  it_e_destroy(itever);

  /*	Compute the largest length of all the edges	*/
  l = 0;
  for(i = 0; i < RT -> num_seq; i ++)	{
    if(l < RT -> len_seq[i])	{
      l = RT -> len_seq[i];
    }
  }
  newreadlist_copy = (READLIST *) ckalloc(l * 2 * sizeof(READLIST));
  newreadlist = (READLIST *) ckalloc(l * 2 * sizeof(READLIST));
  clist = (char *) ckalloc(l * 2 * sizeof(char));

  node = (NODES *) NULL;
  // for each sequence
  for(i = 0; i < RT -> num_seq; i ++)	{
    /*	Compute the list of nodes in graph H	*/
    // ???? char list?
    for(j = 0; j < RT -> len_seq[i]; j ++)	
      clist[j] = 0;
	  
    // for each position
    for(j = 0; j < RT -> len_seq[i] - 1; j ++)	{
      // not sure what k is
      k = 0;
      node = RT -> readlist[i][j].node;
      // subg = subgraph, what is a subgraph?
      if(node -> subg_flag)	{
	// the next node threaded by this read
	node_next = RT -> readlist[i][j + 1].node;
	
	itver = it_ev_new(node, E_OUT_H);
	// I think this is looking for whirls, but I'm not sure.
	while((edge = it_ev_next(itver)) != (EDGE *) NULL)	{
	  if(edge -> end == node_next)	{
	    k = 1;
	  }
	}
	it_ev_destroy(itver);
      }
      if(k == 0)	{
	clist[j] = clist[j + 1] = 1;
      }
    }
    n = 0;
    last = -1;
    gapl = j = 0;
    
    while(j < RT -> len_seq[i])	{
      node = RT -> readlist[i][j].node;
      if(gapl == 0 && clist[j] == 0)	{
	newreadlist[n ++].node = node;
	last = j;
      } else	{
	if(last >= 0)	{
	  for(; j < RT -> len_seq[i]; j ++)	{
	    if(clist[j] == 0)	{
	      break;
	    }
	  }
	  if(j < RT -> len_seq[i])	{
	    node = RT -> readlist[i][last].node;
	    node_next = RT -> readlist[i][j].node;
	    k = fillreadlist(&newreadlist[n],
			     node, node_next, last, j, G);
	    if(k > 0)	{
	      n += k;
	      gapl = 0;
	      last = j;
	    } else	{
	      last_n = last;
	      for(m = 0; m < n; m ++)	{
		newreadlist_copy[m] = newreadlist[m];
	      }

	      /*	Jump back 1 to go forward again until finding the path	*/
	      last -= BACK_GAP;
	      m = n - 1;
	      while(last >= 0)	{
		node = RT -> readlist[i][last].node;
		for(; m >= 0; m --)	{
		  if(newreadlist[m].node == node)	break;
		}
		if(m < 0)	{
		  m = n - 1;
		  last --;
		  continue;
		}
		k = fillreadlist(&newreadlist[m + 1],
				 node, node_next, last, j, G);
		if(k > 0)	{
		  n = m + k + 1;
		  last = j;
		  gapl = 0;
		  break;
		} else	{
		  last --;
		}
	      }
	      if(last < 0)	{
		last = last_n;
		for(m = 0; m < n; m ++)	{
		  newreadlist[m] = newreadlist_copy[m];
		}
		gapl = 1;
	      }
	    }
	  }
	}
      }
      j ++;
    }
    /*	update the readlist of read positions in each supernode and read intevals in each edge	*/
    destroyreadlist(RT -> readlist[i], i, RT -> len_seq[i]);
    destroyreadlist(RT -> readlist[i + RT -> num_seq], i + RT -> num_seq,
		    RT -> len_seq[i + RT -> num_seq]);
    n = makebothreadlist(newreadlist, i, n, RT -> num_seq, RT -> len_seq[i]);
    RT -> new_len_seq[i] = RT -> new_len_seq[i + RT -> num_seq] = n;
  }

  free((void *) clist);
  free((void *) newreadlist);
  free((void *) newreadlist_copy);
}

void destroyreadlist(READLIST *newreadlist, int index, int len)
{
  int	i, j, k, l;
  NODES	*node;

  for(i = 0; i < len; i ++)	{
    node = newreadlist[i].node;
    rem_position(node, index, i);
  }
	
}

int makebothreadlist(READLIST *newreadlist, int index, int len, int num_seq, int len_or)
{
  int	i, j, k, l;
  NODES	*node, *node_next;
  EDGE	*edge1, *edge2;

  for(i = 0; i < len - 1; i ++)	{
    node = newreadlist[i].node;
    node_next = newreadlist[i + 1].node;
    for(j = 0; j < node -> num_nextedge; j ++)      {
      if(node -> nextedge[j] -> subg_flag &&
	 node -> nextedge[j] -> end == node_next)     break;
    }
    if(j == node -> num_nextedge)    {
      printf("Warning: edge not found %d(%d,%d)\n", index, i, i + 1);
      return(0);
    }
  }

  len = min(len, len_or);
  for(i = 0; i < len; i ++)	{
    node = newreadlist[i].node;
    insert_position(node, index, i);
    insert_position(node -> bal_node, index + num_seq, len_or - 1 - i);
    if(i < len - 1)	{
      node_next = newreadlist[i + 1].node;
      for(j = 0; j < node -> num_nextedge; j ++)      {
	if(node -> nextedge[j] -> subg_flag &&
	   node -> nextedge[j] -> end == node_next)     break;
      }
      if(j == node -> num_nextedge)	{
	printf("index %d i %d %d\n", index, i, i + 1);
	exit(-1);
      }
      insert_interval(node -> nextedge[j], index, i, i + 1);
      insert_interval(node -> nextedge[j] -> bal_edge, index + num_seq,	
		      len_or - 2 - i, len_or - 1 - i);
    }
  }
  return(len);
}

int fillreadlist(READLIST *readlist, NODES *node_last, NODES *node, int pos_last, int pos, IGRAPH *G)
{
  int	i, j, k, l;
  char	spath;
  int	len;

  l = pos - pos_last + 1;
  len = spath = 0;
  len = searchnextpath(node_last, node, readlist, l + min(MAX_DIF, WhirlLength), len, &spath);
  return(len);
}

int searchnextpath(NODES *node1, NODES *node2, READLIST *readlist, int max_len, int len, char *spath)
{
  int	i, j, k, l, p, m, plen;
  EDGE	*edge;

  if(max_len < 0)	return(len);

  plen = len;
  for(i = 0; i < node1 -> num_nextedge; i ++)	{
    if(*spath > 1)	return(0);;
    edge = node1 -> nextedge[i];
    if(!edge -> subg_flag)		continue;
    if(edge -> end == node2)	{
      if(*spath >= 1)	{	//found path already?
	*spath = 2;
	return(0);
      } else	{
	readlist[len ++].node = node2;
	*spath = 1;
      }
    } else if(edge -> visit == 0)	{
      k = len;
      p = *spath;
      readlist[len ++].node = edge -> end;
      edge -> visit = 1;
      len = searchnextpath(edge -> end, node2, readlist, max_len - edge -> length + 1, len, spath);
      if(*spath == 0)	{
	len = plen;
      } else if(p == 1 && *spath == 1)	{
	len = k;
      }
      edge -> visit = 0;
    }
  }
  if(*spath > 1)		{
    return(0);
  }
  else			return(len);
}
