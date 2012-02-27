/***************************************************************************
 * Title:          floyd_graph.h
 * Author:         Glenn Tesler
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <values.h>
#define SHP_INF (MAXINT >> 1)

typedef struct shortest_paths {
  int n;                          /* number of vertices */

  /* order the nodes
   * v_order[0..(n-1)] lists the nodes in ascending order
   * node i is v_order[i]
   */
  NODES **v_order;

  /* shortest distance matrix
   * d[i+j*n] = shortest distance from vertex #i to vertex #j
   */
  int *d;

  /* shortest path matrix
   * p[i+j*n] = intermediate vertex on route from i to j
   */
  int *p;
} SHORTEST_PATHS;


int SHP_lookup_dist(SHORTEST_PATHS *sh_paths,
		    NODES *v,
		    NODES *w);

int SHP_lookup_dist_vi(SHORTEST_PATHS *sh_paths,
		       NODES *v,
		       int i_w);

void SHP_destroy(SHORTEST_PATHS *sh_paths);

void SHP_compute_all(IGRAPH *G,
		     SHORTEST_PATHS *sh_paths);

int SHP_cycle(SHORTEST_PATHS *sh_paths,
	      NODES *v);

int SHP_node_num(SHORTEST_PATHS *sh_paths, NODES *v);
