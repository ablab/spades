/***************************************************************************
 * Title:          clean_graph.h
 * Author:         Glenn Tesler
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/* clean graph
 * Glenn Tesler
 * 12/18/02
 */

typedef struct labelvert {
  int depth;                 /* labeling threshold */
  unsigned int dir       :1; /* 0=undirected version, 1=directed */
  unsigned int onepath   :1; /* For erosion:
			      * 1 = choose one max path arbitrarily,
			      * 0 = keep all max paths
			      */
} LABELVERT;

typedef struct {
  int w;               /* window size around indels */
  int max_cov;         /* coverage above this is clipped to this value */

  int C_B;             /* zigzag cycles of length <= C_B are bulges */
  int C_W;             /* directed cycles of length <= C_W are whirls */

  /* depths for classifying internal/external vertices
   * depth <= T from a leaf is considered external, else internal
   */
  LABELVERT *T_c;      /* threshold: chimeric read depth */
  LABELVERT *T_e;      /* threshold: erosion depth */
  //  int T_c;             /* threshold: chimeric read depth */
  //  int T_e;             /* threshold: erosion depth */

  /* low coverage cutoffs:
   * coverage <= L is low coverage
   */
  int L_b;             /* low coverage for bulge detection */
  int L_c;             /* low coverage for chimeric read detection */
  int L_e;             /* low coverage for erosion
			* (reserved) */

  /* phases */
  int do_zz;           /* true to clean up zigzag paths */
} CLEAN_PARAMS;

void clean_graph_once(IGRAPH *G, READTABLE *RT, CLEAN_PARAMS *params);

void graph_grow(IGRAPH *G,
		READTABLE *RT,
		CLEAN_PARAMS *params,
		int treerule);

void clean_graph_adapter(IGRAPH *G, READTABLE *RT);
/*
void clean_graph_adapter(NODES **nodes, int num_nodes,
			 char **src_seq, int num_seq, int *len_seq,
			 LIST **readlist);
*/
