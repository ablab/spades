/***************************************************************************
 * Title:          erosion.c
 * Author:         Glenn Tesler
 * Created:        Jun. 2003
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

/* Erosion
 * Glenn Tesler
 * 6/13/03
 */


#include <stdinc.h>
#include <clean_graph.h>
#include <erosion.h>

//#define DEBUG_OUT
#undef DEBUG_OUT


/* private */
void erode_external_undir(IGRAPH *G, CLEAN_PARAMS *params);
int erode_tree(IGRAPH *G, CLEAN_PARAMS *params, EDGE *e);
void erode_small_components(IGRAPH *G, CLEAN_PARAMS *params);
void mark_path(IGRAPH *G, NODES *v);
int erode_collapse_tree(NODES *v, EDGE *e);
void debug_print_edges(NODES *v, E_TYPE e_type);


void label_external_undir(IGRAPH *G, LABELVERT *T_e);
void label_external_dir(IGRAPH *G, LABELVERT *T_e);
void label_external_dir_v(LABELVERT *T_e, NODES *v);
void erode_external_dir(IGRAPH *G, LABELVERT *T_e);
void erode_edge_dir(IGRAPH *G, LABELVERT *T_e, NODES *v);

/****************************************************************************
 * Erosion and chimeric read detection using depth labelling
 ****************************************************************************
 *
 * undirected labeling (original algorithm):
 * BFS version on nucleotide graph:
 *
 * 1. Label all verts infinity
 * 2. k=1
 * 3. label sources and sinks k
 * 4. delete sources and sinks
 * 5. k = k+1
 * 6. Goto 3
 * 
 ****************************************************************************
 * 
 * DFS on interval graph:
 * 
 * 1. Label all verts infinity
 * 2. label sources and sinks 1
 * 3. for v in source/sinks do {
 *      w = v;
 *      k = label(v)
 *      for e in edges(v) {
 *         x = other vertex of e
 *         new_label(x,v)
 *      }
 *    }
 *
 * new_label(v,z) {
 *    m1 = max(label(x) + length(e))
 *    where e=(x,v) runs over all incoming edges to v
 * 
 *    m2 = max(label(x) + length(e))
 *    where e=(v,x) runs over all outgoing edges from v
 * 
 *    m = min(m1,m2)
 * 
 *    if (m < label(v)) {
 *       set label(v) = m
 *       for x in neighbors(v) - {z} {
 *          new_label(x)
 *       }
 *    }
 * }
 ****************************************************************************/
/****************************************************************************
 *
 * Erosion on nucleotide graph:
 *
 * undirected version on nucleotide graph: 
 * 1. label vertices with undirected labelling algorithm
 * 2. vertices with labels <= T_e are "external", others are "internal"
 * 3. For every edge with vertex labels {T_e+1,T_e} {
 *       Follow a steepest descent path (vertex labels {x,x-1})
 *       down to label 1, marking the vertices along it as "internal"
 *    }
 * 4. Delete all edges except internal-internal
 *
 *
 *
 * Directed version on interval graph:
 *
 * 1. label vertices with directed labelling algorithm
 *
 * 2. vertices with labels <= T_e are "external", others are "internal"
 *
 * 3. For every edge e {
 *       labels are k1, k2, length is L1
 *       if (k1+k2+L1 >= 2T_e+1) {
 *          mark e as internal
 *       } else {
 *          mark e as external
 *       }
 *    }
 *
 * 4. For every edge e with k1+k2+L1 >= 2T_e+1 {
 *       Follow a steepest descent path, stopping if
 *       hit a vertex with label 1 or an internal vertex.
 *       Mark the edges and vertices on the path as internal.
 *    }
 *
 * 5. Delete all external edges
 ****************************************************************************/

/****************************************************************************
 * Labels on interval graph instead of nucleotide graph:
 *
 *    edge length L0 = current method of storing lengths
 *                L1 = L0-1
 * 
 *    So edge [A] --- AGGCT ---> [T]
 *    has L0=5, L1=4
 * 
 *    Suppose labels are
 *    k1 -----> k2
 *    with some length L1
 * 
 *    In nucleotide graph, this would be
 *       k1 -> k1+1 -> ... -> r-1 -> r -> r-1 -> ... -> k2
 *       with k1 + k2 + L1 = 2r
 *    or
 *       k1 -> k1+1 -> ... -> r -> r -> ... -> k2
 *       with k1 + k2 + L1 = 2r+1
 *    depending on parity.
 * 
 *    So r=floor((k1+k2+L1)/2)
 * 
 *    Cases:
 *       |k2 - k1| > L1:
 * 
 *       |k2 - k1| <= L1:
 *          k1 + k2 + L1 >= 2T_e+3:
 *              nucleotide expansion of edge would have int-int
 * 
 *          k1 + k2 + L1 == 2T_e+2:
 *              nucleotide expansion of edge would just reach
 *                  T_e -> T_e+1 -> T_e,
 *              or  T_e -> T_e+1 and end
 *              hence has ext-int
 * 
 *          k1 + k2 + L1 <= 2T_e + 1:
 *              nucleotide expansion of edge is all ext-ext
 *
 *  so      k1 + k2 + L1 >  2T_e + 2: interior edge
 *                       ==           exterior (but "borderline interior")
 *                       <            exterior
 ****************************************************************************/






/****************************************************************************
 * erosion (undirected version):
 *
 * trees hanging off the core should be trimmed or condensed
 ****************************************************************************/

void erode_graph(IGRAPH *G, CLEAN_PARAMS *params)
{
  if (params->T_e->depth <= 0) {
    printf("Erosion skipped.\n");
    return;
  }

  /* classify nodes and edges as internal/external */
  classify_internal_nodes(G, params->T_e);

  /* choose steepest descent paths from core down to ends
   * and delete all the other edges and vertices
   */
  erode_external_dir(G, params->T_e);

#if 0
  /* erosion of trees hanging off an internal core */
  if (params->T_e->dir) {
    erode_external_dir(G, params->T_e);
  } else {
    erode_external_undir(G, params->T_e);

    /* erosion of components of H that are entirely external vertices */
    erode_small_components(G, params);
  }
#endif
}


/* erode trees emanating from internal-external edges */
void erode_external_undir(IGRAPH *G, CLEAN_PARAMS *params)
{
  IEDGE *it_edge;
  EDGE *e;

  int num_eroded_trees = 0;
  int num_eroded_edges = 0;
  int num_collapsed_trees = 0;
  int num_collapsed_edges = 0;
  int num_eroded_edges_cur;
  int num_stable_trees = 0;

  /* locate edges that connect an internal vertex to an external vertex */
  it_edge = it_e_new(G, (E_TYPE) (E_H | E_S));

  while ((e = it_e_next(it_edge)) != (EDGE *) 0) {
    if (e->begin->ext_flag != e->end->ext_flag) {

#ifdef DEBUG_OUT
      fprintf(stderr,"\nDEBUG: Erosion: examining tree on edge %lx\n", e);
#endif

      num_eroded_edges_cur = erode_tree(G, params, e);
      if (num_eroded_edges_cur != 0) {
	if (num_eroded_edges < 0) {
	  /* collapsed a tree onto its largest path */
	  num_collapsed_edges += (-2)*num_eroded_edges_cur;
	  num_collapsed_trees += 2;
	} else {
	  /* deleted the whole tree */
	  num_eroded_edges += 2*num_eroded_edges_cur;
	  num_eroded_trees += 2;
	}
      } else {
	num_stable_trees += 2;
      }
    }
  }

  it_e_destroy(it_edge);

  /***********************************************************************
   * statistics
   ***********************************************************************/

  printf("----------------------------------------------------------------------------------------\n");
  printf("Erosion:\n");
  printf("Removed %d hanging trees with %d edges\n",
	 num_eroded_trees, num_eroded_edges);
  printf("Collapsed %d hanging trees to their max paths, removing %d edges\n",
	 num_collapsed_trees, num_collapsed_edges);
  printf("Left %d trees intact\n",num_stable_trees);

#ifdef DEBUG_OUT
  fprintf(stderr,"DEBUG: Exiting erode_external_undir\n");
#endif
}


void debug_print_edges(NODES *v, E_TYPE e_type)
{
  IVEDGE it_edge_mem, *it_edge;
  EDGE *e;

  it_edge = it_ev_renew(&it_edge_mem, v, e_type);
  while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
    fprintf(stderr," %lx", e);
  }
  fprintf(stderr,"\n");
}

/* erode the tree off of edge e
 * If the tree is deleted and contains k edges, return k
 * If the tree is collapsed onto its largest path and k edges are deleted,
 * return -k
 * If nothing changes, return 0
 *
 * Also erode the complementary tree
 */
int erode_tree(IGRAPH *G, CLEAN_PARAMS *params, EDGE *e)
{
  NODES *v1;            /* internal node */
  NODES *v2;            /* external node */
  NODES *v_farthest;

  int collapse_tree = 0;
  int num_edges_deleted;

  if (e->begin->ext_flag) {
    v1 = e->end;
    v2 = e->begin;
  } else {
    v1 = e->begin;
    v2 = e->end;
  }

#ifdef DEBUG_OUT
  fprintf(stderr,"DEBUG: erode_tree, e=%lx, v1=%lx %s v2=%lx, d(v1)=(%d,%d), d(v2=(%d,%d), ext(v1,v2)=(%d,%d), vval=(%d,%d)\n",
	  e, v1,
	  (e->begin->ext_flag) ? "<-" : "->",
	  v2,
	  gdegree(v1,E_IN_H), gdegree(v1,E_OUT_H),
	  gdegree(v2,E_IN_H), gdegree(v2,E_OUT_H),
	  v1->ext_flag,v2->ext_flag,v1->visit_value,v2->visit_value);
  fprintf(stderr,"DEBUG: v1_in: "); debug_print_edges(v1,E_IN_H);
  fprintf(stderr,"DEBUG: v1_out: "); debug_print_edges(v1,E_OUT_H);
  fprintf(stderr,"DEBUG: v2_in: "); debug_print_edges(v2,E_IN_H);
  fprintf(stderr,"DEBUG: v2_out: "); debug_print_edges(v2,E_OUT_H);
#endif

  if (v2->visit_value + (e->length-1) > params->T_e->depth) {
    /* longest path is at least the threshold, condense the tree onto
     * that path
     */
    collapse_tree = 1; 
    mark_path(G, v2);
  }


  /* delete the tree, except for the longest path */
#ifdef DEBUG_OUT
  fprintf(stderr,"DEBUG: calling erode_collapse_tree\n");
#endif
  num_edges_deleted = erode_collapse_tree(v2,e);
  if (collapse_tree) num_edges_deleted = -num_edges_deleted;
#ifdef DEBUG_OUT
  fprintf(stderr,"DEBUG: exiting erode_tree\n");
#endif

  return num_edges_deleted;
}

/* on the external tree rooted at v,
 * mark the vertices on any longest path as "internal"
 */

void mark_path(IGRAPH *G, NODES *v)
{
  IVEDGE it_edge_mem, *it_edge;
  EDGE *e;
  NODES *v_next, *v2;

  /* examine the vertex labels
   * A path on which they decrease by the edge length on each edge
   * will be a max length path
   */

#ifdef DEBUG_OUT
  fprintf(stderr,"DEBUG: mark_path: ");
#endif
  while (v != (NODES *) 0) {
#ifdef DEBUG_OUT
    fprintf(stderr,"%lx(%d) ", v, v->visit_value);
#endif
    v->ext_flag = 0;
    v_next = (NODES *) 0;
    it_edge = it_ev_renew(&it_edge_mem, v, E_U_H);
    while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
      v2 = (e->begin == v) ? e->end : e->begin;
      if ((v2->visit_value + (e->length-1)) == v->visit_value) {
	v_next = v2;
	break;
      }
    }
    v = v_next;
  }

#ifdef DEBUG_OUT
  fprintf(stderr,"\n");
#endif
}



/* on the tree that sticks out from v away from e,
 * keep internal vertices, delete external vertices and edges on them.
 *
 * return number of deleted edges
 *
 *                ...
 *              /
 * K ----e---- v- ...
 *              \ ...
 *
 * K represents portion of H that v is connected to by e
 *
 * Also performs same operations on complementary tree
 * It's provable that the tree and its complement are disjoint
 *
 */

int erode_collapse_tree(NODES *v, EDGE *e)
{
  NODES *v_next;
  EDGE *e_next;
  IVEDGE it_edge_mem, *it_edge;

  int num_edges = 0;

#ifdef DEBUG_OUT
  int debug_deg;
#endif

#ifdef DEBUG_OUT
  fprintf(stderr,"DEBUG:   erode_collapse_tree  v=%lx e=%lx ext=%d vval=%d\n",
	  v, e, v->ext_flag, v->visit_value);
#endif

  /* delete/collapse subtrees */
  it_edge = it_ev_renew(&it_edge_mem, v, E_U_H);

  while ((e_next = it_ev_next(it_edge)) != (EDGE *) 0) {
    /* don't backtrack */
    if (e_next == e)
      continue;

    /* next vertex */
    v_next = (e_next->begin == v) ? e_next->end : e_next->begin;

    /* recursively collapse/delete tree */
    num_edges += erode_collapse_tree(v_next, e_next);
  }

  /* finally, remove e and v if appropriate */
  if (v->ext_flag) {
    num_edges++;

#ifdef DEBUG_OUT
    debug_deg = gdegree(v,E_U_H);
    if (debug_deg != 1) {
      fprintf(stderr,"***** erode_collapse_tree: v=%lx, e=%lx, has degree %d<>1\n",
	      v, e, debug_deg);
    }
#endif

    e_next = e->bal_edge;
    v_next = v->bal_node;
    subgraph_del_edge(e);
    subgraph_del_edge(e_next);
    subgraph_del_node(v);
    subgraph_del_node(v_next);
  }
  return num_edges;
}


/*****************************************************************************
 * collapse external components onto one path (or small # of paths)
 * emanating from the core
 *****************************************************************************/

void erode_external_dir(IGRAPH *G, LABELVERT *T_e)
{
  EDGE *e;
  IEDGE *it_edge;
  int depth = T_e->depth;
  int e_depth = 2*depth + 2;

#if 0
  IVERT *it_vert;
#endif
  NODES *v1, *v2;
  int d1, d2, s;

  /* statistics */

  int num_vert_del = 0;
  int num_edge_del = 0;
  int num_good_entries = 0;  /* entrances into core that are kept */
  int num_bad_entries = 0;   /* entrances into core that are severed */

  int num_init_internal = 0; /* number of edges initially internal */
  int num_init_external = 0; /* number of remaining edges */

  int num_border_external = 0; /* number of edges borderline external */
  int num_hump_internal = 0; /* neither end is internal but inside is */

  int nedges_sym;


  /**********************************************************************
   * first mark all edges that will be retained as "internal"
   **********************************************************************/

  /* locate internal edges */
  it_edge = it_e_new(G, (E_TYPE) (E_H | E_S));

  while ((e = it_e_next(it_edge)) != (EDGE *) 0) {
    v1 = e->begin;
    v2 = e->end;
    d1 = v1->visit_value;
    d2 = v2->visit_value;

#ifdef DEBUG_OUT
    fflush(stdout);
    fprintf(stderr,"Checking edge e=%lx v1=%lx->v2=%lx d1=%d d2=%d\n",
	    e, v1, v2, d1, d2);
#endif

    s = d1 + d2 + e->length - 1;

    /* is {e, e^c} one edge or two?
     * one edge only happens with some loops.
     */
    //    nedges_sym = (v1 == v2->bal_node && v2 == v1->bal_node) ? 1 : 2;
    nedges_sym = (e == e->bal_edge) ? 1 : 2;

    if (s < e_depth) {
      num_init_external += nedges_sym;

      //      printf("external: e=%lx d1=%d d2=%d len=%d\n", e, d1, d2, e->length - 1);

      if (d1 > depth || d2 > depth)
	num_bad_entries += nedges_sym;

    } else {
      /* edge starts out internal (> e_depth)
       * or nucleotide expansion culminates at internal vertex (== e_depth)
       *
       */
      if (s > e_depth) {
	num_init_internal += nedges_sym;
	if (d1 <= depth && d2 <= depth)
	  num_hump_internal += nedges_sym;
      } else {
	num_border_external += nedges_sym;

	/* change edge to internal for steepest descent path */
	e->ext_flag = e->bal_edge->ext_flag = 0;
      }

      /* mark steepest descent path(s) */
      if (d1 <= depth) {
	erode_edge_dir(G, T_e, v1);
	num_good_entries += nedges_sym;
      }

      if (d2 <= depth) {
	erode_edge_dir(G, T_e, v2);
	num_good_entries += nedges_sym;
      }
    }
  }

  /**********************************************************************
   * delete all external edges and vertices
   * same treatment will occur to each edge/node and its complement,
   * so they are not done concurrently
   **********************************************************************/

#ifdef DEBUG_OUT
  fflush(stdout);
  fprintf(stderr,"Deleting external edges\n");
#endif

  it_edge = it_e_renew(it_edge, G, E_H);
  while ((e = it_e_next(it_edge)) != (EDGE *) 0) {
    if (e->ext_flag) {
      subgraph_del_edge(e);
      num_edge_del++;
    }
  }
  it_e_destroy(it_edge);

  /* isolated vertices */
#ifdef DEBUG_OUT
  fprintf(stderr,"Deleting isolated vertices\n");
#endif
  num_vert_del = mark_deleted_nodes(G);

#if 0
  it_vert = it_v_new(G, V_A_H);
  while ((v1 = it_v_next(it_vert)) != (NODES *) 0) {
    if (gdegree_cmp(v1, E_U_H, 0) == 0) {
      subgraph_del_node(v1);
      num_vert_del++;
    }
  }
  it_v_destroy(it_vert);
#endif

  /**********************************************************************
   * report
   **********************************************************************/

  printf("----------------------------------------------------------------------------------------\n");
  printf("Erosion:\n");
  printf("Graph has %d internal edges and %d external edges.\n",
	 num_init_internal, num_init_external);
  printf("# Path leaders: %d internal, %d external.\n",
	 num_hump_internal, num_border_external);

  //  printf("Graph has %d internal edges (incl. %d hump edges)\n",
  //	 num_init_internal, num_hump_internal);
  //  printf("and %d external edges (incl. %d borderline).\n",
  //	 num_init_external, num_border_external);

  printf("Paths are kept at %d entrances to the core, severed from %d entrances.\n",
	 num_good_entries, num_bad_entries);
  printf("Delete %d edges and %d vertices.\n",
	 num_edge_del, num_vert_del);
  fflush(stdout);

}


/* mark a steepest descent path from vertex v
 */
void erode_edge_dir(IGRAPH *G, LABELVERT *T_e, NODES *v)
{
  EDGE *e;
  IVEDGE it_edge_mem, *it_edge=&it_edge_mem;
  int depth = T_e->depth;

  int d1, d2;
  NODES *w;

#ifdef DEBUG_OUT
  fflush(stdout);
  fprintf(stderr,"erode_edge_dir: steepest descent from v=%lx, label=%d\n",
	  v, v->visit_value);
#endif

  while (v != (NODES *) 0 && v->visit_value > 1) {


    /* reached an internal vertex, hence have
     * joined into a pre-existing marked path.
     */
#ifdef DEBUG_OUT
    fprintf(stderr,"v=%lx label=%d ext_flag=%d\n", v, v->visit_value, v->ext_flag);
#endif

    if (!v->ext_flag) {
#ifdef DEBUG_OUT
      fprintf(stderr,"already internal, aborting\n");
#endif
      return;
    }

    /* set node & complement to internal */
    v->ext_flag = v->bal_node->ext_flag = 0;

    d1 = v->visit_value;

    it_edge = it_ev_renew(it_edge, v, E_U_H);
    while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
      w = (e->begin == v) ? e->end : e->begin;

      d2 = w->visit_value;

#ifdef DEBUG_OUT
      //      if (v == (NODES *) 0x157c772a0
      //	  || v == (NODES *) 0x15bbd7c60
      //	  || v->bal_node == (NODES *) 0x157c772a0
      //	  || v->bal_node == (NODES *) 0x15bbd7c60)
      fprintf(stderr, "examining e=%lx, v=%lx, w=%lx, label(v)=%d, label(w)=%d, len=%d\n",
	      e, v, w, d1, d2, e->length - 1);
#endif

      /* steepest descent   v -> w   or   v <- w
       * requires label(v)-label(w) = length of edge
       */

      if (d1 - d2 == e->length - 1)
	break;
    }

#ifdef DEBUG_OUT
    fprintf(stderr, "chose e=%lx, w=%lx\n", e, w);
#endif

    if (e == (EDGE *) 0) {
      /* shouldn't happen! */
      fprintf(stderr, "erode_external_dir: cannot find steepest descent path. v=%lx\n", v);
      return;
    }

    /* set edge & complement to internal */
    e->ext_flag = e->bal_edge->ext_flag = 0;

    v = w;
  }

  /* set final node & complement to internal */
  v->ext_flag = v->bal_node->ext_flag = 0;

#ifdef DEBUG_OUT
  fprintf(stderr, "path terminated at v=%lx\n", v);
#endif
}

/*****************************************************************************
 * delete components that only contain external vertices
 *
 * some components may be self-complementary, but no special
 * care is required, since everything is deleted eventually
 *****************************************************************************/

void erode_small_components(IGRAPH *G, CLEAN_PARAMS *params)
{
  /* TODO:
   * if a more efficient way of iterating over components is devised,
   * redo this code.
   * Since iterating over components currently involves iterating over
   * all nodes, we will just iterate over all nodes here.
   */

  IVERT *it_vert;
  NODES *v;
  IVEDGE it_edge_mem, *it_edge;
  EDGE *e;

  int num_deleted_nodes = 0;
  int num_deleted_edges = 0;
  int num_deleted_components = 0; 

  classify_internal_components(G);

  it_vert = it_v_new(G, V_A_H);
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    if (comp_id(v)->comp_ext_flag) {
      /* delete all outgoing edges on node */
      it_edge = it_ev_renew(&it_edge_mem, v, E_OUT_H);
      while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
	subgraph_del_edge(e);
	num_deleted_edges++;
      }

      /* delete node */
      subgraph_del_node(v);
      num_deleted_nodes++;

      /* component count */
      if (v->comp == v)
	num_deleted_components++;
    }
  }

  printf("Deleted %d external components, containing %d nodes and %d edges\n",
	 num_deleted_components, num_deleted_nodes, num_deleted_edges);
}

/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * label_external(G, T_e)
 *
 * undirected:
 * For every vertex of depth <= threshold from a leaf,
 * set its visit_value = its depth
 *
 * All other vertices are labeled theshold + 1
 *
 *
 * directed:
 * similar but depth is from sources/sinks instead of leaves
 *
 ****************************************************************************/

void label_external(IGRAPH *G, LABELVERT *T_e)
{

#ifdef DEBUG_OUT
  fprintf(stderr,"label_external: depth=%d dir=%d\n", T_e->depth, T_e->dir);
#endif

  if (T_e->dir) {
    label_external_dir(G, T_e);
  } else {
    label_external_undir(G, T_e);
  }
}

void label_external_undir(IGRAPH *G, LABELVERT *T_e)
{
  IVERT *it_vert;
  NODES *v;             /* current leaf */
  int inf = T_e->depth + 1;

  NODES *w;             /* current vertex */
  NODES *w1, *w2;       /* w1=neighbor vertex w/highest label,
			 * w2=2nd highest
			 */
  NODES *w_next;        /* neighbor currently being explored */
  int d1, d2, d_next;
  EDGE *e1,*e2;

  IVEDGE it_ev_mem;
  IVEDGE *it_edge;
  EDGE *e;


  //  fprintf(stderr,"DEBUG: entering label_external\n");

  /* set all vertices to high depth */
  it_vert = it_v_new(G, V_A_H);
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    v->visit_value = inf;
  }

  /* follow paths in from leaves */
  it_vert = it_v_new(G, (V_TYPE) (V_L_UH | V_S));
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    //    fprintf(stderr,"DEBUG: exploring paths from leaf v=%lx\n", v);
    
    /* follow path from v into graph */
    w = v;
    w->visit_value = 1;
    w->bal_node->visit_value = 1;

    /* get unique edge, if any, on w */
    e = get_first_edge(v, E_U_H);
    if (e == (EDGE *) 0)
      continue;

    /* other vertex on that edge */
    w = (e->begin == v) ? e->end : e->begin;

    /* follow a path into the graph */
    do {
      //      fprintf(stderr,"DEBUG:   w=%lx, label=%d -> ", w, w->visit_value);
      /* compute label(w) based on its top two neighbors
       */
      w1 = w2 = (NODES *) 0;
      d1 = d2 = -1;
      e1 = e2 = (EDGE *) 0;

      it_edge = it_ev_renew(&it_ev_mem, w, E_U_H);
      while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
	w_next = (e->begin == w) ? e->end : e->begin;
	if (w_next == w) {
	  /* vertex loop, hence vertex is in core, abort exploration */
	  d2 = d1 = -1;
	  break;
	}
	d_next = w_next->visit_value + (e->length - 1);

	/* have we displaced one of the top two neighbor labels? */
	if (d_next > d2) {
	  /* yes, update list of top two neighbors */
	  if (d_next > d1) {
	    e2 = e1; 
	    w2 = w1;
	    d2 = d1;
	    e1 = e;
	    w1 = w_next;
	    d1 = d_next;
	  } else {
	    e2 = e;
	    w2 = w_next;
	    d2 = d_next;
	  }
	}
      }

      if (d2 > 0 && d2 < w->visit_value) {
	w->visit_value = d2;     /* update label of w */
	w->bal_node->visit_value = d2;

	if (d2 + (e1->length - 1) < w1->visit_value) {
	  /* if next node label will change, move there */
	  w = w1;
	} else {
	  w = (NODES *) 0;         /* done */
	}
	//      } else if (e1 != (EDGE *) 0 &&
	//		 w->visit_value + (e1->length - 1) < w1->visit_value) {
	//	/* if next node label will change, move there */
	//	w = w1;
      } else {
	w = (NODES *) 0;         /* done */
      }

#if 0
      if (d2 > 0 && d2 < w->visit_value) {
	w->visit_value = d2;     /* update label of w */
	w->bal_node->visit_value = d2;
	//	fprintf(stderr,"%d\n", d2); /* DEBUG */

	if (d1 > d2) {
	  /* advance to next vertex */
	  w = w1;
	} else {
	  /* stop */
	  w = (NODES *) 0;
	}
      } else if (d2 > 0 && d2 == w->visit_value && d1 > d2) {
	  /* advance to next vertex */
	  w = w1;
      } else {
	//	fprintf(stderr,"done\n"); /* DEBUG */
	w = (NODES *) 0;         /* done */
      }
#endif
      
    } while (w != (NODES *) 0);
  }

  //  fprintf(stderr,"DEBUG: exiting label_external\n");

}

void label_external_dir(IGRAPH *G, LABELVERT *T_e)
{
  IVERT *it_vert;
  NODES *v;             /* current leaf */
  int inf = T_e->depth + 1;

  //  fprintf(stderr,"DEBUG: entering label_external_dir\n");

  /* initialize most vertices to depth infinity,
   * sources/sinks to depth _2_ so that if they are encountered
   * early, the branches will be cut off (later they will be set to 1)
   */
  it_vert = it_v_new(G, V_A_H);
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    if (gdegree_cmp(v,E_IN_H,0) == 0 ||
	gdegree_cmp(v,E_OUT_H,0) == 0)
      v->visit_value = 2;
    else
      v->visit_value = inf;
  }

  /* follow paths in from sources/sinks */
  it_vert = it_v_renew(it_vert, G, (V_TYPE) ( V_SOI_H | V_S));
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    //    v->visit_value = v->bal_node->visit_value = 1;
    label_external_dir_v(T_e, v);
  }
  it_v_destroy(it_vert);

  //  fprintf(stderr,"DEBUG: exiting label_external_dir\n");
}

/*
 */
void label_external_dir_v(LABELVERT *T_e, NODES *v)
{
  IVEDGE *it_edge;
  EDGE *e;

  int max_in = 1;      /* label based on incoming neighbors */
  int max_out = 1;     /* label based on outgoing neighbors */
  int new_label;       /* actual new label */

  /* maximum based on incoming neighbors */
  it_edge = it_ev_new(v, E_IN_H);
  while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
    new_label = e->begin->visit_value + e->length - 1;
    if (new_label > max_in)
      max_in = new_label;
  }

  /* maximum based on outgoing neighbors */
  it_edge = it_ev_renew(it_edge, v, E_OUT_H);
  while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
    new_label = e->end->visit_value + e->length - 1;
    if (new_label > max_out)
      max_out = new_label;
  }

  new_label = min(max_in,max_out);

  if (new_label < v->visit_value) {
    /* only assign new label if it decreased */
    v->visit_value = v->bal_node->visit_value = new_label;

    /* explore neighbors that might decrease too */
    if (new_label < max_in) {
      /* only check incoming edges if some of them might decrease */
      it_edge = it_ev_renew(it_edge, v, E_IN_H);
      while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
	if (new_label + e->length - 1 < e->begin->visit_value) {
	  label_external_dir_v(T_e, e->begin);
	}
      }
    }

    if (new_label < max_out) {
      /* only check outgoing edges if some of them might decrease */
      it_edge = it_ev_renew(it_edge, v, E_OUT_H);
      while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
	if (new_label + e->length - 1 < e->end->visit_value) {
	  label_external_dir_v(T_e, e->end);
	}
      }
    }
  }

  it_ev_destroy(it_edge);
}


/****************************************************************************
 * classify_internal_nodes(G, T_e)
 *
 * Sets the internal/external flag for every node and edge
 ****************************************************************************/

void classify_internal_nodes(IGRAPH *G, LABELVERT *T_e)
{
  IVERT *it_vert;
  NODES *v;
  int depth = T_e->depth;

  IEDGE *it_edge;
  EDGE *e;
  int e_depth1 = 2*depth + 3;

  /* label vertices */
  label_external(G, T_e);

  /* mark internal/external vertices */
  it_vert = it_v_new(G, V_A_H);
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    v->ext_flag = (int) v->visit_value <= depth;
  }
  it_v_destroy(it_vert);

  /* mark internal/external edges */
  it_edge = it_e_new(G, E_G);   /* want E_H + E_P */
  while ((e = it_e_next(it_edge)) != (EDGE *) 0) {
    /* should be k1 + k2 + e->length - 1 <= e_depth,
     * but absorb the -1 into e_depth
     */
    e->ext_flag =
      ((int)e->begin->visit_value + (int)e->end->visit_value + e->length)
      <= e_depth1;
  }
  it_e_destroy(it_edge);

}

/****************************************************************************
 * classify_internal_components(G)
 *
 * Sets the internal/external flag for every component
 * Assumes classify_internal_nodes has been done
 * The flag is "internal" if the component has >=1 internal nodes
 * and is external otherwise
 ****************************************************************************/

void classify_internal_components(IGRAPH *G)
{
  ICOMP *it_comp;
  COMP *comp;

  IVERT *it_vert;
  NODES *v;

  /* set all components to "external" */
  it_comp = it_c_new(G, C_A_H);
  while ((comp = it_c_next(it_comp)) != (COMP *) 0) {
    comp->comp_ext_flag = 1;
  }
  it_c_destroy(it_comp);

  /* determine if they have any internal nodes */
  it_vert = it_v_new(G, V_A_H);
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    if (!v->ext_flag) {
      comp_id(v)->comp_ext_flag = 0;
    }
  }

  it_v_destroy(it_vert);
}

/****************************************************************************
 * Determine whether a path is "internal" or "external".
 * One end of the path is defined by a starting vertex v;
 * the other end is determined by following edges in direction e_type
 * until any of these are reached:
 *      1. A branching vertex (or the original vertex)
 *      2. A high coverage edge (stop w/o entering it)
 *      3. An internal vertex or edge
 *
 * The path is "internal" if it stops at an internal vertex or edge
 *
 * Special case: if v is a branching vertex, the path only consists of v.
 *
 * Note: the reverse edge out of e is not yet in the graph, but we
 * are considering it as a candidate edge, so its degree must be counted.
 *
 * TODO: determine if we need a cutoff for the path length as well.
 *
 * input:
 *   v:       starting vertex
 *   e_type:  edge type to follow
 *   L_c:     low coverage edge is <= L_c,  high coverage is > L_c
 ****************************************************************************/

int is_path_internal(NODES *v, int e_type, int L_c)
{
  NODES *w;
  EDGE *e;
  int dir = (e_type & E_DIR_MASK) == E_IN_G;

#if 0
  /* DEBUG: emulate original behavior */
  return !v->ext_flag;
#endif

  /* 0 instead of 1 because we provisionally add in
   * one of the reverse edges from v
   */
  if (gdegree_cmp(v, (E_TYPE) (e_type ^ E_REV_MASK), 0) > 0) {
    /* branching vertex, path is just the vertex */
    return !v->ext_flag;
  }

  w = v;

  /* if hit an internal vertex on the path, immediately return that it's
   * internal
   */
  while (w->ext_flag) {
    /* reverse direction has degree 1
     * direction e_type has unknown degree d
     * if d != 1, path is over
     * if d=1 and w is internal, then path is internal
     * else d=1 and w is external, so extend path
     */

    e = get_unique_edge(w, e_type);

    /* reached branching vertex or self-loop,
     * or coverage is 
     * terminate */
    if (e == (EDGE *) 0             /* reached source/sink/branching vertex */
	|| e->begin == e->end       /* self-loop, returned to start */
#if 0
	|| e->multip > L_c          /* high coverage edge */
#endif
	)
      return !w->ext_flag;

    /* if edge is internal, return */
    if (!e->ext_flag)
      return 1;

    /* advance to next vertex */
    w = dir ? e->begin : e->end;

    /* path is a cycle, terminate */
    if (w == v)
      return !w->ext_flag;
  }

  /* reached internal vertex, terminate */
  return 1;

}


/****************************************************************************
 * Chimeric edge test
 ****************************************************************************/

int is_edge_chimeric(EDGE *e, int L_c)
{
  return !e->ext_flag
    || (is_path_internal(e->begin, E_IN_H, L_c)
	&& is_path_internal(e->end, E_OUT_H, L_c));
}

