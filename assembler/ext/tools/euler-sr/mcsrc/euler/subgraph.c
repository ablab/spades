/***************************************************************************
 * Title:          subgraph.c
 * Author:         Glenn Tesler
 * Created:        Dec. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/* graph iterators
 * Glenn Tesler
 * 12/18/02
 */

#include <stdinc.h>
#include <extvab.h>
//#include <param.h>
#include <extfunc.h>
//#include <subgraph.h>


/* private */

//int explore_path(IGRAPH *G,
//		 NODES *v,
//		 NODES *goal,
//		 EDGE *prev_edge,
//		 int L,
//		 int cur_dist);




/* General routines for the combined graph/subgraph structure
 * Notation:
 * G = graph (directed)
 * H = subgraph (directed)
 * UG, UH = undirected versions
 *
 * vertex sets:   X = G or H
 * V_A_X = all
 * V_L_UX = leaves (total deg = 0 or 1)
 * V_B_UX = branching vertices (total deg != 2), incl. leaves
 * V_SB_UX = strict branching (total deg > 2)
 * V_B_X = branching (not 1-in/1-out)
 * V_BC_X = branching and cycle vertices (in two passes)
 * V_SO_X = sources
 * V_SI_X = sinks
 * V_SOI_X = sources or sinks (indegree=0 or outdegree=0)
 *
 * can OR (|) the above with V_S to only include one of v, v^c in iteration
 *
 * edge sets: for iterating on edges at a vertex, and for gdegree:
 * E_OUT_G = all outgoing edges of G
 * E_OUT_H = all outgoing edges of H
 * E_OUT_P = all outgoing postponed edges
 * E_IN_G, E_IN_H, E_IN_P: same for incoming edges
 * E_U_G, E_U_H, E_U_P: incoming and outgoing
 *
 * for iterating over all edges in a graph:
 * E_G, E_H, E_P  (synonyms for OUT variants above)
 * can OR (|) these with E_S to only include one of e, e^c in iteration
 */


/****************************************************************************
 * Iterators
 *
 * Loop over all vertices of some type
 *    it = it_v_new(G,v_type)       // create iterator
 *    it = it_v_renew(it,G,v_type)  // same but reuse memory from
 *                                  // an existing "it" structure
 *    v = it_v_next(it)             // return ptr to next vertex
 *                                  // returns NULL at end
 *    it_v_destroy(it)              // destroy iterator
 *    it_v_count(G,v_type)          // returns count of # vertices of this type
 *
 *
 * Loop over all edges of some type in (sub)graph
 *    it = it_e_new(G,e_type)
 *    it = it_e_renew(it,G,e_type)
 *    e = it_e_next(it)
 *    it_e_destroy(it)
 *    it_v_count(G,e_type)
 *    
 * Loop over all edges of some type incident with a particular vertex
 *    it = it_ev_new(v,e_type)
 *    it = it_ev_renew(it,v,e_type)
 *    e = it_ev_next(it)
 *    it_ev_destroy(it)
 *    it_ev_count(v,e_type)
 *
 * Loop over all components of H
 *    it = it_c_new(G,c_type)
 *    it = it_c_renew(it,G,c_type)
 *    comp = it_c_next(it)
 *    it_c_destroy(it)
 *    it_c_count(G,c_type)
 *
 ****************************************************************************/

/****************************************************************************
 * iterate over all vertices
 * it_v_*
 ****************************************************************************/

IVERT *it_v_renew(IVERT *it,
		  IGRAPH *G, V_TYPE v_type)
{
  it->G = G;
  it->v_type = v_type & ~V_S;
  it->sym = (v_type & V_S) != 0;

  /* The next vertex will be 0, i.e., the first one */
  it->v_num = -1;

  it->pass = 0;

  if (it->v_type == V_BC_G) {
    clear_visit_flag(G, V_A_G);
    it->pass = 1;
    it->v_type = V_B_G;
  } else if (it->v_type == V_BC_H) {
    clear_visit_flag(G, V_A_H);
    it->pass = 1;
    it->v_type = V_B_H;
  }

  //  fprintf(stderr,"v_type = %d\n", v_type);

  return it;
}

IVERT *it_v_new(IGRAPH *G, V_TYPE v_type)
{
  return it_v_renew((IVERT *) ckalloc(1 * sizeof(IVERT)),
		    G, v_type);
}

void it_v_destroy(IVERT *it)
{
  free ((void *) it);
}

int it_v_count(IGRAPH *G, V_TYPE v_type)
{
  int count = 0;
  IVERT it_v_mem;
  IVERT *it = it_v_renew(&it_v_mem, G, v_type);

  while (it_v_next(it) != (NODES *) 0) {
    count++;
  }

  return count;
}


NODES *it_v_next(IVERT *it)
{
  IGRAPH *G = it->G;
  NODES *v;

  IVEDGE it_edge_mem, *it_edge;
  EDGE *e;
  NODES *v1;


  while (++it->v_num < G->num_nodes) {
    v = G->nodes[it->v_num];
    if (it->v_type >= V_A_H) {
      //      fprintf(stderr,"v=%lx: indegree=%d, outdegree=%d\n",
      //	      v, gdegree(v,E_IN_H), gdegree(v,E_OUT_H));
    }

    /* if only doing one of v, v^c, then do the one with lower pointer value */
    if (it->sym &&
	v > v->bal_node)
      continue;

    /* check if the vertex has suitable type */
    /* TODO: specialized versions of gdegree that
     * will check each degree condition more efficiently
     */
    switch (it->v_type) {
    case V_A_G:
      break;
    case V_L_UG:
      if (gdegree_cmp(v,E_U_G,1) <= 1)
	break;
      continue;  /* not a leaf, continue "while" */
    case V_B_UG:
      if (gdegree_cmp(v,E_U_G,2) != 2)
	break;
      continue;  /* not branching, continue "while" */
    case V_SB_UG:
      if (gdegree_cmp(v,E_U_G,2) > 2)
	break;
      continue;  /* not strict branching, continue "while" */
    case V_B_G:
      if (gdegree_cmp(v,E_IN_G,1) != 1 ||
	  gdegree_cmp(v,E_OUT_G,1) != 1)
	break;
      continue;  /* 1-in/1-out, continue "while" */
    case V_SO_G:
      if (gdegree_cmp(v,E_IN_G,0) == 0 &&
	  gdegree_cmp(v,E_OUT_G,1) == 1)
	break;
      continue;  /* not source, continue "while" */
    case V_SI_G:
      if (gdegree_cmp(v,E_OUT_G,0) == 0 &&
	  gdegree_cmp(v,E_IN_G,1) == 1)
	break;
      continue;  /* not sink, continue "while" */
    case V_SOI_G:
      if (gdegree_cmp(v,E_OUT_G,0) == 0 ||
	  gdegree_cmp(v,E_IN_G,0) == 0)
	break;
      continue;  /* not source or sink, continue "while" */

    case V_A_H:
      if (v->subg_flag)
	break;
      continue;  /* vertex not in subgraph, continue "while" */
    case V_L_UH:
      if (v->subg_flag && gdegree_cmp(v,E_U_H,1) <= 1)
	break;
      continue;  /* not a leaf, continue "while" */
    case V_B_UH:
      if (v->subg_flag && gdegree_cmp(v,E_U_H,2) != 2)
	break;
      continue;  /* not branching, continue "while" */
    case V_SB_UH:
      if (v->subg_flag && gdegree_cmp(v,E_U_H,2) > 2)
	break;
      continue;  /* not strict branching, continue "while" */
    case V_B_H:
      if (v->subg_flag &&
	  (gdegree_cmp(v,E_IN_H,1) != 1 ||
	   gdegree_cmp(v,E_OUT_H,1) != 1)
	  )
	break;
      continue;  /* 1-in/1-out, continue "while" */
    case V_SO_H:
      if (gdegree_cmp(v,E_IN_H,0) == 0 &&
	  gdegree_cmp(v,E_OUT_H,1) == 1)
	break;
      continue;  /* not source, continue "while" */
    case V_SI_H:
      if (gdegree_cmp(v,E_OUT_H,0) == 0 &&
	  gdegree_cmp(v,E_IN_H,1) == 1)
	break;
      continue;  /* not sink, continue "while" */
    case V_SOI_H:
      if (gdegree_cmp(v,E_OUT_H,0) == 0 ||
	  gdegree_cmp(v,E_IN_H,0) == 0)
	break;
      continue;  /* not source or sink, continue "while" */
    }

    
    /* iterator V_BC_X? */
    if (it->pass > 0) {
      /* ignore vertices already visited */
      if (v->visit)
	continue;

      /* mark vertex as visited */
      v->visit = 1;
      if (it->pass == 1) {
	/* Follow each outgoing edge from v1 until another branching vertex
	 * is reached.  Mark the vertices on the path as visited.
	 */
	v1 = v;
	it_edge = it_ev_renew(&it_edge_mem, v1, E_OUT_H);
	while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {

	  /* follow path that starts on edge e, terminating at
	   * a branching vertex.
	   * Mark the intermediate vertices as visited.
	   */
	  while (1) {
	    v1 = e->end;
	    e = find_unique_oedge(v1);
	    if (e == (EDGE *) 0) break;
	    v1->visit = 1;
	  }
	}
      }
    }

    /* vertex type is acceptable, return pointer to the vertex */


    if (it->v_type >= V_A_H) {
      //      fprintf(stderr,".");
    }
    return v;
  }

  /* exhausted all vertices */
  /* in case of multiple calls, adjust this so it will always point
   * to the end in the same manner
   */
  -- it->v_num;

  if (it->pass > 0) {
    /* iterator V_BC_X, check if need another pass */
    it->pass++;
    if (it->pass == 2) {
      /* start pass for cycle vertices */
      it->v_type = (it->v_type == V_B_G) ? V_A_G : V_A_H;
      it->v_num = -1;
      return it_v_next(it);
    }
  }

  /* indicate end of iterations */
  return (NODES *) 0;
}


/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * iterate over all edges in (sub)graph
 * it_e_*
 ****************************************************************************/

/* careful with renew:
 * it assumes the vert_it field points to a valid IVERT that can be reused
 */
IEDGE *it_e_renew(IEDGE *it,
		  IGRAPH *G, E_TYPE e_type)
{
  int v_type;
  int e_in;

  it->G = G;
  it->e_type = e_type & ~E_S;
  it->sym = (e_type & E_S) != 0;

  /* Each edge type requires iterating over only certain vertex types */
  v_type = ((it->e_type & E_TYPE_MASK) == E_OUT_H)
    ? V_A_H    /* e_type = E_*_H */
    : V_A_G;   /* e_type = E_*_G or E_*_P */
#if 0
  switch (it->e_type) {
  case E_OUT_G:
  case E_OUT_P:
  case E_IN_G:
  case E_IN_P:
  case E_U_G:
  case E_U_P:
    v_type = V_A_G; break;
  case E_OUT_H:
  case E_IN_H:
  case E_U_H:
    v_type = V_A_H; break;
  }
#endif

  it_v_renew(it->vert_it, G, (V_TYPE) v_type);

  /* no current vertex */
  it->v = (NODES *) 0;
  it->e_num = it->e_max = 0;

  return it;
}

IEDGE *it_e_new(IGRAPH *G, E_TYPE e_type)
{
  IEDGE *it;

  it = (IEDGE *) ckalloc(1 * sizeof(IEDGE));
  it->vert_it = (IVERT *) ckalloc(1 * sizeof(IVERT));

  return it_e_renew(it, G, e_type);
}

void it_e_destroy(IEDGE *it)
{
  it_v_destroy(it->vert_it);
  free ((void *) it);
}

int it_e_count(IGRAPH *G, E_TYPE e_type)
{
  int count = 0;
  IEDGE *it = it_e_new(G, e_type);

  while (it_e_next(it) != (EDGE *) 0) {
    count++;
  }

  it_e_destroy(it);

  return count;
}

/* NOTE:
 * it_e_next and it_ev_next both share code.
 * If one is changed, be sure to make the corresponding changes
 * in the other.
 * It's done inline instead of as a subroutine for efficiency.
 */
EDGE *it_e_next(IEDGE *it)
{
  IGRAPH *G = it->G;
  EDGE *e;
  NODES *v;

  /* loop through many vertices/edges until find suitable edge */
  while (1) {
    /* should we advance to next vertex? */
    if (++it->e_num >= it->e_max) {
      /* advance to next vertex */
      v = it->v = it_v_next(it->vert_it);

      if (v == (NODES *) 0) {
	/* iterator done */
	it->e_num = it->e_max = 0;
	return (EDGE *) 0;
      }

      /* set edge indices correctly for this vertex */
      switch (it->e_type & E_DIR_MASK) {
      case E_OUT_G:
	it->e_num = 0;
	it->e_max = v->num_nextedge;
	break;
      case E_IN_G:
	it->e_num = -v->num_lastedge;
	it->e_max = 0;
	break;
      case E_U_G:
	it->e_num = -v->num_lastedge;
	it->e_max = v->num_nextedge;
	break;
      }
      if (v->num_nextedge != 1) {
	//	fprintf(stderr,"v=%x, v_num=%d\n",v,it->vert_it->v_num);
	//	fprintf(stderr,"   e_index: [%d , %d)\n", it->e_num, it->e_max);
      }

      /* if vertex has no edges, advance to next vertex */
      if (it->e_num >= it->e_max)
	continue;
    } else {
      v = it->v;
    }


    /* the vertex is in G (or H)
     * and there are still edges to explore
     */

    if (it->e_num >= 0) {
      /* outgoing edge */
      e = v->nextedge[it->e_num];
    } else {
      /* incoming edge */
      e = v->lastedge[-1 - it->e_num];
    }

#if 0
    /* edge may have been deleted from G */
    if (e == (EDGE *) 0)
      continue;
#endif

    /* if only doing one of e, e^c, then do the one with lower pointer value */
    if (it->sym &&
	e > e->bal_edge) {
      //        fprintf(stderr,"rejecting e=%lx, mul %d, has balanced edge ec=%lx\n",
      //		e, e->multip, e->bal_edge);
      //	fprintf(stderr,"  node v=%lx, bal_node=%lx\n", v, v->bal_node);
      continue;
    }

    /* check edge type */
    switch (it->e_type & E_TYPE_MASK) {
    case E_OUT_G:
      break;
    case E_OUT_H:
      if (e->subg_flag == SUBG_IN)
	break;
      //      fprintf(stderr,"rejecting edge, not in H\n");
      continue;    /* edge not in H, continue "while" */
    case E_OUT_P:
      if (e->subg_flag == SUBG_POSTPONE)
	break;
      //      fprintf(stderr,"rejecting edge, not postponed\n");
      continue;    /* edge not postponed, continue "while" */
    }

    /* passed all tests, return the edge */
    //    fprintf(stderr,"using edge\n");
    return e;
  }
}



/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * iterate over all edges incident with a particular vertex
 * it_ev_*
 ****************************************************************************/

IVEDGE *it_ev_renew(IVEDGE *it,
		    NODES *v, E_TYPE e_type)
{
  it->v = v;
  it->e_type = e_type & ~E_S;
  it->sym = (e_type & E_S) != 0;

  /* set edge indices correctly for this vertex */
  switch (e_type & E_DIR_MASK) {
  case E_OUT_G:
    /* point to one before the edge, so that next advances to it */
    it->e_num = 0 - 1;
    it->e_max = v->num_nextedge;
    break;
  case E_IN_G:
    it->e_num = -v->num_lastedge - 1;
    it->e_max = 0;
    break;
  case E_U_G:
    it->e_num = -v->num_lastedge - 1;
    it->e_max = v->num_nextedge;
    break;
  }

  return it;
}

IVEDGE *it_ev_new(NODES *v, E_TYPE e_type)
{
  return it_ev_renew((IVEDGE *) ckalloc(1 * sizeof(IVEDGE)),
		     v, e_type);
}

void it_ev_destroy(IVEDGE *it)
{
  free ((void *) it);
}


int it_ev_count(NODES *v, E_TYPE e_type)
{
  int count = 0;
  IVEDGE it_ev_mem;
  IVEDGE *it = it_ev_renew(&it_ev_mem, v, e_type);

  while (it_ev_next(it) != (EDGE *) 0) {
    count++;
  }

  return count;
}


/* NOTE:
 * it_e_next and it_ev_next both share code.
 * If one is changed, be sure to make the corresponding changes
 * in the other.
 * It's done inline instead of as a subroutine for efficiency.
 */
EDGE *it_ev_next(IVEDGE *it)
{
  EDGE *e;
  NODES *v = it->v;

  while (1) {
    if (++it->e_num >= it->e_max) {
      /* exhausted all edges */
      /* in case of multiple calls, adjust this so it will always point
       * to the end in the same manner
       */
      -- it->e_num;

      /* indicate end of iterations */
      return (EDGE *) 0;
    }


    /* the vertex is in G (or H)
     * and there are still edges to explore
     */

    if (it->e_num >= 0) {
      /* outgoing edge */
      e = v->nextedge[it->e_num];
    } else {
      /* incoming edge */
      e = v->lastedge[-1 - it->e_num];
    }

#if 0
    /* edge may have been deleted from G */
    if (e == (EDGE *) 0)
      continue;
#endif

    /* if only doing one of e, e^c, then do the one with lower pointer value */
    if (it->sym &&
	e > e->bal_edge)
      continue;

    /* check edge type */
    switch (it->e_type & E_TYPE_MASK) {
    case E_OUT_G:
      break;
    case E_OUT_H:
      if (e->subg_flag == SUBG_IN)
	break;
      continue;    /* edge not in H, continue "while" */
    case E_OUT_P:
      if (e->subg_flag == SUBG_POSTPONE)
	break;
      continue;    /* edge not postponed, continue "while" */
    }

    /* passed all tests, return the edge */
    return e;
  }
}

/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * iterate over all components
 * it_c_*
 ****************************************************************************/

ICOMP *it_c_renew(ICOMP *it,
		  IGRAPH *G, C_TYPE c_type)
{
  it->G = G;
  it->c_type = c_type & ~C_S;
  it->sym = (c_type & C_S) != 0;

  it->vert_it = it_v_renew(it->vert_it, G, V_A_H);

  return it;
}

ICOMP *it_c_new(IGRAPH *G, C_TYPE c_type)
{
  ICOMP *it;

  it = (ICOMP *) ckalloc(1 * sizeof(ICOMP));
  it->vert_it = (IVERT *) ckalloc(1 * sizeof(IVERT));

  return it_c_renew(it, G, c_type);
}

void it_c_destroy(ICOMP *it)
{
  it_v_destroy(it->vert_it);
  free ((void *) it);
}

int it_c_count(IGRAPH *G, C_TYPE c_type)
{
  int count = 0;
  ICOMP *it = it_c_new(G, c_type);

  while (it_c_next(it) != (COMP *) 0) {
    count++;
  }

  it_c_destroy(it);

  return count;
}

COMP *it_c_next(ICOMP *it)
{
  NODES *v;

  while ((v = it_v_next(it->vert_it)) != (NODES *) 0) {

    /* node is not root node of a component */
    if (v != v->comp)
      continue;

    /* if only doing one of comp, comp^c, then do the one with lower value */
    if (it->sym &&
	v > comp_id(v->bal_node))
      continue;

    return (COMP *) v;
  }

  /* iterator done */
  return (COMP *) v;
}



/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * degree of a vertex, only counting edges of specified type
 * d = gdegree(v,e_type)
 *
 * partial degree: for testing if the degree is <, =, >  m  
 * d2 = gdegree_cmp(v,e_type,m)
 * The value of d2 is such that
 *     d2<m   iff  d<m
 *     d2==m  iff  d==m
 *     d2>m   iff  d>m
 * However, except in the middle case, there is no guarantee that d2=d.
 ****************************************************************************/

int gdegree(NODES *v, E_TYPE e_type)
{
  int d = 0;
  IVEDGE it_mem, *it;

  if ((e_type & E_TYPE_MASK) == E_OUT_G) {
    /* degree in G, don't need to iterate over edges */
    switch (e_type) {
    case E_OUT_G:
      return v->num_nextedge;
    case E_IN_G:
      return v->num_lastedge;
    case E_U_G:
      return v->num_nextedge + v->num_lastedge;
    default:
      fprintf(stderr, "gdegree: bad e_type=%d\n", e_type);
      exit(-1);
    }
  }

  it = it_ev_renew(&it_mem, v, e_type);
  while (it_ev_next(it) != (EDGE *) 0) {
    d++;
  }

  return d;
}

int gdegree_cmp(NODES *v, E_TYPE e_type, int m)
{
  int d = 0;
  IVEDGE it_mem, *it;

  if ((e_type & E_TYPE_MASK) == E_OUT_G) {
    /* degree in G, don't need to iterate over edges */
    switch (e_type) {
    case E_OUT_G:
      return v->num_nextedge;
    case E_IN_G:
      return v->num_lastedge;
    case E_U_G:
      return v->num_nextedge + v->num_lastedge;
    default:
      fprintf(stderr, "gdegree_cmp: bad e_type=%d\n", e_type);
      exit(-1);
    }
  }

  it = it_ev_renew(&it_mem, v, e_type);
  while (it_ev_next(it) != (EDGE *) 0) {
    d++;
    if (d > m)
      return d;
  }

  return d;
}


/****************************************************************************
 * get_first_edge(v, e_type)
 * returns any edge on v of type e_type
 * or NULL if none
 ****************************************************************************/

EDGE *get_first_edge(NODES *v, int e_type)
{
  IVEDGE it_edge_mem, *it_edge;

  it_edge = it_ev_renew(&it_edge_mem, v, (E_TYPE) e_type);
  return it_ev_next(it_edge);
}

/****************************************************************************
 * get_unique_edge(v, e_type)
 * If v has a unique edge of type e_type, return it;
 * else return NULL
 ****************************************************************************/

EDGE *get_unique_edge(NODES *v, int e_type)
{
  IVEDGE it_edge_mem, *it_edge;
  EDGE *e;

  it_edge = it_ev_renew(&it_edge_mem, v, (E_TYPE) e_type);
  e = it_ev_next(it_edge);

  /* if there is a unique edge, return it */
  if (e != (EDGE *) 0 || it_ev_next(it_edge) == (EDGE *) 0)
    return e;

  /* else return NULL */
  return (EDGE *) 0;
}


/****************************************************************************
 * if v has indegree=1 and outdegree=1 in H,
 *    return the unique outgoing edge
 * else
 *    return null edge
 ****************************************************************************/

EDGE *find_unique_oedge(NODES *v)
{
  int i;
  EDGE *e_in = (EDGE *) 0;           /* set to unique incoming edge */ 
  EDGE *e_out = (EDGE *) 0;          /* set to unique outgoing edge */
  EDGE *e;
  IVEDGE it_edge_mem, *it_edge;

  /* explore incoming edges */
  it_edge = it_ev_renew(&it_edge_mem, v, E_IN_H);
  while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
    if (e_in != (EDGE *) 0) {
      /* already found one incoming edge, hence there's too many */
      return (EDGE *) 0;
    }

    /* first incoming edge encountered, store it */
    e_in = e;
  }
  if (e_in == (EDGE *) 0)
    return e_in;

  /* explore outgoing edges */
  it_edge = it_ev_renew(&it_edge_mem, v, E_OUT_H);
  while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
    if (e_out != (EDGE *) 0) {
      /* already found one outgoing edge, hence there's too many */
      return (EDGE *) 0;
    }

    /* first outgoing edge encountered, store it */
    e_out = e;
  }

  return e_out;
} 

/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * create_full_subgraph(G, edge_subg_flag, init_flag)
 *
 * create a subgraph H of G where V(H) = V(G) and
 * E(H) is determined by setting edge_subg_flag for all edges
 *
 * init_flag=1: does additional initialization at each vertex and edge
 ****************************************************************************/

void create_full_subgraph(IGRAPH *G, int edge_subg_flag, int init_flag)
{
  IVERT it_vert_mem, *it_vert;
  IVEDGE it_edge_mem, *it_edge;
  NODES *v;
  EDGE *e;

  if (init_flag) {
    G->visit_id = 0;
    G->visit_id_next = 1;
  }

  /* iterate over all vertices */
  it_vert = it_v_renew(&it_vert_mem, G, V_A_G);

  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    v->subg_flag = SUBG_IN;

    if (init_flag) {
      v->comp = v;
      v->visit_id = 0;
      v->visit = 0;

      /* fix the read readposition lists so that they are in sorted order */
      sort_nodepos(v);
    }

    it_edge = it_ev_renew(&it_edge_mem, v, E_OUT_G);
    while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
      e->subg_flag = edge_subg_flag;
    }
  }
}

void clear_visit_flag(IGRAPH *G, V_TYPE v_type)
{
  NODES *v;
  IVERT *it_vert = it_v_new(G, v_type);

  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    v->visit = 0;
  }
  it_v_destroy(it_vert);
}

void clear_visit_edge(IGRAPH *G, E_TYPE e_type)
{
  EDGE *e;
  IEDGE *it_edge = it_e_new(G, e_type);

  while ((e = it_e_next(it_edge)) != (EDGE *) 0) {
    e->visit = 0;
  }
  it_e_destroy(it_edge);
}

/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * comp_id(v)
 * compute the component id of a node
 ****************************************************************************/

COMP *comp_id(NODES *v)
{
  NODES *c,*tmp;

  /* follow chain of component ids until it stabilizes */
  c = v;
  while ((tmp = c->comp) != c) {
    c = tmp; 
  }

  /* collapse chain to speed it up for future passes */
  if (c != v->comp) comp_id_set(v, c);

  return (COMP *) c;
}

/****************************************************************************
 * comp_id_set(v,comp)
 * set the component id of v = comp
 * also modifies the component id of some other vertices in v's component
 ****************************************************************************/

void comp_id_set(NODES *v, NODES *comp)
{
  NODES *tmp;
  do {
    tmp = v->comp;
    v->comp = comp;
    v = tmp;
  } while (v != comp);
}

/****************************************************************************
 * form components of G
 * assumes create_full_subgraph(G, ?, 1)
 * was already run
 ****************************************************************************/
void form_components_subgraph(IGRAPH *G)
{
  IEDGE *it_edge;
  EDGE *e;

  it_edge = it_e_new(G, E_G);
  while ((e = it_e_next(it_edge)) != (EDGE *) 0) {
    comp_id_set(e->begin,comp_id(e->end));
  }
  it_e_destroy(it_edge);
}


/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * subgraph_ins_edge(e)
 * insert edge e of G into subgraph H
 * and join the components together
 ****************************************************************************/

void subgraph_ins_edge(EDGE *e)
{
  NODES *v1 = e->begin;
  NODES *v2 = e->end;
  NODES *c1 = comp_id(v1);
  NODES *c2 = comp_id(v2);

  /* edge should be included in H */
  e->subg_flag = SUBG_IN;

  /* Join components.
   * Use convention that smaller component ID becomes the new component ID.
   */
  if (c1 != c2) {
    if (c1 < c2)
      comp_id_set(v2,c1);
    else
      comp_id_set(v1,c2);
  }
}


/****************************************************************************
 * subgraph_del_edge(e)
 *
 * Delete edge e of G from subgraph H.
 * Assumes e was "postponed", so that no adjustments need to be made to
 * components.
 ****************************************************************************/

void subgraph_del_edge(EDGE *e)
{
  /* delete edge in H */
  e->subg_flag = SUBG_DEL;

  /* TODO:
   * do we need to delete the vertices in the event that their H-degree=0?
   */
}

/****************************************************************************
 * subgraph_del_node(v)
 *
 * Delete node v of G from subgraph H.
 ****************************************************************************/

void subgraph_del_node(NODES *v)
{
  /* delete node in H */
  v->subg_flag = SUBG_DEL;
}

/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * label_distances(G, v1, v2, L, e, e_type)
 *
 * Label distances in a neighborhood of radius <=L around vertex v1.
 * Returns a vertex at the end of some path.
 *
 * v1: start dfs here
 * v2: if v2 reached within distance L of v1, then immediately
 *     exit and return v2.  The distance labels may not be correct,
 *     but it is guaranteed v2 is within L of v1.
 *     Set v2=NULL to avoid this.
 * L: only label distances <= L of v1.
 * e: first step from v1 should not use edge e.
 * e_type: permitted edge type to follow at each step.
 *
 * return value:
 * if v2 is within L of v1
 *     return v2
 * else
 *     return a vertex w at the end of a longest path <= L emanating from v1.
 *     In the case of a tree, this is the farthest vertex from v1
 *     (or is tied for farthest).
 *     In other cases, there could be a shorter path to w, so it's not
 *     necessarily the vertex farthest from v1.
 *
 * check_path_length(G, v1, v2, L, e_type)
 *
 * Determine if node v2 is reachable from node v1 within specified length L.
 * Distances are marked on the nodes, so the path can be reconstructed
 * if desired.
 ****************************************************************************/

NODES *label_distances(IGRAPH *G,
		       NODES *v1, NODES *v2,
		       int L,
		       EDGE *e,
		       int e_type)
{
  int visit_id;

  /* The length field will only be valid when the time stamp
   * is correct.
   * Allocate two visit_id's:
   *   visit_id+0: mark nodes whose distances to v1 have been computed.
   *   visit_id+1: reserved to mark some of those nodes further.
   */
  visit_id = get_visit_id(G,2);

  /* Do a depth first search from v1 to mark nodes with
   * their distance from v1.
   * Stop if path length exceeds L, or if we reach v2.
   */

  return
    explore_neighborhood(G,
			 v1, v2,
			 e,
			 L,
			 0,           /* current distance */
			 e_type);
}


int check_path_length(IGRAPH *G,
		      NODES *v1, NODES *v2,
		      int L,
		      int e_type)
{
  return (label_distances(G,v1,v2,L, (EDGE *) 0, e_type) == v2);
}


/****************************************************************************
 * explore_neighborhood(G, v, v_stop, prev_edge, L, cur_dist, e_type)
 *
 * Subroutine of label_distances.
 * Recursively labels the vertices w/in distance L of v1 by their
 * minimum distance to v1.
 *
 * Currently we are at vertex v, a distance cur_dist from v1.
 * The edge we followed to get here is prev_edge, do not backtrack over it.
 *
 * Distances are marked on the nodes, so the path can be reconstructed
 * if desired.
 *
 * If v_stop is reached, return it.
 * Else, return the vertex at the end of the longest path <= L.
 ****************************************************************************/

NODES *explore_neighborhood(IGRAPH *G,
			    NODES *v,
			    NODES *v_stop,
			    EDGE *prev_edge,
			    int L,
			    int cur_dist,
			    int e_type)
{
  IVEDGE it_edge_mem, *it_edge;
  EDGE *e;
  NODES *v_next;
  NODES *v_farthest = v;            /* current farthest node */
  NODES *v_far2;                    /* candidate new farthest node */

  /* explored too far */
  if (cur_dist > L)
    return v;


  /* if (we already encountered this vertex &&
   *     it had a lower distance on the prior encounter)
   *   abort this branch of the search
   */


  if (v->visit_id == G->visit_id &&
      v->visit_value < cur_dist)
    return v;

  v->visit_id = G->visit_id;
  v->visit_value = cur_dist;

  /* return v_stop if we reached it */
  if (v == v_stop) {
    return v_stop;
  }


  /* explore all edges in undirected graph */
  it_edge = it_ev_renew(&it_edge_mem,
			v,
			(E_TYPE) e_type);
			//			E_U_H);

  while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
    /* don't backtrack */
    if (e == prev_edge)
      continue;

    /* next vertex depends on edge direction (in or out) */
    v_next = (e->begin == v) ? e->end : e->begin;

    /* depth first search from v_next */
    v_far2 = explore_neighborhood(G, v_next, v_stop, e, L,
				  cur_dist + e->length - 1,
				  e_type);

    /* exit if reached v_stop */
    if (v_far2 == v_stop)
      return v_far2;

    /* take farther of v_far2 or v_farthest */
    if (v_farthest->visit_value < v_far2->visit_value)
      v_farthest = v_far2;
  }

  return v_farthest;
}


/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * Node visits with timestamps:
 * When computing a value at each node in a sparse subset of the nodes,
 * we do not want to initialize the values in all nodes of G.
 * Instead, mark each node with a time stamp "visit_id".
 * The value at the node is only current if the "visit_id" is current.
 * Since the "visit_id" field may overflow, it will occasionally
 * be necessary to clear the field.
 *
 * ID = get_visit_id(G,r):
 *   Returns a new ID (and stores it in G->visit_id) to use to mark
 *   nodes on this traversal.
 *   The ID's   ID, ID+1, ..., ID+r-1
 *   are allocated for use on this traversal.
 *
 * reset_visit_id(G):
 *   Resets the visit_id field of all nodes.
 ****************************************************************************/

void reset_visit_id(IGRAPH *G)
{
  int i;
  NODES **nodes = G->nodes;
  int num_nodes = G->num_nodes;
  G->visit_id = 0;
  G->visit_id_next = 1;

  for (i=0; i<num_nodes; i++) {
    nodes[i]->visit_id = 0;
  }
}

int get_visit_id(IGRAPH *G, int r)
{
  G->visit_id = G->visit_id_next;
  if (((unsigned int) G->visit_id) + r >= (unsigned int) MAX_VISIT_ID) {
    fprintf(stderr, "Warning: unexpected event -- visit_id overflowed.  Resetting.\n");
    reset_visit_id(G);
    G->visit_id = 1;
  }

  G->visit_id_next = G->visit_id + r;

  return G->visit_id;
}


/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * mark_deleted_nodes(G)
 * locate all nodes in H of degree 0, and mark them deleted
 *
 * return number of nodes deleted
 ****************************************************************************/

int mark_deleted_nodes(IGRAPH *G)
{
  IVERT *it_vert;
  NODES *v;
  int num_del = 0;

  /* Iterate over all nodes of H */
  it_vert = it_v_new(G, V_A_H);

  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    if (gdegree_cmp(v,E_U_H,0) == 0) {
      v->subg_flag = SUBG_DEL;
      num_del ++;
    }
  }

  it_v_destroy(it_vert);

  return num_del;
}



/****************************************************************************/
/****************************************************************************/

/* calc_read_with_edge(e)
 *
 * Assume e has multiplicity 1.
 * Return the index of the unique read that contains this edge.
 */
int calc_read_with_edge(EDGE *e)
{
  NODES *v1 = e->begin;
  READPOSITION *pos1 = v1->readposition;
  int npos1 = v1->npos;

  NODES *v2 = e->end;
  READPOSITION *pos2 = v2->readposition;
  int npos2 = v2->npos;

  int r1,r2;
  int i;

  /* Examine (read,readposition) lists of both nodes till we find a match.
   */
  while (npos1 > 0   &&   npos2 > 0) {
    r1 = pos1->readindex;
    r2 = pos2->readindex;

    if (r1==r2) {
      return r1;
    } else if (r1 < r2) {
      pos1 = pos1->next;
      npos1--;
    } else {
      pos2 = pos2->next;
      npos2--;
    }
  }

  printf("Error: Edge %lx is supposed to have multiplicity 1, but is not contained in any read.\n",e);
  printf("Node %lx readindicies: ", v1);
  pos1 = v1->readposition;
  for (i=0; i<v1->npos; i++) {
    printf(" %d", pos1->readindex);
    pos1 = pos1->next;
  }
  printf("\n");
  printf("Node %lx readindicies: ", v2);
  pos2 = v2->readposition;
  for (i=0; i<v2->npos; i++) {
    printf(" %d", pos2->readindex);
    pos2 = pos2->next;
  }
  printf("\n");
  return -1;
}


/* Locate key in sorted array p[len].
 * If p[i]=key, set found=1 and return i.
 * Else, if p[i-1]<key<p[i], set found=0 and return i.
 * Exceptions: if key<p[0], set found=0 and return 0;
 *             if key>p[len-1], set found=len and return 0.
 */
int bsearch_int(int *p, int len, int key,

		/* return: */
		int *found)
{

  /* look for key in p[i] with i=lo,lo+1,...,hi-1 */
  int hi=len;
  int lo=0;
  int i;

  while (hi > lo) {
    i = (hi+lo)/2;
    if (p[i] == key) {
      *found = 1;
      return i;
    }
    if (p[i] > key) {
      hi = i;
    } else {
      lo = i+1;
    }
  }

  *found = 0;
  return lo;
}



void subgraph_del_read_with_edge(READTABLE *RT, EDGE *e)
{
  int readno = calc_read_with_edge(e);
  int chim_pos, found;
  int i;
  int *chim_old;

  /* ignore error */
  if (readno<0)
    return;

  /* if using the complementary read numbers, replace by original read # */
  if (readno >= RT->num_seq)
    readno -= RT->num_seq;

  /* Determine if readno is in the list of chimeric reads already or not.
   * Also determine its readposition.
   *
   * Usually # chimeric reads is much less than 100, so we will allocate
   * 100 readpositions initially, and use straight insertion sort to maintain
   * the list.
   */
  if (RT->chim == (int *) 0) {
    RT->num_chim = 0;
    RT->num_chim_alloc = 100;
    RT->chim = (int*) ckalloc(100 * sizeof(int *));
  }

  chim_pos = bsearch_int(RT->chim, RT->num_chim, readno, &found);

  if (found)
    return;

  /* insert it into the list */

  /* make sure there's sufficient memory first */
  if (RT->num_chim_alloc <= RT->num_chim) {
    /* allocate more space */
    chim_old = RT->chim;
    RT->num_chim_alloc += 500;
    RT->chim = (int*) ckrealloc((void *) RT->chim,
			 RT->num_chim_alloc * sizeof(int *),
			 RT->num_chim
			 );
  }

  /* list insertion */
  for (i=RT->num_chim; i>chim_pos; i--) {
    RT->chim[i] = RT->chim[i-1];
  }
  RT->chim[chim_pos] = readno;
      
  RT->num_chim++;
}

/****************************************************************************/
/****************************************************************************/


void free_igraph(IGRAPH *G)
{
  free_graph(G->nodes, G->num_nodes);
  free((void *) G->nodes);
}

/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * debug_print_vval(G)
 *   prints out list of edges and visit_value at each vertex
 ****************************************************************************/

void debug_print_vval(IGRAPH *G, READTABLE *RT)
{
  IVERT *it_vert;
  NODES *v;

  IVEDGE it_edge_mem, *it_edge;
  EDGE *e;

  char **src_seq;
  READPOSITION *pos;
  int i;

  fprintf(stderr,"DEBUG: entering debug_print_vval.  # nodes=%d, nodetable=%lx, \n", G->num_nodes, G->nodes);

  src_seq = RT->src_seq;

  it_vert = it_v_new(G, V_A_G);
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    fprintf(stderr,"DEBUG:  v=%lx (%d) [%lx] ->", v, v->visit_value, v->bal_node);
    it_edge = it_ev_renew(&it_edge_mem, v, E_OUT_H);
    while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
      fprintf(stderr," %lx <%lx> (%d);", e->end, e, e->end->visit_value);
    }
    fprintf(stderr," <-");
    it_edge = it_ev_renew(&it_edge_mem, v, E_IN_H);
    while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
      fprintf(stderr," %lx <%lx> (%d);", e->begin, e, e->begin->visit_value);
    }
    fprintf(stderr,"\n");

    /* (read,pos) at this node */
    fprintf(stderr,"DEBUG:  rp");
    pos = v->readposition;
    if (src_seq != (char **) 0) {
      for (i=0 ; i < v->npos ; i++) {
	fprintf(stderr," (%d,%d,%c)",
		pos->readindex, pos->position,
		na_name[src_seq[pos->readindex][pos->position]]
		);
	pos = pos->next;
      }
    } else {
      for (i=0 ; i < v->npos ; i++) {
	fprintf(stderr," (%d,%d)",
		pos->readindex, pos->position
		);
	pos = pos->next;
      }
    }
    fprintf(stderr,"\n");
  }
  it_v_destroy(it_vert);

  fprintf(stderr,"DEBUG: exiting debug_print_vval\n");
}

/****************************************************************************
 * debug_check_sym(G)
 *   verify symmetry of flags in vertices and edges
 ****************************************************************************/

void debug_check_sym(IGRAPH *G)
{
  IVERT *it_vert;
  NODES *v;

  int d1in, d1out, d2in, d2out;

  IEDGE *it_edge;
  EDGE *e;


  fprintf(stderr,"DEBUG: enter debug_check_sym\n");
  fprintf(stderr,"DEBUG:    check vertices\n");

  it_vert = it_v_new(G, V_A_H);
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    if (v->subg_flag != v->bal_node->subg_flag) {
      fprintf(stderr,"subg   %lx: %d   <>   %lx: %d\n",
	      v, v->subg_flag, v->bal_node, v->bal_node->subg_flag);
    }

    d1in = gdegree(v,E_IN_H);
    d1out = gdegree(v,E_OUT_H);
    d2in = gdegree(v->bal_node,E_IN_H);
    d2out = gdegree(v->bal_node,E_OUT_H);

    if (d1in != d2out || d2in != d1out) {
      fprintf(stderr,"deg    %lx: (%d,%d)   %lx: (%d,%d)\n",
	      v, d1in, d1out,
	      v->bal_node, d2in, d2out);
    }

  }
  it_v_destroy(it_vert);

  fprintf(stderr,"DEBUG:    check edges\n");
  it_edge = it_e_new(G,E_H);
  while ((e = it_e_next(it_edge)) != (EDGE *) 0) {
    if (e->subg_flag != e->bal_edge->subg_flag) {
      fprintf(stderr,"subg   %lx: %d   <>   %lx: %d\n",
	      e, e->subg_flag, e->bal_edge, e->bal_edge->subg_flag);
    }
  }
  it_e_destroy(it_edge);
  fprintf(stderr,"DEBUG: exit debug_check_sym\n");

}

/****************************************************************************
 * debug_print_degrees(G)
 *   print table of node degrees
 ****************************************************************************/

#define PD_MAXDEG 6
void debug_print_degrees(IGRAPH *G)
{
  int deg_table[PD_MAXDEG+1][PD_MAXDEG+1];

  IVERT *it_vert;
  NODES *v;

  int i,j;
  int d_in, d_out; 


  for (i=0; i<=PD_MAXDEG; i++)
    for (j=0; j<=PD_MAXDEG; j++)
      deg_table[i][j] = 0;


  it_vert = it_v_new(G, V_A_H);

  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    d_in = gdegree(v,E_IN_H);
    d_out = gdegree(v,E_OUT_H);

    if (d_in > PD_MAXDEG)
      d_in = PD_MAXDEG;

    if (d_out > PD_MAXDEG)
      d_out = PD_MAXDEG;

    deg_table[d_in][d_out]++;

#if 0
    if ((d_in==2 && d_out==0)
	||
	(d_in==0 && d_out==2)) {
      printf("v=%lx has degree (%d,%d)\n", v, d_in, d_out);
    }
#endif
  }


  printf("Degree table\n");
  printf("    \\ out\n");
  printf("  in \\   ");
  for (j=0; j<PD_MAXDEG; j++)
    printf("       %2d", j);

  printf("     >=%2d\n", PD_MAXDEG);

  printf("---------");
  for (j=0; j<=PD_MAXDEG; j++)
    printf("---------", j);
  printf("\n");

  for (i=0; i<=PD_MAXDEG; i++) {
    if (i<PD_MAXDEG)
      printf("  %2d     ", i);
    else
      printf(">=%2d     ", PD_MAXDEG);

    for (j=0; j<=PD_MAXDEG; j++)
      printf(" %8d", deg_table[i][j]);

    printf("\n");
  }
}



/****************************************************************************
 ****************************************************************************/

/****************************************************************************
 * v_ext = augment_graph_ext(G)
 *
 * Augment graph by adding a new vertex v_ext.
 * Add new edges from v_ext to every source and from every sink to v_ext.
 * Assumes node table was malloc'd with enough space for the new vertex
 ****************************************************************************/

/* augmentation edge v->w
 * dir=0:  v=v_ext, w=existing vertex
 * dir=1:  w=existing vertex, v=v_ext
 */

EDGE *new_aug_edge(NODES *v, NODES *w, int dir)
{
  EDGE *e;

  e = (EDGE *) ckalloc(sizeof(EDGE));

  /* fill in EDGE fields */
  e->begin = v;
  e->end = w;
  e->length = 1;  /* After accounting for node widths,
		   * passing through this edge does not affect the
		   * total path length.
		   */
  e->seq = (char *) 0;
  e->start_cover = 0;
  e->multip = 0;
  e->readinterval = (READINTERVAL *) 0;

  /* bal_edge: fill in later
   * visit: unused
   */

  e->subg_flag = SUBG_IN;


  /* Add e to the edge lists of v and w.
   * The malloc for v_ext is already large enough, the malloc for
   * the other vertex is not.
   */
  if (dir == 1) {
    v->nextedge = (EDGE **) ckrealloc(v->nextedge,
				      (v->num_nextedge + 1) * sizeof(EDGE *),
				      v->num_nextedge * sizeof(EDGE *));
  } else {
    w->lastedge = (EDGE **) ckrealloc(w->lastedge,
				      (w->num_lastedge + 1) * sizeof(EDGE *),
				      w->num_lastedge * sizeof(EDGE *));
  }

  v->nextedge[v->num_nextedge++] = e;
  w->lastedge[w->num_lastedge++] = e;

  return e;
}

NODES *augment_graph_ext(IGRAPH *G)
{
  IVERT *it_vert;
  NODES *v;

  EDGE *e, *ec;

  /* create v_ext */
  int num_source = it_v_count(G, V_SO_H);
  NODES *v_ext = (NODES *) ckalloc(sizeof(NODES));

  v_ext->readposition = (READPOSITION *) 0;
  v_ext->npos = 0;
  v_ext->nlinks = 0;
  v_ext->bal_node = v_ext;
  v_ext->num_path = 0;
  v_ext->path_index = (int *) 0;
  v_ext->path_pos = (int *) 0;
  v_ext->num_nextedge = 0;
  v_ext->num_lastedge = 0;
  v_ext->nextedge = (EDGE **) ckalloc(num_source * sizeof(EDGE *));
  v_ext->lastedge = (EDGE **) ckalloc(num_source * sizeof(EDGE *));
  v_ext->subg_flag = SUBG_IN;

  /* not used:
   *  comp, visit_id, visit_value, visit, ext_flag,
   *  comp_ext_flag
   */

  /* connect v_ext to all sources and sinks */

  it_vert = it_v_new(G, V_SO_H);
  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    /* create edges  v->v_ext  and  v_ext->v^c */
    e = new_aug_edge(v_ext, v, 0);
    ec = new_aug_edge(v->bal_node, v_ext, 1);
    e->bal_edge = ec;
    ec->bal_edge = e;
  }
  it_v_destroy(it_vert);


  /* Add v_ext to the node table.
   * It is mandatory that space was preallocated for storing v_ext in
   * the node table.
   */
  G->nodes[G->num_nodes++] = v_ext;
}


/****************************************************************************
 * unaugment_graph_ext(G)
 *
 * Remove augmentation vertex from graph
 ****************************************************************************/

void unaugment_graph_ext(IGRAPH *G)
{
  NODES *v_ext;
  EDGE *e;

  v_ext = G->nodes[G->num_nodes - 1];
  while (v_ext->num_nextedge) {
    /* e = v_ext -> w */
    e = v_ext->nextedge[v_ext->num_nextedge - 1];

    /* delete e from w's incoming list */
    e->end->num_lastedge --;

    /* delete e from v_ext's outgoing list */
    v_ext->num_nextedge --;

    /* delete e */
    free((void *) e);
  }

  while (v_ext->num_lastedge) {
    /* e = w -> v_ext */
    e = v_ext->lastedge[v_ext->num_lastedge - 1];

    /* delete e from w's outgoing list */
    e->begin->num_nextedge --;

    /* delete e from v_ext's incoming list */
    v_ext->num_lastedge --;

    /* delete e */
    free((void *) e);
  }

  free((void *) v_ext->nextedge);
  free((void *) v_ext->lastedge);
  free((void *) v_ext);

  G->num_nodes --;
}
