/***************************************************************************
 * Title:          zigzag.c
 * Author:         Glenn Tesler
 * Created:        Dec. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/* zigzag paths
 * Glenn Tesler
 * 12/26/02
 */

#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>
//#include <subgraph.h>

/* private */

//#define DEBUG_zz
#undef DEBUG_zz

void zzstats(NODES *v1, EDGE *e1,
	     NODES *v_stop,

	     /* return values */
	     NODES **v2,
	     EDGE **e2,
	     int *f,
	     int *r,
	     int *rmax,
	     int *rmin,
	     int *nturns);
EDGE *zz_nextedge(NODES *v, EDGE *e);
int is_zzvert(NODES *v,
	      EDGE **e1, EDGE **e2);
int process_zzpath(IGRAPH *G, READTABLE *RT,
		   NODES *v, EDGE *e_a, EDGE *e_b, int dir);
int process_zzpath_1(IGRAPH *G,
		     READTABLE *RT,
		      NODES *v1, EDGE *e1,
		      NODES *v2, EDGE *e2,
		      int f, int r,
		      int rmax, int rmin);

void zz_print(NODES *v1, EDGE *e1,
	      NODES *v2, EDGE *e2,
	      int f, int r,
	      int rmax, int rmin);


EDGE *zz_nextedge_spath(NODES *v, EDGE *e, E_TYPE e_type);
int zz_mergenode(IGRAPH *G, READTABLE *RT, NODES *v1, NODES *v2);
void zz_mergenode_1(IGRAPH *G, READTABLE *RT, NODES *v1, NODES *v2);
void zz_mergenode_readposition(IGRAPH *G, READTABLE *RT, NODES *v1, NODES *v2);
void zz_mergenode_edges(NODES *v1, NODES *v2, int out_flag);
EDGE *zz_isedgedup(NODES *v1, NODES *v2, int out_flag, EDGE *e2);
void zz_deledge(NODES *v2, int out_flag, EDGE *e2);


EDGE *zz_straighten(IGRAPH *G,
		    READTABLE *RT,
		    NODES *v1, EDGE *e1,
		    NODES *v2, EDGE *e2,
		    int f, int r,
		    int rmax, int rmin);

void zz_createedge(NODES *x, NODES *y);
void zz_createedge_add(NODES *v, EDGE *e, int out_flag);

NODES **zz_allocnodes(IGRAPH *G, int n);


void create_test_zzpath_short(IGRAPH *G, READTABLE *RT,
			      int zz_type, int L, int count0);

void zz_create_zzpath(IGRAPH *G, NODES *x, NODES *y, int sym, char *u, int L);

int zz_split_cycle(NODES *v0, EDGE *e1, int f0, int r0,
		   NODES **w, EDGE **e2);

#ifdef DEBUG_zz
static NODES *GLOBAL_v_final;
#endif



/****************************************************************************
 * zzstats(v1,e1,v_stop,
 *         *v2, *e2, *f, *r, *rmax, *rmin, *nturns)
 *
 * compute measurements of zigzag path emanating from v1 along e1
 * v_stop = vertex at end of path if known, NULL if not known
 *
 * returns:
 *    v2 = terminal vertex
 *    e2 = terminal edge
 *    f = # forward edges
 *    r = # reverse edges
 *    rmax = maximum rank
 *    rmin = minimum rank
 *    nturns = # direction changes
 ****************************************************************************/

void zzstats(NODES *v1, EDGE *e1, NODES *v_stop,

	     /* return values */
	     NODES **v2,
	     EDGE **e2,
	     int *f,
	     int *r,
	     int *rmax,
	     int *rmin,
	     int *nturns)
{
  int rk = 0;
  EDGE *e_next, *e;
  NODES *v;

  int last_dir = 0;    /* direction of last edge:
			*  1 = ->
			* -1 = <-
			* 0  = no previous edge
			*/

  *f = *r = *rmax = *rmin = *nturns = 0;

  v = v1;
  e = e1;

  while (1) {
    /* advance to next vertex */
    if (e->begin == v) {
      v = e->end;
      ++ *f;
      rk++;
      if (rk > *rmax)
	*rmax = rk;

      if (last_dir == -1)
	++*nturns;
      last_dir = 1;
	
    } else {
      v = e->begin;
      ++ *r;
      rk--;
      if (rk < *rmin)
	*rmin = rk;

      if (last_dir == 1)
	++*nturns;
      last_dir = -1;
    }

    /* terminate if path a cycle
     * or if we at the stop vertex */
    if (v == v1 || v == v_stop)
      break;

    /* determine next edge */
    e_next = zz_nextedge(v,e);

    /* terminate if at branching vertex */
    if (e_next == (EDGE *) 0)
      break;
    e = e_next;
  }

  /* in above loop, path terminates with vertex v and edge e */
  *v2 = v;
  *e2 = e;

  if (v_stop != (NODES *) 0
      && v != v_stop) {
    printf("Warning: zzstats unexpected termination. v_stop=%lx v=%lx e=%lx e_next=%lx\n", v_stop, v, e, e_next);
  }
}

/****************************************************************************
 * Find vertex w splitting path into two pieces
 * 
 *        ..p1.... 
 *      /         \
 *    v0           w
 *      \         /
 *        ..p2... 
 *
 *
 * s.t. |d1+d2| is maximized
 * where d1 = f1-r1,   f1 = # forwards edges (v0->w) on p1,  r1=# reverse
 *       d2 similar
 *
 * On input, given v0 and e1=1st edge along p1.
 *  f0 = # forwards edges on cycle in direction v0->p1->w->p2->v0
 *  r0 = # reverse edges
 *
 * On output:
 *  w = vertex maximizing d1+d2
 *  e2 = last clockwise edge into v0
 *
 * return value = |d1+d2| for this path
 *
 * TODO: compute all the other statistics about the paths p1,p2
 ****************************************************************************/

int zz_split_cycle(NODES *v0, EDGE *e1, int f0, int r0,
		   NODES **w, EDGE **e2)
{
  /* first put w at 0 readpositions along clockwise route */
  //  int f1 = 0;
  //  int r1 = 0;
  //  int d1 = f1-r1;

  //  int f2 = r0;
  //  int r2 = f0;
  //  int d2 = f0-r0;

  //  int d  = d1+d2;

  int d = f0 - r0;
  int n = 0;        /* # steps */

  int d_max = d;
  NODES *v_max  = v0;

  NODES *v  = v0;
  EDGE *e = e1;

  do {
    if (e == (EDGE *) 0) {
      printf("Warning: zz_split_cycle: premature cycle end\n");
      break;
    }
    *e2 = e;

    /* advance to next vertex */
    if (e->begin == v) {
      v = e->end;

      //      f1++;
      //      r2--;
      //      d1++;
      //      d2++;
      d += 2;
    } else {
      v = e->begin;

      //      r1++;
      //      f2--;
      //      d1--;
      //      d2--;
      d -= 2;
    }

    if (abs(d) > abs(d_max)) {
      d_max = d;
      v_max = v;
    }

    e = zz_nextedge(v,e);
    n++;
  } while (v != v0);

  *w = v_max;

  if (n != f0 + r0) {
    printf("Warning: zz_split_cycle: expected %d steps, got %d instead.  v0=%lx, v=%lx\n", f0+r0, n, v0, v);
  }

  return abs(d_max);
}

/****************************************************************************
 * Find vertices w1,w2 splitting path into two pieces
 * 
 *        ..p1.... 
 *      /         \
 *    w1           w2
 *      \         /
 *        ..p2... 
 *
 *
 * s.t. |d1+d2| is maximized
 * where d1 = f1-r1,   f1 = # forwards edges (w1->w2) on p1,  r1=# reverse
 *       d2 similar
 *
 * input:
 *    v0 = any source/sink on cycle
 *    e1 = edge on v0
 *    f0 = # forwards edges when traversing cycle in direction v0 through e1
 *    r0 = # reverse edges ...
 *
 * output: w1, w2, and edges e1a, e1b
 *
 * TODO: compute all the other statistics about the paths p1,p2
 ****************************************************************************/

void zz_split_cycle_comp(NODES *v0, EDGE *e1, int f0, int r0,
			 NODES **w1, NODES **w2, EDGE **e1a, EDGE **e1b)
{
  NODES *u1, *u2;

  int d_max = 0;
  int turned = 1;

  EDGE *e;
  EDGE *e2;

  int i;
  int d;

  *w1 = *w2 = (NODES *) 0;
  *e1a = *e1b = (EDGE *) 0;


  /* check if cycle has form x ... xc ... x */
  if (f0 == r0) {
    u1 = v0;
    e = e1;
    d = 0;
    for (i=f0; i>0; i--) {
      if (e->begin == u1) {
	u1 = e->end;
	d += 2;
      } else {
	u1 = e->begin;
	d -= 2;
      }
      e = zz_nextedge(u1,e);
    }

    if (e == e1->bal_edge) {
      /* does have form  x .. xc .. x */
#ifdef DEBUG_zz
      printf("Warning: zz_split_cycle_comp: cycle x..xc..x\n");
#endif

      /* x -> x++  does fwds edge        f1-- r2++
       * x' -> x'++ does backwds edge    r1++ f2--
       *                      net d1 -= 2   d2 -= 2    d -= 4
       *
       * reverse is d += 4
       */

      u1 = v0;
      e = e1;

      d_max = d;
      *w1 = u1;
      *e1a = e1;

      for (i=f0; i>0; i--) {
	if (e->begin == u1) {
	  u1 = e->end;
	  d -= 4;
	} else {
	  u1 = e->begin;
	  d += 4;
	}

	if (abs(d) > abs(d_max)) {
	  d_max = d;
	  *w1 = u1;
	  *e1a = e;
	}

	e = zz_nextedge(u1,e);
      }

      *w2 = (*w1)->bal_node;
      *e1b = zz_nextedge(*w1, *e1a);
      return;
    }
  }



  u1 = v0;
  e = e1;
  for (i=f0+r0; i>0; i--) {
    if (turned) {
      d = zz_split_cycle(u1, e, f0, r0,
			 &u2, &e2);

      if (abs(d) > abs(d_max)) {
	*w1 = u1;
	*w2 = u2;
	*e1a = e;
	*e1b = e2;
	d_max = d;
      }
    }

    /* advance vertex */
    if (e->begin == u1) {
      u1 = e->end;
      turned = 0;
    } else {
      u1 = e->begin;
      turned = 1;
    }

    e = zz_nextedge(u1,e);
    if (e->begin == u1) {
      turned = !turned;
    }
  }
}


/****************************************************************************
 * next edge in zigzag path
 * If v has degree != 2, path is done, return null.
 * Else, return the edge != e.
 ****************************************************************************/

EDGE *zz_nextedge(NODES *v, EDGE *e)
{
  EDGE *e_tmp  = (EDGE *) 0;
  EDGE *e_next = (EDGE *) 0;
  IVEDGE it_edge_mem, *it_edge;

  it_edge = it_ev_renew(&it_edge_mem, v, E_U_H);
  while ((e_tmp = it_ev_next(it_edge)) != (EDGE *) 0) {
    if (e_tmp != e) {
      /* is degree too large? */
      if (e_next != (EDGE *) 0)
	return (EDGE *) 0;

      e_next = e_tmp;
    }
  }

  /* got here, degree is 0, 1, or 2 */
  return e_next;
}

/****************************************************************************
 * next edge in straightened path
 * e_type = E_IN_H, E_OUT_H gives direction
 * e = edge in that direction to avoid
 ****************************************************************************/

EDGE *zz_nextedge_spath(NODES *v, EDGE *e, E_TYPE e_type)
{
  IVEDGE it_edge_mem, *it_edge;
  EDGE *e_tmp;

  it_edge = it_ev_renew(&it_edge_mem, v, e_type);
  while ((e_tmp = it_ev_next(it_edge)) != (EDGE *) 0) {
    if (e_tmp != e)
      return e_tmp;
  } 

  printf("ERROR: zz_nextedge_spath: no edge!\n");
  exit(-1);
  //  return (EDGE *) 0;
}


/****************************************************************************
 * locate and process all zigzag paths
 *
 * TODO:
 * processing zigzag paths involves modifying the nodes of G, not just of H.
 * This code loops over all the vertices of H, which currently is done
 * by looping over the vertices of G as they were on initialization.
 * If the vertex iteration method changes, this code may have to be changed
 * as well.
 ****************************************************************************/

void process_all_zzpaths(IGRAPH *G, READTABLE *RT)
{
  IVERT it_vert_mem, *it_vert;
  int dir;
  EDGE *e1, *e2;
  NODES *v;

  int num_zz = 0;

  /* locate a zigzag path and then process it
   */

  it_vert = it_v_renew(&it_vert_mem, G, V_A_H);

  while ((v = it_v_next(it_vert)) != (NODES *) 0) {
    if (dir = is_zzvert(v, &e1, &e2)) {
      num_zz += process_zzpath(G, RT, v, e1, e2, dir);
    }
  }
  printf("Straightened %d zigzag paths\n", num_zz); 
}

/****************************************************************************
 * check if v has indegree 2, outdegree 0,  or vice-versa.
 * If not:
 *   return 0
 *   e1,e2,dir are garbage
 * If so:
 *   e1,e2 are the edges
 *   return -1 (for indegree 2, outdegree 0)
 *          +1 (for indegree 0, outdegree 2)
 ****************************************************************************/

int is_zzvert(NODES *v,
	      EDGE **e1, EDGE **e2)
{
  int deg = 0;

  IVEDGE it_edge_mem, *it_edge;
  EDGE *e;

  it_edge = it_ev_renew(&it_edge_mem, v, E_IN_H);
  while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
    ++deg;
    if (deg == 1)
      *e1 = e;
    else if (deg == 2)
      *e2 = e;
    else
      return 0;
  }

  /* now indegree = 0, 1, or 2 */
  if (deg == 2) {
    if (gdegree_cmp(v, E_OUT_H, 0) == 0)
      return -1;  /* indegree 2, outdegree 0 */
    return 0; /* indegree 2, outdegree != 0 */
  } else if (deg == 1) {
    return 0; /* indegree 1 */
  }

  /* now indegree = 0 */
  it_edge = it_ev_renew(&it_edge_mem, v, E_OUT_H);
  while ((e = it_ev_next(it_edge)) != (EDGE *) 0) {
    ++deg;
    if (deg == 1)
      *e1 = e;
    else if (deg == 2)
      *e2 = e;
    else
      return 0;
  }

  /* now indegree = 0 and outdegree = 0, 1, or 2 */
  if (deg == 2)
    return +1;

  /* indegree = 0 and outdegree != 2 */
  return 0;
}


#define MIN2(a,b)  ((a) < (b) ? (a) : (b))
#define MAX2(a,b)  ((a) > (b) ? (a) : (b))

/****************************************************************************
 * Process the zigzag path containing vertex v
 * v is either a source or sink of degree 2
 * e_a, e_b are the edges incident with v
 * dir = +1: v is source
 *       -1: v is sink
 *
 * returns # paths straightened (0,1,2,3,4)
 ****************************************************************************/

int process_zzpath(IGRAPH *G, READTABLE *RT,
		   NODES *v, EDGE *e_a, EDGE *e_b, int dir)
{
  int rmax_a, rmax_b;
  int f_a, f_b;
  int r_a, r_b;
  int rmin_a, rmin_b;
  NODES *v2_a, *v2_b;
  EDGE *e2_a, *e2_b;

  int nturns_a, nturns_b;

  int v_flag, w_flag;

  EDGE *e_c;

  NODES *w = (NODES *) 0;  /* Split cycle into two paths.
			    * If non-null, this is the vertex
			    * on which to split.
			    */

  int dual_pair;    /* is cycle = path + path' ? */

  int num_zz = 0;

  /* measure the two branches, a and b, of a zigzag path */
  zzstats(v, e_a, (NODES *) 0,
	  &v2_a, &e2_a,
	  &f_a, &r_a, &rmax_a, &rmin_a, &nturns_a);

  if (v == v2_a) {
    /* zigzag path forms an entire component, and is a cycle */

    if (nturns_a < 2) {
#ifdef DEBUG_zz
      printf("Warning: zigzag path forms cycle component, %d parallel paths, ignoring\n", nturns_a+1);
#endif
      return 0;
    }

#ifdef DEBUG_zz
    printf("Warning: zigzag path forms cycle component\n");
    zz_print(v, e_a, v2_a, e2_a, f_a, r_a, rmax_a, rmin_a);
    printf("\n");
#endif

    /* Find appropriate start, end verts in component.
     * We break it up into several overlapping computations.
     *
     * TODO: Could speed it up by combining them together,
     * but these paths are rare.
     */

    /* split the cycle */
    zz_split_cycle_comp(v, e_a, f_a, r_a,
			&v, &w, &e_a, &e_c);
  }

  /* if not a cycle component, get the other half of the path
   */
  if (w == (NODES *) 0) {
    zzstats(v, e_b, (NODES *) 0,
	    &v2_b, &e2_b,
	    &f_b, &r_b, &rmax_b, &rmin_b, &nturns_b);

    if (v2_a == v2_b) {
      /* zigzag path forms a cycle x ... v ... x */

      if (nturns_a == 0 && nturns_b == 0) {
#ifdef DEBUG_zz
	printf("Warning: zigzag path forms cycle, two parallel paths, ignoring\n");
#endif
	return 0;
      }

      /* Find appropriate end vert in cycle.
       * We break it up into several overlapping computations.
       *
       * TODO: Could speed it up by combining them together,
       * but these paths are rare.
       */

      /* move to start of cycle */
      v = v2_a;
      e_a = e2_a;

      /* get midpoint of cycle, and final edge of cycle */
      zz_split_cycle(v, e_a, f_a+r_b, r_a+f_b,
		     &w, &e_c);

#ifdef DEBUG_zz
      printf("Warning: zigzag path forms cycle, split it.  v=%lx  w=%lx\n", v, w);
#endif
    }
  }


  /* handle cycles and cycle components */
  if (w != (NODES *) 0) {
    /* process two halves of cycle separately */
    zzstats(v, e_a, w,
	    &v2_a, &e2_a,
	    &f_a, &r_a, &rmax_a, &rmin_a, &nturns_a);

    if (v != w) {
      zzstats(v, e_c, w,
	      &v2_b, &e2_b,
	      &f_b, &r_b, &rmax_b, &rmin_b, &nturns_b);

      v_flag = v == e_c->begin;
      w_flag = w == e2_b->begin;
    }

    /* anomaly:
     * cycle could be two dual halves
     *    1->2->...->n->1'->...->n'->1
     * with v=i, w=i', for some i
     */
    dual_pair = e_c == e2_a->bal_edge;


#ifdef DEBUG_zz
    printf("1st half  v=%lx  e_a=%lx  w=%lx\n", v, e_a, w);
#endif
    if (nturns_a > 0) {
      num_zz +=
	process_zzpath_1(G, RT,
			 v, e_a,
			 v2_a, e2_a, f_a, r_a, rmax_a, rmin_a);
    }

    if (v != w && !dual_pair) {
#ifdef DEBUG_zz
      printf("2nd half  v=%lx  e_c=%lx  w=%lx\n", v, e_c, w);
#endif
      if (nturns_b > 0) {
	/* careful: node merging in 1st half may have moved v or w */
	num_zz += 
	  process_zzpath_1(G, RT, 

			   /* v, possibly merged w/another vertex */
			   v_flag ? e_c->begin : e_c->end,
			   e_c,

			   /* w, possibly merged w/another vertex */
			   w_flag ? e2_b->begin : e2_b->end,
			   
			   e2_b, f_b, r_b, rmax_b, rmin_b);
      }
    }

    return num_zz;
  }


  /* Zigzag path is between two distinct branching vertices.
   * Moving from v to v2_a causes the edges on the a-branch to
   * be inverted.
   */

  return
    process_zzpath_1(G, RT,
		   v2_a, e2_a,
		   v2_b, e2_b,
		   r_a + f_b,   /* forward edges from a */
		   f_a + r_b,   /* forward edges from b */
		   MAX2(rmax_a,rmax_b) - (f_a-r_a),
		   MIN2(rmin_a,rmin_b) - (f_a-r_a));

}

/*
 * zigzag path may have been split into pieces.
 * process each piece separately.
 */

int process_zzpath_1(IGRAPH *G,
		     READTABLE *RT,
		      NODES *v1, EDGE *e1,
		      NODES *v2, EDGE *e2,
		      int f, int r,
		      int rmax, int rmin)
{

  EDGE *E0;

  int num_zz = 0;

#ifdef DEBUG_zz
  int rmax_c, rmin_c;
  int f_c, r_c;
  NODES *v2_c;
  EDGE *e2_c;
  int nturns_c;
  int d_expected;

  IVEDGE it_e_mem, *it_edge;
  EDGE *e_tmp;

  if (e1->bal_edge == e2)
    printf("Warning: zigzag path is self-complementary\n");
  else
    printf("Warning: zigzag path good\n");

  zz_print(v1, e1,
	   v2, e2,
	   f, r, rmax, rmin);
  printf("\n");

  printf("Straightening above zigzag path\n");
#endif

  if (e1->bal_edge == e2) {
    num_zz = 1;
  } else {
    num_zz = 2;
  }

  E0 = zz_straighten(G, RT, v1, e1, v2, e2, f, r, rmax, rmin);

#ifdef DEBUG_zz
  /* test its measurements */
  zzstats(E0->begin, E0,
	  GLOBAL_v_final,
	  /* (f<r) ? v1 : v2,*/
	  &v2_c, &e2_c, &f_c, &r_c, &rmax_c, &rmin_c, &nturns_c);
  d_expected = abs(f-r);

  if (f_c - r_c != d_expected) {
    printf("Warning: zigzag path straightening gone awry:\n");
    printf("    v1=%lx e1=%lx v2=%lx e2=%lx f=%d r=%d rmax=%d rmin=%d\n",
	   v1,e1,v2,e2,f,r,rmax,rmin);
    printf("    v2_c=%lx e2_c=%lx f_c=%d r_c=%d rmax_c=%d rmin_c=%d\n",
	   v2_c,e2_c,f_c,r_c,rmax_c,rmin_c);

    zz_print(E0->begin, E0, v2_c, e2_c, f_c, r_c, rmax_c, rmin_c);

    printf("Edges on %lx:\n", v2_c);
    it_edge = it_ev_renew(&it_e_mem, v2_c, E_U_H);
    while ((e_tmp = it_ev_next(it_edge)) != (EDGE *) 0) {
      printf("   e=%lx:   %lx -> %lx\n", e_tmp, e_tmp->begin, e_tmp->end);
    }

  }
#endif

  return num_zz;
}



/****************************************************************************
 * print a zigzag path
 * starts from node v1 on edge e1
 * continues on degree 2 vertices
 * ends at node v2 on entry through e2
 *
 * path should have f forward edges, r reverse edges,
 * and achieve max rank = rmax,  min rank = rmin
 ****************************************************************************/

void zz_print(NODES *v1, EDGE *e1,
	      NODES *v2, EDGE *e2,
	      int f, int r,
	      int rmax, int rmin)
{
  EDGE *e_next, *e;
  NODES *v;

  /* x_ is for verification */
  int f_ = 0;
  int r_ = 0;
  int rmax_ = 0;
  int rmin_ = 0;
  int rk_ = 0;

  /* want to process it in "forwards" direction */
  if (f < r) {
    zz_print(v2, e2, v1, e1, r, f,
	     r-f+rmax, r-f+rmin);
    return;
  }

#if 0
  if ((rmax == f && rmin == 0) ||
      (rmin == -r && rmax == f-r)) {
    printf("Warning: only one dir change\n");
  }
#endif


  /* print path and verify statistics */
  v = v1; e = e1;
  while (1) {
    printf("%lx", v);

    printf(" [%lx]", e); /* DEBUG */

    /* advance to next vertex */
    if (e->begin == v) {
      printf(" -> ");
      v = e->end;
      ++ f_;
      rk_ ++;
      if (rk_ > rmax_ )
	rmax_ = rk_ ;
    } else {
      printf(" <- ");
      v = e->begin;
      ++ r_;
      rk_--;
      if (rk_ < rmin_)
	rmin_ = rk_;
    }


    /* at end of path? */
    if (v == v2) {
      printf("%lx\n", v);
      break;
    }

    /* determine next edge */
    e_next = zz_nextedge(v,e);

    if (e_next == (EDGE *) 0) {
      printf("\nWarning: premature end of zigzag path!\n");
      break;
    }

    e = e_next;

    if (f_ + r_ > f + r) {
      printf("\nWarning: zigzag path overflow, aborting\n");
      break;
    }
  }


  /* verify statistics */
  if (v2 != v ||
      e2 != e ||
      f != f_ ||
      r != r_ ||
      rmax != rmax_ ||
      rmin != rmin_ ) {
    printf("Warning: zigzag path did not have expected parameters!\n");
    printf("Got  v2=%lx, e2=%lx, f=%d, r=%d, rmax=%d, rmin=%d\n",
	   v,e,f_,r_,rmax_,rmin_);
    printf("Expected  v2=%lx, e2=%lx, f=%d, r=%d, rmax=%d, rmin=%d\n",
	   v2,e2,f,r,rmax,rmin);
  }
}



/****************************************************************************
 * straighten zigzag path
 * also straightens dual path
 * handles complication that the path equals its dual path
 *
 * returns the forwards edge that starts the new path at one or other end
 ****************************************************************************/

EDGE *zz_straighten(IGRAPH *G,
		    READTABLE *RT,
		    NODES *v1, EDGE *e1,
		    NODES *v2, EDGE *e2,
		    int f, int r,
		    int rmax, int rmin)
{
  int rk, rk_low, rk_high, rk_max;
  EDGE *e, *E1, *E2;

  int i;

  NODES *v;
  NODES *v_next;
  NODES *v_merge;
  NODES *v_high = v1;   /* The vertex of rank rk_high;
			 * init only needed if rk_high = 0 */
  NODES *v_lo = v1;

  EDGE *e_next;
  EDGE *E0 = (EDGE *) 0;

  /* self-dual paths: */
  int sd_thresh = -1;   /* midpoint of self-dual path */
  int merge_err = 0;    /* self-dual path may encounter anomaly in middle */
  EDGE *E_alt = (EDGE *) 0;  /* parallel to edge to choose */
  int e_alt_type;

  /* traverse the zigzag path, straightening as we go.
   * This splits the zz path into 3 pieces:
   *
   *                             rk
   *      E0                  E1   E2             E3
   *  v1  -> . -> . -> . -> . -> x -> . -> . -> . -> y
   *                           e/?\e
   *                   ...v_next    v_next...
   *  v1 = start node
   *  x  = current node
   *  e  = edge on zigzag path leaving x
   *  v1 -> ... y  = straightening of all nodes up to x
   *  rk = rank
   *  v_next = next vertex on zigzag path
   *           (unique for any given x, but 2 directions are possible)
   *  v_merge = vertex parallel to v_next but on the top path v1->...->y
   *
   *  piece v1 -> ... -> x: straightening of nodes up to x with rank <= rank(x)
   *  piece x -> ... -> y: straightening of nodes up to x with rank >= rank(x)
   *  piece x  :  what's left to process (only one down edge, either / or \ )
   *       /?\
   *      .....
   *
   * The path and the pieces are linked lists
   *
   *
   * Old way:
   *      when e is in range, merge v_next into v_merge and e into E1/E2.
   * New way:
   *      when e is in range, merge v_merge into v_next and E1/E2 into e.
   */

  /* want to process it in "forwards" direction */
  if (f < r) {
    return zz_straighten(G, RT,
			 v2, e2, v1, e1, r, f,
			 r-f+rmax, r-f+rmin);
  }

  if (f == r) {
    /*
     * TODO: erase whole path, merge first and last node
     * That is the net effect of the code below in this case, so
     * just do it faster.
     */
  }

  if (e1 == e2->bal_edge) {
    //    printf("Warning: straightening self-dual path\n");
    sd_thresh = (f+r)/2;
  }


  v = v1; e = e1;

  /* nodes of rank outside this interval will be deleted */
  rk_low = 0;
  rk_high = f-r; 

  rk = 0;         /* current rank */
  rk_max = -1;    /* highest rank straightened so far */


  /* traverse edges, merging nodes/edges as appropriate */
  /* v = current vertex on zz-path
   *     and on straightened path if rank(v) is in range
   */

  v = v1; e = e1;
  for (i=f+r ; i>0 ; i--) {

    if (i == sd_thresh) {
      /* middle of self-dual path, adjust parameters */

      /* edges later on path were combined simultaneously
       * as the ones previously on the path were combined,
       * so there are this many left:
       */
      i = MIN2(rk_max - rk, rk-(rk_high-rk_max));
      if (i<=0) break;

      if (e == (EDGE *) 0 && E_alt != (EDGE *) 0) {
	e = zz_nextedge_spath(v_next, E_alt, (E_TYPE) e_alt_type);
      }
    }

    /* at end of path? */
    if (e == (EDGE *) 0) {
      printf("Warning: zzpath terminating unexpectedly, i=%d\n", i); 
      break;
    }


    /* advance to next vertex */
    if (e->begin == v) {
      /* forwards edge */
      v_next = e->end;
      v_merge = v_next;

      /* get next edge before cur vertex is destroyed */
      e_next = zz_nextedge(v_next,e);

      /* rank of v_next */
      rk++;


      if (rk > rk_high || rk <= 0) {

	/* out of range, must delete the edge we just traversed */
	subgraph_del_edge(e->bal_edge);
	subgraph_del_edge(e);

	if (rk == 0) {

	  /*
	   * at rank 0, must merge v1 into v_next
	   *
	   *          rk
	   *              E0
	   *          v1  -> . -> . -> . -> . -> . -> . -> y
	   *        v_next
	   *      /
	   *     v
	   */
	  if (v_next != v_lo) {
	    v_merge = v_lo;
	    merge_err = zz_mergenode(G, RT, v_merge, v_next);
	    v_lo = v_next;
	  }
	}
      } else {

	/* if E2 exists:
	 *         merge edge with E2 and merge v_next into E2->end
	 * else:
	 *         extend straightened path by setting E2 = e
	 */
	if (rk <= rk_max) {

	  /* E2 already exists,
	   * merge E2 into e and E2->end into v_next
	   *
	   *
	   *                                  rk
	   *      E0                  E1   E2             E3
	   *  v1  -> . -> . -> . -> . -> v -> . -> . -> . -> y
	   *                               \
	   *                                \
	   *                                  v_next...
	   */

	  E2 = (rk == 1) ? E0 : zz_nextedge_spath(v,e,E_OUT_H);
	  v_merge = E2->end;

	  /* midpoint of self-dual path, need to adjust next edge */
	  if (e_next == (EDGE *) 0 &&
	      2*i == f+r+1) {
	    E_alt = zz_nextedge(v_merge, E2);
	    e_alt_type = E_OUT_H;
	  }

	  zz_mergenode(G, RT, v_merge,v_next); /* merge E2->end into v_next */
	} else {

	  /* E2 doesn't yet exist, extend straightened path
	   *      E0                  E1  
	   *  v1  -> . -> . -> . -> . -> v
	   *                              \ (becomes E2)
	   *                               v_next...
	   */

	  rk_max = rk;
	}

	if (rk == rk_high) {
	  /* replace highest vertex so far on straightened path */

	  v_high = v_next;
	}

	if (rk == 1) {
	  /* replace 1st edge on straightened path */
	  E0 = e;
	}
      }

    } else {
      /* reverse edge */
      v_next = e->begin;
      v_merge = v_next;

      /* get next edge before cur vertex is destroyed */
      e_next = zz_nextedge(v_next,e);

      /* rank of v_next */
      rk--;

      if (rk >= rk_high || rk < 0) {
	/* out of range, must delete the edge we just traversed */
	subgraph_del_edge(e->bal_edge);
	subgraph_del_edge(e);

	if (rk == rk_high) {
	  /* merge node only
	   *                                               rk
	   *              E0
	   *          v1  -> . -> . -> . -> . -> . -> . -> y
	   *                                                     v
	   *                                                   /
	   *                                               v_next
	   */

	  v_merge = v_high;
	  merge_err = zz_mergenode(G, RT, v_merge, v_next);

	  /* replace highest vertex on path */
	  v_high = v_next;
	}

      } else {
	/* E1 exists:
	 *         merge E1 into e and merge E1->begin into v_next
	 *                       rk
	 *      E0                  E1   E2             E3
	 *  v1  -> . -> . -> . -> . -> v -> . -> . -> . -> y
	 *                           /
	 *                          /
	 *                        v_next...
	 */
	E1 = zz_nextedge_spath(v,e,E_IN_H);
	v_merge = E1->begin;

	/* midpoint of self-dual path, need to adjust next edge */
	if (e_next == (EDGE *) 0 &&
	    2*i == f+r+1) {
	  E_alt = zz_nextedge(v_merge, E1);
	  e_alt_type = E_IN_H;
	}

	/* merge E1->begin into v_next */
	merge_err =
	  zz_mergenode(G, RT, v_merge, v_next);
      }

      if (rk == 0) {
	/* replaced E0 */
	E0 = e;
	v_lo = v_next;
      }

    }

    if (merge_err) {
      /* corrections */
      printf("zz_straighten: anomaly encountered\n");
    }

#if 0
    /* midpoint of self-dual path */
    if (e_next == (EDGE *) 0 &&
	2*i == f+r+1) {
      /* adjust choice of next edge */
    }
#endif

    //    v = v_merge;
    v = v_next;
    e = e_next;

  }

#ifdef DEBUG_zz
  GLOBAL_v_final = (sd_thresh < 0) ? v : E0->begin->bal_node;
#endif

  return E0;
}




/* merge node v1 into v2
 * and v1c into v2c
 *
 * edges, reads, etc. using v1 (resp. v1c) are changed to use v2 (resp. v2c)
 *
 *     cases expressed as partitions of which of v1,v1c,v2,v2c  are equal:
 *     {v1,v2},{v1c,v2c}: do nothing
 *     {v1,v2,v1c,v2c}:   do nothing
 *
 *     {v1,v1c},{v2},{v2c}:
 *          shouldn't happen on a zigzag path
 *          merge v2 into v1 and v2c into v1c
 *          return 1
 *
 *     {v1},{v1c},{v2},{v2c}:
 *     {v1},{v1c},{v2,v2c}:
 *          merge v1 into v2
 *          and v1c into v2c
 *
 *     {v1,v1c},{v2,v2c}:
 *     {v1,v2c},{v1c,v2}:
 *          merge v1 into v2
 *
 */
int zz_mergenode(IGRAPH *G, READTABLE *RT, NODES *v1, NODES *v2)
{
  NODES *v1c, *v2c;
  int rval = 0;

  if (v1 == v2)
    return rval;

  v1c = v1->bal_node;
  v2c = v2->bal_node;

  if ((v1 == v1c && v2 == v2c) ||
      (v1 == v2c && v1c == v2)) {
    zz_mergenode_1(G, RT, v1, v2);
  } else if (v1 == v1c) {
    zz_mergenode_1(G, RT, v2, v1);
    zz_mergenode_1(G, RT, v2c, v1c);
    printf("zz_mergenode warning, v1=v1c=%lx, v2=%lx, v2c=%lx\n", v1, v2, v2c);
    rval = 1;
  } else {
    zz_mergenode_1(G, RT, v1, v2);
    zz_mergenode_1(G, RT, v1c, v2c);
  }

  return rval;
}


/* merge node v1 into node v2
 * optimize for zigzag paths, where there is some vertex v3 s.t.
 * (v2,v3) & (v1,v3) are directed edges
 * OR (v3,v2) & (v3,v1) are directed edges
 *
 *
 * Node v2 must remain at the same memory address
 * Its contents can be changed, however.
 *
 * Node v1 must remain at the same memory address, but use flags to
 * say it's deleted.
 * However, the set union structure for components may have pointers to
 * this node.  For now, these are not actually used any more, but in the
 * future, if they are, we must update them appropriately.
 */
void zz_mergenode_1(IGRAPH *G, READTABLE *RT, NODES *v1, NODES *v2)
{

#ifdef DEBUG_zz
  printf("zz_mergenode_1 v1=%lx v2=%lx\n", v1, v2);
#endif

  /***********************************************************************
   * update fields
   *    readposition, npos
   ***********************************************************************/

  zz_mergenode_readposition(G, RT, v1, v2);


  /***********************************************************************
   * update field
   *     nlinks
   ***********************************************************************/

  v2->nlinks += v1->nlinks + 1;
  v1->nlinks = 0;

  /***********************************************************************
   * update field
   *     bal_node
   ***********************************************************************/

  /* field is used, but nothing to do */

  /***********************************************************************
   * update fields
   *     num_path, path_index, path_pos
   ***********************************************************************/

  /* fields are not used yet */

  /***********************************************************************
   * update fields
   *     num_nextedge, nextedge
   *     num_lastedge, lastedge
   ***********************************************************************/

  zz_mergenode_edges(v1, v2, 1);    /* outgoing edges */
  zz_mergenode_edges(v1, v2, 0);    /* incoming edges */

  /***********************************************************************
   * update fields
   *     comp
   ***********************************************************************/

  /* We assume they are the same component, and do nothing.
   * If they are in different components, the caller must deal with it.
   */

  /***********************************************************************
   * update fields
   *     visit_id, visit_value
   *     visit
   *     ext_flag
   *     comp_ext_flag
   ***********************************************************************/

  /* values are garbage, caller must recompute them if desired */

  /***********************************************************************
   * update fields
   *     subg_flag
   ***********************************************************************/

  v1->subg_flag = SUBG_DEL;

}


void zz_mergenode_readposition(IGRAPH *G, READTABLE *RT, NODES *v1, NODES *v2)
{
  READLIST **readlist = RT->readlist;
  int i, new_npos;
  READPOSITION *pos1, *pos2;
  READPOSITION **new_poslist;   /* pointer to the readposition pointer */

  new_npos = v1->npos + v2->npos;
  pos1 = v1->readposition;
  pos2 = v2->readposition;

  /* the readposition lists are linked lists in increasing order of readindex.
   * merge the 1st list into the 2nd.
   */

  new_poslist = &v2->readposition;
  for (i=0; i<new_npos; i++) {
    if (pos1 != (READPOSITION *) 0
	&& pos2 != (READPOSITION *) 0
	&& pos1->readindex == pos2->readindex) {
     printf("Warning: merging nodes %lx and %lx creates a fat supernode, with readpositions %d and %d from read %d\n", v1, v2, pos1->position, pos2->position, pos1->readindex);
    }

    if (pos1 == (READPOSITION *) 0
	|| (pos2 != (READPOSITION *) 0 &&
	    pos2->readindex < pos1->readindex)) {
      /* advance one read in pos2 */
      *new_poslist = pos2;
      new_poslist = &pos2->next;
      pos2 = pos2->next;

    } else {
      /* merge one read of pos1 into list, and advance in pos1 */

      /* read has a pointer from this readposition to node v1;
       * change it to point to v2 */
      readlist[pos1->readindex][pos1->position].node = v2;

      /* insert node into list */
      *new_poslist = pos1;
      new_poslist = &pos1->next;
      pos1 = pos1->next;
    }
  }

  /* new linked lists for v1 and v2: terminate and adjust lengths */

  /* v2: terminate and adjust length */
  *new_poslist = (READPOSITION *) 0;
  v2->npos = new_npos;

  /* v1: null list */
  v1->npos = 0;
  v1->readposition = (READPOSITION *) 0;
}


/*
 * This is designed for zigzag paths,
 * in which usually both nodes have (in/out) degree = 1 or 2
 *
 * out_flag = 1: merge together the outgoing edges
 *            0: merge together the incoming edges 
 *
 * v1 is merged into v2, so
 *     if there is an edge (v2,v3) & (v1,v3), keep (v2,v3), discard (v1,v3);
 *     if there is (v2,v3) but no (v1,v3), keep (v2,v3);
 *     if there is (v1,v3) but no (v2,v3), replace by (v2,v3)
 *
 * TODO: should we merge all edges on v2 and v1 in G or in H?
 */

void zz_mergenode_edges(NODES *v1, NODES *v2, int out_flag)
{
  EDGE **edges1;
  EDGE **edges2;
  EDGE **new_edgelist, **new_edgelist_next;

  int e_type;
  int deg1, deg2, ndups, nedges;

  int deg1G, deg2G;
  int i;

  IVEDGE it_e1_mem, *it_e1;
  EDGE *e1, *e2;

  if (out_flag) {
    e_type = E_OUT_H;
    deg1G = v1->num_nextedge;
    deg2G = v2->num_nextedge;
    edges1 = v1->nextedge;
    edges2 = v2->nextedge;
  } else {
    e_type = E_IN_H;
    deg1G = v1->num_lastedge;
    deg2G = v2->num_lastedge;
    edges1 = v1->lastedge;
    edges2 = v2->lastedge;
  }

  /* need to merge together the edge lists,
   * and remove edges that become duplicates
   *
   * For zz paths, usually there are 0 duplicates in one direction
   * and 1 duplicate in the other direction.
   *
   * We will allocate memory for ptrs to deg1 + deg2 - ndups edges.
   *
   * This way of counting duplicates is horribly inefficient in general.
   * For most vertices on zz paths, degrees are one of 0, 1, 2.
   */

  deg1 = gdegree(v1, (E_TYPE) e_type);
  deg2 = gdegree(v2, (E_TYPE) e_type);

  ndups = 0;
  it_e1 = it_ev_renew(&it_e1_mem, v1, (E_TYPE) e_type);

  while ((e1 = it_ev_next(it_e1)) != (EDGE *) 0) {
    if (zz_isedgedup(v1, v2, out_flag, e1))
      ndups++;
  }

  nedges = deg1 + deg2 - ndups;

  new_edgelist = (EDGE **) ckalloc(nedges * sizeof(EDGE *));

  /* copy all pointers from edgelist of v2
   * copy some pointers from edgelist of v1
   */


  /* edgelist of v2 */
  new_edgelist_next = new_edgelist;

  for (i=0 ; i<deg2G; i++) {
    /* cases:
#if 0
     * 1. edge i is null (already deleted in G)
#endif
     * 2. edge i is deleted in H
     * 3. edge i is in H and would collapse into an equivalent edge w/v1
     * 4. edge i is in H and is distinct from edges with v1
     */

    e2 = edges2[i];

#if 0
    /* 1. edge i is null (already deleted in G) */
    if (e2 == (EDGE *) 0)
      continue;
#endif

    /* 2. edge i is deleted (or postponed, treat as deleted) in H */
    if (e2->subg_flag != SUBG_IN) {
      /* change other vertex's pointer to this edge to null */
      zz_deledge(v2, out_flag, e2);

      /* free the edge */
      free((void *) e2);

      continue;
    }



    /* 3. edge i is in H and would collapse into an equivalent edge w/v1.
     *    Then keep this edge, and later delete the one w/v1.
     */
    /* 4. edge i is in H */

    /* append it onto new edge list */
    *new_edgelist_next++ = e2;

    /* no need to update other vertex of e2 */
    
  }


  /* edgelist of v1 */
  /* copy some pointers of edgelist of v1,
   * modifying the ones that we copy,
   * destroying the ones that we don't copy
   */
  for (i=0 ; i<deg1G; i++) {
    /* cases:
#if 0
     * 1. edge i is null (already deleted in G)
#endif
     * 2. edge i is deleted in H
     * 3. edge i is in H and would collapse into an equivalent edge w/v2
     * 4. edge i is in H and is distinct from edges with v2
     */

    e1 = edges1[i];

#if 0
    /* 1. edge i is null (already deleted in G) */
    if (e1 == (EDGE *) 0)
      continue;
#endif

    /* 3. edge i is in H and would collapse into an equivalent edge w/v2 */
    if ((e2 = zz_isedgedup(v1, v2, out_flag, e1)) != (EDGE *) 0) {
      e1->subg_flag = SUBG_DEL;

      /* reduced to case 2; continue with case 2 */
    }

    /* 2. edge i is deleted (or postponed, treat as deleted) in H */
    if (e1->subg_flag != SUBG_IN) {
      /* change other vertex's pointer to this edge to null */
      zz_deledge(v1, out_flag, e1);

      /* free the edge */
      free((void *) e1);

      continue;
    }



    /* 4. edge i is in H and is distinct from edges with v2 */
    /* change v1 into v2 */

    if (out_flag) {
      e1->begin = v2;
    } else {
      e1->end = v2;
    }

    /* append it onto new edge list */
    *new_edgelist_next++ = e1;

    /* no need to update other vertex of e1 */
    
  }

  /* replace edgelists of v2 and v1 */
  if (out_flag) {
    free(v1->nextedge);
    free(v2->nextedge);
    v1->nextedge = (EDGE **) 0;
    v1->num_nextedge = 0;
    v2->nextedge = new_edgelist;
    v2->num_nextedge = nedges;
  } else {
    free(v1->lastedge);
    free(v2->lastedge);
    v1->lastedge = (EDGE **) 0;
    v1->num_lastedge = 0;
    v2->lastedge = new_edgelist;
    v2->num_lastedge = nedges;
  }
}

/*
 * e1 is an edge on v1
 * if edge e1 will be equivalent to an edge e2 on v2 when we merge
 * v2,v1, which are
 *    the beginning vertices  (when out_flag=1)
 *    the ending vertices     (when out_flag=0)
 * then return e2; else return null
 *
 */
EDGE *zz_isedgedup(NODES *v1, NODES *v2, int out_flag, EDGE *e1)
{
  IVEDGE it_e2_mem, *it_e2;
  EDGE *e2;
  int e_type = out_flag ? E_OUT_H : E_IN_H;

  it_e2 = it_ev_renew(&it_e2_mem, v2, (E_TYPE) e_type);
  while ((e2 = it_ev_next(it_e2)) != (EDGE *) 0) {
    if (out_flag) {
      if (e1->end == e2->end)
	return e2;
    } else {
      if (e1->begin == e2->begin)
	return e2;
    }
  }
  return (EDGE *) 0;
}


void zz_deledge(NODES *v1, int out_flag, EDGE *e1)
{
  int i;
  NODES *v3;
  EDGE **edges3;
  int num_edges3;

  if (out_flag) {
    /* e1 is edge  v1  ->  v3
     * Delete it from v3's list
     */
    v3 = e1->end;
    num_edges3 = v3->num_lastedge;
    edges3 = v3->lastedge;

  } else {
    /* e1 is edge  v3  ->  v1
     * Delete it from v3's list
     */
    v3 = e1->begin;
    num_edges3 = v3->num_nextedge;
    edges3 = v3->nextedge;
  }

  for (i=0 ; i<num_edges3 ; i++) {
    if (edges3[i] == e1) {
#if 0
      edges3[i] = (EDGE *) 0;
#endif
      /* move final edge to this readposition, then decrease # edges by 1 */
      edges3[i] = edges3[num_edges3-1];
#if 0
      edges3[num_edges3-1] = (EDGE *) 0;
#endif
      if (out_flag)
	v3->num_lastedge--;
      else
	v3->num_nextedge--;

      return;
    }
  }

  printf("zz_deledge failure.  v1=%lx out_flag=%d e1=%lx\n", v1, out_flag, e1);

  for (i=0; i < v3->num_lastedge; i++)
    printf("  v3 in[%d]: %lx    other vertex: %lx\n", i, v3->lastedge[i], (v3->lastedge[i]) ? v3->lastedge[i]->begin : (NODES *) 0);

  for (i=0; i < v3->num_nextedge; i++)
    printf("  v3 out[%d]: %lx    other vertex: %lx\n", i, v3->nextedge[i], (v3->nextedge[i]) ? v3->nextedge[i]->begin : (NODES *) 0);

  //  exit(-1);
}


/****************************************************************************/
/****************************************************************************/

/****************************************************************************
 * create_test_zzpaths(G, RT, zz_type, u, L, count)
 *
 * create count zigzag paths,
 * where existing vertices x and y are connected by a new path of type u.
 * g is specified as array of 0's (<-) and 1's (->).
 *
 *
 * zz_type:
 * 1.       P--x  y--Q    =>   P  -- x  -u> y  -- Q
 *                             P' -- x' <u- y' -- Q'
 *     P != Q and P != Q'
 *
 * 3.       x--P--y       =>   P  -- x  -u>  y
 *                             P' <- x' <u'- y'
 *
 * 4.       x--P--y'  and  y--P'--x'    =>
 *
 *          x->P->y',    x -u> y,   x' <u- y'
 *
 *          x -> P -> y' -u> x' <- P' <- y <u- x
 *
 * 5.       P--x          =>   P -- x -u> x' -- P'      u self-complementary
 *
 * 6.       P--x          =>   P -> x -u> x
 *
 * 7.       P--x          =>   P -- x -u> x' -- P'
 *                                    <u-
 *
 * 8.       x, y isolated =>   x -u> y       y' <u- x'
 *
 * 9.       x isolated    =>   x -u> x'      u self-complementary
 *
 * 10.      x isolated    =>   x -u> x       x' <u- x'
 *
 * 11.      x isolated    =>   x -u> x'
 *                               <u-
 *
 * and complement in each case.
 ****************************************************************************/

void create_test_zzpaths(IGRAPH *G, READTABLE *RT,
			 int zz_type, char *u, int L, int count0)
{
  IVERT it_v_mem1, *it_v1;
  NODES *v1;

  IVERT it_v_mem2, *it_v2;
  NODES *v2;

  NODES *x, *y, *tmp;
  NODES *xc, *yc;
  COMP *comp_x;

  int count = count0;

  int self_dual;

  NODES **node_table;
  int i;
  NODES *w1,*w1c,*w2,*w2c;

  int sym;

  //  printf("create_test_zzpaths zz_type=%d u=%lx L=%d count=%d\n", zz_type, u, L, count);


  if (zz_type >= 8 && zz_type <= 11) {
    while (count > 0) {
      node_table = zz_allocnodes(G, (zz_type == 8) ? 4 : 2);
      if (node_table == (NODES **) 0) {
	printf("create_test_zzpath: Could not allocate nodes 8\n");
	return;
      }
      sym = 0;
      if (zz_type == 8) {
	x = node_table[0];
	xc = node_table[1];
	y = node_table[2];
	yc = node_table[3];
	free(node_table);

      } else {
	x = node_table[0];
	xc = node_table[1];
	free(node_table);

	if (zz_type == 9) {
	  y = xc; yc = x; sym = 1;
	} else if (zz_type == 10) {
	  y = x; yc = xc;
	} else if (zz_type == 11) {
	  y = xc; yc = x;
	}
      }

      x->bal_node = xc;
      xc->bal_node = x;
      y->bal_node = yc;
      yc->bal_node = y;

      x->subg_flag = xc->subg_flag = y->subg_flag = yc->subg_flag = SUBG_IN;

      zz_create_zzpath(G, x, y, sym, u, L);
      count--;
    }
    //  printf("create_zzpath: Type %d, requested %d, got %d\n", zz_type, count0, count0-count);
    return;
  }



  it_v1 = it_v_renew(&it_v_mem1, G, V_A_H);

  while (count > 0 &&
	 (v1 = it_v_next(it_v1)) != (NODES *) 0) {

    x = v1;

    y = (NODES *) 0;

    if (zz_type == 5 || zz_type == 7) {
      y = x->bal_node;
    }
    else if (zz_type == 6) {
      y = x;
    }
    else {
      /* find another vertex, satisfying appropriate case 1, 3, 4 */
      comp_x = comp_id(x);

      it_v2 = it_v_renew(&it_v_mem2, G, V_A_H);
      while (y == (NODES *) 0
	     && (v2 = it_v_next(it_v2)) != (NODES *) 0) {

	if (zz_type == 1) {
	  if (comp_id(v2) == comp_x ||
	      comp_id(v2->bal_node) == comp_x)
	    continue;
	  y = v2;
	} else if (zz_type == 3) {
	  if ( v2 != x && comp_id(v2) == comp_x)
	    y = v2;
	} else if (zz_type == 4) {
	  if ( v2 != x->bal_node && comp_id(v2->bal_node) == comp_x)
	    y = v2;
	}
      }
    }

    if (y != (NODES *) 0) {
      //      printf("Creating zzpath...\n");
      sym = (zz_type == 5);

      zz_create_zzpath(G, x, y, sym, u, L);
      count--;
    }
  }

  //  printf("create_zzpath: Type %d, requested %d, got %d\n", zz_type, count0, count0-count);
}

/* create zigzag path from x to y of type u
 * also create path from xc to yc
 *
 * sym:
 *    if x != yc, sym=0, and two separate paths are formed.
 *    if x == yc, set sym=0 for two paths or sym=1 for one self-symmetric path.
 *
 *    when sym=1, u must be symmetric.
 */

void zz_create_zzpath(IGRAPH *G, NODES *x, NODES *y, int sym, char *u, int L)
{
  NODES *xc, *yc;

  NODES **node_table;

  int i;
  NODES *w1, *w1c, *w2, *w2c;

  xc = x->bal_node;
  yc = y->bal_node;

  /* L edges on path
   * verts
   *   x -- nt[0] -- nt[1] -- ... nt[L-2] -- y
   *
   * and duals
   *
   * direction of nt[i-1] -- nt[i]
   * is ->  when u[i]=1
   *    <-  when u[i]=0
   */
	     
  if (!sym) {
    node_table = zz_allocnodes(G,2*(L-1));
    if (node_table == (NODES **) 0) {
      printf("create_test_zzpath: Could not allocate nodes 1\n");
      return;
    }

    for (i=0 ; i<L ; i++) {
      /* create ith edge on path and on complementary path */
      if (i==0) {
	w1 = x;
	w1c = x->bal_node;
      } else {
	w1 = node_table[2*(i-1)];
	w1c = node_table[2*(i-1)+1];
	w1->bal_node = w1c;
	w1c->bal_node = w1;
      }
      if (i==L-1) {
	w2 = y;
	w2c = y->bal_node;
      } else {
	w2 = node_table[2*i];
	w2c = node_table[2*i+1];
	w2->bal_node = w2c;
	w2c->bal_node = w2;
      }
      if (u[i])
	zz_createedge(w1,w2);
      else
	zz_createedge(w2,w1);
    }
    free(node_table);
  } else {
    /* symmetric path */

    /* error check input */
    if (x != yc) {
      printf("zigzag path is not symmetric, but is supposed to be: x=%lx yc=%lx\n", x, yc);
    }

    for (i=0 ; i<L ; i++) {
      if (u[i] != u[L-1-i]) {
	printf("zigzag pathtype is not symmetric: u[%d]=%d u[%d]=%d\n", i,u[i], L-1-i, u[L-1-i]);
      }
    }

    node_table = zz_allocnodes(G,L-1);

    if (node_table == (NODES **) 0) {
      printf("create_test_zzpath: Could not allocate nodes 2\n");
      return;
    }

    for (i=0 ; i<(L+1)/2 ; i++) {
      /* create ith edge on path
       * and complementary edge i from the end of the same path
       */
      if (i==0) {
	w1 = x;
	w1c = x->bal_node;
      } else {
	w1 = node_table[i-1];
	w1c = node_table[L-1-i];
	w1->bal_node = w1c;
	w1c->bal_node = w1;
      }

      if (i==L-1) {
	/* can't happen unless L is small */
	w2 = y;
	w2c = y->bal_node;
      } else {
	w2 = node_table[i];
	w2c = node_table[L-2-i];
	w2->bal_node = w2c;
	w2c->bal_node = w2;
      }

      if (u[i])
	zz_createedge(w1,w2);
      else
	zz_createedge(w2,w1);
    }
    free(node_table);
  }

}

/****************************************************************************
 * create_test_zzpath_short(G, RT, zz_type, L, count)
 *
 * create count zigzag paths
 *
 * types:   when L=0:
 *
 * 1.       P->x  y<-Q    =>   P -> {x,y} <- Q
 *     P != Q and P != Q'
 *
 * 2.       P->x  y->Q    =>   P -> {x,y} -> Q
 *     P != Q and P != Q'
 *
 * 3.       x<-P->y       =>   P -> {x,y}
 *
 * 4.       x->P->y'  and  y->P'->x'    =>
 *
 *                    -> P ->
 *           {x,y}             {y',x'}
 *                    -> P' ->
 *
 * 5.       P->x          =>   P ->{x,x'}-> P'
 *
 * and complement in each case.
 *
 *
 *
 * When L>0:
 *    u -L> v  =  connect u and v by L forwards edges
 * When L<0:
 *    u -L> v  =  connect u and v by |L| backwards edges
 * When L=0:
 *    u -L> v  =  merge u and v
 *
 * 1.       P->x  y<-Q    =>   P -> x -L> y <- Q
 *     P != Q and P != Q'
 *
 * 2.       P->x  y->Q    =>   P -> x -L> y -> Q
 *     P != Q and P != Q'
 *
 * 3.       x<-P->y       =>   P -> x -L> y
 *
 * 4.       x->P->y'  and  y->P'->x'    =>
 *
 *          x->P->y',    x -L> y,   x' <L- y'
 *
 * 5.       P->x          =>   P -> x -L> x' -> P'
 *
 * and complement in each case.
 *
 *
 * implement these now for L=0, 1, -1.
 ****************************************************************************/

void create_test_zzpath_short(IGRAPH *G, READTABLE *RT,
			      int zz_type, int L, int count0)
{
  IVERT it_v_mem1, *it_v1;
  NODES *v1;

  IVERT it_v_mem2, *it_v2;
  NODES *v2;

  NODES *x, *y, *tmp;
  COMP *comp_x;

  int count = count0;

  int self_dual;

  NODES **node_table;
  int i;
  NODES *w1,*w1c,*w2,*w2c;

  it_v1 = it_v_renew(&it_v_mem1, G, V_L_UH);

  while (count > 0 &&
	 (v1 = it_v_next(it_v1)) != (NODES *) 0) {

    /* get sink */
    if (gdegree(v1,E_IN_H) == 1)
      x = v1;
    else
      x = v1->bal_node;

    y = (NODES *) 0;

    if (zz_type == 5)
      y = x->bal_node;
    else {
      /* find another leaf, satisfying appropriate case 1-4 */
      comp_x = comp_id(x);

      it_v2 = it_v_renew(&it_v_mem2, G, V_L_UH);
      while (y != (NODES *) 0
	     && (v2 = it_v_next(it_v2)) != (NODES *) 0) {
	if (zz_type == 1 || zz_type == 2) {
	  if (comp_id(v2) == comp_x ||
	      comp_id(v2->bal_node) == comp_x)
	    continue;
	  if ( (gdegree(v2,E_IN_H) == 1) == (zz_type == 1) )
	    y = v2;
	  else
	    y = v2->bal_node;
	} else { /* zz_type == 3 || zz_type == 4 */
	  if ( v2 != x &&
	       ((gdegree(v2,E_IN_H) == 1) == (zz_type == 3)) )
	    y = v2;
	}
      }
    }

    if (y != x && y != (NODES *) 0) {
      printf("Creating zzpath...\n");
      self_dual =  x == y->bal_node;

      if (L == 0) {
	zz_mergenode(G, RT, x, y);
	comp_id_set(x,comp_id(y));
	comp_id_set(x->bal_node,comp_id(y->bal_node));

	if (self_dual)
	  y->bal_node = y;
      } else {
	if (L<0) {
	  tmp = x; x = y; y = tmp; L = -L;
	}
	/* for now, only do single edge */
	//	zz_createedge(x,y);

	if (self_dual) {
	  node_table = zz_allocnodes(G,L-1);
	  if (node_table == (NODES **) 0) {
	    printf("create_test_zzpath: Could not allocate nodes 3\n");
	    return;
	  }

	  /* L edges on path
	   * verts
	   *   x -> nt[0] -> nt[1] -> ... nt[L-2] -> y
	   * and duals
	   */
	     
	  for (i=0 ; i<L/2 ; i++) {
	    /* create ith edge on path */
	    if (i==0) {
	      w1 = x;
	      w1c = x->bal_node;
	    } else {
	      w1 = node_table[i-1];
	      w1c = node_table[L-2 - (i-1)];
	      w1->bal_node = w1c;
	      w1c->bal_node = w1;
	    }
	    if (i==L-1) {
	      w2 = y;
	      w2c = y->bal_node;
	    } else {
	      w2 = node_table[i];
	      w2c = node_table[L-2 - i];
	      w2->bal_node = w2c;
	      w2c->bal_node = w2;
	    }
	    zz_createedge(w1,w2);
	  }
	  free(node_table);
	
	} else {
	  node_table = zz_allocnodes(G,2*(L-1));
	  if (node_table == (NODES **) 0) {
	    printf("create_test_zzpath: Could not allocate nodes 4\n");
	    return;
	  }

	  /* L edges on path
	   * verts
	   *   x -> nt[0] -> nt[1] -> ... nt[L-2] -> y
	   * and duals
	   */
	     
	  for (i=0 ; i<L ; i++) {
	    /* create ith edge on path */
	    if (i==0) {
	      w1 = x;
	      w1c = x->bal_node;
	    } else {
	      w1 = node_table[2*(i-1)];
	      w1c = node_table[2*(i-1)+1];
	      w1->bal_node = w1c;
	      w1c->bal_node = w1;
	    }
	    if (i==L-1) {
	      w2 = y;
	      w2c = y->bal_node;
	    } else {
	      w2 = node_table[2*i];
	      w2c = node_table[2*i+1];
	      w2->bal_node = w2c;
	      w2c->bal_node = w2;
	    }
	    zz_createedge(w1,w2);
	  }
	  free(node_table);
	}

      }

      count--;
    }
  }

  printf("create_zzpath: Type %d, requested %d, got %d\n", zz_type, count0, count0-count);
}

/* create edges x->y and yc->xc */
void zz_createedge(NODES *x, NODES *y)
{
  NODES *xc = x->bal_node;
  NODES *yc = y->bal_node;
  int self_dual = x == yc;

  EDGE *e, *ec;

  /* create edges */
  e = (EDGE *) ckalloc(1 * sizeof(EDGE));
  if (!self_dual)
    ec = (EDGE *) ckalloc(1 * sizeof(EDGE));
  else
    ec = e;

  e->begin = x;
  e->end = y;
  e->seq = 0;
  e->start_cover = 1;
  e->length = 2;
  e->multip = 0;
  e->bal_edge = ec;
  e->readinterval = (READINTERVAL *) 0;
  e->visit = 0;
  e->subg_flag = SUBG_IN;

  /* add e into the edge lists of the vertices */
  zz_createedge_add(x,e,1);
  zz_createedge_add(y,e,0);

  if (!self_dual) {
    ec->begin = yc;
    ec->end = xc;
    ec->seq = 0;
    ec->start_cover = 1;
    ec->length = 2;
    ec->multip = 0;
    ec->bal_edge = e;
    ec->readinterval = (READINTERVAL *) 0;
    ec->visit = 0;
    ec->subg_flag = SUBG_IN;

    /* add ec into the edge lists of the vertices */
    zz_createedge_add(yc,ec,1);
    zz_createedge_add(xc,ec,0);
  }
}

/* add edge e to the edgelist of v
 *    out_flag=0: incoming edge
 *    out_flag=1: outgoing edge
 */
void zz_createedge_add(NODES *v, EDGE *e, int out_flag)
{
  int *num_edges;
  EDGE ***edge_list;

  EDGE **new_edge_list;

  if (out_flag) {
    num_edges = &v->num_nextedge;
    edge_list = &v->nextedge;
  } else {
    num_edges = &v->num_lastedge;
    edge_list = &v->lastedge;
  }

  new_edge_list = (EDGE **) ckrealloc(*edge_list,
				      (*num_edges+1) * sizeof(EDGE *),
				      (*num_edges) * sizeof(EDGE *));
  /* append edge to list of edges */
  new_edge_list[*num_edges] = e;

  ++*num_edges;
  *edge_list = new_edge_list;
}


/* allocate n nodes
 * limitation: must reuse deleted nodes from G
 * returns array of pointers to nodes
 */
NODES **zz_allocnodes(IGRAPH *G, int n)
{
  IVERT it_v_mem, *it_v;
  NODES **node_table;
  NODES *v;
  int i = 0;


  node_table = (NODES **) ckalloc(n * sizeof(NODES *));

  it_v = it_v_renew(&it_v_mem, G, V_A_G);

  while (i < n
	 && ((v = it_v_next(it_v)) != (NODES *) 0)) {
    if (!v->subg_flag) {
      node_table[i] = v;
      i++;

      v->subg_flag = SUBG_IN;
    }
  }

  if (i<n)
    return (NODES **) 0;
  return node_table;
}
