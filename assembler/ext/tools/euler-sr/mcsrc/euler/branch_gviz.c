/***************************************************************************
 * Title:          branch_gviz.c
 * Author:         Glenn Tesler
 * Created:        Dec. 2002
 * Last modified:  Dec. 2002
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

/* Output branching graph of H in .gviz format
 * Leaves G and H intact
 * Glenn Tesler
 * 12/19/02
 */

#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>
//#include <subgraph.h>
#include <branch_gviz.h>


//#define DEBUG_BG
#undef DEBUG_BG



/* The label for a vertex is the 5 least sig hex digits.
 * WARNING: There's a low probability event that this is not unique.
 */
//#define V_ID(v) ((int) ((((unsigned long int) v)/sizeof(NODES)) & 0xFFFFF))
#define V_ID(v) v


void write_bgraph_gviz_file(IGRAPH *G, char *filename)
{
  FILE *fp;

#ifdef DEBUG_BG
  fprintf(stderr,"DEBUG: write_bgraph_gviz_file(G,\"%s\")\n", filename);
#endif

  fp = ckopen(filename, "w");
  write_bgraph_gviz(G, fp);
  fclose(fp);


#ifdef DEBUG_BG
  fprintf(stderr,"Graphviz file %s: G: %d verts, %d nuc-edges; H: %d verts, %d nuc-edges, %d dirbranch verts, %d undir branch verts, %d leaves\n",
	 filename,
	 it_v_count(G,V_A_G), it_e_count(G,E_OUT_G),
	 it_v_count(G,V_A_H), it_e_count(G,E_OUT_H),
	 it_v_count(G,V_B_H), it_v_count(G,V_B_UH), it_v_count(G,V_L_UH));
#endif
}

void write_bgraph_gviz(IGRAPH *G, FILE *fp)
{
  IVERT *it_vert;
  IVEDGE it_edge_mem, *it_edge;

  NODES *v1, *v2;
  EDGE *e1, *e2;
  int num_paths, num_edges, len_edges;

  int cov_min, cov_max, cov_tot;

  int dbpath = 0;

  int pass;


  /* graphviz header */
  fprintf(fp,"digraph G {\n");
  fprintf(fp,"\tsize=\"8,8\";\n");

  /* Detect loops
   * 1. Mark all vertices unvisited
   * 2. Draw branching graph, going through unvisited vertices and
   *    marking the ones encountered as visited.
   * 3. Search for unvisited vertices and draw the cycles.
   */

  clear_visit_flag(G, V_A_H);

  /* Pass 1. loop over all branching vertices
   * Pass 2. loop over all verts to find unvisited ones; this gives
   *         components that are cycles
   */

  for (pass=1; pass<=2; pass++) {
    if (pass == 1)
      it_vert = it_v_new(G, V_B_H);
    else
      it_vert = it_v_new(G, V_A_H);

    while ((v1 = it_v_next(it_vert)) != (NODES *) 0) {
      if (pass==2 && v1->visit)
	continue;

      //    fprintf(stderr,"wbg DEBUG: v1=%lx\n", v1);
      /* follow each outgoing path from v1 until hit another branching vert */
      it_edge = it_ev_renew(&it_edge_mem, v1, E_OUT_H);

      /* # of outgoing paths */
      num_paths = 0;



#ifdef DEBUG_BG
      dbpath = (v1 == (NODES *) 0x148cf0180 ||
		v1 == (NODES *) 0x148cf01e0 ||
		v1 == (NODES *) 0x148cedd00);
#endif

      while ((e1 = it_ev_next(it_edge)) != (EDGE *) 0) {
	//      fprintf(stderr,"wbg DEBUG:   e1=%lx\n", e1);
	e2 = e1;
	v2 = v1;

#ifdef DEBUG_BG
	if (dbpath) {
	  fprintf(stderr,"DEBUG path: e1=%lx ",e1);
	}
#endif


	/* follow the path from e1 along outgoing edges until hit another
	 * branching vertex
	 */

	/* statistics */
	num_edges = 0;            /* number of edges in path */
	len_edges = 0;            /* cumulative length of edges */
	cov_min   = 1000000;      /* minimum coveraage */
	cov_max   = 0;            /* maximum coverage */
	cov_tot   = 0;            /* cumulative coverage */ 
      
	do {
	  v2->visit = 1;
#ifdef DEBUG_BG
	  if (dbpath) {
	    fprintf(stderr,"v=%lx e=%lx ",v2, e2);
	  }
#endif

	  num_edges++;
	  len_edges += e2->length - 1;

	  if (e2->multip < cov_min)
	    cov_min = e2->multip;
	  if (e2->multip > cov_max)
	    cov_max = e2->multip;
	  cov_tot += e2->multip;

	  /* vertex may have self-loop, abort */
	  if (e2->end == v2)
	    break;

	  v2 = e2->end;

	  /* isolated cycles */
	  if (pass==2 && v2==v1)
	    break;
	} while ((e2 = find_unique_oedge(v2)) != (EDGE *) 0);

#ifdef DEBUG_BG
	if (dbpath) {
	  fprintf(stderr,"; v1=%lx v2=%lx\n",v1, v2);
	}
#endif


	/* Draw edge from v1 to v2 */
	fprintf(fp,
		//	      "\t\"%05x\" -> \"%05x\"",
		"\t\"%lx\" -> \"%lx\"",
		V_ID(v1), V_ID(v2));

#if 0
	if (!(num_edges == 1 && len_edges == 1)) {
	  if (num_edges == len_edges) {
	    fprintf(fp,
		    " [label=\"%d\"]",
		    num_edges);
	  } else {
	    fprintf(fp,
		    " [label=\"%d,%d\"]",
		    num_edges, len_edges);
	  }
	}
#endif

	fprintf(fp, " [label=\"");
	if (num_edges == len_edges) {
	  if (num_edges > 1)
	    fprintf(fp, "%d", num_edges);
	} else {
	  fprintf(fp, "%d,%d", num_edges, len_edges);
	}

	fprintf(fp, "(");

	if (cov_min == cov_max)
	  fprintf(fp, "%d", cov_min);
	else
	  fprintf(fp, "%d,%d,%d", cov_min, cov_tot/num_edges, cov_max);

	fprintf(fp, ")\"]");
      

	fprintf(fp, ";\n");
      }


      /* if there were no outgoing paths, verify that there is
       * at least one incoming path.  If not, draw isolated vertex.
       */
      if (num_paths == 0 &&
	  gdegree_cmp(v1,E_IN_H,0) == 0) {
	//      fprintf(fp, "\t\"%05x\";\n", V_ID(v1));
      }
    }
    it_v_destroy(it_vert);
  }

  /* TODO:
   * show isolated loops (v has edge v->v and no other edge)
   * only show one component of each symmetric pair
   */

  /* graphviz trailer */
  fprintf(fp,"}\n");
}


