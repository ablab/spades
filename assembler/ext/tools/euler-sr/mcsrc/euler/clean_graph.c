/***************************************************************************
 * Title:          cleangraph.c
 * Author:             Glenn Tesler
 * Created:        Dec. 2002
 * Last modified:  Dec. 2002
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

/* graph cleaning algorithm
 * Glenn Tesler
 * 12/18/02
 */

#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>
#include <clean_graph.h>
#include <zigzag.h>
#include <erosion.h>

//#define DEBUG_OUT
#undef DEBUG_OUT

#define DEBUG_GVIZ
//#undef DEBUG_GVIZ

#define PARAMS_GROUP 1

#ifdef DEBUG_GVIZ
#define debug_gviz(G,filename) write_bgraph_gviz_file(G,filename)
#else
#define debug_gviz(G,filename)
#endif




/* For now, hard code the parameters.
 * Later, feed them in correctly.
 */
void clean_graph_adapter(IGRAPH *G, READTABLE *RT)
{
  CLEAN_PARAMS params;
  LABELVERT T_e_mem, T_c_mem;

  params.T_e = &T_e_mem;
  params.T_e->onepath = 1;
  params.T_e->dir = 0;


  params.T_c = &T_c_mem;
  params.T_c->onepath = 1;
  params.T_c->dir = 0;

  /* DO NOT MODIFY THESE PARAMETERS HERE.
   * MAKE A COPY OF THIS AND MODIFY THE COPY.
   */
#if PARAMS_GROUP==1
  /* DO NOT MODIFY THESE */
  params.w = 5;
  params.max_cov = 16;
  params.C_B = 60;
  params.C_W = 60;
  params.T_c->depth = 100;
  params.T_e->depth = 100;
  params.L_b = 5;
  params.L_c = 1;
  params.L_e = 1;
#elif PARAMS_GROUP==2
  /* DO NOT MODIFY THESE */
  params.w = 5;
  params.max_cov = 16;
  params.C_B = 150;
  params.C_W = 150;
  params.T_c->depth = 100;
  params.T_e->depth = 100;
  params.L_b = 5;
  params.L_c = 1;
  params.L_e = 1;
#elif PARAMS_GROUP==3
  /* DO NOT MODIFY THESE */
  params.w = 5;
  params.max_cov = 16;
  params.C_B = 60;
  params.C_W = 60;
  params.T_c->depth = 100;
  params.T_e->depth = 100;
  params.L_b = 1;
  params.L_c = 0;
  params.L_e = 1;
#elif PARAMS_GROUP==4
  /* DO NOT MODIFY THESE */
  params.w = 5;
  params.max_cov = 1;
  params.C_B = 30;
  params.C_W = 30;
  params.T_c->depth = 100;
  params.T_e->depth = 100;
  params.L_b = 1;
  params.L_c = 0;
  params.L_e = 1;
#elif PARAMS_GROUP==5
  /* DO NOT MODIFY THESE */
  params.w = 5;
  params.max_cov = 1;
  params.C_B = 30;
  params.C_W = 30;
  params.T_c->depth = 100;
  params.T_e->depth = 100;
  params.L_b = 1;
  params.L_c = 1;
  params.L_e = 1;
#elif PARAMS_GROUP==6
  /* DO NOT MODIFY THESE */
  params.w = 5;
  params.max_cov = 1;
  params.C_B = 30;
  params.C_W = 30;
  params.T_c->depth = 100;
  params.T_e->depth = 0;
  params.L_b = 1;
  params.L_c = 1;
  params.L_e = 1;
#elif PARAMS_GROUP==7
  /* DO NOT MODIFY THESE */
  params.w = 5;
  params.max_cov = 16;
  params.C_B = 60;
  params.C_W = 60;
  params.T_c->depth = 5000;
  params.T_e->depth = 100;
  params.L_b = 5;
  params.L_c = 1;
  params.L_e = 1;
#else
  /* TO MODIFY THESE, PUT IN ANOTHER #elif HERE */
  params.w = 5;
  params.max_cov = 16;
  params.C_B = 60;
  params.C_W = 60;
  params.T_c->depth = 100;
  params.T_e->depth = 100;
  params.L_b = 5;
  params.L_c = 1;
  params.L_e = 1;
#endif

  clean_graph_once(G, RT, &params);
//  if (RT -> chim != (int *) 0)
//    free((void *) RT -> chim);
}

/* debugging */

static char zz1[] = { 1,1,1,1,0,0,1,0,1,1,0,0,0,1,1 };
int L1 = 15;

/* symmetric, middle is odd length backwards */
static char zz2[] = { 1,1,1,1,0,0,1,1,0,0,0,1,1,0,0,1,1,1,1 };
int L2 = 19;

/* symmetric, middle is odd length backwards */
static char zz3[] = { 1,1,1,1,1,0,0,0,1,1,1,1,1 };
int L3 = 13;

/* symmetric, middle is even length forwards */
static char zz4[] = { 1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1 };
int L4 = 30;

/* symmetric, middle is even length backwards */
static char zz5[] = { 1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1 };
int L5 = 18;

/* symmetric, middle is even length backwards */
static char zz6[] = { 1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1 };
int L6 = 34;

static char zz7[] = { 1,1,1,1,0,1,1,1,1 };
int L7 = 9;

static char zz8[] = { 1,1,1,1,0,0,1,1,1,1 };
int L8 = 10;

static char zz9[] = { 1,1,1,1,0,1,0,1,1,1,1 };
int L9 = 11;

static char zz10[] = { 1,1,1,1,0,1,1,0,1,1,1,1 };
int L10 = 12;




void clean_graph_once(IGRAPH *G, READTABLE *RT, CLEAN_PARAMS *params)
{
#if 0
  int i;
#endif

  /* mark chimeric reads */
  RT -> chim = (int *) 0;
  RT -> num_chim = 0;
  RT -> num_chim_alloc = 0;



  create_full_subgraph(G, SUBG_IN, 0);
  debug_gviz(G,"G_init.gvz");

  printf("Degrees in initial graph\n");
  debug_print_degrees(G);
  fflush(stdout);



  /* Form subgraph of G with all vertices but no edges */
  create_full_subgraph(G, SUBG_POSTPONE, 1);

  /* Form maximum spanning forest of G */
  graph_grow(G,RT,params,0);

#ifdef DEBUG_OUT
  //  debug_print_vval(G,RT);
  //  debug_check_sym(G);
#endif

  debug_gviz(G,"G_forest.gvz");

  fflush(stdout);
  printf("Degrees in forest\n");
  debug_print_degrees(G);
  fflush(stdout);

  /* Add additional edges, provided they don't form bulges
   * or come from chimeric reads
   */
  graph_grow(G,RT,params,1);
  fflush(stdout);

#ifdef DEBUG_OUT
  //  debug_print_vval(G,RT);
  //  debug_check_sym(G);
#endif

  debug_gviz(G,"G_graph.gvz");

  printf("Degrees in graph with cycles\n");
  debug_print_degrees(G);
  fflush(stdout);

  /* Erosion of trees hanging off of deep components */
  erode_graph(G, params);
  fflush(stdout);

#ifdef DEBUG_OUT
  fprintf(stderr,"DEBUG: calling debug_gviz(G,\"G_erosion.gvz\")\n");
#endif

#ifdef DEBUG_OUT
  //  debug_print_vval(G,RT);
  //  debug_check_sym(G);
#endif

  debug_gviz(G,"G_erosion.gvz");
  fflush(stdout);

  printf("Degrees after erosion\n");
  debug_print_degrees(G);
  fflush(stdout);

#ifdef DEBUG_OUT
  //  debug_print_vval(G,RT);
  //  debug_check_sym(G);
#endif


#if 0
  /* create synthetic zigzag paths for testing purposes */

  /* Mark degree 0 nodes as deleted so that their memory can be reused */
  mark_deleted_nodes(G);

#define zz_sym_paths(zz,LL) \
  create_test_zzpaths(G, RT, 1, zz, LL, 1); \
  create_test_zzpaths(G, RT, 3, zz, LL, 1); \
  create_test_zzpaths(G, RT, 4, zz, LL, 1); \
  create_test_zzpaths(G, RT, 5, zz, LL, 1); \
  create_test_zzpaths(G, RT, 6, zz, LL, 1); \
  create_test_zzpaths(G, RT, 7, zz, LL, 1);

#define zz_cyc_paths(zz,LL) \
  create_test_zzpaths(G, RT, 8, zz, LL, 1); \
  create_test_zzpaths(G, RT, 9, zz, LL, 1); \
  create_test_zzpaths(G, RT, 10, zz, LL, 1); \
  create_test_zzpaths(G, RT, 11, zz, LL, 1);

  create_test_zzpaths(G, RT, 1, zz1, L1, 1);
  create_test_zzpaths(G, RT, 3, zz1, L1, 1);
  create_test_zzpaths(G, RT, 4, zz1, L1, 1);
  create_test_zzpaths(G, RT, 6, zz1, L1, 1);
  create_test_zzpaths(G, RT, 7, zz1, L1, 1);

  zz_sym_paths(zz2,L2);
  zz_sym_paths(zz3,L3);
  zz_sym_paths(zz4,L4);
  zz_sym_paths(zz5,L5);
  zz_sym_paths(zz6,L6);
  zz_sym_paths(zz7,L7);
  zz_sym_paths(zz8,L8);
  zz_sym_paths(zz9,L9);
  zz_sym_paths(zz10,L10);

  zz_cyc_paths(zz2,L2);
  zz_cyc_paths(zz3,L3);
  zz_cyc_paths(zz4,L4);
  zz_cyc_paths(zz5,L5);
  zz_cyc_paths(zz6,L6);
  zz_cyc_paths(zz7,L7);
  zz_cyc_paths(zz8,L8);
  zz_cyc_paths(zz9,L9);
  zz_cyc_paths(zz10,L10);
#endif

  fflush(stdout);
  if (params->do_zz)
    process_all_zzpaths(G, RT);

#ifdef DEBUG_OUT
  //  debug_print_vval(G,RT);
  //  debug_check_sym(G);
#endif
  debug_gviz(G,"G_zigzag.gvz");

  /* Mark degree 0 nodes as deleted */
  mark_deleted_nodes(G);

  /* Map reads back onto graph, form paths,
   * delete/split reads that don't map well
   *
   * DONE IN over-repeat-new.c INSTEAD
   */

  //  debug_gviz(G,"G_mapped.gvz");



  /* later:
   * do above all over again w/reduced # reads, changed parameters */
  /* TODO */

#if 0
  printf("Chimeric read indices: ");
  for (i=0; i<RT->num_chim; i++) {
    printf(" %d", RT->chim[i]);
  }
  printf("\n");
#endif
  fflush(stdout);
}

/* grow graph by a modified Kruskal algorithm
 * treerule=0: build a forest.  Discard obvious bad edges.  Postpone others.
 * treerule=1: may form cycles.  Discard bad edges, insert others.
 */           
void graph_grow(IGRAPH *G,
		READTABLE *RT,
		CLEAN_PARAMS *params,
		int treerule)
{
  int cov;        /* current coverage being considered */
  int cov2;       /* next highest coverage */
  int mul;
  int max_cov = params->max_cov;

  EDGE *e,*ec;
  NODES *v1,*v2,*v1c,*v2c;
  NODES *comp1, *comp2, *comp1c, *comp2c;
  int internal_flag;

  int has_cycle;

  int computed_external = 0;

  IEDGE *it_edge;

  int is_bulge;

  /* statistics */
  int num_ins = 0;        /* number of edges inserted */
  int num_post = 0;       /* number postponed */
  int num_del = 0;        /* number deleted */
  int num_del_chim = 0;   /* number chimeric read edges deleted */
  int num_del_bulge = 0;  /* number of bulge edges deleted */
  int num_ins_tree = 0;   /* number of tree edges inserted */
  int num_ins_cycle = 0;  /* number of cycle edges inserted */

  int inc_edge;           /* usually above are incremented by 2,
			   * but are inc'd by 1 if e==ec
			   */

  /* subtotals up to prior coverage considered */
  int num_ins0 = 0;
  int num_post0 = 0;
  int num_del0 = 0;
  int num_del_chim0 = 0;
  int num_del_bulge0 = 0;
  int num_ins_tree0 = 0;
  int num_ins_cycle0 = 0;


  /* number of chimeric reads deleted in total */
  int num_del_chim_reads0 = RT->num_chim; 

  int do_del;

  /* at low coverage, do 2 passes, for internal/external vertex combinations:
   * pass 1: examine edges connecting ext->ext, ext->int, int->ext
   * pass 2: examine edges connecting int->int
   *
   * pass "0": not at low coverage,
   *           OR at low coverage and need to initialize pass counter
   */

  int low_pass = 0;


  /***********************************************************************
   * statistics
   ***********************************************************************/

  printf("----------------------------------------------------------------------------------------\n");
  printf("Grow graph pass %d\n", treerule+1);
  printf("Classification of edges at each coverage\n");

    
  printf("                                            deleted......... inserted.........\n");
  printf("Cov  inserted  postponed  deleted    total     chim whrl/blg     tree    cycle\n");




  /***********************************************************************
   * loop over edges in decreasing order of coverage
   ***********************************************************************/

  /* We clip the edge coverage so that we may examine them in
   * decreasing order of coverage by doing a small # of passes through all
   * edges
   */

  for (cov=max_cov, cov2=-1 ; cov>=1 ; cov = cov2, cov2=-1) {
#ifdef DEBUG_OUT
    fprintf(stderr,"Starting cov=%d\n", cov);
#endif

    if (!computed_external &&
	cov <= params->L_c) {
      /* classify nodes as internal/external */
      classify_internal_nodes(G, params->T_c);
      computed_external = 1;

      /* two passes, mark the ones seen on first pass */
      clear_visit_edge(G, E_P);
    }
    if (computed_external) {
      low_pass++;
    }

#ifdef DEBUG_OUT
    for (mul=0; mul<3; mul++) {
      switch (mul) {
      case 0:
	e = (EDGE *) 0x1464280f0; break;
      case 1:
	e = (EDGE *) 0x146428050; break;
      case 2:
	e = (EDGE *) 0x14641c7b0; break;
      }
      fprintf(stderr,"  cov=%d; Monitoring edge e=%lx  %lx(%d,%d)->%lx(%d,%d)  subg=%d  multip=%d\n", cov, e, e->begin, e->begin->visit_value, e->begin->ext_flag, e->end, e->end->visit_value, e->end->ext_flag, e->subg_flag, e->multip);
    }
#endif


    it_edge = it_e_new(G, (E_TYPE) (E_P | E_S));

    while ((e = it_e_next(it_edge)) != (EDGE *) 0) {
      /* check edge multiplicity, and determine next highest multiplicity */
      mul = e->multip;

#ifdef DEBUG_OUT
	//	fprintf(stderr,"checking edge: mul=%d\n", mul);
#endif

#ifdef DEBUG_OUT
      if (mul <= 0)
	fprintf(stderr,"INVALID MULTIPLICITY ON EDGE e=%lx, v1=%lx, v2=%lx, mul=%d\n", e, v1, v2, mul);
#endif

      /* determine next highest multiplicity */
      if (mul < cov && mul > cov2)
	cov2 = mul;

      /* edge doesn't have the multiplicity we are currently examining */
      if (!(cov == mul ||
	    (cov == max_cov && mul >= cov)))
	continue;

      /* vertices in edge e:  v1, v2
       * in complementary edge ec: v1c, v2c
       * components of these: comp1, comp2, comp1c, comp2c
       */
      v1 = e->begin;
      v2 = e->end;

#ifdef DEBUG_OUT
      fprintf(stderr,"PROCESSING EDGE e=%lx, %lx->%lx, mul=%d; ec=%lx, %lx->%lx\n", e, v1, v2, mul, e->bal_edge, e->bal_edge->begin, e->bal_edge->end);
#endif

      /* At low coverage, check internal/external status of verts.
       * Use the same test as for chimeric edges later.
       */
      if (low_pass) {
	if (e->visit)
	  continue;

	internal_flag = is_edge_chimeric(e, params->L_c);

#if 0
	internal_flag =
	  //	  !v1->ext_flag && !v2->ext_flag;
	  is_path_internal(v1, E_IN_H, params->L_c)
	  && is_path_internal(v2, E_OUT_H, params->L_c);
#endif

	if (low_pass == 1 && internal_flag) {
	  /* 1st low pass: avoid internal-internal */
#ifdef DEBUG_OUT
	  fprintf(stderr,"   EDGE DELAYED e=%lx, v1=%lx, v2=%lx, ext1=%d, ext2=%d, low_pass=%d\n", e, v1, v2, v1->ext_flag, v2->ext_flag, low_pass);
#endif
	  continue;
	}
	e->visit = 1;

#if 0
	if (internal_flag) {
//	if (!v1->ext_flag && !v2->ext_flag) {
	  if (low_pass == 1) {
#ifdef DEBUG_OUT
	    fprintf(stderr,"   EDGE DELAYED e=%lx, v1=%lx, v2=%lx, ext1=%d, ext2=%d, low_pass=%d\n", e, v1, v2, v1->ext_flag, v2->ext_flag, low_pass);
#endif
	    continue;
	  }
	} else {
	  if (low_pass != 1) {
#ifdef DEBUG_OUT
	    fprintf(stderr,"   EDGE DELAYED e=%lx, v1=%lx, v2=%lx, ext1=%d, ext2=%d, low_pass=%d\n", e, v1, v2, v1->ext_flag, v2->ext_flag, low_pass);
#endif
	    continue;
	  }
	}
#endif


      }


      /******************************************************************
       * end of loop control statements
       * begin loop body
       ******************************************************************/


      /* complementary edge */
      ec = e->bal_edge;
      inc_edge = (e == ec) ? 1 : 2;
#ifdef DEBUG_OUT
      //      fprintf(stderr,"   ec=%lx\n", ec);
#endif


      /* check if a cycle through e exists in
       * H union {e, ec}
       */

      /* first check if H union {e} has a cycle */
      comp1 = comp_id(v1);
      comp2 = comp_id(v2);
      has_cycle = comp1 == comp2;

      /* it doesn't, but maybe H union {e,e^c} does */
      if (!has_cycle && e != ec) {

	v1c = ec->begin;
	v2c = ec->end;
	comp1c = comp_id(v1c);
	comp2c = comp_id(v2c);

	has_cycle =
	  ((comp1 == comp1c && comp2 == comp2c) ||
	   (comp1 == comp2c && comp2 == comp1c));
      }



      if (has_cycle) {
	if (cov <= params->L_b) {
	  /* if cycle length <= B
	   * then the edge forms a bulge, delete the edge.
	   */

	  /* tentatively put in edge ec but not e,
	   * then check if path length between verts of e
	   * as described below.
	   * In the odd case that e=ec or e=reverse edge of ec,
	   * the graph has a loop and it would be a cycle of length 1,
	   * so we don't need to exclude that possibility.
	   */
	  ec->subg_flag = SUBG_IN;


	  /* new test for bulges/whirls separately:
	   * zigzag cycle length <= C_B is a bulge
	   *   meaning there exists undirected cycle of len <= C_B
	   *   and there does not exist directed cycle of len <= C_B
	   * directed cycle length <= C_W is a whirl
	   *
	   * most efficient test depends on values of C_B, C_W
	   *
	   * TODO: figure out which order to do each test in to be
	   * most efficient
	   */

	  /* simplest but not necessarily efficient: */
#if 0
	  is_bulge =
	    /* check bulge */
	    (check_path_length(G, v1, v2, params->C_B, E_U_H)
	     &&	!check_path_length(G, v1, v2, params->C_B, E_IN_H))
	    ||
	    /* check whirl */
	    check_path_length(G, v1, v2, params->C_W, E_IN_H);
#endif

	  if (params->C_W <= 0) {
	    if (params->C_B > 0) {
	      is_bulge =
		check_path_length(G, v1, v2, params->C_B, E_U_H)
		&& !check_path_length(G, v1, v2, params->C_B, E_IN_H);
	    } else {
	      is_bulge = 0;
	    }

	  } else if (params->C_B <= 0) {
	    /* check whirl */
	    is_bulge = check_path_length(G, v1, v2, params->C_W, E_IN_H);

	  } else if (params->C_B == params->C_W) {
	    is_bulge = check_path_length(G, v1, v2, params->C_B, E_U_H);

	  } else if (params->C_B < params->C_W) {
	    /* TODO: optimize order of tests & eliminate redundancy */
	    is_bulge =
	      /* check bulge */
	      (check_path_length(G, v1, v2, params->C_B, E_U_H)
	       && !check_path_length(G, v1, v2, params->C_B, E_IN_H))
	      ||
	      /* check whirl */
	      check_path_length(G, v1, v2, params->C_W, E_IN_H);

	  } else {  // (params->C_B > params->C_W)
	    /* TODO: optimize order of tests & eliminate redundancy */
	    is_bulge =
	      /* check bulge */
	      (check_path_length(G, v1, v2, params->C_B, E_U_H)
	       && !check_path_length(G, v1, v2, params->C_B, E_IN_H))
	      ||
	      /* check whirl */
	      check_path_length(G, v1, v2, params->C_W, E_IN_H);
	  }
	  

	  //	  if (check_path_length(G, v1, v2, params->C_B)) {
	  if (is_bulge) {
	    /* Is a bulge or a whirl.
	     * Delete the edges e, ec.
	     */
	    subgraph_del_edge(e);
	    if (e != ec)
	      subgraph_del_edge(ec);
	    num_del += inc_edge;
	    num_del_bulge += inc_edge;
	    continue;

	  } else {
	    /* does not have cycle of length <= B */
	    /* undo the tentative insertion of ec;
	     * its fate is still to be determined
	     */
	    ec->subg_flag = SUBG_POSTPONE;
	  }
	}


	/* H union {e,e^c} has a cycle through e, but not a bulge. */
	if (treerule == 0) {
	  /* postpone the edge */
	  /* edge status already is postponed, nothing more to do */
	  num_post += inc_edge;
	  continue;
	}
      }

#ifdef DEBUG_OUT
      //	fprintf(stderr, "test chimeric edge, cov %d  ext flags %d %d\n", cov, v1->ext_flag, v2->ext_flag);
#endif

      /* chimeric read detection
       * A read with this edge is considered chimeric
       * if this edge has low coverage
       * and at least one of its vertices is an internal vertex.
       */
      if (cov <= params->L_c &&
	  internal_flag
//	  !v1->ext_flag && !v2->ext_flag
	  ) {

	do_del = RT->num_chim;

	/* The read containing e is chimeric,
	 * mark it and its complement for deletion
	 */
	subgraph_del_read_with_edge(RT,e);
#ifdef DEBUG_OUT
	fprintf(stderr, "chimeric edge, ext flags %d %d %d\n", v1->ext_flag, v2->ext_flag, RT->num_chim-do_del);
#endif


	/* delete e, ec */
	subgraph_del_edge(e);
	if (e != ec)
	  subgraph_del_edge(ec);

	num_del += inc_edge;
	num_del_chim += inc_edge;

	continue;
      }


      /* insert edges e, e^c, and join components */
      subgraph_ins_edge(e);
      if (e != ec)
	subgraph_ins_edge(ec);
      num_ins += inc_edge;

      if (has_cycle) {
	num_ins_cycle += inc_edge;
      } else {
	num_ins_tree += inc_edge;
      }
    }
    it_e_destroy(it_edge);



    /* report statistics for this coverage */
    if (cov == max_cov) {
      printf(">");
    } else {
      printf(" ");
    }
    printf("%2d  %8d   %8d %8d %8d %8d %8d %8d %8d\n",
	   cov,
	   num_ins-num_ins0, num_post-num_post0, num_del-num_del0,
	   (num_ins-num_ins0)+(num_post-num_post0)+(num_del-num_del0),
	   num_del_chim-num_del_chim0, num_del_bulge-num_del_bulge0,
	   num_ins_tree-num_ins_tree0, num_ins_cycle-num_ins_cycle0);

    /* update subtotals for each statistic */
    num_ins0 = num_ins;
    num_post0 = num_post;
    num_del0 = num_del;
    num_del_chim0 = num_del_chim;
    num_del_bulge0 = num_del_bulge;
    num_ins_tree0 = num_ins_tree;
    num_ins_cycle0 = num_ins_cycle;


    /* at low coverage, set up for next pass */
    if (low_pass) {
      if (low_pass == 1) {
	cov2 = cov;      /* do another pass at this coverage */
	//	if (treerule)
	//	  cov2 = cov;      /* do another pass at this coverage */
	//	else
	//	  low_pass = 0;
      } else {
	low_pass = 0;    /* done with this coverage,
			  * reset pass # for next coverage */
      }
    }
  }


  /* clear visit flag */
  clear_visit_edge(G, E_G);


  /* report statistics */
  printf("Tot  %8d   %8d %8d %8d %8d %8d %8d %8d\n",
	   num_ins, num_post, num_del,
	   num_ins+num_post+num_del,
	   num_del_chim, num_del_bulge,
	   num_ins_tree, num_ins_cycle);
  printf("\n");
  printf("%d chimeric reads detected\n",
	 RT->num_chim - num_del_chim_reads0);
}
