/***************************************************************************
 * Title:          subgraph.h
 * Author:         Glenn Tesler
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/* graph iterators

 * Glenn Tesler
 * 12/18/02
 */

/* see subgraph.c for definitions */
typedef enum _V_TYPE { V_A_G, V_L_UG, V_B_UG, V_SB_UG, V_B_G,
		       V_SO_G, V_SI_G, V_SOI_G,
		       V_BC_G,
		       V_A_H, V_L_UH, V_B_UH, V_SB_UH, V_B_H,
		       V_SO_H, V_SI_H, V_SOI_H,
		       V_BC_H
} V_TYPE;

/* OR with this to include only one of v, v^c */
#define V_S 0x20

typedef enum _E_TYPE { E_OUT_G=0x0, E_OUT_H=0x1, E_OUT_P=0x2,
		       E_IN_G=0x4, E_IN_H=0x5, E_IN_P=0x6,
		       E_U_G=0x8, E_U_H=0x9, E_U_P=0xA
} E_TYPE;
#define E_G E_OUT_G
#define E_H E_OUT_H
#define E_P E_OUT_P

#define E_DIR_MASK 0xC
#define E_TYPE_MASK 0x3

/* XOR with this to reverse direction */
#define E_REV_MASK 0x4

/* OR with this to only include one of e, e^c */
#define E_S 0x10



typedef enum _C_TYPE { C_A_H,
		       C_S=0x20
} C_TYPE;


/* graph that can be iterated over */
typedef struct igraph {
  int num_nodes;              /* number of nodes in G */
  NODES **nodes;              /* array of nodes in G */

  int visit_id;               /* current timestamp for visiting
			       * sparse subsets of nodes.
			       * Coordinate it with
			       * visit_id, visit_value fields in NODES,
			       * and with MAX_VISIT_ID.
			       */
  int visit_id_next;

} IGRAPH;

/* iterator: all vertices of some type in graph */
typedef struct ivert {
  IGRAPH *G;                  /* iterable graph */
  int v_num;                  /* current vertex number in nodes array */
  int v_type;                 /* type of vertex to iterate over */
  int sym;                    /* 1: only include one of v, v^c */

  int pass;                   /* for V_BC_*:
			       * 1=branching vertices
			       * 2=cycle vertices
			       * else
			       * 0
			       */
} IVERT;

/* iterator: all edges of some type in graph */
typedef struct iedge {
  IGRAPH *G;                  /* iterable graph */
  int e_type;                 /* edge types */
  int sym;                    /* 1: only include one of e, e^c */

  IVERT *vert_it;             /* vertex iterator */
  NODES *v;                   /* current vertex (NULL before first one
			       * or after last one) */
  int e_num;                  /* current edge index */
  int e_max;                  /* one past max edge index for this vertex */

  /* when there are d_in, d_out incoming and outgoing edges, the
   * edge indices are:
   * 0,1,...,d_out-1: outgoing edges
   * -d_in,...,-1: incoming edges
   */

} IEDGE;

/* iterator: all edges of some type incident with a specific vertex in graph */
typedef struct ivedge {
  int e_type;                 /* edge types */
  int sym;                    /* 1: only include one of e, e^c */

  NODES *v;                   /* vertex */

  /* see IEDGE for definition of indices */
  int e_num;                  /* current edge index */
  int e_max;                  /* one past max edge index for this vertex */

} IVEDGE;



/* iterator: all components in graph */
typedef struct icomp {
  IGRAPH *G;                  /* iterable graph */
  int c_type;                 /* component types */
  int sym;                    /* 1: only include one of comp, comp^c */

  IVERT *vert_it;             /* vertex iterator */
} ICOMP;

typedef NODES COMP;






IVERT *it_v_renew(IVERT *it,
		  IGRAPH *G, V_TYPE v_type);
IVERT *it_v_new(IGRAPH *G, V_TYPE v_type);
void it_v_destroy(IVERT *it);
NODES *it_v_next(IVERT *it);
int it_v_count(IGRAPH *G, V_TYPE v_type);

IEDGE *it_e_renew(IEDGE *it,
		  IGRAPH *G, E_TYPE e_type);
IEDGE *it_e_new(IGRAPH *G, E_TYPE e_type);
void it_e_destroy(IEDGE *it);
EDGE *it_e_next(IEDGE *it);
int it_e_count(IGRAPH *G, E_TYPE e_type);

IVEDGE *it_ev_renew(IVEDGE *it,
		    NODES *v, E_TYPE e_type);
IVEDGE *it_ev_new(NODES *v, E_TYPE e_type);
void it_ev_destroy(IVEDGE *it);
EDGE *it_ev_next(IVEDGE *it);
int it_ev_count(NODES *v, E_TYPE e_type);

ICOMP *it_c_renew(ICOMP *it,
		  IGRAPH *G, C_TYPE c_type);
ICOMP *it_c_new(IGRAPH *G, C_TYPE c_type);
void it_c_destroy(ICOMP *it);
int it_c_count(IGRAPH *G, C_TYPE c_type);
COMP *it_c_next(ICOMP *it);








int gdegree(NODES *v, E_TYPE e_type);
int gdegree_cmp(NODES *v, E_TYPE e_type, int m);

void create_full_subgraph(IGRAPH *G, int edge_subg_flag, int init_flag);

COMP *comp_id(NODES *v);
void comp_id_set(NODES *v, NODES *comp);
void subgraph_ins_edge(EDGE *e);
void subgraph_del_edge(EDGE *e);
void subgraph_del_node(NODES *v);

void reset_visit_id(IGRAPH *G);
int get_visit_id(IGRAPH *G, int r);


int mark_deleted_nodes(IGRAPH *G);



NODES *label_distances(IGRAPH *G,
		       NODES *v1, NODES *v2,
		       int L,
		       EDGE *e,
		       int e_type);
int check_path_length(IGRAPH *G,
		      NODES *v1, NODES *v2,
		      int L,
		      int e_type);
NODES *explore_neighborhood(IGRAPH *G,
			    NODES *v,
			    NODES *v_stop,
			    EDGE *prev_edge,
			    int L,
			    int cur_dist,
			    int e_type);



int is_path_internal(NODES *v, int e_type, int L_c);
EDGE *get_first_edge(NODES *v, int e_type);
EDGE *get_unique_edge(NODES *v, int e_type);

void subgraph_del_read_with_edge(READTABLE *RT, EDGE *e);

void debug_print_vval(IGRAPH *G, READTABLE *RT);
void debug_check_sym(IGRAPH *G);
void debug_print_degrees(IGRAPH *G);

void free_igraph(IGRAPH *G);

void clear_visit_flag(IGRAPH *G, V_TYPE v_type);
void clear_visit_edge(IGRAPH *G, E_TYPE e_type);

EDGE *find_unique_oedge(NODES *v);

NODES *augment_graph_ext(IGRAPH *G);

void unaugment_graph_ext(IGRAPH *G);

