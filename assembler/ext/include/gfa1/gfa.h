#ifndef __GFA_H__
#define __GFA_H__

#include <stdio.h>
#include <stdint.h>

#define GFA_VERSION "r53"

/*
  A segment is a sequence. A vertex is one side of a segment. In the code,
  segment_id is an integer, and vertex_id=segment_id<<1|orientation. The
  convention is to use variable u, v or w for a vertex, not for a segment. An
  arc is a directed edge between two vertices in the graph. Each arc has a
  complement arc. A link represents an arc and its complement. The following
  diagram shows an arc v->w, and the lengths used in the gfa_arc_t struct:

       |<--- lv --->|<-- ov -->|
    v: ------------------------>
                    ||overlap|||
                 w: -------------------------->
                    |<-- ow -->|<---- lw ---->|

  The graph topology is solely represented by an array of gfa_arc_t objects
  (see gfa_t::arc[]), where both an arc and its complement are present. The
  array is sorted by gfa_arc_t::v_lv and indexed by gfa_t::idx[] most of time.
  gfa_arc_a(g, v), of size gfa_arc_n(g, v), gives the array of arcs that leaves
  a vertex v in the graph g.
*/

typedef struct {
	uint64_t v_lv; // higher 32 bits: vertex_id; lower 32 bits: lv; packed together for sorting
	uint32_t w, lw;
	int32_t ov, ow;
	uint64_t link_id:62, del:1, comp:1;
} gfa_arc_t;

#define gfa_arc_head(a) ((uint32_t)((a).v_lv>>32))
#define gfa_arc_tail(a) ((a).w)
#define gfa_arc_len(a) ((uint32_t)(a).v_lv) // different from the original string graph

#define gfa_arc_n(g, v) ((uint32_t)(g)->idx[(v)])
#define gfa_arc_a(g, v) (&(g)->arc[(g)->idx[(v)]>>32])

typedef struct {
	uint32_t m_aux, l_aux;
	uint8_t *aux;
} gfa_aux_t;

typedef struct {
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t len2, dummy; // len_r: the other length of the unitig
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	char **name;
} gfa_utg_t;

typedef struct {
	int32_t len;
	uint32_t del:16, circ:16;
	char *name, *seq;
	gfa_aux_t aux;
	gfa_utg_t utg;
} gfa_seg_t;

#define gfa_n_vtx(g) ((g)->n_seg << 1)

typedef struct gfa_s {
	// segments
	uint32_t m_seg, n_seg;
	gfa_seg_t *seg;
	void *h_names;
	// links
	uint64_t m_arc, n_arc:62, is_srt:1, is_symm:1;
	gfa_arc_t *arc;
	gfa_aux_t *arc_aux;
	uint64_t *idx;
} gfa_t;

extern int gfa_verbose;

#ifdef __cplusplus
extern "C" {
#endif

gfa_t *gfa_init(void);
int32_t gfa_add_seg(gfa_t *g, const char *name);
void gfa_destroy(gfa_t *g);

gfa_t *gfa_read(const char *fn);

void gfa_print(const gfa_t *g, FILE *fp, int M_only);

void gfa_symm(gfa_t *g); // delete multiple edges and restore skew-symmetry
void gfa_cleanup(gfa_t *g); // permanently delete arcs marked as deleted, sort and then index
int gfa_arc_del_trans(gfa_t *g, int fuzz); // transitive reduction
int gfa_arc_del_short(gfa_t *g, float drop_ratio); // delete short arcs
int gfa_cut_tip(gfa_t *g, int max_ext); // cut tips
int gfa_cut_internal(gfa_t *g, int max_ext); // drop internal segments
int gfa_cut_biloop(gfa_t *g, int max_ext); // Hmm... I forgot... Some type of weird local topology
int gfa_pop_bubble(gfa_t *g, int max_dist); // bubble popping
gfa_t *gfa_ug_gen(const gfa_t *g);

uint8_t *gfa_aux_get(int l_data, const uint8_t *data, const char tag[2]);

int32_t gfa_name2id(const gfa_t *g, const char *name);
void gfa_sub(gfa_t *g, int n, char *const* seg, int step);

#ifdef __cplusplus
}
#endif

static inline void gfa_arc_del(gfa_t *g, uint32_t v, uint32_t w, int del)
{
	uint32_t i, nv = gfa_arc_n(g, v);
	gfa_arc_t *av = gfa_arc_a(g, v);
	for (i = 0; i < nv; ++i)
		if (av[i].w == w) av[i].del = !!del;
}

static inline void gfa_seg_del(gfa_t *g, uint32_t s)
{
	uint32_t k;
	g->seg[s].del = 1;
	for (k = 0; k < 2; ++k) {
		uint32_t i, v = s<<1 | k;
		uint32_t nv = gfa_arc_n(g, v);
		gfa_arc_t *av = gfa_arc_a(g, v);
		for (i = 0; i < nv; ++i) {
			av[i].del = 1;
			gfa_arc_del(g, av[i].w^1, v^1, 1);
		}
	}
}

#endif
