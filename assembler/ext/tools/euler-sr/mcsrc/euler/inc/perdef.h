/***************************************************************************
 * Title:          predef.h
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/*#define pow1(x) ((x) * (x))*/
#define rev(x) ((2 + (x)) % 4)
#define reverse_read(x, y) ((x) >= (y) ? (x - y) : (x + y))
#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) < (y) ? (x) : (y))
#define numc(i, j) ((i) < (j)? (j)*((j)-1)/2+(i): (i)*((i)-1)/2+(j))
#define numl(i, j) ((i) < (j)? (j)*((j)+1)/2+(i): (i)*((i)+1)/2+(j))
#define CODE_GAP 21

#define LINE_LENGTH 60

/* linked list of readintervals used for sorting */
typedef struct index {
	int index;
	struct index *next;
} INDEX;

/* linked list of readintervals of (read #, position in read)
 * used to identify all nodes glued together into a supernode
 */
typedef struct readposition {
	int readindex, position;
	struct readposition *next;
} READPOSITION;

typedef struct POSPOINT {
	READPOSITION *readposition;
} POSPOINT;


//#define MAX_VISIT_ID   ((1L<<20) - 1)
#define MAX_VISIT_ID   ((1L<<28) - 1)

/* supernode */
typedef struct nodes {
	READPOSITION *readposition;
	int	npos;			/* # (read,pos) at this node */
	int	nlinks;			/* # pw overlaps in this supernode */
	struct nodes	*bal_node;	/* complementary node */
	int	num_path;		/* # paths through node */
	int	*path_index;		/* indices of paths */
	int	*path_pos;		/* this node's offset on each path */
	int	num_nextedge;		/* # outgoing edges */
	int	num_lastedge;		/* # incoming edges */
	struct	edge **nextedge;	/* pointers to outgoing edges */
	struct	edge **lastedge;	/* pointers to incoming edges */
  int index;
	/* subgraph */
	struct nodes    *comp;          /* component id (using set union
                                         * structure) */

#if 0
	unsigned int    visit_id  : 20; /* time stamp for computing values
					 * on sparse subsets of nodes
					 */
	unsigned int	visit_value: 12; /* value at node
					  * _if_ visit_id is current time
					  */
#endif
	unsigned int	visit_value;     /* value at node
					  * _if_ visit_id is current time
					  */

	unsigned int    visit_id   : 28; /* time stamp for computing values
					  * on sparse subsets of nodes
					  */


	/* subgraph flags */
	unsigned int	visit		: 1;	/* flag */
	unsigned int	subg_flag	: 1;	/* is vertex in subgraph? */
	unsigned int	ext_flag	: 1;	/* 0: internal; 1: external */
	unsigned int	comp_ext_flag	: 1;	/* component int/ext flag
						 * only valid @ comp root
						 * 0: internal; 1: external */
	int size;
} NODES;

typedef struct vertex {
	int	num_path;		/* # paths through node */
	int	*path_index;		/* indices of paths */
	int	*path_pos;		  /* this node's offset on each path */
	int	num_nextedge;		/* # outgoing edges */
	int	num_lastedge;		/* # incoming edges */
	struct	edge **nextedge;	/* pointers to outgoing edges */
	struct	edge **lastedge;	/* pointers to incoming edges */
} VERTEX;

/* linked list of read overlaps for Phrap overlap detection	*/
typedef struct readoverlap {
	int	read1, read2;		/* read indices	*/
	int	pos1[2], pos2[2];	/* read positions of overlap parts*/
	int	LLR_score;		/* phrap LLR-score */
	struct readoverlap *next;	/* link to next read in readinterval */
} READOVERLAP;

/* array or linked list of readintervals of all reads that intersect a given edge */
typedef struct readinterval {
	int	begin;			/* start pos in the read */
	int	eq_read;		/* read index */
	int	length;			/* length of intersection */
	int	offset;			/* start pos in the edge */
	char	cov;			/* ? */
	struct readinterval *next;	/* link to next read in readinterval */
	int pos;
} READINTERVAL;


/* edges in any kind of graph */
typedef struct edge {
  int index;
	NODES	*begin, *end;		/* directed edge, verts */

	int	length;			/* length in nucleotides,
					 * incl. both vertices */
	char	*seq;			/* reserved: bases on the edge */

	int	start_cover;		/* average coverage of edge */
	int	multip;			/* # reads intersecting this edge */
	READINTERVAL	*readinterval;	/* array of readintervals
					 * intersecting this edge;
					 * length=multip
					 */

	struct edge	*bal_edge;	/* pointer to complementary edge */

	/* subgraph */
	unsigned int	visit		: 1;	/* flag */
	unsigned int	subg_flag	: 2;	/* is edge in subgraph? */
	unsigned int	ext_flag	: 1;	/* 0: internal; 1: external */
	int deleted;
} EDGE;

#define SUBG_DEL      0
#define SUBG_IN       1
#define SUBG_POSTPONE 2


/* reserved */
typedef struct eqreadinterval {
	READINTERVAL	*readinterval1;
	READINTERVAL	*readinterval2;
	struct readinterval *next;
} EQREADINTERVAL;

/* Pairwise overlap alignment.
 * For n reads + n complements, form a sparse 2n x 2n matrix:
 *    ALIGN **eq_class;
 *
 *    For i=0,1,...,2n-1:
 *    eq_class[i] = double-linked list of alignments with read i.
 *                  each member of list has reads[0] = i.
 *                  reads[1] is the other read it's linked with.
 *                  Convention: reads[0] < reads[1].
 *    Only need to put in an alignment in one direction;
 *    if reads i & j are aligned,
 *    the complements (2n-1-j, 2n-1-i) will be inferred.
 */

typedef struct align {
	int	reads[2];		/* indices of reads */
	int	*pos[2];		/* indel @ read 0, position pos[0,i]
					 * &       read 1, position pos[1,i]
                                         * for i=1,...,length-2
                                         * and start is i=0, end is i=length-1
					 */
	int	length;			/* # marks (start + end + indels) */
	int	mis_match;		/* # mismatches
					 *   + # indels (no start/end) */
	char	cov;			/* flag: remove the alignment or
					 * not? */
	struct align *next;		/* doubly-linked list of
					 * read intervals of all alignments
					 * with a read, s.t. the read has
					 * the smaller index */
	struct align *prev;
} ALIGN;

typedef struct linkpos {
	int readindex, position;
	struct linkpos *next;
} LINKPOS;


/* used for locating common 20-mers */
typedef struct hash {
	LINKPOS    *linkpos;
	struct hash *next;
} HASH;

typedef struct table {
	HASH    *prev;
} TABLE;

/* reserved for interval graph */
typedef struct block {
	struct block *bal_block;
	struct block *prev;
	struct block *next;
	int	multiplicity;
	READINTERVAL	*readinterval;
} BLOCK;

/* reserved for interval graph */
typedef struct bintree {
	struct bintree *left, *right;
	int	readindex;
	int	minpos, maxpos;
} BINTREE;

/* no longer used */
typedef struct insert {
	NODES	**node;
	int	num_nodes;
	struct insert *next;
} INSERT;

/* for each read, a list of pointers to consecutive supernodes */
typedef struct readlist {
	NODES	*node;
	INSERT	*insert;
} READLIST;

/* path through graph */
typedef struct path {
	EDGE  **edge;			/* list of pointers to edges */
	int len_path;			/* # of edges */
	int seqLength; /* # of nucleotides.*/
        /* The start position on 1st edge is begin_length from its end */
	int begin_length;

        /* The terminal position on final edge is end_length from its start */
	int end_length;
	/* The read index of the corresponding read paths */
	/* 0 -- means mate-paths; >0 means read paths; x-1 is the read index; */
	int readindex;
	/* 0 -- means read-paths; >0 means mate paths; x-1 is the read index for both [0] and [1]; */
	int pairindex[2];
	int index;
} PATH;


/* table of the reads or other strings */
typedef struct readtable {
  /* reads */
  char **src_seq;             /* symbols in raw reads */
  int  num_seq;               /* number of raw reads, not counting
			       * complements;
			       * reads 0 .. (num_seq-1) are orig reads
			       * reads num_seq .. (2*num_seq-1) are
			       * complements
			       */
  char **src_name;            /* src_name[i] = name of sequence i
			       * i=0..(num_seq-1)
			       */
  int *len_seq;               /* len_seq[i] = length of read i,
			       * i=0..(2*num_seq-1)
			       */
  int *new_len_seq;           /* revised lengths of reads */
  READLIST **readlist;        /* how reads map into original graph */


  /* chimeric reads */
  int *chim;                  /* indices of chimeric reads */
  int num_chim;               /* number of chimeric reads */
  int num_chim_alloc;         /* # entries of space allocated for *chim,
			       * regardless of whether it's all used
			       */
} READTABLE;


/* table of mate-pairs */
typedef struct matepairtable {
  int num_pair;

  /* mate pair i=0,1,...,num_pair-1
   * consists of read pair1[i] and inverse of pair2[i]
   */
  int *pair1;
  int *pair2;

  /* mate pair i segments should have gap between min_dist[i] and max_dist[i]
   */
  int *min_dist;
  int *max_dist;

  /* ntype[i] = type of mate pair i */
  int *ntype;
   /* name of mate-pairs -- haixu */
   char **name;
  /*(Haixu) temporary for EULER-SF */
  int num_sf_pair, *sf_pair;
} MATEPAIRTABLE;


/*****************************************************************************
 * mate-pair rules
 * TODO: this uses fixed size allocations; change to adapt to correct sizes
 *
 * ntyperule[0] = # finishing read types
 * ntyperule[1] = # single-barrelled read types
 * ntypepair    = # double-barrelled read types
 *
 * Finishing reads:                f       ALL     ALL
 * namerule[0][i] = 'f'
 * "ALL"'s not used
 * i=0,1,..., ntyperule[0]-1
 *
 * Single reads:                   s       ALL     ALL
 * namerule[1][i] = 's'
 * "ALL"'s not used
 * i=0,1,..., ntyperule[1]-1
 *
 * Double-barreled reads:          x y     ALL     ALL     1500 3500
 * pairrule[i][0] = 'x' (only uses 1st char)
 * pairrule[i][1] = 'y' (only uses 1st char)
 * librule[i] = "ALL"   (1st ALL)
 * platerule[i][0] = 0            (2nd ALL)
 * platerule[i][1] = 100000       (2nd ALL)
 * pairrange[i][0] = 1500
 * pairrange[i][1] = 3500
 * i=0,1,...,ntypepair-1
 *****************************************************************************/

typedef struct matepairrules {
  char librule[500][100];
  char namerule[2][10];
  char pairrule[500][2];
  int ntyperule[2], ntypepair;
  int pairrange[500][2];
  int platerule[500][2];
  /*   int matetype; */
} MATEPAIRRULES; 

/* table of mate-pair rules, revised
 * INCOMPLETE
 */
typedef struct matepairrules0 {
  int ntype_fin;     /* number of rules for finishing reads */
  int ntype_sing;    /* number of rules for single reads */
  int ntype_pair;    /* number of rules for mate-pairs */

  char **librule;
  char *namerule[2];
  char *pairrule[2];
  int *pairrange[2];
  int *platerule[2];
} MATEPAIRRULES0;


typedef struct stats {
  int min_value;
  int max_value;
  int total;        /* sum of values */
  int total2;       /* sum of squares of values */
  int n;            /* number of data points; this is same as link_cov[i] */
} STATS;

/* table of scaffolding pairs
 * pair i=0,1,...,num_sf
 * is from edge linkedge1[i] to edge linkedge2[i]
 * with coverage link_cov[i]
 *
 * statistics on range [lo,hi]
 *   range_lo[i] = statistics on lo
 *   range_hi[i] = statistics on hi
 */

typedef struct sfpairs {
  EDGE **linkedge1;
  EDGE **linkedge2;
  int *link_cov;

  STATS *range_lo;
  STATS *range_hi;

  int num_sf;

  /* Used for euler-db only 
   * for recording mate-pairs that have
   * been mapped on the single edges;
   * after euler-db transformation,
   * don't map them on the transformed graph again.
   */
  EDGE **all_edge;
  int tot_edge;
} SFPAIRS;
