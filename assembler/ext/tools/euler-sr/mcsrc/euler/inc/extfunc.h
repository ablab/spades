/***************************************************************************
 * Title:          extfunc.h
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <set>
#if !defined( min ) 
#define min(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max
#define max(a,b) ((a) < (b) ? b : a)
#endif

extern ALIGN *AlignReads(char **src_seq, char **score, int *len_seq, int num_seq, READINTERVAL *readinterval, int index);
extern ALIGN *AlignReadsPO(char **src_seq, int *len_seq, int num_seq, READOVERLAP *readoverlap);
extern int overalign(char *seq1, char *seq2, char *score1, char *score2, 
	      int len_seq1, int len_seq2, int offset, int *sapp, int *cutp, int *seql,
	      int index1, int index2, int num_seq);
extern double perc_qual(char *score, int s1, int s2);
extern int cal_identity(char *seq1, int len1, char *seq2, int len2, int *sapp);
extern void compute_score(char *seq1, int len1, char *seq2, int len2, int *sapp, double *score);
extern ALIGN *new_align(int *cutp, int *seql, int *sapp, ALIGN *align, int r1, int r2, int offset, int mis_match);
extern void buildblock(int num_seq, int *len_seq, EQREADINTERVAL **eq_class, READLIST **blocklist, BLOCK *block);
extern char char2intgen(char c);
extern char char2int(char c);
extern void *ckalloc(int amount);
extern void *ckrealloc(void *p, size_t new_size, size_t old_size);
extern FILE *ckopen(char *name, char *mode);
extern int alignlength(ALIGN *align);
extern ALIGN *remove_trans(ALIGN *align1, ALIGN *align2, ALIGN *align3);
extern ALIGN *chkconn(int i1, int i2, ALIGN **eq_class, int num_seq);
extern char chkedgelink(EDGE *edge);
extern char chksingle(NODES *node, EDGE *edge, EDGE *sedge, int offset, int cumoff);
extern int count_edge(NODES **vertex, int num_vertex, int **num_pa);
extern int count_edge_simp(NODES **vertex, int num_vertex, int **num_pa);
extern int forwardtanglelen(EDGE *edge, int *ave, int *multip);
extern int backtanglelen(EDGE *edge, int *ave, int *multip);
extern int count_tangle(NODES **vertex, int num_vertex, int disttangle[][7]);
extern int count_bal(NODES *vertex);
extern void getmaxedge(EDGE *edge, EDGE **maxedge);
extern int count_vertex(EDGE **edge, int num_edge, NODES **vertex);
extern int collect_vertex(NODES *v, NODES **vertex, int num_vertex);
extern EDGE *find_bal_edge(EDGE *edge, int *len_seq, int num_seq, int index);
extern char chk_readinterval(READINTERVAL *readinterval1, READINTERVAL *readinterval2, int n1, int n2, int *len_seq, int num_seq);
extern int singstatnode(NODES *node, int **nstat);
extern void statnode(READLIST **list, int *len_seq, char **src_seq, int num_seq);
extern int countnode(READLIST **list, int *len_seq, int num_seq);
extern int cleannode(READLIST **list, int *len_seq, int num_seq);
extern char del_trans(EDGE *edge1, EDGE *edge2);
extern long trans_seq(char *seq, int len);
extern LINKPOS *ins_linkpos(LINKPOS *linkpos, int i1, int i2);
extern HASH *free_hash(HASH *hash);
extern HASH *ins_hash(HASH *hash, int i1, int i2, char **src_seq, int num_seq);
extern READINTERVAL *ins_readinterval(READINTERVAL *readinterval, int readindex, int offset);
extern int detectoverlap(char **src_seq, char **score, int *len_seq, int num_seq, READINTERVAL **readinterval);
extern int comp_word(char *w1, char *w2, int len);
extern READINTERVAL *new_readinterval(HASH *hash, int i1, int i2, char **src_seq, READINTERVAL *readinterval, int num_seq);
extern void erasedge(EDGE *edge);
extern void erasenext(NODES *vertex, int n);
extern void eraselast(NODES *vertex, int n);
extern int searcherase(EDGE **edge, EDGE *e, int num);
extern void errcorrt(READLIST **list, int *len_seq, char **src_seq, char **score, char **src_name,
	      int num_seq, FILE *fp, FILE *fp1);
extern READINTERVAL *free_readinterval(READINTERVAL *readinterval);
extern READOVERLAP *free_readoverlap(READOVERLAP *readoverlap);
extern int size_readinterval(READINTERVAL *readinterval);
extern void free_graph(NODES **vertex, int num_vertex);
extern int GALIGN(char *A, char *B, char *SA, char *SB, int M, int N, int low, int up, int W[][15],
	  int G, int H, int *S, int MW, int MX);
extern INDEX *insert_index(INDEX *index, int p);
extern INDEX *free_index(INDEX *index);
extern void initialize(READLIST **list, int *len_seq, int num_seq);
extern EDGE *insert_edge(NODES *node, NODES *node_next, int read, int pos1, int pos2);
extern void add_nextedge(NODES *node, EDGE *edge);
extern void add_lastedge(NODES *node, EDGE *edge);
extern int LOCAL_ALIGN(char *A, char *B, char *SA, char *SB, int M, int N, int low, int up, int W[][15],
		int G, int H, int *psi, int *psj, int *pei, int *pej, int MW);
extern int makedge(NODES *node, EDGE **edge, int num_edge, READLIST **list);
extern EDGE *newedge(EDGE **midedge, int num_midedge, READLIST **list);
extern int makedge_new(NODES *node, EDGE **edge, int num_edge, EDGE *edge1, NODES **allnodes);
extern EDGE *newedge_new(EDGE **midedge, int num_midedge, NODES *begin, NODES *end);
extern void sortreadinterval_index(READINTERVAL *readinterval, int multip);
extern int readintervalcompar_index(READINTERVAL *a, READINTERVAL *b);
extern int copyreadinterval(READINTERVAL *readinterval, READINTERVAL *readcov);
extern void sortreadinterval(READINTERVAL *readinterval, int multip);
extern int readintervalcompar(READINTERVAL *a, READINTERVAL *b);
extern void sort_nodepos(NODES *node);
extern int poscompar(POSPOINT *a, POSPOINT *b);
extern char findnode(INSERT *insert, NODES **node0, int num, int read, int pos);
extern void ins_node(READLIST **list, int read1, int read2, int pos1, int pos2, int endpos);
extern INSERT *insert_nodes(INSERT *insert, NODES **nodes, int num, int read, int pos);
extern INSERT *free_insert(INSERT *insert);
extern NODES *chk_merge_node(READLIST **list, int read1, int read2, int pos1, int pos2);
extern void insert_position(NODES *node, int read, int pos);
extern void merge(int num_seq, int *len_seq, ALIGN **eq_class, int num_readinterval, READLIST **list);
extern NODES *combine_nodes(NODES *node1, NODES *node2);
extern NODES *free_nodes(NODES *node);
extern void update_link(NODES *node, int r1, int s1, int r2, int s2);
extern void output_graph(NODES **vertex, int num_vertex, FILE *fp);
extern int  ran_number(int n, int *idum);
extern double random1(int *idum);
extern int readreadinterval(ALIGN **eq_class, int num_seq, FILE *fp);
extern int readoverlap(READINTERVAL **readinterval, FILE *fp);
extern void allocpath(PATH *path, int len);
extern void findpath(PATH *path, EDGE *edge, int *len_seq, int num_seq, char *label);
extern void backpath(PATH *path, NODES *vertex, int reads, int end, int *len_seq);
extern void forpath(PATH *path, NODES *vertex, int reads, int begin, int *len_seq);
extern int readpath(NODES **vertex, int *num, PATH *path, int *len_seq, int num_seq, int *chim, int *num_chim);
extern int readqual(char **val, char **src_seq, int *len, FILE *fp);
extern int readscore(char **val, int *len, FILE *fp);
extern int readseq1by1(char **src_seq, char **src_name, int *len_seq, FILE *fp);
extern int readseq1by1gen(char **src_seq, char **src_name, int *len_seq, FILE *fp);
extern void read_fasta_file(char *seqfile, READTABLE *RT);
extern void free_readtable(READTABLE *RT);
extern READINTERVAL *insert_unitig(READINTERVAL *unitig, NODES *node, EDGE *edge, int *len_seq);
extern int search_unitig(NODES **nodeall, NODES *node, READINTERVAL **unitig, int nunitig, int num_seq, int *len_seq);
extern int merge_graph(NODES **vertex, int num_vertex, PATH *paths, int numPaths);
extern int shave_graph(NODES **vertex, int num_vertex, int *chim, int *num_chim);
extern EDGE *merge_vertex(NODES *vertex);
extern int chklist(int *list, int n, int index);
extern EDGE *new_edge(NODES *vertex, NODES *begin, NODES *end, EDGE *edge1, EDGE *edge2, int *beginlist, int *endlist,	int nbegin, int nend);
extern void combine_readinterval(EDGE *edge1, EDGE *edge2, EDGE *newedge, int *beginlist, int *endlist, int nbegin, int nend);
extern int findposition(READINTERVAL *readinterval, int num, int readindex, int position);
extern int FindIntervalBeginningAtPos(READINTERVAL *readinterval, int num, int readindex, int readposition);
extern int FindIntervalEndingAtPos(READINTERVAL *readinterval, int num, int readindex, int readposition);
extern void visit_comp(NODES **node, NODES *begin_node, int *comp, int num_comp, int num_seq);
extern int countwidth(NODES *node);
extern int countthickness(NODES *node);
extern void countallnode(READLIST **list, int *len_seq, int num_seq);
extern ALIGN *free_align(ALIGN *align);
extern int size_align(ALIGN *align);
extern char chkedge(NODES *node1, NODES *node2);
extern EDGE *merge_vertex_path(NODES *vertex, PATH *path, int num_path);
extern int merge_graph_path(NODES **vertex, int num_vertex, PATH *path, int num_path, int small_cover);
extern int eqtrans_bal(NODES ***vertex, int num_vertex, PATH *path, int num_path, int num_seq, int cut);
extern int rm_edge(NODES *vertex, PATH *path, int num_path);
extern void replace1edge(PATH *path, int num_path, EDGE *edge1, EDGE *edge2);
extern void derivelist(PATH *path, int num_path, NODES *vertex, EDGE *edge1, EDGE *edge2,
		int *edgematch1, int *edgematch2, int *beginlist, int *endlist, int *nbegin, int *nend);
extern void reducepath(PATH *path, int num_path, NODES *vertex, EDGE *edge1, EDGE *edge2, EDGE *newedge, 
		int *edgematch1, int *edgematch2);
extern int numedge(EDGE *edge, PATH *path);
extern EDGE *detach_bal(EDGE *edge1, EDGE *edge2, PATH *path, int num_path, int *f_edge, int num_seq);
extern int searchlast(NODES *vertex, EDGE *edge, PATH *path, int num_path, int *match);
extern int searchnext(NODES *vertex, EDGE *edge, PATH *path, int num_path, int *match);
extern int countmatch(EDGE *edge1, EDGE *edge2, PATH *path, int num_path);
extern int countstartmatch(EDGE *edge, NODES *vertex, PATH *path, int num_path);
extern int countendmatch(EDGE *edge, NODES *vertex, PATH *path, int num_path);
extern int collect2forpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midforpath,
		     PATH *path, int num);
extern int collect2aftpaths(NODES *vertex, EDGE *edge1, EDGE *edge2, PATH *midaftpath,
		     PATH *path, int num);
extern int collectstartpaths(NODES *vertex, EDGE *edge, PATH *startpath, PATH *path, int num);
extern int collectendpaths(NODES *vertex, EDGE *edge, PATH *endpath, PATH *path, int num);
extern void resetpath(PATH *path, int num_path, NODES *vertex, EDGE *lastedge, EDGE *nextedge, EDGE *newedge, int begtag, int endtag);
extern int count_path(PATH *path, int num_path);
extern int cutbegpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *nextedge);
extern int cutendpath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *lastedge);
extern void remove_edge(PATH *path, int path_index, int path_pos);
extern int chk_path(int p, int *vp, int nvp);
extern int splitbeg(NODES ***vertex, int num_vertex, PATH *path, int num_path);
extern NODES *new_vertex(NODES *vertex);
extern void move_path_next(EDGE *nextedge, NODES *vertex, NODES *vertex0, PATH *path);
extern void move_path_last(EDGE *lastedge, NODES *vertex, NODES *vertex0, PATH *path);
extern int chk_consist(PATH *startpath, PATH *midpath, int num_midpath, int *n);
extern void add_path(NODES *vertex, int path_index, int path_pos);
extern void rem_path(NODES *vertex, int path_index, int path_pos);
extern void set_path(NODES **vertex, int num_vertex, PATH *path, int num_path);
extern void merge_gap(char **src_seq, int num_seq, int *len_seq, ALIGN **eq_class, READLIST **list, int gap_k);
extern void errcorrt_pair_mem(ALIGN **eq_class, char **multip, int *len_seq, char **src_seq, char **score, char **src_name,
		   int num_seq, int LOW_COV, int intv, FILE *fp, FILE *fp1);
extern void errcorrt_pair(ALIGN **eq_class, int **multip, int *len_seq, char **src_seq, char **score, char **src_name,
		   int num_seq, int LOW_COV, int intv, FILE *fp, FILE *fp1);
extern void trim_align(ALIGN *align, char **src_seq, int *len_seq);
extern int graph(int num_seq, char **src_name, int *len_seq, READLIST **list, EDGE **edge);
extern void movereadinterval(EDGE *edge1, EDGE *edge);
extern int realmultip1(EDGE *edge, int l);
extern int realmultip2(EDGE *edge, int l);
extern void count_multip(NODES **vertex, int num_vertex);
extern void writeseq(FILE *fp, char *seq, char *name, int length);

extern void initial_edge(NODES **vertex, int num_vertex, char **src_seq, int *len_seq,int num_seq);
extern void initedge(EDGE *edge, int *len_seq, char **src_seq);
extern int readpath_single(NODES *start_node, PATH *path);
extern void output_graph_single(NODES **vertex, int num_vertex, FILE *fp);
extern void singlepath(NODES *start_node, PATH *path, int reads, int begin);
extern READINTERVAL *chk_redund(READINTERVAL *readinterval);
extern void output_align(ALIGN *align, char **src_name, char **src_seq, char **score, int *len_seq, int num_seq);
extern void output_align_ns(ALIGN *align, char **src_name, char **src_seq, int *len_seq, int num_seq);
extern int filter_readinterval(ALIGN **eq_class, READTABLE *rt);
extern double perc_maxnuc(char *seq, int d1, int d2, double *perc);
extern void caln_score(char *seq1, int len1, char *seq2, int len2, double *score, ALIGN *align);
extern int dist_range(int p1, int p2);
extern void update_bound(READINTERVAL *leftreadinterval, READINTERVAL *rightreadinterval, READINTERVAL **readinterval2, int p1, int p2);
extern int comppath(PATH *path1, PATH *path2);
extern int link_edge(int num_seq, int *len_seq, READLIST **list, NODES **nodes);
extern void cleangraph(NODES **nodes, int num_nodes);
extern int branch_graph(int num_seq, int *len_seq, NODES **nodes, EDGE **edge, int num_nodes);
extern int crossedge(EDGE *edge, PATH *path);
extern int newvertex(NODES ***vertex, int num_vertex, EDGE *edge, PATH *path);
extern int ckbeginedge(EDGE *edge, EDGE *edge2, PATH *path);
extern int ckendedge(EDGE *edge, EDGE *edge1, PATH *path);
extern int createvertex(NODES **vertex, int num_vertex, NODES *vertex_now, PATH *path, int num_path);
extern int chkcross(NODES *vertex_now, PATH *path, int num_path);
extern void makenode(READLIST **list, int i, int j, int *len_seq, int num_seq);
extern void separatenode(NODES *node, READLIST **list, int *len_seq, int num_seq);
extern void separatenode0(NODES *node, READLIST **list, int *len_seq, int num_seq, int mthresh);
extern int add_prenode(NODES **prenode, int nnode, int *mnode, NODES *node);
extern void 	splitnode(NODES *node, NODES *prenode, READLIST **list, int *len_seq, int num_seq);
extern void rem_position(NODES *node, int reads, int pos);
extern int searchnextpath(NODES *node1, NODES *node2, READLIST *list, int max_len, int len, char *spath);
extern int fillreadlist(READLIST *list, NODES *node_last, NODES *node, int pos_last, int pos, IGRAPH *G);
extern void mapread(IGRAPH *G, READTABLE *RT);
extern void destroylist(READLIST *newlist, int index, int len);
extern int makebothlist(READLIST *newlist, int index, int len, int num_seq);
extern EDGE *find_unique_oedge(NODES *v2);
extern int newreadpath(NODES **vertex, int num_vertex, PATH *path, int *len_seq, int num_seq);
extern READPOSITION *free_position(READPOSITION *readposition);
extern void findendreads(NODES **vertex, int num_vertex, int *endreads, int *beginreads, int *num_endreads, int num_seq);
extern void findmatereads(NODES **vertex, int num_vertex, int *endreads, int *beginreads, int *num_endreads, int num_seq, MATEPAIRTABLE *MP);
extern int findmate(int k, MATEPAIRTABLE *MP);
extern char chkmate(int n, EDGE *edge, int num_seq);
extern char chk_pair(ALIGN **align, int read1, int read2, int num_seq);
extern int identity(char *seq1, int len1, char *seq2, int len2, int *sapp);
extern int alignend(char **src_name, char **src_seq, char **score, int *len_seq, int num_seq, int *endreads,
	     int *beginreads, int *num_endreads, ALIGN **eq_class);
extern int wholealign(char *seq1, int len_seq1, char *seq2, int len_seq2, int *sapp, int *cutp, int *seql);
extern int trans_pos(char *seq, int len);
extern int chkhash(LINKPOS **hash, int index, char *seq, int len, int word_len, int *pair, int *mount);
extern LINKPOS *inserthash(LINKPOS *hash, int index);
extern int counthash(LINKPOS *hash, int index, int *pair, int *mount, int totp);
extern void buildhash(LINKPOS **hash, int index, char *seq, int len, int word_len);
extern int readblast(ALIGN **align, int *len_seq, int num_seq, FILE *fp, int MIN_LEG, double MIN_ID);
extern void write_gvz_graph(FILE *fp, NODES **vertex, int num_vertex);
extern void write_interval(NODES **vertex, int num_vertex, FILE *fp);
extern int read_interval(NODES **vertex, int num_vertex, FILE *fp);
extern READINTERVAL *insert_readcov(READINTERVAL *readcov, int readindex, int begin, int length, int offset);
extern int erasecycle(NODES **vertex, int num_vertex);
extern void rem_chim(NODES **vertex, int num_vertex, int *chim, int num_chim, int num_seq);
extern int chk_chim2(int *chim, int num_chim, int index, int num_seq);
extern void insert_interval(EDGE *edge, int read, int pos1, int pos2);
extern int filter_path(PATH *path, int num, int num_path);
extern void consensus(NODES **vertex, int num_vertex, char **src_seq, int num_seq, int *len_seq);
extern int revise_consensus(EDGE *edge, char **src_seq, int num_seq);
extern int newsfpair(EDGE *begin_edge, EDGE *end_edge, SFPAIRS *SFP);
extern void insert_chim(READTABLE *RT, int index);

extern
void make_Agraph(READTABLE *RT,
		 ALIGN **eq_class,
		 int split_badnodes,
		 int is_sym,   /* reserved: 1=symmetric graph, 0=not */
		 int use_int,  /* reserved: 1=interval graph, 0=nuc graph */
		 IGRAPH *G
		 );

extern
void read_graph_file(char *edgefile,
		     char *graphfile,
		     int *num_vertex,
		     NODES ***vertex,
		     int *num_edge,
		     EDGE ***edge);

extern
int read_interval_file(char *intvfile,
			int num_vertex,
			NODES **vertex);
extern
void write_interval_file(char *intvfile,
			 int num_vertex,
			 NODES **vertex);

extern
void read_matepair_file(char *pairfile,
			MATEPAIRRULES *MPR,
			READTABLE *RT,
			MATEPAIRTABLE *MP);
extern void free_matepair(MATEPAIRTABLE *MP);


extern void read_matepair_rules_file(char *rulefile, MATEPAIRRULES *MPR);

extern int output_contig_files(char *filestem, int num_vertex, NODES **vertex, READTABLE *RT);

extern void write_sf_gvz_file(char *contig_fname, char *suffix, SFPAIRS *SFP);
extern void free_path(int num_path, PATH *path);

extern int newsfpair(EDGE *begin_edge, EDGE *end_edge, SFPAIRS *SFP);
extern void init_sfpairs(SFPAIRS *SFP, MATEPAIRTABLE *MP, int has_len);
extern void free_sfpairs(SFPAIRS *SFP);
extern int readseqname(char **src_name, FILE *fp);
extern int readseqpar(int *max_leg, int *max_name_leg, FILE *fp);
extern void inputreadname(READTABLE *rt, char *inpfile);
extern void inputread(READTABLE *rt, char *inpfile);
extern char **init_score(READTABLE *RT);
extern void read_score_file(char *qualfile, char **score, READTABLE *RT);
extern int chkforbid(EDGE *edge1, EDGE *edge2, PATH *path);
extern void trimpath(PATH *path, PATH *path_rev, int readindex, int readindex_rev);
extern int findreadinterval(EDGE *edge, int readindex);
extern void remove_readinterval(EDGE *edge, int index);
extern int readphrapovp(READOVERLAP **readoverlap, int threshold, FILE *fp);
extern void print_hl(FILE *fp);
extern void print_header(FILE *fp, char *ptitle);
extern void print_line(FILE *fp, char *line);
extern void print_tailor(FILE *fp);
extern void print_table(FILE *fp, int ncul, int nrow, char ***content, char *caption);
extern char ***allocate_content(int nrow, int ncol, int len);
extern char ***free_content(char ***content, int nrow, int ncol);
extern int defineshift(char *newseq, int len_newseq, char *readseq, int len_readseq, int begin);
extern int defineshift_rev(char *newseq, int len_newseq, char *readseq, int len_readseq, int begin, int len);
extern int findri(READTABLE *RT, READINTERVAL *readinterval, int num_intv, char *seq, int len, char *label);
extern int acculen(char *seq, int len);
extern void chg_intv(READINTERVAL *readinterval, char **newseq, int *len_newseq, READTABLE *RT);
extern int compseq_gap(char *seq1, int len1, char *seq2, int len2);
extern int process_row(int *startpos, char **readseq, int *len_readseq, int num_seq, char *row, int len_row);
extern int addrow(char *readseq, int len1, char *str, int len2);
extern int addseq(char *readseq, char *str, int pos, int length);
extern int nextletter(char *row, int pos, int len);
extern void chg_reads(READTABLE *RT, READINTERVAL *readinterval, int num_intv, int startpos, char *readseq, 
	       int len_readseq, char **newseq, int *len_newseq, int index, char *label);
extern char comp_qual(int n1, int n2);
extern int readalnseq(FILE *fp, char **readseq, int *len_readseq, int *startpos, int num_intv, int length,
		char *contigseq, char *contigqual, int *len_contigseq);

extern void print_text_line(FILE *fp, int length);
extern void print_emptyline(FILE *fp);
extern void print_chimtable(FILE *fp, READTABLE *RT);
extern void print_hl(FILE *fp);
extern void print_header(FILE *fp, char *ptitle);
extern void print_line(FILE *fp, char *line);
extern void print_tailor(FILE *fp);
extern void print_table(FILE *fp, int nrow, int ncol, char ***content, char *caption);
extern void write_gvz_graph(FILE *fp, NODES **vertex, int num_vertex);
extern void print_section_head(FILE *fp, char *line);
extern void write_gvz_file(char *fname, int num_vertex, NODES **vertex, int verbose);
extern void statspath(PATH *path, int num_path);
extern int filter_path_mate(PATH *path, int num, int num_path);
extern int shave_graph(NODES **vertex, int num_vertex, int *chim, int *num_chim);
extern int shave_graph_new(NODES **vertex, int num_vertex, READTABLE *RT, int MIN_LENGTH, int MIN_MULTIP);
extern int locpath(EDGE *edge1, EDGE *edge2, int min_leg, int max_leg, int len, PATH *path, 
									 std::set<EDGE*> &traversedEdges);
extern void ResetTraversedEdges(std::set<EDGE*> &traversedEdges);
extern double random1(int *idum);
extern int vertex_run(NODES **vertex, int num_vertex, PATH *path, int num_path, int small_cover, int num_seq);
extern void write_graph(NODES **vertex, int num_vertex, FILE *fp, FILE *fp1);
void write_branching_graph(NODES **vertex, int num_vertex, FILE *fp);
extern int updateedge(char **edgeprof, EDGE *edge);
extern int LOCAL_ALIGN0(char *A, char *B, int M, int N, int low, int up, int W[][15],
		int G, int H, int *psi, int *psj, int *pei, int *pej, int MW);
extern int ALIGN0(char *A, char *B, int M, int N, int low, int up, int W[][15],
		  int G, int H, int *S, int MW, int MX);
extern void readpar();
extern int filter_path_read(PATH *path, int num, int num_path, int filterMate);
extern void output_contig_cons(NODES **vertex, int num_vertex, char **src_name, char **src_seq,
		 int num_seq, FILE *fp, FILE *fp1, FILE *fp2);
extern int output_contig(NODES **vertex, int num_vertex, 
												 PATH *paths, int num_path,
												 char **src_name, char **src_seq,
												 int num_seq, FILE *fp, FILE *fp1, FILE *fp2, FILE *fp3, 
												 FILE *braphFP, FILE *pathFP);
extern void writeseq(FILE *fp, char *seq, char *name, int length);
extern int output_contig_files(char *filestem, int num_vertex, NODES **vertex,
															 PATH * paths, int numPaths,
															 READTABLE *RT);
extern void set_score();
extern void print_gaps(FILE *fp, int *dist);

extern int readclass(ALIGN **eq_class, int num_seq, FILE *fp);
extern void makegraph(IGRAPH *G, READTABLE *RT, ALIGN **eq_class);
extern void destroygraph(IGRAPH *G);
extern void free_readlist(READTABLE *RT);
extern void write_bgraph_gviz_file(IGRAPH *G, char *filename);
extern int readrules(FILE *fp);
extern int readrules2(FILE *fp, MATEPAIRRULES *MPR);
extern char splitname(char *src_name, int *plate, char *library);
extern char extractname(char *src_name);
extern int loc_pair(char *name, char ext, char **src_name, char *ext_name, int num_seq);
extern char pairrules(char ext_name, int plate, char *library, int *range);
extern void free_label();
int FindRead(READINTERVAL *intervals, int numIntervals, int readIndex, int *intervalPos);
void write_paths(NODES **vertex, int numVertices, 
								 PATH *path, int numPaths, FILE *fp);

extern void RemoveUninformativePaths(PATH *path, int num_path);
extern int GetMaxPathLength(PATH *paths, int numPaths);
extern void TrimShortEnds(PATH *paths, int numPaths, int minPathLength, int maxEdgeLength);
extern void UpdatePathIndices(NODES *vertex, int path, int offset);
extern void replacepath(PATH *path, int num_path, NODES *vertex, EDGE *edge, EDGE *newedge);
extern int RemoveIndexedReadinterval(EDGE *edge, int pos);
extern int  LookupReadintervalIndex(EDGE *edge, int index, int pos);
extern void RemovePathsFromVertex(NODES *vertex, int *pathIndices, int *pathPos, int numPaths);
extern void EraseEdgeAndPassingPaths(EDGE *edge, path *paths, int numPaths);
extern int RemovePathIntervalFromVertex(NODES *vertex, int path, int pos);
extern void CountCoverage(NODES **vertices, int numVertices, int *totalLengthParam, int *totalReadsParam);
extern void CheckEdgeMultiplicities(NODES **vertices, int numVertices, int totalLength, int totalReadStarts);
extern int RemoveLowCoverageEdges(NODES **vertices, int numVertices,
																	PATH *paths,int numPaths,
																	int totalLength, int totalReadStarts, 
																	float maxStddev, int minCoverage);
extern int EdgeConnectsPaths(EDGE *edge);

extern void MarkEdgeForRemoval(EDGE *edge);
extern void RemovePathsWithRemovedEdges(PATH *paths, int numPaths );
extern void EraseEdgesMarkedForDeletion(NODES **vertices, int numVertices);
extern void RemovePathFromVertex(NODES *vertex, int path);
extern void EraseMarkedEdgesAndPaths(NODES **vertices, int numVertices,
																		 PATH *paths, int numPaths);
extern int CountShortDirectedCycles(NODES **vertices, int numVertices, int maxCycleSize);
extern int StraightenShortDirectedCycle(NODES ***verticesP, int numVertices,
																				PATH *paths,
																				EDGE *enter, EDGE *repeat, EDGE *link, EDGE *exitEdge,
																				EDGE **repeatCopy);
extern int FindIndexBoundByEdge(PATH *paths, int pathIndex, int pathPos,
													EDGE *exit);
extern int FindReadPos(PATH *paths, int pathIndex, int pathPos);
extern int FindIntervalIndex(EDGE *edge, int readIndex, int readPos);
extern void RemovePathsBetweenEdges(PATH *paths,
														 int  *pathIndices, int *pathPositions,
														 int numPathIndices, 
														 EDGE *enter, EDGE *exit,
														 READINTERVAL ***removedIntervals,
														 int **numPathIntervals);

extern int CountRepeatCopy(PATH *paths, int *pathIndices, int *pathPositions, int numIndices,
										EDGE *enterEdge, EDGE *repeat, EDGE *exitEdge);

extern int RemoveDuplicatePathIndices(int **pathIndices, int numIndices);
extern int CountPathsStartingAtVertex(PATH *paths,
															 VERTEX *firstVertex, EDGE *firstEdge, 
															 int *pathIndices);
extern int StorePathsStartingAtVertex(PATH *paths,
															 VERTEX *firstVertex, EDGE *firstEdge, int **pathIndices);
extern int CountPathsThroughVertex(PATH *paths, EDGE *enter, NODES *vertex, EDGE *exit,
														int *pathIndices, int *pathPositions);
extern int StorePathsThroughVertex(PATH *paths,
														EDGE *enter, NODES *vertex, EDGE *exit, 
														int **pathIndices,
														int **pathPositions);
extern int EdgeIsShortSelfCycle(EDGE *edge, int maxLength);

extern int EdgeConnectsShortTandemRepeats(EDGE *edge, int maxLength, 
																					EDGE **enterEdge, EDGE **repeatEdge, EDGE **exitEdge);
extern int FindShortTandemRepeats(NODES **vertices, int numVertices,
																	PATH *paths, int cycleLength);

extern int StraightenShortTandemRepeats(NODES ***vertices, int numVertex, 
																				PATH *paths, int cycleLength);
extern int MarkIndexedReadIntervalForRemoval(EDGE *edge, int pos);
extern int RemoveMarkedReadIntervals(EDGE *edge);
extern void RemovePath(PATH *paths, int pathIndex);
extern void CountN50(NODES **vertices, int numVertices);
