/***************************************************************************
 * Title:          euler_db.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <rule.h>
#include <param.h>
#include <extfunc.h>
#include <stdlib.h>
char    noshave;
int DoXCut, DoStraighten;
char	htmlout, caption[2000], ***content;
FILE	*flog;
extern void readpar();
extern void print_text_line(FILE *fp, int length);


int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;

void initenv(int argc, char **argv);
//int newsfpair(EDGE *begin_edge, EDGE *end_edge, SFPAIRS *SFP);
//void free_sfpairs(SFPAIRS *SFP);


void build_mp_paths(READTABLE *RT,
										MATEPAIRTABLE *MP,
										MATEPAIRRULES *MPR,
										int verbose,
										int *num_path,
										PATH *path,
										SFPAIRS *SFP, int mpType);

int	votethresh;

char outfile[100], seqfile[100], edgefile[100], graphfile[100],
  intvfile[100], rulefile[100], pairfile[100];

int mpType = -1;

main(int argc, char **argv)
{
  int i,j,k, l, m, n;
  int	*chim, num_chim;

  int num_vertex;
  NODES	**vertex;

  int	num_edge;
  EDGE	**edge, *begin_edge, *end_edge, *edge1, *edge2;

  int	num_path;
  PATH	*path;

  char temp[1000];

  SFPAIRS SFP_mem, *SFP=&SFP_mem;

  READTABLE RT_mem, *RT=&RT_mem;
  MATEPAIRTABLE MP_mem, *MP=&MP_mem;
  MATEPAIRRULES MPR_mem, *MPR=&MPR_mem;

  DoXCut = 0;
  DoStraighten = 0;

  /**********************************************************************
   * Get inputs and parameters
   **********************************************************************/

  readpar();
  initenv(argc, argv);

  if(htmlout)	{
    flog = ckopen("EULER-report.html", "a");
  } else	{
    flog = ckopen("EULER-report.txt", "a");
  }
  if(!htmlout)	{
    print_text_line(flog, LINE_LENGTH);
    fprintf(flog, "EULER-DB Equivalent transformation with mate-pairs:\n\n");
  } else	{
    /* Commented out print_hl */  
    /* print_hl(flog); */
    /* new function for printing section header */
    /* print_line(flog, "Summary of EULER-DB"); */
    print_section_head(flog, "Summary of EULER-DB");
  }

  /**********************************************************************
   * Read mate-pair naming rules file
   **********************************************************************/

  read_matepair_rules_file(rulefile, MPR);

  /**********************************************************************
   * Input the reads, their lengths, and names
   **********************************************************************/

  read_fasta_file(seqfile, RT);

  /**********************************************************************
   * Input the graph
   **********************************************************************/

  read_graph_file(edgefile, graphfile,
									&num_vertex, &vertex,
									&num_edge, &edge);

  /**********************************************************************
   * Input the read intervals in each edge
   **********************************************************************/

  read_interval_file(intvfile, num_vertex, vertex);

  /**********************************************************************
   * Input mate-pair information
   **********************************************************************/

  read_matepair_file(pairfile, MPR, RT, MP);

  /**********************************************************************
   * Build read paths
   **********************************************************************/

  num_chim = 0;
  chim = (int *) ckalloc(RT->num_seq * sizeof(int));
  path = (PATH *) ckalloc(2 * (RT->num_seq + MP->num_pair) * sizeof(PATH));
  num_path = readpath(vertex, &num_vertex, path, RT->len_seq, RT->num_seq, chim, &num_chim);
  if(!htmlout)	{
    fprintf(flog, "num_path %d num_seq %d\n", num_path, RT -> num_seq);
  } else	{
    sprintf(temp, "Summary of input: # path %d, # reads %d", num_path, RT -> num_seq);
    print_line(flog, temp);
  }
  free((void *) chim);

  /* Trim read paths such that the length of every read path is at most 2 */ 
  k = 0;
  for(i = 0; i < RT -> num_seq; i ++)	{
    if(path[i].len_path > 1)	{
      trimpath(&path[i], &path[i + RT -> num_seq], i, i + RT -> num_seq);
      k += 2;
    }
  }

  /* remove 0 multip edges */
  k = 0;
  for (i = 0; i < num_vertex; i ++) {
    for (j = 0; j < vertex[i] -> num_nextedge; j ++) {
      edge1 = vertex[i] -> nextedge[j];
      if(edge1 -> multip == 0)	{
				edge2 = edge1 -> bal_edge;
				erasedge(edge1);
				k ++;
				if(edge2 && edge2 != edge1)	{
					erasedge(edge2);
					k ++;
				}
      }
    }
  }
  statspath(path, num_path);

  /**********************************************************************
   * Build mate-pair paths
   **********************************************************************/
  int numReadPaths = num_path;

  build_mp_paths(RT, MP, MPR, 1,
								 &num_path, path, SFP, mpType);

  num_path = filter_path_read(path, numReadPaths, num_path, DoXCut);
  set_path(vertex, num_vertex, path, num_path);

  /**********************************************************************
   * equivalent transformation of the reads
   **********************************************************************/

  statspath(path, num_path);
  num_vertex = eqtrans_bal(&vertex, num_vertex, path, num_path, RT -> num_seq, 1);
  statspath(path, num_path);
  free_path(num_path, path);

  /**********************************************************************
   * Build read paths on the transformed graph
   **********************************************************************/

  num_chim = 0;
  chim = (int *) ckalloc(RT->num_seq * sizeof(int));
  path = (PATH *) ckalloc(2 * RT->num_seq * sizeof(PATH));
  num_path = readpath(vertex, &num_vertex, path, RT->len_seq, RT->num_seq, chim, &num_chim);
  free((void *) chim);

  /*********************************************************************
   * link EULER-SF edges -- must do this after equivalent transformation
   *********************************************************************/
  for(j = 0; j < MP->num_sf_pair; j ++)	{
    i=MP->sf_pair[j];
    if (path[MP->pair1[i]].len_path > 0) {
      begin_edge = path[MP->pair1[i]].edge[path[MP->pair1[i]].len_path - 1];
    } else {
      begin_edge = (EDGE *) NULL;
    }
    if (path[MP->pair2[i]].len_path > 0) {
      end_edge = path[RT->num_seq + MP->pair2[i]].edge[0];
    } else {
      end_edge = (EDGE *) NULL;
    }
    if(begin_edge && end_edge && begin_edge != end_edge)
      newsfpair(begin_edge, end_edge, SFP);
  }
  free_matepair(MP);

  if(!noshave)	{
    RT -> num_chim = 0;
    RT -> chim = (int *) NULL;
    num_vertex = shave_graph_new(vertex, num_vertex, RT, EndLength, SecondChimericCoverage);
    /*	Remove the read intervals of the chimeric reads (optional 
				if skip building the graph for the second time	*/
    /*	Skip this to keep the partial read intervals of reads	*/

    print_chimtable(flog, RT);
    rem_chim(vertex, num_vertex, RT -> chim, RT -> num_chim, RT -> num_seq);
  }

  /* remove 0 multip edges */
  k = 0;
  for (i = 0; i < num_vertex; i ++) {
    for (j = 0; j < vertex[i] -> num_nextedge; j ++) {
      edge1 = vertex[i] -> nextedge[j];
      if(edge1 -> multip == 0)	{
				edge2 = edge1 -> bal_edge;
				erasedge(edge1);
				k ++;
				if(edge2 && edge2 != edge1)	{
					k ++;
					erasedge(edge2);
				}
      }
    }
  }
  num_vertex = merge_graph(vertex, num_vertex, NULL, 0);

  m = l = 0;
  for (i = 0; i < num_vertex; i ++) {
    for (j = 0; j < vertex[i] -> num_nextedge; j ++) {
      edge1 = vertex[i] -> nextedge[j];
      l += edge1 -> length;
      if(edge1 -> length > m)	{
				m = edge1 -> length;
      }
      sortreadinterval(edge1 -> readinterval, edge1 -> multip);
      n = 0;
      for(k = 0; k < edge1 -> multip; k ++)	{
				n += edge1 -> readinterval[k].length;
      }
    }
  }

  /**********************************************************************
   * Make consensus of edges
   **********************************************************************/

  initial_edge(vertex, num_vertex, RT->src_seq, RT->len_seq,RT->num_seq);

  /**********************************************************************
   * Output graph & contigs
   **********************************************************************/

  sprintf(temp, "%s.db", seqfile);
  SFP->tot_edge = output_contig_files(temp, num_vertex, vertex, path, num_path, RT);
  SFP->all_edge = (EDGE **) ckalloc(2 * SFP->tot_edge * sizeof(EDGE *));
  for(i = 0; i < num_vertex; i ++)	{
    for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
      k = vertex[i] -> nextedge[j] -> start_cover;
      SFP->all_edge[k] = vertex[i] -> nextedge[j];
    }
  }

  /**********************************************************************
   * output scaffolding graph
   **********************************************************************/
  sprintf(temp, "%s.db.contig", seqfile);
  //  write_sf_gvz_file(temp, "", SFP);
  //  free_sfpairs(SFP);

  /**********************************************************************
   * Output intervals
   **********************************************************************/

  sprintf(temp, "%s.db.intv", seqfile);
  write_interval_file(temp, num_vertex, vertex);

  /**********************************************************************
   * Output graphviz format graph
   **********************************************************************/

  write_gvz_file(outfile, num_vertex, vertex, 1);

  if(!htmlout)	{
    print_text_line(flog, LINE_LENGTH);
  } else	{
    print_hl(flog);
  }
  fclose(flog);

  /**********************************************************************
   * free memory
   **********************************************************************/

  free_graph(vertex, num_vertex);
  free((void **) vertex);
  free((void **) edge);
  free_readtable(RT);

  return(0);
}

/*****************************************************************************
 * Build mate-pair paths
 *****************************************************************************/

void build_mp_paths(READTABLE *RT,
										MATEPAIRTABLE *MP,
										MATEPAIRRULES *MPR,
										int verbose,

										/* array of mate-pair paths */
										/* input & output: index to current path number */
										int *num_path,
										PATH *path,
										SFPAIRS *SFP, int mpType
										)
{
  int i,j;
  int m;
  int n1,n2,n;
  int l;
  double	r;
  int el[9];
  int	num_zero_pair, num_bad_pair, num_one_pair, num_more_pair, num_good_pair, num_rep_pair;
  int	b1, b2, e1, e2, clone_len, delta;
  EDGE	*begin_edge, *end_edge;
  int dist[500], n_pair[500][10], num_pair_tot[500];
  int len;
  char mark, temp[100];

  int	num_path_temp;
  PATH path_temp;


  /**********************************************************************
   * Initialize statistics on how well mate-pairs map into graph
   **********************************************************************/

  for (i = 0; i < MPR->ntypepair; i ++) {
    num_pair_tot[i] = dist[i] = 0;
    for (j = 0; j < 10; j ++) {
      n_pair[i][j] = 0;
    }
  }


  /**********************************************************************
   * el[L] = # paths of length L+2 (when L<5), or length >= 7 (when L=5)
   **********************************************************************/

  for (i = 0; i < 9; i ++)	el[i] = 0;

  /**********************************************************************
   * init list of scaffolding pairs for EULER-SF
   **********************************************************************/

  init_sfpairs(SFP, MP, 0);

  /**********************************************************************
   * build list of mate-pair paths
   **********************************************************************/

  path_temp.edge = (EDGE **) ckalloc(MAX_TMP_LEG * sizeof(EDGE *));
  path_temp.readindex = -1;
  num_good_pair = num_one_pair = num_bad_pair = num_more_pair = num_zero_pair = num_rep_pair = 0;
  std::set<EDGE*> traversedEdges;

  for (i = 0; i < MP->num_pair; i ++) {
    if (i % 1000 == 0) {
      printf(".");
      fflush(stdout);
    }
    if (mpType != -1 && MP->ntype[i] != mpType) {
      continue;
    }
    /*	allow 20% error in the clone length measure	*/
    delta = (int) (MP -> max_dist[i] * 0.0);

    int forw, revs;
    forw = MP->pair1[i];
    revs = MP->pair2[i];

    if (path[MP->pair1[i]].len_path > 0) {
      // begin_edge is the last edge mapped to the read MP->pair1
      begin_edge = path[MP->pair1[i]].edge[path[MP->pair1[i]].len_path - 1];
      b1 = path[MP->pair1[i]].begin_length;
      e1 = path[MP->pair1[i]].end_length;
    } else {
      begin_edge = (EDGE *) NULL;
    }
    if (path[MP->pair2[i]].len_path > 0) {
      // end_edge is the first edge mapped to the read MP->pair2
      end_edge = path[RT->num_seq + MP->pair2[i]].edge[0];
      b2 = path[RT -> num_seq + MP->pair2[i]].begin_length;
      e2 = path[RT -> num_seq + MP->pair2[i]].end_length;
    } else {
      end_edge = (EDGE *) NULL;
    }
    if (!begin_edge || !end_edge) {
      // num_rep_pair is output in statistics
      // cannot map both ends of the mate-pair, so it is not used to resolve anything.
      num_rep_pair ++;
    } else if (begin_edge == end_edge && path[MP->pair1[i]].len_path == 1 && path[MP->pair2[i]].len_path == 1) {
      // the mate-pairs map to the same edge
      // use minimum to reflect which mate-pair is first on the edge
      clone_len = begin_edge -> length - min(b1 + e2, b2 + e1);
      /*
				printf("pair %d is on unique edge %d  b1: %d e1: %d b2: %d e2: %d\n",
				i, begin_edge, b1, e1, b2, e2);
      */
      num_pair_tot[MP->ntype[i]] ++;
      // store distance between clones for bookkeeping purposes
      dist[MP->ntype[i]] += clone_len;

      // Check if the distance between the mate-pairs is a valid distance
      if (clone_len >= MP->min_dist[i] && clone_len <= MP->max_dist[i]) {
				n_pair[MP->ntype[i]][0] ++;
      } else {
				// invalid distance (either too long or too short), set l to the error, and do 
				// some book keeping to record the types of error found.
				l = min(abs(clone_len - MP->min_dist[i]), abs(clone_len - MP->max_dist[i]));
				// r is the ratio of error to maximum distance
				// record this for classifying the error
				r = ((double) l) / MP->max_dist[i];
				if (r < 0.1) {
					if (clone_len < MP->min_dist[i]) {
						n_pair[MP->ntype[i]][1] ++;
					} else {
						n_pair[MP->ntype[i]][2] ++;
					}
				} else if (r < 0.2) {
					if (clone_len < MP->min_dist[i]) {
						n_pair[MP->ntype[i]][3] ++;
					} else {
						n_pair[MP->ntype[i]][4] ++;
					}
				} else {
					if (clone_len < MP->min_dist[i]) {
						n_pair[MP->ntype[i]][5] ++;
					} else {
						n_pair[MP->ntype[i]][6] ++;
					}
				}
      }
      // Allow some flexibility in the distance between mate-pairs.
      if (clone_len > MP->min_dist[i] - delta && clone_len < MP->max_dist[i] + delta) {
				num_one_pair ++;
      } else {
				num_bad_pair ++;
      }
    } else if(begin_edge == end_edge)	{
      // The begin edge is the same as the end edge, but the reads before the mate-pairs spans more 
      // than one edge.
      clone_len = e1 + RT -> len_seq[MP->pair1[i]] + b2 + RT -> len_seq[MP->pair2[i]] - begin_edge -> length;
      if (clone_len > MP->min_dist[i] - delta && clone_len < MP->max_dist[i] + delta) {
				len = e1 + RT -> len_seq[MP->pair1[i]] + b2 + RT -> len_seq[MP->pair2[i]];
				traversedEdges.clear();
				num_path_temp = locpath(begin_edge, end_edge, 
																MP->min_dist[i] - delta, MP->max_dist[i] + delta, len, 
																&path_temp, traversedEdges);
				ResetTraversedEdges(traversedEdges);
				if(num_path_temp == 0)	{
					n = path[MP->pair1[i]].len_path + path[MP->pair2[i]].len_path - 1;
					path[*num_path].edge = (EDGE **) ckalloc(n * sizeof(EDGE *));
					path[*num_path].pairindex[0] = MP->pair1[i] + 1;
					path[*num_path].pairindex[1] = RT -> num_seq + MP->pair2[i] + 1;
					n = 0;
					for (j = 0; j < path[MP->pair1[i]].len_path; j ++) {
						path[*num_path].edge[n ++] = path[MP->pair1[i]].edge[j];
					}
					for (j = 1; j < path[RT -> num_seq + MP->pair2[i]].len_path; j ++) {
						path[*num_path].edge[n ++] = path[RT -> num_seq + MP->pair2[i]].edge[j];
					}
					path[*num_path].len_path = n;
					(*num_path) ++;
					path[*num_path].edge = (EDGE **) ckalloc(n * sizeof(EDGE *));
					path[*num_path].pairindex[0] = MP->pair2[i] + 1;
					path[*num_path].pairindex[1] = RT -> num_seq + MP->pair1[i] + 1;
					for (j = 0; j < path[*num_path - 1].len_path; j ++) {
						path[*num_path].edge[j] = path[*num_path - 1].edge[path[*num_path - 1].len_path - 1 - j] -> bal_edge;
					}
					path[*num_path].len_path = path[*num_path - 1].len_path;
					if (path[*num_path].len_path > 6)	el[5] ++;
					else el[path[*num_path].len_path - 2] ++;
					(*num_path) ++;
					num_good_pair ++;
				} else	{
					num_more_pair ++;
				}
      } else {
				num_bad_pair ++;
      }
    } else {
      // There are several edges between the ends of the mate-pair. Look to see if they 
      // are 
      len = e1 + RT -> len_seq[MP->pair1[i]] + b2 + RT -> len_seq[MP->pair2[i]] - VERTEX_SIZE;
      path_temp.len_path = 0;
      traversedEdges.clear();
      num_path_temp = locpath(begin_edge, end_edge, 
															MP->min_dist[i] - delta, MP->max_dist[i] + delta, len, 
															&path_temp, traversedEdges);
      ResetTraversedEdges(traversedEdges);
      if (num_path_temp == 0) {
        MP->sf_pair[MP->num_sf_pair ++] = i;
				num_zero_pair ++;
      } else if (num_path_temp > 1) {
				num_more_pair ++;
      } else {
				n1 = path[MP->pair1[i]].len_path;
				n2 = path[MP->pair2[i]].len_path;
				n = n1 + n2;
				for (j = 0; j < path_temp.len_path - 1; j ++) {
					n ++;
          len += path_temp.edge[j] -> length - VERTEX_SIZE;
				}
				if (len <= MP->max_dist[i] + delta) {
					path[*num_path].edge = (EDGE **) ckalloc(n * sizeof(EDGE *));
					path[*num_path].pairindex[0] = MP->pair1[i] + 1;
					// Determine the entire mate-path.  This is the path for 
					// read1-> gap -> read2.

					// Store the entire path in 'path'
					path[*num_path].pairindex[1] = RT -> num_seq + MP->pair2[i] + 1;
					n = 0;

					for (j = 0; j < path[MP->pair1[i]].len_path; j ++) {
						path[*num_path].edge[n ++] = path[MP->pair1[i]].edge[j];
					}
					for (j = 0; j < path_temp.len_path - 1; j ++) {
						path[*num_path].edge[n ++] = path_temp.edge[j];
					}
					for (j = 0; j < path[RT -> num_seq + MP->pair2[i]].len_path; j ++) {
						path[*num_path].edge[n ++] = path[RT -> num_seq + MP->pair2[i]].edge[j];
					}

					// Now store the reverse complement of this path.
					path[*num_path].readindex = 0;
					path[*num_path].len_path = n;
					(*num_path) ++;
					path[*num_path].edge = (EDGE **) ckalloc(n * sizeof(EDGE *));
					path[*num_path].pairindex[0] = MP->pair2[i] + 1;
					path[*num_path].pairindex[1] = RT -> num_seq + MP->pair1[i] + 1;
					for (j = 0; j < path[*num_path - 1].len_path; j ++) {
						path[*num_path].edge[j] = path[*num_path - 1].edge[path[*num_path - 1].len_path - 1 - j] -> bal_edge;
					}
					path[*num_path].len_path = path[*num_path - 1].len_path;
					if (path[*num_path].len_path > 6)	el[5] ++;
					else el[path[*num_path].len_path - 2] ++;
					(*num_path) ++;
					num_good_pair ++;
				} else {
					num_zero_pair ++;
				}
      }
    }
  }

  free((void **) path_temp.edge);

  if (!verbose) return;

  /**********************************************************************
   * Report on mate-pair paths
   **********************************************************************/

  if(!htmlout)	{
    fprintf(flog, "Summary of EULER-DB mate-pair mapping:\n%d mate-pairs = (%d+%d+%d+%d+%d+%d)\n", MP->num_pair,
						num_one_pair, num_good_pair, num_more_pair, num_bad_pair, num_zero_pair, num_rep_pair);
    fprintf(flog, "%d mate-pairs are mapped to the same edge with acceptable distance;\n", num_one_pair);
    fprintf(flog, "%d mate-pairs are mapped to different edges with unique acceptable distance -- used to define mate-paths;\n", num_good_pair);
    fprintf(flog, "Distribution of the length of mate-paths (number of edges it contains):\n");
    print_text_line(flog, LINE_LENGTH);
    fprintf(flog, "Length             2        3        4        5        6        >6\n");
    fprintf(flog, "# paths      %4d %8d %8d %8d %8d %8d\n", el[0], el[1], el[2], el[3], el[4], el[5]);
    print_text_line(flog, LINE_LENGTH);
    fprintf(flog, "%d mate-pairs are mapped to different edges with non-unique acceptable distance (ignored);\n", num_more_pair);
    fprintf(flog, "%d mate-pairs are mapped to the same edge with unacceptable distance;\n", num_bad_pair);
    fprintf(flog, "Warning: a large number of such pairs indicates errors in naming rules or insert size estimate.\n");
    fprintf(flog, "%d mate-pairs are mapped to different edges with no path or path with unacceptable distance -- some will be used in EULER-SF.\n", num_zero_pair);
    fprintf(flog, "%d mate-pairs have at least one read that is in the repeat region or cannot be mapped to the contigs;\n", num_rep_pair);
    print_text_line(flog, LINE_LENGTH);
    fprintf(flog, "EULER-DB insert length statistics based on number of mate-pairs mapped to the same edge:\n");
    print_text_line(flog, LINE_LENGTH);
    fprintf(flog, "               ");
    for (i = 0; i < MPR->ntypepair; i ++) {
      fprintf(flog, " %c-%c (plates %d-%d)  ",
							MPR->pairrule[i][0], MPR->pairrule[i][1], MPR->platerule[i][0], MPR->platerule[i][1]);
    }
    fprintf(flog, "\n");
    fprintf(flog, "Distance range ");
    for (i = 0; i < MPR->ntypepair; i ++) {
      fprintf(flog, "    %5d-%5d    ", MPR->pairrange[i][0], MPR->pairrange[i][1]);
    }
    fprintf(flog, "\n");
    fprintf(flog, "Average        ");
    for (i = 0; i < MPR->ntypepair; i ++) {
      if (num_pair_tot[i] == 0)	m = 0;
      else m = dist[i] / num_pair_tot[i];
      fprintf(flog, "     %6d       ", m);
    }
    fprintf(flog, "\n");
    fprintf(flog, ">20%% above/below ");
    for (i = 0; i < MPR->ntypepair; i ++) {
      fprintf(flog, " %6d/%-6d  ", n_pair[i][6], n_pair[i][5]);
    }
    fprintf(flog, "\n");
    fprintf(flog, "10-20%% above/below ");
    for (i = 0; i < MPR->ntypepair; i ++) {
      fprintf(flog, " %6d/%-6d ", n_pair[i][4], n_pair[i][3]);
    }
    fprintf(flog, "\n");
    fprintf(flog, "0-10%% above/below ");
    for (i = 0; i < MPR->ntypepair; i ++) {
      fprintf(flog, " %6d/%-6d  ", n_pair[i][2], n_pair[i][1]);
    }
    fprintf(flog, "\n");
    fprintf(flog, "Within the range");
    for (i = 0; i < MPR->ntypepair; i ++) {
      fprintf(flog, "     %6d       ", n_pair[i][0]);
    }
    fprintf(flog, "\n");
    fprintf(flog, "Total           ");
    for (i = 0; i < MPR->ntypepair; i ++) {
      fprintf(flog, "     %6d       ", num_pair_tot[i]);
    }
    fprintf(flog, "\n");
    print_text_line(flog, LINE_LENGTH);
    fprintf(flog, "Warning: In case of a skewed distribution, you may want to reconsider your estimate of insert length.\n");
  } else  {
    print_line(flog, "Summary of EULER-DB mate-pair mapping");
    sprintf(temp, "Totally %d mate-pairs = (%d+%d+%d+%d+%d+%d)", MP->num_pair,
						num_one_pair, num_good_pair, num_more_pair, num_bad_pair, num_zero_pair, num_rep_pair);
    print_line(flog, temp);
    sprintf(temp, "%d mate-pairs are mapped to the same edge with acceptable distance.", num_one_pair);
    print_line(flog, temp);
    sprintf(temp, "%d mate-pairs are mapped to different edges with unique acceptable distance -- used to define mate-paths.", num_good_pair);
    print_line(flog, temp);
    content = allocate_content(2, 7, 50);
    strcpy(content[0][0], "Length");
    strcpy(content[1][0], "# paths");
    for(i = 1; i < 7; i ++)	{
      if(i == 6)	{
				sprintf(content[0][i], ">6");
      } else	{
				sprintf(content[0][i], "%d", i + 1);
      }
      sprintf(content[1][i], "%d", el[i - 1]);
    }
    sprintf(caption, "Distribution of the length of mate-paths (number of edges it contains)");
    print_table(flog, 2, 7, content, caption);
    content = free_content(content, 2, 7);
    sprintf(temp, "%d mate-pairs are mapped to different edges with non-unique acceptable distance (ignored).", num_more_pair);
    print_line(flog, temp);
    sprintf(temp, "%d mate-pairs are mapped to the same edge with unacceptable distance.\n", num_bad_pair);
    print_line(flog, temp);
    print_line(flog, "Warning: a large number of such pairs indicates errors in naming rules or insert size estimate.");
    sprintf(temp, "%d mate-pairs are mapped to different edges with no path or path with unacceptable distance -- some will be used in EULER-SF.", num_zero_pair);
    print_line(flog, temp);
    sprintf(temp, "%d mate-pairs have at least one read that is in the repeat region or cannot be mapped to the contigs.", num_rep_pair);
    print_line(flog, temp);

    /*	Print table of statistics of insert length	*/
    content = allocate_content(8, MPR->ntypepair + 1, 50);
    for (i = 0; i < MPR->ntypepair; i ++) {
      sprintf(content[0][i + 1], "%c-%c (plates %d-%d)  ",
							MPR->pairrule[i][0], MPR->pairrule[i][1], MPR->platerule[i][0], MPR->platerule[i][1]);
    }
    strcpy(content[1][0], "Distance range");
    for (i = 0; i < MPR->ntypepair; i ++) {
      sprintf(content[1][i + 1], "%d-%d", MPR->pairrange[i][0], MPR->pairrange[i][1]);
    }
    strcpy(content[2][0], "Average");
    for (i = 0; i < MPR->ntypepair; i ++) {
      if (num_pair_tot[i] == 0)	m = 0;
      else m = dist[i] / num_pair_tot[i];
      sprintf(content[2][i + 1], "%d", m);
    }
    strcpy(content[3][0], ">20% above/below");
    for (i = 0; i < MPR->ntypepair; i ++) {
      sprintf(content[3][i + 1], " %d/%-d  ", n_pair[i][6], n_pair[i][5]);
    }
    strcpy(content[4][0], "10-20% above/below");
    for (i = 0; i < MPR->ntypepair; i ++) {
      sprintf(content[4][i + 1], " %d/%-d ", n_pair[i][4], n_pair[i][3]);
    }
    strcpy(content[5][0], "0-10% above/below");
    for (i = 0; i < MPR->ntypepair; i ++) {
      sprintf(content[5][i + 1], " %d/%-d  ", n_pair[i][2], n_pair[i][1]);
    }
    strcpy(content[6][0], "Within the range");
    for (i = 0; i < MPR->ntypepair; i ++) {
      sprintf(content[6][i + 1], "%d", n_pair[i][0]);
    }
    strcpy(content[7][0], "Total");
    for (i = 0; i < MPR->ntypepair; i ++) {
      sprintf(content[7][i + 1], "%d", num_pair_tot[i]);
    }
    sprintf(caption, "EULER-DB insert length statistics based on number of mate-pairs mapped to the same edge.");
    print_table(flog, 8, MPR->ntypepair + 1, content, caption);
    content = free_content(content, 8, MPR->ntypepair + 1);
    print_line(flog, "Warning: In case of a skewed distribution, you may want to reconsider your estimate of insert length.");
    print_emptyline(flog);
  }
}

/*****************************************************************************
 * Parse command line switches
 *****************************************************************************/


void initenv(int argc, char **argv)
{
  int copt;
  int inpseq, outseq, inppair;
  extern char *optarg;

  noshave = 1;
	DoXCut  = 0;
  inpseq = outseq = inppair = 0;
  overlaplen = 20;
  votethresh = 2;
  htmlout = 0;
  strcpy(rulefile, "name.rul");

  while ((copt=getopt(argc,argv,"i:e:g:m:k:vr:o:c:w:l:x:p:t:E:Q:H:M:XS")) != EOF) {
    switch(copt) {
    case 'i':
      inpseq = 1;
      sscanf(optarg,"%s", seqfile);
      sprintf(edgefile, "%s.et.edge", seqfile);
      sprintf(graphfile, "%s.et.graph", seqfile);
      sprintf(intvfile, "%s.et.intv", seqfile);
      continue;
		case 'X':
			DoXCut = 1;
			continue;
    case 'S':
      DoStraighten = 1;
      continue;
    case 'E':
      sscanf(optarg,"%d", &EndLength);
      continue;
    case 'e':
      sscanf(optarg,"%s", edgefile);
      continue;
    case 'p':
      sscanf(optarg,"%s", intvfile);
      continue;
    case 'w':
      sscanf(optarg,"%d", &LOW_COV_PATH);
      continue;
    case 'c':
      sscanf(optarg,"%d", &LOW_COV);
      continue;
    case 't':
      sscanf(optarg,"%d", &ChimericTerm);
      continue;
    case 'g':
      sscanf(optarg,"%s", graphfile);
      continue;
    case 'M':
      sscanf(optarg, "%d", &mpType);
      continue;
    case 'm':
      inppair = 1;
      sscanf(optarg,"%s", pairfile);
      continue;
    case 'o':
      outseq = 1;
      sscanf(optarg,"%s", outfile);
      continue;
    case 'r':
      sscanf(optarg,"%s", rulefile);
      continue;
    case 'k':
      sscanf(optarg,"%d", &overlaplen);
      continue;
    case 'Q':
      sscanf(optarg,"%d", &SecondChimericCoverage);
      continue;
    case 'l':
      sscanf(optarg,"%d", &SMALL_EDGE);
      continue;
    case 'v':
      noshave = 0;
      continue;
    case 'x':
      sscanf(optarg,"%d", &VERTEX_SIZE);
      continue;
    case 'H':
      htmlout = 1;
      continue;
    default:
      if(!htmlout)	{
				printf("euler_db -i SeqFile [-e EdgeFile -g GraphFile -p IntvFile] -m Matefile [-k overlaplen] [-v votethresh -x VERTEX_SIZE] \n");
				printf("-i SourceSeqFile: The input file name of reads\n");
				printf("-e EdgeFile (optional): input edge file name\n");
				printf("-g GraphFile (optional): input graph file name\n");
				printf("-p IntvFile (optional): input path file name\n");
				printf("-m MateFile (optional): input mate file name\n");
				printf("-k k-tuple (optional): length of k-tuple\n");
				exit(-1);
      } else {
				print_line(flog, "euler_db -i SeqFile [-e EdgeFile -g GraphFile -p IntvFile] -m Matefile [-k overlaplen] [-v votethresh -x VERTEX_SIZE]");
				print_line(flog, "-i SourceSeqFile: The input file name of reads");
				print_line(flog, "-e EdgeFile (optional): input edge file name");
				print_line(flog, "-g GraphFile (optional): input graph file name");
				print_line(flog, "-p IntvFile (optional): input path file name");
				print_line(flog, "-m MateFile (optional): input mate file name");
				print_line(flog, "-k k-tuple (optional): length of k-tuple");
				exit(-1);
      }
    }
    optind--;
  }

  if (inpseq == 0 || outseq == 0 || inppair == 0) {
    if(!htmlout)	{
      printf("euler_db -i SeqFile [-e EdgeFile -g GraphFile -p IntvFile] -m Matefile [-k overlaplen] [-v votethresh -x VERTEX_SIZE] \n");
      printf("-i SourceSeqFile: The input file name of reads\n");
      printf("-e EdgeFile (optional): input edge file name\n");
      printf("-g GraphFile (optional): input graph file name\n");
      printf("-p IntvFile (optional): input path file name\n");
      printf("-m MateFile (optional): input mate file name\n");
      printf("-k k-tuple (optional): length of k-tuple\n");
    } else {
      print_line(flog, "euler_db -i SeqFile [-e EdgeFile -g GraphFile -p IntvFile] -m Matefile [-k overlaplen] [-v votethresh -x VERTEX_SIZE]");
      print_line(flog, "-i SourceSeqFile: The input file name of reads");
      print_line(flog, "-e EdgeFile (optional): input edge file name");
      print_line(flog, "-g GraphFile (optional): input graph file name");
      print_line(flog, "-p IntvFile (optional): input path file name");
      print_line(flog, "-m MateFile (optional): input mate file name");
      print_line(flog, "-k k-tuple (optional): length of k-tuple");
    }
    exit(-1);
  }
}
