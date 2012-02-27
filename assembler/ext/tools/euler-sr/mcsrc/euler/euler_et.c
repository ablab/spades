/***************************************************************************
 * Title:          euler_et.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "stdinc.h"
#include "param.h"
#include "extfunc.h"
#include <assert.h>
char  noshave;
int	  qualinp;
int   DoStraighten, DoXCut;
char	htmlout, ***content, caption[2000];
char	inpfile[100], outfile[100], seqfile[100], qualfile[100],
  edgefile[100], graphfile[100], intvfile[100];
FILE	*flog;

void initenv(int argc, char **argv);
int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;
main(int argc, char **argv)
{
  int num_vertex;
  NODES **vertex;

  int num_edge;
  EDGE **edge;
  EDGE *edge1, *edge2;

  int **num_pa;
  int i, m, l, j, k, n;

  char temp[100];

  int num_path;
  PATH *path;

  int	*chim, num_chim;
	int v, e;

  READTABLE RT_mem, *RT=&RT_mem;
	
	DoStraighten = 0;
	DoXCut = 0;
  /**********************************************************************
   * Get inputs and parameters
   **********************************************************************/

  readpar();
  //  random1(&idum);
  initenv(argc, argv);
  if(htmlout)	{
		flog = ckopen("EULER-report.html", "a");
  } else	{
		flog = ckopen("EULER-report.txt", "a");
  }
  if(!htmlout)	{
  	print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "EULER-ET Equivalent transformation with reads:\n\n");
  } else	{
		print_section_head(flog, "Summary of EULER-ET");
	 
  }

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
   * Build some statistics on the number of edges in the graph.
   * Nothing is done with this...
   **********************************************************************/

  num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
  for (i = 0; i < MAX_BRA; i ++) {
    num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
  }
  num_edge = count_edge_simp(vertex, num_vertex, num_pa);

  /**********************************************************************
   * Build read paths
   **********************************************************************/

  num_chim = 0;
  chim = (int *) ckalloc(RT->num_seq * sizeof(int));
  path = (PATH *) ckalloc(2 * RT->num_seq * sizeof(PATH));
  num_path = readpath(vertex, &num_vertex, path,
											RT->len_seq, RT->num_seq, chim, &num_chim);
  free((void *) chim);

  /**********************************************************************
   * remove 0-coverage edges
   **********************************************************************/

	/* First mark edges for removal.
	 */
  k = 0;
  for (i = 0; i < num_vertex; i ++) {
    for (j = 0; j < vertex[i] -> num_nextedge; j ++) {
      edge1 = vertex[i] -> nextedge[j];
      if(edge1 -> multip == 0)	{
				edge2 = edge1 -> bal_edge;
				MarkEdgeForRemoval(edge1);
				//				edge1->deleted = 1;
				//				EraseEdgeAndPassingPaths(edge1, path, num_path);
				//				erasedge(edge1);
				k ++;

				if(edge2 && edge2 != edge1)	{
					/*					EraseEdgeAndPassingPaths(edge2, path, num_path);*/
					//					erasedge(edge2);
					k ++;
				}
      }
    }
  }
	/*
		Now actually remove the edges.
	*/
	EraseMarkedEdgesAndPaths(vertex, num_vertex, path, num_path);

  /* 
     Remove paths that are mapped to a single edge, they will not help with the equivalent 
     transformation.
  */
  printf("erased: %d paths\n", k);

  num_path = filter_path_read(path, RT->num_seq * 2, RT-> num_seq * 2 , 0);

  /* 
     Currently all paths are stored in edge lists. To enable equivalent transformation, 
     the paths are mapped to the corresponding vertices using 'set_path'.
  */
  set_path(vertex, num_vertex, path, num_path);

  /*
    Merge any paths that are simple.  A --> B --> C  _> A-->C
  */
  int prevNumVertex = num_vertex;
  num_vertex = merge_graph_path(vertex, num_vertex, path, num_path, LOW_COV_PATH);
  printf("prev n vertices: %d new n vertices: %d\n", prevNumVertex, num_vertex);
  if(!htmlout)	{
	  fprintf(flog, "num_path %d num_seq %d\n", num_path, RT -> num_seq);
  } else	{
	  sprintf(temp, "Summary of input: # path %d, # reads %d", num_path, RT -> num_seq);
	  print_line(flog, temp);
  }

  /**********************************************************************
   * equivalent transformation of the reads
   **********************************************************************/

  statspath(path, num_path);
	/*
	for (v= 0; v < num_vertex; v++) {
		for (e= 0 ; e < vertex[v]->num_nextedge; e++) {
			if (vertex[v]->nextedge[e]->multip != 
					vertex[v]->nextedge[e]->bal_edge->multip) {
				printf("edge: %d, multip: %d len: %d  bal: %d multip: %d\n", vertex[v]->nextedge[e]->index,
							 vertex[v]->nextedge[e]->multip,
							 vertex[v]->nextedge[e]->length,
							 vertex[v]->nextedge[e]->bal_edge->multip,
							 vertex[v]->nextedge[e]->bal_edge->length);
			}
		}
	}
	*/
  num_vertex = eqtrans_bal(&vertex, num_vertex, path, num_path, RT -> num_seq, 0);
	/*	
	for (v= 0; v < num_vertex; v++) {
		for (e= 0 ; e < vertex[v]->num_nextedge; e++) {
			if (vertex[v]->nextedge[e]->multip !=
					vertex[v]->nextedge[e]->bal_edge->multip) {
				printf("edge: %d, multip: %d %d bal: %d multip: %d\n", vertex[v]->nextedge[e]->index,
							 vertex[v]->nextedge[e]->multip,
							 vertex[v]->nextedge[e]->bal_edge->index,
							 vertex[v]->nextedge[e]->bal_edge->multip,
							 vertex[v]->nextedge[e]->bal_edge->length);
			}
		}
	}
	*/
  statspath(path, num_path);
	int read, ri;

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

  /**********************************************************************
   * remove 0-coverage edges
   **********************************************************************/

												 
  for (i = 0; i < num_vertex; i ++) {
    for (j = 0; j < vertex[i] -> num_nextedge; j ++) {
      edge1 = vertex[i] -> nextedge[j];
      if(edge1 -> multip == 0)	{
				edge2 = edge1 -> bal_edge;
				/*				printf("removing edge %d len: %d\n", edge1->index, edge1->length);*/
				erasedge(edge1);
				if(edge2 && edge2 != edge1)	{
					/*					printf("removing bal edge %d\n", edge2->index);*/
					erasedge(edge2);
				}
      }
    }
  }

  num_vertex = merge_graph_path(vertex, num_vertex, path, num_path, 0);

	int totalReads, totalLength, totalReadStarts;
	CountCoverage(vertex, num_vertex, &totalLength, &totalReadStarts);
	printf("removing low coverage again, %d %d\n", totalLength, totalReadStarts);
	RemoveLowCoverageEdges(vertex, num_vertex, path, num_path,
												 totalLength, totalReadStarts,
												 4, 5);


	if (DoStraighten) {
		num_vertex = StraightenShortTandemRepeats(&vertex, num_vertex, path, 40000);
		num_vertex = merge_graph_path(vertex, num_vertex, path, num_path, LOW_COV_PATH);
		num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	}

	/*	num_vertex = merge_graph_path(vertex, num_vertex, path, num_path, 0);*/

  num_vertex = eqtrans_bal(&vertex, num_vertex, path, num_path, RT -> num_seq, DoXCut);
  num_vertex = merge_graph_path(vertex, num_vertex, path, num_path, LOW_COV_PATH);

	CountN50(vertex, num_vertex);

  /**********************************************************************
   * Make consensus of edges
   **********************************************************************/
  /*
    consensus(vertex, num_vertex, RT->src_seq, RT->num_seq, RT->len_seq);
  */
	//	CountCoverage(vertex, num_vertex, &totalLength, &totalReads);
	//	CheckEdgeMultiplicities(vertex, num_vertex, totalLength, totalReads);
	//	num_vertex = RemoveLowCoverageEdges(vertex, num_vertex, path, num_path, totalLength, totalReads, 3.0, 5);
												 
	//	printf("total reads: %d  length: %d  average reads/nucleotide: %f\n", 
	//				 totalReads, totalLength, totalReads / (1.0*totalLength));

  initial_edge(vertex, num_vertex, RT->src_seq, RT->len_seq,RT->num_seq);

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

  num_edge = count_edge_simp(vertex, num_vertex, num_pa);
	/*
		printf("after sorting and stuff, num edge: %d\n", num_edge);
		CountN50(vertex, num_vertex);
	*/
  numtangle = (int*) ckalloc(num_edge * sizeof(int));
  maxmultip = (int*) ckalloc(num_edge * sizeof(int));
  avelength = (int*) ckalloc(num_edge * sizeof(int));
  mlength = (int*) ckalloc(num_edge * sizeof(int));
  /**********************************************************************
   * Output graph & contigs
   **********************************************************************/

  sprintf(temp, "%s.et", seqfile);
  output_contig_files(temp, num_vertex, vertex, path, num_path, RT);

  /**********************************************************************
   * Output graphviz format graph
   **********************************************************************/

  write_gvz_file(outfile, num_vertex, vertex, 1);

  /**********************************************************************
   * Output intervals
   **********************************************************************/

  sprintf(temp, "%s.et.intv", seqfile);
  write_interval_file(temp, num_vertex, vertex);

  if(!htmlout)	{
  	print_text_line(flog, LINE_LENGTH);
  } else	{
		print_hl(flog);
  }
  fclose(flog);

  /**********************************************************************
   * free memory
   **********************************************************************/

  for (i = 0; i < MAX_BRA; i ++) {
    free((void *) num_pa[i]);
  }
  free((void **) num_pa);

  free_path(num_path, path);
	free_graph(vertex, num_vertex);
  free((void **) vertex);
  free((void **) edge);
  free_readtable(RT);

  return(0);
}

void initenv(int argc, char **argv)
{
  int copt;
  int inpseq, outseq;
  extern char *optarg;

  noshave = 1;
  inpseq = outseq = qualinp = 0;
  MIN_OVERLAP = 100;
  htmlout = 0;

  while ((copt=getopt(argc,argv,"o:s:c:w:p:l:x:vE:Q:t:HSX")) != EOF) {
    switch(copt) {
		case 'S':
			DoStraighten = 1;
			continue;
		case 'X':
			DoXCut = 1;
			continue;
    case 's':
      inpseq = 1;
      sscanf(optarg,"%s", seqfile);
      continue;
    case 'p':
      sscanf(optarg,"%d", &MIN_OVERLAP);
      continue;
    case 'c':
      sscanf(optarg,"%d", &LOW_COV);
      continue;
    case 'w':
      sscanf(optarg,"%d", &LOW_COV_PATH);
      continue;
    case 'o':
      outseq = 1;
      sscanf(optarg,"%s", outfile);
      continue;
    case 'l':
      sscanf(optarg,"%d", &SMALL_EDGE);
      continue;
    case 'x':
      sscanf(optarg,"%d", &VERTEX_SIZE);
      continue;
    case 'E':
      sscanf(optarg,"%d", &EndLength);
      continue;
    case 't':
      sscanf(optarg,"%d", &ChimericTerm);
      continue;
    case 'Q':
      sscanf(optarg,"%d", &SecondChimericCoverage);
      continue;
    case 'v':
      noshave = 0;
      continue;
    case 'H':
      htmlout = 1;
      continue;
    default:
      if(!htmlout)	{
	      printf("euler_et -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP -x VERTEX_SIZE]\n");
	      printf("-s SeqFile: The input file name of reads\n");
	      printf("-o Outfile: The output file name of contigs\n");
	      printf("-c LOW_COV (%d): minimal coverage\n");
	      printf("-w LOW_COV_PATH (%d): remove paths if there are less than LOW_COV_PATH\n"
		     "            copies of the path present in the set of reads.\n");
	      printf("-p MIN_OVERLAP (optional): minimal overlap length\n");
	      exit(-1);
      } else	{
	      print_line(flog, "euler_et -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP -x VERTEX_SIZE]");
	      print_line(flog, "-s SeqFile: The input file name of reads");
	      print_line(flog, "-o Outfile: The output file name of contigs");
	      print_line(flog, "-c LOW_COV (optional): minimal coverage");
	      print_line(flog, "-p MIN_OVERLAP (optional): minimal overlap length");
	      exit(-1);
      }
    }
    optind--;
  }
  sprintf(edgefile, "%s.edge", seqfile);
  sprintf(graphfile, "%s.graph", seqfile);
  sprintf(intvfile, "%s.intv", seqfile);

  if (inpseq == 0 || outseq == 0) {
		if(!htmlout)	{
			printf("euler_et -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP -x VERTEX_SIZE]\n");
			printf("-s SeqFile: The input file name of reads\n");
			printf("-o Outfile: The output file name of contigs\n");
			printf("-c LOW_COV (optional): minimal coverage\n");
			printf("-p MIN_OVERLAP (optional): minimal overlap length\n");
		} else	{
			print_line(flog, "euler_et -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP -x VERTEX_SIZE]");
			print_line(flog, "-s SeqFile: The input file name of reads");
			print_line(flog, "-o Outfile: The output file name of contigs");
			print_line(flog, "-c LOW_COV (optional): minimal coverage");
		}
		exit(-1);
  }
}
