/***************************************************************************
 * Title:          over-repeat-new.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <param.h>
#include <zigzag.h>
#include <extfunc.h>
#include <clean_graph.h>

int	noshave, filterinp;
int	qualinp;
char	inpfile[100], outfile[100], seqfile[100], qualfile[100];
char	htmlout, ***content, caption[2000];
FILE	*flog;

void initenv(int argc, char **argv);
int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;
main(int argc, char **argv)
{
	int	i, j, k, l, m, n, m0, n0, n1, n2;
	int	reads, multip;
	int	num_class, num_edge;
	int	**num_pa;
	char	temp[100], temp1[100];
	ALIGN	**eq_class, *align;
	NODES	**badnode;
	EDGE	**edge, *edge1, *edge2, *edge3, *bal_edge1, *bal_edge2;
	int	nbad;
	int	num_bv, num_bredge;
	NODES	**vertex, *begin, *node, *node_next;
	NODES	**branch_vertex;
	IGRAPH	G, GB;
	READTABLE RT;
	IVERT	*it_vert;
	IEDGE	*it_edge;
	READINTERVAL	*readinterval;
	READPOSITION	*readposition;
	FILE	*fp, *fp1, *fp2, *fp3, *branchFP;

	CLEAN_PARAMS clean_params;
	LABELVERT T_e_mem, T_c_mem;

	readpar();
	random1(&idum);
	initenv(argc, argv);
	if(htmlout)	{
		flog = ckopen("EULER-report.html", "a");
	} else	{
		flog = ckopen("EULER-report.txt", "a");
	}

	/* parameters for cleaning the graph */
	clean_params.w = 5;
	clean_params.max_cov = 16;
	clean_params.C_B = BulgeLength;
	clean_params.L_b = BulgeCoverage;
	clean_params.C_W = WhirlLength;
	clean_params.L_c = ChimericCoverage;
	clean_params.L_e = 1;
	clean_params.do_zz = 1;              /* straighten zigzag paths */

	//	clean_params.T_c = ChimericTerm;
	//	clean_params.T_e = ErosionLength;

	clean_params.T_c = &T_c_mem;
	if (ChimericTerm >= 0) {
	  clean_params.T_c->depth = ChimericTerm;
	  clean_params.T_c->dir = 0;
	  clean_params.T_c->onepath = 1;     /* not used */
	} else {
	  clean_params.T_c->depth = -ChimericTerm;
	  clean_params.T_c->dir = 1;
	  clean_params.T_c->onepath = 0;     /* not used */
	}
	clean_params.T_e = &T_e_mem;

	if (ErosionLength >= 0) {
	  clean_params.T_e->depth = ErosionLength;
	  clean_params.T_e->dir = 0;
	  clean_params.T_e->onepath = 1;
	} else {
	  clean_params.T_e->depth = -ErosionLength;
	  clean_params.T_e->dir = 1;
	  clean_params.T_e->onepath = 0;
	}

/*	Input the length of the reads (required) */

	read_fasta_file(seqfile, &RT);

/*	Input equivalent readintervals between reads --
	see the format of the equivalent readinterval files	*/

	eq_class = (ALIGN **) ckalloc(2 * RT.num_seq * sizeof(ALIGN *));
	fp = ckopen(inpfile, "r");
	num_class = readclass(eq_class, RT.num_seq, fp);
	fclose(fp);

/*	Filter redundunt overlaps	*/
	if(filterinp)	{
		num_class = filter_readinterval(eq_class, &RT);
	}

/*	Make nucleotide graph with one procedure	*/

	makegraph(&G, &RT, eq_class);

/*	Clean up the graph
	The following program only mark every node and
	edge as "IN" --
  	Replace this program by true cleaning later	*/

/*
	cleangraph(vertex, num_vertex);
*/

/*	Glenn's code for graph cleaning	*/

        clean_graph_once(&G, &RT, &clean_params);

	n = m = 0;
	it_edge = it_e_new(&G, E_G);
	while((edge1 = it_e_next(it_edge)) != (EDGE *) NULL)	{
		if(edge1 -> subg_flag == 2)	{
			edge1 -> subg_flag = 1;
			n ++;
		}
		if(edge1 -> subg_flag > 0)	{
			m ++;
		}
	}

/*	map the read to graph H:
	(1) change the pointers (readlist) to the new nodes it corresponds to;
	(2) update the array of read positions of each node;
	(3) update the array of read intervals of each edge;
	readlist is not useful after this procedure.
	To build branching graph,   */

	mapread(&G, &RT);

/*	Free the space of the pointers from read positions to nodes	*/
	free_readlist(&RT);

/*	Free spaces for the pairwise alignments	*/
	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		while(eq_class[i])	{
			eq_class[i] = free_align(eq_class[i]);
		}
	}
	free((void **) eq_class);

/*	Remove edges with no multiplicity in nucleotide graph	*/
	n = 0;
	it_edge = it_e_new(&G, E_H);
	while((edge1 = it_e_next(it_edge)) != (EDGE *) NULL)	{
		if(edge1 -> multip == 0)	{
			edge1 -> subg_flag = 0;
			n ++;
		}
	}

/*	Building branching graph from nucleotide graph		*/
/*	(1) creating branching vertices;	*/

	num_bv = it_v_count(&G, V_B_H);
	branch_vertex = (NODES **) ckalloc(num_bv * sizeof(NODES *));
	it_vert = it_v_new(&G, V_B_H);
	i = n = 0;
	while((node = it_v_next(it_vert)) != (NODES *) NULL)	{
		branch_vertex[i] = (NODES *) ckalloc(1 * sizeof(NODES));
		branch_vertex[i] -> bal_node = node -> bal_node;
		branch_vertex[i] -> num_lastedge = 0;
		branch_vertex[i] -> nextedge = (EDGE **) ckalloc(node -> num_nextedge * sizeof(EDGE *));
		m = 0;
		for(j = 0; j < node -> num_nextedge; j ++)	{
			if(node -> nextedge[j] -> subg_flag == 1)	{
				branch_vertex[i] -> nextedge[m ++] = node -> nextedge[j];
			}
		}
		branch_vertex[i] -> num_nextedge = m;
/*	node -> link eq. index of branching vertices	*/
		node -> nlinks = i;
		n += branch_vertex[i] -> num_nextedge;
		i ++;
	}
	it_v_destroy(it_vert);
	num_bv = i;

/*	(2) connecting branching vertices by long edges;	*/

	edge = (EDGE **) ckalloc(n * sizeof(EDGE *));
	num_bredge = branch_graph(RT.num_seq, RT.new_len_seq, branch_vertex, edge, num_bv); 

/*	Free the space of nucleotide graph	*/
	destroygraph(&G);

/*	Nucleotide graph are stored in (vertex, num_vertex) while
	branching graph are stored in (branch_graph, num_bv)	*/

/*	Assign the complementary edges of each edge	*/

	for(i = 0; i < num_bredge; i ++)	{
		edge[i] -> bal_edge = find_bal_edge(edge[i], RT.len_seq, RT.num_seq, i);
	}

/*	Shave short ends linking sources & sinks	*/

	if(!noshave)	{
		num_bv = shave_graph_new(branch_vertex, num_bv, &RT, EndLength, SecondChimericCoverage);
	}

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}
	num_bredge = count_edge_simp(branch_vertex, num_bv, num_pa);

	sprintf(temp, "%s.single", seqfile);
	fp = ckopen(temp, "w");
	m = n = 0;
	for(i = 0; i < num_bv; i ++)	{
		for(j = 0; j < branch_vertex[i] -> num_nextedge; j ++)	{
			edge1 = branch_vertex[i] -> nextedge[j];
			if(edge1 -> end -> num_nextedge == 0 && branch_vertex[i] -> num_lastedge == 0 &&
			   edge1 -> multip == 1)	{
/*	Output single read edges	*/
				reads = edge1 -> readinterval -> eq_read;
				if(reads < RT.num_seq)	{
					writeseq(fp, RT.src_seq[reads], RT.src_name[reads], RT.len_seq[reads]);
					n ++;
					m += edge1 -> length;
				}
				erasedge(edge1);
			}
		}
	}
	fclose(fp);
	if(!htmlout)	{
		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "Summary after FragmentGluer\n");
		fprintf(flog, "Removed %d single read with total length: %d. They are stored in file: %s.\n", n, m, temp);
	} else	{
		print_section_head(flog, "Summary by FragmentGluer");
		sprintf(caption, "Removed %d single read with total length: %d. They are stored in file: %s.",
			 n, m, temp);
		print_line(flog, caption);
		fflush(flog);
	}

	num_bv = merge_graph(branch_vertex, num_bv, NULL, 0);
	num_bredge = count_edge_simp(branch_vertex, num_bv, num_pa);
	m = l = 0;
	for(i = 0; i < num_bv; i ++)	{
		for(j = 0; j < branch_vertex[i] -> num_nextedge; j ++)	{
			l += branch_vertex[i] -> nextedge[j] -> length;
			edge1 = branch_vertex[i] -> nextedge[j];
			if(branch_vertex[i] -> nextedge[j] -> length > m)	{
				m = branch_vertex[i] -> nextedge[j] -> length;
			}
		}
	}

/*	Remove the read intervals of the chimeric reads (optional 
	if skip building the graph for the second time	*/
/*	Skip this to keep the partial read intervals of reads	*/

	print_chimtable(flog, &RT);
	rem_chim(branch_vertex, num_bv, RT.chim, RT.num_chim, RT.num_seq);

/**********************************************************************
 * remove 0-coverage edges
 **********************************************************************/

	k = 0;
	for (i = 0; i < num_bv; i ++) {
		for (j = 0; j < branch_vertex[i] -> num_nextedge; j ++) {
			edge1 = branch_vertex[i] -> nextedge[j];
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
	num_bv = merge_graph(branch_vertex, num_bv, NULL, 0);

/* free list of chimeric reads	*/
	if (RT.chim)	free((void *) RT.chim);

/*	Make consensus of edges	*/
	initial_edge(branch_vertex, num_bv, RT.src_seq, RT.len_seq,RT.num_seq);

	m = l = 0;
	for(i = 0; i < num_bv; i ++)	{
		for(j = 0; j < branch_vertex[i] -> num_nextedge; j ++)	{
			l += branch_vertex[i] -> nextedge[j] -> length;
			if(branch_vertex[i] -> nextedge[j] -> length > m)	{
				m = branch_vertex[i] -> nextedge[j] -> length;
			}
		}
	}
	num_bredge = count_edge_simp(branch_vertex, num_bv, num_pa);
	if(!htmlout)	{
		fprintf(flog, "Summary of graph\n");
		fprintf(flog, "%d vertices %d edges (%d source %d sinks) remained.\n",
			 num_bv, num_bredge, num_pa[0][1], num_pa[1][0]);
		fprintf(flog, "total length %d (maximal %d).\n", l, m);
	}
	fflush(flog);

/*	Output graph & contigs	*/
	sprintf(temp, "%s.edge", seqfile);
	fp = ckopen(temp, "w");
	sprintf(temp, "%s.graph", seqfile);
	fp1 = ckopen(temp, "w");
	sprintf(temp, "%s.contig", seqfile);
	fp2 = ckopen(temp, "w");
	sprintf(temp, "%s.contig.ace", seqfile);
	fp3 = ckopen(temp, "w");
	sprintf(temp, "%s.bgraph", seqfile);
	branchFP = ckopen(temp, "w");
	
	output_contig(branch_vertex, num_bv, NULL, 0, RT.src_name, RT.src_seq, RT.num_seq, 
								fp, fp1, fp2, fp3, branchFP, NULL);
	fclose(fp);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);

/*	Output graphviz format graph	*/

	sprintf(temp, "%s", outfile);
	fp = ckopen(temp, "w");
	output_graph(branch_vertex, num_bv, fp);
	fclose(fp);

/*	Output read intervals in each edge	*/
	sprintf(temp, "%s.intv", seqfile);
	fp = ckopen(temp, "w");
	write_interval(branch_vertex, num_bv, fp);
	fclose(fp);

	if(!htmlout)	{
		print_text_line(flog, LINE_LENGTH);
	} else	{
		print_hl(flog);
	}
	fclose(flog);

	GB.nodes = branch_vertex;
	GB.num_nodes = num_bv;
	destroygraph(&GB);
	for(i = 0; i < MAX_BRA; i ++)	{
		free((void *) num_pa[i]);
	}
	free((void **) num_pa);
	free_readtable(&RT);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;	
	extern char *optarg;

	inpseq = outseq = qualinp = 0;
	noshave = 0;
	gap_k = 0;
	END_MERGE = 1;
	MIN_OVERLAP = 100;
	htmlout = 0;
	filterinp = 0;

	while ((copt=getopt(argc,argv,"i:o:s:c:p:q:k:C:b:B:w:T:t:e:E:m:Q:Sf:F:M:H")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'S':
			  noshave = 1;
			  continue;
			case 'b':
			  sscanf(optarg,"%d", &BulgeLength);
			  continue;
			case 'w':
			  sscanf(optarg,"%d", &WhirlLength);
			  continue;
			case 'H':
			  htmlout = 1;
			  continue;
			case 'E':
			  sscanf(optarg,"%d", &EndLength);
			  continue;
			case 'B':
			  sscanf(optarg,"%d", &BulgeCoverage);
			  continue;
			case 't':
			  sscanf(optarg,"%d", &ChimericTerm);
			  continue;
			case 'T':
			  sscanf(optarg,"%d", &ChimericCoverage);
			  continue;
			case 'Q':
			  sscanf(optarg,"%d", &SecondChimericCoverage);
			  continue;
			case 'e':
			  sscanf(optarg,"%d", &ErosionLength);
			  continue;
			case 's':
			  sscanf(optarg,"%s", seqfile);
			  continue;
			case 'q':
			  qualinp = 1;
			  sscanf(optarg,"%s", qualfile);
			  continue;
			case 'p':
			  sscanf(optarg,"%d", &MIN_OVERLAP);
			  continue;
			case 'k':
			  sscanf(optarg,"%d", &gap_k);
			  continue;
			case 'm':
			  sscanf(optarg,"%d", &END_LEG);
			  continue;
			case 'c':
			  sscanf(optarg,"%d", &LOW_COV);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'f':
			  filterinp = 1;
			  sscanf(optarg,"%lf", &FILTER_THRESH_PERC);
			  continue;
			case 'F':
			  sscanf(optarg,"%d", &FILTER_THRESH_ABSOLUTE);
			  continue;
			case 'M':
			  sscanf(optarg,"%d", &MIS_SCORE);
			  continue;
			default:
			  if(!htmlout)	{
				  printf("over-repeat-new -i InpFile -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP]\n");
				  printf("-i InpFile: The input file name of alignments\n");
				  printf("-s SeqFile: The input file name of reads\n");
			  	  printf("-o Outfile: The output file name of contigs\n");
				  printf("-c LOW_COV (optional): minimal coverage\n");
				  printf("-p MIN_OVERLAP (optional): minimal overlap length\n");
				  exit(-1);
			  } else	{
				  print_line(flog, "over-repeat-new -i InpFile -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP]");
				  print_line(flog, "-i InpFile: The input file name of alignments");
				  print_line(flog, "-s SeqFile: The input file name of reads");
			  	  print_line(flog, "-o Outfile: The output file name of contigs");
				  print_line(flog, "-c LOW_COV (optional): minimal coverage");
				  print_line(flog, "-p MIN_OVERLAP (optional): minimal overlap length");
				  exit(-1);
			  }
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		if(!htmlout)	{
			  printf("over-repeat-new -i InpFile -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP]\n");
			  printf("-i InpFile: The input file name of alignments\n");
			  printf("-s SeqFile: The input file name of reads\n");
		  	  printf("-o Outfile: The output file name of contigs\n");
			  printf("-c LOW_COV (optional): minimal coverage\n");
			  printf("-p MIN_OVERLAP (optional): minimal overlap length\n");
			  exit(-1);
		} else	{
			  print_line(flog, "over-repeat-new -i InpFile -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP]");
			  print_line(flog, "-i InpFile: The input file name of alignments");
			  print_line(flog, "-s SeqFile: The input file name of reads");
		  	  print_line(flog, "-o Outfile: The output file name of contigs");
			  print_line(flog, "-c LOW_COV (optional): minimal coverage");
			  print_line(flog, "-p MIN_OVERLAP (optional): minimal overlap length");
			  exit(-1);
		}
	}
}
