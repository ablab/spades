/***************************************************************************
 * Title:          euler_cons.c
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
#include <extfunc.h>

int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;
int	qualinp;
char	inpfile[400], outfile[400], seqfile[400], qualfile[400],
  edgefile[400], graphfile[400], intvfile[400], acefile[400], contigfile[400];
char	htmlout, caption[2000], ***content;
FILE	*flog;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
  int num_vertex;
  NODES **vertex;

  int num_edge;
  EDGE **edge;
  EDGE *edge1;

  int **num_pa;
  int i, m, l, j;

  char temp[400];

  int num_path;
  PATH *path;

  int	*chim, num_chim;

  READTABLE RT_mem, *RT=&RT_mem;

  FILE *fp, *fp1, *fp2;


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
	fprintf(flog, "EULER-Cons Consense step\n\n");
  } else	{
	print_section_head(flog, "Summary of EULER-Cons");
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
   * get the topology of the graph
   **********************************************************************/

  num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
  for (i = 0; i < MAX_BRA; i ++) {
    num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
  }
  num_edge = count_edge_simp(vertex, num_vertex, num_pa);

  /**********************************************************************
   * Make consensus of edges
   **********************************************************************/
  /*
    consensus(vertex, num_vertex, RT->src_seq, RT->num_seq, RT->len_seq);
  */
  initial_edge(vertex, num_vertex, RT->src_seq, RT->len_seq,RT->num_seq);

/*	Output contig files	*/

  fp = ckopen(edgefile, "w");
  fp1 = ckopen(contigfile, "w");
  fp2 = ckopen(acefile, "w");
  output_contig_cons(vertex, num_vertex, RT->src_name, RT->src_seq, RT->num_seq, fp, fp1, fp2);
  fclose(fp);
  fclose(fp1);
  fclose(fp2);

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

  inpseq = outseq = qualinp = 0;
  MIN_OVERLAP = 100;
  htmlout = 0;

  sprintf(edgefile, "%s.edge", seqfile);
  sprintf(graphfile, "%s.graph", seqfile);
  sprintf(intvfile, "%s.intv", seqfile);
  sprintf(acefile, "%s.contig.ace", seqfile);
  sprintf(contigfile, "%s.contig", seqfile);

  while ((copt=getopt(argc,argv,"s:g:i:e:c:a:o:H")) != EOF) {
    switch(copt) {
    case 's':
      inpseq = 1;
      sscanf(optarg,"%s", seqfile);
      continue;
    case 'g':
      sscanf(optarg,"%s", graphfile);
      continue;
    case 'o':
      sscanf(optarg,"%s", outfile);
      continue;
    case 'i':
      sscanf(optarg,"%s", intvfile);
      continue;
    case 'e':
      outseq = 1;
      sscanf(optarg,"%s", edgefile);
      continue;
    case 'c':
      sscanf(optarg,"%s", contigfile);
      continue;
    case 'a':
      sscanf(optarg,"%s", acefile);
      continue;
    case 'H':
      htmlout = 1;
      continue;
    default:
      if(!htmlout)	{
   	      printf("euler_cons -s SeqFile -g GraphFile -i IntvFile -e EdgeFile -c ContigFle -a AceFile \n");
	      printf("-s SeqFile: The input file name of reads\n");
	      printf("-g GraphFile: The input file name of graph\n");
	      printf("-i IntvFile: The input file name of intervals\n");
	      printf("-e EdgeFile: The (input)output file name of edge sequences\n");
	      printf("-c ContigFile: The output file name of contig sequences\n");
	      printf("-a AceFile: The output .ace file name\n");
	      printf("-o GVZFile: The output file name of graph\n");
       } else	{
   	      print_line(flog, "euler_cons -s SeqFile -g GraphFile -i IntvFile -e EdgeFile -c ContigFle -a AceFile");
	      print_line(flog, "-s SeqFile: The input file name of reads");
	      print_line(flog, "-g GraphFile: The input file name of graph");
	      print_line(flog, "-i IntvFile: The input file name of intervals");
	      print_line(flog, "-e EdgeFile: The (input)output file name of edge sequences");
	      print_line(flog, "-c ContigFile: The output file name of contig sequences");
	      print_line(flog, "-a AceFile: The output .ace file name");
	      print_line(flog, "-o GVZFile: The output file name of graph");
       }
       exit(-1);
    }
    optind--;
  }
  if (inpseq == 0 || outseq == 0) {
      if(!htmlout)	{
   	      printf("euler_cons -s SeqFile -g GraphFile -i IntvFile -e EdgeFile -c ContigFle -a AceFile \n");
	      printf("-s SeqFile: The input file name of reads\n");
	      printf("-g GraphFile: The input file name of graph\n");
	      printf("-i IntvFile: The input file name of intervals\n");
	      printf("-e EdgeFile: The (input)output file name of edge sequences\n");
	      printf("-c ContigFile: The output file name of contig sequences\n");
	      printf("-a AceFile: The output .ace file name\n");
	      printf("-o GVZFile: The output file name of graph\n");
       } else	{
   	      print_line(flog, "euler_cons -s SeqFile -g GraphFile -i IntvFile -e EdgeFile -c ContigFle -a AceFile");
	      print_line(flog, "-s SeqFile: The input file name of reads");
	      print_line(flog, "-g GraphFile: The input file name of graph");
	      print_line(flog, "-i IntvFile: The input file name of intervals");
	      print_line(flog, "-e EdgeFile: The (input)output file name of edge sequences");
	      print_line(flog, "-c ContigFile: The output file name of contig sequences");
	      print_line(flog, "-a AceFile: The output .ace file name");
	      print_line(flog, "-o GVZFile: The output file name of graph");
       }
       exit(-1);
  }
}
