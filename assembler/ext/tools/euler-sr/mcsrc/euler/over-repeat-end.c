/***************************************************************************
 * Title:          over-repeat-end.c
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
#include <rule.h>
#include <extfunc.h>

int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;
int	qualinp, mateinp;
char	inpfile[100], outfile[100], seqfile[100], qualfile[100], rulefile[100],
  	edgefile[100], graphfile[100], intvfile[100], matefile[100];
char	htmlout, caption[2000], ***content;
FILE	*flog;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	reads, multip;
	int	num_vertex, num_class, num_edge;
	int	*beginreads, *endreads, num_endreads[2];
	ALIGN	**eq_class, *align;
	EDGE	**edge;
	char	**score;
	char	temp[100];
	READTABLE RT;
	MATEPAIRTABLE MP;
	MATEPAIRRULES MPR;
	NODES	**vertex;
	FILE	*fp;

	readpar();
	random1(&idum);
	initenv(argc, argv);
	set_score();
	if(htmlout)	{
		flog = ckopen("EULER-report.html", "a");
	} else	{
		flog = ckopen("EULER-report.txt", "a");
	}
	if(!htmlout)	{
		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "Overlap detection between end reads.\n\n");
	} else	{
		print_section_head(flog, "Overlap detection between end reads.");
	}

/*	Input the length of the reads (required) */

	read_fasta_file(seqfile, &RT);

/*	Input the quality values (optional): 
	score = 2,       if quality >= 30
	        1,       if 30 > quality > 20
		0,	 if quality <= 20
	if no quality values are input, all scores are set up
	as 1				*/

	score = init_score(&RT);
	if(qualinp)	{
		read_score_file(qualfile, score, &RT);
	}

/*	Input equivalent readintervals between reads --
	see the format of the equivalent readinterval files	*/

	eq_class = (ALIGN **) ckalloc(2 * RT.num_seq * sizeof(ALIGN *));
	fp = ckopen(inpfile, "r");
	num_class = readclass(eq_class, RT.num_seq, fp);
	fclose(fp);
	if(!htmlout)	{
		fprintf(flog, "# equivalent readintervals input: %d\n", num_class);
	} else	{
		sprintf(temp, "# equivalent readintervals input: %d", num_class);
		print_line(flog, temp);
	}

/*	Input the graph	*/

  	read_graph_file(edgefile, graphfile, &num_vertex, &vertex, &num_edge, &edge);

/**********************************************************************
 * Input the read intervals in each edge
 **********************************************************************/

	read_interval_file(intvfile, num_vertex, vertex);

  /**********************************************************************
   * Input mate-pair information
   **********************************************************************/

	MP.num_pair = 0;
	if(mateinp)	{
		read_matepair_rules_file(rulefile, &MPR);
		read_matepair_file(matefile, &MPR, &RT, &MP);
  	}

/*	Find end reads	*/
	beginreads = (int *) ckalloc(2 * RT.num_seq * sizeof(int));
	endreads = (int *) ckalloc(2 * RT.num_seq * sizeof(int));
	findendreads(vertex, num_vertex, endreads, beginreads, num_endreads, RT.num_seq);
	if(mateinp)	{
		findmatereads(vertex, num_vertex, endreads, beginreads, num_endreads, RT.num_seq, &MP);
	}
	free_label();

/*	Detect overlaps between endreads	*/
	n = alignend(RT.src_name, RT.src_seq, score, RT.len_seq, RT.num_seq, endreads, beginreads, num_endreads, eq_class);
	if(!htmlout)	{
		fprintf(flog, "# of overlaps detected between end reads: %d\n", n);
	} else	{
		sprintf(temp, "# of overlaps detected between end reads: %d", n);
		print_line(flog, temp);
	}
	free((void *) beginreads);
	free((void *) endreads);

/*	Write alignment 	*/
	fp = ckopen(outfile, "w");
	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		n = size_align(eq_class[i]);
		fwrite(&n, sizeof(int), 1, fp);
		align = eq_class[i];
		while(align)	{
			if(align -> reads[1] >= RT.num_seq)	{
				k = align -> reads[1] - RT.num_seq;
			} else	{
				k = align -> reads[1];
			}
			fwrite(&(align -> reads[1]), sizeof(int), 1, fp);
			fwrite(&(align -> mis_match), sizeof(int), 1, fp);
			fwrite(&(align -> length), sizeof(int), 1, fp);
			fwrite(align -> pos[0], sizeof(int), align -> length, fp);
			fwrite(align -> pos[1], sizeof(int), align -> length, fp);
			align = free_align(align);
		}
	}
	fclose(fp);

	if(!htmlout)	{
  		print_text_line(flog, LINE_LENGTH);
	} else	{
		print_hl(flog);
	}
	fclose(flog);

	free_graph(vertex, num_vertex);
	free((void **) vertex);
	free((void **) edge);
	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		free((void *) score[i]);
	}
	free((void **) score);
	free_readtable(&RT);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;	
	extern char *optarg;

	strcpy(rulefile, "name.rul");
	inpseq = outseq = qualinp = mateinp = 0;
	MIN_OVERLAP = 100;
	htmlout = 0;

	while ((copt=getopt(argc,argv,"i:o:s:c:v:q:d:m:M:w:g:h:u:e:R:H")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'd':
			  sscanf(optarg,"%lf", &MIN_IDENTITY);
			  continue;
			case 'q':
			  sscanf(optarg,"%s", qualfile);
			  qualinp = 1;
			  continue;
			case 'e':
			  sscanf(optarg,"%d", &END_LEG);
			  continue;
			case 'g':
			  sscanf(optarg,"%d", &g);
			  continue;
			case 'h':
			  sscanf(optarg,"%d", &h);
			  continue;
			case 'R':
			  sscanf(optarg,"%s", rulefile);
			  continue;
			case 'w':
			  sscanf(optarg,"%d", &LPAT);
			  continue;
			case 's':
			  sscanf(optarg,"%s", seqfile);
			  continue;
			case 'v':
			  sscanf(optarg,"%d", &MIN_OVERLAP);
			  continue;
			case 'c':
			  sscanf(optarg,"%d", &LOW_COV);
			  continue;
			case 'u':
			  sscanf(optarg,"%d", &match_score);
			  continue;
			case 'm':
			  mateinp = 1;
			  sscanf(optarg,"%s", matefile);
			  continue;
			case 'M':
			  sscanf(optarg,"%d", &MIS_SCORE);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'H':
			  htmlout = 1;
			  continue;
			default:
			  if(!htmlout)	{
				  printf("over-repeat-end -i InpFile -s SeqFile -o outfile [-q QualFile -c LOW_COV -v MIN_OVERLAP -d MIN_ID -m MateFile]\n");
				  printf("-i InpFile: The input file name of alignments\n");
				  printf("-s SeqFile: The input file name of reads\n");
				  printf("-o Outfile: The output file name of contigs\n");
				  printf("-m (optional): The input file name of mate-pairs\n");
				  printf("-c LOW_COV (optional): minimal coverage\n");
				  printf("-v MIN_OVERLAP (optional): minimal overlap length\n");
			  } else	{
				  print_line(flog, "over-repeat-end -i InpFile -s SeqFile -o outfile [-q QualFile -c LOW_COV -v MIN_OVERLAP -d MIN_ID -m MateFile]");
				  print_line(flog, "-i InpFile: The input file name of alignments");
				  print_line(flog, "-s SeqFile: The input file name of reads");
				  print_line(flog, "-o Outfile: The output file name of contigs");
				  print_line(flog, "-m (optional): The input file name of mate-pairs");
				  print_line(flog, "-c LOW_COV (optional): minimal coverage");
				  print_line(flog, "-v MIN_OVERLAP (optional): minimal overlap length");
			  }
			  exit(-1);
		}
		optind--;
	}
	sprintf(edgefile, "%s.edge", seqfile);
	sprintf(graphfile, "%s.graph", seqfile);
	sprintf(intvfile, "%s.intv", seqfile);

	if(inpseq == 0 || outseq == 0)	{
		if(!htmlout)	{
			printf("over-repeat-end -i InpFile -s SeqFile -o outfile [-q QualFile -c LOW_COV -v MIN_OVERLAP -d MIN_ID -m MateFile]\n");
			printf("-i InpFile: The input file name of alignments\n");
			printf("-s SeqFile: The input file name of reads\n");
			printf("-o Outfile: The output file name of contigs\n");
			printf("-m (optional): The input file name of mate-pairs\n");
			printf("-c LOW_COV (optional): minimal coverage\n");
			printf("-v MIN_OVERLAP (optional): minimal overlap length\n");
		} else	{
			print_line(flog, "over-repeat-end -i InpFile -s SeqFile -o outfile [-q QualFile -c LOW_COV -v MIN_OVERLAP -d MIN_ID -m MateFile]");
			print_line(flog, "-i InpFile: The input file name of alignments");
			print_line(flog, "-s SeqFile: The input file name of reads");
			print_line(flog, "-o Outfile: The output file name of contigs");
			print_line(flog, "-m (optional): The input file name of mate-pairs");
			print_line(flog, "-c LOW_COV (optional): minimal coverage");
			print_line(flog, "-v MIN_OVERLAP (optional): minimal overlap length");
		}
		exit(-1);
	}
}
