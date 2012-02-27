/***************************************************************************
 * Title:          over-repeat-end-po.c
 * Author:         Haixu Tang
 * Created:        Jun. 2003
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

int	LLR_score;
int	qualinp;
char	inpfile[100], outfile[100], seqfile[100], qualfile[100], rulefile[100],
  	edgefile[100], graphfile[100], intvfile[100], ovpfile[100];
char	htmlout, caption[2000], ***content;
FILE	*flog;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	reads, multip;
	int	num_vertex, num_class, num_edge;
	int	*beginreads, *endreads, num_endreads[2];
	int	read1, read2;
	int	seql[2], cutp[2], *sapp;
	ALIGN	**eq_class, *align;
	EDGE	**edge;
	char	**score;
	READTABLE RT;
	READOVERLAP **readoverlap, *rovp;
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
		fprintf(flog, "Transforming Phrap overlaps between end reads with LLR score %d+.\n\n", LLR_score);
	} else	{
		/* Commented out print_hl -Vagisha */
	       	/* print_hl(flog); */
	       	sprintf(caption, "Transforming Phrap overlaps between end reads with LLR score %d+", LLR_score);
	       	print_line(flog, caption);
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

	printf("Read equivalent readintervals...\n");
	eq_class = (ALIGN **) ckalloc(2 * RT.num_seq * sizeof(ALIGN *));
	fp = ckopen(inpfile, "r");
	num_class = readclass(eq_class, RT.num_seq, fp);
	fclose(fp);

/*	Input the graph	*/

  	read_graph_file(edgefile, graphfile, &num_vertex, &vertex, &num_edge, &edge);

/**********************************************************************
 * Input the read intervals in each edge
 **********************************************************************/

	read_interval_file(intvfile, num_vertex, vertex);

/*	Find end reads	*/
	beginreads = (int *) ckalloc(2 * RT.num_seq * sizeof(int));
	endreads = (int *) ckalloc(2 * RT.num_seq * sizeof(int));
	findendreads(vertex, num_vertex, endreads, beginreads, num_endreads, RT.num_seq);
	free_label();

/*	read into possible overlaps between reads	*/

	readoverlap = (READOVERLAP **) ckalloc(2 * RT.num_seq * sizeof(READOVERLAP *));
	fp = ckopen(ovpfile, "r");
	n = readphrapovp(readoverlap, LLR_score, fp);
	fclose(fp);

	sapp = (int *) ckalloc(MAX_TMP_LEG * sizeof(int));
/*	Detect overlaps between endreads	*/
	n = 0;
	for(m = 0; m < num_endreads[0]; m ++)	{
		if(beginreads[m] >= RT.num_seq)	{
			i = beginreads[m] - RT.num_seq;
		} else	{
			i = beginreads[m];
		}
		rovp = readoverlap[i];
		while(rovp)	{
			align = eq_class[i];
			while(align)	{
				if(rovp -> read2 == align -> reads[1])	{
					break;
				}
				align = align -> next;
			}
			if(!align)	{
				for(k = 0; k < num_endreads[1]; k ++)	{
					if(i == beginreads[m] && rovp -> read2 == endreads[k] ||
					   i == beginreads[m] - RT.num_seq && rovp -> read2
					     == reverse_read(endreads[k], RT.num_seq))	{
						cutp[0] = rovp -> pos1[0];
						cutp[1] = rovp -> pos2[0];
						seql[0] = rovp -> pos1[1] - rovp -> pos1[0] + 1;
						seql[1] = rovp -> pos2[1] - rovp -> pos2[0] + 1;
						read1 = rovp -> read1;
						read2 = rovp -> read2;
						k = ALIGN0(&(RT.src_seq[read1][cutp[0] - 1]),
							&(RT.src_seq[read2][cutp[1] - 1]), 
							seql[0], seql[1], -band, band, W, g, h, sapp,
							RT.len_seq[read1] + RT.len_seq[read2],
							RT.len_seq[read1] + RT.len_seq[read2]);
/*	transfer the LLR_score to align -> mis_match	*/
						eq_class[i] = new_align(cutp, seql, sapp, eq_class[i],
							read1, read2, 1, rovp -> LLR_score);
						n ++;
						break;
					}
				}
			}
			rovp = rovp -> next;
		}
	}
	if(!htmlout)	{
		fprintf(flog, "Adding # of overlaps between end reads: %d\n", n);
	} else	{
		sprintf(caption, "Adding # of overlaps between end reads: %d", n);
		print_line(flog, caption);
	}

	free((void *) sapp);
	free((void *) beginreads);
	free((void *) endreads);

/*	Write alignment 	*/
	fp = ckopen(outfile, "w");
	for(i = 0; i < RT.num_seq; i ++)	{
		n = size_align(eq_class[i]);
		fwrite(&n, sizeof(int), 1, fp);
		align = eq_class[i];
		while(readoverlap[i])	{
			readoverlap[i] = free_readoverlap(readoverlap[i]);
		}
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

/*	Free memory	*/

	free((void **) readoverlap);
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

	strcpy(ovpfile, "phrap.ovp");
	inpseq = outseq = qualinp = 0;
	LLR_score = 10;
	htmlout = 0;

	while ((copt=getopt(argc,argv,"i:o:s:q:p:r:H")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'q':
			  sscanf(optarg,"%s", qualfile);
			  qualinp = 1;
			  continue;
			case 's':
			  sscanf(optarg,"%s", seqfile);
			  continue;
			case 'p':
			  sscanf(optarg,"%s", ovpfile);
			  continue;
			case 'r':
			  sscanf(optarg,"%d", &LLR_score);
			  continue;
			case 'H':
			  htmlout = 1;
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			default:
			  if(!htmlout)	{
				  printf("over-repeat-end-po -i InpFile -s SeqFile -o outfile -p Ovpfile [-q QualFile -r LLR]\n");
				  printf("-i InpFile: The input file name of alignments\n");
				  printf("-s SeqFile: The input file name of reads\n");
				  printf("-o Outfile: The output file name of contigs\n");
				  printf("-p Ovpfile: The overlap file produced by Phrap\n");
				  printf("-r: mininal LLR score.\n");
			  } else	{
				  print_line(flog, "over-repeat-end-po -i InpFile -s SeqFile -o outfile -p Ovpfile [-q QualFile -r LLR]");
				  print_line(flog, "-i InpFile: The input file name of alignments");
				  print_line(flog, "-s SeqFile: The input file name of reads");
				  print_line(flog, "-o Outfile: The output file name of contigs");
				  print_line(flog, "-p Ovpfile: The overlap file produced by Phrap");
				  print_line(flog, "-r: mininal LLR score.");
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
			printf("over-repeat-end-po -i InpFile -s SeqFile -o outfile -p Ovpfile [-q QualFile -r LLR]\n");
			printf("-i InpFile: The input file name of alignments\n");
			printf("-s SeqFile: The input file name of reads\n");
			printf("-o Outfile: The output file name of contigs\n");
			printf("-p Ovpfile: The overlap file produced by Phrap\n");
			printf("-r: mininal LLR score.\n");
		} else	{
			print_line(flog, "over-repeat-end-po -i InpFile -s SeqFile -o outfile -p Ovpfile [-q QualFile -r LLR]");
			print_line(flog, "-i InpFile: The input file name of alignments");
			print_line(flog, "-s SeqFile: The input file name of reads");
			print_line(flog, "-o Outfile: The output file name of contigs");
			print_line(flog, "-p Ovpfile: The overlap file produced by Phrap");
			print_line(flog, "-r: mininal LLR score.");
		}
		exit(-1);
	}
}
