/***************************************************************************
 * Title:          errcorr_pair_mem.c
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

char	inpfile[100], outfile[100], seqfile[100], qualfile[100];
int	qualinp;
int	gap_k;
char	htmlout, caption[2000], ***content;
FILE	*flog;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n;
	int	num_vertex, num_class, num_edge;
	int	k1, k2;
	int	read, pos;
	READTABLE RT;
	int	read1, read2, pos1, pos2;
	int	read_rev1, read_rev2, pos_rev1, pos_rev2;
	char	**multip;
	int	dist[10][10];
	int	intv;
	char	**score;
	char	temp[100];
	int	n_match, n_mis, n_indel;
	int 	reads[2];
	ALIGN	**eq_class, *align, *align0;
	FILE	*fp, *fp1;

	readpar();
	random1(&idum);
	initenv(argc, argv);
	if(htmlout)	{
		flog = ckopen("EULER-report.html", "a");
	} else	{
		flog = ckopen("EULER-report.txt", "a");
	}
	if(!htmlout)	{
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "Error correction\n\n");
	} else	{
		print_section_head(flog, "Error correction.");
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

/*	Input equivalent classes between reads --
	see the format of the equivalent class files	*/

	eq_class = (ALIGN **) ckalloc(2 * RT.num_seq * sizeof(ALIGN *));
	fp = ckopen(inpfile, "r");
	num_class = readclass(eq_class, RT.num_seq, fp);
	fclose(fp);

/*	Do statistics of match/mismatch, insertion/deletion in pairwise alignments	*/

	n = n_match = n_mis = n_indel = 0;
	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			for(j = 0; j < align -> length - 1; j ++)	{
				pos1 = align -> pos[0][j];
				pos2 = align -> pos[1][j];
				while(pos1 < align -> pos[0][j + 1] && pos2 < align -> pos[1][j + 1])	{
					if(RT.src_seq[read1][pos1] == RT.src_seq[read2][pos2])	{
						n_match ++;
					} else	{
						n_mis ++;
					}
					pos1 ++;
					pos2 ++;
				}
				n_indel += (align -> pos[0][j + 1] - pos1 + align -> pos[1][j + 1] - pos2);
			}
			n ++;
			align = align -> next;
		}
	}
	if(!htmlout)	{
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "Summary of pairwise alignments\n");
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "# alignments:              %d.\n", n);
		fprintf(flog, "# matched positions:       %d.\n", n_match);
		fprintf(flog, "# mismatched positions:    %d.\n", n_mis);
		fprintf(flog, "# indels:                  %d.\n", n_indel);
		fprintf(flog, "# total positions:         %d.\n", n_indel + n_mis + n_match);
  		print_text_line(flog, LINE_LENGTH);
	} else	{
		sprintf(caption, "Summary of pairwise alignments");
		content = allocate_content(2, 5, 50);
		strcpy(content[0][0], "# alignments");
		strcpy(content[0][1], "# matched positions");
		strcpy(content[0][2], "# mismatched positions");
		strcpy(content[0][3], "# indels");
		strcpy(content[0][4], "# total positions");
		sprintf(content[1][0], "%d", n);
		sprintf(content[1][1], "%d", n_match);
		sprintf(content[1][2], "%d", n_mis);
		sprintf(content[1][3], "%d", n_indel);
		sprintf(content[1][4], "%d", n_indel + n_mis + n_match);
		print_table(flog, 2, 5, content, caption);
		content = free_content(content, 2, 5);
	}

/*	Computer multiplicity of each position in each read	*/

	multip = (char **) ckalloc(2 * RT.num_seq * sizeof(char *));
	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		multip[i] = (char *) ckalloc(RT.len_seq[i] * sizeof(char));
		for(j = 0; j < RT.len_seq[i]; j ++)	{
			multip[i][j] = 1;
		}
	}
	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			read_rev1 = reverse_read(read1, RT.num_seq);
			read_rev2 = reverse_read(read2, RT.num_seq);
			for(j = 0; j < align -> length - 1; j ++)	{
				if(j == 0)	{
					pos1 = align -> pos[0][j];
					pos2 = align -> pos[1][j];
				} else	{
					pos1 = align -> pos[0][j] + gap_k;
					pos2 = align -> pos[1][j] + gap_k;
				}
				if(j == align -> length - 2)	{
					k1 = align -> pos[0][j + 1];
					k2 = align -> pos[1][j + 1];
				} else	{
					k1 = align -> pos[0][j + 1] - gap_k;
					k2 = align -> pos[1][j + 1] - gap_k;
				}
				while(pos1 < k1 && pos2 < k2)	{
					pos_rev1 = RT.len_seq[read1] - pos1 - 1;
					pos_rev2 = RT.len_seq[read2] - pos2 - 1;
					if(RT.src_seq[read1][pos1] == RT.src_seq[read2][pos2])	{
						if(multip[read1][pos1] < 9)	multip[read1][pos1] ++;
						if(multip[read2][pos2] < 9)	multip[read2][pos2] ++;
						if(multip[read_rev1][pos_rev1] < 9)	multip[read_rev1][pos_rev1] ++;
						if(multip[read_rev2][pos_rev2] < 9)	multip[read_rev2][pos_rev2] ++;
					}
					pos1 ++;
					pos2 ++;
				}
			}
			align = align -> next;
		}
	}

/*	Simple statistics of read coverage	*/

	for(i = 0; i < 10; i ++)
		for(j = 0; j < 10; j ++)
			dist[i][j] = 0;

	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		for(j = 0; j < RT.len_seq[i]; j ++)	{
			if(multip[i][j] >= 7)	{
				dist[score[i][j]][7] ++;
			} else	{
				dist[score[i][j]][multip[i][j]] ++;
			}
		}
	}
	if(!htmlout)	{
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "Distribution of multiplicity of read overlaps.\n");
		fprintf(flog, "---Multiplicity k means one positions perfectly overlaped with k-1 other reads\n");
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "Quality  1        2        3        4        5        6        >=7\n");
		fprintf(flog, "<20      ");
		for(i = 1; i <= 7; i ++)	{
			fprintf(flog, "%-9d", dist[0][i]);
		}
		fprintf(flog, "\n");
		fprintf(flog, "20-30    ");
		for(i = 1; i <= 7; i ++)	{
			fprintf(flog, "%-9d", dist[1][i]);
		}
		fprintf(flog, "\n");
		fprintf(flog, ">30      ");
		for(i = 1; i <= 7; i ++)	{
			fprintf(flog, "%-9d", dist[2][i]);
		}
		fprintf(flog, "\n");
  		print_text_line(flog, LINE_LENGTH);
	} else	{
		sprintf(caption, "Distribution of multiplicity of read overlaps.");
		content = allocate_content(4, 8, 50);
		strcpy(content[0][0], "Quality");
		for(i = 1; i < 7; i ++)	{
			sprintf(content[0][7], "%d", i);
		}
		strcpy(content[0][7], ">6");
		strcpy(content[1][0], "<20");
		strcpy(content[2][0], "20-30");
		strcpy(content[3][0], ">30");
		for(i = 1; i < 4; i ++)	{
			for(j = 1; j < 8; j ++)	{
				sprintf(content[i][j], "%d", dist[i - 1][j]);
			}
		}
		print_table(flog, 4, 8, content, caption);
		content = free_content(content, 4, 8);
	}

/*	Fix errors	*/

	fp = ckopen(outfile, "w");
	sprintf(temp, "%s.score", outfile);
	intv = gap_k;
	fp1 = ckopen(temp, "w");
	errcorrt_pair_mem(eq_class, multip, RT.len_seq, RT.src_seq, score, RT.src_name, RT.num_seq, LOW_COV, intv, fp, fp1);
	fclose(fp);
	fclose(fp1);

	if(!htmlout)	{
  		print_text_line(flog, LINE_LENGTH);
	} else	{
		print_hl(flog);
	}
	fclose(flog);

/*	Free memory	*/

	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		while(eq_class[i])	{
			eq_class[i] = free_align(eq_class[i]);
		}
	}
	free((void **) eq_class);

	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		free((void *) multip[i]);
	}
	free((void **) multip);

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

	inpseq = qualinp = outseq = 0;
	gap_k = 5;
	END_MERGE = 0;
	mov_qual = 0;
	LOW_COV = 1;
	htmlout = 0;

	while ((copt=getopt(argc,argv,"i:s:c:o:q:k:emH")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 's':
			  sscanf(optarg,"%s", seqfile);
			  continue;
			case 'c':
			  sscanf(optarg,"%d", &LOW_COV);
			  continue;
			case 'k':
			  sscanf(optarg,"%d", &gap_k);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'm':
			  mov_qual = 1;
			  continue;
			case 'e':
			  END_MERGE = 1;
			  continue;
			case 'q':
			  qualinp = 1;
			  sscanf(optarg,"%s", qualfile);
			  continue;
			case 'H':
			  htmlout = 1;
			  continue;
			default:
			  if(!htmlout)	{
				  printf("errcorr_pair_mem -i InpFile -s SeqFile -o outfile -k gap_k [-c LOW_COV -q qualinp]\n");
				  printf("-i InpFile: The input file name of alignments\n");
				  printf("-s SeqFile: The input file name of reads\n");
				  printf("-o OutFile: The output file name of reads\n");
				  printf("-k gap_k: Distance from the gap\n");
				  printf("-c LOW_COV(optional): low coverage\n");
				  printf("-q qualinp(optional): The input quality files\n");
			  } else	{
				  print_line(flog, "errcorr_pair_mem -i InpFile -s SeqFile -o outfile -k gap_k [-c LOW_COV -q qualinp]");
				  print_line(flog, "-i InpFile: The input file name of alignments");
				  print_line(flog, "-s SeqFile: The input file name of reads");
				  print_line(flog, "-o OutFile: The output file name of reads");
				  print_line(flog, "-k gap_k: Distance from the gap");
				  print_line(flog, "-c LOW_COV(optional): low coverage");
				  print_line(flog, "-q qualinp(optional): The input quality files");
			  }
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		if(!htmlout)	{
			printf("errcorr_pair_mem -i InpFile -s SeqFile -o outfile -k gap_k [-c LOW_COV -q qualinp]\n");
			printf("-i InpFile: The input file name of alignments\n");
			printf("-s SeqFile: The input file name of reads\n");
			printf("-o OutFile: The output file name of reads\n");
			printf("-k gap_k: Distance from the gap\n");
			printf("-c LOW_COV(optional): low coverage\n");
			printf("-q qualinp(optional): The input quality files\n");
		} else	{
			print_line(flog, "errcorr_pair_mem -i InpFile -s SeqFile -o outfile -k gap_k [-c LOW_COV -q qualinp]");
			print_line(flog, "-i InpFile: The input file name of alignments");
			print_line(flog, "-s SeqFile: The input file name of reads");
			print_line(flog, "-o OutFile: The output file name of reads");
			print_line(flog, "-k gap_k: Distance from the gap");
			print_line(flog, "-c LOW_COV(optional): low coverage");
			print_line(flog, "-q qualinp(optional): The input quality files");
		}
		exit(-1);
	}
}
