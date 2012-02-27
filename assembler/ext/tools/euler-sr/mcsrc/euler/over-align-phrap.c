/***************************************************************************
 * Title:          over-align-phrap.c
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
#include <extfunc.h>

char	inpfile[100], outfile[100], overlapfile[100];
int	qualinp;
int 	llrthresh;
char	htmlout, caption[2000], ***content;
FILE	*flog;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, nconf, nk;
	int	pos1, pos2;
	READTABLE RT;
	int	num;
	int	dist[10];
	int	*minc, *maxc;
	char	**score;
	char	temp[100];
	READOVERLAP	**readoverlap;
	ALIGN	**align, *aln, *aln0;
	FILE	*fp, *fp1, *fp2;

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
		fprintf(flog, "Transforming overlap alignment produced by Phrap.\n\n");
	} else	{
		print_section_head(flog, "Transforming overlap alignment produced by Phrap.");
	}

/*	Input the reads (required) */

	read_fasta_file(inpfile, &RT);

/*	read into possible overlaps between reads	*/

	readoverlap = (READOVERLAP **) ckalloc(2 * RT.num_seq * sizeof(READOVERLAP *));
	fp = ckopen(overlapfile, "r");
	n = readphrapovp(readoverlap, llrthresh, fp);
	fclose(fp);
	if(!htmlout)	{
		fprintf(flog, "# Phrap overlaps with LLR score %d+: %d on %d reads.\n", llrthresh, n, RT.num_seq);
	} else	{
		sprintf(caption, "# Phrap overlaps with LLR score %d+: %d on %d reads.\n", llrthresh, n, RT.num_seq);
		print_line(flog, caption);
	}

/*	Generate alignment from overlaps	*/

	nconf = 0;
	align = (ALIGN **) ckalloc(2 * RT.num_seq * sizeof(ALIGN *));
	for(i = 0; i < RT.num_seq; i ++)	{
		align[i] = AlignReadsPO(RT.src_seq, RT.len_seq, RT.num_seq, readoverlap[i]);
		while(readoverlap[i])	{
			readoverlap[i] = free_readoverlap(readoverlap[i]);
		}
		aln = align[i];
		while(aln)	{
			nconf ++;
			aln = aln -> next;
		}
	}
	free((void **) readoverlap);
	
	for(i = 0; i < 7; i ++)	{
		dist[i] = 0;
	}
	num = 0;
	fp = ckopen(outfile, "w");
	for(i = 0; i < RT.num_seq; i ++)	{
		n = size_align(align[i]);
		fwrite(&n, sizeof(int), 1, fp);
		aln = align[i];
		while(aln)	{
			if(aln -> reads[1] >= RT.num_seq)	{
				k = aln -> reads[1] - RT.num_seq;
			} else	{
				k = aln -> reads[1];
			}
			fwrite(&(aln -> reads[1]), sizeof(int), 1, fp);
			fwrite(&(aln -> mis_match), sizeof(int), 1, fp);
			fwrite(&(aln -> length), sizeof(int), 1, fp);
			fwrite(aln -> pos[0], sizeof(int), aln -> length, fp);
			fwrite(aln -> pos[1], sizeof(int), aln -> length, fp);
			if(aln -> length < 6)	{
				dist[aln -> length] ++;
			} else	{
				dist[6] ++;
			}
			if(aln -> length == 2)	{
				l = 0;
				pos1 = aln -> pos[0][0];
				pos2 = aln -> pos[1][0];
				while(pos1 < aln -> pos[0][1] && pos2 < aln -> pos[1][1])	{
					if(RT.src_seq[aln -> reads[0]][pos1] != RT.src_seq[aln -> reads[1]][pos2])	{
						l = 1;
						break;
					}
					pos1 ++;
					pos2 ++;
				}
				if(l == 0)	{
					num ++;
				}
			}
			aln = free_align(aln);
		}
	}
	fclose(fp);
	print_gaps(flog, dist);

	if(!htmlout)	{
		print_text_line(flog, LINE_LENGTH);
	} else	{
		print_hl(flog);
	}
	fclose(flog);

/*	Free memory	*/
	free((void **) align);
	free_readtable(&RT);

	return 1;
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	llrthresh = 0;
	inpseq = outseq = 0;
	htmlout = 0;

	while ((copt=getopt(argc,argv,"i:p:r:o:u:H")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'p':
			  sscanf(optarg,"%s", overlapfile);
			  continue;
			case 'r':
			  sscanf(optarg,"%d", &llrthresh);
			  continue;
			case 'u':
			  sscanf(optarg,"%d", &match_score);
			  continue;
			case 'H':
			  htmlout = 1;
			  continue;
			default:
			  if(!htmlout)	{
				  printf("over-align-phrap -i InpFile -p Overlapfile -o OutFile [-r cut-off for LLR-score]\n");
				  printf("-i InpFile: The input file name of reads\n");
				  printf("-o OutFile: The out file name of overlaps\n");
				  printf("-p OverlapFile: The input overlap file name of reads\n");
				  printf("-r cut-off (optional): cut-off for LLR-score\n");
			  } else	{
				  print_line(flog, "over-align-phrap -i InpFile -p Overlapfile -o OutFile [-r cut-off for LLR-score]");
				  print_line(flog, "-i InpFile: The input file name of reads");
				  print_line(flog, "-o OutFile: The out file name of overlaps");
				  print_line(flog, "-p OverlapFile: The input overlap file name of reads");
				  print_line(flog, "-r cut-off (optional): cut-off for LLR-score");
			  }
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		if(!htmlout)	{
			printf("over-align-phrap -i InpFile -p Overlapfile -o OutFile [-r cut-off for LLR-score]\n");
			printf("-i InpFile: The input file name of reads\n");
			printf("-o OutFile: The out file name of overlaps\n");
			printf("-p OverlapFile: The input overlap file name of reads\n");
			printf("-r cut-off (optional): cut-off for LLR-score\n");
		} else	{
			print_line(flog, "over-align-phrap -i InpFile -p Overlapfile -o OutFile [-r cut-off for LLR-score]");
			print_line(flog, "-i InpFile: The input file name of reads");
			print_line(flog, "-o OutFile: The out file name of overlaps");
			print_line(flog, "-p OverlapFile: The input overlap file name of reads");
			print_line(flog, "-r cut-off (optional): cut-off for LLR-score");
		}
		exit(-1);
	}
}
