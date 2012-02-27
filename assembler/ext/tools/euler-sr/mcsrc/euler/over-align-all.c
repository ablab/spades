/***************************************************************************
 * Title:          over-align-all.c
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

char	inpfile[100], outfile[100], qualfile[100], overlapfile[100], outseqfile[100];
int	qualinp;
char	htmlout, caption[2000], ***content;
FILE	*flog;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, nconf, nk;
	int	num_vertex, num_edge;
	int	num_comp;
	READTABLE RT;
	char	**score;
	int	dist[10];
	char	temp[100];
	READINTERVAL	**readinterval;
	ALIGN	**align, *aln, *aln0;
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
		fprintf(flog, "Verify overlaps from other program for overlap detection; Mininum identity: %2f, Minimum overlap length %d\n",
			 MIN_IDENTITY, MIN_OVERLAP);
		if(qualinp)	{
			fprintf(flog, "Quality file considered.\n");
		} else	{
			fprintf(flog, "Quality file not considered.\n");
		}
	} else	{
	       	print_hl(flog);
		sprintf(caption, "Verify overlaps from other program for overlap detection; Mininum identity: %2f, Minimum overlap length %d\n",
			 MIN_IDENTITY, MIN_OVERLAP);
	       	print_line(flog, caption);
		if(qualinp)	{
			print_line(flog, "Quality file considered.\n");
		} else	{
			print_line(flog, "Quality file not considered.\n");
		}
	}

/*	Input the reads (required) */

	read_fasta_file(inpfile, &RT);

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

/*	read into possible overlaps between reads	*/

	readinterval = (READINTERVAL **) ckalloc(2 * RT.num_seq * sizeof(READINTERVAL *));
	fp = ckopen(overlapfile, "r");
	n = readoverlap(readinterval, fp);
	fclose(fp);
	if(!htmlout)	{
		fprintf(flog, "# overlaps input: %d on %d reads.\n", n, RT.num_seq);
	} else	{
		sprintf(caption, "# overlaps input: %d on %d reads.", n, RT.num_seq);
		print_line(flog, caption);
	}

/*	Generate alignment from overlaps	*/

	align = (ALIGN **) ckalloc(2 * RT.num_seq * sizeof(ALIGN *));
	nconf = 0;
	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		align[i] = AlignReads(RT.src_seq, score, RT.len_seq, RT.num_seq, readinterval[i], i);
		while(readinterval[i])	{
			readinterval[i] = free_readinterval(readinterval[i]);
		}

/*	Determine the beginning and ending positions of a particular read covered by overlaps	*/

		aln = align[i];
		while(aln)	{
			aln = aln -> next;
			nconf ++;
		}
	}
	free((void **) readinterval);
	if(!htmlout)	{
		fprintf(flog, "Input %d overlap candidates, %d verified by alignment.\n", n, nconf);
	} else	{
		sprintf(caption, "Input %d overlap candidates, %d verified by alignment.\n", n, nconf);
		print_line(flog, caption);
	}

	for(i = 0; i < 7; i ++)	{
		dist[i] = 0;
	}
	fp = ckopen(outfile, "w");
	for(i = 0; i < 2 * RT.num_seq; i ++)	{
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
	for(i = 0; i < 2 * RT.num_seq; i ++)	{
		free((void *) score[i]);
	}
	free((void **) score);
	free_readtable(&RT);
	return(num_comp);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	inpseq = qualinp = outseq = 0;
	htmlout = 0;

	while ((copt=getopt(argc,argv,"i:p:q:o:d:u:l:s:H")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'p':
			  sscanf(optarg,"%s", overlapfile);
			  continue;
			case 'l':
			  sscanf(optarg,"%d", &MIN_OVERLAP);
			  continue;
			case 'u':
			  sscanf(optarg,"%d", &match_score);
			  continue;
			case 'd':
			  sscanf(optarg,"%lf", &MIN_IDENTITY);
			  continue;
			case 'q':
			  qualinp = 1;
			  sscanf(optarg,"%s", qualfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 's':
			  sscanf(optarg,"%s", outseqfile);
			  continue;
			case 'H':
			  htmlout = 1;
			  continue;
			default:
			  if(!htmlout)	{
				  printf("over-align-all -i InpFile -p Overlapfile [-q qualfile -d MIN_IDENTITY -l MIN_OVERLAP] -o outfile\n");
				  printf("-i InpFile: The input file name of reads\n");
				  printf("-p OverlapFile: The input overlap file name of reads\n");
				  printf("-d MIN_IDENTITY(optional): minimum identity for overlaps\n");
				  printf("-l MIN_OVERLAP(optional): minimum length for overlaps\n");
				  printf("-q qualfile (optional): Quality file\n");
				  printf("-o outfile: Output overlap file\n");
			  } else	{
				  print_line(flog, "over-align-all -i InpFile -p Overlapfile [-q qualfile -d MIN_IDENTITY -l MIN_OVERLAP] -o outfile");
				  print_line(flog, "-i InpFile: The input file name of reads");
				  print_line(flog, "-p OverlapFile: The input overlap file name of reads");
				  print_line(flog, "-d MIN_IDENTITY(optional): minimum identity for overlaps");
				  print_line(flog, "-l MIN_OVERLAP(optional): minimum length for overlaps");
				  print_line(flog, "-q qualfile (optional): Quality file");
				  print_line(flog, "-o outfile: Output overlap file");
			  }
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		if(!htmlout)	{
			printf("over-align-all -i InpFile -p Overlapfile [-q qualfile -d MIN_IDENTITY -l MIN_OVERLAP] -o outfile\n");
			printf("-i InpFile: The input file name of reads\n");
			printf("-p OverlapFile: The input overlap file name of reads\n");
			printf("-d MIN_IDENTITY(optional): minimum identity for overlaps\n");
			printf("-l MIN_OVERLAP(optional): minimum length for overlaps\n");
			printf("-q qualfile (optional): Quality file\n");
			printf("-o outfile: Output overlap file\n");
		} else	{
			print_line(flog, "over-align-all -i InpFile -p Overlapfile [-q qualfile -d MIN_IDENTITY -l MIN_OVERLAP] -o outfile");
			print_line(flog, "-i InpFile: The input file name of reads");
			print_line(flog, "-p OverlapFile: The input overlap file name of reads");
			print_line(flog, "-d MIN_IDENTITY(optional): minimum identity for overlaps");
			print_line(flog, "-l MIN_OVERLAP(optional): minimum length for overlaps");
			print_line(flog, "-q qualfile (optional): Quality file");
			print_line(flog, "-o outfile: Output overlap file");
		}
		exit(-1);
	}
}
