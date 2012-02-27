/***************************************************************************
 * Title:          trans_qual.c
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

char	inpfile[100], outfile[100], qualfile[100], overlapfile[100];
int	qualinp;
char	htmlout, caption[2000], ***content;
FILE	*flog;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, m, n, hq, hq2;
	char	*index, c;
	char	temp[200];
	int	num_seq;
	char	**src_seq, **score, **src_name;
	int	*len_seq;
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
		print_text(flog, "Removing duplicated reads.");
	} else	{
		print_hl(flog);
	       	print_line(flog, "Removing duplicated reads");
	}

/*	Input the reads (required) */

	len_seq = (int *) ckalloc(2 * max_seq * sizeof(int));
	src_seq = (char **) ckalloc(2 * max_seq * sizeof(char *));
	src_name = (char **) ckalloc(max_seq * sizeof(char *));
	for(i = 0; i < max_seq; i ++)	{
		src_name[i] = (char *) ckalloc(100 * sizeof(char));
	}
	fp = ckopen(inpfile, "r");
	num_seq = readseq1by1gen(src_seq, src_name, len_seq, fp);
	fclose(fp);

/*	Input the quality values (optional): 
	score = 2,       if quality >= 30
	        1,       if 30 > quality > 20
		0,	 if quality <= 20
	if no quality values are input, all scores are set up
	as 1				*/

	score = (char **) ckalloc(2 * num_seq * sizeof(char *));
	for(i = 0; i < num_seq; i ++)	{
		score[i] = (char *) ckalloc(len_seq[i] * sizeof(char));
		for(j = 0; j < len_seq[i]; j ++)	{
			score[i][j] = 1;
		}
	}
	if(qualinp)	{
		fp = ckopen(qualfile, "r");
		n = readqual(score, src_seq, len_seq, fp);
		if(n != num_seq)	{
			printf("The quality file (%d sequence) and the reads (%d sequence) file are not consistent!\n", n,
				num_seq);
			exit(0);
		}
		fclose(fp);
	}

/*	Delete redundent reads	*/

	fp = ckopen("deg.list", "w");
	fprintf(fp, "redundent reads:\n");
	index = (char *) ckalloc(num_seq * sizeof(char));
	for(i = 0; i < num_seq; i ++)	{
		for(j = i + 1; j < num_seq; j ++)	{
			if(len_seq[i] == len_seq[j])	{
				c = comp_word(src_seq[i], src_seq[j], len_seq[i]);
				if(c)	{
					fprintf(fp, "%d %s %d -- %d %s %d\n", i, src_name[i], len_seq[i],
						j, src_name[j], len_seq[j]);
					index[j] = 1;
				}
			}
		}
	}

/*	Delete very low quality reads	*/

	k = hq = hq2 = 0;
	fprintf(fp, "low quality reads:\n");
	for(i = 0; i < num_seq; i ++)	{
		n = 0;
		m = l = 0;
		for(j = 0; j < len_seq[i]; j ++)	{
			if(score[i][j] == 2)	hq ++;
			if(score[i][j] >= 1)	hq2 ++;
			if(score[i][j] > 0)	{
				n ++;
				l ++;
			} else	{
				if(l > m)	m = l;
				l = 0;
			}
		}
	}
	fclose(fp);

	if(!htmlout)	{
		printf("%d low quality reads.\n", k);
		printf("%d(>=30) and %d(>=20) high quality positions.\n", hq, hq2);
	}

	n = 0;
	fp = ckopen(outfile, "w");
	sprintf(temp, "%s.score", outfile);
	fp1 = ckopen(temp, "w");
	for(i = 0; i < num_seq; i ++)	{
		if(index[i] == 1)	continue;
		fprintf(fp, ">%s\n", src_name[i]);
		for(j = 0; j < len_seq[i]; j ++)	{
			if(src_seq[i][j] >= 4)	src_seq[i][j] = ran_number(4, &idum);
			fprintf(fp, "%c", na_name[src_seq[i][j]]);
			if(j % 50 == 49)	{
				fprintf(fp, "\n");
			}
		}
		if(j % 50 != 0)	{
			fprintf(fp, "\n");
		}
		fprintf(fp1, ">%s\n", src_name[i]);
		for(j = 0; j < len_seq[i]; j ++)	{
			fprintf(fp1, "%d", score[i][j]);
			if(j % 50 == 49)	{
				fprintf(fp1, "\n");
			}
		}
		if(j % 50 != 0)	{
			fprintf(fp1, "\n");
		}
		n ++;
	}
	fclose(fp);
	fclose(fp1);
	if(!htmlout)	{
		fprintf(flog, "Reads output: %d\n", n);
	} else	{
		sprintf(caption, "Reads remained: %d\n", n);
		print_line(flog, caption);
	}
	if(!htmlout)	{
		print_text_line(flog, LINE_LENGTH);
	} else	{
		print_hl(flog);
	}
	fclose(flog);

/*	Free memory	*/

	free((void *) index);
	for(i = 0; i < max_seq; i ++)	{
		free((void *) src_name[i]);
	}
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void *) src_seq[i]);
		free((void *) score[i]);
	}
	free((void **) score);
	free((void **) src_seq);
	free((void **) src_name);
	free((void *) len_seq);
}

void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;
	extern char *optarg;

	inpseq = qualinp = outseq = 0;
	htmlout = 0;

	while ((copt=getopt(argc,argv,"i:q:o:H")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'q':
			  qualinp = 1;
			  sscanf(optarg,"%s", qualfile);
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
				  printf("trans_qual -i InpFile [-q qualfile] -o outfile\n");
				  printf("-i InpFile: The input file name of reads\n");
				  printf("-q qualfile (optional): Quality file\n");
			  } else	{
				  print_line(flog, "trans_qual -i InpFile [-q qualfile] -o outfile");
				  print_line(flog, "-i InpFile: The input file name of reads");
				  print_line(flog, "-q qualfile (optional): Quality file");
			  }
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0)	{
		if(!htmlout)	{
		  printf("trans_qual -i InpFile [-q qualfile] -o outfile\n");
		  printf("-i InpFile: The input file name of reads\n");
		  printf("-q qualfile (optional): Quality file\n");
		} else	{
		  print_line(flog, "trans_qual -i InpFile [-q qualfile] -o outfile");
		  print_line(flog, "-i InpFile: The input file name of reads");
		  print_line(flog, "-q qualfile (optional): Quality file");
		}
		exit(-1);
	}
}
