/***************************************************************************
 * Title:          trimseq_qual.c
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
double	cutoff;
int	qualinp;
char	htmlout, caption[2000], ***content;
FILE	*flog;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, n, m, q;
	char	temp[100];
	int	max_l, max_j, el[5], ek[5], nl[5], al[20], onl[20], nr;
	double	elf[5];
	char	**src_name;
	int	num_seq, *len_seq;
	int	startpos, endpos;
	int	length, max_length, min_length;
	int	region_from, region_to;
	int	**region;
	char 	**src_seq;
	char	**qualdist;
	int	LN;
	int	max_leg;
	double	v;
	FILE	*fp, *fp1;

	readpar();
	initenv(argc, argv);

	len_seq = (int *) ckalloc(max_seq * sizeof(int));
	src_seq = (char **) ckalloc(max_seq * sizeof(char *));
	src_name = (char **) ckalloc(max_seq * sizeof(char *));
	qualdist = (char **) ckalloc(max_seq * sizeof(char *));
	for(i = 0; i < max_seq; i ++)	{
		src_name[i] = (char *) ckalloc(100 * sizeof(char));
	}

	fp = ckopen(inpfile, "r");
	num_seq = readseq1by1(src_seq, src_name, len_seq, fp);
	fclose(fp);

	max_leg = 0;
	for(i = 0; i < num_seq; i ++)	{
		if(len_seq[i] > max_leg)	max_leg = len_seq[i];
		qualdist[i] = (char *) ckalloc((1 + len_seq[i]) * sizeof(char));
	}

	region = (int **) ckalloc(max_leg * sizeof(int *));
	for(i = 0; i < max_leg; i ++)	{
		region[i] = (int *) ckalloc(2 * sizeof(int));
	}

	LN = 20;

	nr = 0;
	for(i = 0; i < 5; i ++)		{
		nl[i] = 0;
		el[i] = 0;
		ek[i] = 0;
	}
	for(i = 0; i < 8; i ++)	{
		onl[i] = al[i] = 0;
	}

	m = l = q = 0;
	min_length = max_length = 0;
	fp = ckopen(qualfile, "r");
	n = readqualraw(qualdist, src_seq, len_seq, fp);
	if(n != num_seq)	{
		printf("The quality file and the reads file are not consistent!\n");
		exit(-1);
	}
	fclose(fp);
	fp = ckopen(outfile, "w");
	sprintf(temp, "%s.qual", outfile);
	fp1 = ckopen(temp, "w");
	for(i = 0; i < num_seq; i ++)	{
		l += len_seq[i];
		n = 0;
		for(j = LN; j < len_seq[i] - LN; j ++)	{
			v = 0;
			for(k = j - LN; k < j + LN; k ++)	{
				v += qualdist[i][k];
			}
			v /= (LN * 2);
			if(region[n][0] == 0 && v >= cutoff)	{
				region[n][0] = j;
			} else if(region[n][0] > 0)	{
				if(j == len_seq[i] - LN - 1)	{
					region[n][1] = j;
					n ++;
				} else if(v < cutoff)	{
					region[n][1] = j - 1;
					n ++;
				}
			}
		}
		max_l = 0;
		for(j = 0; j < n; j ++)	{
			if((k = region[j][1] - region[j][0] + 1) > max_l)	{
				max_l = k;
				max_j = j;
			}
		}
		if(max_l == 0)	{
			startpos = endpos = 0;
		} else	{
			startpos = region[max_j][0];;
			endpos = region[max_j][1];;
		}
		length = endpos - startpos + 1;
		if(length < 100)	{
			continue;
		}
		if(min_length == 0 || length < min_length)	{
			min_length = length;
		}
		if(max_length == 0 || length > max_length)	{
			max_length = length;
		}
		if(length <= 200)	nl[0] ++;
		else if(length <= 300)	nl[1] ++;
		else if(length <= 400)	nl[2] ++;
		else if(length <= 500)	nl[3] ++;
		else			nl[4] ++;
		nr ++;
		fprintf(fp, ">%s\n", src_name[i]);
		fprintf(fp1, ">%s\n", src_name[i]);
		for(j = startpos; j <= endpos; j ++)	{
			fprintf(fp, "%c", na_name[src_seq[i][j]]);
			fprintf(fp1, "%d ", qualdist[i][j]);
			k = j - startpos;
			if(k <= 100)	{
				al[0] += qualdist[i][k];
				onl[0] ++;
			} else if(k <= 200)	{
				al[1] += qualdist[i][k];
				onl[1] ++;
			} else if(k <= 300)	{
				al[2] += qualdist[i][k];
				onl[2] ++;
			} else if(k <= 400)	{
				al[3] += qualdist[i][k];
				onl[3] ++;
			} else if(k <= 500)	{
				al[4] += qualdist[i][k];
				onl[4] ++;
			} else if(k <= 600)	{
				al[5] += qualdist[i][k];
				onl[5] ++;
			} else if(k <= 700)	{
				al[6] += qualdist[i][k];
				onl[6] ++;
			} else	{
				al[7] += qualdist[i][k];
				onl[7] ++;
			}
			if(src_seq[i][j] > 3)		{
				if(qualdist[i][j] <= 10)	{
					el[0] ++;
				} else if(qualdist[i][j] <= 15)	{
					el[1] ++;
				} else if(qualdist[i][j] <= 20)	{
					el[2] ++;
				} else if(qualdist[i][j] <= 25)	{
					el[3] ++;
				} else	{
					el[4] ++;
				}
			} else	{
				if(qualdist[i][j] <= 10)	{
					ek[0] ++;
				} else if(qualdist[i][j] <= 15)	{
					ek[1] ++;
				} else if(qualdist[i][j] <= 20)	{
					ek[2] ++;
				} else if(qualdist[i][j] <= 25)	{
					ek[3] ++;
				} else	{
					ek[4] ++;
				}
			}
			if((j - startpos) % 60 == 59)	{
				fprintf(fp, "\n");
			}
			if((j - startpos) % 30 == 29)	{
				fprintf(fp1, "\n");
			}
			q ++;
		}
		if((j - startpos) % 60 != 0)	{
			fprintf(fp, "\n");
		}
		if((j - startpos) % 30 != 0)	{
			fprintf(fp1, "\n");
		}
		m ++;
	}
	fclose(fp);
	fclose(fp1);

	printf("------------------------------------------------------------------------------------\n");
	printf("EULER-TRIM based on the Phred quality value: \n\n");
	if(cutoff == 15)	{
		printf("Ends of reads with Phred quality value above %.2f (default) are trimmed\n", cutoff);
	} else	{
		printf("Ends of reads with error rate above %.2f are trimmed\n", cutoff);
	}
	printf("Trimmed reads are stored in \"%s\";\n", outfile);
	printf("Trimming deletes %d out of %d reads;\n", num_seq - m, num_seq);
	printf("%d out of %d characters left after trimming -- %.2f%% characters are trimmed out.\n",
		 q, l, (1 - ((double) q) / l) * 100);
	if(m == 0)	m = 1;
	printf("After trimming, read length varies from %d to %d.\nThe average read length is %d.\n", min_length, max_length, (q / m));
	printf("Phred quality value distribution of ambiguous positions (N or X):\n");
	printf("------------------------------------------------------------------------------------\n");
	printf("                    <10     10-15   16-20   21-25   >25     total\n");
	k = el[0] + el[1] + el[2] + el[3] + el[4];
	for(i = 0; i < 5; i ++)	{
		elf[i] = 100 * ((double) el[i]) / (k + 1);
	}
	printf(" # of positions     %-6d  %-6d  %-6d  %-6d  %-6d  %-6d\n",
		 el[0], el[1], el[2], el[3], el[4], k);
	printf(" Percentage      %5.1f%%  %5.1f%%  %5.1f%%  %5.1f%%  %5.1f%%\n",
		 elf[0], elf[1], elf[2], elf[3], elf[4]);
	printf("------------------------------------------------------------------------------------\n");
	printf("Phred quality value distribution of non-ambiguous positions:\n");
	printf("------------------------------------------------------------------------------------\n");
	printf("                    <10     10-15   16-20   21-25   >25     total\n");
	k = ek[0] + ek[1] + ek[2] + ek[3] + ek[4];
	for(i = 0; i < 5; i ++)	{
		elf[i] = 100 * ((double) ek[i]) / (k + 1);
	}
	printf(" # of positions     %-6d  %-6d  %-6d  %-6d  %-6d  %-6d\n",
		 ek[0], ek[1], ek[2], ek[3], ek[4], k);
	printf(" Percentage         %-5.1f%%  %-5.1f%%  %-5.1f%%  %-5.1f%%  %-5.1f%%\n",
		 elf[0], elf[1], elf[2], elf[3], elf[4]);
	printf("------------------------------------------------------------------------------\n");
	printf("Distribution of read length after trimming:\n");
	printf("------------------------------------------------------------------------------\n");
	printf("Length       100-200  201-300  301-400  401-500  >501     total\n");
	printf("# of reads   %-7d  %-7d  %-7d  %-7d  %-7d  %-7d\n", nl[0], nl[1], nl[2], nl[3], nl[4], nr);
	for(i = 0; i < 5; i ++)	{
		elf[i] = 100 * ((double) nl[i]) / nr;
	}
	printf("Percentage %6.1f%%  %6.1f%%  %6.1f%%  %6.1f%%  %7.1f%%\n",
		 elf[0], elf[1], elf[2], elf[3], elf[4]);
	printf("------------------------------------------------------------------------------\n");
	for(i = 0; i < 8; i ++)	{
		if(onl[i] > 0)	al[i] /= onl[i];
		else		al[i] = 0;
	}
	printf("Distribution of Phred quality values along the length of the reads.\n");
	printf("------------------------------------------------------------------------------\n");
	printf("Position     0-100  101-200  201-300  301-400  401-500  501-600  601-600  >700\n");
	printf("Phred values %-5d  %-7d  %-7d  %-7d  %-7d  %-7d  %-7d  %-4d\n", al[0], al[1], al[2],
		al[3], al[4], al[5], al[6], al[7]);
	printf("------------------------------------------------------------------------------\n");

	for(i = 0; i < max_leg; i ++)	{
		free((void *) region[i]);
	}
	free((void **) region);

	free((void *) len_seq);
	for(i = 0; i < max_seq; i ++)   {
		free((void *) src_name[i]);
	}
	for(i = 0; i < num_seq; i ++)   {
		free((void *) src_seq[i]);
		free((void *) qualdist[i]);
	}
	free((void **) src_seq);
	free((void **) src_name);
	free((void **) qualdist);
	return(0);
}


void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;	
	int qualinp;
	extern char *optarg;

	inpseq = outseq = 0;
	qualinp = 0;
	cutoff = 15;
	htmlout = 0;

	while ((copt=getopt(argc,argv,"i:o:c:q:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", outfile);
			  continue;
			case 'c':
			  sscanf(optarg,"%lf", &cutoff);
			  continue;
			case 'q':
			  qualinp = 1;
			  sscanf(optarg,"%s", qualfile);
			  continue;
			default:
			  printf("trimseq_qual -i SeqFile -o Outseq -c Cutoff -q QualityFile \n\n");
			  printf("-i SeqFile: The file name of the inputting reads\n");
			  printf("-o Outseq: The file name of the outputing reads\n");
			  printf("-c Cutoff: Cutoff value of trimming; it will be applied to quality-value-based trimming\n");
			  printf("-q QualityFile: the file name of quality values\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0 || qualinp == 0)	{
		printf("trimseq_qual -i SeqFile -o Outseq -c Cutoff -q QualityFile \n\n");
		printf("-i SeqFile: The file name of the inputting reads\n");
		printf("-o Outseq: The file name of the outputing reads\n");
		printf("-c Cutoff: Cutoff value of trimming; it will be applied to quality-value-based trimming\n");
		printf("-q QualityFile: the file name of quality values\n");
		exit(-1);
	}
}
