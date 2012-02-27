/***************************************************************************
 * Title:          makealn.c
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

#define MAX_LEG 2500
#define LINELEN	50

char	seqfile[100], contigfile[100], alnfile[100], intvfile[100], outfile[100], outintvfile[100],
	newcontigfile[100], newcontigqualfile[100];
char    htmlout, caption[2000], ***content;

void initenv(int argc, char **argv);
void read_one_inteval(FILE *fp, READINTERVAL *readinterval, int num_intv);

main(int argc, char **argv)
{
  int i, m, l, j, k, n;
  int max_leg, max_name_leg;
  int ii, pbegin, pend, id, ll, kk;

  char temp[400], str[200];
  char *label;

  READTABLE RT_mem, *RT=&RT_mem;
  char	**contigseq, **contigname, **readseq, c;
  int   *len_contigs, num_contigs, *len_readseq;
  int   *num_intv, *startpos;
  char  **newcontigseq, **newcontigqual, *tmpseq;
  int   *len_newcontig;
  int	**index;

  READINTERVAL **readinterval;

  FILE *fp, *fp1, *fp2;

  /**********************************************************************
   * Get inputs and parameters
   **********************************************************************/

  readpar();
  //  random1(&idum);
  initenv(argc, argv);

  /**********************************************************************
   * Input the reads, their lengths, and names
   **********************************************************************/

  read_fasta_file(seqfile, RT);

  /**********************************************************************
   * Input contigs
   **********************************************************************/
  fp = ckopen(contigfile, "r");
  num_contigs = readseqpar(&max_leg, &max_name_leg, fp);
  fclose(fp);
  contigseq = (char **) ckalloc(2 * RT -> num_seq * sizeof(char *));
  len_contigs = (int *) ckalloc(2 * RT -> num_seq * sizeof(int));
  contigname = (char **) ckalloc(RT -> num_seq * sizeof(char *));
  for (i = 0; i < num_contigs; i ++) {
    contigname[i] = (char *) ckalloc((max_name_leg + 1) * sizeof(char));
  }
  fp = ckopen(contigfile, "r");
  num_contigs = readseq1by1(contigseq, contigname, len_contigs, fp);
  fclose(fp);
  printf("# contigs input: %d\n", num_contigs);

  len_newcontig = (int *) ckalloc(num_contigs * sizeof(int));
  newcontigseq = (char **) ckalloc(num_contigs * sizeof(char *));
  newcontigqual = (char **) ckalloc(num_contigs * sizeof(char *));
  for(i = 0; i < num_contigs; i ++)	{
	newcontigseq[i] = (char *) ckalloc((len_contigs[i] + MAX_LEG) * sizeof(char));
	newcontigqual[i] = (char *) ckalloc((len_contigs[i] + MAX_LEG) * sizeof(char));
  }

  /**********************************************************************
   * Input the multiple alignment of the reads in each contig
   **********************************************************************/
   readinterval = (READINTERVAL **) ckalloc(num_contigs * sizeof(READINTERVAL *));
   index = (int **) ckalloc(num_contigs * sizeof(int *));
   num_intv = (int *) ckalloc(num_contigs * sizeof(int));
 
   n = 0;
   fp = ckopen(intvfile, "r");
   fp1 = ckopen(alnfile, "r");
   fp2 = ckopen(outfile, "w");
   while(fgets(str, 190, fp1))	{
	if(!strncmp(str, "%%%", 3))	continue;
	if(!strncmp(str, "% number", 8))	{
/*	Skip a line with %%%%%%%%%%%%%%%%%%%% */
		sscanf(str, "%*s%*s%*s%*s%*s%*s%d", &num_intv[n]);
   		fgets(str, 190, fp1);
		startpos = (int *) ckalloc(num_intv[n] * sizeof(int));
		len_readseq = (int *) ckalloc((100 + num_intv[n]) * sizeof(int));
		readseq = (char **) ckalloc((100 + num_intv[n]) * sizeof(char *));
		index[n] = (int *) ckalloc(num_intv[n] * sizeof(int));
		label = (char *) ckalloc(num_intv[n] * sizeof(char));
		for(i = 0; i < num_intv[n] + 100; i ++)	{
			readseq[i] = (char *) ckalloc(MAX_LEG * sizeof(char));
		}
		m = readalnseq(fp1, readseq, len_readseq, startpos, num_intv[n], len_contigs[n] + MAX_LEG,
			newcontigseq[n], newcontigqual[n], &len_newcontig[n]);
		readinterval[n] = (READINTERVAL *) ckalloc(num_intv[n] * sizeof(READINTERVAL));
		read_one_inteval(fp, readinterval[n], num_intv[n]);
		for(i = 0; i < num_intv[n]; i ++)	{
			index[n][i] = findri(RT, readinterval[n], num_intv[n], readseq[i], len_readseq[i], label);
		}
		m = len_newcontig[n] / LINELEN + 1;
		for(j = 0; j < m; j ++)	{
			if(j == m - 1)	{
				l = len_newcontig[n];
			} else	{
				l = (j + 1) * LINELEN;
			}
/*	Output the contig sequence	*/
			/*fprintf(fp2, "Contig%-9d ", n + 1);*/
			fprintf(fp2, "Contig%-9d ", n);
			for(k = j * LINELEN; k < l; k ++)	{
				if(newcontigseq[n][k] == CODE_GAP)	{
					fprintf(fp2, "-");
				} else	{
					fprintf(fp2, "%c", na_name[newcontigseq[n][k]] - 'a' + 'A');
				}
			}
			fprintf(fp2, "\n");
/*	Output the read alignment	*/
			for(i = 0; i < num_intv[n]; i ++)	{
				if(l - 1 < startpos[i] || j * LINELEN >= len_readseq[i] + startpos[i])	{
					continue;
				}
				k = index[n][i];
				if(readinterval[n][k].eq_read < RT -> num_seq)	{
					fprintf(fp2, "%-14sF ", RT -> src_name[readinterval[n][k].eq_read]);
				} else	{
					id = readinterval[n][k].eq_read - RT -> num_seq;
					fprintf(fp2, "%-14sR ", RT -> src_name[id]);
				}
				for(k = j * LINELEN; k < l; k ++)	{
					if(k >= len_readseq[i] + startpos[i])	{
						break;
					} else if(k < startpos[i])	{
						fprintf(fp2, " ");
					} else if(readseq[i][k - startpos[i]] == CODE_GAP)	{
						fprintf(fp2, "-");
					} else	{
						if(readseq[i][k - startpos[i]] == newcontigseq[n][k])	{
							fprintf(fp2, "%c", na_name[readseq[i][k - startpos[i]]] - 'a' + 'A');
						} else	{
							fprintf(fp2, "%c", na_name[readseq[i][k - startpos[i]]]);
						}
					}
				}
				fprintf(fp2, "\n");
			}
			fprintf(fp2, "\n");
		}

		for(i = 0; i < num_intv[n]; i ++)	{
			free((void *) readseq[i]);
		}
		free((void *) label);
		free((void **) readseq);
		free((void *) len_readseq);
		free((void *) startpos);
		n ++;
	}
   }
   fclose(fp);
   fclose(fp1);
   fclose(fp2);
   if(n != num_contigs)	{
	printf("# contigs are not consistent: in contig file %d; in alignment file %d.\n", num_contigs, n);
	exit(-1);
   }

  /**********************************************************************
   * free memory
   **********************************************************************/

  free((void *) len_newcontig);
  free((void *) len_contigs);
  for(i = 0; i < num_contigs; i ++)	{
  	free((void *) newcontigseq[i]);
  	free((void *) newcontigqual[i]);
  	free((void *) contigseq[i]);
  	free((void *) contigname[i]);
  	free((void *) readinterval[i]);
	free((void *) index[i]);
  }
  free((void **) readinterval);
  free((void **) newcontigseq);
  free((void **) newcontigqual);
  free((void **) contigseq);
  free((void **) contigname);
  free((void *) num_intv);
  free((void **) index);
  free_readtable(RT);

  return(0);
}

void read_one_inteval(FILE *fp, READINTERVAL *readinterval, int num_intv)
{
	int	i, j, k, l;
	char	str[200];

	fgets(str, 190, fp);
	if(strncmp(str, "Contig", 6))	{
		printf("Inteval file format error!\n");
		exit(-1);
	}
	for(i = 0; i < num_intv; i ++)	{
		fgets(str, 190, fp);
		if(!strncmp(str, "Contig", 6))	{
			printf("Number of intevals are not equal: %d %d\n", i, num_intv);
			exit(0);
		}
		sscanf(str, "%*s%d%d%d%d", &(readinterval[i].eq_read), &(readinterval[i].begin),
			&(readinterval[i].length), &(readinterval[i].offset));
	}
}

void initenv(int argc, char **argv)
{
  int copt;
  int inpseq, outseq;
  extern char *optarg;

  inpseq = outseq = 0;
  MIN_OVERLAP = 100;
  htmlout = 0;

  sprintf(newcontigfile, "contig.new");
  sprintf(newcontigqualfile, "contig.qual.new");
  sprintf(intvfile, "%s.contig.ace", seqfile);
  sprintf(contigfile, "%s.contig", seqfile);

  while ((copt=getopt(argc,argv,"s:i:c:a:o:v:")) != EOF) {
    switch(copt) {
    case 's':
      inpseq = 1;
      sscanf(optarg,"%s", seqfile);
      continue;
    case 'o':
      outseq = 1;
      sscanf(optarg,"%s", outfile);
      continue;
    case 'i':
      sscanf(optarg,"%s", intvfile);
      continue;
    case 'a':
      sscanf(optarg,"%s", alnfile);
      continue;
    case 'c':
      sscanf(optarg,"%s", contigfile);
      continue;
    default:
      printf("makealn -s SeqFile -i IntvFile -a AlignFile -c ContigFle -o OutFile\n");
      printf("-s SeqFile: The input file name of reads\n");
      printf("-i IntvFile: The input file name of intervals\n");
      printf("-a AlignFile: The input file name of multiple alignments\n");
      printf("-c ContigFile: The output file name of contig sequences\n");
      printf("-o OutFile: The output file name modified reads\n");
      exit(-1);
    }
    optind--;
  }
  if (inpseq == 0 || outseq == 0) {
      printf("makeintv -s SeqFile -i IntvFile -a AlignFile -c ContigFle -o OutFile\n");
      printf("-s SeqFile: The input file name of reads\n");
      printf("-i IntvFile: The input file name of intervals\n");
      printf("-a AlignFile: The input file name of multiple alignments\n");
      printf("-c ContigFile: The output file name of contig sequences\n");
      printf("-o OutFile: The output file name modified reads\n");
      exit(-1);
  }
}
