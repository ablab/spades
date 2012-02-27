/***************************************************************************
 * Title:          makeintv.c
 * Author:         Haixu Tang
 * Created:        May. 2004
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

char	seqfile[100], contigfile[100], alnfile[100], intvfile[100], outfile[100], outintvfile[100],
	newcontigfile[100], newcontigqualfile[100];
char    htmlout, caption[2000], ***content;

void initenv(int argc, char **argv);
void read_one_inteval(FILE *fp, READINTERVAL *readinterval, int num_intv);

main(int argc, char **argv)
{
  int i, m, l, j, k, n;
  int max_leg, max_name_leg;

  char temp[400], str[200];
  char *label;

  READTABLE RT_mem, *RT=&RT_mem;
  char	**contigseq, **contigname, **newseq, **readseq;
  int   *len_contigs, num_contigs, *len_newseq, *len_readseq;
  int   *num_intv, *startpos;
  char  **newcontigseq, **newcontigqual;
  int   *len_newcontig;

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

  len_newseq = (int *) ckalloc(RT -> num_seq * sizeof(int));
  newseq = (char **) ckalloc(RT -> num_seq * sizeof(char *));
  for(i = 0; i < RT -> num_seq; i ++)	{
  	newseq[i] = (char *) ckalloc(RT -> len_seq[i] * sizeof(char));
  	len_newseq[i] = RT -> len_seq[i];
  }

  
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
   num_intv = (int *) ckalloc(num_contigs * sizeof(int));
 
   n = 0;
   fp = ckopen(intvfile, "r");
   fp1 = ckopen(alnfile, "r");
   while(fgets(str, 190, fp1))	{
	if(!strncmp(str, "%%%", 3))	continue;
	if(!strncmp(str, "% number", 8))	{
/*	Skip a line with %%%%%%%%%%%%%%%%%%%% */
		sscanf(str, "%*s%*s%*s%*s%*s%*s%d", &num_intv[n]);
   		fgets(str, 190, fp1);
		startpos = (int *) ckalloc(num_intv[n] * sizeof(int));
		len_readseq = (int *) ckalloc(num_intv[n] * sizeof(int));
		readseq = (char **) ckalloc(num_intv[n] * sizeof(char *));
		label = (char *) ckalloc(num_intv[n] * sizeof(char));
		for(i = 0; i < num_intv[n]; i ++)	{
			readseq[i] = (char *) ckalloc(MAX_LEG * sizeof(char));
		}
		m = readalnseq(fp1, readseq, len_readseq, startpos, num_intv[n], len_contigs[n] + MAX_LEG,
			newcontigseq[n], newcontigqual[n], &len_newcontig[n]);
		readinterval[n] = (READINTERVAL *) ckalloc(num_intv[n] * sizeof(READINTERVAL));
		read_one_inteval(fp, readinterval[n], num_intv[n]);
		for(i = 0; i < num_intv[n]; i ++)	{
			chg_reads(RT, readinterval[n], num_intv[n], startpos[i], readseq[i], len_readseq[i],
				 newseq, len_newseq, i, label);
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
   if(n != num_contigs)	{
	printf("# contigs are not consistent: in contig file %d; in alignment file %d.\n", num_contigs, n);
	exit(-1);
   }
printf("num_contigs %d\n", num_contigs);

  /**********************************************************************
   * Output temporary inteval file
   **********************************************************************/

/*
for(i = 0; i < RT -> num_seq; i ++)	{
	printf("i %d %d %d\n", i, RT -> len_seq[i], len_newseq[i]);
	getchar();
}
*/

   fp = ckopen(outintvfile, "w");
   for(i = 0; i < num_contigs; i ++)	{
	fprintf(fp, "Contig%d length %d # of reads %d\n", i + 1, len_contigs[i], num_intv[i]);
	fprintf(fp, "Index    Read-name                direction length   Position(contig)  Postion(read)\n");
	for(j = 0; j < num_intv[i]; j ++)	{
		chg_intv(&readinterval[i][j], newseq, len_newseq, RT);
	}
	for(j = 0; j < num_intv[i]; j ++)	{
		if(readinterval[i][j].eq_read < RT -> num_seq)	{
			fprintf(fp, "%-9d%-25s      for %-8d %-8d-%-8d %-8d-%-8d\n",
				readinterval[i][j].eq_read + 1, RT -> src_name[readinterval[i][j].eq_read],
				len_newseq[readinterval[i][j].eq_read], readinterval[i][j].offset + 1,
				readinterval[i][j].offset + readinterval[i][j].length,
				readinterval[i][j].begin + 1, readinterval[i][j].begin + readinterval[i][j].length);
		} else	{
			fprintf(fp, "%-9d%-25s      rev %-8d %-8d-%-8d %-8d-%-8d\n",
				readinterval[i][j].eq_read - RT -> num_seq + 1,
				RT -> src_name[readinterval[i][j].eq_read - RT -> num_seq],
				len_newseq[readinterval[i][j].eq_read - RT -> num_seq], readinterval[i][j].offset + 1,
				readinterval[i][j].offset + readinterval[i][j].length,
				readinterval[i][j].begin + 1, readinterval[i][j].begin + readinterval[i][j].length);
		}
	}
   }
   fclose(fp);

  /**********************************************************************
   * Output revised sequence
   **********************************************************************/

   fp = ckopen(outfile, "w");
   for(i = 0; i < RT -> num_seq; i ++)	{
	k = 0;
	fprintf(fp, ">%s\n", RT -> src_name[i]);
	for(j = 0; j < RT -> len_seq[i]; j ++)	{
		fprintf(fp, "%c", na_name[RT -> src_seq[i][j]]);
		if(k % 60 == 59)	{
			fprintf(fp, "\n");
		}
		k ++;
		for(m = 0; m < newseq[i][j]; m ++)	{
			fprintf(fp, "-");
			if(k % 60 == 59)	{
				fprintf(fp, "\n");
			}
			k ++;
		}
	}
	if(k % 60 != 0)	{
		fprintf(fp, "\n");
	}
   }
   fclose(fp);

  /**********************************************************************
   * output new contigs and quality files
   **********************************************************************/

   fp = ckopen(newcontigfile, "w");
   fp1 = ckopen(newcontigqualfile, "w");
   for(i = 0; i < num_contigs; i ++)	{
	k = 0;
	fprintf(fp, ">Contig%d %d\n", i + 1, len_newcontig[i]);
	fprintf(fp1, ">Contig%d\n", i + 1);
	for(j = 0; j < len_newcontig[i]; j ++)	{
		if(newcontigseq[i][j] == CODE_GAP)	{
			fprintf(fp, "*");
		} else	{
			fprintf(fp, "%c", na_name[newcontigseq[i][j]]);
		}
		if(newcontigseq[i][j] != CODE_GAP)	{
			fprintf(fp1, "%d ", newcontigqual[i][j]);
		}
		if(k % 60 == 59)	{
			fprintf(fp, "\n");
		}
		if(k % 20 == 19)	{
			fprintf(fp1, "\n");
		}
		k ++;
	}
	if(k % 60 != 0)	{
		fprintf(fp, "\n");
	}
	if(k % 20 != 0)	{
		fprintf(fp1, "\n");
	}
   }
   fclose(fp);
   fclose(fp1);

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
  }
  free((void **) readinterval);
  free((void **) newcontigseq);
  free((void **) newcontigqual);
  free((void **) contigseq);
  free((void **) contigname);
  for(i = 0; i < RT -> num_seq; i ++)	{
	free((void *) newseq[i]);
  }
  free((void **) newseq);
  free((void *) len_newseq);
  free((void *) num_intv);
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

  while ((copt=getopt(argc,argv,"s:i:c:a:o:v:q:n:")) != EOF) {
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
    case 'n':
      sscanf(optarg,"%s", newcontigfile);
      continue;
    case 'q':
      sscanf(optarg,"%s", newcontigqualfile);
      continue;
    case 'v':
      sscanf(optarg,"%s", outintvfile);
      continue;
    case 'c':
      sscanf(optarg,"%s", contigfile);
      continue;
    default:
      printf("makeintv -s SeqFile -i IntvFile -a AlignFile -c ContigFle -o OutFile -v OutIntvFile -n NewContigFile -q NewContigQualFile\n");
      printf("-s SeqFile: The input file name of reads\n");
      printf("-i IntvFile: The input file name of intervals\n");
      printf("-a AlignFile: The input file name of multiple alignments\n");
      printf("-c ContigFile: The output file name of contig sequences\n");
      printf("-o OutFile: The output file name modified reads\n");
      printf("-v OutIntvFile: The output file name modified read intervals\n");
      exit(-1);
    }
    optind--;
  }
  if (inpseq == 0 || outseq == 0) {
      printf("makeintv -s SeqFile -i IntvFile -a AlignFile -c ContigFle -o OutFile -v OutIntvFile -n NewContigFile -q NewContigQualFile\n");
      printf("-s SeqFile: The input file name of reads\n");
      printf("-i IntvFile: The input file name of intervals\n");
      printf("-a AlignFile: The input file name of multiple alignments\n");
      printf("-c ContigFile: The output file name of contig sequences\n");
      printf("-o OutFile: The output file name modified reads\n");
      printf("-v OutIntvFile: The output file name modified read intervals\n");
      exit(-1);
  }
}
