/***************************************************************************
 * Title:          po_prod.c
 * Author:         Haixu Tang
 * Created:        Jul. 2004
 * Last modified:  Jul. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <stdinc.h>
#include <param.h>
#include <extfunc.h>

char	inpfile[100], outfile[100], seqfile[100];

void initenv(int argc, char **argv);
int findname(char *name, READTABLE *RT);

main(int argc, char **argv)
{
  int i, m, l, j, k, n;

  int k1, k2, pos1, pos2, pos3, pos4, score, len;
  char temp[100], str[500], name1[100], name2[100];
  char rev;
  FILE *fp, *fp1;

  READTABLE RT_mem, *RT=&RT_mem;


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
   * Input the overlaps from cross_match
   **********************************************************************/

  fp = ckopen(inpfile, "r");
  fp1 = ckopen(outfile, "w");
  while(fgets(str, 190, fp))	{
	if(!strncmp(str, "Maximal single", 14))	{
		break;
	}
  }
  while(fgets(str, 190, fp))	{
	if(!strncmp(str, "Discrepancy", 11))	{
		break;
	}
	len = strlen(str);
	if(len < 50)	continue;
	rev = 0;
	for(i = 1; i < len - 1; i ++)	{
		if(str[i - 1] == ' ' && str[i] == 'C' && str[i + 1] == ' ')	{
			str[i] = ' ';
			rev = 1;
			break;
		}
	}
	if(rev)	{
		sscanf(str, "%d%*s%*s%*s%s%d%d%*s%s%*s%d%d", &score, name1, &pos1, &pos2,
			name2, &pos3, &pos4);
	} else	{
		sscanf(str, "%d%*s%*s%*s%s%d%d%*s%s%d%d", &score, name1, &pos1, &pos2,
			name2, &pos3, &pos4);
	}
	pos1 --;
	pos2 --;
	pos3 --;
	pos4 --;
	k1 = findname(name1, RT);
	k2 = findname(name2, RT);
	if(rev)	{
		pos4 = RT -> len_seq[k2] - 1 - pos4;
		pos3 = RT -> len_seq[k2] - 1 - pos3;
		k2 += RT -> num_seq;
	}
	fprintf(fp1, "%d %d %d %d %d %d %d\n", k1, k2, pos1, pos2, pos3, pos4, score);
  }
  fclose(fp);
  fclose(fp1);

  free_readtable(RT);
  return(0);
}

int findname(char *name, READTABLE *RT)
{
	int	i;

	for(i = 0; i < RT -> num_seq; i ++)	{
		if(!strcmp(RT -> src_name[i], name))	{
			return(i);
		}
	}
	printf("Read %s not found.\n", name);
	exit(-1);
}

void initenv(int argc, char **argv)
{
  int copt;
  int inpseq, outseq;
  extern char *optarg;

  inpseq = outseq = 0;
  MIN_OVERLAP = 100;

  while ((copt=getopt(argc,argv,"o:s:i:")) != EOF) {
    switch(copt) {
    case 'i':
      inpseq = 1;
      sscanf(optarg,"%s", inpfile);
      continue;
    case 's':
      sscanf(optarg,"%s", seqfile);
      continue;
    case 'o':
      outseq = 1;
      sscanf(optarg,"%s", outfile);
      continue;
    default:
	     printf("po-prod -i InpFile -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP -x VERTEX_SIZE]\n");
	     printf("-i InpFile: The input file name of overlaps\n");
	     printf("-s SeqFile: The input file name of reads\n");
	     printf("-o Outfile: The output file name of contigs\n");
	     printf("-c LOW_COV (optional): minimal coverage\n");
	     printf("-p MIN_OVERLAP (optional): minimal overlap length\n");
	     exit(-1);
    }
    optind--;
  }

  if (inpseq == 0 || outseq == 0) {
      printf("po-prod -i InpFile -s SeqFile -o outfile [-c LOW_COV -p MIN_OVERLAP -x VERTEX_SIZE]\n");
      printf("-i InpFile: The input file name of overlaps\n");
      printf("-s SeqFile: The input file name of reads\n");
      printf("-o Outfile: The output file name of contigs\n");
      printf("-c LOW_COV (optional): minimal coverage\n");
      printf("-p MIN_OVERLAP (optional): minimal overlap length\n");
      exit(-1);
  }
}
