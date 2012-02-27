/***************************************************************************
 * Title:          pairreads_file.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <rule.h>
#include <param.h>
#include <extfunc.h>

char rulefile[200], inpfile[200], pairfile[200];
char	htmlout, caption[2000], ***content;
FILE	*flog;

void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, n, m, q, n1, n2;
	char	clone1[100], clone2[100], cname[100];
	int	max_leg, max_name_leg;
	char	str[500];
	READTABLE	rt;
	FILE	*fp, *fp1;

	readpar();
	initenv(argc, argv);

	inputreadname(&rt, inpfile);

	fp = ckopen(rulefile, "r");
	fp1 = ckopen(pairfile, "w");
	n = 0;
	while(fgets(str, 200, fp))	{
		sscanf(str, "%s%s", clone1, clone2);
		n1 = n2 = -1;
		for(i = 0; i < rt.num_seq; i ++)	{
			if(!strcmp(clone1, rt.src_name[i]))	{
				n1 = i;
			} else if(!strcmp(clone2, rt.src_name[i]))	{
				n2 = i;
			}
		}
		if(n1 >= 0 && n2 >= 0)	{
			fprintf(fp1, "%s %d %d %d\n", clone1, n1, n2, matetype);
			n ++;
		}
	}
	fclose(fp);
	fclose(fp1);

	printf("There are %d of type%d mates are found.\n", n, matetype);
	return(0);
}


void initenv(int argc, char **argv)
{
	int copt;
	int inpseq, outseq;	
	extern char *optarg;

	inpseq = outseq = 0;
	htmlout = 0;

	strcpy(rulefile, "name.rul");

	while ((copt=getopt(argc,argv,"i:o:r:t:")) != EOF)	{
		switch(copt) {
			case 'i':
			  inpseq = 1;
			  sscanf(optarg,"%s", inpfile);
			  continue;
			case 'o':
			  outseq = 1;
			  sscanf(optarg,"%s", pairfile);
			  continue;
			case 'r':
			  sscanf(optarg,"%s", rulefile);
			  continue;
			case 't':
			  sscanf(optarg,"%d", &matetype);
			  continue;
			default:
			  printf("pairreads_file -i SeqFile -o PairFile -r RuleFile -t MateType\n\n");
			  printf("-i SeqFile: The file name of the inputting reads\n");
			  printf("-o PairFile: The file name of the output pairs\n");
			  printf("-r RuleFile: name of rule file\n");
			  printf("-t MateType: mate type index\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("pairreads_file -i SeqFile -o PairFile -r RuleFile -t MateType\n\n");
		printf("-i SeqFile: The file name of the inputting reads\n");
		printf("-o PairFile: The file name of the output pairs\n");
		printf("-r RuleFile: name of rule file\n");
		printf("-t MateType: mate type index\n");
		exit(-1);
	}
}
