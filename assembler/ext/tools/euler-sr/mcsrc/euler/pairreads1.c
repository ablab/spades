/***************************************************************************
 * Title:          pairreads1.c
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
	int	i, j, k, l, n, m, q;
	char	c;
	int	range[2];
	int	max_leg, max_name_leg;
	int	el[10];
	READTABLE	rt;
	FILE	*fp;

	readpar();
	initenv(argc, argv);

	inputreadname(&rt, inpfile);
	fp = ckopen(rulefile, "r");
	readrules(fp);
	fclose(fp);

	el[0] = el[1] = 0;
	fp = ckopen(pairfile, "w");
	for(i = 0; i < rt.num_seq; i ++)	{
		l = strlen(rt.src_name[i]);
		if(rt.src_name[i][l - 2] == '.' && rt.src_name[i][l - 1] == 'f')	{
			for(j = 0; j < rt.num_seq; j ++)	{
				k = strlen(rt.src_name[j]);
				if(k == l - 2 && !strncmp(rt.src_name[i], rt.src_name[j], l - 2))	{
					if(!strncmp(rt.src_name[i], "met1", 4))	{
						fprintf(fp, "%s %d %d %d\n", rt.src_name[j], j, i, 0);
						el[0] ++;
					} else	{
						fprintf(fp, "%s %d %d %d\n", rt.src_name[j], j, i, 1);
						el[1] ++;
					}
					break;
				}
			}
		}
	}
	fclose(fp);

	printf("-------------------------------------------------------------------------\n");
	printf("EULER mate-pair INPUT:\n\n");
	printf("Totally %d sequences and %d pairs are found.\n", rt.num_seq, el[0] + el[1]);
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

	while ((copt=getopt(argc,argv,"i:o:r:")) != EOF)	{
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
			default:
			  printf("pairreads -i SeqFile -o PairFile [-r RuleFile]\n\n");
			  printf("-i SeqFile: The file name of the inputting reads\n");
			  printf("-o PairFile: The file name of the output pairs\n");
			  printf("-r RuleFile(optional): name of rule file\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("pairreads -i SeqFile -o PairFile [-r RuleFile]\n\n");
		printf("-i SeqFile: The file name of the inputting reads\n");
		printf("-o PairFile: The file name of the output pairs\n");
		printf("-r RuleFile(optional): name of rule file\n");
		exit(-1);
	}
}
