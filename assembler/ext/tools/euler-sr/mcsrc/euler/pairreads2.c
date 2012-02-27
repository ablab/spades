/***************************************************************************
 * Title:          pairreads2.c
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

void cutname(char *part_name, char *src_name);
void initenv(int argc, char **argv);

main(int argc, char **argv)
{
	int	i, j, k, l, n, m, q, n1, n2;
	char	str[500];
	char	part_name[2][100];
	READTABLE	rt;
	FILE	*fp, *fp1;

	readpar();
	initenv(argc, argv);

	inputreadname(&rt, inpfile);
	fp = ckopen(rulefile, "r");
	readrules(fp);
	fclose(fp);

	n = 0;
	fp = ckopen(pairfile, "w");
	for(i = 0; i < rt.num_seq - 1; i ++)	{
		cutname(part_name[0], rt.src_name[i]);
		cutname(part_name[1], rt.src_name[i]);
		if(!strcmp(part_name[0], part_name[1]))	{
			l = strlen(part_name[0]);
			matetype = 0;
			for(j = l - 1; j >= 0; j --)	{
				if(part_name[0][j] == '_')	{
					break;
				}
				if(part_name[0][j] == 'p')	{
					matetype = 1;
				}
			}
			fprintf(fp, "%s %d %d %d\n", part_name[0], i, i + 1, matetype);
			n ++;
		}
	}
	fclose(fp);

	printf("-------------------------------------------------------------------------\n");
	printf("EULER mate-pair INPUT:\n\n");
	for(i = 0; i < ntypepair; i ++)	{
		printf("%c-%c(plates %d-%d): distance range from %d to %d.\n",
			pairrule[i][0], pairrule[i][1], 
			platerule[i][0], platerule[i][1], pairrange[i][0], pairrange[i][1]);
	}
	printf("There are %d of mates are found.\n", n);
	return(0);
}

void cutname(char *part_name, char *src_name)
{
	int	i, j, k, l;

	l = strlen(src_name);
	for(i = l - 1; i >= 0; i --)	{
		if(src_name[i] == '_')	{
			break;
		}
	}
	for(j = 0; j < i; j ++)	{
		part_name[j] = src_name[j];
	}
	part_name[j] = '\0';
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
			  printf("pairreads -i SeqFile -o PairFile -r RuleFile -t MateType\n\n");
			  printf("-i SeqFile: The file name of the inputting reads\n");
			  printf("-o PairFile: The file name of the output pairs\n");
			  printf("-r RuleFile: name of rule file\n");
			  printf("-t MateType: mate type index\n");
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		printf("pairreads -i SeqFile -o PairFile -r RuleFile -t MateType\n\n");
		printf("-i SeqFile: The file name of the inputting reads\n");
		printf("-o PairFile: The file name of the output pairs\n");
		printf("-r RuleFile: name of rule file\n");
		printf("-t MateType: mate type index\n");
		exit(-1);
	}
}
