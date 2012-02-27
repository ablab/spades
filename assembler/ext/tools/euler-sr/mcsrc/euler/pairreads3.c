/***************************************************************************
 * Title:          pairreads3.c
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
	int	max_leg, max_name_leg;
	char	c;
	int	range[2];
	char	*ext_name;
	int	*plate;
	char	**library;
	char	str[100], temp[100];
	READTABLE	rt;
	int	*el;
	FILE	*fp;

	readpar();
	initenv(argc, argv);

	fp = ckopen(inpfile, "r");
	rt.num_seq = 0;
	while(fgets(str, 90, fp))	{
		if(str[0] == '>')	rt.num_seq ++;
	}
	fclose(fp);
	ext_name = (char *) ckalloc(rt.num_seq * sizeof(char));
	plate = (int *) ckalloc(rt.num_seq * sizeof(int));
	library = (char **) ckalloc(rt.num_seq * sizeof(char *));
	rt.src_name = (char **) ckalloc(rt.num_seq * sizeof(char *));
	for(i = 0; i < rt.num_seq; i ++)	{
		library[i] = (char *) ckalloc(100 * sizeof(char));
		rt.src_name[i] = (char *) ckalloc(100 * sizeof(char));
	}
	i = 0;
	fp = ckopen(inpfile, "r");
	while(fgets(str, 90, fp))	{
		if(str[0] == '>')	{
			sscanf(str, "%*s%s%s%s", library[i], rt.src_name[i], temp);
			ext_name[i] = temp[0];
			i ++;
		}
	}
	fclose(fp);
	printf("num_seq: %d\n", rt.num_seq);
	fp = ckopen(rulefile, "r");
	readrules(fp);
	fclose(fp);

	el = (int *) ckalloc(ntypepair * sizeof(int));
	for(i = 0; i < ntypepair; i ++)	{
		el[i] = 0;
	}

	n = 0;
	fp = ckopen(pairfile, "w");
	for(i = 0; i < rt.num_seq; i ++)	{
		c = pairrules(ext_name[i], plate[i], library[i], range);
		if(c != ext_name[i])	{
			k = loc_pair(rt.src_name[i], c, rt.src_name, ext_name, rt.num_seq);
			if(k >= 0)	{
				el[range[0]] ++;
				fprintf(fp, "%s %d %d %d\n", rt.src_name[i], i, k, range[0]);
				n ++;
			}
		}
	}
	fclose(fp);
	free((void *) plate);
	for(i = 0; i < rt.num_seq; i ++)	{
		free((void *) library[i]);
	}
	free((void **) library);

	printf("-------------------------------------------------------------------------\n");
	printf("EULER mate-pair INPUT:\n\n");
	for(i = 0; i < ntypepair; i ++)	{
		printf("%c-%c(library %s plates %d-%d): distance range from %d to %d, %d pairs.\n",
			pairrule[i][0], pairrule[i][1], librule[i],
			platerule[i][0], platerule[i][1], pairrange[i][0], pairrange[i][1], el[i]);
	}
	printf("Total %d sequences and %d pairs are found.\n", rt.num_seq, n);

	free((void *) ext_name);
	free((void *) el);
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
