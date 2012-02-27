/***************************************************************************
 * Title:          pairreads.c
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
	READTABLE	rt;
	int	*el;
	FILE	*fp;

	readpar();
	initenv(argc, argv);
	if(htmlout)	{
		flog = ckopen("EULER-report.html", "a");
	} else	{
		flog = ckopen("EULER-report.txt", "a");
	}

	if(!htmlout)	{
  	        print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "Defining mate-pairs based on read names.\n\n");
	} else	{
	       	print_section_head(flog, "Defining mate-pairs based on read names.");
	}

/*	Input the reads (required) */

	inputreadname(&rt, inpfile);
	fp = ckopen(rulefile, "r");
	readrules(fp);
	fclose(fp);

	ext_name = (char *) ckalloc(rt.num_seq * sizeof(char));
	plate = (int *) ckalloc(rt.num_seq * sizeof(int));
	library = (char **) ckalloc(rt.num_seq * sizeof(char *));
	for(i = 0; i < rt.num_seq; i ++)	{
		library[i] = (char *) ckalloc(100 * sizeof(char));
	}
	fp = ckopen(rulefile, "r");
	for(i = 0; i < rt.num_seq; i ++)	{
		ext_name[i] = splitname(rt.src_name[i], &plate[i], library[i]);
	}
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

	if(!htmlout)	{
		printf("-------------------------------------------------------------------------\n");
		printf("EULER mate-pair INPUT:\n\n");
		for(i = 0; i < ntypepair; i ++)	{
			printf("%c-%c(library %s plates %d-%d): distance range from %d to %d, %d pairs.\n",
				pairrule[i][0], pairrule[i][1], librule[i],
				platerule[i][0], platerule[i][1], pairrange[i][0], pairrange[i][1], el[i]);
		}
		printf("Total %d sequences and %d pairs are found.\n", rt.num_seq, n);
	} else	{
		sprintf(caption, "Total %d sequences and %d pairs are found.\n", rt.num_seq, n);
		content = allocate_content(ntypepair + 1, 3, 50);
		sprintf(content[0][0], "%s", "library and plates");
		sprintf(content[0][1], "%s", "distance range");
		sprintf(content[0][2], "%s", "# pairs");
		for(i = 0; i < ntypepair; i ++)	{
			sprintf(content[i + 1][0], "%c-%c(library %s plates %d-%d)",
				pairrule[i][0], pairrule[i][1], librule[i],
				platerule[i][0], platerule[i][1]);
			sprintf(content[i + 1][1], "%d-%d", pairrange[i][0], pairrange[i][1]);
			sprintf(content[i + 1][2], "%d", el[i]);
		}
		print_table(flog, ntypepair + 1, 3, content, caption);
		content = free_content(content, ntypepair + 1, 3);
	}

	if(!htmlout)	{
  		print_text_line(flog, LINE_LENGTH);
	} else	{
		print_hl(flog);
	}
	fclose(flog);

/*	Free Memory	*/

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

	while ((copt=getopt(argc,argv,"i:o:r:H")) != EOF)	{
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
			case 'H':
			  htmlout = 1;
			  continue;
			default:
			  if(!htmlout)	{
				  printf("pairreads -i SeqFile -o PairFile [-r RuleFile]\n\n");
				  printf("-i SeqFile: The file name of the inputting reads\n");
				  printf("-o PairFile: The file name of the output pairs\n");
				  printf("-r RuleFile(optional): name of rule file\n");
			  } else	{
				  print_line(flog, "pairreads -i SeqFile -o PairFile [-r RuleFile]");
				  print_line(flog, "-i SeqFile: The file name of the inputting reads");
				  print_line(flog, "-o PairFile: The file name of the output pairs");
				  print_line(flog, "-r RuleFile(optional): name of rule file");
			  }
			  exit(-1);
		}
		optind--;
	}

	if(inpseq == 0 || outseq == 0)	{
		if(!htmlout)	{
			  printf("pairreads -i SeqFile -o PairFile [-r RuleFile]\n\n");
			  printf("-i SeqFile: The file name of the inputting reads\n");
			  printf("-o PairFile: The file name of the output pairs\n");
			  printf("-r RuleFile(optional): name of rule file\n");
		} else	{
			  print_line(flog, "pairreads -i SeqFile -o PairFile [-r RuleFile]");
			  print_line(flog, "-i SeqFile: The file name of the inputting reads");
			  print_line(flog, "-o PairFile: The file name of the output pairs");
			  print_line(flog, "-r RuleFile(optional): name of rule file");
		}
		exit(-1);
	}
}
