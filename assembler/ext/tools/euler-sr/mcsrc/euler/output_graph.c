/***************************************************************************
 * Title:          output_graph.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

extern int nsuper, *numtangle, *maxmultip, *maxlength, *avelength, *mlength;
extern char	htmlout, ***content, caption[2000];
extern FILE	*flog;
extern void print_text_line(FILE *fp, int length);
/* TODO:
 * 1. Use IGRAPH structure
 * 2. Use consistent parameter order
 */

void output_graph(NODES **vertex, int num_vertex, FILE *fp);
void write_gvz_graph(FILE *fp, NODES **vertex, int num_vertex);

void write_gvz_file(char *fname, int num_vertex, NODES **vertex, int verbose);

void write_sf_gvz_file(char *contig_fname,
		       char *suffix,
		       SFPAIRS *SFP);

void output_graph(NODES **vertex, int num_vertex, FILE *fp)
{
	int	i, j, k, l, m, n, len;
	int	tot_edge;
	int	**num_pa, disttangle[8][7], num_tangle;
	int	*label;
	char	temp[1000];
	EDGE	*edge0;
	NODES	*v;

	num_pa = (int **) ckalloc(MAX_BRA * sizeof(int *));
	for(i = 0; i < MAX_BRA; i ++)	{
		num_pa[i] = (int *) ckalloc(MAX_BRA * sizeof(int));
	}

	for(i = 0; i < num_vertex; i ++)	{
		vertex[i] -> num_path = i + 1;
	}

	tot_edge = 0;
	for(i = 0; i < num_vertex; i ++)	{
		tot_edge += vertex[i] -> num_nextedge;
	}
	write_gvz_graph(fp, vertex, num_vertex);
	numtangle = (int *) ckalloc(tot_edge * sizeof(int));;
	maxmultip = (int *) ckalloc(tot_edge * sizeof(int));;
	mlength = (int *) ckalloc(tot_edge * sizeof(int));;
	maxlength = (int *) ckalloc(tot_edge * sizeof(int));;
	avelength = (int *) ckalloc(tot_edge * sizeof(int));;
	nsuper = 0;
	count_multip(vertex, num_vertex);
	num_tangle = count_tangle(vertex, num_vertex, disttangle);
	tot_edge = count_edge(vertex, num_vertex, num_pa);
	len = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			len += vertex[i] -> nextedge[j] -> length - 1;
		}
	}
	if(!htmlout)	{
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "Vertices Edge Source Sink Tangles Super-tangles Overall-length\n");
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "%-8d %-4d %-6d %-4d %-7d %-13d %-14d\n",
			num_vertex, tot_edge, num_pa[0][1], num_pa[1][0], num_tangle, nsuper, len);
  		print_text_line(flog, LINE_LENGTH);
		if(nsuper > 0)	{
			fprintf(flog, "\nStatistics of the super-tangles:\n");
  			print_text_line(flog, LINE_LENGTH);
			fprintf(flog, "Supertangle #tangles max-multip length-of-the-central-edge max-length ave-length\n");
  			print_text_line(flog, LINE_LENGTH);
		}
		for(j = 0; j < nsuper; j ++)	{
			fprintf(flog, "%-11d %-8d %-9d  %-26d %10d %11d\n", j + 1, numtangle[j], maxmultip[j],
				 mlength[j], maxlength[j], avelength[j]);
			numtangle[j] = maxmultip[j] = mlength[j] = 0;
		}
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "Distribution of vertex degrees:\n");
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "         \\ Indegree     0    1    2    3    4    5    >5\n");
		fprintf(flog, "Outdegree \\       \n");
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "   0                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[0][0],
			num_pa[0][1], num_pa[0][2], num_pa[0][3], num_pa[0][4], num_pa[0][5], num_pa[0][6]);
		fprintf(flog, "   1                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[1][0],
			num_pa[1][1], num_pa[1][2], num_pa[1][3], num_pa[1][4], num_pa[1][5], num_pa[1][6]);
		fprintf(flog, "   2                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[2][0],
			num_pa[2][1], num_pa[2][2], num_pa[2][3], num_pa[2][4], num_pa[2][5], num_pa[2][6]);
		fprintf(flog, "   3                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[3][0],
			num_pa[3][1], num_pa[3][2], num_pa[3][3], num_pa[3][4], num_pa[3][5], num_pa[3][6]);
		fprintf(flog, "   4                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[4][0],
			num_pa[4][1], num_pa[4][2], num_pa[4][3], num_pa[4][4], num_pa[4][5], num_pa[4][6]);
		fprintf(flog, "   5                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[5][0],
			num_pa[5][1], num_pa[5][2], num_pa[5][3], num_pa[5][4], num_pa[5][5], num_pa[5][6]);
		fprintf(flog, "  >5                 %4d %4d %4d %4d %4d %4d %4d\n", num_pa[6][0],
			num_pa[6][1], num_pa[6][2], num_pa[6][3], num_pa[6][4], num_pa[6][5], num_pa[6][6]);
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "\nNumber of tangles (repeat edges): %d.\n", num_tangle);
		fprintf(flog, "Distribution of tangle multiplicities:\n");
  		print_text_line(flog, LINE_LENGTH);
       	        fprintf(flog, "            \\ Length     <500     <1000    <5000   <10000   >10000   Total\n");
		fprintf(flog, "Multiplicity \\   \n");
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "    1 (non-repeated) %8d %8d %8d %8d %8d %8d\n", disttangle[0][0], disttangle[0][1],
			disttangle[0][2], disttangle[0][3], disttangle[0][4], disttangle[0][5]);
		fprintf(flog, "    2                %8d %8d %8d %8d %8d %8d\n", disttangle[1][0], disttangle[1][1],
			disttangle[1][2], disttangle[1][3], disttangle[1][4], disttangle[1][5]);
		fprintf(flog, "    3                %8d %8d %8d %8d %8d %8d\n", disttangle[2][0], disttangle[2][1],
			disttangle[2][2], disttangle[2][3], disttangle[2][4], disttangle[2][5]);
		fprintf(flog, "    4                %8d %8d %8d %8d %8d %8d\n", disttangle[3][0], disttangle[3][1],
			disttangle[3][2], disttangle[3][3], disttangle[3][4], disttangle[3][5]);
		fprintf(flog, "    5                %8d %8d %8d %8d %8d %8d\n", disttangle[4][0], disttangle[4][1],
			disttangle[4][2], disttangle[4][3], disttangle[4][4], disttangle[4][5]);
		fprintf(flog, "    6                %8d %8d %8d %8d %8d %8d\n", disttangle[5][0], disttangle[5][1],
			disttangle[5][2], disttangle[5][3], disttangle[5][4], disttangle[5][5]);
		fprintf(flog, "   >6                %8d %8d %8d %8d %8d %8d\n", disttangle[6][0], disttangle[6][1],
			disttangle[6][2], disttangle[6][3], disttangle[6][4], disttangle[6][5]);
		fprintf(flog, "Total                %8d %8d %8d %8d %8d %8d\n", disttangle[7][0], disttangle[7][1],
			disttangle[7][2], disttangle[7][3], disttangle[7][4], disttangle[7][5]);
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "The overall length of the edges is %d.\n", len);
	} else	{
/*	Print summary Table	*/
		sprintf(temp, "The overall length of the edges is %d.", len);
		print_line(flog, temp);
		content = allocate_content(2, 7, 50);
		strcpy(content[0][0], "Vertices");
		strcpy(content[0][1], "Edges");
		strcpy(content[0][2], "Sources");
		strcpy(content[0][3], "Sinks");
		strcpy(content[0][4], "Tangles");
		strcpy(content[0][5], "Super-tangles");
		strcpy(content[0][6], "Overall-length");
		sprintf(content[1][0], "%d", num_vertex);
		sprintf(content[1][1], "%d", tot_edge);
		sprintf(content[1][2], "%d", num_pa[0][1]);
		sprintf(content[1][3], "%d", num_pa[1][0]);
		sprintf(content[1][4], "%d", num_tangle);
		sprintf(content[1][5], "%d", nsuper);
		sprintf(content[1][6], "%d", len);
		sprintf(caption, "Summary of the repeat graph");
		print_table(flog, 2, 7, content, caption);
		content = free_content(content, 2, 7);
/*	Print Super-tangle Table	*/
		if(nsuper > 0)	{
			content = allocate_content(nsuper + 1, 6, 50);
			strcpy(content[0][0], "Super-tangle");
			strcpy(content[0][1], "# tangles");
			strcpy(content[0][2], "Maximum multiplicity");
			strcpy(content[0][3], "Length of the central edge");
			strcpy(content[0][4], "Maximum edge length");
			strcpy(content[0][5], "Average length");
			for(j = 1; j <= nsuper; j ++)	{
				sprintf(content[j][0], "%d", j);;
				sprintf(content[j][1], "%d", numtangle[j - 1]);;
				sprintf(content[j][2], "%d", maxmultip[j - 1]);;
				sprintf(content[j][3], "%d", mlength[j - 1]);;
				sprintf(content[j][4], "%d", maxlength[j - 1]);;
				sprintf(content[j][5], "%d", avelength[j - 1]);;
				numtangle[j - 1] = maxmultip[j - 1] = mlength[j - 1] = 0;
			}
			sprintf(caption, "Statistics of the super-tangles");
			print_table(flog, nsuper + 1, 6, content, caption);
			content = free_content(content, nsuper + 1, 6);
		}
/*	Output degree distribution table	*/
		content = allocate_content(8, 8, 50);
		for(i = 0; i < 8; i ++)	{
			for(j = 0; j < 8; j ++)	{
				if(i == 0)	{
					if(j > 0 && j < 7)	{
						sprintf(content[i][j], "%d", j - 1);
					} else if(j == 7)	{
						sprintf(content[i][j], ">5");
					}
				} else {
					if(j == 0)	{
						if(i == 7)	{
							sprintf(content[i][j], ">5");
						} else	{
							sprintf(content[i][j], "%d", i - 1);
						}
					} else	{
						sprintf(content[i][j], "%d", num_pa[i - 1][j - 1]);
					}
				}
			}
		}
		sprintf(caption, "Distribution of vertex degrees");
		print_table(flog, 8, 8, content, caption);
		content = free_content(content, 8, 8);
/*	Output tangle description table	*/
  		print_hl(flog);
		sprintf(temp, "Number of tangles (repeat edges): %d.", num_tangle);
		print_line(flog, temp);
		content = allocate_content(9, 7, 50);
		sprintf(content[0][1], "Length<500");
		sprintf(content[0][2], "500-1000");
		sprintf(content[0][3], "1001-5000");
		sprintf(content[0][4], "5001-10000");
		sprintf(content[0][5], ">10000");
		sprintf(content[0][6], "Total");
		for(i = 1; i < 9; i ++)	{
			if(i == 1)	{
				sprintf(content[i][0], "1 (non-repeated)");
			} else if(i == 7)	{
				sprintf(content[i][0], ">6");
			} else if(i == 8)	{
				sprintf(content[i][0], "Total");
			} else	{
				sprintf(content[i][0], "%d", i);
			}
			for(j = 1; j < 7; j ++)	{
				sprintf(content[i][j], "%d", disttangle[i - 1][j - 1]);
			}
		}
		sprintf(caption, "Distribution of tangle multiplicities");
		print_table(flog, 9, 7, content, caption);
		content = free_content(content, 9, 7);
	}

	free((void *) numtangle);
	free((void *) maxmultip);
	free((void *) mlength);
	free((void *) maxlength);
	free((void *) avelength);
	for(i = 0; i < MAX_BRA; i ++)	{
		free((void *) num_pa[i]);
	}
	free((void **) num_pa);
}

void write_gvz_graph(FILE *fp, NODES **vertex, int num_vertex)
{
	int	i, j, k, l, m;
	EDGE	*edge0;
	fprintf(fp, "digraph G {\n");
	fprintf(fp, "\tsize=\"8,8\";\n");
	for(m = 0; m < num_vertex; m ++)	{
		for(j = 0; j < vertex[m] -> num_nextedge; j ++)	{
			edge0 = vertex[m] -> nextedge[j];
			if(edge0 -> end -> num_nextedge == 0 && edge0 -> begin -> num_lastedge == 0)	continue;
			fprintf(fp, "\t%d -> %d [", edge0 -> begin , edge0 -> end );
			fprintf(fp, "label = \"L:%d M:%d I:%d", edge0 -> length, edge0 -> multip, edge0->index);
			fprintf(fp, ")\"];\n");
			fprintf(fp, "%d [label = \"%d\"]\n", edge0->begin, edge0->begin->index);
			fprintf(fp, "%d [label = \"%d\"]\n", edge0->end, edge0->end->index);
			
		}
	}
	fprintf(fp, "}\n");
}

void write_gvz_file(char *fname, int num_vertex, NODES **vertex, int verbose)
{
  FILE *fp;

  fp = ckopen(fname, "w");

  if (verbose) {
    output_graph(vertex, num_vertex, fp);
  } else {
    write_gvz_graph(fp, vertex, num_vertex);
  }

  fclose(fp);
}


void write_sf_gvz_file(char *contig_fname,
		       char *suffix,
		       SFPAIRS *SFP)
{
  char temp[1000];
  FILE *fp;
  int i;

  sprintf(temp, "%s.sf%s.gvz", contig_fname, suffix);
  fp = ckopen(temp, "w");
  fprintf(fp, "digraph G {\n");
  fprintf(fp, "\tsize=\"6,6\";\n");
  fprintf(fp, "\tlabel=\"EULER-SF graph\\nContig file: %s\\nNodes: 10 means contig10, R10 means its reverse complement.\\nEdges: 5 means supported by 5 reads.\";\n", contig_fname);
  for (i = 0; i < SFP->tot_edge; i ++) {
	if(SFP->all_edge[i] -> length >= LINK_MIN_LEN)	{
      		fprintf(fp, "\t%dL%d;\n", SFP->all_edge[i] -> start_cover,
			SFP->all_edge[i] -> length);
	}
  }
  for (i = 0; i < SFP->num_sf; i ++) {
    if(SFP->link_cov[i] <= LINK_COV || SFP -> linkedge1[i] -> length < LINK_MIN_LEN
       || SFP -> linkedge2[i] -> length < LINK_MIN_LEN)		continue;
    if (SFP->linkedge1[i] -> start_cover > 0) {
      fprintf(fp, "\t%dL%d ", SFP->linkedge1[i] -> start_cover,
		SFP->linkedge1[i] -> length);
    } else {
      fprintf(fp, "\tR%dL%d ", SFP->linkedge1[i] -> bal_edge -> start_cover,
		SFP->linkedge2[i] -> bal_edge -> length);
    }
    fprintf(fp, "-> ");
    if (SFP->linkedge2[i] -> start_cover > 0) {
      fprintf(fp, "%dL%d [", SFP->linkedge2[i] -> start_cover,
		SFP->linkedge1[i] -> length);
    } else {
      fprintf(fp, "R%dL%d [", SFP->linkedge2[i] -> bal_edge -> start_cover,
		SFP->linkedge2[i] -> bal_edge -> length);
    }
    if (SFP->link_cov[i] > 3) {
      fprintf(fp, "style=bold,weight=8,");
    }
    fprintf(fp, "label = \"%d\"];\n", SFP->link_cov[i]);
  }
  fprintf(fp, "}\n");
  fclose(fp);
}
