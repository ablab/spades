/***************************************************************************
 * Title:          output_contig.c
 * Author:         Haixu Tang
 * Created:        Jun. 2003
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

void output_contig_cons(NODES **vertex, int num_vertex, char **src_name, char **src_seq,
		 int num_seq, FILE *fp, FILE *fp1, FILE *fp2);

void writeseq(FILE *fp, char *seq, char *name, int length);
int output_contig_files(char *filestem, int num_vertex, NODES **vertex, 
												PATH * paths, int numPaths,
												READTABLE *RT);

int output_contig(NODES **vertex, int num_vertex, 
									PATH *path, int numPaths,
									char **src_name, char **src_seq,
									int num_seq, FILE *fp, FILE *fp1, FILE *fp2, FILE *fp3, 
									FILE *bgraphFP, FILE *pathFP)
{
	int	i, j, k, l, n, m;
	char	temp[100];
	EDGE	*edge;

	write_graph(vertex, num_vertex, fp, fp1);
	write_branching_graph(vertex, num_vertex, bgraphFP);
	/*
	if (pathFP != NULL) 
		write_paths(vertex, num_vertex, path, numPaths, pathFP);
	*/
	fclose(bgraphFP);
	fflush(fp);
	fflush(fp1);
	n = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			edge -> visit = 0;
			n ++;
		}
	}
	m = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			if(edge -> begin -> num_lastedge == 0 && edge -> end -> num_nextedge == 0 &&
			   edge -> multip <= 1)	continue;
			if(edge -> visit == 0)	{
				sprintf(temp, "Contig%d", m + 1);
				writeseq(fp2, edge -> seq, temp, edge -> length);
/*	Temporatory: output intervals to ace file now	*/
				fprintf(fp3, "Contig %d Length %d Multiplicity %d\n", m + 1, edge -> length, edge -> multip);
				for(k = 0; k < edge -> multip; k ++)	{
					fprintf(fp3, "INTV %d %d %d %d\n", edge -> readinterval[k].eq_read, edge -> readinterval[k].begin,
						edge -> readinterval[k].length, edge -> readinterval[k].offset);
				}
/*
				writeace(fp3, edge -> seq, edge -> class, edge -> multip, edge -> length, src_name,
					 src_seq, num_seq);
*/
				edge -> visit = edge -> bal_edge -> visit = 1;
				m ++;
			}
		}
	}
	return(n);
}

void output_contig_cons(NODES **vertex, int num_vertex, char **src_name, char **src_seq,
		 int num_seq, FILE *fp, FILE *fp1, FILE *fp2)
{
	int	i, j, k, l, n;
	char	temp[100];
	EDGE	*edge;

	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			edge -> visit = 0;
		}
	}
	n = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			writeseq(fp, edge -> seq, temp, edge -> length);
			if(edge -> visit == 0)	{
				sprintf(temp, "Contig%d", n);
				writeseq(fp1, edge -> seq, temp, edge -> length);
/*
				writeace(fp2, edge -> seq, edge -> class, edge -> multip, edge -> length, src_name,
					 src_seq, num_seq);
*/
				n ++;
				edge -> visit = edge -> bal_edge -> visit = 1;
			}
		}
	}
}

void writeseq(FILE *fp, char *seq, char *name, int length)
{
	int	i, j, k;

	fprintf(fp, ">%s %d\n", name, length);
	for(i = 0; i < length; i ++)	{
		fprintf(fp, "%c", na_name[seq[i]]);
		if(i % 50 == 49)	{
			fprintf(fp, "\n");
		}
	}
	if(i % 50 != 0)	{
		fprintf(fp, "\n");
	}
}


int output_contig_files(char *filestem,
												int num_vertex,
												NODES **vertex,
												PATH *paths, int numPaths,
												READTABLE *RT)
{
  int n;
  char *fname;
	fname = (char*) ckalloc(strlen(filestem) + 12);
  FILE *fp, *fp1, *fp2, *fp3, *bgraphFP, *pathFP;

  sprintf(fname, "%s.edge", filestem);
  fp = ckopen(fname, "w");
  sprintf(fname, "%s.graph", filestem);
  fp1 = ckopen(fname, "w");
  sprintf(fname, "%s.contig", filestem);
  fp2 = ckopen(fname, "w");
  sprintf(fname, "%s.contig.ace", filestem);
  fp3 = ckopen(fname, "w");
	sprintf(fname, "%s.bgraph", filestem);
	bgraphFP = ckopen(fname, "w");
	sprintf(fname, "%s.path", filestem);
	pathFP = ckopen(fname, "w");
  n = output_contig(vertex, num_vertex,	paths, numPaths, 
										RT->src_name, RT->src_seq, RT->num_seq, 
										fp, fp1, fp2, fp3, bgraphFP, pathFP);
  fclose(fp);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
	free(fname);
  return(n);

}

