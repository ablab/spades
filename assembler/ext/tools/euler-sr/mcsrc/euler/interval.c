/***************************************************************************
 * Title:          interval.c
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

void write_interval(NODES **vertex, int num_vertex, FILE *fp);
int read_interval(NODES **vertex, int num_vertex, FILE *fp);

void write_interval_file(char *fname, int num_vertex, NODES **vertex);

void write_interval(NODES **vertex, int num_vertex, FILE *fp)
{
	int	i, j, k, l, m, n;
	EDGE	*edge;

	n = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			fprintf(fp, "EDGE %d Length %d Multiplicity %d\n", n + 1, edge -> length, edge -> multip);
			for(k = 0; k < edge -> multip; k ++)	{
				fprintf(fp, "INTV %d %d %d %d\n", edge -> readinterval[k].eq_read, edge -> readinterval[k].begin,
					edge -> readinterval[k].length, edge -> readinterval[k].offset);
			}
			n ++;
		}
	}
}

int read_interval(NODES **vertex, int num_vertex, FILE *fp)
{
	int	i, j, k, l, m, n;
	char	str[500];
	EDGE	*edge;
	int num_intv;

	num_intv = 0;

	n = 0;
	for(i = 0; i < num_vertex; i ++)	{
		for(j = 0; j < vertex[i] -> num_nextedge; j ++)	{
			edge = vertex[i] -> nextedge[j];
			fgets(str, 400, fp);
			sscanf(str, "%*s%*d%*s%*d%*s%d\n", &(edge -> multip));
			num_intv += edge->multip;
			edge -> readinterval = (READINTERVAL *) ckalloc(edge -> multip * sizeof(READINTERVAL));
			for(k = 0; k < edge -> multip; k ++)	{
				fgets(str, 400, fp);
				sscanf(str, "%*s%d%d%d%d\n", &(edge -> readinterval[k].eq_read), &(edge -> readinterval[k].begin),
					&(edge -> readinterval[k].length), &(edge -> readinterval[k].offset));
			}
			n ++;
		}
	}
	printf("Input %d read-intervals on %d edges\n", num_intv, n);
	return num_intv;
}


/**********************************************************************
 * Input/output the read intervals in each edge
 **********************************************************************/

int read_interval_file(char *intvfile,
			int num_vertex,
			NODES **vertex)
{
  FILE *fp;
  int num_intv;

  fp = ckopen(intvfile, "r");
  num_intv = read_interval(vertex, num_vertex, fp);
  fclose(fp);

  return num_intv;
}


void write_interval_file(char *intvfile,
			 int num_vertex,
			 NODES **vertex)
{
  FILE *fp;

  fp = ckopen(intvfile, "w");
  write_interval(vertex, num_vertex, fp);
  fclose(fp);
}

