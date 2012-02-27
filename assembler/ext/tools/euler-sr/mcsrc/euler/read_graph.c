/***************************************************************************
 * Title:          read_graph.c
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

int read_graph(NODES **vertex, EDGE **edge, int *num_edge, FILE *fp, FILE *fp1);

int read_graph(NODES **vertex, EDGE **edge, int *num_edge, FILE *fp, FILE *fp1)
{
	int	i, j, k, l, n, m, d, num_vertex, p1, p2;
	char	targ_name[100], *targ_seq, str[500];

	targ_seq = (char *) ckalloc(MAX_TMP_LEG * sizeof(char));

	k = -1;
	while(fgets(str, 490, fp))	{
		if(str[0] == '>')	{
			sscanf(&str[2], "%s", targ_name);
			if(k >= 0)	{
				edge[k] = (EDGE *) ckalloc(1 * sizeof(EDGE));
				edge[k] -> seq = (char *) ckalloc(n * sizeof(char));
				for(i = 0; i < n; i ++)	{
					edge[k] -> seq[i] = targ_seq[i];
				}
				edge[k] -> length = n;
				edge[k]->index = k;
				edge[k]->deleted=0;
			}
			k ++;
			n = 0;
		} else {
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					targ_seq[n ++] = char2int(str[i]);
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					targ_seq[n ++] = char2int(str[i] - 'A' + 'a');
				}
			}
		}
	}
	/* process the last edge */
	edge[k] = (EDGE *) ckalloc(1 * sizeof(EDGE));
	edge[k] -> seq = (char *) ckalloc(n * sizeof(char));
	for(i = 0; i < n; i ++)	{
		edge[k] -> seq[i] = targ_seq[i];
	}
	edge[k] -> length = n;
	edge[k]->index = k;
	edge[k]->deleted = 0;
	k ++;

	*num_edge = k;
	printf("created %d edges\n", *num_edge);
	p1 = p2 = 0;

	n = -1;
	while(fgets(str, 490, fp1))	{
		if(!strncmp(str, "Number_of_Vertex", 16))	{
			sscanf(str, "%*s%d", &num_vertex);
		} else if(!strncmp(str, "Vertex", 6))	{
			n ++;
			vertex[n] = (NODES *) ckalloc(1 * sizeof(NODES));
			sscanf(str, "%*s%*d%d%d", &(vertex[n] -> num_nextedge), &(vertex[n] -> num_lastedge));
			vertex[n]->index = n;
			vertex[n] -> lastedge = (EDGE **) ckalloc(vertex[n] -> num_lastedge * sizeof(EDGE *));
			vertex[n] -> nextedge = (EDGE **) ckalloc(vertex[n] -> num_nextedge * sizeof(EDGE *));
			p1 += vertex[n] -> num_lastedge;
			p2 += vertex[n] -> num_nextedge;
		} else if(!strncmp(str, "Last_edge", 9))	{
			l = 9;
			for(j = 0; j < vertex[n] -> num_lastedge; j ++)	{
				sscanf(&str[l + 1], "%d %d", &m, &d);
				vertex[n] -> lastedge[j] = edge[m];
				vertex[n] -> lastedge[j] -> bal_edge = edge[d];
				//				printf("storing last edge %d vertex %d for pos %d\n", m, n, j);
				edge[m] -> end = vertex[n];
				for(l ++; l < strlen(str) - 1; l ++)	{
					if(str[l] == ' ')	{
						break;
					}
				}
				for(l ++; l < strlen(str) - 1; l ++)	{
					if(str[l] == ' ')	{
						break;
					}
				}
			}
		} else if(!strncmp(str, "Next_edge", 9))	{
			l = 9;
			for(j = 0; j < vertex[n] -> num_nextedge; j ++)	{
				sscanf(&str[l + 1], "%d %d", &m, &d);
				//				printf("storing nextedge %d vertex %d for pos %d\n", m, n, j);
				vertex[n] -> nextedge[j] = edge[m];
				vertex[n] -> nextedge[j] -> bal_edge = edge[d];
				edge[m] -> begin = vertex[n];
				for(l ++; l < strlen(str) - 1; l ++)	{
					if(str[l] == ' ')	{
						break;
					}
				}
				for(l ++; l < strlen(str) - 1; l ++)	{
					if(str[l] == ' ')	{
						break;
					}
				}
			}
		}
	}
	//	printf("done storing edges\n");
	num_vertex = n + 1;
	for(i = 0; i < num_vertex; i ++)	{
		if(vertex[i] -> num_lastedge > 0)	{
			vertex[i] -> bal_node = vertex[i] -> lastedge[0] -> bal_edge -> begin;
		} else	{
			vertex[i] -> bal_node = vertex[i] -> nextedge[0] -> bal_edge -> end;
		}
	}

	free((void *) targ_seq);
	return(num_vertex);
}


/**********************************************************************
 * Read the graph files 
 *
 * TODO:
 * 1. adjust allocation sizes to input, instead of assuming fixed sizes
 * 2. Return a graph structure instead of individual parameters
 **********************************************************************/

int CountEdges(FILE *edgeFile, int*numEdges) {
  char c;
  if (! edgeFile)
    return 0;

  while((c = fgetc(edgeFile)) != EOF)
    if (c == '>') {
      (*numEdges)++;
    }

  return 1;
}

int CountVertices(FILE *vertexFile, int *numVertices) {
  char str[490];
  if (fgets(str, 490, vertexFile))	{
    if(!strncmp(str, "Number_of_Vertex", 16))	{
      sscanf(str, "%*s%d", numVertices);
      return 1;
    }
    else {
      printf("Vertex file should start with 'Number_of_Vertex'\n");
      exit(1);
    }
  }
  return 0;
}

void read_graph_file(char *edgefile,
		     char *graphfile,
		     int *num_vertex,
		     NODES ***vertex,
		     int *num_edge,
		     EDGE ***edge)
{
  FILE *fp;
  FILE *fp1;
  int numEdges, numVertices;
  fp = ckopen(edgefile, "r");
  if (!CountEdges(fp, &numEdges)) {
    printf("Error counting edges\n");
    exit(1);
  }
  fclose(fp);
  fp = ckopen(edgefile, "r");


  fp1 = ckopen(graphfile, "r");
  if (!CountVertices(fp1, &numVertices)) {
    printf("Error counting vertices\n");
    exit(1);
  }
  fclose(fp1);
  fp1 = ckopen(graphfile, "r");
  *vertex = (NODES **) ckalloc(numVertices * sizeof(NODES *));
  *edge = (EDGE **) ckalloc(numEdges * sizeof(EDGE *));
  *num_vertex = read_graph(*vertex, *edge, num_edge, fp, fp1);
  printf("Input graph ... done. %d vertices and %d edges.\n",
	 *num_vertex, *num_edge);
  fclose(fp1);
  fclose(fp);
}
