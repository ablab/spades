/***************************************************************************
 * Title:          countnodes.c
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

void statspath(PATH *path, int num_path);
int singstatnode(NODES *node, int **nstat);
void statnode(READLIST **readlist, int *len_seq, char **src_seq, int num_seq);
int countnode(READLIST **readlist, int *len_seq, int num_seq);
int cleannode(READLIST **readlist, int *len_seq, int num_seq);

int countnode(READLIST **readlist, int *len_seq, int num_seq)
{
	int	i, j, k, l, n;
	READPOSITION *readposition;
	/* 
	   Count the number of unique nodes in the graph.  readlist[i][j].node
	   and readlist[k][l].node may point to the same node.
	*/
	n = 0;
	for(i = 0; i < num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			if(readlist[i][j].node && readlist[i][j].node -> visit == 0)	{
				n ++;
				readlist[i][j].node -> visit = 1;
			}
		}
	}
	cleannode(readlist, len_seq, num_seq/2);
	return(n);
}

void statnode(READLIST **readlist, int *len_seq, char **src_seq, int num_seq)
{
	int	i, j, k, l, n, max_l;
	int	**stat, **nstat, *lax;
	int	width;
	READPOSITION *readposition;

	stat = (int **) ckalloc(500 * sizeof(int *));
	nstat = (int **) ckalloc(500 * sizeof(int *));
	lax = (int *) ckalloc(500 * sizeof(int));
	for(i = 0; i < 500; i ++)	{
		stat[i] = (int *) ckalloc(500 * sizeof(int));
		nstat[i] = (int *) ckalloc(500 * sizeof(int));
	}

	for(i = 0; i < num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			for(k = 0; k < 5; k ++)	{
				lax[k] = 0;
			}
			if(readlist[i][j].node && readlist[i][j].node -> visit == 0)	{
				readlist[i][j].node -> visit = 1;
				width = singstatnode(readlist[i][j].node, nstat);
				if(width > 1)	continue;
				readposition = readlist[i][j].node -> readposition;
				while(readposition)	{
					if(readposition -> position >= 0)	{
						lax[src_seq[readposition -> readindex][readposition -> position]] ++;
					} else	{
						lax[4] ++;
					}
					readposition = readposition -> next;
				}
				max_l = lax[0];
				for(k = 1; k < 5; k ++)	{
					if(lax[k] > max_l)	{
						max_l = lax[k];
					}
				}
				if(max_l < 10 && readlist[i][j].node -> npos < 30)	{
					stat[readlist[i][j].node -> npos][max_l] ++;
				} else if(max_l < 10)	{
					stat[30][max_l] ++;
				} else if(readlist[i][j].node -> npos < 30)	{
					stat[readlist[i][j].node -> npos][10] ++;
				} else	{
					stat[30][10] ++;
				}
			}
		}
	}
	printf("----------------------------------------------------------------------------------------\n");
	printf("# of supernodes that contains n nodes and the major node has multiplicity k\n");
	printf("----------------------------------------------------------------------------------------\n");
	printf("major   1      2      3      4      5      6      7      8      9      10+     Total\n");
	printf("----------------------------------------------------------------------------------------\n");
	for(i = 1; i <= 30; i ++)	{
		if(i < 30)	{
			printf("%-7d ", i);
		} else	{
			printf("30+     ", i);
		}
		k = 0;
		for(j = 1; j <= 10; j ++)	{
			printf("%-6d ", stat[i][j]);
			k += stat[i][j];
		}
		printf("%-6d", k);
		printf("\n");
	}
	printf("----------------------------------------------------------------------------------------\n");
	printf("# of supernodes with width k and thickness n\n");
	printf("----------------------------------------------------------------------------------------\n");
	printf("          2         3         4         5         6        7+\n");
	printf("----------------------------------------------------------------------------------------\n");
	for(i = 2; i <= 50; i ++)	{
		if(i < 50)	{
			printf("%-9d ", i);
		} else	{
			printf("50+       ");
		}
		for(j = 2; j <= 7; j ++)	{
			printf("%-9d ", nstat[i][j]);
		}
		printf("\n");
	}
	printf("----------------------------------------------------------------------------------------\n");
	cleannode(readlist, len_seq, num_seq/2);
	for(i = 0; i < 500; i ++)	{
		free((void *) stat[i]);
		free((void *) nstat[i]);
	}
	free((void **) stat);
	free((void **) nstat);
	free((void *) lax);
}

int singstatnode(NODES *node, int **nstat)
{
	int	i, j, k, l, n, width, thickness;

	width = countwidth(node);
	thickness = countthickness(node);
	if(width < 7)	{
		if(thickness <= 50)	{
			nstat[thickness][width] ++;
		} else	{
			nstat[50][width] ++;
		}
	} else	{
		if(thickness <= 50)	{
			nstat[thickness][7] ++;
		} else	{
			nstat[50][7] ++;
		}
	}
	return(width);
}

int cleannode(READLIST **readlist, int *len_seq, int num_seq)
{
	int	i, j, k, l, n;

	n = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		for(j = 0; j < len_seq[i]; j ++)	{
			if(readlist[i][j].node)	{
				readlist[i][j].node -> visit = 0;
			}
		}
	}
	return(n);
}

void statspath(PATH *path, int num_path)
{
	int	i, j;
	int	dist[20];

	for(i = 0; i < 10; i ++)	dist[i] = 0;
	for(i = 0; i < num_path; i ++)	{
		if(path[i].len_path >= 7)	{
			dist[7] ++;
		} else	{
			dist[path[i].len_path] ++;
		}
		dist[8] ++;
	}
	printf("# paths: %d.\n", num_path);
	printf("-----------------------------------------------------------------------------\n");
	printf("        1        2        3        4        5        6       7+    total\n");
	printf("-----------------------------------------------------------------------------\n");
	for(i = 1; i < 9; i ++)	{
		printf("%9d", dist[i]);
	}
	printf("\n");
	printf("-----------------------------------------------------------------------------\n");
}
