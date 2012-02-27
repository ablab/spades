/***************************************************************************
 * Title:          errcorrt_pair.c
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

void errcorrt_pair(ALIGN **eq_class, int **multip, int *len_seq, char **src_seq, char **score, char **src_name,
		   int num_seq, int LOW_COV, int intv, FILE *fp, FILE *fp1);

void errcorrt_pair(ALIGN **eq_class, int **multip, int *len_seq, char **src_seq, char **score, char **src_name,
		   int num_seq, int LOW_COV, int intv, FILE *fp, FILE *fp1)
{
	int	i, j, k, l, m, n, na1;
	int	j1, j2, k1, k2;
	int	rr1, rr2, rp1, rp2;
	int	m1, m2, c1, c2;
	int	**newmultip;
	int	nerr, n1, n2;
	int	dist[10];
	int	read1, read2, pos1, pos2;
	ALIGN	*align;
	signed char	*newseq, **repseq, **insseq;
	char	*score_new;

	for(k = 0; k < 10; k ++)		dist[k] = 0;
	repseq = (signed char **) ckalloc(2 * num_seq * sizeof(signed char *));
	insseq = (signed char **) ckalloc(2 * num_seq * sizeof(signed char *));
	newmultip = (int **) ckalloc(2 * num_seq * sizeof(int *));
	for(i = 0; i < 2 * num_seq; i ++)	{
		repseq[i] = (signed char *) ckalloc(len_seq[i] * sizeof(signed char));
		insseq[i] = (signed char *) ckalloc(len_seq[i] * sizeof(signed char));
		newmultip[i] = (int *) ckalloc(len_seq[i] * sizeof(int));
	}
	nerr = 0;
	for(i = 0; i < 2 * num_seq; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			rr1 = reverse_read(read1, num_seq);
			rr2 = reverse_read(read2, num_seq);
			for(j = 0; j < align -> length - 1; j ++)	{
				pos1 = align -> pos[0][j];
				pos2 = align -> pos[1][j];
				while(pos1 < align -> pos[0][j + 1] && pos2 < align -> pos[1][j + 1])	{
					if(src_seq[read1][pos1] != src_seq[read2][pos2] &&
					   !repseq[read1][pos1] && !repseq[read2][pos2] &&
					   !insseq[read1][pos1] && !insseq[read2][pos2])	{
						if(multip[read1][pos1] <= LOW_COV &&
						   (multip[read1][pos1] < multip[read2][pos2]) ||
						   (multip[read1][pos1] == multip[read2][pos2] &&
						    read1 < read2) && newmultip[read1][pos1] < multip[read2][pos2])	{
							repseq[read1][pos1] = src_seq[read2][pos2] + 1;
							rp1 = len_seq[read1] - pos1 - 1;
							rp2 = len_seq[read2] - pos2 - 1;
							repseq[rr1][rp1] = src_seq[rr2][rp2] + 1;
							newmultip[read1][pos1] = multip[read2][pos2];
							newmultip[rr1][rp1] = multip[rr2][rp2];
						} else if(multip[read2][pos2] <= LOW_COV &&
							  newmultip[read2][pos2] < multip[read1][pos1])	{
							repseq[read2][pos2] = src_seq[read1][pos1] + 1;
							rp1 = len_seq[read1] - pos1 - 1;
							rp2 = len_seq[read2] - pos2 - 1;
							repseq[rr2][rp2] = src_seq[rr1][rp1] + 1;
							newmultip[read2][pos2] = multip[read1][pos1];
							newmultip[rr2][rp2] = multip[rr1][rp1];
						}
					}
					pos1 ++;
					pos2 ++;
				}
				if(pos1 < align -> pos[0][j + 1] && !insseq[read1][pos1] &&
				   !insseq[read2][pos2 - 1])	{
					k1 = pos1 - 1;
					k2 = align -> pos[0][j + 1];
					j1 = pos2 - 1;
/*	pos2 == align -> pos[1][j + 1]	*/
					j2 = pos2;
					c1 = (multip[read1][k1] + multip[read1][k2]) / 2;
					c2 = (multip[read2][j1] + multip[read2][j2]) / 2;
					n1 = (newmultip[read1][k1] + newmultip[read1][k2]) / 2;
					n2 = (newmultip[read2][j1] + newmultip[read2][j2]) / 2;
					if(c2 <= LOW_COV && (c2 < c1 || c2 == c1 && read2 < read1) && n2 < c1)	{
/*	insert (read1,pos1) into (read2, j1)	*/
						m1 = max(0, j1 - intv);
						m2 = min(len_seq[read2] - 1, j2 + intv);
						for(m = m1; m <= m2; m ++)	{
							rp2 = len_seq[read2] - m - 1;
							insseq[read2][m] = 0;
							insseq[rr2][rp2 - 1] = 0;
							newmultip[read2][m] = c1;
							newmultip[rr2][rp2] = c1;
						}
						insseq[read2][j1] = src_seq[read1][pos1] + 1;
						rp1 = len_seq[read1] - pos1 - 1;
						rp2 = len_seq[read2] - j1 - 1;
						insseq[rr2][rp2 - 1] = src_seq[rr1][rp1] + 1;
					} else if(c1 <= LOW_COV && n1 < c2)	{
/*	delete (read1, pos1-> k2 - 1)	*/
						m1 = max(0, k1 - intv);
						m2 = min(len_seq[read1] - 1, k2 + intv);
						for(m = m1; m <= m2; m ++)	{
							rp1 = len_seq[read1] - m - 1;
							insseq[read1][m] = 0;
							insseq[rr1][rp1 - 1] = 0;
							newmultip[read1][m] = c2;
							newmultip[rr1][rp1] = c2;
						}
						for(m = pos1; m < k2; m ++)	{
							insseq[read1][m] = -1;
							rp1 = len_seq[read1] - m - 1;
							insseq[rr1][rp1] = -1;
						}
					}
				} else if(pos2 < align -> pos[1][j + 1] && !insseq[read1][pos1 - 1] &&
					  !insseq[read2][pos2])	{
					j1 = pos2 - 1;
					j2 = align -> pos[1][j + 1];
					k1 = pos1 - 1;
/*	pos1 == align -> pos[0][j + 1]	*/
					k2 = pos1;
					c1 = (multip[read1][k1] + multip[read1][k2]) / 2;
					c2 = (multip[read2][j1] + multip[read2][j2]) / 2;
					n1 = (newmultip[read1][k1] + newmultip[read1][k2]) / 2;
					n2 = (newmultip[read2][j1] + newmultip[read2][j2]) / 2;
					if(c1 <= LOW_COV && (c1 < c2 || c1 == c2 && read1 < read2) && n1 < c2)	{
/*	insert (read2,pos2) into (read1, k1)	*/
						m1 = max(0, k1 - intv);
						m2 = min(len_seq[read2] - 1, k2 + intv);
						for(m = m1; m <= m2; m ++)	{
							rp1 = len_seq[read1] - m - 1;
							insseq[read1][m] = 0;
							insseq[rr1][rp1 - 1] = 0;
							newmultip[read1][m] = c2;
							newmultip[rr1][rp1] = c2;
						}
						insseq[read1][k1] = src_seq[read2][pos2] + 1;
						rp1 = len_seq[read1] - k1 - 1;
						rp2 = len_seq[read2] - pos2 - 1;
						insseq[rr1][rp1 - 1] = src_seq[rr2][rp2] + 1;
					} else if(c2 <= LOW_COV && n2 < c1)	{
/*	delete (read2, pos2->j2 - 1)	*/
						m1 = max(0, j1 - intv);
						m2 = min(len_seq[read2] - 1, j2 + intv);
						for(m = m1; m <= m2; m ++)	{
							rp2 = len_seq[read2] - m - 1;
							insseq[read2][m] = 0;
							insseq[rr2][rp2 - 1] = 0;
							newmultip[read2][m] = c1;
							newmultip[rr2][rp2] = c1;
						}
						for(m = pos2; m < j2; m ++)	{
							insseq[read2][m] = -1;
							rp2 = len_seq[read2] - m - 1;
							insseq[rr2][rp2] = -1;
						}
					}
				}
			}
			align = align -> next;
		}
	}
	newseq = (signed char *) ckalloc(MAX_TMP_LEG * sizeof(char));
	score_new = (char *) ckalloc(MAX_TMP_LEG * sizeof(char));
	for(i = 0; i < num_seq; i ++)	{
		k = 0;
		for(j = 0; j < len_seq[i]; j ++)	{
			if(insseq[i][j] >= 0)	{
				if(repseq[i][j] == 0)	{
					score_new[k] = score[i][j];
					newseq[k ++] = src_seq[i][j];
				} else	{
					score_new[k] = 1;
					newseq[k ++] = repseq[i][j] - 1;
					dist[score[i][j]] ++;
					nerr ++;
				}
				if(insseq[i][j] > 0)	{
					score_new[k] = 1;
					newseq[k ++] = insseq[i][j] - 1;
					dist[score[i][j] + 3] ++;
					nerr ++;
				}
			} else	{
				dist[score[i][j] + 3] ++;
				nerr ++;
			}
		}

/*	Output revised sequences	*/

		fprintf(fp, ">%s\n", src_name[i]);
		for(j = 0; j < k; j ++)	{
			fprintf(fp, "%c", na_name[newseq[j]]);
			if(j % 60 == 59)	{
				fprintf(fp, "\n");
			}
		}
		if(j % 60 != 0)	{
			fprintf(fp, "\n");
		}
		fprintf(fp1, ">%s\n", src_name[i]);
		for(j = 0; j < k; j ++)	{
			fprintf(fp1, "%d", score_new[j]);
			if(j % 60 == 59)	{
				fprintf(fp1, "\n");
			}
		}
		if(j % 60 != 0)	{
			fprintf(fp1, "\n");
		}
	}
	free((void *) newseq);
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void *) repseq[i]);
		free((void *) insseq[i]);
		free((void *) newmultip[i]);
	}
	free((void **) repseq);
	free((void **) insseq);
	free((void **) newmultip);
	free((void *) score_new);
	printf("---------------------------------------------------------------------\n");
	printf("The distribution of the quality values at error corrected positions.\n");
	printf("---------------------------------------------------------------------\n");
	printf("# errors        <20            20-30          >30            Total\n");
	printf("Substitutions   %-14d %-14d %-14d %-14d\n", dist[0], dist[1], dist[2], dist[0] + dist[1] + dist[2]);
	printf("Indels          %-14d %-14d %-14d %-14d\n", dist[3], dist[4], dist[5], dist[3] + dist[4] + dist[5]);
	printf("Total           %-14d %-14d %-14d %-14d\n", dist[0] + dist[3], dist[1] + dist[4], dist[2] + dist[5],
		 nerr);
	printf("---------------------------------------------------------------------\n");
}
