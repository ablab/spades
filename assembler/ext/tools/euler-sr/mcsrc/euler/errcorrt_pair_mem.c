/***************************************************************************
 * Title:          errcorrt_pair_mem.c
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

extern char htmlout, caption[2000], ***content;
extern FILE *flog;

void errcorrt_pair_mem(ALIGN **eq_class, char **multip, int *len_seq, char **src_seq, char **score, char **src_name,
		   int num_seq, int LOW_COV, int intv, FILE *fp, FILE *fp1);

void errcorrt_pair_mem(ALIGN **eq_class, char **multip, int *len_seq, char **src_seq, char **score, char **src_name,
		   int num_seq, int LOW_COV, int intv, FILE *fp, FILE *fp1)
{
	int	i, j, k, l, m, n, na1;
	int	j1, j2, k1, k2;
	int	rr1, rr2, rp1, rp2;
	int	m1, m2, c1, c2, s1, s2;
	int	nerr;
	int	dist[10];
	int	read1, read2, pos1, pos2;
	int	end1, end2;
	ALIGN	*align;
	signed char	*newseq, **repseq, **insseq;
	char	*score_new, *label;

	for(k = 0; k < 10; k ++)		dist[k] = 0;
	repseq = (signed char **) ckalloc(2 * num_seq * sizeof(signed char *));
	insseq = (signed char **) ckalloc(2 * num_seq * sizeof(signed char *));
	for(i = 0; i < 2 * num_seq; i ++)	{
		repseq[i] = (signed char *) ckalloc((len_seq[i] + 20) * sizeof(signed char));
		insseq[i] = (signed char *) ckalloc((len_seq[i] + 20) * sizeof(signed char));
	}
	label = (char *) ckalloc(2 * num_seq * sizeof(char));
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
				end1 = align -> pos[0][j + 1];
				end2 = align -> pos[1][j + 1];
				while(pos1 < end1 && pos2 < end2)	{
					if(src_seq[read1][pos1] == src_seq[read2][pos2])	{
						if(mov_qual)	{
							score[read1][pos1] = score[read2][pos2] =
								max(1, max(score[read1][pos1], score[read2][pos2]));
						}
					} else if(!repseq[read1][pos1] && !repseq[read2][pos2] &&
					   !insseq[read1][pos1] && !insseq[read2][pos2])	{
						if(multip[read1][pos1] <= LOW_COV && (score[read2][pos2] >= 1 &&
						   multip[read2][pos2] >= LOW_COV * 2))	{
							score[read2][pos2] = max(1, score[read2][pos2]);
							repseq[read1][pos1] = src_seq[read2][pos2] + 1;
							rp1 = len_seq[read1] - pos1 - 1;
							rp2 = len_seq[read2] - pos2 - 1;
							repseq[rr1][rp1] = src_seq[rr2][rp2] + 1;
						}
						if(multip[read2][pos2] <= LOW_COV && (score[read1][pos1] >= 1 &&
						   multip[read1][pos1] >= LOW_COV * 2))	{
							score[read1][pos1] = max(1, score[read1][pos1]);
							repseq[read2][pos2] = src_seq[read1][pos1] + 1;
							rp1 = len_seq[read1] - pos1 - 1;
							rp2 = len_seq[read2] - pos2 - 1;
							repseq[rr2][rp2] = src_seq[rr1][rp1] + 1;
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
					c1 = min(multip[read1][k1], multip[read1][k2]);
					c2 = min(multip[read2][j1], multip[read2][j2]);
					s1 = min(score[read1][k1], score[read1][k2]);
					s2 = min(score[read1][j1], score[read1][j2]);
					if(c2 <= LOW_COV && (s1 == 2 || c1 >= LOW_COV * 2))	{
/*	insert (read1,pos1) into (read2, j1)	*/
						rp1 = len_seq[read1] - pos1 - 1;
						rp2 = len_seq[read2] - j1 - 1;
						insseq[read2][j1] = src_seq[read1][pos1] + 1;
						insseq[rr2][rp2 - 1] = src_seq[rr1][rp1] + 1;
					} else if(c1 <= LOW_COV && (s2 == 2 || c2 >= LOW_COV * 2))	{
/*	delete (read1, pos1-> k2 - 1)	*/
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
					c1 = min(multip[read1][k1], multip[read1][k2]);
					c2 = min(multip[read2][j1], multip[read2][j2]);
					s1 = min(score[read1][k1], score[read1][k2]);
					s2 = min(score[read2][j1], score[read2][j2]);
					if(c1 <= LOW_COV && (s2 == 2 || c2 >= LOW_COV * 2))	{
/*	insert (read2,pos2) into (read1, k1)	*/
						insseq[read1][k1] = src_seq[read2][pos2] + 1;
						rp1 = len_seq[read1] - k1 - 1;
						rp2 = len_seq[read2] - pos2 - 1;
						insseq[rr1][rp1 - 1] = src_seq[rr2][rp2] + 1;
					} else if(c2 <= LOW_COV && (s1 == 2 || c1 >= LOW_COV * 2))	{
/*	delete (read2, pos2->j2 - 1)	*/
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
	free((void *) label);
	newseq = (signed char *) ckalloc(MAX_TMP_LEG * sizeof(signed char));
	score_new = (char *) ckalloc(MAX_TMP_LEG * sizeof(char));
	for(i = 0; i < num_seq; i ++)	{
		j = 0;
		while(j < len_seq[i])	{
			if(insseq[i][j])	{
				for(m = 1; m < 10; m ++)	{
					insseq[i][j + m] = 0;
					repseq[i][j + m] = 0;
				}
				j += 10;
			} else	{
				j ++;
			}
		}
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
	for(i = 0; i < 2 * num_seq; i ++)	{
		free((void *) repseq[i]);
		free((void *) insseq[i]);
	}
	free((void **) repseq);
	free((void **) insseq);
	free((void *) score_new);
	free((void *) newseq);
	if(!htmlout)	{
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "The distribution of the quality values at error corrected positions.\n");
  		print_text_line(flog, LINE_LENGTH);
		fprintf(flog, "# errors        <20            20-30          >30            Total\n");
		fprintf(flog, "Substitutions   %-14d %-14d %-14d %-14d\n", dist[0], dist[1], dist[2], dist[0] + dist[1] + dist[2]);
		fprintf(flog, "Indels          %-14d %-14d %-14d %-14d\n", dist[3], dist[4], dist[5], dist[3] + dist[4] + dist[5]);
		fprintf(flog, "Total           %-14d %-14d %-14d %-14d\n", dist[0] + dist[3], dist[1] + dist[4], dist[2] + dist[5],
			 nerr);
  		print_text_line(flog, LINE_LENGTH);
	} else	{
		sprintf(caption, "The distribution of the quality values at error corrected positions.\n");
		content = allocate_content(4, 5, 50);
		sprintf(content[0][0], "# errors");
		sprintf(content[0][1], "< 20");
		sprintf(content[0][2], "< 20-30");
		sprintf(content[0][3], "> 30");
		sprintf(content[0][4], "Total");
		sprintf(content[1][0], "Substitution");
		sprintf(content[1][1], "%d", dist[0]);
		sprintf(content[1][2], "%d", dist[1]);
		sprintf(content[1][3], "%d", dist[2]);
		sprintf(content[1][4], "%d", dist[0] + dist[1] + dist[2]);
		sprintf(content[2][0], "Indel");
		sprintf(content[2][1], "%d", dist[3]);
		sprintf(content[2][2], "%d", dist[4]);
		sprintf(content[2][3], "%d", dist[5]);
		sprintf(content[2][4], "%d", dist[3] + dist[4] + dist[5]);
		sprintf(content[2][0], "Total");
		sprintf(content[3][1], "%d", dist[0] + dist[3]);
		sprintf(content[3][2], "%d", dist[1] + dist[4]);
		sprintf(content[3][3], "%d", dist[2] + dist[5]);
		sprintf(content[3][4], "%d", nerr);
		print_table(flog, 4, 5, content, caption);
		content = free_content(content, 4, 5);
	}
}
