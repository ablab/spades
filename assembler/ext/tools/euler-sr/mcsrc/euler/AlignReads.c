/***************************************************************************
 * Title:          AlignReads.c
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

#define HIGHEND 50

extern int qualinp;

ALIGN *AlignReadsPO(char **src_seq, int *len_seq, int num_seq, READOVERLAP *readoverlap);
ALIGN *AlignReads(char **src_seq, char **score, int *len_seq, int num_seq, READINTERVAL *readinterval, int index);
int overalign(char *seq1, char *seq2, char *score1, char *score2, 
	      int len_seq1, int len_seq2, int offset, int *sapp, int *cutp, int *seql,
	      int index1, int index2, int num_seq);
double perc_qual(char *score, int s1, int s2);
int cal_identity(char *seq1, int len1, char *seq2, int len2, int *sapp);
void compute_score(char *seq1, int len1, char *seq2, int len2, int *sapp, double *score);
ALIGN *new_align(int *cutp, int *seql, int *sapp, ALIGN *align, int r1, int r2, int offset, int mis_match);

ALIGN *AlignReadsPO(char **src_seq, int *len_seq, int num_seq, READOVERLAP *readoverlap)
{
	int	i, j, k, l, m, n;
	int	offset;
	ALIGN	*align;
	int	read1, read2;
	int	*sapp;
	int	c;
	int	cutp[2], seql[2];

	sapp = (int *) ckalloc(MAX_TMP_LEG * sizeof(int));
	align = (ALIGN *) NULL;

	while(readoverlap)	{
		cutp[0] = readoverlap -> pos1[0];
		cutp[1] = readoverlap -> pos2[0];
		seql[0] = readoverlap -> pos1[1] - readoverlap -> pos1[0] + 1;
		seql[1] = readoverlap -> pos2[1] - readoverlap -> pos2[0] + 1;
		read1 = readoverlap -> read1;
		read2 = readoverlap -> read2;
		k = ALIGN0(&src_seq[read1][cutp[0] - 1], &src_seq[read2][cutp[1] - 1], 
			seql[0], seql[1], -band, band, W, g, h, sapp, len_seq[read1] + len_seq[read2],
			len_seq[read1] + len_seq[read2]);
/*	transfer the LLR_score to align -> mis_match	*/
		align = new_align(cutp, seql, sapp, align, readoverlap -> read1, readoverlap -> read2,
			 1, readoverlap -> LLR_score);
		readoverlap = readoverlap -> next;
	}
	free((void *) sapp);
	return(align);
}

ALIGN *AlignReads(char **src_seq, char **score, int *len_seq, int num_seq, READINTERVAL *readinterval, int index)
{
	int	i, j, k, l, m, n;
	int	offset;
	ALIGN	*align;
	int	*sapp;
	int	c;
	int	cutp[2], seql[2];

	sapp = (int *) ckalloc(MAX_TMP_LEG * sizeof(int));
	align = (ALIGN *) NULL;

	while(readinterval)	{
		offset = readinterval -> offset;
		k = readinterval -> eq_read;
		if(offset >= 0)	{
			c = overalign(src_seq[index], src_seq[k], score[index], score[k],
				len_seq[index], len_seq[k], offset, sapp, cutp, seql, index, k, num_seq);
		} else	{
			c = overalign(src_seq[k], src_seq[index], score[k], score[index],
				len_seq[k], len_seq[index], -offset, sapp, cutp, seql, k, index, num_seq);
		}
/*
c = 1;
*/
		if(c >= 0)	{
		  //		  printf("found alignment %d to %d \n", index, k);
			align = new_align(cutp, seql, sapp, align, index, k, offset, c);
		}
		readinterval = readinterval -> next;
	}
	free((void *) sapp);
	return(align);
}

int overalign(char *seq1, char *seq2, char *score1, char *score2, 
	      int len_seq1, int len_seq2, int offset, int *sapp, int *cutp, int *seql,
	      int index1, int index2, int num_seq)
{
	int	i, j, k, l, e1, e2;
	int	st, ln;
	int	pre, suf;
	double	p, p1, p2, p3, p4;
	int	beg1, beg2, ent1, ent2, range[8];
	char	c1, c2, c3, c4;

	st = max(0, offset - band);
	ln = len_seq1 - st;

	k = LOCAL_ALIGN(&seq1[st - 1], &seq2[-1], &score1[st - 1], &score2[-1], ln, len_seq2, -2 * band, 2 * band,
		 W, g, h, &cutp[0], &cutp[1], &e1, &e2, len_seq1 + len_seq2);

/*	Use regular alignment	*/

/*
	k = LOCAL_ALIGN0(&seq1[st - 1], &seq2[-1], ln, len_seq2, band * 2, 0,
		 W, g, h, &cutp[0], &cutp[1], &e1, &e2, len_seq1 + len_seq2);
*/

	if(k < 0)	{
		return(-1);
	}
	//	printf("results: %d %d \n", cutp[0], cutp[1]);
	seql[0] = e1 - cutp[0] + 1;
	seql[1] = e2 - cutp[1] + 1;
	if(seql[0] < overlaplen || seql[1] < overlaplen)	{
		return(-1);
	}
	cutp[0] += st - 1;
	e1 += st - 1;
	cutp[1] -= 1;
	e2 -= 1;
	pre = min(cutp[0], cutp[1]);
	suf = min(len_seq1 - e1 - 1, len_seq2 - e2 - 1);
	c1 = c2 = c3 = c4 = 0;
	if(index1 < num_seq)	{
		beg1 = cutp[0];
		ent1 = e1;
	} else	{
		beg1 = len_seq1 - e1 - 1;
		ent1 = len_seq1 - cutp[0] - 1;
	}
	if(index2 < num_seq)	{
		beg2 = cutp[1];
		ent2 = e2;
	} else	{
		beg2 = len_seq2 - e2 - 1;
		ent2 = len_seq2 - cutp[1] - 1;
	}
	range[0] = cutp[0] - pre;
	range[1] = cutp[0] - 1;
	range[2] = e1 + 1;
	range[3] = e1 + suf;
	range[4] = cutp[1] - pre;
	range[5] = cutp[1] - 1;
	range[6] = e2 + 1;
	range[7] = e2 + suf;
	if(qualinp)	{
		p1 = perc_qual(score1, range[0], range[1]);
		l = (1 - p1) * (range[1] - range[0]);
		if(l < HIGHEND && p1 > MIN_PERC || range[1] - range[0] < SHORT_D)	{
			c1 = 1;
		}
	} else	{
		if(beg1 < start_ent)	{
			c1 = 1;
		}
	}
	if(qualinp)	{
		p2 = perc_qual(score2, range[4], range[5]);
		l = (1 - p2) * (range[5] - range[4]);
		if(l < HIGHEND && p2 > MIN_PERC || range[5] - range[4] < SHORT_D)	{
			c2 = 1;
		}
	} else	{
		if(beg2 < start_ent)	{
			c2 = 1;
		}
	}
	if(qualinp)	{
		p3 = perc_qual(score1, range[2], range[3]);
		l = (1 - p3) * (range[3] - range[2]);
		if(l < HIGHEND && p3 > MIN_PERC || range[3] - range[2] < SHORT_D)	{
			c3 = 1;
		}
	} else	{
		if(ent1 > end_ent || len_seq1 - ent1 - 1 < MIN_END)	{
			c3 = 1;
		}
	}
	if(qualinp)	{
		p4 = perc_qual(score2, range[6], range[7]);
		l = (1 - p4) * (range[7] - range[6]);
		if(l < HIGHEND && p4 > MIN_PERC || range[7] - range[6] < SHORT_D)	{
			c4 = 1;
		}
	} else	{
		if(ent2 > end_ent || len_seq2 - ent2 - 1 < MIN_END)	{
			c4 = 1;
		}
	}
	if((c1 || c2) && (c3 || c4))	{
		k = GALIGN(&seq1[cutp[0] - 1], &seq2[cutp[1] - 1], &score1[cutp[0] - 1], &score2[cutp[1] - 1],
			seql[0], seql[1], -2 * band, 2 * band, W, g, h, sapp, len_seq1 + len_seq2,
			len_seq1 + len_seq2);

/*	Use regular alignment	*/

/*
		k = ALIGN0(&seq1[cutp[0] - 1], &seq2[cutp[1] - 1], 
			seql[0], seql[1], -band, band, W, g, h, sapp, len_seq1 + len_seq2,
			len_seq1 + len_seq2);
*/
		k = cal_identity(&seq1[cutp[0]], seql[0], &seq2[cutp[1]], seql[1], sapp);
		if(k >= 0)	{
			return(k);
		} else {
			return(-1);
		}
	} else	{
		return(-1);
	}
}

int cal_identity(char *seq1, int len1, char *seq2, int len2, int *sapp)
{
	double	s[2];

	compute_score(seq1, len1, seq2, len2, sapp, s);
	if(s[0] > MIN_OVERLAP && s[1] >= MIN_IDENTITY)	{
/*
	if(s[0] > MIN_OVERLAP)	{
*/
		return(1);
	} else	{
		return(-1);
	}
}

void compute_score(char *seq1, int len1, char *seq2, int len2, int *sapp, double *score)
{
	int	i, j, k, l, m, n;
	int	op;
	double	p, s;

	s = 0;
	i = j = k = l = 0;
	while(i < len1 || j < len2)	{
		op = *sapp ++;
		if(op == 0)	{
			if(seq1[i] == seq2[j])	{
				s += 1.0;
				k ++;
			} else	{
				s -= MIS_SCORE;
				l ++;
			}
			n ++;
			i ++;
			j ++;
		} else if(op > 0)	{
			for(m = 0; m < op; m ++)	{
				s -= MIS_SCORE;
				l ++;
			}
			j += op;
		} else 	{
			for(m = 0; m < -op; m ++)	{
				s -= MIS_SCORE;
				l ++;
			}
			i -= op;
		}
	}
	if(k + l > 0)	{
		score[1] = p = ((double) k) / (k + l);
	}
	score[0] = s;
}

double perc_qual(char *score, int s1, int s2)
{
	int	i, j, k;

	k = 0;
	for(i = s1; i <= s2; i ++)	{
		if(score[i] == 0)	{
			k ++;
		}
	}
	if(s1 > s2)	{
		return(1.0);
	} else	{
		return(((double) k) / (s2 - s1 + 1));
	}
}

ALIGN *new_align(int *cutp, int *seql, int *sapp, ALIGN *align, int r1, int r2, int offset, int mis_match)
{
	int	i, j, k, l, m, n, len, mark;
	int	op;
	int	*pos[2];
	ALIGN	*newalign, *aln, *aln_last;

	newalign = (ALIGN *) ckalloc(1 * sizeof(ALIGN));
	newalign -> reads[0] = r1;
	newalign -> reads[1] = r2;
	len = seql[0] + seql[1];
	for(i = 0; i < 2; i ++)	{
		pos[i] = (int *) ckalloc(len * sizeof(int));
	}
	i = j = k = mark = 0;
	while(i < seql[0] && j < seql[1])	{
		op = *sapp ++;
		if(op == 0)	{
			if(mark == 0)	{
				mark = 1;
				pos[0][k] = cutp[0] + i;
				pos[1][k] = cutp[1] + j;
				k ++;
			}
			i ++;
			j ++;
		} else if(op > 0)	{
			mark = 0;
			j += op;
		} else 	{
			mark = 0;
			i -= op;
		}
	}
	pos[0][k] = cutp[0] + i;
	pos[1][k] = cutp[1] + j;
	k ++;
	newalign -> length = k;
	newalign -> mis_match = mis_match;
	newalign -> pos[0] = (int *) ckalloc(k * sizeof(int));
	newalign -> pos[1] = (int *) ckalloc(k * sizeof(int));
	if(offset >= 0)	{
		for(i = 0; i < k; i ++)	{
			newalign -> pos[0][i] = pos[0][i];
			newalign -> pos[1][i] = pos[1][i];
		}
	} else	{
		for(i = 0; i < k; i ++)	{
			newalign -> pos[0][i] = pos[1][i];
			newalign -> pos[1][i] = pos[0][i];
		}
	}

/*	Remove multiple overlapers between two reads	*/
/*
	aln = align;
	aln_last = (ALIGN *) NULL;
	while(aln)	{
		if(aln -> reads[1] == newalign -> reads[1])	{
			if(newalign -> pos[1][newalign -> length - 1] - newalign -> pos[1][0] <
			   aln -> pos[1][aln -> length - 1] - aln -> pos[1][0])	{
				break;
			} else	{
				if(!aln_last)	{
					align = aln -> next;
				} else	{
					aln_last -> next = aln -> next;
				}
				free_align(aln);
			}
		}
		aln_last = aln;
		aln = aln -> next;
	}
	if(!aln)	{
		newalign -> next = align;
		align = newalign;
	}
*/
	newalign -> next = align;
	align = newalign;

	free((void *) pos[0]);
	free((void *) pos[1]);
	return(align);
}
