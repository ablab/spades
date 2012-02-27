/***************************************************************************
 * Title:          filter_class.c
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

#define MIN_GAP 50
#define MAX_PERC1 0.80
#define MAX_PERC2 0.95

int dist_range(int p1, int p2);
void update_bound(READINTERVAL *leftreadinterval, READINTERVAL *rightreadinterval, READINTERVAL **readinterval2, int p1, int p2);
int filter_readinterval(ALIGN **eq_class, READTABLE *RT);
double perc_maxnuc(char *seq, int d1, int d2, double *perc);
void caln_score(char *seq1, int len1, char *seq2, int len2, double *score, ALIGN *align);

int filter_readinterval(ALIGN **eq_class, READTABLE *RT)
{
	int	i, j, k, l, m, n;
	int	pos1, pos2, pos3, pos4, read1, read2;
	int	l1, l2, l3, l4, d1, d2, d3, d4;
	int	length;
	double	score[2], *leftbestscore, *rightbestscore;
	int	*leftmost, *rightmost;
	double	perc1[2],perc2[2];
	ALIGN	*align, *align0, *align1;

/*	Determine the leftmost and rightmost positions for each read	*/
	leftmost = (int *) ckalloc(RT -> num_seq * 2 * sizeof(int));
	rightmost = (int *) ckalloc(RT -> num_seq * 2 * sizeof(int));
	for(i = 0; i < RT -> num_seq * 2; i ++)	{
		leftmost[i] = RT -> len_seq[i] - 1;
	}
	for(i = 0; i < RT -> num_seq * 2; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			d1 = align -> pos[0][0];
			d2 = align -> pos[0][align -> length - 1] - 1;
			d3 = align -> pos[1][0];
			d4 = align -> pos[1][align -> length - 1] - 1;
			l1 = reverse_read(read1, RT -> num_seq);
			l2 = reverse_read(read2, RT -> num_seq);
			if(d1 < leftmost[read1])	{
				leftmost[read1] = d1;
				rightmost[l1] = RT -> len_seq[read1] - 1 - d1;
			}
			if(d2 > rightmost[read1])	{
				rightmost[read1] = d2;
				leftmost[l1] = RT -> len_seq[read1] - 1 - d2;
			}
			if(d3 < leftmost[read2])	{
				leftmost[read2] = d3;
				rightmost[l2] = RT -> len_seq[read2] - 1 - d3;
			}
			if(d4 > rightmost[read2])	{
				rightmost[read2] = d4;
				leftmost[l2] = RT -> len_seq[read2] - 1 - d4;
			}
			align = align -> next;
		}
	}

/*	Find the best alignment score of each read	*/

	leftbestscore = (double *) ckalloc(RT -> num_seq * 2 * sizeof(double));
	rightbestscore = (double *) ckalloc(RT -> num_seq * 2 * sizeof(double));
	for(i = 0; i < RT -> num_seq * 2; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			caln_score(RT -> src_seq[read1], RT -> len_seq[read1],
				RT -> src_seq[read2], RT -> len_seq[read2], score, align);
			d1 = align -> pos[0][0];
			d2 = align -> pos[0][align -> length - 1] - 1;
			d3 = align -> pos[1][0];
			d4 = align -> pos[1][align -> length - 1] - 1;
/* r1	    -------->
   r2	  ------------>	*/
			if(d1 - leftmost[read1] < MIN_GAP && rightmost[read1] - d2 < MIN_GAP ||
			   d3 - leftmost[read2] < MIN_GAP && rightmost[read2] - d4 < MIN_GAP)	{
/* r1	-------->
   r2	    ---------->	*/
			} else if(d1 - leftmost[read1] > MIN_GAP && rightmost[read2] - d4 > MIN_GAP)	{
				if(score[0] > rightbestscore[read1])	{
					rightbestscore[read1] = score[0];
					m = reverse_read(read1, RT -> num_seq);
					leftbestscore[m] = score[0];
				}
				if(score[0] > leftbestscore[read2])	{
					leftbestscore[read2] = score[0];
					m = reverse_read(read2, RT -> num_seq);
					rightbestscore[m] = score[0];
				}
			}
/* r1	           -------->
   r2	    ---------->	*/
			if(d3 - leftmost[read2] > MIN_GAP && rightmost[read1] - d2 > MIN_GAP)	{
				if(score[0] > rightbestscore[read2])	{
					rightbestscore[read2] = score[0];
					m = reverse_read(read2, RT -> num_seq);
					leftbestscore[m] = score[0];
				}
				if(score[0] > leftbestscore[read1])	{
					leftbestscore[read1] = score[0];
					m = reverse_read(read1, RT -> num_seq);
					rightbestscore[m] = score[0];
				}
			}
			align = align -> next;
		}
	}

/*	Determine the threshold of the overlap score between one read with the others.
	Score_threshold = Best_score * FILTER_THRESH_PERC.
	default: FILTER_THRESH_PERC = 0.6
	The overlaps between read i and j are not filtered if and only of
	Score >= max(Score_threshold(i), Score_threshold(j), min_good_overlap);
*/

	for(i = 0; i < RT -> num_seq * 2; i ++)	{
		leftbestscore[i] = leftbestscore[i] * FILTER_THRESH_PERC;
		rightbestscore[i] = rightbestscore[i] * FILTER_THRESH_PERC;
	}
	

/*	Mark redundunt overlaps	*/

	for(i = 0; i < RT -> num_seq * 2; i ++)	{
		align = eq_class[i];
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			d1 = align -> pos[0][0];
			d2 = align -> pos[0][align -> length - 1] - 1;
			d3 = align -> pos[1][0];
			d4 = align -> pos[1][align -> length - 1] - 1;
/*	remove overlaps in the low complexity regions */
			perc_maxnuc(RT -> src_seq[read1], d1, d2, perc1);
			perc_maxnuc(RT -> src_seq[read2], d3, d4, perc2);
			if(perc1[0] > MAX_PERC1 || perc2[0] > MAX_PERC1 ||
			   perc1[1] > MAX_PERC2 || perc2[1] > MAX_PERC2)	{
				align -> cov = 1;
				align = align -> next;
				continue;
			}
			caln_score(RT -> src_seq[read1], RT -> len_seq[read1],
				RT -> src_seq[read2], RT -> len_seq[read2], score, align);
/* r1	    -------->
   r2	  ------------>	*/
			if(d1 - leftmost[read1] < MIN_GAP && rightmost[read1] - d2 < MIN_GAP ||
			   d3 - leftmost[read2] < MIN_GAP && rightmost[read2] - d4 < MIN_GAP)	{
				align -> cov = 0;
			} else if(d1 - leftmost[read1] > MIN_GAP && score[0] < rightbestscore[read1] ||
/* r1	-------->
   r2	    ---------->	*/
			   rightmost[read2] - d4 > MIN_GAP && score[0] < leftbestscore[read2] ||
/* r1	           -------->
   r2	    ---------->	*/
			   d3 - leftmost[read2] > MIN_GAP && score[0] < rightbestscore[read2] ||
			   rightmost[read1] - d2 > MIN_GAP && score[0] < leftbestscore[read1])	{
/*
if(read1 == 2072 || read2 == 2072 || read1 == reverse_read(2072, RT -> num_seq) || read2 == reverse_read(2072, RT -> num_seq))	{
	printf("score %f read1 %d %d %f %f read2 %d %d %f %f\n", score[0], read1, reverse_read(read1, RT -> num_seq), leftbestscore[read1], rightbestscore[read1],
		read2, reverse_read(read2, RT -> num_seq), leftbestscore[read2], rightbestscore[read2]);
	printf("removed\n");
	getchar();
}
*/
				align -> cov = 1;
			}
			align = align -> next;
		}
	}
	free((void *) leftmost);
	free((void *) rightmost);

/*	Remove redundunt overlaps	*/

	n = 0;
	for(i = 0; i < RT -> num_seq * 2; i ++)	{
		align = eq_class[i];
		align0 = (ALIGN *) NULL;
		while(align)	{
			read1 = align -> reads[0];
			read2 = align -> reads[1];
			if(align -> cov == 1)	{
				align1 = align -> next;
				if(!align0)	{
					eq_class[i] = align1;
				} else	{
					align0 -> next = align1;
				}
				free_align(align);
				align = align1;
			} else	{
				n ++;
				align0 = align;
				align = align -> next;
			}
		}
	}
	free((void *) leftbestscore);
	free((void *) rightbestscore);
	return(n);
}

double perc_maxnuc(char *seq, int d1, int d2, double *perc)
{
	int	i, j, k, l;
	int	count[10];

	for(i = 0; i < 10; i ++)	{
		count[i] = 0;
	}
	for(i = d1; i <= d2; i ++)	{
		count[seq[i]] ++;
	}
	j = k = 0;
	for(i = 0; i < 10; i ++)	{
		if(count[i] > k)	{
			j = k;
			k = count[i];
		} else if(count[i] > j)	{
			j = count[i];
		}
	}
	if(d2 > d1)	{
		perc[0] = ((double) k) / (d2 - d1 + 1);
		perc[1] = ((double) j + k) / (d2 - d1 + 1);
	} else {
		perc[0] = perc[1] = 0;
	}
}

void caln_score(char *seq1, int len1, char *seq2, int len2, double *score, ALIGN *align)
{
	int	i, j, k, l, m, n;
	double	p, s;

	k = l = 0;
	s = 0;
	for(m = 0; m < align -> length - 1; m ++)	{
		i = align -> pos[0][m];
		j = align -> pos[1][m];
		while(i < align -> pos[0][m + 1] && j < align -> pos[1][m + 1])	{
			if(seq1[i] == seq2[j])	{
				s += 1.0;
				k ++;
			} else	{
				s -= 1.0;
				l ++;
			}
			i ++;
			j ++;
		}
		n = align -> pos[1][m + 1] - j + align -> pos[0][m + 1] - i;
		l += n;
		s -= n;
		
	}
	if(k + l > 0)	{
		score[1] = ((double) k) / (k + l);
	}
	score[0] = s;
/*
printf("s %f p %f k %d l %d\n", s, p, k, l);
*/
}

void update_bound(READINTERVAL *leftreadinterval, READINTERVAL *rightreadinterval, READINTERVAL **readinterval2, int p1, int p2)
{
	int	i, j, k, l;
	int	d;
	READINTERVAL	*readinterval;

	readinterval2[0] = leftreadinterval;
	readinterval2[1] = rightreadinterval;
	readinterval = leftreadinterval;
	while(readinterval)	{
		d = dist_range(p1, readinterval -> begin);
		if(d == 0)	{
			readinterval -> offset ++;
			if(p2 > readinterval -> length)	readinterval -> length = p2;
			break;
		}
		readinterval = readinterval -> next;
	}
	if(!readinterval)	{
		readinterval2[0] = (READINTERVAL *) ckalloc(1 * sizeof(READINTERVAL));
		readinterval2[0] -> begin = p1;
		readinterval2[0] -> offset = 1;
		readinterval2[0] -> length = p2;
		readinterval2[0] -> next = leftreadinterval;
	}
	readinterval = rightreadinterval;
	while(readinterval)	{
		d = dist_range(p2, readinterval -> begin);
		if(d == 0)	{
			readinterval -> offset ++;
			if(p1 < readinterval -> length)	readinterval -> length = p1;
			break;
		}
		readinterval = readinterval -> next;
	}
	if(!readinterval)	{
		readinterval2[1] = (READINTERVAL *) ckalloc(1 * sizeof(READINTERVAL));
		readinterval2[1] -> begin = p2;
		readinterval2[1] -> offset = 1;
		readinterval2[1] -> length = p1;
		readinterval2[1] -> next = rightreadinterval;
	}
}

int dist_range(int p1, int p2)
{
	int	d;

	d = p1 - p2;
	if(abs(d) < SHORT_D)	d = 0;
	return(d);
}
