/***************************************************************************
 * Title:          alignend.c
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

int trans_pos(char *seq, int len);
int chkhash(LINKPOS **hash, int index, char *seq, int len, int word_len, int *pair, int *mount);
LINKPOS *inserthash(LINKPOS *hash, int index);
int counthash(LINKPOS *hash, int index, int *pair, int *mount, int totp);
void buildhash(LINKPOS **hash, int index, char *seq, int len, int word_len);
char chk_pair(ALIGN **align, int read1, int read2, int num_seq);
int identity(char *seq1, int len1, char *seq2, int len2, int *sapp);
int alignend(char **src_name, char **src_seq, char **score, int *len_seq, int num_seq, int *endreads,
	     int *beginreads, int *num_endreads, ALIGN **eq_class);
int wholealign(char *seq1, int len_seq1, char *seq2, int len_seq2, int *sapp, int *cutp, int *seql);

int alignend(char **src_name, char **src_seq, char **score, int *len_seq, int num_seq, int *endreads,
	     int *beginreads, int *num_endreads, ALIGN **eq_class)
{
	int	i, j, k, l, m, n, i1, i2;
	int	read1, read2;
	char	c, c1;
	int	*pair, *mount;
	int	*sapp;
	LINKPOS 	**hash, *hashtmp;
	int	cutp[2], seql[2];

	sapp = (int *) ckalloc(2 * MAX_TMP_LEG * sizeof(int));
	pair = (int *) ckalloc(2 * num_endreads[0] * sizeof(int));
	mount = (int *) ckalloc(2 * num_endreads[0] * sizeof(int));
	hash = (LINKPOS **) ckalloc(n_ban * sizeof(LINKPOS *));

	n = 0;
	for(i = 0; i < num_endreads[0]; i ++)	{
		read1 = beginreads[i];
		buildhash(hash, read1, src_seq[read1], len_seq[read1], word_len);
	}
	for(j = 0; j < num_endreads[1]; j ++)	{
		read2 = endreads[j];
		k = chkhash(hash, read2, src_seq[read2], len_seq[read2], word_len, pair, mount);
		for(m = 0; m < k; m ++)	{
			read1 = pair[m];
			if(mount[m] < overlaplen - word_len)	continue;
			c = chk_pair(eq_class, read1, read2, num_seq);
			if(c)	{
				i1 = max(read1, read2);
				i2 = min(read1, read2);
				band = len_seq[i1] + len_seq[i2];
				c1 = overalign(src_seq[i1], src_seq[i2], score[i1], score[i2],
					len_seq[i1], len_seq[i2], 0, sapp, cutp, seql, i1, i2, num_seq);
				if(c1 >= 0)	{
					eq_class[i1] = new_align(cutp, seql, sapp, eq_class[i1], i1, i2, 1, c);
					
					n ++;
				}
			}
		}
	}
	free((void *) mount);
	free((void *) pair);
	free((void *) sapp);
	for(i = 0; i < n_ban; i ++)	{
		while(hash[i])	{
			hashtmp = hash[i] -> next;
			free((void *) hash[i]);
			hash[i] = hashtmp;
		}
	}
	free((void **) hash);
	return(n);
}

int trans_pos(char *seq, int len)
{
	int	i, j, k;

	k = seq[0];
	for(i = 1; i < len; i ++)	{
		if(seq[i] > 4)	{
			k = k * 4;
		} else	{
			k = k * 4 + seq[i];
		}
	}

	return(k);
}

void buildhash(LINKPOS **hash, int index, char *seq, int len, int word_len)
{
	int	i, j, k;

	if(len < overlaplen)	return;
	k = trans_pos(seq, word_len);
	hash[k] = inserthash(hash[k], index);
	for(i = 1; i < len - overlaplen + 1; i ++)	{
		k = (k * 4 + seq[i + word_len - 1]) % n_ban;
		hash[k] = inserthash(hash[k], index);
	}
}

LINKPOS *inserthash(LINKPOS *hash, int index)
{
	int	i, j, k;
	LINKPOS	*hashtmp;

	hashtmp = hash;
	while(hashtmp)	{
		if(hashtmp -> readindex == index)	{
			hashtmp -> position ++; 
			return(hash);
		}
		hashtmp = hashtmp -> next;
	}
	hashtmp = (LINKPOS *) ckalloc(1 * sizeof(LINKPOS));
	hashtmp -> readindex = index;
	hashtmp -> position = 1;
	hashtmp -> next = hash;
	return(hashtmp);
}

int chkhash(LINKPOS **hash, int index, char *seq, int len, int word_len, int *pair, int *mount)
{
	int	i, j, k, totp;

	totp = 0;
	if(len < word_len)	return(0);
	k = trans_pos(seq, word_len);
	totp = counthash(hash[k], index, pair, mount, totp);
	for(i = 1; i < len - word_len + 1; i ++)	{
		k = (k * 4 + seq[i + word_len - 1]) % n_ban;
		totp = counthash(hash[k], index, pair, mount, totp);
	}
	return(totp);
}

int counthash(LINKPOS *hash, int index, int *pair, int *mount, int totp)
{
	int	i, j, k;
	LINKPOS	*hashtmp;

	hashtmp = hash;
	while(hashtmp)	{
		if(hashtmp -> readindex != index)	{
			for(i = 0; i < totp; i ++)	{
				if(pair[i] == hashtmp -> readindex)	{
					mount[i] += hashtmp -> position;
					break;
				}
			}
			if(i == totp)	{
				pair[totp] = hashtmp -> readindex;
				mount[totp ++] = 1;
			}
		}
		hashtmp = hashtmp -> next;
	}
	return(totp);
}

char chk_pair(ALIGN **align, int read1, int read2, int num_seq)
{
	int	i, j, k, m1, m2;
	ALIGN	*aln;

	if(read1 == read2 || read1 == reverse_read(read2, num_seq))	{
		return(0);
	}

	aln = align[read1];
	while(aln)	{
		if(aln -> reads[1] == read2)	{
			return(0);
		}
		aln = aln -> next;
	}
	aln = align[read2];
	while(aln)	{
		if(aln -> reads[1] == read1)	{
			return(0);
		}
		aln = aln -> next;
	}
	read1 = reverse_read(read1, num_seq);
	read2 = reverse_read(read2, num_seq);
	aln = align[read1];
	while(aln)	{
		if(aln -> reads[1] == read2)	{
			return(0);
		}
		aln = aln -> next;
	}
	aln = align[read2];
	while(aln)	{
		if(aln -> reads[1] == read1)	{
			return(0);
		}
		aln = aln -> next;
	}
	return(1);
}

int wholealign(char *seq1, int len_seq1, char *seq2, int len_seq2, int *sapp, int *cutp, int *seql)
{
	int	i, j, k, l, e1, e2;
	int	st, ln;
	double	p;

/*	Use regular alignment	*/

	ln = max(len_seq1, len_seq2);
	k = LOCAL_ALIGN0(&seq1[-1], &seq2[-1], len_seq1, len_seq2, -ln, ln,
		 W, g, h, &cutp[0], &cutp[1], &e1, &e2, len_seq1 + len_seq2);
	if(k < 0)	{
		return(-1);
	}
	seql[0] = e1 - cutp[0] + 1;
	seql[1] = e2 - cutp[1] + 1;
	cutp[0] --;
	cutp[1] --;
	if(seql[0] < MIN_OVERLAP2 || seql[1] < MIN_OVERLAP2)	{
		return(-1);
	}

/*	Use regular alignment	*/

	k = ALIGN0(&seq1[cutp[0] - 1], &seq2[cutp[1] - 1], 
		seql[0], seql[1], -band, band, W, g, h, sapp, len_seq1 + len_seq2,
		len_seq1 + len_seq2);
	k = identity(&seq1[cutp[0]], seql[0], &seq2[cutp[1]], seql[1], sapp);
	p = 1 - ((double) k) / min(seql[0], seql[1]);
	if(p >= MIN_IDENTITY)	{
		return(k);
	} else {
		return(-1);
	}
}

int identity(char *seq1, int len1, char *seq2, int len2, int *sapp)
{
	int	i, j, k, l, m, n;
	int	op;
	double	p, s;

	i = j = k = l = 0;
	while(i < len1 || j < len2)	{
		op = *sapp ++;
		if(op == 0)	{
			if(seq1[i] != seq1[j])	{
				k ++;
			}
			i ++;
			j ++;
		} else if(op > 0)	{
			k += op;
			j += op;
		} else 	{
			k -= op;
			i -= op;
		}
	}
	return(k);
}
