/***************************************************************************
 * Title:          chg_reads.c
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

int defineshift(char *newseq, int len_newseq, char *readseq, int len_readseq, int begin);
int defineshift_rev(char *newseq, int len_newseq, char *readseq, int len_readseq, int begin, int len);
int findri(READTABLE *RT, READINTERVAL *readinterval, int num_intv, char *seq, int len, char *label);
int acculen(char *seq, int len);
void chg_intv(READINTERVAL *readinterval, char **newseq, int *len_newseq, READTABLE *RT);
int compseq_gap(char *seq1, int len1, char *seq2, int len2);
void chg_reads(READTABLE *RT, READINTERVAL *readinterval, int num_intv, int startpos, char *readseq, 
	       int len_readseq, char **newseq, int *len_newseq, int index, char *label);

void chg_reads(READTABLE *RT, READINTERVAL *readinterval, int num_intv, int startpos, char *readseq, 
	       int len_readseq, char **newseq, int *len_newseq, int index, char *label)
{
	int	i, j, k, l, m, n;

	k = findri(RT, readinterval, num_intv, readseq, len_readseq, label);
	label[k] = 1;
	j = readinterval[k].eq_read;
	if(j < RT -> num_seq)	{
		len_newseq[j] = defineshift(newseq[j], len_newseq[j],
			 readseq, len_readseq, readinterval[k].begin);
	} else	{
		len_newseq[j - RT -> num_seq] = defineshift_rev(newseq[j - RT -> num_seq], len_newseq[j - RT -> num_seq],
			 readseq, len_readseq, readinterval[k].begin, RT -> len_seq[j - RT -> num_seq]);
	}
	readinterval[k].offset = startpos;
}

int findri(READTABLE *RT, READINTERVAL *readinterval, int num_intv, char *seq, int len, char *label)
{
	int	i, j, k, l, c, sml;

	for(i = 0; i < num_intv; i ++)	{
		if(label[i] == 1)	continue;
		k = readinterval[i].eq_read;
		l = readinterval[i].begin;
		sml = min(RT -> len_seq[k] - l, readinterval[i].length);
		c = compseq_gap(&(RT -> src_seq[k][l]), sml, seq, len);
		if(c)	{
			return(i);
		}
	}
	printf("Sequence not found. \n");
	exit(-1);
}

int compseq_gap(char *seq1, int len1, char *seq2, int len2)
{
	int i, j, k;

	j = 0;
	for(i = 0; i < len2; i ++)	{
		if(seq2[i] == CODE_GAP)	{
			continue;
		} else if(seq2[i] != seq1[j])	{
			return(0);
		} else	{
			j ++;
		}
	}
	if(j == len1)	{
		return(1);
	} else	{
		return(0);
	}
}

int defineshift(char *newseq, int len_newseq, char *readseq, int len_readseq, int begin)
{
	int	i, j, k, l;

	i = 0;
	k = begin - 1;
	while(i < len_readseq)	{
		if(readseq[i] == CODE_GAP)	{
			newseq[k] ++;
			len_newseq ++;
		} else	{
			k ++;
		}
		i ++;
	}
	return(len_newseq);
}

int defineshift_rev(char *newseq, int len_newseq, char *readseq, int len_readseq, int begin, int len)
{
	int	i, j, k, l;

	i = 0;
	k = begin;
	while(i < len_readseq)	{
		if(readseq[i] == CODE_GAP)	{
			newseq[len - 1 - k] ++;
			len_newseq ++;
		} else	{
			k ++;
		}
		i ++;
	}
	return(len_newseq);
}

void chg_intv(READINTERVAL *readinterval, char **newseq, int *len_newseq, READTABLE *RT)
{
	int	i, j, k, l, m;

	if(readinterval -> eq_read < RT -> num_seq)	{
		j = acculen(newseq[readinterval -> eq_read], readinterval -> begin);
		k = acculen(newseq[readinterval -> eq_read], readinterval -> length + readinterval -> begin - 1);
		readinterval -> begin = j;
		readinterval -> length = k - j + 1;
	} else	{
		m = readinterval -> eq_read - RT -> num_seq;
		l = RT -> len_seq[readinterval -> eq_read];
		j = acculen(newseq[m], l - 1 - readinterval -> begin);
		k = acculen(newseq[m], l - readinterval -> length - readinterval -> begin);
		readinterval -> begin = len_newseq[m] - j - 1;
		readinterval -> length = j - k + 1;
	}
}

int acculen(char *seq, int len)
{
	int	i, n;

	n = len;
	for(i = 0; i < len; i ++)	{
		n += seq[i];
	}
	return(n);
}
