/***************************************************************************
 * Title:          detectoverlap.c
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

long trans_seq(char *seq, int len);
LINKPOS *ins_linkpos(LINKPOS *linkpos, int i1, int i2);
HASH *free_hash(HASH *hash);
HASH *ins_hash(HASH *hash, int i1, int i2, char **src_seq, int num_seq);
READINTERVAL *ins_readinterval(READINTERVAL *readinterval, int readindex, int offset);
int detectoverlap(char **src_seq, char **score, int *len_seq, int num_seq, READINTERVAL **readinterval);
int comp_word(char *w1, char *w2, int len);
READINTERVAL *new_readinterval(HASH *hash, int i1, int i2, char **src_seq, READINTERVAL *readinterval, int num_seq);
READINTERVAL *chk_redund(READINTERVAL *readinterval);

int detectoverlap(char **src_seq, char **score, int *len_seq, int num_seq, READINTERVAL **readinterval)
{
	int	i, j, k, l;
	int	num;
	long    hash_ban;
	READINTERVAL	*readinterval0, *readinterval1;
	TABLE	*hash_table;

	hash_table = (TABLE *) ckalloc(n_ban * sizeof(TABLE));
	for(i = 0; i < num_seq; i ++)	{
		if(len_seq[i] < overlaplen)     continue;
		hash_ban = trans_seq(&src_seq[i][0], word_len);
		readinterval[i] = new_readinterval(hash_table[hash_ban].prev, i, 0, src_seq, readinterval[i], num_seq);
		hash_table[hash_ban].prev = ins_hash(hash_table[hash_ban].prev, i, 0, src_seq, num_seq);
		for(j = 1; j <= len_seq[i] - word_len; j ++)	{
			hash_ban = (hash_ban * 4 + src_seq[i][j + word_len - 1]) % n_ban;
			readinterval[i] = new_readinterval(hash_table[hash_ban].prev, i, j, src_seq, readinterval[i], num_seq);
			hash_table[hash_ban].prev = ins_hash(hash_table[hash_ban].prev, i, j,
						 src_seq, num_seq);
		}
		readinterval[i] = chk_redund(readinterval[i]);
	}

	for(i = num_seq; i < 2 * num_seq; i ++)	{
		if(len_seq[i] < overlaplen)     continue;
		hash_ban = trans_seq(&src_seq[i][0], word_len);
		readinterval[i] = new_readinterval(hash_table[hash_ban].prev, i, 0, src_seq, readinterval[i], num_seq);
		for(j = 1; j <= len_seq[i] - word_len; j ++)	{
			hash_ban = (hash_ban * 4 + src_seq[i][j + word_len - 1]) % n_ban;
			readinterval[i] = new_readinterval(hash_table[hash_ban].prev, i, j, src_seq, readinterval[i], num_seq);
		}
		readinterval[i] = chk_redund(readinterval[i]);
	}
	for(i = 0; i < n_ban; i ++)	{
		while(hash_table[i].prev)	{
			hash_table[i].prev = free_hash(hash_table[i].prev);
		}
	}
	free((void *) hash_table);
	num = 0;
	for(i = 0; i < num_seq * 2; i ++)	{
		num += size_readinterval(readinterval[i]);
	}
	return(num);
}

READINTERVAL *chk_redund(READINTERVAL *readinterval)
{
	int	i, j, k, l;
	READINTERVAL	*readinterval0, *readinterval1, *readinterval2;

	readinterval0 = readinterval;
	while(readinterval0)	{
		readinterval1 = readinterval0 -> next;
		while(readinterval1)	{
			if(readinterval0 -> eq_read == readinterval1 -> eq_read)	{
				if(readinterval0 -> length >= readinterval1 -> length)	{
					readinterval1 -> cov = 1;
				} else	{
					readinterval0 -> cov = 1;
				}
			}
			readinterval1 = readinterval1 -> next;
		}
		readinterval0 = readinterval0 -> next;
	}
	readinterval2 = (READINTERVAL *) NULL;
	readinterval0 = readinterval;
	while(readinterval0)	{
		if(readinterval0 -> cov == 1)	{
			readinterval1 = readinterval0 -> next;
			free((void *) readinterval0);
			readinterval0 = readinterval1;
			if(readinterval2)	{
				readinterval2 -> next = readinterval1;
			} else	{
				readinterval = readinterval1;
			}
		} else	{
			readinterval2 = readinterval0;
			readinterval0 = readinterval0 -> next;
		}
	}
	return(readinterval);
}

HASH *free_hash(HASH *hash)
{
	HASH	*hash0;
	LINKPOS	*linkpos;

	hash0 = hash -> next;
	while(hash -> linkpos)	{
		linkpos = hash -> linkpos -> next;
		free((void *) hash -> linkpos);
		hash -> linkpos = linkpos;
	}
	free((void *) hash);
	return(hash0);
}

long trans_seq(char *seq, int len)
{
        int     i, j, k;
        long    res;

        res = 0;
        for(i = 0; i < len; i ++)       {
                res = res * 4 + seq[i];
        }

        return(res);
}

READINTERVAL *new_readinterval(HASH *hash, int i1, int i2, char **src_seq, READINTERVAL *readinterval, int num_seq)
{
        int     i, j, c;
        HASH    *h0;
	LINKPOS *linkpos;

        h0 = hash;
        while(h0)       {
		linkpos = h0 -> linkpos;
		while(linkpos)	{
			if(i1 >= num_seq && i1 - num_seq < linkpos -> readindex)	{
				linkpos = linkpos -> next;
				continue;
			}
                	c = comp_word(&src_seq[linkpos -> readindex][linkpos -> position + word_len],
			      &src_seq[i1][i2 + word_len], overlaplen - word_len);
			if(c && linkpos -> readindex != i1 && linkpos -> readindex !=
				 reverse_read(i1, num_seq))	{
				readinterval = ins_readinterval(readinterval, linkpos -> readindex,
					 i2 - linkpos -> position);
			}
			linkpos = linkpos -> next;
		}
                h0 = h0 -> next;
        }

        return(readinterval);
}

HASH *ins_hash(HASH *hash, int i1, int i2, char **src_seq, int num_seq)
{
        int     i, j, c;
        HASH    *h0;
	LINKPOS *linkpos;

        h0 = hash;
        while(h0)       {
		linkpos = h0 -> linkpos;
                c = comp_word(&src_seq[linkpos -> readindex][linkpos -> position + word_len],
			      &src_seq[i1][i2 + word_len], overlaplen - word_len);
                if(c)   {
			h0 -> linkpos = ins_linkpos(linkpos, i1, i2);
			return(hash);
		}
                h0 = h0 -> next;
        }

        h0 = (HASH *) ckalloc(1 * sizeof(HASH));
	h0 -> linkpos = (LINKPOS *) ckalloc(1 * sizeof(LINKPOS));
	h0 -> linkpos -> readindex = i1;
	h0 -> linkpos -> position = i2;
	h0 -> next = hash;
        return(h0);
}

LINKPOS *ins_linkpos(LINKPOS *linkpos, int i1, int i2)
{
	LINKPOS *newlinkpos;

	newlinkpos = (LINKPOS *) ckalloc(1 * sizeof(LINKPOS));
	newlinkpos -> readindex = i1;
	newlinkpos -> position = i2;
	newlinkpos -> next = linkpos;
	return(newlinkpos);
}

READINTERVAL *ins_readinterval(READINTERVAL *readinterval, int readindex, int offset)
{
	READINTERVAL *newreadinterval, *readinterval1;

	readinterval1 = readinterval;
	while(readinterval1)	{
		if(readinterval1 -> eq_read == readindex && abs(readinterval1 -> offset - offset) < band / 2)
			break;
		readinterval1 = readinterval1 -> next;
	}

	if(!readinterval1)	{
		newreadinterval = (READINTERVAL *) ckalloc(1 * sizeof(READINTERVAL));
		newreadinterval -> eq_read = readindex;
		newreadinterval -> offset = offset;
		newreadinterval -> next = readinterval;
		newreadinterval -> length = 1;
		return(newreadinterval);
	} else	{
		readinterval1 -> length ++;
		return(readinterval);
	}
}

int comp_word(char *w1, char *w2, int len)
{
	int	i, j, k, j1, j2;

	for(i = 0; i < len; i ++)	{
		if(w1[i] != w2[i])	{
			return(0);
		}
	}

	return(1);
}
