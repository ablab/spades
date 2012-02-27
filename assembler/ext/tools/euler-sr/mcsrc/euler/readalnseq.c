/***************************************************************************
 * Title:          readalnseq.c
 * Author:         Haixu Tang
 * Created:        Jun. 2004
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extvab.h>
#include <extfunc.h>

#define SEC_LEN 80
#define MAX_LEG 500

int readalnseq(FILE *fp, char **readseq, int *len_readseq, int *startpos, int num_intv, int length,
		char *contigseq, char *contigqual, int *len_contigseq);
int process_row(int *startpos, char **readseq, int *len_readseq, int num_seq, char *row, int len_row);
int addrow(char *readseq, int len1, char *str, int len2);
int addseq(char *readseq, char *str, int pos, int length);
int nextletter(char *row, int pos, int len);
char comp_qual(int n1, int n2);

int readalnseq(FILE *fp, char **readseq, int *len_readseq, int *startpos, int num_intv, int length,
		char *contigseq, char *contigqual, int *len_contigseq)
{
	int	i, j, k, l, m, num_sec, ln, lc, n1, n2;
	int	num_seq;
	int	*len_row;
	int	lastrow, currrow, tot_row, nonl;
	char	**rowstr;
	char	str[200];

	lc = tot_row = 0;
	lastrow = -1;
	nonl = num_sec = currrow = 0;
	len_row = (int *) ckalloc(num_intv * sizeof(int));
	rowstr = (char **) ckalloc(num_intv * sizeof(char *));
	while(fgets(str, 190, fp))	{
		if(!strncmp(str, "-------------------------------------------------------------------", 50))	{
			if(!fgets(str, 190, fp) || !strncmp(str, "%%%", 3))	{
				break;
			}
			if(tot_row < currrow)		{
				tot_row = currrow;
			}
			ln = strlen(str);
			lc = addrow(contigseq, lc, str, ln);
			for(m = currrow; m < tot_row; m ++)	{
				len_row[m] += SEC_LEN;
			}
			lastrow = currrow;
/*	Skip 3 empty lines	*/
			fgets(str, 190, fp);
/*			if(str[0] == '\n')	continue; */
			fgets(str, 190, fp);
			fgets(str, 190, fp);
			nonl = currrow = 0;
			num_sec ++;
		} else 	{
			if(len_row[currrow] == 0)	{
				rowstr[currrow] = (char *) ckalloc((length + MAX_LEG) * sizeof(char));
				len_row[currrow] = num_sec * SEC_LEN;
			}
			ln = strlen(str);
			len_row[currrow] = addrow(rowstr[currrow], len_row[currrow], str, ln);
			currrow ++;
		}
	}

/*	Derive contig quality	*/

	for(i = 0; i < lc; i ++)	{
		if(contigseq[i] >= 'A' && contigseq[i] <= 'Z')	{
			contigseq[i] += 'a' - 'A';
		}
		n1 = n2 = 0;
		if(contigseq[i] == '-')	{
			contigseq[i] = CODE_GAP;
		} else	{
			for(j = 0; j < tot_row; j ++)	{
				if(i < len_row[j] && rowstr[j][i] != ' ' && rowstr[j][i] != 0)	{
					if(rowstr[j][i] == contigseq[i])	{
						n1 ++;
					} else	{
						n2 ++;
					}
				}
			}
			contigqual[i] = comp_qual(n1, n2);
			contigseq[i] = char2int(contigseq[i]);
		}
	}

/*	Extract reads from rows	*/

	num_seq = 0;
	for(i = 0; i < tot_row; i ++)	{
		num_seq = process_row(startpos, readseq, len_readseq, num_seq, rowstr[i], len_row[i]);
		free((void *) rowstr[i]);
	}
	free((void **) rowstr);
	free((void *) len_row);
	*len_contigseq = lc;
	return(num_seq);
}

char comp_qual(int n1, int n2)
{
	int	i, j, n;
	n = n1 + n2;
	j = n1 - n2;
	if(j <= 0)	{
		return(0);
	} else	{
		if(j <= 6)	{
			return(10 * j);
		} else	{
			return(60);
		}
	}
}

int process_row(int *startpos, char **readseq, int *len_readseq, int num_seq, char *row, int len_row)
{
	int	i, j, k, l;

/*
for(i = 0; i < len_row; i ++)	{
	printf("%c", row[i]);
}
printf("\n");
*/

	i = 0;
	while(i < len_row)	{
		i = nextletter(row, i, len_row); 
		if(i < len_row)	{
			startpos[num_seq] = i;
			len_readseq[num_seq] = addseq(readseq[num_seq], row, i, len_row);
			i += len_readseq[num_seq];
			num_seq ++;
			i = nextletter(row, i, len_row); 
		}
	}
	return(num_seq);
}

int nextletter(char *row, int pos, int len)
{
	int	i, j;
	for(i = pos; i < len; i ++)	{
		if(row[i] >= 'a' && row[i] <= 'z' || row[i] == '-')	{
			break;
		}
	}
	return(i);
}

int addrow(char *readseq, int len1, char *str, int len2)
{
	int	i, j, k;

	for(i = 0; i < len2; i ++)	{
		if((str[i] >= 'A' && str[i] <= 'Z') ||
		   (str[i] >= 'a' && str[i] <= 'z') ||
		   str[i] == ' ' || str[i] == '-')	{
			readseq[len1 ++] = str[i];
		}
	}
	return(len1);
}

int addseq(char *readseq, char *str, int pos, int length)
{
	int	i, j, k;

	k = 0;
	for(i = pos; i < length; i ++)	{
		if(str[i] >= 'a' && str[i] <= 'z')	{
			readseq[k ++] = char2int(str[i]);
		} else if(str[i] == '-')	{
			readseq[k ++] = CODE_GAP;
		} else	{
			return(k);
		}
	}
	return(k);
}
