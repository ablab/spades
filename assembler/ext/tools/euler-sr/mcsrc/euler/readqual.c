/***************************************************************************
 * Title:          readqual.c
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

char **src_name;

int readqualraw(char **val, char **src_seq, int *len, FILE *fp);
int readqual(char **val, char **src_seq, int *len, FILE *fp);
int readscore(char **val, int *len, FILE *fp);

int readqualraw(char **val, char **src_seq, int *len, FILE *fp)
{
	int	i, j, k, l, n, num_seq, s;
	char	str[800], temp[100];

	j = -1;
	n = 0;
	while(fgets(str, 490, fp))	{
		if(str[0] != '>')	{
			k = 0;
			do	{
				if(n >= len[j])	{
					break;
				}
				sscanf(&str[k], "%d", &s);
				if(src_seq[j][n] >= 4 || s > 90)	{
					val[j][n] = 0;
				} else 	{
					val[j][n] = s;
				}
				n ++;
				for(k ++; k < strlen(str); k ++)	{
					if(str[k - 1] == ' ')	{
						break;
					}
				}
			} while(k < strlen(str) - 1);
		} else	{
			if(j >= 0 && len[j] > n)	{
				len[j] = n;
			}
			j ++;
			n = 0;
		}
	}
	len[j] = n;
	num_seq = j + 1;

	return(num_seq);
}

int readqual(char **val, char **src_seq, int *len, FILE *fp)
{
	int	i, j, k, l, n, num_seq, s;
	char	str[800], temp[100];

	j = -1;
	n = 0;
	while(fgets(str, 490, fp))	{
		if(str[0] != '>')	{
			k = 0;
			do	{
				if(n >= len[j])	{
					break;
				}
				sscanf(&str[k], "%d", &s);
				if(src_seq[j][n] >= 4 || s > 90)	{
					val[j][n] = 0;
				} else if(s <= 20)	{
					val[j][n] = 0;
				} else if(s >= 30)	{
					val[j][n] = 2;
				} else 	{
					val[j][n] = 1;
				}
				n ++;
				for(k ++; k < strlen(str); k ++)	{
					if(str[k - 1] == ' ')	{
						break;
					}
				}
			} while(k < strlen(str) - 1);
		} else	{
			if(j >= 0 && len[j] > n)	{
				len[j] = n;
			}
			j ++;
			n = 0;
		}
	}
	len[j] = n;
	num_seq = j + 1;

	return(num_seq);
}

int readscore(char **val, int *len, FILE *fp)
{
	int	i, j, k, l, n, num_seq, s;
	char	str[800], temp[100];

	j = -1;
	n = 0;
	while(fgets(str, 490, fp))	{
		if(str[0] != '>')	{
			k = 0;
			do	{
				if(n >= len[j])	{
					break;
				}
				val[j][n ++] = str[k ++] - '0';
			} while(k < strlen(str) - 1);
		} else	{
			if(j >= 0 && len[j] > n)	{
				len[j] = n;
			}
			j ++;
			n = 0;
		}
	}
	len[j] = n;
	num_seq = j + 1;

	return(num_seq);
}
