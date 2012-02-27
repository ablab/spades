/***************************************************************************
 * Title:          readseq1by1.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>
#include <extfunc.h>
#include <extvab.h>

#define MAXNUM 100000

void inputreadname(READTABLE *rt, char *inpfile);
void inputread(READTABLE *rt, char *inpfile);
int readseqname(char **src_name, FILE *fp);
int readseqpar(int *max_leg, int *max_name_leg, FILE *fp);
int readseq1by1(char **src_seq, char **src_name, int *len_seq, FILE *fp);
int readseq1by1gen(char **src_seq, char **src_name, int *len_seq, FILE *fp);
char **init_score(READTABLE *RT);
void read_score_file(char *qualfile, char **score, READTABLE *RT);


void inputreadname(READTABLE *rt, char *inpfile)
{
	int	i;
	int	max_leg, max_name_leg;
	FILE	*fp;

	fp = ckopen(inpfile, "r");
	rt -> num_seq = readseqpar(&max_leg, &max_name_leg, fp);
	fclose(fp);
	rt -> src_name = (char **) ckalloc(rt -> num_seq * sizeof(char *));
	for(i = 0; i < rt -> num_seq; i ++)	{
		rt -> src_name[i] = (char *) ckalloc((max_name_leg + 1)* sizeof(char));
	}
	fp = ckopen(inpfile, "r");
	rt -> num_seq = readseqname(rt -> src_name, fp);
	fclose(fp);
}

void inputread(READTABLE *rt, char *inpfile)
{
	int	i;
	int	max_leg, max_name_leg;
	FILE	*fp;

	fp = ckopen(inpfile, "r");
	rt -> num_seq = readseqpar(&max_leg, &max_name_leg, fp);
	fclose(fp);
	rt -> src_name = (char **) ckalloc(rt -> num_seq * sizeof(char *));
	rt -> src_seq = (char **) ckalloc(rt -> num_seq * sizeof(char *));
	for(i = 0; i < rt -> num_seq; i ++)	{
		rt -> src_name[i] = (char *) ckalloc((max_name_leg + 1) * sizeof(char));
	}
	fp = ckopen(inpfile, "r");
	rt -> num_seq = readseq1by1(rt -> src_seq, rt -> src_name, rt -> len_seq, fp);
	fclose(fp);
}

int readseqname(char **src_name, FILE *fp)
{
	int	i, j, k, l, n;
	char	*seq, c;
	char	str[5000];

	k = 0;
	while(fgets(str, 4950, fp))	{
		if(str[0] == '>')	{
			sscanf(&str[1], "%s", src_name[k]);
			k ++;
		}
	}
	return(k);
}

int readseqpar(int *max_leg, int *max_name_leg, FILE *fp)
{
	int	i, j, k, l, n;
	char	*seq, c;
	char	str[5000], src_name[5000];

	*max_name_leg = *max_leg = 1;
	n = 0;
	k = -1;
	while(fgets(str, 4950, fp))	{
		if(str[0] == '>')	{
			if(k >= 0)	{
				if(n > *max_leg)	{
					*max_leg = n;
				}
			}
			n = 0;
			k ++;
			sscanf(&str[1], "%s", src_name);
			if((l = strlen(src_name)) > *max_name_leg)		*max_name_leg = l;
		} else {
			n += strlen(str);
		}
	}
	if(n > *max_leg)	{
		*max_leg = n;
	}
	k ++;
	return(k);
}

int readseq1by1(char **src_seq, char **src_name, int *len_seq, FILE *fp)
{
	int	i, j, k, l, n;
	char	*seq, c;
	char	str[5000];

	seq = (char *) ckalloc(MAXNUM * sizeof(char));

	n = 0;
	k = -1;
	while(fgets(str, 4950, fp))	{
		if(str[0] == '#')	continue;
		if(str[0] == '>')	{
			if(k >= 0)	{
				len_seq[k] = n;
				src_seq[k] = (char *) ckalloc((n + 1) * sizeof(char));
				for(i = 0; i < n; i ++)		src_seq[k][i] = seq[i];
			}
			n = 0;
			k ++;
			sscanf(&str[1], "%s", src_name[k]);
		} else {
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					c = char2int(str[i]);
					seq[n ++] = c;
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					c = char2int(str[i] - 'A' + 'a');
					seq[n ++] = c;
				}
			}
		}
	}

	if(k >= 0)	{
		len_seq[k] = n;
		src_seq[k] = (char *) ckalloc((n + 1) * sizeof(char));
		for(i = 0; i < n; i ++)		src_seq[k][i] = seq[i];
	}
	k ++;

	free((void *) seq);
	return(k);
}

int readseq1by1gen(char **src_seq, char **src_name, int *len_seq, FILE *fp)
{
	int	i, j, k, l, n;
	char	*seq, c;
	char	str[500];

	seq = (char *) ckalloc(MAXNUM * sizeof(char));

	n = 0;
	k = -1;
	while(fgets(str, 450, fp))	{
		if(str[0] == '#')	continue;
		if(str[0] == '>')	{
			if(k >= 0)	{
				len_seq[k] = n;
				src_seq[k] = (char *) ckalloc((n + 1) * sizeof(char));
				for(i = 0; i < n; i ++)		src_seq[k][i] = seq[i];
			}
			n = 0;
			k ++;
			sscanf(&str[1], "%s", src_name[k]);
		} else {
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] >= 'a' && str[i] <= 'z') {
					c = char2intgen(str[i]);
					seq[n ++] = c;
				} else if(str[i] >= 'A' && str[i] <= 'Z') {
					c = char2intgen(str[i] - 'A' + 'a');
					seq[n ++] = c;
				}
			}
		}
	}

	if(k >= 0)	{
		len_seq[k] = n;
		src_seq[k] = (char *) ckalloc((n + 1) * sizeof(char));
		for(i = 0; i < n; i ++)		src_seq[k][i] = seq[i];
	}
	k ++;

	free((void *) seq);
	return(k);
}

char **init_score(READTABLE *RT)
{
	int	i, j;
	char	**score;

	score = (char **) ckalloc(2 * RT -> num_seq * sizeof(char *));
	for(i = 0; i < RT -> num_seq; i ++)	{
		score[i] = (char *) ckalloc(RT -> len_seq[i] * sizeof(char));
		for(j = 0; j < RT -> len_seq[i]; j ++)	{
			score[i][j] = 1;
		}
	}
	for(i = RT -> num_seq; i < RT -> num_seq * 2; i ++)	{
		score[i] = (char *) ckalloc(RT -> len_seq[i] * sizeof(char));
		for(j = 0; j < RT -> len_seq[i]; j ++)	{
			score[i][j] = score[i - RT -> num_seq][RT -> len_seq[i] - 1 - j];
		}
	}
	return(score);
}

void read_score_file(char *qualfile, char **score, READTABLE *RT)
{
	int	i, j, n;
	FILE	*fp;

	fp = ckopen(qualfile, "r");
	n = readscore(score, RT -> len_seq, fp);
	if(n != RT -> num_seq)	{
		printf("The quality file (%d sequence) and the reads (%d sequence) file are not consistent!\n", n,
			RT -> num_seq);
		exit(0);
	}
	for(i = RT -> num_seq; i < RT -> num_seq * 2; i ++)	{
		for(j = 0; j < RT -> len_seq[i]; j ++)	{
			score[i][j] = score[i - RT -> num_seq][RT -> len_seq[i] - 1 - j];
		}
	}
	fclose(fp);
}

/****************************************************************************
 * Read fasta file
 * Create reverse reads too
 *
 * TODO:
 * 1. This assumes a fixed size for the read names, a fixed max. # reads,
 *    etc.; adapt it to figure out these values itself.
 *    --by Haixu, it is done. 4/9/2003
 * 2. Infer compressed versions of the data:
 *    just the read lengths,
 *    just the reads but not the complements,
 *    2-bits per character,
 *    etc., and only do the ones needed for a particular program.
 ****************************************************************************/

void read_fasta_file(char *seqfile, READTABLE *RT)
{
  int i, j;
  int num_seq;
  int max_leg, max_name_leg;

  FILE *fp;

  fp = ckopen(seqfile, "r");
  num_seq = readseqpar(&max_leg, &max_name_leg, fp);
  fclose(fp);
  /**********************************************************************
   * Input the reads, their lengths, and names
   **********************************************************************/
// added by haixu
  RT->new_len_seq = (int *) ckalloc(2 * num_seq * sizeof(int));
  RT->len_seq = (int *) ckalloc(2 * num_seq * sizeof(int));
  RT->src_seq = (char **) ckalloc(2 * num_seq * sizeof(char *));
  RT->src_name = (char **) ckalloc(num_seq * sizeof(char *));
  for (i = 0; i < num_seq; i ++) {
    RT->src_name[i] = (char *) ckalloc((max_name_leg + 1) * sizeof(char));
  }
  fp = ckopen(seqfile, "r");
  RT->num_seq = readseq1by1(RT->src_seq, RT->src_name, RT->len_seq,
				      fp);
  fclose(fp);
  printf("# reads input: %d\n", num_seq);

/* set up max_leg, i.e. the maximal length of all the reads
 * and max_seq i.e. the total number of reads.
 * Both are external variables.	*/

//  max_seq = num_seq;

  /**********************************************************************
   * Make reverse complements of input sequences rev(i) --> i + num_seq
   **********************************************************************/

  for (i = 0; i < num_seq; i ++) {
    RT->len_seq[num_seq + i] = RT->len_seq[i];
    RT->src_seq[num_seq + i] = (char *) ckalloc(RT->len_seq[i] * sizeof(char));
    for (j = 0; j < RT->len_seq[i]; j ++) {
      RT->src_seq[num_seq + i][j] =
	rev(RT->src_seq[i][RT->len_seq[i] - j - 1]);
    }
  }
}

void free_readlist(READTABLE *RT)
{
	int i,j;
	
	for(i = 0; i < 2 * RT -> num_seq; i ++)	{
		free((void *) RT -> readlist[i]);
	}
	free((void **) RT -> readlist);
}

void free_readtable(READTABLE *RT)
{
  int i;

  for (i = 0; i < 2 * RT->num_seq; i ++) {
    free((void *) RT->src_seq[i]);
  }
  for (i = 0; i < RT -> num_seq; i ++) {
    free((void *) RT->src_name[i]);
  }
  free((void **) RT->src_seq);
  free((void **) RT->src_name);
  free((void *) RT->len_seq);

}
