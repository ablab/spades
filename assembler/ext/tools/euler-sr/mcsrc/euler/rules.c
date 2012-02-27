/***************************************************************************
 * Title:          rules.c
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

extern int overlaplen;

int readrules(FILE *fp);
int readrules2(FILE *fp, MATEPAIRRULES *MPR);

int readrules(FILE *fp)
{
  MATEPAIRRULES MPR_mem, *MPR = &MPR_mem;
  int ntypepair1 = readrules2(fp,MPR);

  /* populate global variables */
  /* TODO: Convert all calling routines to use MATEPAIRRULES structure
   * instead
   */

  memcpy(&librule, &MPR->librule, sizeof(MPR->librule));
  memcpy(&namerule, &MPR->namerule, sizeof(MPR->namerule));
  memcpy(&pairrule, &MPR->pairrule, sizeof(MPR->pairrule));
  memcpy(&ntyperule, &MPR->ntyperule, sizeof(MPR->ntyperule));
  memcpy(&pairrange, &MPR->pairrange, sizeof(MPR->pairrange));
  memcpy(&platerule, &MPR->platerule, sizeof(MPR->platerule));
  memcpy(&ntypepair, &MPR->ntypepair, sizeof(MPR->ntypepair));
  return ntypepair1;
}

/* TODO:
 * 1. Convert all pgms to use readrules2 instead of global vars + readrules.
 * 2. This uses fixed size allocations; change to adapt to correct sizes.
 */

int readrules2(FILE *fp, MATEPAIRRULES *MPR)
{
	int	i, j, k, l;
	char	temp[100];
	char	c, s1[10], s2[10], str[500];

	MPR->ntyperule[0] = MPR->ntyperule[1] = 0;
	MPR->ntypepair = 0;
	while(fgets(str, 190, fp))	{
		if(!strncmp(str, "Finish", 6))	{
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] == ':')	{
					break;
				}
			}
			for(i ++; i < strlen(str); i ++)	{
				if((c = str[i]) != ' ')	{
					MPR->namerule[0][MPR->ntyperule[0] ++] = c;
				}
			}
		} else if(!strncmp(str, "Single", 6))	{
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] == ':')	{
					break;
				}
			}
			for(i ++; i < strlen(str); i ++)	{
				if((c = str[i]) != ' ')	{
					MPR->namerule[1][MPR->ntyperule[1] ++] = c;
				}
			}
		} else if(!strncmp(str, "Double", 6))	{
			for(i = 0; i < strlen(str); i ++)	{
				if(str[i] == ':')	{
					break;
				}
			}
			sscanf(str + i + 1, "%s%s%s%s%d%d", s1, s2, MPR->librule[MPR->ntypepair], temp, &MPR->pairrange[MPR->ntypepair][0], &MPR->pairrange[MPR->ntypepair][1]);
			for(i = 0; i < strlen(temp); i ++)	{
				if(temp[i] == '-')	temp[i] = ' ';
			}
			if(!strcmp(temp, "ALL"))	{
				MPR->platerule[MPR->ntypepair][0] = 0;
				MPR->platerule[MPR->ntypepair][1] = 100000000;
			} else	{
				sscanf(temp, "%d%d", &MPR->platerule[MPR->ntypepair][0], &MPR->platerule[MPR->ntypepair][1]);
			}
			MPR->pairrule[MPR->ntypepair][0] = s1[0];
			MPR->pairrule[MPR->ntypepair][1] = s2[0];
			MPR->ntypepair ++;
		}
	}
	return(MPR->ntypepair);
}


void read_matepair_rules_file(char *rulefile, MATEPAIRRULES *MPR)
{
  FILE *fp;
  fp = ckopen(rulefile, "r");
  readrules2(fp, MPR);
  fclose(fp);
}
