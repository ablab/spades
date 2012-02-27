/***************************************************************************
 * Title:          readpair.c
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

/* TODO:
 * 1. Use MATEPAIRRULES structure as input instead of global variables
 * 2. Change fixed-size allocations to correct sizes
 * 3. Create new structure for names, pair1, pair2, etc.
 */

int readpair(char **names, int *pair1, int *pair2, int *min_dist, int *max_dist, int *type, FILE *fp);

int readpair(char **names, int *pair1, int *pair2, int *min_dist, int *max_dist, int *type, FILE *fp)
{
	int	i, j, k, l;
	char	str[100];

	i = 0;
	while(fgets(str, 90, fp))	{
		sscanf(str, "%s%d%d%d", names[i], &pair1[i], &pair2[i], &type[i]);
		min_dist[i] = pairrange[type[i]][0];
		max_dist[i] = pairrange[type[i]][1];
		i ++;
	}

	return(i);
}


/**********************************************************************
 * Input mate-pair information
 *
 * TODO: adapt allocation sizes to input
 **********************************************************************/

void read_matepair_file(char *pairfile,
			MATEPAIRRULES *MPR,
			READTABLE *RT,
			MATEPAIRTABLE *MP)
{
  FILE *fp;
  int i;
  char str[100];

  fp = ckopen(pairfile, "r");

  i = 0;
  while(fgets(str, 90, fp))	{
	i++;
  }
  MP->num_pair = i;
  MP->pair1 = (int *) ckalloc(MP->num_pair* sizeof(int));
  MP->pair2 = (int *) ckalloc(MP->num_pair* sizeof(int));
  MP->min_dist = (int *) ckalloc(MP->num_pair* sizeof(int));
  MP->max_dist = (int *) ckalloc(MP->num_pair* sizeof(int));
  MP->ntype = (int *) ckalloc(MP->num_pair* sizeof(int));
  MP->name = (char **) ckalloc(MP->num_pair* sizeof(char *));
  i = 0;
  rewind(fp);
  while(fgets(str, 90, fp))	{
    MP->name[i] = (char *) ckalloc(100* sizeof(char));
    sscanf(str, "%s%d%d%d", MP->name[i], &(MP->pair1[i]), &(MP->pair2[i]), &(MP->ntype[i]));
    MP->min_dist[i] = MPR->pairrange[MP->ntype[i]][0];
    MP->max_dist[i] = MPR->pairrange[MP->ntype[i]][1];
    i ++;
  }
  printf("Input pairs ...done; %d pairs are found\n", MP->num_pair);
  MP->num_sf_pair=0;
  MP->sf_pair = (int *) ckalloc(i * sizeof(int));
  fclose(fp);
}

/**********************************************************************
 * deallocate matepair table memory
 **********************************************************************/

void free_matepair(MATEPAIRTABLE *MP)
{
  int i;
  for(i = 0; i < MP->num_pair; i ++)	free((void *) MP->name[i]);
  free((void **) MP->name);
  free((void *) MP->pair1);
  free((void *) MP->pair2);
  free((void *) MP->sf_pair);
  free((void *) MP->ntype);
  free((void *) MP->min_dist);
  free((void *) MP->max_dist);
}
