/***************************************************************************
 * Title:          sf_sub.c
 * Author:         Glenn Tesler
 * Created:        Mar. 2003
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/*
 * Subroutines for maintaining SFPAIRS structure
 * Haixu Tang, Glenn Tesler
 * 3/27/03
 *
 * TODO if the # SF pairs gets to be large:
 *  1. instead of storing (e1,e2) and (e2^c,e1^c),
 *     just store one of them (min. of e1 and e2^c)
 *  2. it's in no particular order now, store it in
 *     a sorted order
 */

#include <stdinc.h>
#include <extfunc.h>
#include <stats.h>

int newsfpair(EDGE *begin_edge, EDGE *end_edge, SFPAIRS *SFP)
{
  int	i;

  for (i = 0; i < SFP->num_sf; i ++) {
    if (SFP->linkedge1[i] == begin_edge && SFP->linkedge2[i] == end_edge) {
      SFP->link_cov[i] ++;
      if (i % 2 == 1) {
	SFP->link_cov[i - 1] ++;
      } else {
	SFP->link_cov[i + 1] ++;
      }
      return(SFP->num_sf);
    }
  }
  SFP->linkedge1[SFP->num_sf] = begin_edge;
  SFP->linkedge2[SFP->num_sf] = end_edge;
  SFP->link_cov[SFP->num_sf ++] = 1;
  SFP->linkedge1[SFP->num_sf] = end_edge -> bal_edge;
  SFP->linkedge2[SFP->num_sf] = begin_edge -> bal_edge;
  SFP->link_cov[SFP->num_sf ++] = 1;
  return(SFP->num_sf);
}



int newsfpair_w_stats(EDGE *begin_edge, EDGE *end_edge,
		      int lo, int hi,
		      SFPAIRS *SFP)
{
  int	i;

  for (i = 0; i < SFP->num_sf; i ++) {
    if (SFP->linkedge1[i] == begin_edge && SFP->linkedge2[i] == end_edge) {
      break;
    }
  }

  if (i == SFP->num_sf) {
    /* new pair */
    SFP->linkedge1[i] = begin_edge;
    SFP->linkedge2[i] = end_edge;
    SFP->link_cov[i] = 1;
    SFP->linkedge1[i+1] = end_edge -> bal_edge;
    SFP->linkedge2[i+1] = begin_edge -> bal_edge;
    SFP->link_cov[i+1] = 1;
    SFP->num_sf += 2;
  } else {
    /* old pair */
    i = i & ~1;
    SFP->link_cov[i]++;
    SFP->link_cov[i+1]++;
  }

  /* update additional statistics */
  STATS_addpt(&SFP->range_lo[i], lo);
  STATS_addpt(&SFP->range_hi[i], hi);
  STATS_addpt(&SFP->range_lo[i+1], lo);
  STATS_addpt(&SFP->range_hi[i+1], hi);

  return(SFP->num_sf);
}




/* init_sfpairs(SFP, MP, 0)
 *     initialize an SFP structure for storing scaffolding edge pairs
 *     but no length statistics
 * init_sfpairs(SFP, MP, 1)
 *     initialize an SFP structure for storing scaffolding edge pairs
 *     and also length statistics
 */
void init_sfpairs(SFPAIRS *SFP, MATEPAIRTABLE *MP, int has_len)
{
  SFP->num_sf = 0;
  SFP->tot_edge = 0;
  SFP->linkedge1 = (EDGE **) ckalloc(2 * MP->num_pair * sizeof(EDGE *));
  SFP->linkedge2 = (EDGE **) ckalloc(2 * MP->num_pair * sizeof(EDGE *));
  SFP->link_cov = (int *) ckalloc(2 * MP->num_pair * sizeof(int));

  if (has_len) {
    SFP->range_lo = (STATS *) ckalloc(2 * MP->num_pair * sizeof(STATS));
    SFP->range_hi = (STATS *) ckalloc(2 * MP->num_pair * sizeof(STATS));
  } else {
    SFP->range_lo = (STATS *) 0;
    SFP->range_hi = (STATS *) 0;
  }
}

void free_sfpairs(SFPAIRS *SFP)
{
  if(SFP->tot_edge>0) free((void *) SFP->all_edge);
  free((void *) SFP->linkedge1);
  free((void *) SFP->linkedge2);
  free((void *) SFP->link_cov);
  free((void *) SFP->range_lo);
  free((void *) SFP->range_hi);
}
