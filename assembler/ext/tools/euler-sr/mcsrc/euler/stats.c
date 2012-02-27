/***************************************************************************
 * Title:          stats.c
 * Author:         Haixu Tang
 * Created:        Jun. 2002
 * Last modified:  May. 2004
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdinc.h>

void STATS_addpt(STATS *s, int x)
{
  if (s->n == 0) {
    s->min_value = s->max_value = x;
  } else {
    s->min_value = min(s->min_value,x);
    s->max_value = max(s->max_value,x);
  }

  s->total += x;
  s->total2 += x*x;
  s->n ++;
}

/* should we use doubles instead? */
void STATS_compute(STATS *s,
		   int *val_min, int *val_max,
		   int *val_mean, int *val_sig)
{
  if (s->n < 2) {
    if (s->n < 1) {
      *val_min = *val_max = *val_mean = *val_sig = 0;
      return;
    }

    *val_min = s->min_value;
    *val_max = s->max_value;
    *val_mean = s->total;
    *val_sig = 0;
    return;
  }

    *val_min = s->min_value;
    *val_max = s->max_value;
    *val_mean = ((2 * s->total) / s->n + 1) / 2;  /* round */
    *val_sig = (int)(0.5 +
		     sqrt(((double) (s->total2 - ((double) (s->total * s->total) / (double) s->n)))/(s->n - 1)));
}
