/***************************************************************************
 * Title:          stats.h
 * Author:         Glenn Tesler
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
void STATS_addpt(STATS *s, int x);
void STATS_compute(STATS *s,
		   int *val_min, int *val_max,
		   int *val_mean, int *val_sig);


