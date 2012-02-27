/***************************************************************************
 * Title:          zigzag.h
 * Author:         Glenn Tesler
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
/* zigzag paths
 * Glenn Tesler
 * 12/26/02
 */

void process_all_zzpaths(IGRAPH *G, READTABLE *RT);

void create_test_zzpaths(IGRAPH *G, READTABLE *RT,
			 int zz_type, char *u, int L, int count0);

