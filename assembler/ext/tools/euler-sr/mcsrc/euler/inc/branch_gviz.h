/***************************************************************************
 * Title:          branch_gviz.h
 * Author:         Glenn Tesler
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

/* Output branching graph of H in .gviz format
 * Leaves G and H intact
 * Glenn Tesler
 * 12/19/02
 */

void write_bgraph_gviz_file(IGRAPH *G, char *filename);
void write_bgraph_gviz(IGRAPH *G, FILE *fp);

