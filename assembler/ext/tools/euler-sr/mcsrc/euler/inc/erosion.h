/***************************************************************************
 * Title:          erosion.h
 * Author:         Glenn Tesler
 * Created:        Jun. 2002
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2001-2004 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _EROSION_H_
#define _EROSION_H_



void erode_graph(IGRAPH *G, CLEAN_PARAMS *params);
//void label_external(IGRAPH *G, int T_e);
//void classify_internal_nodes(IGRAPH *G, int depth);
void label_external(IGRAPH *G, LABELVERT *T_e);
void classify_internal_nodes(IGRAPH *G, LABELVERT *T_e);
void classify_internal_components(IGRAPH *G);
int is_path_internal(NODES *v, int e_type, int L_c);
int is_edge_chimeric(EDGE *e, int L_c);






#endif /* _EROSION_H_ */

