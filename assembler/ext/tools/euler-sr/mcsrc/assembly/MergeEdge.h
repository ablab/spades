/***************************************************************************
 * Title:          MergeEdge.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  04/13/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MERGE_EDGE_H_
#define MERGE_EDGE_H_

#include "IntervalGraph.h"
#include "align/alignutils.h"

void MergeEdge(IntervalGraph &g,
							 TVertexList &vertices,
							 TEdgeList   &edges,
							 ssize_t vertexIndex,
							 ssize_t edge1Index, ssize_t edge2Index);



#endif
