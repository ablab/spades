/***************************************************************************
 * Title:          RepeatSearch.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef FSDRTR_H_
#define FSDRTR_H_

void FindShortestDistanceBetweenTwoRepeats(TVertexList &vertices, TEdgeList &edges);
ssize_t SearchForRepeat(TVertexList &vertices, TEdgeList &edges, 
		    ssize_t curEdge, ssize_t curSearchLength, ssize_t maxSearchLength,
		    std::map<ssize_t, ssize_t> &distMap,
		    ssize_t &repeatEdge, ssize_t &len1, ssize_t &len2);


#endif
