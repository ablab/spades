/***************************************************************************
 * Title:          RepeatSearch.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "DeBruijnGraph.h"
#include "RepeatSearch.h"
#include <map>
#include <algorithm>

void FindShortestDistanceBetweenTwoRepeats(TVertexList &vertices, TEdgeList &edges) {

  ssize_t e;
  std::map<ssize_t, ssize_t> lengthToEdge;
  ssize_t repeatEdge, pathLen1, pathLen2;
  std::vector<ssize_t> distToRepeat, repeatEdgeList;
  distToRepeat.resize(edges.size());
	repeatEdgeList.resize(edges.size());
	
  std::fill(distToRepeat.begin(), distToRepeat.end(), 0);
  std::fill(repeatEdgeList.begin(), repeatEdgeList.end(), 0);

  for (e = 0; e < edges.size(); e++ ) {
    if (vertices[edges[e].src].InDegree() > 1) {
      lengthToEdge.clear();
      if (SearchForRepeat(vertices, edges,
			  e, 0, 10000,
			  lengthToEdge, repeatEdge, pathLen1, pathLen2)) {
				distToRepeat[e] = std::max(pathLen1, pathLen2) + edges[repeatEdge].length;
				repeatEdgeList[e] = repeatEdge;
      }
    }
  }

  std::cout << "distances that were found: " << std::endl;
  for (e = 0; e < distToRepeat.size(); e++ ) {
    if (distToRepeat[e] > 0)
      std::cout << "e: " << e << " " << repeatEdgeList[e] 
								<< " " << distToRepeat[e] << std::endl;
		
  }
}


ssize_t SearchForRepeat(TVertexList &vertices, TEdgeList &edges, 
		    ssize_t curEdge, ssize_t curSearchLength, ssize_t maxSearchLength,
		    std::map<ssize_t, ssize_t> &distMap,
		    ssize_t &repeatEdge, ssize_t &len1, ssize_t &len2) {
  if (distMap.find(curEdge) != distMap.end()) {
    // found a repeat
    len1 = distMap[curEdge];
    len2 = curSearchLength;
    repeatEdge = curEdge;
    return 1;
  }
  else {
    // record how long it took to get here
    distMap[curEdge] = curSearchLength;

    // traverse this edge
    curSearchLength += edges[curEdge].length;
    if (curSearchLength > maxSearchLength) 
      return 0;

    ssize_t dest = edges[curEdge].dest;

    //UNUSED+// ssize_t  outEdge;
    ssize_t outEdgeIndex;
    for (outEdgeIndex = vertices[dest].FirstOut();
	 outEdgeIndex < vertices[dest].EndOut();
	 outEdgeIndex = vertices[dest].NextOut(outEdgeIndex)) {
      if (SearchForRepeat(vertices, edges, 
			  vertices[dest].out[outEdgeIndex],
			  curSearchLength, maxSearchLength, 
			  distMap, repeatEdge, len1, len2))
	return 1;
    }

    // no repeat found
    return 0;
  }

}
