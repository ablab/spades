/***************************************************************************
 * Title:          RepeatResolution.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <map>
#include <set>
#include "RepeatResolution.h"


ssize_t IsRepeatEdgeBalanced(IntervalGraph &g, ssize_t edge) {
	ssize_t srcVertex, destVertex;
	srcVertex = g.edges[edge].src;
	destVertex = g.edges[edge].dest;
	
	return g.vertices[srcVertex].InDegree() == g.vertices[destVertex].OutDegree();
}



ssize_t ConsistentRepeatEdgePaths(IntervalGraph &g, ssize_t edge) {
	
	ssize_t i;
	ssize_t src, dest;
	src = g.edges[e].src;
	dest = g.edges[e].dest;
	ReadIntervalList *intervals;
	intervals = g.edges[edge].intervals;
	ssize_t p, pi;

	std::vector<ssize_t> outPair, uniqueOutPair;
	outPair.resize(g.vertices[dest].OutDegree());
	uniqueOutPair.resize(g.vertices[dest].OutDegree());
	
	ssize_t outEdge;
	for (outEdge = 0; outEdge < uniqueOutPair.size(); outEdge++) {

	}
	
	if (g.vertices[src].InDegree() > 1) {

		
		ssize_t inEdge, inEdgeIndex;

		for (inEdgeIndex = g.vertices[src].FirstIn();
				 inEdgeIndex != g.vertices[src].EndIn();
				 inEdgeIndex = g.vertices[src].NextIn(inEdgeIndex)) {


			inEdge = g.vertices[src].in[inEdgeIndex];

			for (i = 0; i < intervals->size(); i++ ) {
				p  = (*intervals)[i].readIndex;
				pi = (*intervals)[i].pathPos;
				if (pi > 0 and graph.paths[p][pi].edge == inEdge) {
					if (pi < g.pathLengths[p] - 1) {
						
		
		
	}
	



}
