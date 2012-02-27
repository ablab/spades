/***************************************************************************
 * Title:          PrintGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "PrintGraph.h"

void PrintGraphSubset(IntervalGraph &g, 
											std::set<ssize_t> &vertexSet,
											ssize_t minEdgeLength,
											std::ostream &graphOut) {
	
	graphOut << "digraph G {"<<std::endl
					 << "\tsize=\"8,8\";"<<std::endl;

	std::set<ssize_t>::iterator vIt, vEnd;
	vIt = vertexSet.begin(); vEnd = vertexSet.end();
	ssize_t inEdgeIndex, inEdge, outEdgeIndex, outEdge;

	for (; vIt != vEnd; ++vIt) {
		for (inEdgeIndex = g.vertices[*vIt].FirstIn();
				 inEdgeIndex != g.vertices[*vIt].EndIn();
				 inEdgeIndex =  g.vertices[*vIt].NextIn(inEdgeIndex)) {
			inEdge = g.vertices[*vIt].in[inEdgeIndex];
			if (vertexSet.find(g.edges[inEdge].src) == vertexSet.end() and
					g.edges[inEdge].length > minEdgeLength) {
				graphOut << "   " << g.edges[inEdge].src << " -> " << *vIt << " [style=bold, color=red, label=\""
								 << g.edges[inEdge].index << "  " << g.edges[inEdge].length  << "\"];" << std::endl;
			}
		}
		for (outEdgeIndex = g.vertices[*vIt].FirstOut();
				 outEdgeIndex != g.vertices[*vIt].EndOut();
				 outEdgeIndex =  g.vertices[*vIt].NextOut(outEdgeIndex)) {
			outEdge = g.vertices[*vIt].out[outEdgeIndex];
			if (vertexSet.find(g.edges[outEdge].dest) != vEnd or 
					g.edges[outEdge].length > minEdgeLength) {
				graphOut << "   " << *vIt << " -> " << g.edges[outEdge].dest;
				
				if (vertexSet.find(g.edges[outEdge].dest) == vEnd) {
					graphOut << " [style=bold,color=red ";
				}
				else {
					graphOut << " [";
				}
				graphOut << " label=\""
								 << g.edges[outEdge].index << "  " << g.edges[outEdge].length  << "\"];" << std::endl;
			}
		}
	}
	graphOut << "}";
}
