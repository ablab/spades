/***************************************************************************
 * Title:          Printervals.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/03/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "IntervalGraph.h"
#include "ReadIntervals.h"
#include <iostream>
#include <vector>
#include "IntegralTupleStatic.h"

int main(int argc, char* argv[]) {

  std::string baseName;
	std::string graphFileName, intvFileName;
  IntervalGraph graph;
  //UNUSED// ssize_t numEdges;
  int argi = 1;
  int vertexSize;
  if (argc < 3) {
    std::cout << "usage: printervals graphBase vertexSize [edgesToPrint]" << std::endl;
    return 1;
  }
  baseName = argv[argi++];
	graphFileName = baseName + ".bgraph";
	intvFileName   = baseName + ".intv";
  vertexSize  = atoi(argv[argi++]);

	graph.vertexSize = vertexSize;

  graph.ReadIntervalGraph(graphFileName, intvFileName);
  
  std::vector<ssize_t> edgesToPrint;
  
  std::vector<ssize_t> revPath, forPath, fullPath;

  while (argi < argc) {
    edgesToPrint.push_back(atosz(argv[argi]));
    argi++;
  }

  ssize_t e, i, p;
  for (e = 0; e < edgesToPrint.size(); e++ ) {
    std::cout << "Intervals for edge " << edgesToPrint[e] << std::endl;
    assert(edgesToPrint[e] < graph.edges.size());

    for ( i = 0; i < graph.edges[edgesToPrint[e]].intervals->size(); i++ ) {
      revPath.clear();
      forPath.clear();
      graph.edges[edgesToPrint[e]].index = edgesToPrint[e];
      StoreIntervalPathReverse(graph.edges[edgesToPrint[e]], i,
															 graph.vertices, graph.edges,
															 vertexSize, revPath);
      
      StoreIntervalPathForward(graph.edges[edgesToPrint[e]], i,
															 graph.vertices, graph.edges,
															 vertexSize, forPath);
      std::cout << (*graph.edges[edgesToPrint[e]].intervals)[i].read << " ";
      for (p = revPath.size()-1; p >= 0; p-- )
				std::cout << " " << revPath[p];
      std::cout << " " << edgesToPrint[e];
      for (p = 0; p < forPath.size(); p++ )
				std::cout << " " << forPath[p];
      std::cout << std::endl;
    }
  }
}
