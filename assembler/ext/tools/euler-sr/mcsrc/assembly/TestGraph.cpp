/***************************************************************************
 * Title:          TestGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include "graph/GraphAlgo.h"
#include "graph/MSTAlgo.h"
#include <string>
#include "IntegralTupleStatic.h"

int main(int argc, char* argv[] ) {
  std::string base;
  std::string graphFileName, intervalFileName, edgeFileName;
  if (argc != 2) {
    std::cout << "usage: testGraph graphBase " << std::endl;
    return 1;
  }
  base = argv[1];
  graphFileName = base + ".bgraph";
  intervalFileName = base + ".intv";
  edgeFileName = base + ".edge";
  IntervalGraph graph;
  graph.ReadIntervalGraph(graphFileName, intervalFileName);
  ReadSequences(edgeFileName, graph.edges);
  graph.vertexSize = 20;
  graph.MarkSuspectEdges(4);
	//  graph.FindAlternatePaths();
  return 0;
}
