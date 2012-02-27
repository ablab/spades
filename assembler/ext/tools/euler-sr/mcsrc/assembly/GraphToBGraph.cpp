/***************************************************************************
 * Title:          GraphToBGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "GraphReader.h"
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "DeBruijnGraph.h"
#include <iostream>
#include "IntegralTupleStatic.h"


int main(int argc, char* argv[]) {

  std::string graphName, bgraphName, edgeFileName;
  IntervalGraph graph;
  if (argc != 4) {
    std::cout << "usage: graphToBGraph graphfile edgefile bgraphfile" << std::endl;
    return 1;
  }

  graphName = argv[1];
  edgeFileName = argv[2];
  bgraphName = argv[3];
  ReadGraph(graphName, graph);
  SimpleSequenceList simpleSequences;
  ReadSimpleSequences(edgeFileName, simpleSequences);
  if (simpleSequences.size() != graph.edges.size()) {
    std::cout << "error, there should be the same number of edges " << std::endl
	      << " in the graph as in the edge file " << std::endl;
    exit(1);
  }

  ssize_t e;
  for (e = 0; e < graph.edges.size(); e++ ) {
    graph.edges[e].length = simpleSequences[e].length;
  }
  PrintBGraph(graph.vertices, graph.edges, graph.vertexSize, bgraphName);

  return 0;
}
