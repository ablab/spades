/***************************************************************************
 * Title:          testShortCycle.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "IntervalGraph.h"
#include "DeBruijnGraph.h"
#include "ReadIntervals.h"
#include <string>

int main(int argc, char* argv[]) {

  std::string graphFileName, intvFileName, edgeFileName;
  if (argc < 6) {
    std::cout << "usage: checkGraph b_graph intv edge source cyclesize" << std::endl << std::endl;
    std::cout << "  This reads in a graph, and checks for consistency " << std::endl;
    std::cout << "  using the following checks: " << std::endl
	      << "  - balanced graph (make sure each edge is the balance of it's balanced edge" << std::endl;
    return 1;
  }
  ssize_t source, cycleSize;
  intvFileName  = "";
  graphFileName = argv[1];
  intvFileName  = argv[2];
  edgeFileName  = argv[3];
  source = atosz(argv[4]);
  cycleSize = atoi(argv[5]);
  IntervalGraph graph;
  graph.vertexSize = 20;
  graph.ReadIntervalGraph(graphFileName, intvFileName);
  ReadSequences(edgeFileName, graph.edges);
  std::set<ssize_t> dagVertices, dagEdges;
  std::vector<ssize_t> optimalPath;
  ssize_t out;
  graph.Unmark();
  graph.MarkSuspectEdges(5);
	/*
  if (graph.FindShortContainedDAG(source, cycleSize, dagVertices, dagEdges, optimalPath) ) {
    std::cout << "found a short bulge from " << source << " out: " << out << std::endl;

    ssize_t sinkVertex = graph.edges[out].src;
    std::vector<ssize_t> path;
    //    if (graph.FindMaximumPath(source, sinkVertex, path)) {
    std::cout << "found a path of len: " << path.size() << std::endl;
    ssize_t i;
    for (i =0; i < optimalPath.size(); i++) {
      std::cout << optimalPath[i] << std::endl;
    }
  }
	
  else {
    std::cout << "no bulgex found " << std::endl;
  }
	*/
  std::string graphOutName = graphFileName + ".test.dot";
  std::cout << "printing gvz graph asdf"<< std::endl; 
  GVZPrintBGraph(graph.vertices, graph.edges, graphOutName);
  return 0;
}
