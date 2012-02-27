/***************************************************************************
 * Title:          PrintDMST.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "graph/GraphAlgo.h"
#include "graph/MSTAlgo.h"
#include "IntegralTupleStatic.h"



int main(int argc, char* argv[]) {

  std::string baseName, graphName, intvName, outBase, graphOut, intvOut, gvzOutName;
  if (argc <2) {
    std::cout << "usage: printDMST graphBase outName graphVizOutName " << std::endl;
    exit(0);
  }
  baseName = argv[1];
  outBase  = argv[2];
  gvzOutName  = argv[3];
  graphName = baseName + ".bgraph";
  intvName  = baseName + ".intv";
  graphOut  = outBase + ".bgraph";
  intvOut   = outBase + ".intv";
  std::string edgeName  = baseName + ".edge";
  IntervalGraph graph;

  graph.ReadIntervalGraph(graphName, intvName);
  ReadSequences(edgeName, graph.edges);

  graph.CalcDMST();
  //  std::vector<std::vector<int> > sccs;
  //  FindStronglyConnectedComponents(graph.vertices, graph.edges, sccs);
  graph.PrintIntervalGraph(graphOut, intvOut);

  GVZPrintBGraph(graph.vertices, graph.edges, gvzOutName);
  return 0;
}
