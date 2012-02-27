/***************************************************************************
 * Title:          PrintSourceSinks.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "SeqReader.h"
#include "SimpleSequence.h"
#include <iostream>
#include <sstream>
#include "IntegralTupleStatic.h"

int main(int argc, char *argv[]) {

  std::string baseName, graphFileName, intvFileName, edgeFileName;
  std::string danglingEdgeFileName, danglingIntvFileName;
  if (argc != 5) {
    std::cout << "usage: printSourceSinks grap_base  danglingEdge danglingIntv vertexSize"
							<< std::endl;
		std::cout << "    Prints the source/sink edges to danglingEdge, and the corresponding intervals to danglingIntv." << endl;
    return 1;
  }
  baseName      = argv[1];
  danglingEdgeFileName = argv[2];
  danglingIntvFileName = argv[3];
	int vertexSize = atoi(argv[4]);
  graphFileName = baseName + ".bgraph";
  intvFileName  = baseName + ".intv";
  edgeFileName  = baseName + ".edge";
  IntervalGraph graph;
  graph.vertexSize = vertexSize;
  graph.ReadIntervalGraph(graphFileName, intvFileName);
  
  SimpleSequenceList edges;
  ReadSimpleSequences(edgeFileName, edges);
  
  std::ofstream dangEdgeOut, dangIntvOut;

  openck(danglingEdgeFileName, dangEdgeOut, std::ios::out);
  openck(danglingIntvFileName, dangIntvOut, std::ios::out);

  ssize_t e;
  std::stringstream titlestrm;
  std::string title;
  for (e = 0; e < graph.edges.size(); e++ ) {
    if (graph.vertices[graph.edges[e].src].InDegree() == 0 or
	graph.vertices[graph.edges[e].dest].OutDegree() == 0) {
      // found a dangling edge.
      title = "";
      titlestrm.str(title);
      titlestrm << e << " " << graph.edges[e].length;
      title = titlestrm.str();
      edges[e].PrintSeq(dangEdgeOut, title);

      // Now print the intervals that correspond to this edge
      dangIntvOut << ">" << title << std::endl;
      ssize_t intv;
      for (intv = 0; intv < graph.edges[e].intervals->size(); intv++) {
	dangIntvOut << (*graph.edges[e].intervals)[intv].read << " "
		    << (*graph.edges[e].intervals)[intv].readPos << std::endl;
      }
    }
  }
  dangEdgeOut.close();
  dangIntvOut.close();
  return 0;
}
