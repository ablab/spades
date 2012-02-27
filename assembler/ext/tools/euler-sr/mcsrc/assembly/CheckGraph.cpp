/***************************************************************************
 * Title:          CheckGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "BVertex.h"
#include "DeBruijnGraph.h"
#include "ReadIntervals.h"
#include <string>
#include "IntegralTupleStatic.h"


ssize_t CheckBalance(IntervalEdgeList &edges);
ssize_t CheckBalanceMultiplicity(IntervalEdgeList &edges);

int main(int argc, char* argv[]) {

  std::string graphFileName, intvFileName;
  if (argc < 2) {
    std::cout << "usage: checkGraph b_graphfile" << std::endl << std::endl;
    std::cout << "  This reads in a graph, and checks for consistency " << std::endl;
    std::cout << "  using the following checks: " << std::endl
	      << "  - balanced graph (make sure each edge is the balance of it's balanced edge" << std::endl;
    return 1;
  }
  intvFileName  = "";
  graphFileName = argv[1];
  int argi = 2;
  if (argi < argc)
    intvFileName = argv[argi++];
  
  BVertexList vertices;
  IntervalEdgeList   edges;

  ReadBGraph(graphFileName, vertices, edges);
  CheckBalance(edges);


  if (intvFileName != "") {
    ReadReadIntervals(intvFileName, edges);
    if (CheckBalanceMultiplicity(edges)) {
      std::cout << "invalid balanced multiplicity " << std::endl;
    }
  }
  return 0;
}

ssize_t CheckBalanceMultiplicity(IntervalEdgeList &edges) {

  ssize_t e;
  ssize_t balanced;
  ssize_t retval = 0;
  for (e = 0; e < edges.size(); e++ ) {
    balanced = edges[e].balancedEdge;
    std::cout << edges[balanced].multiplicity << " " <<  edges[e].multiplicity << std::endl;
    if (edges[balanced].multiplicity != edges[e].multiplicity) {
      std::cout << "error: edges " << e << " (" << edges[e].multiplicity << ") "
		<< balanced  << " (" << edges[balanced].multiplicity << ") " << std::endl;
      retval = 1;
    }
  }
  return retval;
}

ssize_t CheckBalance(IntervalEdgeList &edges) {
  ssize_t e;
  ssize_t balancedEdge;
  for (e = 0; e < edges.size(); e++) {
    balancedEdge = edges[e].balancedEdge;
    if (edges[balancedEdge].balancedEdge != e) {
      std::cout << e << " should be balanced with " << balancedEdge 
		<< " except it is " << edges[balancedEdge].balancedEdge 
		<< std::endl;
      return 0;
    }
  }
  return 1;
}
