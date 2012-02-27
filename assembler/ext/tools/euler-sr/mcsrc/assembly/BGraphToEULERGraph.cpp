/***************************************************************************
 * Title:          BGraphToEULERGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "BVertex.h"
#include "BEdge.h"
#include <string>
#include "IntegralTupleStatic.h"


int main(int argc, char* argv[]) {
  std::string bGraphFileName, eulerGraphFileName;
  std::string eulerBaseFileName, eulerEdgeFileName;
  if (argc < 3) {
    std::cout << "usage: bGraphToEuler bgraphIn eulerOut(base) " << std::endl;
    std::cout << "  take a Bgraph and output in the format that euler_et reads in " << std::endl;
    return 1;
  }

  bGraphFileName = argv[1];
  eulerBaseFileName = argv[2];
  eulerGraphFileName = eulerBaseFileName + ".graph";
  eulerEdgeFileName = eulerBaseFileName + ".edge";
  BVertexList vertices;
  BEdgeList edges;
  ReadBGraph(bGraphFileName, vertices, edges);
  
  PrintGraph(vertices, edges, eulerGraphFileName);
  PrintEdges(vertices, edges, eulerEdgeFileName);

  return 0;
}

