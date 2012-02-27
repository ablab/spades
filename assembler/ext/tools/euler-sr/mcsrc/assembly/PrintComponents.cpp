/***************************************************************************
 * Title:          PrintComponents.cpp 
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
#include <string>


int main(int argc, char* argv[]) {
  std::string graphFile;
  
  if (argc < 2) {
    std::cout << "usage: printComponents file.bgraph" << std::endl;
    exit(1);
  }
  graphFile = argv[1];

  TVertexList vertices;
  TEdgeList edges;
  ReadBGraph(graphFile, vertices, edges);

  PrintComponents(vertices, edges);

  return 0;
}
