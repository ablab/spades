/***************************************************************************
 * Title:          PrintReadPaths.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "ReadIntervals.h"
#include "IntervalGraph.h"
#include "IntegralTupleStatic.h"


int main(int argc, char* argv[]) {
  std::string inGraph, inIntv, inReads, outPaths;
  if (argc != 5 ){
    std::cout << "usage: printReadPaths graph intv reads pathsOut" << std::endl;
    return 1;
  }
  inGraph = argv[1];
  inIntv  = argv[2];
  inReads = argv[3];
  outPaths = argv[4];
  IntervalGraph graph;
  graph.ReadIntervalGraph(inGraph, inIntv);
  graph.vertexSize = 20;
  SimpleSequenceList reads;
  ReadSimpleSequences(inReads, reads);
  

  std::vector<std::vector<ssize_t > > pathEdges, pathIndices;
  pathEdges.resize(reads.size()*2);
  pathIndices.resize(reads.size()*2);

  // cycle through all edges adding paths

  ssize_t e,i;
  //UNUSED// ssize_t prevEdge, curEdge;
  //UNUSED// ssize_t prevIndex, curIndex;
  //UNUSED// ssize_t nextEdge, nextIndex;
  ssize_t read;
  for (e = 0; e < graph.edges.size(); e++ ) {
    for (i= 0; i < graph.edges[e].intervals->size(); i++ ) {
      read = (*graph.edges[e].intervals)[i].read;
      if (pathEdges[read].size() > 0) 
	continue;
      // Look backwards for a read interval
    }
  }
  std::ofstream out;
  openck(outPaths, out, std::ios::out);
  ssize_t p;
  for (p = 0 ; p < reads.size()*2; p++ ) {
    out << p << " ";
    for (i = 0; i < pathEdges[p].size(); i++ ) {
      out << pathEdges[p][i] << " " << pathIndices[p][i] << ", ";
    }
    out << std::endl;
  }
  return 0;
}

