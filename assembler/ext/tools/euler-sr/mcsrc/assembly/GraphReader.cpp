/***************************************************************************
 * Title:          GraphReader.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "GraphReader.h"


ssize_t ReadGraph(std::string &graphFileName, IntervalGraph &graph) {
  //UNUSED+// ssize_t numEdges;
  ssize_t numVertices ;
  //  GetDimensions(graphFileName, numVertices, numEdges);
  //  std::cout << "resizing graph to: " << numVertices << " " << numEdges << std::endl;
  std::ifstream in;
  openck(graphFileName, in, std::ios::in);
  std::string word;

  ssize_t vertexIndex = -1;
  ssize_t balEdgeIndex;
  ssize_t edgeIndex;
  ssize_t nextEdge, lastEdge;
  ssize_t i;
  std::vector<ssize_t> balancedEdges;

  // Read in the graph
  while(in) {
    if (! (in >> word)) {
      break;
    }
    if (word == "Number_of_Vertex") {
      in >> numVertices;
      graph.vertices.resize(numVertices);
    }
    else if (word == "Vertex") {
      in >> vertexIndex >> nextEdge >> lastEdge;
    }
    else if (word == "Last_edge") {
      for (i = 0; i < lastEdge; i++) {
				in >> edgeIndex;
				graph.vertices[vertexIndex].in[i] = edgeIndex;
				in >> balEdgeIndex;
				balancedEdges.push_back(balEdgeIndex);
      }
    }
    else if (word == "Next_edge") {
      for (i = 0; i < nextEdge; i++ ) {
				in >> edgeIndex;
				graph.vertices[vertexIndex].out[i] = edgeIndex;
				in >> balEdgeIndex;
				balancedEdges.push_back(balEdgeIndex);
      }
    }
    else {
      std::cout << "Bad input: value: " << word << " not correct " << std::endl;
      exit(1);
    }
  }
  in.close();

  graph.edges.resize(balancedEdges.size()/2);
  std::cout << "allocating " << graph.edges.size() << " edges" << std::endl;
  // Now transform the edges 
  ssize_t v;
  ssize_t curBalancedEdge = 0;
  ssize_t balancedEdgeIndex;
  for (v = 0; v < graph.vertices.size(); v++ ) {
    // last edge is first
    for (i = 0; i < 4 and graph.vertices[v].in[i] >= 0; i++ ) {
      edgeIndex = graph.vertices[v].in[i];
      graph.edges[edgeIndex].src = v;
      balancedEdgeIndex = balancedEdges[curBalancedEdge];
      curBalancedEdge++;
      graph.edges[edgeIndex].balancedEdge = balancedEdgeIndex;
    }
    for (i = 0; i < 4 and graph.vertices[v].out[i] >= 0; i++ ) {
      edgeIndex = graph.vertices[v].out[i];
      graph.edges[edgeIndex].dest = v;
      balancedEdgeIndex = balancedEdges[curBalancedEdge];
      curBalancedEdge++;
      graph.edges[edgeIndex].balancedEdge = balancedEdgeIndex;
    }
  }
  // Do a sanity check on the edges
  ssize_t e;
  for (e = 0; e < graph.edges.size(); e++ ) {
    if (graph.edges[e].src == -1 or
				graph.edges[e].dest == -1) {
      std::cout << "error: edge " << e << " not assigned vertices " << std::endl;
      exit(1);
    }
    if (graph.edges[e].balancedEdge < 0 or 
				graph.edges[graph.edges[e].balancedEdge].balancedEdge != e) {
      std::cout << "error: edge " << e << " has a bad balanced edge " 
								<< graph.edges[e].balancedEdge << " and reflected: " 
                << graph.edges[graph.edges[e].balancedEdge].balancedEdge 
                << std::endl;
      exit(1);
    }
  }
  return numVertices;
}
