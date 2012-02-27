/***************************************************************************
 * Title:          CalcGraphStatistics.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "graph/GraphAlgo.h"
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include "IntegralTupleStatic.h"


void PrintUsage() {
  std::cout << "usage: calcGraphStatistics bgraph_file interval_file " << std::endl;
  std::cout << "calcGraphStatistics - read in a condensed de Bruijn graph " << std::endl
	    << "  and output: directed/ indirected cycles statistics, " << std::endl
	    << "  and source/sink edges " << std::endl;
}

ssize_t CountZeroDegreeVertices(TVertexList &vertices) {
  ssize_t v;
  ssize_t numZero = 0;
  for (v = 0; v < vertices.size(); v++ ) {
    if (vertices[v].OutDegree() == 0 and vertices[v].InDegree() == 0)
      numZero++;
  }
  return numZero;
}

ssize_t CountSimpleBulges(TVertexList &vertices, TEdgeList &edges) {
  ssize_t v;
  ssize_t numSimpleBulges = 0;
  for (v = 0; v < vertices.size(); v++ ) {
    ssize_t a, b;

    for (a = 0; a < 3; a++ ) {
      for (b = a + 1; b < 4; b++ ) {
	if (vertices[v].out[a] != -1 and
	    vertices[v].out[b] != -1 and
	    edges[vertices[v].out[a]].dest == 
	    edges[vertices[v].out[b]].dest) {
	  /*	  std::cout << v << " " << edges[vertices[v].out[a]].dest << " "
		    << edges[vertices[v].out[a]].length << " "
		    << edges[vertices[v].out[b]].length << " "
		    << edges[vertices[v].out[a]].multiplicity << " "
		    << edges[vertices[v].out[b]].multiplicity << std::endl;
	  */
	  numSimpleBulges++;
	}
      }
    }
  }
  return numSimpleBulges;
}

template<typename V, typename E>
class CountSize {
public:
  ssize_t size;
  ssize_t length;
  ssize_t dest, src;
  std::vector<E> *edges;
  std::vector<V> *vertices;
  void operator()(ssize_t vertexIndex) {
    ssize_t e;
    size++;
    for (e = (*vertices)[vertexIndex].FirstIn(); 
	 e < (*vertices)[vertexIndex].EndIn(); 
	 e = (*vertices)[vertexIndex].NextIn(e)) {
      src = (*edges)[(*vertices)[vertexIndex].in[e]].src;
      if (!(*edges)[(*vertices)[vertexIndex].in[e]].IsMarked())
	length += (*edges)[(*vertices)[vertexIndex].in[e]].length;
    }
    for (e = (*vertices)[vertexIndex].FirstOut(); 
	 e < (*vertices)[vertexIndex].EndOut(); 
	 e = (*vertices)[vertexIndex].NextOut(e)) {
      dest = (*edges)[(*vertices)[vertexIndex].out[e]].dest;
      if (!(*edges)[(*vertices)[vertexIndex].out[e]].IsMarked()) 
	length += (*edges)[(*vertices)[vertexIndex].out[e]].length;
    }
  }
};


template<typename V, typename E>
ssize_t CountConnectedComponents(std::vector<V> &vertices, std::vector<E> &edges, std::vector<ssize_t> &sizes, std::vector<ssize_t> &lengths) {
  ssize_t numComponents = 0;
  ssize_t v;
  ssize_t size, length;
  Unmark(vertices);
  Unmark(edges);
  CountSize<V, E> countSize;
  countSize.edges = &edges;
  countSize.vertices = &vertices;
  for (v = 0; v < vertices.size(); v++ ) {
    if (! vertices[v].IsMarked()) {
      size = 0;
      length = 0;
      countSize.size = 0;
      countSize.length = 0;
      TraverseDFS(vertices, edges, v, countSize);
      std::cout << "got size " << countSize.size << " " << countSize.length << std::endl;
      if (countSize.size == 1) {
	std::cout << "vertex is size 1: ";
	WriteVertex(std::cout, vertices[v]);
	std::cout << std::endl;
      }
      sizes.push_back(countSize.size);
      lengths.push_back(countSize.length);
      numComponents++;
    }
  }
  return numComponents;
}


int main(int argc, char* argv[]) {

  IntervalGraph graph;
  std::string graphFileName, intervalFileName;
  int argi = 1;
  if (argc < 3) {
    PrintUsage();
    exit(1);
  }
  graphFileName = argv[argi++];
  intervalFileName = argv[argi++];
  graph.ReadIntervalGraph(graphFileName, intervalFileName);

  std::cout << "done reading (and checking) the graph " << std::endl;
  ssize_t numSimpleBulges;
  numSimpleBulges = CountSimpleBulges(graph.vertices, graph.edges);
  std::cout << "got " << numSimpleBulges << " bulges " << std::endl;

  ssize_t numComponents = 0;
  std::vector<ssize_t> sizes, lengths;
  numComponents = CountConnectedComponents(graph.vertices, graph.edges, sizes, lengths);
  std::cout << "counted " << numComponents << std::endl;

  ssize_t numZero;
  numZero = CountZeroDegreeVertices(graph.vertices);
  std::cout << "num zero: " << numZero << std::endl;
  /*
    ssize_t s;
    for (s = 0; s < sizes.size(); s++ ) {
    std::cout << sizes[s] << " " << lengths[s] << std::endl;
    }
  */
  
  return 0;
}
