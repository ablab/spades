#include "graph.hpp"

void CVertex::AddEdge(CEdge e) {
//	edgeArray.push_back(e);
}

void CGraph::AddVertex(CVertex v) {
//	vertexArray.push_back(v);
}

extern ostream & operator << ( ostream & os, const CGraph & g )
{
  os << "\nGraph: ";
  return os;
}
