#include <cassert>
#include <limits>
#include <iostream>

#include "mph_index/trigraph.h"

using std::cerr;
using std::endl;
using std::vector;

namespace {
static const uint32_t kInvalidEdge = std::numeric_limits<uint32_t>::max();
} 

namespace cxxmph {

TriGraph::TriGraph(uint32_t nvertices, uint32_t nedges)
      : nedges_(0),
        edges_(nedges),
        next_edge_(nedges),
        first_edge_(nvertices, kInvalidEdge),
        vertex_degree_(nvertices, 0) { }
TriGraph::~TriGraph() {}

void TriGraph::ExtractEdgesAndClear(vector<Edge>* edges) {
  vector<Edge>().swap(next_edge_);
  vector<uint32_t>().swap(first_edge_);
  vector<uint8_t>().swap(vertex_degree_);
  nedges_ = 0;
  edges->swap(edges_);
}

void TriGraph::AddEdge(const Edge& edge) { 
  edges_[nedges_] = edge; 
  assert(first_edge_.size() > edge[0]);
  assert(first_edge_.size() > edge[1]);
  assert(first_edge_.size() > edge[0]);
  assert(first_edge_.size() > edge[1]);
  assert(first_edge_.size() > edge[2]);
  assert(next_edge_.size() > nedges_);
  next_edge_[nedges_] = Edge(first_edge_[edge[0]],
                             first_edge_[edge[1]],
                             first_edge_[edge[2]]);
  first_edge_[edge[0]] = first_edge_[edge[1]] = first_edge_[edge[2]] = nedges_;
  ++vertex_degree_[edge[0]];
  ++vertex_degree_[edge[1]];
  ++vertex_degree_[edge[2]];
  ++nedges_;
}

void TriGraph::RemoveEdge(uint32_t current_edge) {
  // cerr << "Removing edge " << current_edge << " from " << nedges_ << " existing edges " << endl;
  for (uint8_t i = 0; i < 3; ++i) {
    uint32_t vertex = edges_[current_edge][i];
    uint32_t edge1 = first_edge_[vertex];
    uint32_t edge2 = kInvalidEdge;
    uint8_t j = 0;
    while (edge1 != current_edge && edge1 != kInvalidEdge) {
      edge2 = edge1;
      if (edges_[edge1][0] == vertex) j = 0;
      else if (edges_[edge1][1] == vertex) j = 1;
      else j = 2;
      edge1 = next_edge_[edge1][j];
    }
    assert(edge1 != kInvalidEdge);
    if (edge2 != kInvalidEdge) next_edge_[edge2][j] = next_edge_[edge1][i];
    else first_edge_[vertex] = next_edge_[edge1][i];
    --vertex_degree_[vertex];
  }
}

void TriGraph::DebugGraph() const {
  uint32_t i;
  for(i = 0; i < edges_.size(); i++){
    cerr << i << "  " << edges_[i][0] << " " << edges_[i][1] << " " << edges_[i][2]
         << " nexts " << next_edge_[i][0] << " " << next_edge_[i][1] << " " << next_edge_[i][2] << endl;
  }
  for(i = 0; i < first_edge_.size();i++){
    cerr << "first for vertice " <<i << " " << first_edge_[i] << endl;
  }
}

     
}  // namespace cxxmph
