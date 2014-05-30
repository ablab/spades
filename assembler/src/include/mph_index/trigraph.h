//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef __CXXMPH_TRIGRAPH_H__
#define __CXXMPH_TRIGRAPH_H__
// Build a trigraph using a memory efficient representation.
//
// Prior knowledge of the number of edges and vertices for the graph is
// required. For each vertex, we store how many edges touch it (degree) and the
// index of the first edge in the vector of triples representing the edges.

#include <stdint.h>  // for uint32_t and friends

#include <vector>

namespace cxxmph {

class TriGraph {
 public:
  struct Edge {
    Edge() { }
    Edge(uint32_t v0, uint32_t v1, uint32_t v2) { 
      vertices[0] = v0;
      vertices[1] = v1;
      vertices[2] = v2;
    }
    uint32_t& operator[](uint8_t v) { return vertices[v]; }
    const uint32_t& operator[](uint8_t v) const { return vertices[v]; }
    uint32_t vertices[3];
  };
  TriGraph(uint32_t nedges, uint32_t nvertices);
  ~TriGraph();

  void AddEdge(const Edge& edge);
  void RemoveEdge(uint32_t edge_id);
  void ExtractEdgesAndClear(std::vector<Edge>* edges);
  void DebugGraph() const;

  const std::vector<Edge>& edges() const { return edges_; }
  uint8_t vertex_degree(size_t idx) const {
    return vertex_degree_[idx];
  }
  const std::vector<uint32_t>& first_edge() const { return first_edge_; }

 private:
  uint32_t nedges_;  // total number of edges
  std::vector<Edge> edges_;
  std::vector<Edge> next_edge_;  // for implementing removal
  std::vector<uint32_t> first_edge_;  // the first edge for this vertex
  std::vector<uint8_t> vertex_degree_;  // number of edges for this vertex
};

}  // namespace cxxmph

#endif  // __CXXMPH_TRIGRAPH_H__
