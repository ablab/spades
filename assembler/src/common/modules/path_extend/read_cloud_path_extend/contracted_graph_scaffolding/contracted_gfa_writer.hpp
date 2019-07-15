#pragma once

#include <common/assembly_graph/contracted_graph/contracted_graph.hpp>
#include <common/io/utils/edge_namer.hpp>
namespace path_extend {
namespace read_cloud {
//todo remove duplication with GFAWriter
class ContractedGFAWriter {
 public:
  typedef contracted_graph::ContractedGraph ContractedGraph;
  typedef ContractedGraph::EdgeId EdgeId;

  ContractedGFAWriter(const ContractedGraph &graph, std::ostream &os,
                      io::EdgeNamingF<ContractedGraph> naming_f = io::IdNamingF<ContractedGraph>())
      : graph_(graph),
        edge_namer_(graph_, naming_f),
        os_(os) {
  }

  void WriteSegmentsAndLinks();

 private:
  void WriteSegments();
  void WriteLinks();
  void WriteLink(ContractedGFAWriter::EdgeId edge,
                 VertexId vertex,
                 size_t overlap_size,
                 std::ostream &os,
                 io::CanonicalEdgeHelper<ContractedGFAWriter::ContractedGraph> &namer);
  void WriteLink(VertexId vertex,
                 ContractedGFAWriter::EdgeId edge,
                 size_t overlap_size,
                 std::ostream &os,
                 io::CanonicalEdgeHelper<ContractedGFAWriter::ContractedGraph> &namer);
  void WriteSegment(const std::string &edge_id, const std::string &seq, double cov, std::ostream &os);

  const ContractedGraph &graph_;
  io::CanonicalEdgeHelper<ContractedGraph> edge_namer_;
  std::ostream &os_;
};

//class ContractedCSVWriter {
// public:
//    typedef contracted_graph::ContractedGraph ContractedGraph;
//    typedef ContractedGraph::EdgeId EdgeId;
//
//    void WriteProperties();
//
// private:
//    void WriteProperty(ContractedGraph::VertexId vertex);
//    void WriteProperty(ContractedGraph::EdgeId edge);
//
//    const ContractedGraph &graph_;
//    std::ostream &os_;
//};
}
}